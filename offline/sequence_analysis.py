#!/usr/bin/env python3
"""Sequence-only conservation and SCA scores for comparison with residue graphs.

Implements the unregularized Frobenius SCA matrix described by Rivoire,
Reynolds & Ranganathan (2016), DOI:10.1371/journal.pcbi.1004817. The matrix
is evaluated in blocks to avoid allocating a (20L)^2 tensor. Continuous
eigenvector scores are exploratory; they are NOT inferred protein sectors.
No structural score is used to fit the sequence model.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse, stats

AA = 'ACDEFGHIKLMNPQRSTVWY'
BACKGROUND = np.array([.073,.025,.050,.061,.042,.072,.023,.053,.064,.089,
                       .023,.043,.052,.040,.052,.073,.056,.063,.013,.033])


def read_fasta(path):
    records = {}
    name = None
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            name = line[1:].split()[0]
            if name in records:
                raise ValueError(f'Duplicate FASTA identifier: {name}')
            records[name] = ''
        elif line.strip():
            if name is None:
                raise ValueError('Sequence before FASTA header')
            records[name] += line.strip().upper().replace('.', '-')
    if not records or len({len(v) for v in records.values()}) != 1:
        raise ValueError('An existing rectangular alignment is required')
    return records


def encode(sequences):
    lookup = np.full(256, 20, dtype=np.uint8)
    for i, letter in enumerate(AA):
        lookup[ord(letter)] = i
    a = np.frombuffer(''.join(sequences).encode('ascii'), dtype=np.uint8)
    return lookup[a].reshape(len(sequences), len(sequences[0]))


def one_hot(a, states=20):
    n, length = a.shape
    row, pos = np.nonzero(a < states)
    col = pos * states + a[row, pos]
    return sparse.csr_matrix((np.ones(len(row)), (row, col)),
                             shape=(n, length * states))


def sequence_weights(a, identity=.8, block=256):
    """Inverse neighbor count at strictly >identity, including gap matches.

    This matches pySCA's default convention. It reduces redundancy but does
    not remove all effects of shared ancestry. Exact duplicates are collapsed
    by the caller before weighting; this policy is recorded in the output.
    """
    if not 0 <= identity < 1 or not a.size:
        raise ValueError('Require a nonempty alignment and 0 <= identity < 1')
    x = one_hot(a, 21)
    counts = np.empty(len(a), dtype=int)
    for start in range(0, len(a), block):
        similarity = (x[start:start+block] @ x.T).toarray() / a.shape[1]
        counts[start:start+block] = (similarity > identity).sum(axis=1)
    if np.any(counts < 1):
        raise ValueError('Every sequence must be its own neighbor')
    return 1 / counts


def sca_matrix(a, weights, block=24):
    """Return Frobenius SCA matrix and single-site frequencies (no pseudocount)."""
    n, length = a.shape
    w = np.asarray(weights, dtype=float)
    if w.shape != (n,) or np.any(w <= 0) or not np.isfinite(w).all():
        raise ValueError('One finite positive weight per sequence is required')
    w = w / w.sum()
    x = one_hot(a)
    f = np.asarray(x.T @ w).reshape(length, 20)
    f = np.clip(f, 0, 1)
    q = np.broadcast_to(BACKGROUND, f.shape)
    position_weight = np.zeros_like(f)
    valid = (f > 0) & (f < 1)
    position_weight[valid] = abs(np.log(f[valid] * (1-q[valid]) /
                                       ((1-f[valid]) * q[valid])))
    xf = x.multiply(np.sqrt(w)[:, None]).tocsr()
    matrix = np.empty((length, length))
    for start in range(0, length, block):
        end = min(length, start+block)
        joint = (xf[:,start*20:end*20].T @ xf).toarray()
        covariance = joint - np.outer(f[start:end].ravel(), f.ravel())
        covariance *= np.outer(position_weight[start:end].ravel(),
                               position_weight.ravel())
        matrix[start:end] = np.sqrt(np.sum(
            covariance.reshape(end-start,20,length,20)**2, axis=(1,3)))
    if not np.allclose(matrix, matrix.T, atol=1e-10):
        raise ValueError('SCA matrix is not symmetric')
    return (matrix+matrix.T)/2, f


def single_site_scores(f):
    occupied = f.sum(axis=1)
    p = np.divide(f, occupied[:,None], out=np.zeros_like(f),
                  where=occupied[:,None] > 0)
    h = -np.sum(np.where(p > 0, p*np.log(np.maximum(p,1e-300)), 0), axis=1)
    theta = max(0., 1 - occupied.mean())
    q = (1-theta)*BACKGROUND
    kl = np.sum(np.where(f > 0, f*np.log(np.maximum(f,1e-300)/q), 0), axis=1)
    gaps = np.maximum(0.,1-occupied)
    if theta > 1e-12:
        kl += np.where(gaps > 1e-12, gaps*np.log(np.maximum(gaps,1e-300)/theta), 0)
    return dict(conservation=1-h/np.log(20), entropy_nats=h,
                relative_entropy=kl, weighted_gap_fraction=gaps)


def analyze_alignment(path, out, reference_ids=(), identity=.8,
                      max_sequence_gap=.2, max_column_gap=.2):
    """Fit sequence scores; references map columns but are excluded from fitting."""
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    records = read_fasta(path)
    missing = set(reference_ids)-records.keys()
    if missing:
        raise ValueError(f'Reference identifiers not found: {missing}')
    all_ids = list(records)
    all_a = encode(list(records.values()))
    ref_indices = [all_ids.index(k) for k in reference_ids]
    columns = np.flatnonzero(np.any(all_a[ref_indices] < 20, axis=0)) if ref_indices else np.arange(all_a.shape[1])
    train = [i for i,k in enumerate(all_ids) if k not in reference_ids]
    a = all_a[np.ix_(train,columns)]
    ids = np.array(all_ids)[train]
    valid = (a == 20).mean(axis=1) <= max_sequence_gap
    sequence_rejected = int((~valid).sum())
    a, ids = a[valid], ids[valid]
    if len(a) < 3:
        raise ValueError('Fewer than three sequences after sequence-gap filtering')
    _, unique = np.unique(a, axis=0, return_index=True)
    unique = np.sort(unique)
    duplicates = len(a)-len(unique)
    a, ids = a[unique], ids[unique]
    weights = sequence_weights(a, identity)
    keep = np.average(a==20,axis=0,weights=weights) <= max_column_gap
    a, columns = a[:,keep], columns[keep]
    if a.shape[1] < 2:
        raise ValueError('Fewer than two columns after gap filtering')
    weights = sequence_weights(a, identity)
    print(f'{out.name}: {len(a)} unique sequences, {len(columns)} columns, N_eff={weights.sum():.2f}', flush=True)
    matrix, f = sca_matrix(a, weights)
    vals, vecs = np.linalg.eigh(matrix)
    order = np.argsort(vals)[::-1]
    vals, vecs = vals[order], vecs[:,order]
    k = min(6,len(columns))
    scores = pd.DataFrame(dict(msa_column=columns+1, **single_site_scores(f),
        sca_leading_loading=np.abs(vecs[:,0]),
        sca_top6_loading=np.sqrt(np.sum(vecs[:,:k]**2,axis=1)),
        sca_strength=matrix.sum(axis=1)-matrix.diagonal()))
    scores.to_csv(out/'sequence_scores.csv',index=False)
    np.savez_compressed(out/'sca_matrix.npz',matrix=matrix,msa_columns=columns+1,
                        eigenvalues=vals, eigenvectors=vecs,weights=weights)
    pd.DataFrame(dict(sequence_id=ids,weight=weights)).to_csv(out/'sequence_weights.csv',index=False)
    mapping=[]
    included=set(columns)
    for ref in reference_ids:
        pos=0
        for col,aa in enumerate(records[ref]):
            if aa in AA:
                pos += 1
                mapping.append(dict(reference_id=ref,reference_position=pos,
                    residue=aa,msa_column=col+1,scored=col in included))
    pd.DataFrame(mapping).to_csv(out/'reference_mapping.csv',index=False)
    metadata=dict(input_sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest(),
        raw_sequences=len(records),reference_ids=list(reference_ids),
        references_excluded_from_fit=True,sequence_gap_rejected=sequence_rejected,
        exact_duplicates_removed=duplicates,fit_sequences=len(a),scored_columns=len(columns),
        effective_sequences=float(weights.sum()),effective_sequences_per_column=float(weights.sum()/len(columns)),
        max_sequence_gap=max_sequence_gap,column_gap_filter_threshold=max_column_gap,
        column_filter_policy='One pass using pre-column-filter weights; weights are recomputed after filtering.',
        final_max_weighted_column_gap=float(single_site_scores(f)['weighted_gap_fraction'].max()),
        identity_threshold=identity,weight_rule='1/count(identity > threshold); gaps count as matching state',
        gap_fraction=float(np.average((a==20).mean(axis=1),weights=weights)),
        method='unregularized Frobenius SCA; pySCA background; no structure in fitting',
        score_note='Leading and top-six eigenvector loadings are exploratory continuous scores, not inferred sectors; k=6 fixed without using network outcomes.',
        depth_note='Effective depth is a redundancy diagnostic, not a guarantee of independence or adequate phylogenetic sampling.',
        low_depth_flag=bool(weights.sum()<len(columns)))
    (out/'summary.json').write_text(json.dumps(metadata,indent=2)+'\n')
    return scores,metadata


def partial_spearman(x,y,controls):
    """Correlation of rank residuals; no independence-based p-value asserted."""
    a=np.column_stack([x,y,np.asarray(controls)])
    a=a[np.isfinite(a).all(axis=1)]
    if len(a)<5: return float('nan')
    ranks=np.apply_along_axis(stats.rankdata,0,a)
    design=np.column_stack([np.ones(len(a)),ranks[:,2:]])
    residual=ranks[:,:2]-design@np.linalg.lstsq(design,ranks[:,:2],rcond=None)[0]
    return float(np.corrcoef(residual.T)[0,1])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('alignment',type=Path); p.add_argument('--out',type=Path,required=True)
    p.add_argument('--reference',action='append',default=[])
    p.add_argument('--identity',type=float,default=.8)
    args=p.parse_args()
    analyze_alignment(args.alignment,args.out,args.reference,args.identity)
