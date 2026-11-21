# Comparing contact networks with sequence-based scores

`sequence_analysis.py` takes an existing protein multiple-sequence alignment
and computes conservation and continuous statistical coupling analysis (SCA)
scores without using structures or network rankings in model fitting.

```bash
python -m pip install -r offline/requirements-sequence.txt
python offline/sequence_analysis.py homologs_with_reference.fasta \
  --reference REF_target --out results/sequence
```

Use `--reference` once per mapping-only reference sequence. Those sequences
are excluded from fitting. Keep references in the alignment, including gaps;
do not supply an unaligned FASTA. Output `reference_mapping.csv` relates each
reference's ungapped, one-based position to the original one-based alignment
column. It explicitly marks columns not scored after filtering. Structure
author residue numbers and insertion codes require a separate validated map;
they must not be assumed to equal alignment indices.

## Method and outputs

- Select the union of reference amino-acid columns, if references are supplied.
- Remove sequences with more than 20% gaps/unknown states, then collapse exact
  duplicates. Unrecognized amino-acid symbols are treated as missing states.
- Weight sequences by the inverse number of neighbors at strictly greater than
  80% identity, including matching gaps. `--identity .9` gives a useful
  sensitivity analysis. Effective sequence count is the sum of these weights.
- Remove columns with more than 20% weighted gaps in one pass, then recompute
  weights. Because the weights change, the final gap fraction can differ
  slightly from the filtering threshold; both policy and final maximum are
  recorded. This weighting reduces redundancy, not all phylogenetic effects.
- Compute normalized Shannon conservation (gaps excluded from amino-acid
  frequencies), relative entropy using the pySCA background, and the
  unregularized Frobenius SCA matrix in memory-bounded blocks.
- Save the leading eigenvector's absolute loading and the norm of loadings
  across the top six eigenvectors. Six is fixed before any structural
  comparison. These are **exploratory continuous scores, not inferred sectors**.
  Sector inference requires additional null, component, phylogeny and
  stability analyses. DCA/Potts couplings are not computed by this module.

`sequence_scores.csv`, `sequence_weights.csv`, `reference_mapping.csv`,
`sca_matrix.npz` and `summary.json` retain the scores, matrix, eigenvectors,
input hash, thresholds and depth diagnostics. A low-depth flag (`N_eff < L`)
is a heuristic warning, not a universal SCA validity threshold. No sequence
alignment, structure repair, downloads or experimental validation are implied
by running this module.

The SCA equations follow Rivoire, Reynolds & Ranganathan (2016),
[Evolution-Based Functional Decomposition of Proteins](https://doi.org/10.1371/journal.pcbi.1004817).
The implementation was compared with the authors' pySCA at commit
`d848ee4d1d7c5270a87e0fdc7362b7016b3dbe02` on a 667-sequence, 157-position DHFR
alignment: maximum matrix difference 2.67e-15 and weight difference 6.94e-17.
This verifies numerical agreement, not biological predictive validity.

## Comparing with Pampano

Join by original alignment column and a structure-specific chain/residue map.
For full capsids, average equivalent symmetry copies before residue-level
statistics. Do not treat sixty copies as sixty independent sequences or
replicates. Missing structural residues and filtered sequence positions stay
missing. A modeled residue absent from an exhaustive shortest-path ranking
can receive a zero path score; an unmodeled residue cannot.

Compare the global and partner-contact scores separately, and examine direct
contacts and noncontact residues separately. Report conservation and degree
controls alongside SCA correlations. Packing, distance to a functional site,
sequence coverage and domain can also confound overlap. Rank correlations
among residues are descriptive unless an appropriate dependent-data null is
specified. Agreement with a sequence model alone does not establish allostery.

For capsid interchain coevolution, keep VP4–VP2–VP3–VP1 sequences paired from
the same viral polyprotein. For NS5 in multiple complexes, use one shared
alignment and map each chain independently. Comparing NS5 with RNA or STAT2
does not justify arbitrarily concatenating their sequences into a paired
coevolution alignment.

## DHFR benchmark interpretation

Reynolds, McLaughlin & Ranganathan (2011),
[Hot Spots for Allosteric Regulation on Protein Surfaces](https://doi.org/10.1016/j.cell.2011.10.049),
provides fixed sector labels for E. coli DHFR (1RX2). The lab's later tutorial
alignment is different from the 418-sequence alignment underlying those labels;
comparison is not an exact rerun of the 2011 analysis.

McCormick et al. (2021),
[Structurally distributed surface sites tune allosteric regulation](https://doi.org/10.7554/eLife.68346),
provides an independent experimental readout in a DHFR–LOV2 fusion. It tests a
different object from unmodified 1RX2: preserve published missing-data masks,
separate loss of activity from changes in allostery, and flag sequence
differences (1RX2 residue 37 is D, while the assayed DHFR background has N).
Published sectors, sequence conservation, global centrality and site-specific
paths need not predict the same biological property.
