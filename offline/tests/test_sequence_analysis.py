import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from sequence_analysis import (BACKGROUND, encode, sequence_weights, sca_matrix,
                               single_site_scores, analyze_alignment)


def test_weights_match_pairwise_enumeration():
    a=encode(['ACD--','ACD--','ACDEF','ACDEY','WYDEF'])
    expected=1/((a[:,None,:]==a[None,:,:]).mean(axis=2)>.8).sum(axis=1)
    np.testing.assert_allclose(sequence_weights(a,block=2),expected)


def test_block_sca_matches_direct_frequency_tensor():
    a=encode(['ACDE','ACDY','FCEY','AC-E','FCDY'])
    w=np.array([.5,.5,1,1,1]); wn=w/w.sum()
    got,f=sca_matrix(a,w,block=2)
    expected=np.zeros_like(got)
    for i in range(a.shape[1]):
        for j in range(a.shape[1]):
            joint=np.zeros((20,20))
            for s in range(len(a)):
                if a[s,i]<20 and a[s,j]<20:joint[a[s,i],a[s,j]]+=wn[s]
            z=[]
            for pos in [i,j]:
                valid=(f[pos]>0)&(f[pos]<1)
                u=np.zeros(20)
                u[valid]=abs(np.log(f[pos,valid]*(1-BACKGROUND[valid])/
                       ((1-f[pos,valid])*BACKGROUND[valid])))
                z.append(u)
            expected[i,j]=np.linalg.norm((joint-np.outer(f[i],f[j]))*np.outer(*z))
    np.testing.assert_allclose(got,expected,atol=1e-12)


def test_reference_excluded_and_gap_column_unscored(tmp_path):
    p=tmp_path/'a.fasta'
    p.write_text('>a\nACD-\n>b\nACE-\n>c\nAFD-\n>REF\nACDE\n')
    scores,meta=analyze_alignment(p,tmp_path/'out',['REF'],max_sequence_gap=.3)
    assert meta['fit_sequences']==3
    assert scores.msa_column.tolist()==[1,2,3]
    mapping=pd.read_csv(tmp_path/'out/reference_mapping.csv')
    assert mapping.scored.tolist()==[True,True,True,False]
    with pytest.raises(ValueError,match='not found'):
        analyze_alignment(p,tmp_path/'bad',['missing'])


def test_entropy_extremes():
    f=np.zeros((2,20));f[0,0]=1;f[1,:]=1/20
    s=single_site_scores(f)
    np.testing.assert_allclose(s['conservation'],[1,0],atol=1e-12)
