import math
import sys
from pathlib import Path
import networkx as nx
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from network_analysis import exact_site_scores
from structure_view import residue_key, export_score_cif
from interface_analysis import exact_cross_scores

def test_all_tied_paths_match_full_networkx():
    g=nx.Graph([('a','x'),('a','y'),('x','b'),('y','b'),('x','z'),('z','c'),('b','c')])
    terminals=['a','b','c']
    local,scores,stats=exact_site_scores(g,terminals)
    expected=nx.betweenness_centrality_subset(g,terminals,terminals,normalized=False)
    assert stats['pairs']==3
    for n in g: assert scores.get(n,0)==pytest.approx(expected[n]/3)

def test_disconnected_pairs_are_reported_and_zero():
    g=nx.path_graph(['a','x','b']);g.add_node('c')
    local,scores,stats=exact_site_scores(g,['a','b','c'])
    assert stats['reachable_pairs']==1 and stats['unreachable_pairs']==2
    assert scores['x']==pytest.approx(1/3)

def test_single_terminal_is_not_claimed_as_an_analysis():
    g=nx.path_graph(3)
    local,scores,stats=exact_site_scores(g,[0])
    assert scores=={} and stats['pairs']==0

def test_cross_interface_scores_match_enumerated_unordered_paths():
    g=nx.Graph([('s1','x'),('s1','y'),('s2','y'),('x','t1'),('y','t1'),('y','z'),('z','t2')])
    sources=['s1','s2'];targets=['t1','t2']
    local,scores,stats=exact_cross_scores(g,sources,targets)
    expected={n:0. for n in g}
    for s in sources:
        for t in targets:
            paths=list(nx.all_shortest_paths(g,s,t))
            for path in paths:
                for n in path[1:-1]:expected[n]+=1/(len(paths)*4)
    assert stats['pairs']==4
    for n in g:assert scores.get(n,0)==pytest.approx(expected[n])
    with pytest.raises(ValueError):exact_cross_scores(g,['s1'],['s1'])

def test_rna_atom_records_are_explicitly_selectable_as_one_partner(tmp_path):
    from network_analysis import ligand_analysis
    g=nx.Graph(cutoff=5.)
    for n,chain,res,num,het in [('0','A','ALA',1,False),('1','A','GLY',2,False),('2','A','VAL',3,False),('3','B','A',1,False),('4','B','U',2,False),('5','A','ZN',99,True)]:
        g.add_node(n,chain=chain,base_chain=chain,res_name=res,res_seq=num,icode='',hetatm=het,operator=0,copy=0)
    g.add_edges_from([('0','1'),('1','2'),('0','3'),('2','4'),('1','5')],distance=4.)
    frame,rank,sites,summary=ligand_analysis(g,tmp_path,partner_chains=['B'],group='chain',contact_cutoff=5.)
    assert len(sites)==1 and sites[0]['ligand_nodes']==['3','4']
    assert sites[0]['scores']['1']==1. and summary['path_graph']['nodes']==3
    assert set(sites[0]['terminals'])=={'0','2'}

def test_residue_mapping_preserves_insertion_copy_and_chain():
    keys={residue_key('A',10,'','ALA',0),residue_key('A',10,'A','ALA',0),residue_key('A_4',10,'','ALA',0),residue_key('A',10,'','ALA',1)}
    assert len(keys)==4

def test_score_cif_retains_coordinates_insertion_codes_and_source_b(tmp_path):
    import gemmi
    records=[]
    for i,(chain,number,icode) in enumerate([('AA',10,''),('AA',10,'A'),('AA_4',10,''),('AA',301,'')]):
        records.append(dict(chain=chain,original_chain=chain,resi=number,icode=icode,resn='ALA',copy=0,
            atom='CA',elem='C',x=4.*i,y=1.,z=2.,hetflag=False,score=.123,occupancy=1.,original_b=32.,node=str(i),contact=False,ligand=False))
    path=tmp_path/'scores.cif';export_score_cif(records,path,'test')
    st=gemmi.read_structure(str(path));atoms=[a for c in st[0] for r in c for a in r]
    assert len(atoms)==4 and all(a.b_iso==pytest.approx(12.3) for a in atoms)
    assert st[0]['AA'][1].seqid.icode=='A'
    doc=gemmi.cif.read_file(str(path)).sole_block()
    assert list(doc.find_values('_pampano_atom_mapping.original_b_iso'))==['32.0']*4
    assert list(doc.find_values('_pampano_atom_mapping.node_id'))==['0','1','3','2']

def test_batch_reuse_and_input_preservation(tmp_path,monkeypatch):
    import json
    from analyze_structure import main, run_job
    records=[]
    for i,x in enumerate([0.,3.5,7.]):
        records.append(dict(chain='A',original_chain='A',resi=i+1,icode='',resn='ALA',copy=0,atom='CA',elem='C',x=x,y=0.,z=0.,hetflag=False,score=0.,occupancy=1.,original_b=20.,node=str(i),contact=False,ligand=False))
    source=tmp_path/'tiny.cif';export_score_cif(records,source,'fixture')
    batch=tmp_path/'jobs.json'
    jobs=[dict(id=n,input='tiny.cif',assembly='asym',no_views=True,betweenness_sources=10) for n in ['first','second']]
    batch.write_text(json.dumps(jobs));out=tmp_path/'results'
    argv=['analyze_structure.py','--batch',str(batch),'--out',str(out)]
    monkeypatch.setattr(sys,'argv',argv);main()
    assert all(r['status']=='ok' for r in json.loads((out/'batch_results.json').read_text()))
    monkeypatch.setattr(sys,'argv',argv+['--reuse']);main()
    assert all(r['status']=='ok' for r in json.loads((out/'batch_results.json').read_text()))
    saved=out/'first/inputs/tiny.cif';original=saved.read_bytes()
    records[0]['x']=1.;export_score_cif(records,source,'modified fixture')
    with pytest.raises(ValueError,match='Output already contains'):
        run_job(dict(jobs[0],input=str(source)),out/'first',reuse=True)
    assert saved.read_bytes()==original
