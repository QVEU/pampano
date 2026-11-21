"""Exact path participation across a specified protein-protein interface."""
import gzip, json, math
from collections import defaultdict
from pathlib import Path
import networkx as nx
import pandas as pd
from network_analysis import protein_node, describe

def exact_cross_scores(g,sources,targets):
    sources=list(dict.fromkeys(sources));targets=list(dict.fromkeys(targets))
    if set(sources)&set(targets):raise ValueError('Cross-interface endpoint sets must be disjoint')
    if any(n not in g for n in sources+targets):raise ValueError('Endpoint absent from path graph')
    pairs=len(sources)*len(targets);keep=set(sources+targets);reachable=0;internal=0;maximum=0
    ds={n:nx.single_source_shortest_path_length(g,n) for n in sources+targets}
    for s in sources:
        for t in targets:
            if t not in ds[s]:continue
            length=ds[s][t];reachable+=1;internal+=length-1;maximum=max(maximum,length)
            keep.update(v for v in ds[s].keys()&ds[t].keys() if ds[s][v]+ds[t][v]==length)
    local=g.subgraph(keep).copy()
    # NetworkX assigns half a path for disjoint source/target sets on an
    # undirected graph. Restore one unit per unordered cross-interface pair.
    raw={n:2*v for n,v in nx.betweenness_centrality_subset(local,sources,targets,normalized=False,weight=None).items()}
    if not math.isclose(sum(raw.values()),internal,abs_tol=1e-7):raise ValueError('Cross-interface path credit failed validation')
    return local,{n:v/pairs if pairs else 0 for n,v in raw.items()},dict(pairs=pairs,reachable_pairs=reachable,unreachable_pairs=pairs-reachable,max_path_length=maximum)

def interface_analysis(g,out,left_chains,right_chains,contact_cutoff=None):
    if not left_chains or not right_chains or set(left_chains)&set(right_chains):raise ValueError('Specify two disjoint sets of author protein chains')
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    p=nx.Graph(g.subgraph([n for n,d in g.nodes(data=True) if protein_node(d)]))
    groups=defaultdict(list);contacts=[]
    for u,v,d in p.edges(data=True):
        if contact_cutoff is not None and d['distance']>contact_cutoff:continue
        if p.nodes[v]['base_chain'] in left_chains and p.nodes[u]['base_chain'] in right_chains:u,v=v,u
        a,b=p.nodes[u],p.nodes[v]
        if a['base_chain'] not in left_chains or b['base_chain'] not in right_chains:continue
        groups[(a['chain'],b['chain'])].append((u,v,d))
    if not groups:raise ValueError('No contacts connect the requested protein-chain groups')
    paths=[];rows=[];site_rows=[]
    for site,((left,right),edges) in enumerate(sorted(groups.items())):
        sources=sorted({u for u,v,d in edges},key=int);targets=sorted({v for u,v,d in edges},key=int)
        local,scores,stats=exact_cross_scores(p,sources,targets);terminals=set(sources+targets)
        for u,v,d in edges:
            a,b=p.nodes[u],p.nodes[v]
            contacts.append(dict(site=site,left_node=u,right_node=v,left_chain=a['chain'],left_res_name=a['res_name'],left_res_seq=a['res_seq'],left_icode=a['icode'],right_chain=b['chain'],right_res_name=b['res_name'],right_res_seq=b['res_seq'],right_icode=b['icode'],distance=d['distance'],interaction=d['interaction']))
        state=dict(site=site,ligand=f'{left} / {right}',ligand_nodes=[],analysis_kind='interface',left_chain=left,right_chain=right,
            left_contacts=len(sources),right_contacts=len(targets),contact_edges=len(edges),contacts=len(terminals),
            status='ok',geodesic_nodes=len(local),geodesic_edges=local.number_of_edges(),**stats)
        site_rows.append(state)
        for n,score in scores.items():
            d=p.nodes[n];rows.append(dict(site=site,node_id=n,chain=d['chain'],base_chain=d['base_chain'],res_name=d['res_name'],res_seq=d['res_seq'],icode=d['icode'],direct_interface_contact=n in terminals,interface_degree=sum(1 for u,v,e in edges if n in (u,v)),pair_fraction=score))
        paths.append(dict(**state,sources=sources,targets=targets,terminals=sorted(terminals,key=int),scores=scores,nodes=list(local),edges=list(local.edges())))
    frame=pd.DataFrame(rows).sort_values(['site','pair_fraction'],ascending=[True,False])
    frame.to_csv(out/'site_residue_scores.csv',index=False)
    frame[~frame.direct_interface_contact].to_csv(out/'noncontact_intermediates.csv',index=False)
    pd.DataFrame(contacts).to_csv(out/'direct_contacts.csv',index=False);pd.DataFrame(site_rows).to_csv(out/'sites.csv',index=False)
    with gzip.open(out/'site_paths.json.gz','wt') as f:json.dump(paths,f)
    summary=dict(site_count=len(paths),path_graph=describe(p),left_chains=left_chains,right_chains=right_chains,
        normalization='number of cross-interface contact-residue pairs; unreachable pairs contribute zero',
        endpoint_exclusion=True,path_weights=None,contact_cutoff=contact_cutoff or g.graph['cutoff'])
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    return paths,summary
