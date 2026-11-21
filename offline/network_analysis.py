"""Graph analysis shared by the single-structure and batch CLI."""
import csv, gzip, hashlib, json, math, random
from collections import Counter
from itertools import combinations
from pathlib import Path

import gemmi
import networkx as nx
import numpy as np
import pandas as pd

WATERS = {'HOH', 'WAT', 'DOD'}

def load_graph(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as f: data = json.load(f)
    g = nx.Graph(**data['data'])
    for item in data['elements']['nodes']:
        d = item['data']; n = str(d['id'])
        if n in g: raise ValueError('Duplicate node ID')
        g.add_node(n, **d)
    for item in data['elements']['edges']:
        d = item['data']; u, v = str(d['source']), str(d['target'])
        if u not in g or v not in g or u == v or g.has_edge(u, v):
            raise ValueError('Invalid or duplicate graph edge')
        g.add_edge(u, v, **d)
    if any(g.degree(n) != d['degree'] for n,d in g.nodes(data=True)):
        raise ValueError('Exported degrees disagree with edges')
    return g

def protein_node(d):
    return gemmi.find_tabulated_residue(d['res_name']).is_amino_acid()

def describe(g):
    sizes = sorted(map(len, nx.connected_components(g)), reverse=True)
    return dict(nodes=len(g), edges=g.number_of_edges(), components=len(sizes),
                largest_component=sizes[0] if sizes else 0, isolates=nx.number_of_isolates(g),
                mean_degree=2*g.number_of_edges()/len(g) if len(g) else 0, density=nx.density(g))

def topology(g):
    h = nx.Graph(); h.add_nodes_from(g); h.add_edges_from(g.edges()); return h

def total_analysis(g, out, k=256, seeds=(17,29)):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    h = topology(g)
    nodes = pd.DataFrame([dict(d, id=n) for n,d in g.nodes(data=True)]).set_index('id', drop=False)
    nodes['degree'] = pd.Series(dict(h.degree()))
    nodes['clustering'] = pd.Series(nx.clustering(h, weight=None))
    nodes['core_number'] = pd.Series(nx.core_number(h))
    for seed in seeds:
        print(f'Total betweenness: {min(k,len(h))} sources, seed {seed}', flush=True)
        bc = nx.betweenness_centrality(h,k=min(k,len(h)),normalized=True,weight=None,seed=seed)
        nodes[f'betweenness_seed{seed}'] = pd.Series(bc)
        with gzip.open(out/f'betweenness_seed{seed}.json.gz','wt') as f: json.dump(bc,f)
    nodes['betweenness_mean'] = nodes[[f'betweenness_seed{s}' for s in seeds]].mean(axis=1)
    community_runs = []
    for seed in seeds:
        print(f'Louvain communities: seed {seed}', flush=True)
        if h.number_of_edges():
            communities = nx.community.louvain_communities(h,weight=None,resolution=1,seed=seed)
            q = nx.community.modularity(h,communities,weight=None)
        else: communities=[{n} for n in h]; q=None
        communities = sorted(communities,key=lambda c:(-len(c),min(c)))
        community_runs.append(dict(seed=seed,communities=len(communities),modularity=q))
        if seed == seeds[0]: nodes['community'] = pd.Series({n:i for i,c in enumerate(communities) for n in c})
    interactions = Counter(); interfaces = Counter(); pairs = Counter(); scopes = Counter(); idegree=Counter()
    for u,v,d in g.edges(data=True):
        a,b=g.nodes[u],g.nodes[v]
        cls = 'between_symmetry_copies' if a['operator']!=b['operator'] else ('between_chains_same_copy' if a['base_chain']!=b['base_chain'] else 'within_chain_copy')
        interactions[d['interaction']]+=1; interfaces[cls]+=1; scopes[d['scope']]+=1
        if cls!='within_chain_copy':
            idegree[u]+=1; idegree[v]+=1
            pairs[(*sorted((a['base_chain'],b['base_chain'])),cls)]+=1
    nodes['interface_degree'] = pd.Series({n:idegree[n] for n in g})
    keys=['base_chain','res_name','res_seq','icode','hetatm']
    residues=nodes.groupby(keys,dropna=False).agg(copies=('id','size'),degree_mean=('degree','mean'),
        degree_min=('degree','min'),degree_max=('degree','max'),clustering_mean=('clustering','mean'),
        core_number_mean=('core_number','mean'),betweenness_mean=('betweenness_mean','mean'),
        interface_degree_mean=('interface_degree','mean'),**{f'betweenness_seed{s}':(f'betweenness_seed{s}','mean') for s in seeds}).reset_index()
    residues=residues.sort_values('betweenness_mean',ascending=False)
    pp=[n for n,d in g.nodes(data=True) if protein_node(d)]
    graphs={'full_with_ligands':describe(g),'protein_only':describe(g.subgraph(pp))}
    no_peptide=nx.Graph(); no_peptide.add_nodes_from(g)
    no_peptide.add_edges_from((u,v) for u,v,d in g.edges(data=True) if d['interaction']!='peptide')
    graphs['peptide_label_edges_removed']=describe(no_peptide)
    largest=h.subgraph(max(nx.connected_components(h),key=len)).copy() if len(h) else h
    sources=random.Random(101).sample(list(largest),min(64,len(largest)))
    path_means=[]; diameter_bound=0
    for source in sources:
        ds=nx.single_source_shortest_path_length(largest,source)
        path_means.append(sum(ds.values())/max(1,len(largest)-1)); diameter_bound=max(diameter_bound,max(ds.values()))
    correlation=residues[f'betweenness_seed{seeds[0]}'].corr(residues[f'betweenness_seed{seeds[1]}'],method='spearman') if len(seeds)>1 else None
    summary=dict(graphs=graphs,mean_clustering=float(nodes.clustering.mean()),max_core_number=int(nodes.core_number.max()),
        community_runs=community_runs,betweenness_sources_per_seed=min(k,len(h)),betweenness_seeds=list(seeds),
        betweenness_exact=k>=len(h),betweenness_residue_rank_correlation=None if correlation is None or not np.isfinite(correlation) else float(correlation),
        sampled_mean_path_length_largest_component=float(np.mean(path_means)) if path_means else None,
        path_sample_sources=len(sources),diameter_lower_bound=diameter_bound,
        interactions=dict(interactions),interfaces=dict(interfaces),networkx_version=nx.__version__)
    nodes.to_csv(out/'node_metrics.csv.gz',index=False); residues.to_csv(out/'residue_metrics.csv',index=False)
    pd.DataFrame(community_runs).to_csv(out/'community_runs.csv',index=False)
    pd.DataFrame([dict(chain_1=a,chain_2=b,interface_class=c,edges=n) for (a,b,c),n in pairs.items()]).to_csv(out/'interface_chain_pairs.csv',index=False)
    for name,counts,col in [('interaction_counts',interactions,'interaction'),('interface_counts',interfaces,'interface_class'),('scope_counts',scopes,'pampano_scope')]:
        pd.DataFrame(counts.items(),columns=[col,'edges']).to_csv(out/f'{name}.csv',index=False)
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    return nodes, residues, summary

def exact_site_scores(g, terminals):
    """Exact undirected subset betweenness; disconnected pairs contribute zero.

    Keep all vertices on any shortest terminal-terminal path. An induced graph
    on these vertices preserves every tied shortest path for the endpoint set.
    """
    terminals=list(dict.fromkeys(terminals))
    if any(n not in g for n in terminals): raise ValueError('Terminal absent from path graph')
    pair_count=math.comb(len(terminals),2)
    if not pair_count: return g.subgraph(terminals).copy(),{},dict(pairs=0,reachable_pairs=0,unreachable_pairs=0)
    component={n:i for i,c in enumerate(nx.connected_components(g)) for n in c}
    by_component={}
    for n in terminals: by_component.setdefault(component[n],[]).append(n)
    keep=set(terminals); reachable=0; internal_total=0; max_length=0
    for ts in by_component.values():
        distances={}
        for source in ts:
            cutoff=4
            while True:
                ds=nx.single_source_shortest_path_length(g,source,cutoff=cutoff)
                if all(t in ds for t in ts): break
                cutoff*=2
            distances[source]=ds
        for a,b in combinations(ts,2):
            da,db=distances[a],distances[b]; length=da[b]
            keep.update(v for v in da.keys()&db.keys() if da[v]+db[v]==length)
            reachable+=1;internal_total+=length-1;max_length=max(max_length,length)
    local=g.subgraph(keep).copy()
    raw=nx.betweenness_centrality_subset(local,terminals,terminals,normalized=False,weight=None)
    if not math.isclose(sum(raw.values()),internal_total,abs_tol=1e-7): raise ValueError('Path-credit validation failed')
    return local,{n:v/pair_count for n,v in raw.items()},dict(pairs=pair_count,reachable_pairs=reachable,
        unreachable_pairs=pair_count-reachable,max_path_length=max_length)

def ligand_analysis(g, out, ligands=None, ligand_chain=None, ligand_residue=None, group='residue', contact_cutoff=None, partner_chains=None):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    p=nx.Graph();p.add_nodes_from((n,d) for n,d in g.nodes(data=True) if protein_node(d));p.add_edges_from((u,v) for u,v in g.edges() if u in p and v in p)
    selected=[]
    for n,d in g.nodes(data=True):
        if protein_node(d) or d['res_name'] in WATERS: continue
        if partner_chains:
            if d['base_chain'] not in partner_chains:continue
        elif not d['hetatm']:continue
        if ligands and d['res_name'] not in ligands: continue
        if ligand_chain is not None and d['base_chain']!=ligand_chain: continue
        if ligand_residue is not None and d['res_seq']!=ligand_residue: continue
        selected.append(n)
    if (ligands or partner_chains or ligand_chain is not None or ligand_residue is not None) and not selected: raise ValueError('No ligand/partner matched the requested selector')
    groups={}
    for n in sorted(selected,key=lambda n:(g.nodes[n]['operator'],g.nodes[n]['chain'],g.nodes[n]['res_seq'] or 0,g.nodes[n]['res_name'])):
        d=g.nodes[n];key=(d['chain'],d['operator']) if group=='chain' else (n,)
        groups.setdefault(key,[]).append(n)
    rows=[];contact_rows=[];sites=[];paths=[]
    for site,lids in enumerate(groups.values()):
        name=('partner chain '+g.nodes[lids[0]]['base_chain']) if partner_chains else '+'.join(sorted({g.nodes[n]['res_name'] for n in lids}))
        terminals=set()
        for ligand in lids:
            for n in g[ligand]:
                if n in p and (contact_cutoff is None or g[ligand][n]['distance']<=contact_cutoff):
                    terminals.add(n);d=g.nodes[n]
                    contact_rows.append(dict(site=site,ligand=name,ligand_node=ligand,partner_chain=g.nodes[ligand]['chain'],partner_res_name=g.nodes[ligand]['res_name'],partner_res_seq=g.nodes[ligand]['res_seq'],partner_icode=g.nodes[ligand]['icode'],node_id=n,chain=d['chain'],base_chain=d['base_chain'],res_name=d['res_name'],res_seq=d['res_seq'],icode=d['icode'],operator=d['operator'],distance=g[ligand][n]['distance']))
        terminals=sorted(terminals,key=lambda n:int(n))
        local,scores,stats=exact_site_scores(p,terminals)
        state=dict(site=site,ligand=name,ligand_nodes=lids,contacts=len(terminals),status='ok' if stats['pairs'] else 'insufficient_contacts',
            geodesic_nodes=len(local),geodesic_edges=local.number_of_edges(),**stats)
        sites.append(state)
        for n,score in scores.items():
            d=g.nodes[n];rows.append(dict(site=site,ligand=name,node_id=n,chain=d['chain'],base_chain=d['base_chain'],res_name=d['res_name'],res_seq=d['res_seq'],icode=d['icode'],operator=d['operator'],copy=d['copy'],direct_ligand_contact=n in terminals,pair_fraction=score))
        paths.append(dict(**state,terminals=terminals,scores=scores,nodes=list(local),edges=list(local.edges())))
    columns=['site','ligand','node_id','chain','base_chain','res_name','res_seq','icode','operator','copy','direct_ligand_contact','pair_fraction']
    frame=pd.DataFrame(rows,columns=columns)
    frame.to_csv(out/'site_residue_scores.csv',index=False)
    pd.DataFrame(contact_rows,columns=['site','ligand','ligand_node','partner_chain','partner_res_name','partner_res_seq','partner_icode','node_id','chain','base_chain','res_name','res_seq','icode','operator','distance']).to_csv(out/'direct_contacts.csv',index=False)
    pd.DataFrame(sites,columns=['site','ligand','ligand_nodes','contacts','status','geodesic_nodes','geodesic_edges','pairs','reachable_pairs','unreachable_pairs','max_path_length']).to_csv(out/'sites.csv',index=False)
    summaries=[]
    for lig,part in frame.groupby('ligand'):
        count=sum(s['status']=='ok' and s['ligand']==lig for s in sites)
        for key,sub in part.groupby(['base_chain','res_name','res_seq','icode'],dropna=False):
            summaries.append(dict(ligand=lig,base_chain=key[0],res_name=key[1],res_seq=key[2],icode=key[3],
                valid_sites=count,sites_with_path_residue=sub.site.nunique(),mean_pair_fraction=sub.pair_fraction.sum()/count,
                max_physical_residue_score=sub.pair_fraction.max(),direct_ligand_contact=bool(sub.direct_ligand_contact.any()),
                max_copies_per_site=int(sub.groupby('site').size().max())))
    ranking=pd.DataFrame(summaries,columns=['ligand','base_chain','res_name','res_seq','icode','valid_sites','sites_with_path_residue','mean_pair_fraction','max_physical_residue_score','direct_ligand_contact','max_copies_per_site'])
    ranking=ranking.sort_values('mean_pair_fraction',ascending=False)
    ranking.to_csv(out/'residue_ranking.csv',index=False)
    ranking[~ranking.direct_ligand_contact.astype(bool)].to_csv(out/'noncontact_intermediates.csv',index=False)
    with gzip.open(out/'site_paths.json.gz','wt') as f:json.dump(paths,f)
    summary=dict(site_count=len(sites),valid_sites=sum(s['status']=='ok' for s in sites),ligands=sorted({s['ligand'] for s in sites}),
        path_graph=describe(p),path_weights=None,endpoint_exclusion=True,normalization='all unordered contact pairs per site; unreachable pairs contribute zero',
        aggregation='sum physical-copy contributions per base residue identity within each site, then average over all valid sites of the same ligand label',
        contact_cutoff=contact_cutoff or g.graph['cutoff'],group=group,partner_chains=partner_chains)
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    return frame,ranking,paths,summary
