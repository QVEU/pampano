"""Copy-resolved structure painting and visualization-only mmCIF exports."""
import csv, hashlib, json, math
from pathlib import Path
import gemmi
import numpy as np
from scipy.spatial import cKDTree
from network_analysis import protein_node

def residue_key(chain, number, icode, name, copy):
    return (str(chain),None if number is None else int(number),str(icode or ''),str(name),int(copy))

def map_atoms(g, atoms):
    lookup={}
    for n,d in g.nodes(data=True):
        key=residue_key(d['chain'],d['res_seq'],d['icode'],d['res_name'],d['copy'])
        if key in lookup: raise ValueError('Ambiguous residue identity while mapping scores')
        lookup[key]=n
    groups={n:[] for n in g}
    for i in range(len(atoms)):
        if atoms.element[i] in ('H','D'):continue
        key=residue_key(atoms.chain[i],None if atoms.resSeqNaN[i] else atoms.resSeq[i],atoms.iCode[i],atoms.resName[i],atoms.copyId[i])
        if key not in lookup:raise ValueError(f'Atom does not map to a graph residue: {key}')
        groups[lookup[key]].append(i)
    if any(not v for v in groups.values()):raise ValueError('A scored graph residue has no coordinate atoms')
    return groups

def export_score_cif(records,path,title):
    """B_iso carries score percent in a NEW visualization-only file.

    Atom positions and residue identities are unchanged. Original B and node
    identity are retained in a custom loop; the input structure is never edited.
    """
    # Build children before appending parents: Gemmi's vector-backed references
    # can be invalidated when another chain/residue is appended. A ligand may
    # share an author chain with protein yet occur after other chains in input.
    grouped={}
    for a in records:
        key=(a['resi'],a['icode'],a['resn'],a['copy'])
        grouped.setdefault(a['chain'],{}).setdefault(key,[]).append(a)
    st=gemmi.Structure();st.name='pampano_scores';model=gemmi.Model('1')
    ordered=[]
    for ch,rs in grouped.items():
        chain=gemmi.Chain(ch)
        for key,items in rs.items():
            a=items[0];r=gemmi.Residue();r.name=a['resn'];r.seqid=gemmi.SeqId(a['resi'],a['icode'] or ' ')
            r.het_flag='H' if a['hetflag'] else 'A';r.subchain=ch
            for a in items:
                atom=gemmi.Atom();atom.name=a['atom'];atom.element=gemmi.Element(a['elem']);atom.pos=gemmi.Position(a['x'],a['y'],a['z'])
                atom.b_iso=100*a['score'];atom.occ=a['occupancy'];atom.serial=len(ordered)+1
                r.add_atom(atom);ordered.append(a)
            chain.add_residue(r)
        model.add_chain(chain)
    st.add_model(model);records=ordered
    st.setup_entities();doc=st.make_mmcif_document();block=doc.sole_block()
    block.set_pair('_pampano_visualization.description',gemmi.cif.quote(title+'; B_iso_or_equiv is 100 times network score, NOT an experimental B-factor'))
    loop=block.init_loop('_pampano_atom_mapping.',['ordinal','node_id','original_chain','copy','original_b_iso','score_fraction','direct_ligand_contact','selected_ligand'])
    for i,a in enumerate(records,1):
        loop.add_row([str(i),gemmi.cif.quote(a['node']),gemmi.cif.quote(a['original_chain']),str(a['copy']),str(a['original_b']),str(a['score']),str(int(a['contact'])),str(int(a['ligand']))])
    doc.write_file(str(path))
    check=gemmi.read_structure(str(path))
    back=[a for m in check for c in m for r in c for a in r]
    if len(back)!=len(records):raise ValueError('Coordinate export lost atoms')
    for a,b in zip(back,records):
        if np.linalg.norm(np.array([a.pos.x,a.pos.y,a.pos.z])-np.array([b['x'],b['y'],b['z']]))>0.002:
            raise ValueError('Coordinate export changed atom identity/order')

def make_atom(g,atoms,raw,groups,n,i,score,contact=False,ligand=False):
    d=g.nodes[n]
    if d['res_seq'] is None:raise ValueError('Structure painting requires author residue numbers')
    # Duplicate same-ID residues become distinct display chains, with a mapping.
    chain=d['chain'] if not d['copy'] else f"{d['chain']}__copy{d['copy']}"
    original=raw[int(atoms.src[i])]
    return dict(x=float(atoms.x[i]),y=float(atoms.y[i]),z=float(atoms.z[i]),elem=atoms.element[i],atom=atoms.name[i],
        chain=chain,original_chain=d['chain'],resi=int(d['res_seq']),icode=d['icode'],resn=d['res_name'],copy=d['copy'],node=n,
        hetflag=bool(d['hetatm']),protein=protein_node(d),score=float(score),contact=bool(contact),ligand=bool(ligand),
        original_b=float(original.get('b_iso',0)),occupancy=float(original.get('occupancy',1)))

def packed_atoms(records, bonds=False):
    adjacency=[[] for _ in records]
    if bonds and records:
        points=np.array([[a['x'],a['y'],a['z']] for a in records])
        radii={'C':.76,'N':.71,'O':.66,'S':1.05,'P':1.07,'F':.57,'CL':1.02,'BR':1.2,'I':1.39}
        for i,j in cKDTree(points).query_pairs(2.8):
            cutoff=radii.get(records[i]['elem'].upper(),.77)+radii.get(records[j]['elem'].upper(),.77)+.4
            distance=np.linalg.norm(points[i]-points[j])
            if .4<distance<=cutoff:adjacency[i].append(j);adjacency[j].append(i)
    return [[round(a['x'],4),round(a['y'],4),round(a['z'],4),a['elem'],a['atom'],a['chain'],a['resi'],a['icode'],a['resn'],a['node'],int(a['protein']),a['score'],int(a['contact']),int(a['ligand']),adjacency[i]] for i,a in enumerate(records)]

def export_views(g,atoms,raw,total,sites,out,title,viewer_js=None,radius=12,max_view_sites=120,global_exact=False):
    if radius<=0 or max_view_sites<1:raise ValueError('Viewer radius and site limit must be positive')
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    groups=map_atoms(g,atoms)
    scores={str(n):float(x) for n,x in total['betweenness_mean'].items()}
    global_atoms=[]
    for n,d in g.nodes(data=True):
        if not protein_node(d):continue
        indices=[i for i in groups[n] if atoms.name[i]=='CA']
        if not indices:continue
        global_atoms.append(make_atom(g,atoms,raw,groups,n,indices[0],scores[n]))
    xyz=np.column_stack((atoms.x,atoms.y,atoms.z));tree=cKDTree(xyz)
    node_for_atom={i:n for n,indices in groups.items() for i in indices}
    payload_sites=[];first_atoms=None;first_site=None
    for site in sites[:max_view_sites]:
        ligand_ids=site['ligand_nodes'];li=[i for n in ligand_ids for i in groups[n]]
        context=li or [i for n in site['terminals'] for i in groups[n] if atoms.name[i]=='CA']
        near=set(i for around in tree.query_ball_point(xyz[context],radius) for i in around)
        residue_ids={node_for_atom[i] for i in near if i in node_for_atom and protein_node(g.nodes[node_for_atom[i]])}
        residue_ids.update(site['nodes']);residue_ids.update(ligand_ids)
        # Add sequential backbone context only when physically adjacent.
        expanded=set(residue_ids)
        for n in residue_ids:
            for neighbor in g[n]:
                if g[n][neighbor]['interaction']=='peptide':expanded.add(neighbor)
        records=[];terminals=set(site['terminals'])
        for n in sorted(expanded,key=lambda x:int(x)):
            for i in groups[n]:records.append(make_atom(g,atoms,raw,groups,n,i,site['scores'].get(n,0),n in terminals,n in ligand_ids))
        labels=[dict(node=n,name=g.nodes[n]['name'],score=v,contact=n in terminals) for n,v in sorted(site['scores'].items(),key=lambda x:x[1],reverse=True)]
        site_name=site['ligand']
        if 0<len(ligand_ids)<=3:site_name+=' | '+', '.join(g.nodes[n]['name'] for n in ligand_ids)
        elif ligand_ids:site_name+=f' | {len(ligand_ids)} residues'
        payload_sites.append(dict(site=site['site'],name=site_name,kind=site.get('analysis_kind','ligand'),
            contacts=site['contacts'],pairs=site['pairs'],status=site['status'],atoms=packed_atoms(records,bonds=True),ranking=labels))
        if first_atoms is None and site['status']=='ok':first_atoms=records;first_site=site
    payload=dict(title=title,global_atoms=packed_atoms(global_atoms),sites=payload_sites,site_count=len(sites),global_exact=global_exact,
        displayed_sites=len(payload_sites),network_mode=g.graph['assembly_mode'],network_cutoff=g.graph['cutoff'],
        global_ranking=[dict(node=n,name=g.nodes[n]['name'],score=scores[n],contact=False) for n in sorted(scores,key=scores.get,reverse=True)[:30] if protein_node(g.nodes[n])])
    template=(Path(__file__).parent/'structure_view.html').read_text()
    script=f'<script>{Path(viewer_js).read_text()}</script>' if viewer_js else '<script src="https://unpkg.com/3dmol@2.5.5/build/3Dmol-min.js"></script>'
    data=json.dumps(payload,separators=(',',':'),allow_nan=False).replace('<','\\u003c')
    (out/'structure_view.html').write_text(template.replace('<!-- RENDERER -->',script).replace('/* PAYLOAD */',data))
    export_score_cif(global_atoms,out/'total_CA_score_percent.cif','Global betweenness: alpha-carbon representation')
    if first_atoms is not None:
        export_score_cif(first_atoms,out/'ligand_site_000_score_percent.cif','Ligand contact-restricted betweenness: first valid site')
        (out/'ligand_site_000_scores.json').write_text(json.dumps(first_site,indent=2))
        (out/'paint_ligand.pml').write_text('''# Visualization only: the supplied CIF B column contains score percent.
load ligand_site_000_score_percent.cif, ligand_network
hide everything
show cartoon, polymer.protein
spectrum b, white_yellow_red, polymer.protein, minimum=0, maximum=25
show sticks, b > 0 or not polymer.protein
color forest, not polymer.protein
bg_color white
orient not polymer.protein
zoom all
''')
    (out/'paint_total.pml').write_text('''# Visualization only: CA-only model; B column contains score percent.
load total_CA_score_percent.cif, total_network
hide everything
show spheres
set sphere_scale, 0.7
spectrum b, white_yellow_red, all, minimum=0, maximum=2
bg_color white
orient
''')
    metadata=dict(atom_mapping='physical chain, author residue number, insertion code, residue name, duplicate-copy index',
        source_coordinates_unchanged=True,visualization_b_factor='100 times network score; original B values retained in custom mmCIF atom-mapping loop',
        global_representation='one alpha carbon per protein residue; atoms lacking CA are omitted only from this view',
        global_CA_nodes=len(global_atoms),protein_nodes=sum(protein_node(d) for _,d in g.nodes(data=True)),
        ligand_context_radius=radius,view_sites=len(payload_sites),total_sites=len(sites),
        renderer_version='3Dmol.js 2.5.5',embedded_renderer=bool(viewer_js),renderer_sha256=hashlib.sha256(Path(viewer_js).read_bytes()).hexdigest() if viewer_js else None)
    (out/'mapping_validation.json').write_text(json.dumps(metadata,indent=2)+'\n')
