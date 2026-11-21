#!/usr/bin/env python3
"""Pampano -> total/ligand NetworkX analysis -> colored molecular structures."""
import argparse, csv, gzip, hashlib, json, re, shutil, sys, time, tempfile
from pathlib import Path
import numpy as np
import pandas as pd
import pampano_offline as p
from network_analysis import WATERS, load_graph, total_analysis, ligand_analysis

SCHEMA_VERSION=1

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def read_input(path):
    text=p.read_structure_text(str(path))
    cif=p.is_mmcif(text,str(path))
    cats=p.read_cif_categories(text,('_atom_site.','_entity.','_entity_poly.','_pdbx_unobs_or_zero_occ_residues.','_pdbx_struct_assembly_gen.')) if cif else {}
    if cif:
        site=cats['_atom_site'];models=set(site.get('pdbx_PDB_model_num',['1']))
        if len(models)>1: raise ValueError('Multiple coordinate models found; extract one model before analysis')
        alt=site.get('label_alt_id',['.']*len(site['Cartn_x']))
        selected=[i for i,a in enumerate(alt) if a in ('.','?','A','')]
        meta=[(float(site.get('B_iso_or_equiv',['0']*len(alt))[i]),float(site.get('occupancy',['1']*len(alt))[i])) for i in selected]
    else:
        if sum(l.startswith('MODEL ') for l in text.splitlines())>1: raise ValueError('Multiple coordinate models found; extract one model before analysis')
        lines=[l for l in text.splitlines() if l[:6] in ('ATOM  ','HETATM') and l[16:17].strip() in ('','A')]
        meta=[(float(l[60:66].strip() or 0),float(l[54:60].strip() or 1)) for l in lines]
    parsed=p.parse_structure(text,str(path))
    if len(meta)!=len(parsed['atoms']): raise ValueError('Source atom metadata does not match parsed atoms')
    for a,(b,occ) in zip(parsed['atoms'],meta):a['b_iso']=b;a['occupancy']=occ
    waters=sum(a['resName'] in WATERS for a in parsed['atoms'])
    parsed['atoms']=[a for a in parsed['atoms'] if a['resName'] not in WATERS]
    missing=cats.get('_pdbx_unobs_or_zero_occ_residues',{})
    missing_rows=[dict(chain=c,residue_name=r,residue_number=int(n)) for c,r,n in zip(missing.get('auth_asym_id',[]),missing.get('auth_comp_id',[]),missing.get('auth_seq_id',[]))]
    entities=cats.get('_entity',{});descriptions=dict(zip(entities.get('id',[]),entities.get('pdbx_description',[])))
    ep=cats.get('_entity_poly',{});chains=[]
    for i,eid in enumerate(ep.get('entity_id',[])):
        for chain in ep.get('pdbx_strand_id',[])[i].split(','):
            seq=ep.get('pdbx_seq_one_letter_code_can',ep.get('pdbx_seq_one_letter_code',[]))[i]
            chains.append(dict(chain=chain,description=descriptions.get(eid,''),deposited_sequence_length=len(re.sub(r'\s','',seq))))
    return parsed,dict(water_atoms_removed=waters,missing_residues=missing_rows,chains=chains,
        assembly_operators=len(parsed['biomt']),assembly_id=(cats.get('_pdbx_struct_assembly_gen',{}).get('assembly_id') or ['1'])[0],
        alternate_location_policy='Pampano: blank or A; other alternatives excluded',model_policy='single coordinate model required')

def run_job(job,out_root,viewer_js=None,reuse=False):
    if not job.get('input'):
        text,filename=p.fetch_entry(job['pdb_id'])
        with tempfile.TemporaryDirectory(prefix='pampano_input_') as tmp:
            source=Path(tmp)/Path(filename).name;source.write_text(text)
            return run_job(dict(job,input=str(source)),out_root,viewer_js,reuse)
    name=job.get('id') or job.get('pdb_id') or Path(job['input']).stem
    if not re.fullmatch(r'[A-Za-z0-9_.-]+',name) or name in {'.','..'}:raise ValueError('Job id must be a safe filename')
    out=Path(out_root);out.mkdir(parents=True,exist_ok=True)
    for sub in ['inputs','networks','total','ligand','structure']: (out/sub).mkdir(exist_ok=True)
    if job.get('input'):
        source=Path(job['input']).resolve()
        saved=out/'inputs'/source.name
    parsed,metadata=read_input(source)
    assembly=job.get('assembly','bioassembly')
    if assembly not in ('bioassembly','asym'):raise ValueError('Use bioassembly or asym; projected biocontacts cannot provide copy-resolved structure painting')
    if assembly=='bioassembly' and not parsed['biomt']: raise ValueError('No biological assembly operators found; explicitly select assembly=asym if appropriate')
    cutoff=float(job.get('cutoff',5.0));contact_cutoff=float(job.get('contact_cutoff',cutoff))
    if not (0<contact_cutoff<=cutoff):raise ValueError('Require 0 < contact_cutoff <= network cutoff')
    opts=p.default_opts();opts.update(assemblyMode=assembly,granularity='residue',atomScope='heavy',includeHetatm=True,cutoff=cutoff)
    if job.get('ligand_group','residue') not in ('residue','chain'):raise ValueError('ligand_group must be residue or chain')
    settings=dict(input_sha256=sha(source),assembly=assembly,cutoff=cutoff,contact_cutoff=contact_cutoff,
        ligands=job.get('ligands'),ligand_chain=job.get('ligand_chain'),ligand_residue=job.get('ligand_residue'),
        ligand_group=job.get('ligand_group','residue'),betweenness_sources=int(job.get('betweenness_sources',256)),seeds=[17,29],
        core_sha256=sha(Path(p.__file__)),schema_version=SCHEMA_VERSION)
    if job.get('partner_chains'):settings['partner_chains']=job['partner_chains']
    if job.get('interface_left') or job.get('interface_right'):
        if not job.get('interface_left') or not job.get('interface_right'):raise ValueError('Both interface sides are required')
        settings['interface_left']=job['interface_left'];settings['interface_right']=job['interface_right']
    if settings['betweenness_sources']<1:raise ValueError('betweenness_sources must be positive')
    previous=out/'manifest.json'
    cached=reuse and previous.exists() and json.loads(previous.read_text()).get('settings')==settings
    if previous.exists() and not cached:raise ValueError('Output already contains a run; use a new directory, or --reuse with identical settings')
    if source!=saved.resolve():shutil.copy2(source,saved)
    atoms=p.prepare_atoms(parsed,opts)
    if not cached:
        print(f'{name}: building {assembly}, {len(atoms)} atoms',flush=True)
        graph=p.build_graph(atoms,opts)
        paths=p.write_outputs(graph,opts,str(out/'networks/contact_network'),saved.name+'; waters excluded')
        del graph
        for path in paths:
            with open(path,'rb') as src,gzip.open(path+'.gz','wb') as dst:shutil.copyfileobj(src,dst)
            Path(path).unlink()
    g=load_graph(out/'networks/contact_network.json.gz')
    if not len(g):raise ValueError('No graph residues survived the selection')
    if not cached:
        ligand_scores,ligand_rank,site_paths,ligand_summary=ligand_analysis(g,out/'ligand',ligands=settings['ligands'],
            ligand_chain=settings['ligand_chain'],ligand_residue=settings['ligand_residue'],group=settings['ligand_group'],contact_cutoff=contact_cutoff,partner_chains=settings.get('partner_chains'))
        total_nodes,total_res,total_summary=total_analysis(g,out/'total',k=settings['betweenness_sources'])
    else:
        total_nodes=pd.read_csv(out/'total/node_metrics.csv.gz',dtype={'id':str},keep_default_na=False).set_index('id',drop=False)
        total_res=pd.read_csv(out/'total/residue_metrics.csv',keep_default_na=False)
        total_summary=json.loads((out/'total/summary.json').read_text())
        ligand_scores=pd.read_csv(out/'ligand/site_residue_scores.csv',dtype={'node_id':str},keep_default_na=False)
        ligand_rank=pd.read_csv(out/'ligand/residue_ranking.csv',keep_default_na=False)
        ligand_summary=json.loads((out/'ligand/summary.json').read_text())
        with gzip.open(out/'ligand/site_paths.json.gz','rt') as f:site_paths=json.load(f)
    pd.DataFrame(metadata['missing_residues'],columns=['chain','residue_name','residue_number']).to_csv(out/'missing_residues.csv',index=False)
    pd.DataFrame(metadata['chains']).to_csv(out/'chain_metadata.csv',index=False)
    manifest=dict(structure=name,settings=settings,metadata=metadata,software={'pampano':p.__version__},
        source_files={f:sha(Path(__file__).parent/f) for f in ['analyze_structure.py','network_analysis.py','pampano_offline.py']})
    if cached:manifest['source_files']=json.loads(previous.read_text())['source_files']
    if settings.get('interface_left'):
        from interface_analysis import interface_analysis
        site_paths,interface_summary=interface_analysis(g,out/'interface',settings['interface_left'],settings['interface_right'],contact_cutoff)
        manifest['interface_summary']=interface_summary
        manifest['interface_source_sha256']=sha(Path(__file__).parent/'interface_analysis.py')
    previous.write_text(json.dumps(manifest,indent=2)+'\n')
    if not job.get('no_views',False):
        from structure_view import export_views
        export_views(g,atoms,parsed['atoms'],total_nodes,site_paths,out/'structure',name,viewer_js,
            radius=float(job.get('view_radius',12)),max_view_sites=int(job.get('max_view_sites',120)),global_exact=total_summary['betweenness_exact'])
    print(name, json.dumps(dict(total=total_summary['graphs'],ligand=ligand_summary)),flush=True)
    return dict(id=name,status='ok',output=str(out),nodes=len(g),edges=g.number_of_edges(),ligand_sites=ligand_summary['site_count'])

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('input',nargs='?');ap.add_argument('--pdb-id');ap.add_argument('--batch',type=Path)
    ap.add_argument('--out',type=Path,required=True);ap.add_argument('--id')
    ap.add_argument('--ligand',action='append',help='Residue name; repeat for alternatives. Default: every non-water, non-protein HETATM.')
    ap.add_argument('--ligand-chain');ap.add_argument('--ligand-residue',type=int)
    ap.add_argument('--partner-chain',dest='partner_chains',action='append',help='Explicit non-protein author chain, including RNA/DNA ATOM records. Use --ligand-group chain for whole-chain endpoints.')
    ap.add_argument('--interface-left',action='append',help='Author protein chain on one interface side; repeat to select several')
    ap.add_argument('--interface-right',action='append',help='Author protein chain on the other interface side')
    ap.add_argument('--ligand-group',choices=['residue','chain'],default='residue')
    ap.add_argument('--assembly',choices=['bioassembly','asym'],default='bioassembly')
    ap.add_argument('--cutoff',type=float,default=5);ap.add_argument('--contact-cutoff',type=float)
    ap.add_argument('--betweenness-sources',type=int,default=256)
    ap.add_argument('--viewer-js',type=Path,help='Embed a local 3Dmol.js bundle; otherwise use pinned CDN')
    ap.add_argument('--view-radius',type=float,default=12);ap.add_argument('--max-view-sites',type=int,default=120)
    ap.add_argument('--no-views',action='store_true');ap.add_argument('--reuse',action='store_true')
    a=ap.parse_args()
    if a.batch:
        if a.input or a.pdb_id:ap.error('--batch cannot be combined with a single input')
        data=json.loads(a.batch.read_text());jobs=data['jobs'] if isinstance(data,dict) else data
        ids=[j.get('id') or j.get('pdb_id') or Path(j['input']).stem for j in jobs]
        if len(ids)!=len(set(ids)):ap.error('Batch job ids must be unique')
        results=[]
        for job,name in zip(jobs,ids):
            try:
                if not re.fullmatch(r'[A-Za-z0-9_.-]+',name) or name in {'.','..'}:raise ValueError('Unsafe batch job id')
                if job.get('input'):job=dict(job,input=str((a.batch.parent/job['input']).resolve()))
                results.append(run_job(job,a.out/name,a.viewer_js,a.reuse))
            except Exception as e:results.append(dict(id=name,status='error',error=str(e)))
            a.out.mkdir(parents=True,exist_ok=True);(a.out/'batch_results.json').write_text(json.dumps(results,indent=2)+'\n')
        if any(r['status']=='error' for r in results):sys.exit(1)
    else:
        if bool(a.input)==bool(a.pdb_id):ap.error('Provide exactly one input path or --pdb-id')
        job={k:v for k,v in vars(a).items() if v is not None and k not in {'batch','out','viewer_js','reuse','ligand'}}
        job['ligands']=a.ligand
        run_job(job,a.out,a.viewer_js,a.reuse)

if __name__=='__main__':main()
