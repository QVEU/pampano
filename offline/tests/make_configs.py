"""Write the parity-test configurations used for the offline tool.
usage: python make_configs.py PDB_DIR configs.json
PDB_DIR holds pdbXXXX.ent files (see README: Testing)."""
import itertools, json, os, sys

pdb_dir, out = sys.argv[1], sys.argv[2]
P = lambda i: os.path.join(os.path.abspath(pdb_dir), 'pdb%s.ent' % i)
base = {'cutoff': 5.0, 'atomScope': 'heavy'}
cfgs = []
small = [i for i in ['1crn', '4e43', '1hvr', '1a28', '4hhb', '1stp', '1osm'] if os.path.exists(P(i))]
for pdb in small:
    for mode, gran, scope in itertools.product(['monomer', 'asym', 'bioassembly', 'biocontacts'],
                                               ['residue', 'atom'], ['heavy', 'all', 'ca']):
        if gran == 'atom' and pdb in ('1osm', '1a28', '4hhb') and scope == 'all':
            continue
        cfgs.append({'pdb': P(pdb), 'opts': {'assemblyMode': mode, 'granularity': gran, 'atomScope': scope, 'cutoff': 5.0}})
    for extra in [{'includeHetatm': True}, {'excludeSeqAdj': True}, {'covalentOnly': True}, {'interChain': False},
                  {'intraChain': False}, {'interProtomer': False},
                  {'cutoff': 4.0, 'includeHetatm': True, 'excludeSeqAdj': True}]:
        for mode in ['asym', 'biocontacts']:
            cfgs.append({'pdb': P(pdb), 'opts': dict(base, assemblyMode=mode, granularity='residue', **extra)})
    cfgs.append({'pdb': P(pdb), 'opts': dict(base, assemblyMode='biocontacts', granularity='residue', includeHetatm=True), 'sasa': True})
if os.path.exists(P('4hhb')):
    cfgs.append({'pdb': P('4hhb'), 'opts': dict(base, assemblyMode='asym', granularity='residue'), 'groups': 'A,B ; C,D'})
if os.path.exists(P('1stp')):
    cfgs.append({'pdb': P('1stp'), 'opts': dict(base, assemblyMode='bioassembly', granularity='residue', includeHetatm=True), 'sasa': True})
for pdb in [i for i in ['3vbs', '2plv'] if os.path.exists(P(i))]:
    for mode in ['asym', 'biocontacts']:
        for gran in ['residue', 'atom']:
            cfgs.append({'pdb': P(pdb), 'opts': dict(base, assemblyMode=mode, granularity=gran)})
    cfgs.append({'pdb': P(pdb), 'opts': dict(base, assemblyMode='biocontacts', granularity='residue', includeHetatm=True), 'sasa': True})
    cfgs.append({'pdb': P(pdb), 'opts': dict(base, assemblyMode='bioassembly', granularity='residue')})
json.dump(cfgs, open(out, 'w'))
print(len(cfgs), 'configs')
