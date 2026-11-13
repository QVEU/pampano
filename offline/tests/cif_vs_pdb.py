"""Check that the PDB and mmCIF files of the same entries give the same network.
usage: python cif_vs_pdb.py 1crn,4hhb,3vbs PDB_DIR/   (needs pdbXXXX.ent and XXXX.cif in PDB_DIR)"""
import sys, itertools
sys.path.insert(0, __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), ".."))
import pampano_offline as po
D = sys.argv[2].rstrip("/") + "/" if len(sys.argv) > 2 else "pdb/"

def graph(parsed, o):
    g = po.build_graph(po.prepare_atoms(parsed, o), o, verbose=False)
    return po.browser_style_graph(g)

def keyed(g):
    n = g['nodes']
    ed = {}
    for l in g['links']:
        a, b = n[l['source']]['title'], n[l['target']]['title']
        flip = a > b
        k = (b, a) if flip else (a, b)
        sp, tp = (l['targetPart'], l['sourcePart']) if flip else (l['sourcePart'], l['targetPart'])
        pc = {}
        for kk, v in l['partCounts'].items():
            x, y = kk.split('-'); pc[(y + '-' + x) if flip else kk] = v
        ed[k] = (l['type'], l['scope'], l['count'], l['bond'], sp, tp, tuple(sorted(pc.items())), round(l['dist'], 3))
    nodes = {x['title']: (x['label'], x['chain'], x['baseChain'], x['hetatm'], None if x['sasa'] is None else round(x['sasa'], 1)) for x in n}
    return nodes, ed

ids = sys.argv[1].split(',')
tot = bad = 0
for i in ids:
    P = po.parse_pdb(open(D + 'pdb%s.ent' % i).read())
    C = po.parse_mmcif(open(D + '%s.cif' % i).read())
    for mode, gran, het in itertools.product(['asym', 'bioassembly', 'biocontacts'], ['residue', 'atom'], [False, True]):
        if gran == 'atom' and mode == 'bioassembly' and i in ('3vbs', '2plv', '9mxc', '1aon', '6vxx'): continue
        o = dict(po.default_opts(), assemblyMode=mode, granularity=gran, includeHetatm=het)
        gp, gc = graph(P, o), graph(C, o)
        tot += 1
        if gp == gc:
            print('IDENTICAL', i, mode, gran, 'het' if het else '', len(gp['nodes']), len(gp['links'])); continue
        np_, ep = keyed(gp); nc, ec = keyed(gc)
        same = np_ == nc and ep == ec
        if not same: bad += 1
        print('SAME-SET ' if same else 'DIFF     ', i, mode, gran, 'het' if het else '', len(gp['nodes']), len(gp['links']), len(gc['nodes']), len(gc['links']))
        if not same:
            only_p = set(ep) - set(ec); only_c = set(ec) - set(ep)
            chg = [k for k in set(ep) & set(ec) if ep[k] != ec[k]]
            print('    edges only PDB', len(only_p), list(only_p)[:2], ' only CIF', len(only_c), list(only_c)[:2], ' changed', len(chg))
            for k in chg[:2]: print('     ', k, ep[k], ec[k])
            nd = [k for k in set(np_) | set(nc) if np_.get(k) != nc.get(k)]
            print('    node diffs', len(nd), [(k, np_.get(k), nc.get(k)) for k in nd[:2]])
print('%d/%d configs same network' % (tot - bad, tot))
