"""Compare the offline tool with the browser dump.
usage: python parity.py configs.json browser.json  (prints OK/DIFF per config)"""
import json, sys, itertools, math, time
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import pampano_offline as po

def py_dump(cfgs):
    out, cache = [], {}
    for c in cfgs:
        if c['pdb'] not in cache:
            cache[c['pdb']] = po.parse_pdb(open(c['pdb']).read())
        o = dict(po.default_opts(), **c['opts'])
        o['protomerGroups'] = po.parse_protomer_groups(c.get('groups', ''))
        o['sasa'] = bool(c.get('sasa'))
        t = time.time()
        g = po.build_graph(po.prepare_atoms(cache[c['pdb']], o), o, verbose=False)
        r = po.browser_style_graph(g); r['ms'] = (time.time()-t)*1000
        out.append(r)
    return out

def same(a, b, tol):
    if isinstance(a, float) or isinstance(b, float):
        if a is None or b is None: return a is b
        return abs(a-b) <= tol*max(1, abs(a))
    return a == b

def compare(br, py):
    bad = []
    if br['meta'] != py['meta']: bad.append(('meta', br['meta'], py['meta']))
    if len(br['nodes']) != len(py['nodes']): bad.append(('nnodes', len(br['nodes']), len(py['nodes'])))
    if len(br['links']) != len(py['links']): bad.append(('nlinks', len(br['links']), len(py['links'])))
    for a, b in zip(br['nodes'], py['nodes']):
        for k in a:
            if not same(a[k], b[k], 1e-9 if k in ('sasa','relsasa') else 0):
                bad.append(('node', a['id'], k, a[k], b[k])); break
        if len(bad) > 5: break
    for e, (a, b) in enumerate(zip(br['links'], py['links'])):
        for k in a:
            if not same(a[k], b[k], 0):
                bad.append(('link', e, k, a[k], b[k])); break
        if len(bad) > 10: break
    return bad

if __name__ == '__main__':
    cfgs = json.load(open(sys.argv[1])); br = json.load(open(sys.argv[2]))
    py = py_dump(cfgs)
    nbad = 0
    for c, b, p in zip(cfgs, br, py):
        bad = compare(b, p)
        tag = 'OK  ' if not bad else 'DIFF'
        if bad: nbad += 1
        print(tag, c['pdb'].split('/')[-1], json.dumps(c['opts']), c.get('groups',''), 'sasa' if c.get('sasa') else '',
              '%dn %de' % (len(b['nodes']), len(b['links'])), 'browser %.0fms py %.0fms' % (b['ms'], p['ms']))
        for x in bad[:6]: print('     ', x)
    print('%d/%d configs identical' % (len(cfgs)-nbad, len(cfgs)))
