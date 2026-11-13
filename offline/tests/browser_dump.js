// usage: node browser_dump.js ../../index.html configs.json browser.json
// Loads index.html headless (d3/three served from node_modules), runs parsePDB/prepareAtoms/buildGraph
// for each config and dumps the graphs for parity.py.
const fs = require('fs'), path = require('path');
const { chromium } = require('playwright');
(async () => {
  const [html, cfgPath, outPath] = process.argv.slice(2);
  const cfgs = JSON.parse(fs.readFileSync(cfgPath, 'utf8'));
  const browser = await chromium.launch();
  const page = await browser.newPage();
  await page.route('**/d3.min.js', r => r.fulfill({ path: path.join(__dirname, 'node_modules/d3/dist/d3.min.js'), contentType: 'text/javascript' }));
  await page.route('**/three.min.js', r => r.fulfill({ path: path.join(__dirname, 'node_modules/three/build/three.min.js'), contentType: 'text/javascript' }));
  await page.goto('file://' + path.resolve(html));
  const results = [];
  let lastPdb = null;
  for (const c of cfgs) {
    if (c.pdb !== lastPdb) {
      const text = fs.readFileSync(c.pdb, 'utf8');
      await page.evaluate(t => { parsed = parsePDB(t); }, text);
      lastPdb = c.pdb;
    }
    const r = await page.evaluate(o => {
      const opts = Object.assign({colorBy: o.sasa ? 'sasa' : 'chain', sizeBy: 'uniform', edgeFilter: 'all', edgeColorBy: 'type', assemblyMode:'asym', granularity:'residue', cutoff:5, atomScope:'heavy', covalentOnly:false, intraChain:true, interChain:true, interProtomer:true, excludeSeqAdj:false, includeHetatm:false}, o.opts);
      opts.protomerGroups = parseProtomerGroups(o.groups || '');
      const t0 = performance.now();
      const atoms = prepareAtoms(opts);
      const g = buildGraph(atoms, opts);
      const ms = performance.now() - t0;
      return { ms,
        nodes: g.nodes.map(n => ({id:n.id, label:n.label, chain:n.chain, baseChain:n.baseChain, opId:n.opId, copyId:n.copyId,
          hetatm:n.hetatm, resName:n.resName, resSeq:Number.isNaN(n.resSeq)?null:n.resSeq, iCode:n.iCode, hydro:n.hydro, vol:n.vol,
          sasa:n.sasa, relsasa:n.relsasa, resclass:n.resclass, title:n.title})),
        links: g.links.map(l => ({source:l.source, target:l.target, dist:l.dist, count:l.count||1, type:l.type, scope:l.scope,
          bond:!!l.bond, sourcePart:l.sourcePart, targetPart:l.targetPart, partCounts:l.partCounts})),
        meta: {mateAtoms: g.meta.mateAtoms||0, selfSymPairs: g.meta.selfSymPairs||0} };
    }, c);
    results.push(r);
    process.stderr.write(`${path.basename(c.pdb)} ${JSON.stringify(c.opts)} -> ${r.nodes.length}n ${r.links.length}e ${r.ms.toFixed(0)}ms\n`);
  }
  fs.writeFileSync(outPath, JSON.stringify(results));
  await browser.close();
})();
