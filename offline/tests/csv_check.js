// usage: node csv_check.js ../../index.html file.pdb <assembly> <granularity> out.csv
// builds the network through the real UI (#buildBtn) and captures the CSV the Export button downloads
const fs = require('fs'), path = require('path');
const { chromium } = require('playwright');
(async () => {
  const [html, pdb, mode, gran, out] = process.argv.slice(2);
  const browser = await chromium.launch();
  const page = await browser.newPage({ acceptDownloads: true });
  await page.route('**/d3.min.js', r => r.fulfill({ path: path.join(__dirname, 'node_modules/d3/dist/d3.min.js') }));
  await page.route('**/three.min.js', r => r.fulfill({ path: path.join(__dirname, 'node_modules/three/build/three.min.js') }));
  await page.goto('file://' + path.resolve(html));
  await page.setInputFiles('#fileInput', pdb);
  await page.selectOption('#assemblyMode', mode);
  await page.selectOption('#granularity', gran);
  await page.waitForFunction(() => !document.getElementById('buildBtn').disabled);
  await page.click('#buildBtn');
  await page.waitForFunction(() => !document.getElementById('exportBtn').disabled, null, { timeout: 120000 });
  const [dl] = await Promise.all([page.waitForEvent('download'), page.click('#exportBtn')]);
  await dl.saveAs(out);
  await browser.close();
})();
