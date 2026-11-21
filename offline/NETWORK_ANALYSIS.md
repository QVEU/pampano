# Total networks, binding partners, and structural score maps

`analyze_structure.py` adds NetworkX analysis and interactive structure coloring
to the existing Pampano offline contact builder. The contact builder and its
interaction definitions are unchanged. Python 3.10+ is recommended.

```bash
python -m pip install -r offline/requirements-analysis.txt
python offline/analyze_structure.py --pdb-id 3VBS --ligand SPH --out results/3VBS
python offline/analyze_structure.py --pdb-id 9DTT --partner-chain B --ligand-group chain --betweenness-sources 100000 --out results/9DTT
python offline/analyze_structure.py --pdb-id 8T12 --interface-left B --interface-right A --betweenness-sources 100000 --out results/8T12
python offline/analyze_structure.py --batch offline/analysis_batch.example.json --out results
```

A positional PDB/mmCIF path can replace `--pdb-id`. Batch `input` paths are
relative to the JSON file. Each job supplies its own analysis options; the
batch-level flags `--viewer-js` and `--reuse` apply to every job. Jobs run
sequentially, and `batch_results.json` records success or failure for each job.
The process exits nonzero if any job fails. Use separate output directories for
different selections or cutoffs. `--reuse` requires identical analysis settings
and input/core hashes and can regenerate visualizations from saved metrics.
It does not revalidate changed analysis implementation: after upgrading analysis
code, use a new directory to recompute results.

## Graph and coordinate policies

- Default: first deposited biological assembly, residue nodes, heavy atoms,
  5 Å minimum atom-pair cutoff, all non-water residues retained. Explicitly use
  `--assembly asym` when asymmetric-unit analysis is intended.
- Water names HOH/WAT/DOD are removed. Pampano's blank/A alternate-location rule
  is retained; multiple coordinate models are rejected. No atoms, loops, missing
  residues, protonation states, or side chains are modeled or refined.
- `--cutoff` controls the graph. `--contact-cutoff` can select a smaller contact
  endpoint set, but cannot exceed the graph cutoff. All path calculations are
  undirected and unweighted; atom distances are not path weights.
- `biocontacts` collapses physical symmetry copies. It can alter shortest paths
  and is unsuitable as a substitute for copy-resolved global centrality or
  structure painting. Use the original Pampano CLI for compact projections and
  compare them explicitly if needed.
- Modified amino acids recognized by Gemmi enter the protein graph. Other
  nonprotein residues, including nucleic acids and metals, do not enter the
  protein-only path calculations below. This is a topology model, not an
  energetic model or a demonstration of biological signaling.

## Total-network analysis

Outputs include degree, local clustering, core number, interface degree,
connected components, interaction counts, and Louvain communities. Global
betweenness uses 256 sampled sources at each of seeds 17 and 29 and their mean.
Set `--betweenness-sources` at least as large as the graph to compute exact
global betweenness. `summary.json` records which was used. Symmetry-equivalent
residues are averaged by author chain, residue number, insertion code, and name;
physical-copy tables remain available. Symmetry copies are not independent
experimental replicates. Mean path length samples up to 64 sources in the
largest component; the reported diameter is a lower bound.

The total graph includes ligand/RNA/metal nodes, which can form paths in this
analysis. Protein-only topology and topology after removing peptide-labeled
edges are reported separately; the latter is not a reconstruction of every
noncovalent atom interaction.

## Ligands and RNA/DNA partners

Without selectors, every non-water, nonprotein HETATM residue becomes a site.
`--ligand` can be repeated; `--ligand-chain` and `--ligand-residue` filter author
identities. Use `--partner-chain B --ligand-group chain` to select an explicit
nonprotein polymer chain, including RNA/DNA ATOM records. Chain grouping unites
the selected residues of each physical chain/operator into a single partner.
Check entity identities first: author chain letters have no universal meaning.

For each site, C is the set of protein residues contacting the selected partner.
The score is the fraction of shortest paths that use residue v internally,
averaged over all unordered pairs in C. Tied paths share credit. Endpoints do
not receive credit for their own pair, but contact residues may be internal to
other pairs. All nonprotein nodes are removed before finding paths, preventing
RNA/ligand shortcuts. Unreachable pairs contribute zero and are counted;
fewer than two contacts is reported as insufficient. Nodes on every tied
shortest path are retained, and the local induced graph is checked by the total
internal path-credit identity. The denominator is choose(|C|, 2), not all
protein pairs. This score is not numerically comparable with global centrality.

`ligand/direct_contacts.csv` includes both protein and partner residue
identities, distance, and site. `site_residue_scores.csv` is copy-resolved;
`residue_ranking.csv` sums contributions from physical copies of each base
residue within a site, then averages over valid sites with the same partner
label. For heterogeneous binding sites, use site-level tables. The aggregate
`direct_ligand_contact` flag means contact in any contributing site.

## Protein-protein interfaces

Repeat `--interface-left` and `--interface-right` to select two disjoint sets of
author protein chains. Each contacting physical chain pair is analyzed
separately. The endpoints are its contacting residues on the two sides, and
paths use the full protein-only complex. The score averages internal path
participation over |left contacts| × |right contacts| cross-interface pairs.
The implementation corrects NetworkX's half-path convention for disjoint
source/target sets in an undirected graph. Direct contact edges, per-residue
interface degree, and exact centralities are in `interface/`. The selected
interface is painted in the viewer; small-molecule analyses remain in `ligand/`.

## Structure coloring and exports

Open `structure/structure_view.html` in a WebGL-capable browser. Select a contact
site or total network, click residues for values, and adjust the numeric color
maximum. Light gray → yellow → red represents increasing score. Selected RNA or
small ligands are green, and direct protein-contact CAs have cyan markers. For
interfaces both protein partners use the score colors. The total view shows
protein alpha carbons only; RNA and metals remain part of its underlying graph.
The contact view includes a 12 Å neighborhood and every path residue; its
visual crop never restricts the graph used for calculations. At most 120 sites
are embedded by default (`--max-view-sites`); tables cover every site.

The viewer loads pinned 3Dmol.js 2.5.5 from its CDN unless `--viewer-js PATH`
embeds a downloaded bundle. Embedded bundles make the HTML work offline;
distribute the renderer's license with any bundled version. Visual covalent
bonds are inferred from element radii solely for display; graph edges and
scores always come from Pampano/NetworkX. `--no-views` skips rendering exports.

New visualization-only mmCIF files carry **100 × score in B_iso_or_equiv**:

- `total_CA_score_percent.cif`: protein alpha-carbon global score map.
- `ligand_site_000_score_percent.cif`: all atoms for the first valid contact
  site (the historical filename also applies to RNA and interface analyses).
  `ligand_site_000_scores.json` identifies the actual site and definition.
- `paint_total.pml` and `paint_ligand.pml`: starting PyMOL coloring commands.

These B values are not experimental displacement factors. Original B values,
graph node IDs and physical copies are retained in a custom mmCIF mapping loop.
Source coordinates are unchanged, and exports are read back to check atom
counts and coordinate identity. Mapping preserves author chain, number,
insertion code, residue name and copy. Duplicate IDs get explicit display-chain
suffixes. Repeated runs on different proteins do not establish homologous
residue correspondence; sequence/structural alignment is a separate step.

```bash
python -m pip install pytest
python -m pytest offline/tests/test_network_analysis.py
```
