# Pampano offline

A command-line version of Pampano for structures that are too big for the browser: virus capsids, large biological assemblies, atom-level networks of big complexes. It builds **the same network as the web app** (same PDB parsing, BIOMT handling, assembly modes, atom selection, contact search, interaction types, scope labels, backbone/sidechain partners and SASA) and writes it to files instead of drawing it:

| File | What it is | Open it with |
|---|---|---|
| `PREFIX.json` | the network in Cytoscape.js JSON | Cytoscape desktop, Cytoscape.js, or any JSON reader (Python, R, D3, Gephi via plugin…) |
| `PREFIX_edges.csv` | one row per edge (contact or bond) | R, pandas, Excel |
| `PREFIX_nodes.csv` | one row per node (residue or atom) | R, pandas, Excel |

For example, the full 3VBS capsid (60 copies, 50,520 residues, 312,120 edges) takes about 15 seconds and 1.2 GB of memory at residue level.

## Install

You need Python 3.8 or newer with **numpy** and **scipy**. If you already have them, there is nothing to install: just run the script.

Otherwise, use the conda environment in this folder:

```bash
conda env create -f environment.yml
conda activate pampano
```

or pip: `pip install -r requirements.txt`.

## Quick start

```bash
cd offline

# residue network of the asymmetric unit
python pampano_offline.py 1abc.pdb

# full biological assembly (every symmetry copy gets its own nodes)
python pampano_offline.py 3vbs.pdb -a bioassembly -o 3vbs_capsid

# one copy of each chain, with every contact to neighbouring copies mapped back
python pampano_offline.py 3vbs.pdb -a biocontacts -o 3vbs_asu

# fetch from RCSB instead of a local file (needs internet)
python pampano_offline.py --pdb-id 4hhb

# a small form instead of the command line
python pampano_offline.py --gui
```

Progress goes to the terminal; `-q` silences it. Input can be `.pdb`, `.ent` or either gzipped (`.pdb.gz`, `.ent.gz`, as downloaded from the PDB archive). Like the web app, only PDB format is read, not mmCIF.

## Options

Each option matches a control in the web app.

| Option | Web app control | Default |
|---|---|---|
| `-a, --assembly monomer\|asym\|bioassembly\|biocontacts` | Assembly: monomer (first chain), full asymmetric unit, biological assembly (applies REMARK 350 BIOMT of assembly 1), bio-assembly contacts on asym unit | `asym` |
| `-g, --granularity residue\|atom` | Granularity | `residue` |
| `-c, --cutoff Å` | Contact cutoff | `5.0` |
| `--atoms all\|heavy\|ca` | Atoms: all, heavy only, CA/representative | `heavy` |
| `--covalent-only` | Covalent bonds only | off |
| `--no-intra-chain` | untick Include intra-chain contacts | included |
| `--no-inter-chain` | untick Include inter-chain contacts | included |
| `--no-inter-protomer` | untick Include inter-protomer contacts | included |
| `--protomer-groups "A,B ; C,D"` | Protomer groups | none (classify by sequence) |
| `--exclude-seq-adjacent` | Exclude sequential neighbors | off |
| `--hetatm` | Include HETATM | off |
| `--sasa` | Color/size by SASA (computes SASA) | off |
| `--edge-filter FILTER` | Show edges: `all`, `bonds`, `contacts`, `peptide`, `disulfide`, `salt_bridge`, `hydrogen_bond`, `hydrophobic`, `scope_intra`, `scope_intra_protomer`, `scope_homotypic`, `scope_heterotypic`, `scope_inter_protomer` | `all` |
| `-o, --out PREFIX` | output file prefix | `<file>_<assembly>_<granularity>` |

`--edge-filter` writes only the matching edges, and node degree is counted on those edges, just as the web app's degree follows the shown edges. Nodes are kept either way. For analysis it is usually simpler to write all edges and filter in R.

SASA is the slowest step: about 25 seconds for the 390,000 heavy atoms of the 3VBS capsid. In `biocontacts` mode SASA includes the neighbouring symmetry copies, so interface burial is real.

### Which assembly mode for a capsid?

- **`bioassembly`** gives every copy its own nodes: 60 × the asymmetric unit for an icosahedral capsid. Use it when you want the whole particle as one network (e.g. paths across the shell).
- **`biocontacts`** gives the size of one asymmetric unit, with each contact to a neighbouring copy credited back to the residue it touches (scope `inter_protomer`, or `homotypic`/`heterotypic` for different chains). Use it to ask "which residues make which interfaces". Contacts between a residue and its own symmetry copy have no edge (it would be a self-loop); their count is stored as `self_symmetry_atom_pairs` in the JSON.

## Output formats

### JSON (Cytoscape.js)

```json
{
  "format_version": "1.0",
  "generated_by": "pampano_offline 1.0",
  "data": { "name": "Pampano 3vbs.pdb", "assembly_mode": "biocontacts", "cutoff": 5.0, ... },
  "elements": {
    "nodes": [ { "data": { "id": "0", "name": "A GLY 1", "chain": "A", "res_name": "GLY", "degree": 7, ... },
                 "position": { "x": 1262.3, "y": 2237.4 } }, ... ],
    "edges": [ { "data": { "id": "e0", "source": "0", "target": "11", "interaction": "hydrogen_bond",
                           "scope": "intra", "distance": 2.87, "count": 9, ... } }, ... ]
  }
}
```

`data` at the top records every setting used, so a file documents how it was made.

**Node attributes:** `id`, `name` (e.g. `A GLY 1`, the label shown in Cytoscape), `label` (`GLY1`), `chain` (copy label: `A`, or `A_5` for BIOMT operator 5), `base_chain` (deposited chain), `operator` (BIOMT index), `copy` (distinguishes repeated residues with the same ID), `res_name`, `res_seq`, `icode`, `atom_name` (atom level only), `hetatm`, `res_class` (neg/pos/polar/hydrophobic/special/unknown), `hydropathy` (Kyte–Doolittle), `volume` (Å³), `sasa` and `rel_sasa` (with `--sasa`), `degree`, and `x`, `y`, `z` (residue centroid or atom position, Å).

**Edge attributes:** `source`, `target`, `name`, `interaction` (peptide, disulfide, covalent, salt_bridge, hydrogen_bond, hydrophobic, vdw), `scope` (intra, intra_protomer, homotypic, heterotypic, inter_protomer), `distance` (closest atom pair, Å), `count` (atom pairs), `is_bond`, `source_part` / `target_part` (backbone, sidechain or hetatm: the parts of the atom pair that defines the interaction type), `part_counts` (e.g. `sc-sc:12;bb-sc:3`) and one `pc_<source>_<target>` count per part combination present (e.g. `pc_sc_sc`). See the web app's Help for definitions.

Missing values are left out rather than written as `null`, so Cytoscape does not create empty columns.

**Opening in Cytoscape desktop:** File → Import → Network from File…, pick the `.json`. Nodes are placed at their structure x/y coordinates (×10), so the network starts as a projection of the structure; apply any layout from the Layout menu to rearrange. Every attribute appears in the Node Table and Edge Table and can drive a Style (e.g. map `interaction` to edge color, `degree` to node size). For networks over a few hundred thousand edges, give Cytoscape more memory (Cytoscape → Preferences → Memory, or `-Xmx` in `Cytoscape.vmoptions`) and turn on the graphics-detail threshold.

**Other tools:** in Cytoscape.js, `cytoscape({ elements: json.elements })`. In Python:

```python
import json, networkx as nx
j = json.load(open("3vbs_asu.json"))
g = nx.Graph()
g.add_nodes_from((n["data"]["id"], n["data"]) for n in j["elements"]["nodes"])
g.add_edges_from((e["data"]["source"], e["data"]["target"], e["data"]) for e in j["elements"]["edges"])
```

### CSV edge list

`PREFIX_edges.csv` has one row per edge. The first 12 columns are exactly the web app's CSV export:

| Column | Meaning |
|---|---|
| `source`, `target` | node ids (0-based; match `id` in the nodes CSV) |
| `source_label`, `target_label` | e.g. `A GLY 1` |
| `distance` | closest atom pair, Å (2 decimals) |
| `count` | number of atom pairs in contact |
| `interaction` | strongest interaction type |
| `scope` | structural scope |
| `is_bond` | 1 for covalent bonds (peptide, disulfide, covalent), else 0 |
| `source_part`, `target_part` | backbone / sidechain / hetatm |
| `part_counts` | atom pairs by part, largest first |

followed by:

- `pc_bb_bb`, `pc_bb_sc`, … `pc_het_het`: the same part counts as nine numeric columns (zero when absent), source part first.
- `source_chain`, `source_base_chain`, `source_operator`, `source_res_name`, `source_res_seq`, `source_icode`, `source_hetatm`, `source_res_class`, and the same for `target_` (plus `atom_name` at atom level), so most analyses need no join.

`PREFIX_nodes.csv` has one row per node with the node attributes listed above. Missing values are `NA` and booleans are `TRUE`/`FALSE`, which R reads natively.

### Reading the CSVs in R

```r
library(readr); library(dplyr)

edges <- read_csv("3vbs_asu_edges.csv", col_types = cols(source_icode = "c", target_icode = "c"))
nodes <- read_csv("3vbs_asu_nodes.csv", col_types = cols(icode = "c"))

# contacts between chains (and symmetry copies), by chain pair and type
edges |>
  filter(scope != "intra") |>
  count(source_chain, target_chain, interaction, sort = TRUE)

# hydrogen bonds made by side chains, per residue
edges |>
  filter(interaction == "hydrogen_bond", source_part == "sidechain") |>
  count(source_label, sort = TRUE)

# as a graph
library(igraph)
g <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)
```

## How it stays identical to the web app

The script is a line-by-line port of `buildGraph` and its helpers in `../index.html`, vectorised with numpy and a scipy KD-tree. It even replays the browser's order of atom pairs, so ties (two atom pairs at exactly the same distance) resolve the same way and edges come out in the same order.

It was checked against the web app running in headless Chromium (`tests/`): 275 configurations across 1CRN, 4E43, 1HVR, 1A28, 4HHB, 1STP, 1OSM, 2PLV and 3VBS (every assembly mode × residue/atom × all/heavy/CA atoms, plus HETATM, sequential-neighbour exclusion, covalent-only, each scope switched off, protomer groups and SASA, and the full 3VBS and 2PLV assemblies of about 50,000 residues and 310,000 edges each) produced identical nodes and edges: same ids, labels, types, scopes, parts, part counts and bit-identical distances, with SASA equal to within 1e-9 Å². The edge CSV also matches the web app's Export CSV download byte for byte in its 12 columns.

**If you change `index.html`'s network logic, make the same change here** and rerun the tests.

### Running the tests

```bash
cd offline/tests
npm install                    # d3, three, playwright (for headless Chromium)
npx playwright install chromium   # if you don't already have it
# PDB files named pdbXXXX.ent in some folder, e.g. from https://files.rcsb.org/download/4HHB.pdb
python make_configs.py ~/pdb configs.json
node browser_dump.js ../../index.html configs.json browser.json
python parity.py configs.json browser.json     # prints OK/DIFF per configuration
```

`csv_check.js` drives the real page (file upload, Build, Export CSV) and saves the downloaded CSV for comparison with `PREFIX_edges.csv`.

## Limits

- PDB format only (no mmCIF), like the web app. Entries too large for PDB format (over 62 chains or 99,999 atoms in the asymmetric unit) are distributed only as mmCIF and can't be read yet. Capsids like 3VBS are fine, because the PDB file holds one asymmetric unit and the assembly is built from its BIOMT records.
- Memory grows with the number of atom pairs: about 1.2 GB for the 3VBS capsid at residue level, 1.8 GB at atom level (390,000 nodes, 4.8 million edges; 1.5 GB JSON).
