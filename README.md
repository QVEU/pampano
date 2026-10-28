# Pampano

**Pampano** generates a network of contacts and bonds from a protein structure (PDB file). Each node is a residue (or atom), and each edge is a physical interaction between them — a covalent bond, a hydrogen bond, a salt bridge, a hydrophobic contact, and so on. It runs entirely in your browser; nothing is uploaded to a server.

---

## Quick start

1. **Load a structure.** Either upload a `.pdb` file, or type a 4-character PDB ID (e.g. `1CRN`, `4HHB`, `1HHO`) and click **Fetch** to pull it directly from RCSB.
2. **Choose an assembly** (monomer, asymmetric unit, or biological assembly).
3. Click **Build Network.**
4. Explore: drag nodes, zoom/pan, change how nodes and edges are colored, filter edges, and **click any node** to open a 3D structural view of that residue and its contacts.
5. **Export** the result as a CSV edge list or a JSON graph for use in Cytoscape, Gephi, NetworkX, etc.

---

## Controls

### Loading & assembly

- **PDB file / Fetch by ID** — Provide the structure either way. Fetch requires an internet connection and a valid RCSB entry.
- **Assembly**
  - *Monomer* — the first protein chain only.
  - *Full asymmetric unit* — every chain exactly as deposited.
  - *Biological assembly* — applies the `REMARK 350` BIOMT symmetry operators of biological assembly 1 (each to the chains it lists) to rebuild the biologically relevant multimer. If the file has no BIOMT records, this falls back to the asymmetric unit.

### Network definition

- **Granularity**
  - *Residue-level* — one node per residue; an edge aggregates all atom-atom contacts between two residues (storing the closest distance and the count of contributing atom pairs).
  - *Atom-level* — one node per atom; edges are individual atom-atom contacts. More detailed and heavier for large structures.
- **Contact cutoff (Å)** — Maximum distance for two atoms to count as being in contact (default 5.0 Å).
- **Atoms** — Which atoms to consider: all atoms, heavy atoms only (no hydrogens), or CA/representative only.

### Visualization

- **Color nodes by** — chain, hydrophobicity, surface accessibility (SASA), residue type (charge/polarity), or degree (number of contacts).
- **Size nodes by** — uniform, residue volume, SASA, or degree.
- **Show edges** — filter to a single category (all, bonds only, contacts only, or one specific interaction/scope type).
- **Color edges by** — interaction type, contact scope, or a simple intra/inter split.

### Advanced

- **Covalent bonds only** — build a network of covalent linkages alone (≤ 1.95 Å; S–S ≤ 2.4 Å; X–H ≤ 1.3 Å).
- **Include intra-chain / inter-chain / inter-protomer contacts** — toggle each structural scope on or off.
- **Exclude sequential neighbors** — when on, removes contacts between consecutive polymer residues (HETATM groups are unaffected) so you see only tertiary/quaternary contacts. Off by default (the backbone chain is included). Covalent backbone bonds are never removed by this option when off.
- **Protomer groups** — optionally group chains into protomers (e.g. `A,B ; C,D`) to distinguish within-protomer from between-protomer interfaces. Leave blank to classify inter-chain contacts automatically by sequence.
- **Include HETATM** — include ligands, cofactors, ions, and water.

### 3D inset

Click any node to open a structural view (bottom-left) of that residue and the residues it contacts. The selected residue is shown with full side-chain atoms and highlighted; context and sequential residues are shown as backbone only, for orientation. Drag to rotate, scroll to zoom. Clicking the background or the × closes it.

---

## Exports

- **CSV edge list** — one row per edge: node indices and labels, minimum distance, atom-pair count, interaction type, scope, whether it is a covalent bond, `source_part` / `target_part` (backbone, sidechain or hetatm for each partner) and `part_counts`.
- **JSON graph** — nodes (with per-node properties: residue name, hydropathy, volume, SASA, degree, chain, base chain, operator) and links (with distance, count, interaction, scope, bond flag, source/target part, part counts).

---

## Glossary

### Structural scope (where a contact sits)

- **Intra-chain** — Both atoms belong to the same physical chain copy (same chain ID, same symmetry operator). Internal to one subunit.
- **Inter-chain** — The two atoms come from different deposited chains (e.g. chain A and chain B). Automatically subdivided by sequence into:
  - **Homotypic** — the two chains have the *same* sequence (the same protein — typical of a homo-oligomer protomer interface).
  - **Heterotypic** — the two chains have *different* sequences (different proteins — a true hetero-complex interface).
- **Intra-protomer** — The two chains belong to the *same* manually defined protomer group (set via **Protomer groups**). Represents an interface internal to a multi-chain protomer.
- **Inter-protomer** — The two atoms are the same original chain under *different* BIOMT symmetry operators (e.g. A and A′). These interfaces exist only once the biological assembly is built. Only appears in biological-assembly mode.

> **Why both homotypic/heterotypic and intra/inter-protomer?** Chain labels (A, B, C…) are assigned by whoever deposited the structure and don't inherently tell you which chains form one functional unit. Pampano classifies inter-chain contacts automatically by sequence identity (which is unambiguous within a single structure), and lets you override that with explicit protomer groups when the biological grouping is something only you know.

### Interaction types (what kind of contact)

Listed from strongest/most-specific to weakest. When several atom-pair interactions contribute to one residue edge, the strongest one labels the edge.

- **Peptide bond** — The covalent backbone linkage between consecutive residues, C(i)–N(i+1). Detected geometrically from the actual backbone atoms, not from residue numbering.
- **Disulfide** — A covalent S–S bond between two cysteines (~2.05 Å).
- **Covalent (other)** — Any other covalent linkage (≤ ~1.95 Å; ≤ 1.3 Å when a hydrogen is involved), such as ligand attachments.
- **Salt bridge** — An electrostatic interaction between oppositely charged groups: a basic nitrogen (Arg, Lys, His) and an acidic oxygen (Asp, Glu) within 4.0 Å.
- **Hydrogen bond** — A donor N/O and an acceptor N/O within 3.5 Å that also satisfy an angular geometry check (see below).
- **Hydrophobic** — A carbon–carbon contact within 4.0 Å involving a hydrophobic residue.
- **Van der Waals (vdW)** — A close contact within the cutoff that doesn't meet any of the more specific criteria above.

### Hydrogen-bond geometry

Pampano doesn't accept an H-bond on distance alone; it also checks the angle.

- **With explicit hydrogens** (all-atom scope on a structure that has them): the true donor–H···acceptor angle must be ≥ 120° (near-linear).
- **Without hydrogens** (the usual case for crystal structures): it uses the heavy-atom surrogate — the angle at the acceptor between the donor and the acceptor's covalently bonded antecedent (D···A–AA) must be ≥ 90°, so the donor approaches from the lone-pair side rather than through a bonded atom.

### Interaction partners (backbone / side chain / HETATM)
Each edge records which part of each partner is involved: **backbone** (N, CA, C, O, OXT and their hydrogens), **sidechain** (any other atom of a polymer residue), or **hetatm** (ligand, ion, water).
- At atom level this is simply the part of each of the two atoms.
- At residue level an edge combines many atom pairs. `source_part` / `target_part` give the parts of the atom pair that defines the edge's interaction type (the closest pair of the strongest type) — so a hydrogen-bond edge with `backbone` / `sidechain` means that H-bond runs from the source residue's backbone to the target residue's side chain. `part_counts` tallies every contributing atom pair, e.g. `sc-sc:12;bb-sc:3` (bb = backbone, sc = side chain, het = HETATM; the first part belongs to the source node).

### Node properties

- **Hydrophobicity** — Kyte–Doolittle hydropathy; higher = more hydrophobic.
- **Residue volume** — Approximate side-chain/residue volume (Zamyatnin values).
- **SASA (solvent-accessible surface area)** — How exposed a residue is to solvent, computed with a Shrake–Rupley approximation (1.4 Å water probe, Fibonacci sphere sampling). At residue level it's also expressed as a *relative* value — the fraction of the maximum possible SASA for that residue type (Tien et al. reference values), i.e. buried (0) to fully exposed (1).
- **Degree** — The number of edges connected to a node; high-degree nodes are network hubs.
- **Base chain / operator** — The original chain ID and the BIOMT operator index that produced a given atom copy (relevant in biological-assembly mode).

### Other terms

- **BIOMT** — The transformation matrices in a PDB file's `REMARK 350` records that specify how to generate the biological assembly from the deposited coordinates.
- **Asymmetric unit** — The set of coordinates actually deposited in the file; may be a fraction of, or a multiple of, the biological unit.
- **Protomer** — The repeating structural unit of an oligomer. May be a single chain or several chains together.
- **HETATM** — Non-standard residues in a PDB file: ligands, cofactors, ions, modified residues, and water.

---

## Caveats & notes

- **SASA is computed on the current selection.** A residue buried in the full complex will look artificially exposed if you compute SASA on an isolated monomer. Build the biological assembly for meaningful burial values. SASA is also the most computationally intensive step on large all-atom selections.
- **Hydrogen bonds are approximate.** Without explicit hydrogens the geometry check is a heavy-atom surrogate; His is treated as potentially positive regardless of protonation state. This is standard for heavy-atom contact networks but will differ somewhat from dedicated tools (DSSP, HBPLUS).
- **CA-only mode cannot resolve bond chemistry.** With only alpha carbons present, true covalent geometry is unavailable; the backbone is approximated from CA(i)–CA(i+1) adjacency.
- **Biological assemblies multiply atom count** by the number of symmetry operators, so large assemblies are best viewed at residue-level granularity.
- **Sequence classification uses exact identity.** Within one structure, copies of the same protein share an identical sequence; any difference marks a different protein. There's no fuzzy matching, so it cannot misclassify.

---

*Pampano runs locally in your browser. Structures fetched by ID are retrieved from RCSB; uploaded files never leave your machine.*
