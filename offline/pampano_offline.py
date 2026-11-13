#!/usr/bin/env python3
"""Pampano offline: build Pampano's contact/bond network from a PDB file on the
command line, for structures too large for the browser (capsids, big
bio-assemblies).

It reproduces index.html's network exactly (same parsing, BIOMT handling,
assembly modes, atom selection, contact search, interaction classification,
scope labels, backbone/sidechain partner parts and SASA) and writes:

  <out>.json        Cytoscape.js JSON (.cyjs layout). Opens directly in
                    Cytoscape desktop (File > Import > Network from File) and
                    in Cytoscape.js; plain {nodes, edges} for anything else.
  <out>_edges.csv   one row per edge, same columns as the browser CSV export
                    plus node attributes for both ends (tidy, for R).
  <out>_nodes.csv   one row per node.

Requires Python 3.8+, numpy and scipy.  Run with --help for options, or --gui
for a small form.
"""

import argparse
import csv
import gzip
import json
import math
import os
import re
import sys
import time
import urllib.request

import numpy as np
from scipy.spatial import cKDTree

__version__ = "1.0"

# ---------------------------------------------------------------------------
# property tables (copied from index.html)
# ---------------------------------------------------------------------------
KD = {"ILE": 4.5, "VAL": 4.2, "LEU": 3.8, "PHE": 2.8, "CYS": 2.5, "MET": 1.9, "ALA": 1.8, "GLY": -0.4,
      "THR": -0.7, "SER": -0.8, "TRP": -0.9, "TYR": -1.3, "PRO": -1.6, "HIS": -3.2, "GLU": -3.5, "GLN": -3.5,
      "ASP": -3.5, "ASN": -3.5, "LYS": -3.9, "ARG": -4.5}
VOL = {"GLY": 60.1, "ALA": 88.6, "SER": 89.0, "CYS": 108.5, "ASP": 111.1, "PRO": 112.7, "ASN": 114.1,
       "THR": 116.1, "GLU": 138.4, "VAL": 140.0, "GLN": 143.8, "HIS": 153.2, "MET": 162.9, "ILE": 166.7,
       "LEU": 166.7, "LYS": 168.6, "ARG": 173.4, "PHE": 189.9, "TYR": 193.6, "TRP": 227.8}
MAXASA = {"ALA": 129, "ARG": 274, "ASN": 195, "ASP": 193, "CYS": 167, "GLU": 223, "GLN": 225, "GLY": 104,
          "HIS": 224, "ILE": 197, "LEU": 201, "LYS": 236, "MET": 224, "PHE": 240, "PRO": 159, "SER": 155,
          "THR": 172, "TRP": 285, "TYR": 263, "VAL": 174}
RESCLASS = {"ASP": "neg", "GLU": "neg", "LYS": "pos", "ARG": "pos", "HIS": "pos",
            "SER": "polar", "THR": "polar", "ASN": "polar", "GLN": "polar", "TYR": "polar", "CYS": "polar",
            "ALA": "hydrophobic", "VAL": "hydrophobic", "LEU": "hydrophobic", "ILE": "hydrophobic",
            "MET": "hydrophobic", "PHE": "hydrophobic", "TRP": "hydrophobic", "PRO": "special", "GLY": "special"}
VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80, "H": 1.20, "F": 1.47,
       "CL": 1.75, "BR": 1.85, "I": 1.98, "FE": 2.0, "ZN": 1.39, "MG": 1.73, "CA": 2.31, "NA": 2.27}
THREE_TO_ONE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
                "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
                "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "U", "PYL": "O"}
POS_ATOMS = {"ARG": ("NH1", "NH2", "NE"), "LYS": ("NZ",), "HIS": ("ND1", "NE2")}
NEG_ATOMS = {"ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2")}
HYDROPHOBIC_RES = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "CYS"}
BACKBONE_ATOMS = {"N", "CA", "C", "O", "OXT", "H", "H1", "H2", "H3", "HN", "HA", "HA2", "HA3"}
PARTS = ["backbone", "sidechain", "hetatm"]
PART_ABBR = ["bb", "sc", "het"]
INTERACTION_PRIORITY = ["disulfide", "peptide", "covalent", "salt_bridge", "hydrogen_bond", "hydrophobic", "vdw"]
SCOPE_RANK = ["heterotypic", "homotypic", "inter_protomer", "intra_protomer", "intra"]
COVALENT_SEARCH = 2.4
HBOND_DHA_MIN = 120.0
HBOND_DAAA_MIN = 90.0

# interaction codes = index in INTERACTION_PRIORITY (lower = stronger)
T_DISULFIDE, T_PEPTIDE, T_COVALENT, T_SALT, T_HBOND, T_HYDROPHOBIC, T_VDW = range(7)
# scope codes = index in SCOPE_RANK (lower = more notable)
S_HETERO, S_HOMO, S_INTERPROT, S_INTRAPROT, S_INTRA = range(5)

ASSEMBLY_MODES = ["monomer", "asym", "bioassembly", "biocontacts"]
EDGE_FILTERS = ["all", "bonds", "contacts", "peptide", "disulfide", "salt_bridge", "hydrogen_bond",
                "hydrophobic", "scope_intra", "scope_intra_protomer", "scope_homotypic",
                "scope_heterotypic", "scope_inter_protomer"]


def log(msg):
    print(msg, file=sys.stderr, flush=True)


# ---------------------------------------------------------------------------
# PDB parsing (mirrors parsePDB)
# ---------------------------------------------------------------------------
_INT_RE = re.compile(r"^\s*([+-]?\d+)")
_FLOAT_RE = re.compile(r"^\s*([+-]?(?:\d+\.?\d*(?:[eE][+-]?\d+)?|\.\d+(?:[eE][+-]?\d+)?))")


def js_parse_int(s):
    m = _INT_RE.match(s)
    return int(m.group(1)) if m else None      # None plays the role of NaN


def js_parse_float(s):
    m = _FLOAT_RE.match(s)
    return float(m.group(1)) if m else float("nan")


def parse_pdb(text):
    """Return {'atoms': list of dicts, 'biomt': list of ops} like parsePDB."""
    atoms = []
    biomt = []
    cur_op = None
    mol = None
    apply_chains = None
    re_bm = re.compile(r"BIOMOLECULE:\s*(\d+)")
    re_ap = re.compile(r"(APPLY THE FOLLOWING TO CHAINS|AND CHAINS):(.*)")
    re_bt = re.compile(r"BIOMT([123])\s+\d+\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)")
    for line in text.splitlines():
        rec = line[0:6]
        if rec == "ATOM  " or rec == "HETATM":
            alt = line[16:17].strip()
            if alt not in ("", "A"):
                continue
            el = line[76:78].strip() or re.sub(r"[0-9]", "", line[12:14].strip())
            atoms.append({
                "hetatm": rec == "HETATM",
                "name": line[12:16].strip(),
                "resName": line[17:20].strip(),
                "chain": line[21:22].strip() or "A",
                "resSeq": js_parse_int(line[22:26]),
                "iCode": line[26:27].strip(),
                "x": js_parse_float(line[30:38]),
                "y": js_parse_float(line[38:46]),
                "z": js_parse_float(line[46:54]),
                "element": el.upper(),
            })
        elif line.startswith("REMARK 350"):
            body = line[10:]
            bm = re_bm.search(body)
            if bm:
                mol = int(bm.group(1))
                apply_chains = None
                continue
            if mol is not None and mol != 1:
                continue
            ap = re_ap.search(body)
            if ap:
                lst = [c.strip() for c in ap.group(2).split(",") if c.strip()]
                if ap.group(1) == "AND CHAINS" and apply_chains:
                    apply_chains = apply_chains + lst
                else:
                    apply_chains = lst
                continue
            bt = re_bt.search(body)
            if bt:
                row = int(bt.group(1)) - 1
                nums = []
                for v in bt.groups()[1:5]:
                    try:
                        nums.append(float(v))
                    except ValueError:
                        nums.append(float("nan"))
                if row == 0:
                    cur_op = {"matrix": [[0.0] * 3 for _ in range(3)], "trans": [0.0] * 3,
                              "chains": set(apply_chains) if apply_chains else None}
                if cur_op is not None:
                    cur_op["matrix"][row] = nums[0:3]
                    cur_op["trans"][row] = nums[3]
                    if row == 2:
                        biomt.append(cur_op)
                        cur_op = None
    return {"atoms": atoms, "biomt": biomt}


def read_structure_text(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", errors="replace") as fh:
        return fh.read()


def fetch_pdb(pdb_id):
    url = "https://files.rcsb.org/download/%s.pdb" % pdb_id.lower()
    with urllib.request.urlopen(url, timeout=60) as r:
        return r.read().decode("utf-8", "replace")


# ---------------------------------------------------------------------------
# Atom table: columnar numpy arrays for speed
# ---------------------------------------------------------------------------
class Atoms:
    """Columnar atom table. String columns hold Python lists for output and
    integer codes for vectorised comparisons."""

    def __init__(self, **cols):
        self.__dict__.update(cols)

    def __len__(self):
        return len(self.x)

    def take(self, idx):
        idx = np.asarray(idx, dtype=np.int64)
        out = {}
        for k, v in self.__dict__.items():
            if isinstance(v, np.ndarray):
                out[k] = v[idx]
            elif isinstance(v, list):
                out[k] = [v[i] for i in idx]
            else:
                out[k] = v
        return Atoms(**out)


def atoms_from_parsed(raw):
    n = len(raw)
    res_seq = [a["resSeq"] for a in raw]
    return Atoms(
        x=np.array([a["x"] for a in raw], dtype=np.float64).reshape(n),
        y=np.array([a["y"] for a in raw], dtype=np.float64).reshape(n),
        z=np.array([a["z"] for a in raw], dtype=np.float64).reshape(n),
        hetatm=np.array([a["hetatm"] for a in raw], dtype=bool).reshape(n),
        name=[a["name"] for a in raw],
        resName=[a["resName"] for a in raw],
        chain=[a["chain"] for a in raw],
        baseChain=[a["chain"] for a in raw],
        resSeq=np.array([r if r is not None else 0 for r in res_seq], dtype=np.int64).reshape(n),
        resSeqNaN=np.array([r is None for r in res_seq], dtype=bool).reshape(n),
        iCode=[a["iCode"] for a in raw],
        element=[a["element"] for a in raw],
        opId=np.zeros(n, dtype=np.int64),
        src=np.arange(n, dtype=np.int64),
        mate=np.zeros(n, dtype=bool),
    )


def concat_atoms(parts):
    out = {}
    for k, v in parts[0].__dict__.items():
        if isinstance(v, np.ndarray):
            out[k] = np.concatenate([p.__dict__[k] for p in parts])
        else:
            out[k] = [e for p in parts for e in p.__dict__[k]]
    return Atoms(**out)


def apply_biomt(atoms, biomt):
    if not biomt:
        return atoms
    parts = []
    for oi, op in enumerate(biomt):
        if op["chains"] is not None:
            keep = np.array([c in op["chains"] for c in atoms.chain], dtype=bool)
            idx = np.nonzero(keep)[0]
        else:
            idx = np.arange(len(atoms))
        if len(idx) == 0:
            continue
        p = atoms.take(idx)
        (r0, r1, r2), t = op["matrix"], op["trans"]
        x, y, z = p.x, p.y, p.z
        # same evaluation order as the JS (left to right), so coordinates are bit-identical
        p.x = r0[0] * x + r0[1] * y + r0[2] * z + t[0]
        p.y = r1[0] * x + r1[1] * y + r1[2] * z + t[1]
        p.z = r2[0] * x + r2[1] * y + r2[2] * z + t[2]
        p.opId = np.full(len(idx), oi, dtype=np.int64)
        p.baseChain = list(p.chain)
        if oi != 0:
            p.chain = ["%s_%d" % (c, oi) for c in p.chain]
        parts.append(p)
    return concat_atoms(parts)


def _cell_keys(ix, iy, iz):
    # pack three cell indices into one int64 (cells are small integers)
    return ((ix + 2 ** 20) << 42) | ((iy + 2 ** 20) << 21) | (iz + 2 ** 20)


def symmetry_mate_atoms(atoms, cutoff):
    """Keep each chain's first-operator copy, plus symmetry-mate atoms in a
    grid cell adjacent to a primary atom (mirrors symmetryMateAtoms)."""
    primary_op = {}
    for bc, op in zip(atoms.baseChain, atoms.opId.tolist()):
        if bc not in primary_op:
            primary_op[bc] = op
    is_primary = np.array([primary_op[bc] == op for bc, op in zip(atoms.baseChain, atoms.opId.tolist())],
                          dtype=bool)
    cell = max(cutoff or 0.0, COVALENT_SEARCH)
    ci = np.floor(atoms.x / cell).astype(np.int64)
    cj = np.floor(atoms.y / cell).astype(np.int64)
    ck = np.floor(atoms.z / cell).astype(np.int64)
    pkeys = np.unique(_cell_keys(ci[is_primary], cj[is_primary], ck[is_primary]))
    others = np.nonzero(~is_primary)[0]
    near = np.zeros(len(others), dtype=bool)
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            for dk in (-1, 0, 1):
                k = _cell_keys(ci[others] + di, cj[others] + dj, ck[others] + dk)
                near |= np.isin(k, pkeys)
    keep = np.concatenate([np.nonzero(is_primary)[0], others[near]])
    out = atoms.take(keep)
    out.mate = np.concatenate([np.zeros(int(is_primary.sum()), dtype=bool), np.ones(int(near.sum()), dtype=bool)])
    return out


def assign_residue_copies(atoms):
    seen = {}
    copy = np.zeros(len(atoms), dtype=np.int64)
    for idx, (ch, rs, nan, ic, rn, op, nm) in enumerate(zip(
            atoms.chain, atoms.resSeq.tolist(), atoms.resSeqNaN.tolist(), atoms.iCode,
            atoms.resName, atoms.opId.tolist(), atoms.name)):
        key = (ch, "NaN" if nan else rs, ic, rn, op)
        rec = seen.get(key)
        if rec is None:
            rec = [0, set()]
            seen[key] = rec
        if nm in rec[1]:
            rec[0] += 1
            rec[1] = set()
        rec[1].add(nm)
        copy[idx] = rec[0]
    atoms.copyId = copy
    return atoms


def prepare_atoms(parsed, opts):
    atoms = atoms_from_parsed(parsed["atoms"])
    mode = opts["assemblyMode"]
    if mode == "biocontacts":
        atoms = symmetry_mate_atoms(apply_biomt(atoms, parsed["biomt"]), opts["cutoff"])
    elif mode == "monomer":
        first = next((c for c, h in zip(atoms.chain, atoms.hetatm.tolist()) if not h), None)
        atoms = atoms.take([i for i, c in enumerate(atoms.chain) if c == first])
    elif mode == "bioassembly":
        atoms = apply_biomt(atoms, parsed["biomt"])
    assign_residue_copies(atoms)
    return atoms


def select_atoms(atoms, opts):
    keep = np.ones(len(atoms), dtype=bool)
    if not opts["includeHetatm"]:
        keep &= ~atoms.hetatm
    if opts["atomScope"] == "heavy":
        keep &= np.array([e not in ("H", "D") for e in atoms.element], dtype=bool)
    if opts["atomScope"] == "ca":
        keep &= np.array([n == "CA" or n == "P" or (h and e != "H")
                          for n, h, e in zip(atoms.name, atoms.hetatm.tolist(), atoms.element)], dtype=bool)
    return atoms.take(np.nonzero(keep)[0])


# ---------------------------------------------------------------------------
# chain sequences / clustering
# ---------------------------------------------------------------------------
def chain_sequences(atoms):
    by_chain = {}
    for bc, rs, nan, ic, rn, h in zip(atoms.baseChain, atoms.resSeq.tolist(), atoms.resSeqNaN.tolist(),
                                      atoms.iCode, atoms.resName, atoms.hetatm.tolist()):
        if h:
            continue
        d = by_chain.setdefault(bc, {})
        k = ("NaN" if nan else rs, ic)
        if k not in d:
            d[k] = (rs, ic, rn)
    seqs = {}
    for bc, d in by_chain.items():
        res = sorted(d.values(), key=lambda r: (r[0], r[1] or ""))
        seqs[bc] = "".join(THREE_TO_ONE.get(r[2], "X") for r in res)
    return seqs


def cluster_chains(seqs):
    seq_to_id, cluster = {}, {}
    for ch, s in seqs.items():
        if s not in seq_to_id:
            seq_to_id[s] = len(seq_to_id)
        cluster[ch] = seq_to_id[s]
    return cluster


def parse_protomer_groups(s):
    if not s or not s.strip():
        return None
    groups = {}
    for gi, grp in enumerate(s.split(";")):
        for ch in [c.strip() for c in grp.split(",") if c.strip()]:
            groups[ch] = gi
    return groups or None


# ---------------------------------------------------------------------------
# neighbour searches
# ---------------------------------------------------------------------------
def find_pairs(xyz, cutoff):
    """All atom pairs i<j with d<=cutoff, in the browser's buildContacts order
    (by i, then grid-cell offset, then j) so tie-breaks match exactly."""
    tree = cKDTree(xyz)
    pairs = tree.query_pairs(cutoff * (1 + 1e-9) + 1e-12, output_type="ndarray")
    if len(pairs) == 0:
        return np.zeros(0, np.int64), np.zeros(0, np.int64), np.zeros(0)
    i = pairs[:, 0].astype(np.int64)
    j = pairs[:, 1].astype(np.int64)
    i, j = np.minimum(i, j), np.maximum(i, j)
    del pairs
    dx = xyz[i, 0] - xyz[j, 0]
    dy = xyz[i, 1] - xyz[j, 1]
    dz = xyz[i, 2] - xyz[j, 2]
    d2 = dx * dx + dy * dy + dz * dz
    ok = d2 <= cutoff * cutoff
    i, j, d2 = i[ok], j[ok], d2[ok]
    cells = np.floor(xyz / cutoff).astype(np.int64)
    dc = cells[j] - cells[i] + 1
    off = dc[:, 0] * 9 + dc[:, 1] * 3 + dc[:, 2]
    order = np.lexsort((j, off, i))
    return i[order], j[order], np.sqrt(d2[order])


def build_antecedents(sel, is_h):
    """Bonded heavy atoms (<=1.95 A) and bonded hydrogens (<=1.3 A) of every
    heavy atom, as padded index arrays (-1 = none)."""
    n = len(sel)
    xyz = np.column_stack([sel.x, sel.y, sel.z])
    tree = cKDTree(xyz)
    p = tree.query_pairs(1.95 * (1 + 1e-9), output_type="ndarray")
    if len(p) == 0:
        return np.full((n, 1), -1, np.int64), np.full((n, 1), -1, np.int64)
    a = np.concatenate([p[:, 0], p[:, 1]]).astype(np.int64)
    b = np.concatenate([p[:, 1], p[:, 0]]).astype(np.int64)
    dx = xyz[a, 0] - xyz[b, 0]
    dy = xyz[a, 1] - xyz[b, 1]
    dz = xyz[a, 2] - xyz[b, 2]
    d2 = dx * dx + dy * dy + dz * dz
    heavy_a = ~is_h[a]
    ante_m = heavy_a & ~is_h[b] & (d2 <= 1.95 * 1.95)
    hyd_m = heavy_a & is_h[b] & (d2 <= 1.3 * 1.3)

    def pad(src, dst):
        if len(src) == 0:
            return np.full((n, 1), -1, np.int64)
        order = np.argsort(src, kind="stable")
        src, dst = src[order], dst[order]
        counts = np.bincount(src, minlength=n)
        k = max(1, int(counts.max()))
        out = np.full((n, k), -1, np.int64)
        starts = np.concatenate([[0], np.cumsum(counts)[:-1]])
        pos = np.arange(len(src)) - starts[src]
        out[src, pos] = dst
        return out

    return pad(a[ante_m], b[ante_m]), pad(a[hyd_m], b[hyd_m])


def _angle_ge(V, P, Q, thresh):
    """angleAt(V, P, Q) >= thresh, vectorised (rows)."""
    u = P - V
    w = Q - V
    lu = np.sqrt(u[:, 0] * u[:, 0] + u[:, 1] * u[:, 1] + u[:, 2] * u[:, 2])
    lw = np.sqrt(w[:, 0] * w[:, 0] + w[:, 1] * w[:, 1] + w[:, 2] * w[:, 2])
    denom = lu * lw
    dot = u[:, 0] * w[:, 0] + u[:, 1] * w[:, 1] + u[:, 2] * w[:, 2]
    with np.errstate(invalid="ignore", divide="ignore"):
        c = np.where(denom == 0, 1.0, dot / np.where(denom == 0, 1.0, denom))
    c = np.clip(c, -1.0, 1.0)
    ang = np.where(denom == 0, 0.0, np.arccos(c) * 180.0 / math.pi)
    return ang >= thresh


def hbond_geometry_ok(xyz, ante, hydro, ia, ib):
    """Vectorised hbondGeometryOK for pair arrays ia, ib."""
    m = len(ia)
    if m == 0:
        return np.zeros(0, bool)

    def dha(don, acc):
        hs = hydro[don]                       # (m, k)
        has = (hs >= 0).any(axis=1)
        ok = np.zeros(m, bool)
        for col in range(hs.shape[1]):
            h = hs[:, col]
            v = h >= 0
            if not v.any():
                continue
            r = np.zeros(m, bool)
            idx = np.nonzero(v)[0]
            r[idx] = _angle_ge(xyz[h[idx]], xyz[don[idx]], xyz[acc[idx]], HBOND_DHA_MIN)
            ok |= r
        return has, ok

    has1, ok1 = dha(ia, ib)
    has2, ok2 = dha(ib, ia)
    use_h = has1 | has2
    res_h = (has1 & ok1) | (has2 & ok2)

    def as_acceptor(acc, don):
        an = ante[acc]
        has = (an >= 0).any(axis=1)
        ok = np.zeros(m, bool)
        for col in range(an.shape[1]):
            aa = an[:, col]
            v = aa >= 0
            if not v.any():
                continue
            idx = np.nonzero(v)[0]
            r = np.zeros(m, bool)
            r[idx] = _angle_ge(xyz[acc[idx]], xyz[don[idx]], xyz[aa[idx]], HBOND_DAAA_MIN)
            ok |= r
        return ~has | ok                       # no antecedent: can't test, don't reject

    res_heavy = as_acceptor(ib, ia) | as_acceptor(ia, ib)
    return np.where(use_h, res_h, res_heavy)


def compute_sasa(xyz, elements, probe=1.4, n_sphere=92):
    """Shrake-Rupley SASA per atom, same sampling and neighbour grid as the browser."""
    n = len(xyz)
    if n == 0:
        return np.zeros(0)
    gr = math.pi * (3 - math.sqrt(5))
    pts = []
    for i in range(n_sphere):
        y = 1 - (i / (n_sphere - 1)) * 2
        r = math.sqrt(1 - y * y)
        th = gr * i
        pts.append([math.cos(th) * r, y, math.sin(th) * r])
    pts = np.array(pts)
    rad = np.array([VDW.get(e, 1.70) for e in elements]) + probe
    cell = 2.0 + probe * 2 + 2.0
    cells = np.floor(xyz / cell).astype(np.int64)
    tree = cKDTree(xyz)
    rmax = float(rad.max())
    neigh_lists = tree.query_ball_point(xyz, rad + rmax + 1e-9)
    sasa = np.zeros(n)
    for idx in range(n):
        nb = np.asarray(neigh_lists[idx], dtype=np.int64)
        nb = nb[nb != idx]
        if len(nb):
            # the browser only looks in the 27 grid cells around the atom
            adj = np.all(np.abs(cells[nb] - cells[idx]) <= 1, axis=1)
            nb = nb[adj]
        ri = rad[idx]
        P = xyz[idx] + ri * pts                       # (92, 3)
        if len(nb) == 0:
            sasa[idx] = 4 * math.pi * ri * ri
            continue
        B = xyz[nb]
        dx = P[:, None, 0] - B[None, :, 0]
        dy = P[:, None, 1] - B[None, :, 1]
        dz = P[:, None, 2] - B[None, :, 2]
        buried = (dx * dx + dy * dy + dz * dz < (rad[nb] ** 2)[None, :]).any(axis=1)
        accessible = int(n_sphere - buried.sum())
        sasa[idx] = (accessible / n_sphere) * 4 * math.pi * ri * ri
    return sasa


# ---------------------------------------------------------------------------
# graph construction (mirrors buildGraph)
# ---------------------------------------------------------------------------
def _codes(values):
    """Map a list of hashables to int codes (first-seen order) and the lookup list."""
    lut, out = {}, np.empty(len(values), dtype=np.int64)
    for k, v in enumerate(values):
        c = lut.get(v)
        if c is None:
            c = lut[v] = len(lut)
        out[k] = c
    return out, list(lut)


def _first_per_group(group, n_groups, mask=None):
    """Index of the first row (lowest index) of each group, -1 if none."""
    rows = np.arange(len(group))
    if mask is not None:
        rows = rows[mask]
    first = np.full(n_groups, np.iinfo(np.int64).max, dtype=np.int64)
    np.minimum.at(first, group[rows], rows)
    first[first == np.iinfo(np.int64).max] = -1
    return first


def build_graph(atoms, opts, verbose=True):
    t0 = time.time()
    eff_cutoff = COVALENT_SEARCH if opts["covalentOnly"] else opts["cutoff"]
    sel = select_atoms(atoms, opts)
    n = len(sel)
    xyz = np.column_stack([sel.x, sel.y, sel.z])
    ca_only = opts["atomScope"] == "ca"
    if verbose:
        log("  %d atoms selected (%d symmetry-mate atoms)" % (n, int(sel.mate.sum())))

    pi, pj, pd = find_pairs(xyz, eff_cutoff) if n else (np.zeros(0, np.int64),) * 2 + (np.zeros(0),)
    if verbose:
        log("  %d atom pairs within %.2f A (%.1fs)" % (len(pi), eff_cutoff, time.time() - t0))

    is_h = np.array([e in ("H", "D") for e in sel.element], dtype=bool)
    need_geom = not opts["covalentOnly"] and not ca_only
    if need_geom and n:
        ante, hydro = build_antecedents(sel, is_h)
    else:
        ante = np.full((n, 1), -1, np.int64)
        hydro = np.full((n, 1), -1, np.int64)

    seqs = chain_sequences(atoms)
    cluster = cluster_chains(seqs)
    groups = opts.get("protomerGroups")

    need_sasa = bool(opts.get("sasa"))
    atom_sasa = None
    if need_sasa:
        heavy_idx = np.nonzero(~is_h)[0]
        t1 = time.time()
        s = compute_sasa(xyz[heavy_idx], [sel.element[k] for k in heavy_idx])
        atom_sasa = np.zeros(n)
        atom_sasa[heavy_idx] = s
        if verbose:
            log("  SASA on %d heavy atoms (%.1fs)" % (len(heavy_idx), time.time() - t1))

    has_mates = bool(sel.mate.any())
    home = np.arange(n, dtype=np.int64)
    if has_mates:
        by_src = {}
        for k in np.nonzero(~sel.mate)[0].tolist():
            by_src[int(sel.src[k])] = k
        for k in np.nonzero(sel.mate)[0].tolist():
            home[k] = by_src[int(sel.src[k])]

    # ---- per-atom codes ----
    chain_c, _ = _codes(sel.chain)
    base_c, base_lut = _codes(sel.baseChain)
    clus = np.array([cluster.get(b, -1 - k) for k, b in enumerate(base_lut)], dtype=np.int64)
    grp = np.array([groups[b] if groups and b in groups else -1 for b in base_lut], dtype=np.int64)
    icode_c, _ = _codes(sel.iCode)
    names = sel.name
    resn = sel.resName
    elem = sel.element
    pos_q = np.array([nm in POS_ATOMS.get(r, ()) for nm, r in zip(names, resn)], dtype=bool)
    neg_q = np.array([nm in NEG_ATOMS.get(r, ()) for nm, r in zip(names, resn)], dtype=bool)
    is_c = np.array([e == "C" for e in elem], dtype=bool)
    is_s = np.array([e == "S" for e in elem], dtype=bool)
    is_no = np.array([e in ("N", "O") for e in elem], dtype=bool)
    hyd_res = np.array([r in HYDROPHOBIC_RES for r in resn], dtype=bool)
    nm_C = np.array([nm == "C" for nm in names], dtype=bool)
    nm_N = np.array([nm == "N" for nm in names], dtype=bool)
    nm_CA = np.array([nm == "CA" for nm in names], dtype=bool)
    part = np.array([2 if h else (0 if nm in BACKBONE_ATOMS else 1)
                     for nm, h in zip(names, sel.hetatm.tolist())], dtype=np.int64)
    het = sel.hetatm
    op = sel.opId
    rs = sel.resSeq
    rsnan = sel.resSeqNaN
    mate = sel.mate

    # drop mate-mate pairs
    keep = ~(mate[pi] & mate[pj])
    pi, pj, pd = pi[keep], pj[keep], pd[keep]

    # ---- scope (contactScope) ----
    same_op = op[pi] == op[pj]
    same_base = base_c[pi] == base_c[pj]
    ga, gb = grp[base_c[pi]], grp[base_c[pj]]
    ca_, cb_ = clus[base_c[pi]], clus[base_c[pj]]
    scope = np.where(ca_ == cb_, S_HOMO, S_HETERO)
    scope = np.where(same_op & (ga >= 0) & (ga == gb), S_INTRAPROT, scope)
    scope = np.where(same_base, np.where(same_op, S_INTRA, S_INTERPROT), scope)
    inc = np.ones(len(pi), bool)
    if not opts["intraChain"]:
        inc &= scope != S_INTRA
    if not opts["interChain"]:
        inc &= ~np.isin(scope, [S_INTRAPROT, S_HOMO, S_HETERO])
    if not opts["interProtomer"]:
        inc &= scope != S_INTERPROT
    pi, pj, pd, scope = pi[inc], pj[inc], pd[inc], scope[inc]

    # ---- classification (classifyPair) ----
    m = len(pi)
    typ = np.full(m, T_VDW, dtype=np.int64)
    same_res = (rs[pi] == rs[pj]) & ~rsnan[pi] & ~rsnan[pj] & (icode_c[pi] == icode_c[pj])
    pep_base = ~het[pi] & ~het[pj] & (chain_c[pi] == chain_c[pj]) & ~same_res
    pep = pep_base & (((nm_C[pi] & nm_N[pj]) | (nm_N[pi] & nm_C[pj])) & (pd <= 1.95))
    if ca_only:
        pep |= pep_base & nm_CA[pi] & nm_CA[pj] & ~((nm_C[pi] & nm_N[pj]) | (nm_N[pi] & nm_C[pj])) & \
            (np.abs(rs[pi] - rs[pj]) == 1) & ~rsnan[pi] & ~rsnan[pj] & (pd <= 4.2)
    ss = is_s[pi] & is_s[pj]
    anyh = is_h[pi] | is_h[pj]
    cov = np.where(ss, pd <= 2.4, np.where(anyh, pd <= 1.3, pd <= 1.95))
    bond = pep | cov
    unresolved = ~bond
    # salt bridge
    salt = unresolved & (pd <= 4.0) & ((pos_q[pi] & neg_q[pj]) | (neg_q[pi] & pos_q[pj]))
    typ[salt] = T_SALT
    unresolved &= ~salt
    # hydrogen bond
    cand = unresolved & (pd <= 3.5) & is_no[pi] & is_no[pj]
    ci = np.nonzero(cand)[0]
    if len(ci):
        okh = hbond_geometry_ok(xyz, ante, hydro, pi[ci], pj[ci])
        hb = np.zeros(m, bool)
        hb[ci[okh]] = True
        typ[hb] = T_HBOND
        unresolved &= ~hb
    hph = unresolved & (pd <= 4.0) & is_c[pi] & is_c[pj] & (hyd_res[pi] | hyd_res[pj])
    typ[hph] = T_HYDROPHOBIC
    typ[cov] = np.where(ss[cov], T_DISULFIDE, T_COVALENT)
    typ[pep] = T_PEPTIDE
    if opts["covalentOnly"]:
        pi, pj, pd, scope, typ, bond = pi[bond], pj[bond], pd[bond], scope[bond], typ[bond], bond[bond]
    if verbose:
        log("  %d classified atom contacts (%.1fs)" % (len(pi), time.time() - t0))

    self_sym = 0
    prim_idx = np.nonzero(~mate)[0]

    if opts["granularity"] == "atom":
        node_of = np.full(n, -1, np.int64)
        node_of[prim_idx] = np.arange(len(prim_idx))
        s = node_of[home[pi]]
        t = node_of[home[pj]]
        selfm = s == t
        self_sym = int(selfm.sum())
        k = ~selfm
        s, t, pd2, sc2, ty2, bd2 = s[k], t[k], pd[k], scope[k], typ[k], bond[k]
        pa, pb = part[pi[k]], part[pj[k]]
        lo, hi = np.minimum(s, t), np.maximum(s, t)
        key = lo * len(prim_idx) + hi
        uk, first, inv = np.unique(key, return_index=True, return_inverse=True)
        order = np.argsort(first, kind="stable")           # edges in first-seen order
        rank = np.empty_like(order)
        rank[order] = np.arange(len(order))
        g = rank[inv]
        E = len(uk)
        best_t = np.full(E, 99, np.int64)
        np.minimum.at(best_t, g, ty2)
        is_bt = ty2 == best_t[g]
        mind_bt = np.full(E, np.inf)
        np.minimum.at(mind_bt, g, np.where(is_bt, pd2, np.inf))
        win = _first_per_group(g, E, is_bt & (pd2 == mind_bt[g]))
        # The browser folds duplicates one at a time and compares a same-type
        # newcomer against the running minimum distance of *all* pairs seen so
        # far, so replay that fold exactly for atom pairs seen more than once.
        multi = np.bincount(g, minlength=E) > 1
        if multi.any():
            rows = np.nonzero(multi[g])[0]
            rows = rows[np.argsort(g[rows], kind="stable")]
            cur_g, kept, run_d = -1, -1, 0.0
            for r in rows.tolist():
                gr = int(g[r])
                if gr != cur_g:
                    if cur_g >= 0:
                        win[cur_g] = kept
                    cur_g, kept, run_d = gr, r, float(pd2[r])
                    continue
                tp, tc = ty2[kept], ty2[r]
                if tc < tp or (tc == tp and not run_d <= pd2[r]):
                    kept = r
                run_d = min(run_d, float(pd2[r]))
            win[cur_g] = kept
        dist = np.full(E, np.inf)
        np.minimum.at(dist, g, pd2)
        scp = np.full(E, 99, np.int64)
        np.minimum.at(scp, g, sc2)
        bnd = np.zeros(E, bool)
        np.logical_or.at(bnd, g, bd2)
        links = {
            "source": s[win], "target": t[win], "dist": dist, "count": np.ones(E, np.int64),
            "type": best_t, "scope": scp, "bond": bnd, "sourcePart": pa[win], "targetPart": pb[win],
            "partCounts": None, "partPair": pa[win] * 3 + pb[win],
        }
        nodes = _atom_nodes(sel, prim_idx, atom_sasa)
    else:
        rkey = list(zip(sel.chain, ["NaN" if q else v for q, v in zip(rsnan.tolist(), rs.tolist())],
                        sel.iCode, sel.resName, sel.copyId.tolist()))
        res_of_atom = np.full(n, -1, np.int64)
        res_lut = {}
        res_first = []
        for k in prim_idx.tolist():
            r = res_lut.get(rkey[k])
            if r is None:
                r = res_lut[rkey[k]] = len(res_first)
                res_first.append(k)
            res_of_atom[k] = r
        R = len(res_first)
        res_sasa = np.zeros(R)
        if atom_sasa is not None:
            np.add.at(res_sasa, res_of_atom[prim_idx], atom_sasa[prim_idx])
        ri = res_of_atom[home[pi]]
        rj = res_of_atom[home[pj]]
        same = ri == rj
        self_sym = int((same & ((pi != home[pi]) | (pj != home[pj]))).sum())
        k = ~same
        if opts["excludeSeqAdj"]:
            seqn = ~het[pi] & ~het[pj] & (chain_c[pi] == chain_c[pj]) & \
                (np.abs(rs[pi] - rs[pj]) <= 1) & ~rsnan[pi] & ~rsnan[pj]
            k &= ~seqn
        ri, rj, pd2, sc2, ty2, bd2 = ri[k], rj[k], pd[k], scope[k], typ[k], bond[k]
        pa, pb = part[pi[k]], part[pj[k]]
        lo, hi = np.minimum(ri, rj), np.maximum(ri, rj)
        ps = np.where(ri == lo, pa, pb)
        pt = np.where(ri == lo, pb, pa)
        key = lo * max(R, 1) + hi
        uk, first, inv = np.unique(key, return_index=True, return_inverse=True)
        order = np.argsort(first, kind="stable")
        rank = np.empty_like(order)
        rank[order] = np.arange(len(order))
        g = rank[inv]
        E = len(uk)
        src = (uk // max(R, 1))[order]
        tgt = (uk % max(R, 1))[order]
        cnt = np.bincount(g, minlength=E)
        dist = np.full(E, np.inf)
        np.minimum.at(dist, g, pd2)
        best_t = np.full(E, 99, np.int64)
        np.minimum.at(best_t, g, ty2)
        is_bt = ty2 == best_t[g]
        mind_bt = np.full(E, np.inf)
        np.minimum.at(mind_bt, g, np.where(is_bt, pd2, np.inf))
        win = _first_per_group(g, E, is_bt & (pd2 == mind_bt[g]))
        scp = np.full(E, 99, np.int64)
        np.minimum.at(scp, g, sc2)
        bnd = np.zeros(E, bool)
        np.logical_or.at(bnd, g, bd2)
        pk = ps * 3 + pt
        pc = np.zeros((E, 9), np.int64)
        np.add.at(pc, (g, pk), 1)
        pc_first = np.full((E, 9), np.iinfo(np.int64).max, np.int64)
        np.minimum.at(pc_first, (g, pk), np.arange(len(g)))
        # keep only residues with an edge, renumbered in residue order
        used = np.zeros(R, bool)
        used[src] = True
        used[tgt] = True
        remap = np.cumsum(used) - 1
        links = {
            "source": remap[src], "target": remap[tgt], "dist": dist, "count": cnt,
            "type": best_t, "scope": scp, "bond": bnd, "sourcePart": ps[win], "targetPart": pt[win],
            "partCounts": pc, "partFirst": pc_first,
        }
        nodes = _residue_nodes(sel, [res_first[r] for r in np.nonzero(used)[0]],
                               res_sasa[used] if need_sasa else None,
                               res_of_atom, prim_idx, np.nonzero(used)[0])
    if verbose:
        log("  graph: %d nodes, %d edges (%.1fs)" % (len(nodes["id"]), len(links["source"]), time.time() - t0))
    meta = {"seqs": seqs, "cluster": cluster,
            "mateAtoms": int(mate.sum()) if has_mates else 0, "selfSymPairs": self_sym,
            "selectedAtoms": n}
    return {"nodes": nodes, "links": links, "meta": meta}


def _residue_centroids(sel, res_of_atom, prim_idx, used_res):
    R = int(res_of_atom.max()) + 1 if len(prim_idx) else 0
    sums = np.zeros((R, 3))
    cnt = np.zeros(R)
    r = res_of_atom[prim_idx]
    np.add.at(sums, r, np.column_stack([sel.x[prim_idx], sel.y[prim_idx], sel.z[prim_idx]]))
    np.add.at(cnt, r, 1)
    c = sums / np.maximum(cnt, 1)[:, None]
    return c[used_res]


def _residue_nodes(sel, first_atoms, res_sasa, res_of_atom, prim_idx, used_res):
    cen = _residue_centroids(sel, res_of_atom, prim_idx, used_res)
    nodes = {"id": list(range(len(first_atoms)))}
    ch = [sel.chain[k] for k in first_atoms]
    rn = [sel.resName[k] for k in first_atoms]
    rsq = [None if sel.resSeqNaN[k] else int(sel.resSeq[k]) for k in first_atoms]
    ic = [sel.iCode[k] for k in first_atoms]
    nodes.update({
        "label": ["%s%s" % (r, "NaN" if s is None else s) for r, s in zip(rn, rsq)],
        "chain": ch,
        "baseChain": [sel.baseChain[k] for k in first_atoms],
        "opId": [int(sel.opId[k]) for k in first_atoms],
        "copyId": [int(sel.copyId[k]) for k in first_atoms],
        "hetatm": [bool(sel.hetatm[k]) for k in first_atoms],
        "resName": rn, "resSeq": rsq, "iCode": ic,
        "atomName": [None] * len(first_atoms),
        "hydro": [KD.get(r) for r in rn],
        "vol": [VOL.get(r) for r in rn],
        "sasa": None if res_sasa is None else [float(v) for v in res_sasa],
        "relsasa": None if res_sasa is None else
        [min(v / MAXASA[r], 1) if r in MAXASA else None for v, r in zip(res_sasa.tolist(), rn)],
        "resclass": [RESCLASS.get(r, "unknown") for r in rn],
        "title": ["%s %s %s%s" % (c, r, "NaN" if s is None else s, i) for c, r, s, i in zip(ch, rn, rsq, ic)],
        "xyz": cen,
    })
    return nodes


def _atom_nodes(sel, prim_idx, atom_sasa):
    idx = prim_idx.tolist()
    ch = [sel.chain[k] for k in idx]
    rn = [sel.resName[k] for k in idx]
    rsq = [None if sel.resSeqNaN[k] else int(sel.resSeq[k]) for k in idx]
    ic = [sel.iCode[k] for k in idx]
    nm = [sel.name[k] for k in idx]
    return {
        "id": list(range(len(idx))),
        "label": nm, "chain": ch,
        "baseChain": [sel.baseChain[k] for k in idx],
        "opId": [int(sel.opId[k]) for k in idx],
        "copyId": [int(sel.copyId[k]) for k in idx],
        "hetatm": [bool(sel.hetatm[k]) for k in idx],
        "resName": rn, "resSeq": rsq, "iCode": ic, "atomName": nm,
        "hydro": [KD.get(r) for r in rn],
        "vol": [VOL.get(r) for r in rn],
        "sasa": None if atom_sasa is None else [float(atom_sasa[k]) for k in idx],
        "relsasa": None,
        "resclass": [RESCLASS.get(r, "unknown") for r in rn],
        "title": ["%s %s%s%s %s" % (c, r, "NaN" if s is None else s, i, a)
                  for c, r, s, i, a in zip(ch, rn, rsq, ic, nm)],
        "xyz": np.column_stack([sel.x[prim_idx], sel.y[prim_idx], sel.z[prim_idx]]),
    }


# ---------------------------------------------------------------------------
# edge filter / degree / part counts
# ---------------------------------------------------------------------------
def filter_mask(links, ef):
    E = len(links["source"])
    if ef == "all":
        return np.ones(E, bool)
    if ef == "bonds":
        return links["bond"].copy()
    if ef == "contacts":
        return ~links["bond"]
    if ef.startswith("scope_"):
        return links["scope"] == SCOPE_RANK.index(ef[6:])
    return links["type"] == INTERACTION_PRIORITY.index(ef)


def part_counts_list(links, e):
    """[(key, count)] in the browser object's insertion order."""
    pc = links.get("partCounts")
    if pc is None:
        p = int(links["partPair"][e])
        return [("%s-%s" % (PART_ABBR[p // 3], PART_ABBR[p % 3]), 1)]
    row, first = pc[e], links["partFirst"][e]
    ks = [k for k in np.argsort(first, kind="stable") if row[k] > 0]
    return [("%s-%s" % (PART_ABBR[k // 3], PART_ABBR[k % 3]), int(row[k])) for k in ks]


def format_part_counts(lst):
    # largest first; ties keep insertion order (JS sort is stable)
    return ";".join("%s:%d" % kv for kv in sorted(lst, key=lambda kv: -kv[1]))


# ---------------------------------------------------------------------------
# writers
# ---------------------------------------------------------------------------
def _num(v):
    if v is None:
        return None
    if isinstance(v, (np.floating, float)):
        f = float(v)
        return None if math.isnan(f) else f
    if isinstance(v, (np.integer,)):
        return int(v)
    return v


def _edge_part_strings(links, eidx):
    """part_counts strings ("sc-sc:12;bb-sc:3") and per-combination count matrix for edges eidx."""
    if links.get("partCounts") is None:
        pp = links["partPair"][eidx]
        lut = np.array(["%s-%s:1" % (PART_ABBR[k // 3], PART_ABBR[k % 3]) for k in range(9)], dtype=object)
        mat = np.zeros((len(eidx), 9), np.int64)
        mat[np.arange(len(eidx)), pp] = 1
        return lut[pp], mat
    strs = np.empty(len(eidx), dtype=object)
    for n, e in enumerate(eidx.tolist()):
        strs[n] = format_part_counts(part_counts_list(links, e))
    return strs, links["partCounts"][eidx]


def write_outputs(g, opts, out_prefix, source_name, edge_filter="all", position_scale=10.0):
    """Write PREFIX.json (Cytoscape.js), PREFIX_edges.csv and PREFIX_nodes.csv, streaming rows."""
    nodes, links = g["nodes"], g["links"]
    atom_mode = opts["granularity"] == "atom"
    mask = filter_mask(links, edge_filter)
    eidx = np.nonzero(mask)[0]
    N = len(nodes["id"])
    src = links["source"][eidx]
    tgt = links["target"][eidx]
    deg = np.bincount(src, minlength=N) + np.bincount(tgt, minlength=N)
    tnames = np.array(INTERACTION_PRIORITY, dtype=object)[links["type"][eidx]]
    snames = np.array(SCOPE_RANK, dtype=object)[links["scope"][eidx]]
    pnames = np.array(PARTS, dtype=object)
    spart = pnames[links["sourcePart"][eidx]]
    tpart = pnames[links["targetPart"][eidx]]
    dist = links["dist"][eidx]
    cnt = links["count"][eidx]
    bond = links["bond"][eidx]
    pc_str, pc_mat = _edge_part_strings(links, eidx)
    pc_keys = ["%s_%s" % (a, b) for a in PART_ABBR for b in PART_ABBR]
    titles = nodes["title"]
    xyz = nodes["xyz"]

    node_keys = ["chain", "base_chain", "operator", "res_name", "res_seq", "icode"] + \
        (["atom_name"] if atom_mode else []) + ["hetatm", "res_class"]

    def node_row(k):
        return {
            "id": str(k), "name": titles[k], "label": nodes["label"][k],
            "chain": nodes["chain"][k], "base_chain": nodes["baseChain"][k],
            "operator": nodes["opId"][k], "copy": nodes["copyId"][k],
            "res_name": nodes["resName"][k], "res_seq": nodes["resSeq"][k], "icode": nodes["iCode"][k],
            "atom_name": nodes["atomName"][k], "hetatm": nodes["hetatm"][k],
            "res_class": nodes["resclass"][k], "hydropathy": nodes["hydro"][k], "volume": nodes["vol"][k],
            "sasa": None if nodes["sasa"] is None else nodes["sasa"][k],
            "rel_sasa": None if nodes["relsasa"] is None else nodes["relsasa"][k],
            "degree": int(deg[k]),
            "x": float(xyz[k, 0]), "y": float(xyz[k, 1]), "z": float(xyz[k, 2]),
        }

    def clean(d):
        out = {}
        for k, v in d.items():
            v = _num(v)
            if v is not None:
                out[k] = v
        return out

    # ---- nodes CSV + collect node JSON lines ----
    nodes_path = out_prefix + "_nodes.csv"
    ncols = ["id", "name", "label", "chain", "base_chain", "operator", "copy", "res_name", "res_seq", "icode"] + \
        (["atom_name"] if atom_mode else []) + \
        ["hetatm", "res_class", "hydropathy", "volume", "sasa", "rel_sasa", "degree", "x", "y", "z"]
    node_attr_cache = []
    json_path = out_prefix + ".json"
    header = {
        "format_version": "1.0",
        "generated_by": "pampano_offline %s" % __version__,
        "target_cytoscapejs_version": "~2.1",
        "data": {
            "name": "Pampano %s" % source_name, "shared_name": "Pampano %s" % source_name,
            "source": source_name,
            "assembly_mode": opts["assemblyMode"], "granularity": opts["granularity"],
            "cutoff": opts["cutoff"], "atom_scope": opts["atomScope"],
            "covalent_only": opts["covalentOnly"], "intra_chain": opts["intraChain"],
            "inter_chain": opts["interChain"], "inter_protomer": opts["interProtomer"],
            "exclude_seq_adjacent": opts["excludeSeqAdj"], "include_hetatm": opts["includeHetatm"],
            "protomer_groups": opts.get("protomerGroupsText") or "",
            "sasa": bool(opts.get("sasa")), "edge_filter": edge_filter,
            "symmetry_mate_atoms": g["meta"]["mateAtoms"],
            "self_symmetry_atom_pairs": g["meta"]["selfSymPairs"],
        },
    }
    with open(json_path, "w") as jf, open(nodes_path, "w", newline="") as nf:
        nw = csv.writer(nf)
        nw.writerow(ncols)
        jf.write(json.dumps(header, separators=(",", ":"))[:-1])
        jf.write(',"elements":{"nodes":[')
        for k in range(N):
            r = node_row(k)
            nw.writerow([_csv(r[c]) for c in ncols])
            node_attr_cache.append([_csv(r[c]) for c in node_keys])
            if k:
                jf.write(",")
            jf.write(json.dumps({"data": clean(r), "position": {"x": r["x"] * position_scale,
                                                                 "y": r["y"] * position_scale}},
                                separators=(",", ":")))
        jf.write('],"edges":[')
        # ---- edges: JSON + CSV in one pass ----
        edges_path = out_prefix + "_edges.csv"
        ecols = ["source", "target", "source_label", "target_label", "distance", "count", "interaction", "scope",
                 "is_bond", "source_part", "target_part", "part_counts"]
        with open(edges_path, "w", newline="") as ef:
            ew = csv.writer(ef)
            ew.writerow(ecols + ["pc_" + k for k in pc_keys] +
                        ["source_" + c for c in node_keys] + ["target_" + c for c in node_keys])
            dumps = json.JSONEncoder(separators=(",", ":")).encode
            for n in range(len(eidx)):
                s, t = int(src[n]), int(tgt[n])
                typ = tnames[n]
                d = {"id": "e%d" % n, "source": str(s), "target": str(t),
                     "name": "%s (%s) %s" % (titles[s], typ, titles[t]),
                     "interaction": typ, "scope": snames[n], "distance": float(dist[n]), "count": int(cnt[n]),
                     "is_bond": bool(bond[n]), "source_part": spart[n], "target_part": tpart[n],
                     "part_counts": pc_str[n]}
                row = pc_mat[n]
                for q in range(9):
                    if row[q]:
                        d["pc_" + pc_keys[q]] = int(row[q])
                if n:
                    jf.write(",")
                jf.write('{"data":')
                jf.write(dumps(d))
                jf.write("}")
                ew.writerow([s, t, titles[s], titles[t], "%.2f" % dist[n], int(cnt[n]), typ, snames[n],
                             1 if bond[n] else 0, spart[n], tpart[n], pc_str[n]] + row.tolist() +
                            node_attr_cache[s] + node_attr_cache[t])
        jf.write("]}}")
    return json_path, edges_path, nodes_path


def _csv(v):
    if v is None:
        return "NA"
    if isinstance(v, bool):
        return "TRUE" if v else "FALSE"
    if isinstance(v, float):
        return "NA" if math.isnan(v) else ("%.4f" % v).rstrip("0").rstrip(".") if v != int(v) else "%g" % v
    return v


def browser_style_graph(g):
    """The graph as the browser's buildGraph returns it (for parity tests)."""
    nodes, links = g["nodes"], g["links"]
    out_nodes = []
    for k in range(len(nodes["id"])):
        out_nodes.append({
            "id": k, "label": nodes["label"][k], "chain": nodes["chain"][k], "baseChain": nodes["baseChain"][k],
            "opId": nodes["opId"][k], "copyId": nodes["copyId"][k], "hetatm": nodes["hetatm"][k],
            "resName": nodes["resName"][k], "resSeq": nodes["resSeq"][k], "iCode": nodes["iCode"][k],
            "hydro": nodes["hydro"][k], "vol": nodes["vol"][k],
            "sasa": None if nodes["sasa"] is None else nodes["sasa"][k],
            "relsasa": None if nodes["relsasa"] is None else nodes["relsasa"][k],
            "resclass": nodes["resclass"][k], "title": nodes["title"][k],
        })
    out_links = []
    for e in range(len(links["source"])):
        out_links.append({
            "source": int(links["source"][e]), "target": int(links["target"][e]),
            "dist": float(links["dist"][e]), "count": int(links["count"][e]),
            "type": INTERACTION_PRIORITY[int(links["type"][e])], "scope": SCOPE_RANK[int(links["scope"][e])],
            "bond": bool(links["bond"][e]),
            "sourcePart": PARTS[int(links["sourcePart"][e])], "targetPart": PARTS[int(links["targetPart"][e])],
            "partCounts": dict(part_counts_list(links, e)),
        })
    return {"nodes": out_nodes, "links": out_links,
            "meta": {"mateAtoms": g["meta"]["mateAtoms"], "selfSymPairs": g["meta"]["selfSymPairs"]}}


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------
def default_opts():
    return {"assemblyMode": "asym", "granularity": "residue", "cutoff": 5.0, "atomScope": "heavy",
            "covalentOnly": False, "intraChain": True, "interChain": True, "interProtomer": True,
            "protomerGroups": None, "protomerGroupsText": "", "excludeSeqAdj": False,
            "includeHetatm": False, "sasa": False}


def run(input_path=None, pdb_id=None, out_prefix=None, opts=None, edge_filter="all", verbose=True):
    opts = dict(default_opts(), **(opts or {}))
    t0 = time.time()
    if pdb_id:
        name = pdb_id.upper()
        if verbose:
            log("Fetching %s from RCSB" % name)
        text = fetch_pdb(pdb_id)
    else:
        name = os.path.basename(input_path)
        text = read_structure_text(input_path)
    parsed = parse_pdb(text)
    if verbose:
        log("Parsed %s: %d atoms, %d BIOMT operators" % (name, len(parsed["atoms"]), len(parsed["biomt"])))
    if opts["assemblyMode"] in ("bioassembly", "biocontacts") and not parsed["biomt"] and verbose:
        log("  No BIOMT records; using the asymmetric unit.")
    atoms = prepare_atoms(parsed, opts)
    if verbose:
        log("  %d atoms after assembly (%s)" % (len(atoms), opts["assemblyMode"]))
    g = build_graph(atoms, opts, verbose=verbose)
    if out_prefix is None:
        stem = re.sub(r"(\.gz)?$", "", name)
        stem = os.path.splitext(stem)[0]
        out_prefix = "%s_%s_%s" % (stem, opts["assemblyMode"], opts["granularity"])
    paths = write_outputs(g, opts, out_prefix, name, edge_filter)
    if verbose:
        if g["meta"]["mateAtoms"]:
            log("  %d symmetry-mate atoms; %d atom pairs between a residue and its own symmetry copy (no edge)"
                % (g["meta"]["mateAtoms"], g["meta"]["selfSymPairs"]))
        log("Wrote %s (%.1fs total)" % (", ".join(paths), time.time() - t0))
    return g, paths


def build_parser():
    p = argparse.ArgumentParser(
        prog="pampano_offline.py",
        description="Build Pampano's residue/atom contact network from a PDB file and write Cytoscape JSON "
                    "plus CSV edge and node tables. Same network as the browser app; no size limit.")
    src = p.add_mutually_exclusive_group()
    src.add_argument("input", nargs="?", help="PDB file (.pdb, .ent, optionally .gz)")
    src.add_argument("--pdb-id", help="fetch this 4-character entry from RCSB instead of reading a file")
    p.add_argument("-o", "--out", help="output prefix (default: <file>_<assembly>_<granularity>); writes "
                                       "PREFIX.json, PREFIX_edges.csv, PREFIX_nodes.csv")
    p.add_argument("-a", "--assembly", choices=ASSEMBLY_MODES, default="asym",
                   help="monomer (first chain), asym (all chains as deposited, default), bioassembly "
                        "(apply BIOMT, every copy gets nodes), biocontacts (one copy, symmetry-mate "
                        "contacts mapped back)")
    p.add_argument("-g", "--granularity", choices=["residue", "atom"], default="residue")
    p.add_argument("-c", "--cutoff", type=float, default=5.0, help="contact cutoff in A (default 5.0)")
    p.add_argument("--atoms", choices=["all", "heavy", "ca"], default="heavy", help="atom scope (default heavy)")
    p.add_argument("--covalent-only", action="store_true", help="covalent bonds only")
    p.add_argument("--no-intra-chain", action="store_true", help="drop intra-chain contacts")
    p.add_argument("--no-inter-chain", action="store_true", help="drop inter-chain contacts")
    p.add_argument("--no-inter-protomer", action="store_true", help="drop inter-protomer (A-A') contacts")
    p.add_argument("--protomer-groups", default="", help='manual protomer groups, e.g. "A,B ; C,D"')
    p.add_argument("--exclude-seq-adjacent", action="store_true",
                   help="drop contacts between consecutive residues (residue mode)")
    p.add_argument("--hetatm", action="store_true", help="include HETATM (ligands, ions, water)")
    p.add_argument("--sasa", action="store_true", help="compute SASA (slow on very large selections)")
    p.add_argument("--edge-filter", choices=EDGE_FILTERS, default="all", metavar="FILTER",
                   help="write only these edges: " + ", ".join(EDGE_FILTERS) +
                        " (default all); degree is counted on written edges")
    p.add_argument("--gui", action="store_true", help="open a small form instead of using the command line")
    p.add_argument("-q", "--quiet", action="store_true")
    p.add_argument("--version", action="version", version=__version__)
    return p


def opts_from_args(a):
    return {
        "assemblyMode": a.assembly, "granularity": a.granularity, "cutoff": a.cutoff, "atomScope": a.atoms,
        "covalentOnly": a.covalent_only, "intraChain": not a.no_intra_chain,
        "interChain": not a.no_inter_chain, "interProtomer": not a.no_inter_protomer,
        "protomerGroups": parse_protomer_groups(a.protomer_groups), "protomerGroupsText": a.protomer_groups,
        "excludeSeqAdj": a.exclude_seq_adjacent, "includeHetatm": a.hetatm, "sasa": a.sasa,
    }


def main(argv=None):
    p = build_parser()
    a = p.parse_args(argv)
    if a.gui:
        return run_gui()
    if not a.input and not a.pdb_id:
        p.error("give a PDB file or --pdb-id (or --gui)")
    run(a.input, a.pdb_id, a.out, opts_from_args(a), a.edge_filter, verbose=not a.quiet)
    return 0


# ---------------------------------------------------------------------------
# minimal Tk GUI
# ---------------------------------------------------------------------------
def run_gui():
    import threading
    import tkinter as tk
    from tkinter import filedialog, messagebox, ttk

    root = tk.Tk()
    root.title("Pampano offline")
    v = {
        "input": tk.StringVar(), "out": tk.StringVar(), "assembly": tk.StringVar(value="biocontacts"),
        "granularity": tk.StringVar(value="residue"), "cutoff": tk.StringVar(value="5.0"),
        "atoms": tk.StringVar(value="heavy"), "groups": tk.StringVar(), "edge_filter": tk.StringVar(value="all"),
    }
    flags = {k: tk.BooleanVar(value=d) for k, d in [
        ("covalentOnly", False), ("intraChain", True), ("interChain", True), ("interProtomer", True),
        ("excludeSeqAdj", False), ("includeHetatm", False), ("sasa", False)]}
    frm = ttk.Frame(root, padding=12)
    frm.grid(sticky="nsew")
    row = 0

    def browse_in():
        f = filedialog.askopenfilename(filetypes=[("PDB", "*.pdb *.ent *.pdb.gz *.ent.gz *.txt"), ("All", "*")])
        if f:
            v["input"].set(f)
            if not v["out"].get():
                stem = re.sub(r"(\.gz)?$", "", f)
                v["out"].set(os.path.splitext(stem)[0] + "_pampano")

    def browse_out():
        f = filedialog.asksaveasfilename(defaultextension="")
        if f:
            v["out"].set(re.sub(r"\.json$", "", f))

    for label, key, btn in [("PDB file", "input", browse_in), ("Output prefix", "out", browse_out)]:
        ttk.Label(frm, text=label).grid(row=row, column=0, sticky="w")
        ttk.Entry(frm, textvariable=v[key], width=50).grid(row=row, column=1, sticky="we")
        ttk.Button(frm, text="Browse…", command=btn).grid(row=row, column=2)
        row += 1
    for label, key, values in [("Assembly", "assembly", ASSEMBLY_MODES), ("Granularity", "granularity",
                                                                          ["residue", "atom"]),
                               ("Atoms", "atoms", ["all", "heavy", "ca"]), ("Edges", "edge_filter", EDGE_FILTERS)]:
        ttk.Label(frm, text=label).grid(row=row, column=0, sticky="w")
        ttk.Combobox(frm, textvariable=v[key], values=values, state="readonly").grid(row=row, column=1, sticky="w")
        row += 1
    for label, key in [("Cutoff (Å)", "cutoff"), ("Protomer groups", "groups")]:
        ttk.Label(frm, text=label).grid(row=row, column=0, sticky="w")
        ttk.Entry(frm, textvariable=v[key], width=20).grid(row=row, column=1, sticky="w")
        row += 1
    for key, label in [("covalentOnly", "Covalent bonds only"), ("intraChain", "Intra-chain contacts"),
                       ("interChain", "Inter-chain contacts"), ("interProtomer", "Inter-protomer contacts"),
                       ("excludeSeqAdj", "Exclude sequential neighbours"), ("includeHetatm", "Include HETATM"),
                       ("sasa", "Compute SASA (slow on huge assemblies)")]:
        ttk.Checkbutton(frm, text=label, variable=flags[key]).grid(row=row, column=1, sticky="w")
        row += 1
    status = tk.StringVar(value="Pick a PDB file and press Build.")
    ttk.Label(frm, textvariable=status, wraplength=480).grid(row=row + 1, column=0, columnspan=3, sticky="w")
    btn = ttk.Button(frm, text="Build network")
    btn.grid(row=row, column=1, sticky="w", pady=8)

    def work():
        try:
            opts = {k: b.get() for k, b in flags.items()}
            opts.update(assemblyMode=v["assembly"].get(), granularity=v["granularity"].get(),
                        cutoff=float(v["cutoff"].get()), atomScope=v["atoms"].get(),
                        protomerGroups=parse_protomer_groups(v["groups"].get()),
                        protomerGroupsText=v["groups"].get())
            g, paths = run(v["input"].get(), None, v["out"].get() or None, opts, v["edge_filter"].get())
            msg = "Done: %d nodes, %d edges.\n%s" % (len(g["nodes"]["id"]), len(g["links"]["source"]),
                                                     "\n".join(paths))
            root.after(0, lambda: (status.set(msg), btn.state(["!disabled"])))
        except Exception as exc:  # show any failure in the window
            err = str(exc)
            root.after(0, lambda: (status.set("Failed: " + err), btn.state(["!disabled"]),
                                   messagebox.showerror("Pampano", err)))

    def start():
        if not v["input"].get():
            messagebox.showwarning("Pampano", "Choose a PDB file first.")
            return
        btn.state(["disabled"])
        status.set("Building… (large assemblies can take a few minutes)")
        threading.Thread(target=work, daemon=True).start()

    btn.configure(command=start)
    root.mainloop()
    return 0


if __name__ == "__main__":
    sys.exit(main())
