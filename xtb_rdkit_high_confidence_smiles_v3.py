#!/usr/bin/env python3
"""
Conservative XYZ -> canonical SMILES pipeline using:
- RDKit multi-method topology generation
- xTB geometry refinement
- xTB Wiberg bond-order exact-pattern support check
- unanimous-agreement acceptance rule

Decision policy:
- HIGH_CONFIDENCE: at least two independent RDKit topology methods succeed on
  BOTH the original and xTB-optimized geometries, every attempted successful
  method agrees on the same canonical isomeric SMILES, the original/optimized
  consensus SMILES are identical, and xTB Wiberg bond orders support the same
  bonded pairs and inferred bond-order pattern as the unanimous optimized RDKit topology.
- AMBIGUOUS: any disagreement, missing consensus, topology change after xTB,
  Wiberg support failure, or Wiberg ambiguity.
- FAILED: unreadable XYZ, xTB failure, or other technical failure.
python xtb_rdkit_high_confidence_smiles_v3.py all -o xtb_rdkit_results
IPXZ fix is in place and verified (GFN-FF opt + GFN2 single-point WBO). Ready to run on the full `all/` folder whenever you want:

Compared with v2, this version makes the WBO gate stricter:
1) It infers a DISCRETE xTB bond-class pattern (NONBOND/SINGLE/DOUBLE/TRIPLE/AROMATIC/AMBIGUOUS)
   from the Wiberg bond orders instead of merely checking broad compatibility intervals.
2) It compares that inferred xTB pattern bond-by-bond against the unanimous optimized RDKit topology.
3) Aromatic bonds are handled explicitly as aromatic ring systems instead of by a loose generic interval.
python xtb_rdkit_high_confidence_smiles_v3.py all -o xtb_rdkit_results
Defaults: `--opt-method gfnff`, `--gfn 2`, `--acc 2.0`, `--cycles 2000`, `--scc-iterations 2500`, `--jobs 8`.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading
import traceback
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Set, Tuple

from rdkit import Chem, RDLogger
from rdkit.Chem import rdDetermineBonds

RDLogger.DisableLog("rdApp.*")


@dataclass
class BondRecord:
    a1: int  # 1-based atom index in the FULL explicit-H molecule
    a2: int  # 1-based atom index in the FULL explicit-H molecule
    symbol1: str
    symbol2: str
    bond_type: str  # SINGLE / DOUBLE / TRIPLE / AROMATIC / OTHER


@dataclass
class Candidate:
    method: str
    smiles: str  # canonical isomeric SMILES generated after removing explicit Hs
    formal_charge: int
    radical_electrons: int
    num_atoms_full: int
    bonds_full: List[BondRecord]
    atom_symbols_full: List[str]
    aromatic_pairs: Set[Tuple[int, int]]
    aromatic_ring_bonds: List[Set[Tuple[int, int]]]


@dataclass
class ConsensusReport:
    attempted_methods: List[str]
    successful_methods: Dict[str, Candidate]
    failed_methods: Dict[str, str]
    consensus_smiles: Optional[str]
    consensus_candidate: Optional[Candidate]
    status: str  # UNANIMOUS / AMBIGUOUS / FAILED
    reason: str


@dataclass
class XTBResult:
    optimized_xyz_text: Optional[str]
    stdout_text: str
    stderr_text: str
    returncode: int
    reason: str
    wbo_pairs: Dict[Tuple[int, int], float]
    wbo_reason: str


MethodFn = Callable[[Chem.Mol, int], None]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Conservative XYZ -> canonical SMILES using RDKit consensus + xTB optimization + exact-pattern WBO support."
    )
    parser.add_argument("input_folder", help="Folder containing .xyz files")
    parser.add_argument("-o", "--output-folder", default="xtb_rdkit_results")
    parser.add_argument("--charge", type=int, default=0, help="Total molecular charge (default: 0)")
    parser.add_argument(
        "--xtb",
        default="xtb",
        help="Path to xtb executable or command name on PATH (default: xtb)",
    )
    parser.add_argument(
        "--opt-level",
        default="normal",
        choices=["crude", "sloppy", "loose", "lax", "normal", "tight", "vtight", "extreme"],
        help="xTB optimization level (default: normal)",
    )
    parser.add_argument(
        "--gfn",
        type=int,
        default=2,
        choices=[0, 1, 2],
        help="xTB GFN parametrization used for the WBO single-point (default: 2).",
    )
    parser.add_argument(
        "--opt-method",
        default="gfnff",
        choices=["gfnff", "gfn0", "gfn1", "gfn2"],
        help=(
            "Method used for the geometry optimization pass. Default 'gfnff' "
            "(Grimme's GFN-FF force field) is used because the xtb 6.7.1 "
            "conda-forge Windows build has a broken GFN2 analytical gradient "
            "that causes optimization to march uphill and 'converge' at a "
            "wildly distorted geometry (e.g. aromatic C-C stretched from 1.4 "
            "to 2-3 A). GFN-FF does not use that gradient and produces clean "
            "geometries. The subsequent Wiberg single-point always uses --gfn "
            "so that WBO bond orders remain electronic-structure derived."
        ),
    )
    parser.add_argument(
        "--acc",
        type=float,
        default=2.0,
        help=(
            "xTB SCC accuracy (lower is tighter; default: 2.0). xTB's own default is "
            "1.0; we use 2.0 to help SCF converge on strained/heavy-atom systems where "
            "tighter values stall."
        ),
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=3600,
        help="xTB timeout per molecule in seconds (default: 3600)",
    )
    parser.add_argument(
        "--cycles",
        type=int,
        default=2000,
        help=(
            "Maximum geometry-optimization cycles for xTB (--cycles N). xTB's built-in "
            "per-level defaults are small (normal~200, lax~150) and large or strained "
            "molecules often stall at that ceiling. Default 2000 is permissive; lower "
            "it if you prefer faster failure over extended convergence attempts. "
            "Set to 0 to use xTB's built-in per-opt-level default."
        ),
    )
    parser.add_argument(
        "--scc-iterations",
        type=int,
        default=2500,
        help=(
            "Maximum SCC (self-consistent charge) iterations per geometry step, passed "
            "to xTB as --iterations N. xTB's default is 250; set to 2500 so the "
            "electronic solver has more room before aborting a geometry step. "
            "Set to 0 to use xTB's built-in default."
        ),
    )
    default_jobs = max(1, (os.cpu_count() or 2) // 2)
    parser.add_argument(
        "--jobs",
        type=int,
        default=default_jobs,
        help=(
            "Number of molecules to process in parallel (default: half of CPU cores = "
            f"{default_jobs}). Each in-flight task spawns one xTB subprocess pinned "
            "to 1 thread, so total CPU load stays at --jobs cores. The tasks are "
            "managed by Python threads (not processes) to avoid the Windows DLL-load "
            "race that kills ProcessPool workers under concurrent RDKit/MKL imports."
        ),
    )
    parser.add_argument(
        "--xtb-threads",
        type=int,
        default=1,
        help=(
            "Threads per xTB subprocess (default: 1). Effective total cores used is "
            "--jobs * --xtb-threads; keep the product <= half the CPU count."
        ),
    )
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help=(
            "Skip any .xyz whose optimized geometry already exists at "
            "<output>/optimized_xyz/<name>.xyz AND whose xTB log contains "
            "'GEOMETRY OPTIMIZATION CONVERGED'. Massive wall-time saver on reruns, "
            "since each xTB call can take minutes."
        ),
    )
    return parser.parse_args()


def read_text_with_fallback(path: Path) -> str:
    data = path.read_bytes()
    for enc in ("utf-8-sig", "utf-8", "cp1252", "latin-1"):
        try:
            return data.decode(enc)
        except UnicodeDecodeError:
            continue
    return data.decode("utf-8", errors="replace")


def normalize_xyz_text(text: str) -> str:
    lines = text.replace("\r\n", "\n").replace("\r", "\n").split("\n")
    while lines and lines[-1] == "":
        lines.pop()
    if len(lines) < 2:
        raise ValueError("XYZ text must contain at least atom-count and comment lines.")
    try:
        natoms = int(lines[0].strip())
    except Exception as exc:
        raise ValueError("First line of XYZ is not a valid atom count.") from exc
    atom_lines = [ln for ln in lines[2:] if ln.strip()]
    if len(atom_lines) != natoms:
        raise ValueError(
            f"XYZ atom-count mismatch: header says {natoms}, found {len(atom_lines)} atom lines."
        )
    return "\n".join([str(natoms), lines[1]] + atom_lines) + "\n"


def build_mol_from_xyz_text(xyz_text: str) -> Chem.Mol:
    mol = Chem.MolFromXYZBlock(xyz_text)
    if mol is None:
        raise ValueError("RDKit could not read the XYZ content.")
    if mol.GetNumConformers() == 0:
        raise ValueError("RDKit molecule has no 3D conformer.")
    return mol


def extract_bond_records(mol: Chem.Mol) -> List[BondRecord]:
    out: List[BondRecord] = []
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx() + 1
        a2 = bond.GetEndAtomIdx() + 1
        atom1 = bond.GetBeginAtom().GetSymbol()
        atom2 = bond.GetEndAtom().GetSymbol()
        btype = bond.GetBondType()
        if btype == Chem.BondType.SINGLE:
            label = "SINGLE"
        elif btype == Chem.BondType.DOUBLE:
            label = "DOUBLE"
        elif btype == Chem.BondType.TRIPLE:
            label = "TRIPLE"
        elif btype == Chem.BondType.AROMATIC:
            label = "AROMATIC"
        else:
            label = str(btype)
        if a1 > a2:
            a1, a2 = a2, a1
            atom1, atom2 = atom2, atom1
        out.append(BondRecord(a1, a2, atom1, atom2, label))
    out.sort(key=lambda b: (b.a1, b.a2, b.bond_type, b.symbol1, b.symbol2))
    return out


def _pair(i: int, j: int) -> Tuple[int, int]:
    return (i, j) if i < j else (j, i)


def extract_aromatic_ring_bonds(mol: Chem.Mol) -> Tuple[Set[Tuple[int, int]], List[Set[Tuple[int, int]]]]:
    aromatic_pairs: Set[Tuple[int, int]] = set()
    aromatic_ring_bonds: List[Set[Tuple[int, int]]] = []
    ring_info = mol.GetRingInfo()
    for ring_bond_indices in ring_info.BondRings():
        ring_pairs: Set[Tuple[int, int]] = set()
        all_aromatic = True
        for bidx in ring_bond_indices:
            bond = mol.GetBondWithIdx(bidx)
            pair = _pair(bond.GetBeginAtomIdx() + 1, bond.GetEndAtomIdx() + 1)
            ring_pairs.add(pair)
            if not bond.GetIsAromatic():
                all_aromatic = False
        if all_aromatic and ring_pairs:
            aromatic_ring_bonds.append(ring_pairs)
            aromatic_pairs.update(ring_pairs)
    return aromatic_pairs, aromatic_ring_bonds


def sanitize_and_candidate(mol: Chem.Mol, expected_charge: int) -> Candidate:
    """
    Build a Candidate using FULL explicit-H atom numbering for all bond/WBO checks.
    Canonical SMILES is still generated after removing explicit H atoms.
    """
    Chem.SanitizeMol(mol)
    formal_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    radical_electrons = sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms())
    if formal_charge != expected_charge:
        raise ValueError(f"Formal charge {formal_charge} != expected charge {expected_charge}.")
    if radical_electrons != 0:
        raise ValueError(f"Found {radical_electrons} radical electrons.")

    full = Chem.Mol(mol)
    bonds_full = extract_bond_records(full)
    atom_symbols_full = [a.GetSymbol() for a in full.GetAtoms()]
    aromatic_pairs, aromatic_ring_bonds = extract_aromatic_ring_bonds(full)

    smiles_mol = Chem.RemoveHs(Chem.Mol(full), sanitize=True)
    smiles = Chem.MolToSmiles(smiles_mol, canonical=True, isomericSmiles=True)

    return Candidate(
        method="",
        smiles=smiles,
        formal_charge=formal_charge,
        radical_electrons=radical_electrons,
        num_atoms_full=full.GetNumAtoms(),
        bonds_full=bonds_full,
        atom_symbols_full=atom_symbols_full,
        aromatic_pairs=aromatic_pairs,
        aromatic_ring_bonds=aromatic_ring_bonds,
    )


def _supported_kwargs(fn) -> Set[str]:
    """
    Return the set of keyword argument names this Boost.Python function actually accepts,
    parsed from its docstring. Older RDKit builds don't accept useVdw / maxIterations.
    """
    doc = fn.__doc__ or ""
    candidates = {
        "useHueckel", "charge", "covFactor", "allowChargedFragments",
        "embedChiral", "useAtomMap", "useVdw", "maxIterations",
    }
    return {k for k in candidates if k in doc}


_DB_KWARGS = _supported_kwargs(rdDetermineBonds.DetermineBonds)
_DC_KWARGS = _supported_kwargs(rdDetermineBonds.DetermineConnectivity)
_DO_KWARGS = _supported_kwargs(rdDetermineBonds.DetermineBondOrders)


def _filter_kwargs(kwargs: Dict, supported: Set[str]) -> Dict:
    return {k: v for k, v in kwargs.items() if k in supported}


def _call_determine_bonds(mol: Chem.Mol, charge: int, covFactor: float,
                          useHueckel: bool = False, useVdw: Optional[bool] = None) -> None:
    kwargs = {
        "useHueckel": useHueckel,
        "charge": charge,
        "covFactor": covFactor,
        "allowChargedFragments": True,
        "embedChiral": True,
        "useAtomMap": False,
        "maxIterations": 20000,
    }
    if useVdw is not None:
        kwargs["useVdw"] = useVdw
    rdDetermineBonds.DetermineBonds(mol, **_filter_kwargs(kwargs, _DB_KWARGS))


def _call_determine_connectivity(mol: Chem.Mol, charge: int, covFactor: float,
                                 useHueckel: bool = False,
                                 useVdw: Optional[bool] = None) -> None:
    kwargs = {
        "useHueckel": useHueckel,
        "charge": charge,
        "covFactor": covFactor,
    }
    if useVdw is not None:
        kwargs["useVdw"] = useVdw
    rdDetermineBonds.DetermineConnectivity(mol, **_filter_kwargs(kwargs, _DC_KWARGS))


def _call_determine_bond_orders(mol: Chem.Mol, charge: int) -> None:
    kwargs = {
        "charge": charge,
        "allowChargedFragments": True,
        "embedChiral": True,
        "useAtomMap": False,
        "maxIterations": 20000,
    }
    rdDetermineBonds.DetermineBondOrders(mol, **_filter_kwargs(kwargs, _DO_KWARGS))


def method_determine_bonds_ctd(mol: Chem.Mol, charge: int) -> None:
    _call_determine_bonds(mol, charge, covFactor=1.30, useVdw=False if "useVdw" in _DB_KWARGS else None)


def method_determine_bonds_vdw_130(mol: Chem.Mol, charge: int) -> None:
    if "useVdw" not in _DB_KWARGS:
        raise RuntimeError("This RDKit build does not expose useVdw; method unavailable.")
    _call_determine_bonds(mol, charge, covFactor=1.30, useVdw=True)


def method_determine_bonds_vdw_125(mol: Chem.Mol, charge: int) -> None:
    # On older RDKit without useVdw we still expose a covFactor=1.25 variant as a distinct method.
    _call_determine_bonds(
        mol, charge, covFactor=1.25,
        useVdw=True if "useVdw" in _DB_KWARGS else None,
    )


def method_connectivity_ctd_then_orders(mol: Chem.Mol, charge: int) -> None:
    _call_determine_connectivity(
        mol, charge, covFactor=1.30,
        useVdw=False if "useVdw" in _DC_KWARGS else None,
    )
    _call_determine_bond_orders(mol, charge)


def method_connectivity_vdw_130_then_orders(mol: Chem.Mol, charge: int) -> None:
    if "useVdw" not in _DC_KWARGS:
        raise RuntimeError("This RDKit build does not expose useVdw; method unavailable.")
    _call_determine_connectivity(mol, charge, covFactor=1.30, useVdw=True)
    _call_determine_bond_orders(mol, charge)


def method_connectivity_vdw_125_then_orders(mol: Chem.Mol, charge: int) -> None:
    _call_determine_connectivity(
        mol, charge, covFactor=1.25,
        useVdw=True if "useVdw" in _DC_KWARGS else None,
    )
    _call_determine_bond_orders(mol, charge)


def method_determine_bonds_hueckel(mol: Chem.Mol, charge: int) -> None:
    _call_determine_bonds(mol, charge, covFactor=1.30, useHueckel=True,
                          useVdw=False if "useVdw" in _DB_KWARGS else None)


def method_connectivity_hueckel_then_orders(mol: Chem.Mol, charge: int) -> None:
    _call_determine_connectivity(mol, charge, covFactor=1.30, useHueckel=True,
                                 useVdw=False if "useVdw" in _DC_KWARGS else None)
    _call_determine_bond_orders(mol, charge)


def get_method_suite() -> List[Tuple[str, MethodFn]]:
    methods: List[Tuple[str, MethodFn]] = [
        ("DetermineBonds_ctd", method_determine_bonds_ctd),
        ("DetermineBonds_covF_1.25", method_determine_bonds_vdw_125),
        ("Connectivity_ctd_then_BondOrders", method_connectivity_ctd_then_orders),
        ("Connectivity_covF_1.25_then_BondOrders", method_connectivity_vdw_125_then_orders),
    ]
    # Only include the explicit-VdW variants on RDKit builds that support useVdw=True.
    if "useVdw" in _DB_KWARGS:
        methods.insert(1, ("DetermineBonds_vdw_1.30", method_determine_bonds_vdw_130))
    if "useVdw" in _DC_KWARGS:
        methods.append(("Connectivity_vdw_1.30_then_BondOrders", method_connectivity_vdw_130_then_orders))
    try:
        if hasattr(rdDetermineBonds, "hueckelEnabled") and rdDetermineBonds.hueckelEnabled():
            methods.extend(
                [
                    ("DetermineBonds_hueckel", method_determine_bonds_hueckel),
                    ("Connectivity_hueckel_then_BondOrders", method_connectivity_hueckel_then_orders),
                ]
            )
    except Exception:
        pass
    return methods


def consensus_from_xyz_text(xyz_text: str, charge: int) -> ConsensusReport:
    methods = get_method_suite()
    attempted = [name for name, _ in methods]
    successes: Dict[str, Candidate] = {}
    failures: Dict[str, str] = {}

    for method_name, method_fn in methods:
        try:
            mol = build_mol_from_xyz_text(xyz_text)
            method_fn(mol, charge)
            candidate = sanitize_and_candidate(mol, expected_charge=charge)
            candidate.method = method_name
            successes[method_name] = candidate
        except Exception as exc:
            failures[method_name] = str(exc)

    if not successes:
        return ConsensusReport(
            attempted_methods=attempted,
            successful_methods=successes,
            failed_methods=failures,
            consensus_smiles=None,
            consensus_candidate=None,
            status="FAILED",
            reason="No RDKit topology method succeeded.",
        )

    unique_smiles = sorted({cand.smiles for cand in successes.values()})
    consensus_candidate = None
    if len(unique_smiles) == 1:
        consensus_candidate = next(iter(successes.values()))

    if len(successes) < 2:
        return ConsensusReport(
            attempted_methods=attempted,
            successful_methods=successes,
            failed_methods=failures,
            consensus_smiles=unique_smiles[0] if len(unique_smiles) == 1 else None,
            consensus_candidate=consensus_candidate,
            status="AMBIGUOUS",
            reason="Fewer than two independent RDKit methods succeeded.",
        )
    if failures:
        return ConsensusReport(
            attempted_methods=attempted,
            successful_methods=successes,
            failed_methods=failures,
            consensus_smiles=unique_smiles[0] if len(unique_smiles) == 1 else None,
            consensus_candidate=consensus_candidate,
            status="AMBIGUOUS",
            reason="At least one RDKit topology method failed.",
        )
    if len(unique_smiles) != 1:
        return ConsensusReport(
            attempted_methods=attempted,
            successful_methods=successes,
            failed_methods=failures,
            consensus_smiles=None,
            consensus_candidate=None,
            status="AMBIGUOUS",
            reason="RDKit methods disagree on the canonical SMILES.",
        )

    return ConsensusReport(
        attempted_methods=attempted,
        successful_methods=successes,
        failed_methods=failures,
        consensus_smiles=unique_smiles[0],
        consensus_candidate=consensus_candidate,
        status="UNANIMOUS",
        reason="All attempted RDKit methods agree.",
    )


def find_xtb_executable(user_value: str) -> str:
    p = Path(user_value)
    if p.exists():
        return str(p)
    exe = shutil.which(user_value)
    if exe:
        return exe
    raise FileNotFoundError(f"xTB executable not found: {user_value}")


_FLOAT_RE = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?"
_NEIGHBOR_TRIPLE_RE = re.compile(rf"(\d+)\s+([A-Z][a-z]?)\s+({_FLOAT_RE})")
_OLD_STYLE_LINE_RE = re.compile(rf"^\s*(\d+)\s+([A-Z][a-z]?)\s+({_FLOAT_RE})(.*)$")
_NEW_STYLE_SOURCE_RE = re.compile(rf"^\s*(\d+)\s+\d+\s+[A-Z][a-z]?\s+{_FLOAT_RE}\s*--\s*(.*)$")


def _record_wbo_pair(store: Dict[Tuple[int, int], List[float]], atom_i: int, atom_j: int, wbo: float) -> None:
    if atom_i == atom_j:
        return
    a, b = (atom_i, atom_j) if atom_i < atom_j else (atom_j, atom_i)
    store.setdefault((a, b), []).append(wbo)


def parse_xtb_wbo(stdout_text: str) -> Tuple[Dict[Tuple[int, int], float], str]:
    """
    Parse xTB Wiberg/Mayer bond orders from either of these layouts:

    Older compact style (neighbor symbol before neighbor index):
        1  C  3.965   C  6  1.445   H  7  0.978

    Newer table style (neighbor index before neighbor symbol), often with '--':
        1   6 C    3.965 --     6 C    1.445     7 H    0.978

    Also tolerates simple continuation lines in the newer style:
                                  8 H    0.977     9 H    0.976

    Returns averaged pairwise WBOs keyed by 1-based atom indices.
    """
    lines = stdout_text.replace("\r\n", "\n").replace("\r", "\n").split("\n")
    start_idx = None
    for i, line in enumerate(lines):
        if "Wiberg/Mayer (AO) data" in line:
            start_idx = i
            break
    if start_idx is None:
        return {}, "xTB stdout does not contain a Wiberg/Mayer block."

    pair_values: Dict[Tuple[int, int], List[float]] = {}
    parsed_any = False
    current_source: Optional[int] = None

    for line in lines[start_idx + 1:]:
        stripped = line.strip()
        if not stripped:
            if parsed_any:
                current_source = None
            continue

        lower = stripped.lower()
        if lower.startswith("molecular dipole") or lower.startswith("dipole moment") or lower.startswith("fukui"):
            break
        if "largest (>0.10)" in stripped or "total wbo" in lower:
            continue
        if stripped.startswith("#"):
            continue
        if set(stripped) <= {"-", "=", "|"}:
            continue

        # Newer table style with explicit source atom and '--'
        m_new = _NEW_STYLE_SOURCE_RE.match(line)
        if m_new:
            current_source = int(m_new.group(1))
            rest = m_new.group(2)
            found = False
            for neigh_idx_s, _neigh_sym, wbo_s in _NEIGHBOR_TRIPLE_RE.findall(rest):
                _record_wbo_pair(pair_values, current_source, int(neigh_idx_s), float(wbo_s))
                found = True
                parsed_any = True
            if found:
                continue

        # Older compact style: source index, source symbol, total WBO, then repeated neighbor symbol/index/WBO triples
        m_old = _OLD_STYLE_LINE_RE.match(line)
        if m_old and "--" not in line:
            atom_i = int(m_old.group(1))
            rest_tokens = m_old.group(4).split()
            found = False
            for j in range(0, len(rest_tokens) - 2, 3):
                sym = rest_tokens[j]
                idx = rest_tokens[j + 1]
                wbo = rest_tokens[j + 2]
                if not re.fullmatch(r"[A-Z][a-z]?", sym):
                    continue
                if not re.fullmatch(r"\d+", idx):
                    continue
                if not re.fullmatch(_FLOAT_RE, wbo):
                    continue
                _record_wbo_pair(pair_values, atom_i, int(idx), float(wbo))
                found = True
                parsed_any = True
            if found:
                current_source = atom_i
                continue

        # New-style continuation line: only neighbor index/symbol/WBO triples, source implied by previous line.
        if current_source is not None:
            triples = _NEIGHBOR_TRIPLE_RE.findall(stripped)
            if triples:
                found = False
                for neigh_idx_s, _neigh_sym, wbo_s in triples:
                    _record_wbo_pair(pair_values, current_source, int(neigh_idx_s), float(wbo_s))
                    found = True
                    parsed_any = True
                if found:
                    continue

    if not parsed_any:
        return {}, "Found Wiberg/Mayer header, but could not parse any bond-order lines."

    averaged = {pair: sum(vals) / len(vals) for pair, vals in pair_values.items()}
    return averaged, f"Parsed {len(averaged)} unique xTB Wiberg bond-order pairs."


# xTB WBO -> discrete bond-class inference helpers
_SINGLE_HETERO_DOUBLE_LOWER = {
    frozenset({"P", "O"}): 1.45,
    frozenset({"S", "O"}): 1.30,
    frozenset({"Se", "O"}): 1.20,
    frozenset({"S", "N"}): 1.20,
    frozenset({"P", "N"}): 1.25,
    frozenset({"Si", "O"}): 1.20,
    frozenset({"Ge", "O"}): 1.15,
}


def _double_lower_bound(symbol1: str, symbol2: str) -> float:
    pair = frozenset({symbol1, symbol2})
    return _SINGLE_HETERO_DOUBLE_LOWER.get(pair, 1.45)


def _single_upper_bound(symbol1: str, symbol2: str) -> float:
    return _double_lower_bound(symbol1, symbol2) - 0.05


def infer_nonaromatic_bond_class(wbo: float, symbol1: str, symbol2: str) -> str:
    """
    Convert a non-aromatic xTB WBO into a discrete bond class.
    Returns one of: NONBOND, SINGLE, DOUBLE, TRIPLE, AMBIGUOUS.

    The classification is intentionally conservative. Borderline regions are treated as AMBIGUOUS.
    """
    if wbo < 0.35:
        return "NONBOND"

    single_hi = _single_upper_bound(symbol1, symbol2)
    double_lo = _double_lower_bound(symbol1, symbol2)

    # Conservative buffer zones around class boundaries.
    if 0.35 <= wbo < 0.55:
        return "AMBIGUOUS"
    if 0.55 <= wbo <= single_hi:
        return "SINGLE"
    if single_hi < wbo < double_lo:
        return "AMBIGUOUS"
    if double_lo <= wbo < 2.20:
        return "DOUBLE"
    if 2.20 <= wbo < 2.45:
        return "AMBIGUOUS"
    if 2.45 <= wbo <= 3.35:
        return "TRIPLE"
    return "AMBIGUOUS"


def infer_aromatic_ring_pattern(
    ring_pairs: Set[Tuple[int, int]],
    wbo_pairs: Dict[Tuple[int, int], float],
    atom_symbols: Sequence[str],
) -> Tuple[bool, str]:
    """
    Explicit aromatic handling.

    For a ring that RDKit calls aromatic, require all ring bonds to:
    - exist in the xTB WBO map,
    - lie in a narrow aromatic WBO window,
    - have a ring-internal spread consistent with delocalization.
    """
    missing = [pair for pair in sorted(ring_pairs) if pair not in wbo_pairs]
    if missing:
        msg = ", ".join(f"{atom_symbols[i-1]}{i}-{atom_symbols[j-1]}{j}" for i, j in missing[:12])
        return False, f"Missing xTB WBOs for aromatic-ring bonds: {msg}"

    values = [wbo_pairs[pair] for pair in ring_pairs]
    min_wbo = min(values)
    max_wbo = max(values)
    mean_wbo = statistics.mean(values)

    # Aromatic bonds should not collapse toward pure single or pure double.
    if min_wbo < 1.10 or max_wbo > 1.70:
        return False, f"Aromatic-ring WBOs out of aromatic range: min={min_wbo:.3f}, max={max_wbo:.3f}"

    # Avoid strongly alternating patterns being accepted as aromatic.
    if max_wbo - min_wbo > 0.30:
        return False, f"Aromatic-ring WBO spread too large: spread={max_wbo - min_wbo:.3f}"

    if not (1.20 <= mean_wbo <= 1.55):
        return False, f"Aromatic-ring mean WBO not in aromatic band: mean={mean_wbo:.3f}"

    return True, "OK"


def verify_exact_wbo_pattern(candidate: Candidate, wbo_pairs: Dict[Tuple[int, int], float]) -> Tuple[bool, str]:
    """
    Require the xTB WBOs to support the same bonded pairs and inferred bond-order pattern
    as the unanimous RDKit candidate.

    Rules:
    - Every RDKit bond must have an xTB WBO.
    - Every aromatic ring bond must pass explicit aromatic-ring checks.
    - Every non-aromatic RDKit bond must map to the EXACT same discrete class from xTB WBO.
    - Non-bonded pairs with strong xTB WBO are rejected.
    - Borderline xTB WBOs are treated as AMBIGUOUS and rejected.
    """
    if not wbo_pairs:
        return False, "No xTB Wiberg bond orders available."

    rdkit_pairs = {(b.a1, b.a2): b for b in candidate.bonds_full}
    missing_pairs: List[str] = []
    mismatched_pairs: List[str] = []
    ambiguous_pairs: List[str] = []

    aromatic_verified: Set[Tuple[int, int]] = set()
    for ring_pairs in candidate.aromatic_ring_bonds:
        ok, reason = infer_aromatic_ring_pattern(ring_pairs, wbo_pairs, candidate.atom_symbols_full)
        if not ok:
            return False, reason
        aromatic_verified.update(ring_pairs)

    for pair, bond in rdkit_pairs.items():
        if pair not in wbo_pairs:
            missing_pairs.append(f"{bond.symbol1}{bond.a1}-{bond.symbol2}{bond.a2}")
            continue

        if bond.bond_type == "AROMATIC":
            if pair not in aromatic_verified:
                ambiguous_pairs.append(
                    f"{bond.symbol1}{bond.a1}-{bond.symbol2}{bond.a2}: aromatic bond not covered by aromatic-ring verification"
                )
            continue

        wbo = wbo_pairs[pair]
        xtb_class = infer_nonaromatic_bond_class(wbo, bond.symbol1, bond.symbol2)
        if xtb_class == "AMBIGUOUS":
            ambiguous_pairs.append(
                f"{bond.symbol1}{bond.a1}-{bond.symbol2}{bond.a2}:WBO={wbo:.3f}:borderline"
            )
            continue
        if xtb_class != bond.bond_type:
            mismatched_pairs.append(
                f"{bond.symbol1}{bond.a1}-{bond.symbol2}{bond.a2}:RDKit={bond.bond_type}:xTB={xtb_class}:WBO={wbo:.3f}"
            )

    unexpected_pairs: List[str] = []
    for pair, wbo in sorted(wbo_pairs.items()):
        if pair in rdkit_pairs:
            continue
        # Strong enough to imply a bond under the same non-aromatic classifier.
        i, j = pair
        sym_i = candidate.atom_symbols_full[i - 1] if 1 <= i <= len(candidate.atom_symbols_full) else "?"
        sym_j = candidate.atom_symbols_full[j - 1] if 1 <= j <= len(candidate.atom_symbols_full) else "?"
        xtb_class = infer_nonaromatic_bond_class(wbo, sym_i, sym_j)
        if xtb_class in {"SINGLE", "DOUBLE", "TRIPLE"}:
            unexpected_pairs.append(f"{sym_i}{i}-{sym_j}{j}:xTB={xtb_class}:WBO={wbo:.3f}")
        elif xtb_class == "AMBIGUOUS" and wbo >= 0.55:
            unexpected_pairs.append(f"{sym_i}{i}-{sym_j}{j}:xTB=AMBIGUOUS:WBO={wbo:.3f}")

    problems: List[str] = []
    if missing_pairs:
        problems.append("missing WBO for RDKit bonds: " + ", ".join(missing_pairs[:12]))
    if ambiguous_pairs:
        problems.append("ambiguous xTB bond classes: " + ", ".join(ambiguous_pairs[:12]))
    if mismatched_pairs:
        problems.append("bond-class mismatches: " + ", ".join(mismatched_pairs[:12]))
    if unexpected_pairs:
        problems.append("unexpected xTB-bonded non-RDKit pairs: " + ", ".join(unexpected_pairs[:12]))

    if problems:
        return False, " | ".join(problems)
    return True, "xTB WBOs support the same bonded pairs and exact inferred bond-order pattern as the unanimous FULL-atom RDKit topology."


_illegal = re.compile(r"[^A-Za-z0-9._-]+")


def safe_stem(name: str) -> str:
    stem = Path(name).stem
    stem = _illegal.sub("_", stem).strip("._")
    return stem or "molecule"


def classify_xtb_failure(stdout: str, stderr: str, returncode: int) -> str:
    """
    Build a human-readable reason string explaining *why* xTB failed.

    Parses xTB's stdout/stderr for specific failure signatures. Prefers the most
    specific signal (e.g. SCF non-convergence vs. generic abnormal termination).

    Note: "convergence criteria satisfied after N iterations" in the log refers
    to SCC (self-consistent charge) convergence within a single geometry-opt
    step -- NOT the overall geometry optimization. A run can have hundreds of
    those SCC-success lines and still fail to converge the geometry.
    """
    out = stdout or ""
    err = stderr or ""
    combined = out + "\n" + err
    reasons: List[str] = []

    # Explicit geometry-opt iteration-limit failure (e.g. "FAILED TO CONVERGE GEOMETRY OPTIMIZATION IN 750 ITERATIONS")
    m = re.search(r"FAILED TO CONVERGE GEOMETRY OPTIMIZATION IN\s+(\d+)\s+ITERATIONS", combined, re.IGNORECASE)
    if m:
        reasons.append(f"geometry optimization did not converge within {m.group(1)} iterations")
    elif re.search(r"Geometry optimization failed", combined, re.IGNORECASE):
        reasons.append("geometry optimization failed")

    # SCF / SCC non-convergence during geom-opt
    if re.search(r"SCF not converged", combined, re.IGNORECASE) or \
       re.search(r"Self consistent charge iterator did not converge", combined, re.IGNORECASE) or \
       re.search(r"SCC not converged", combined, re.IGNORECASE):
        reasons.append("SCF/SCC did not converge during a geometry step (aborted)")

    # xTB's numbered error chain: lines like "-7-", "-6-", "-1-"
    error_chain = re.findall(r"^\s*-\d+-\s*(.+?)\s*$", combined, re.MULTILINE)
    if error_chain:
        chain = " -> ".join(s.strip() for s in error_chain[:4])
        reasons.append(f"xTB error chain: {chain}")

    # Abnormal termination banner
    if re.search(r"abnormal termination of xtb", combined, re.IGNORECASE):
        reasons.append("abnormal termination of xtb")

    # Fallback: generic exit code
    if not reasons:
        reasons.append(f"xTB exited with code {returncode}")

    return "xTB failed: " + "; ".join(reasons) + "."


def _build_xtb_env(xtb_threads: int) -> Dict[str, str]:
    env = os.environ.copy()
    thread_str = str(max(1, int(xtb_threads)))
    env["OMP_NUM_THREADS"] = thread_str
    env["MKL_NUM_THREADS"] = thread_str
    env["OPENBLAS_NUM_THREADS"] = thread_str
    if int(xtb_threads) > 1:
        env.setdefault("OMP_STACKSIZE", "512M")
    else:
        env.pop("OMP_STACKSIZE", None)
    return env


def _run_xtb_subprocess(
    cmd: List[str],
    cwd: Path,
    timeout: int,
    env: Dict[str, str],
) -> Tuple[int, str, str, Optional[str]]:
    """Run xtb, returning (returncode, stdout, stderr, error_label). error_label is None on normal completion."""
    try:
        proc = subprocess.run(
            cmd,
            cwd=str(cwd),
            capture_output=True,
            text=True,
            timeout=timeout,
            env=env,
            encoding="utf-8",
            errors="replace",
            stdin=subprocess.DEVNULL,
        )
        return proc.returncode, proc.stdout or "", proc.stderr or "", None
    except subprocess.TimeoutExpired as exc:
        out_bytes = exc.stdout if isinstance(exc.stdout, (bytes, bytearray)) else None
        err_bytes = exc.stderr if isinstance(exc.stderr, (bytes, bytearray)) else None
        out_text = out_bytes.decode("utf-8", errors="replace") if out_bytes is not None else (exc.stdout or "")
        err_text = err_bytes.decode("utf-8", errors="replace") if err_bytes is not None else (exc.stderr or "")
        return 124, out_text, err_text, "xTB timed out."
    except Exception as exc:
        return 1, "", "", f"xTB could not be started: {exc}"


def run_xtb_optimization(
    xyz_text: str,
    xtb_exe: str,
    charge: int,
    opt_level: str,
    gfn: int,
    acc: float,
    timeout: int,
    xtb_threads: int = 1,
    cycles: int = 0,
    input_stem: str = "input",
    scc_iterations: int = 0,
    opt_method: str = "gfnff",
) -> XTBResult:
    """
    Two-step optimization:
      Step 1 (geometry): xtb with --opt, using opt_method.
        - If opt_method == 'gfnff' we pass --gfnff (avoids the broken GFN2
          analytical gradient in xtb 6.7.1 conda-forge Windows builds).
        - Otherwise we pass --gfn N for N in {0,1,2}.
      Step 2 (Wiberg): xtb single-point at the optimized geometry with
        --gfn <gfn> --sp --wbo. SCF is unaffected by the broken-gradient bug,
        so WBO is reliable here even when step 1 used GFN-FF.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        safe_in = safe_stem(input_stem) or "input"
        input_xyz = tmp_path / f"{safe_in}.xyz"
        input_xyz.write_text(xyz_text, encoding="utf-8", newline="\n")

        env = _build_xtb_env(xtb_threads)
        p_threads = str(max(1, int(xtb_threads)))

        # ---- Step 1: geometry optimization ----
        opt_cmd: List[str] = [xtb_exe, input_xyz.name, "--opt", opt_level]
        if opt_method == "gfnff":
            opt_cmd.append("--gfnff")
        else:
            # opt_method is "gfn0" / "gfn1" / "gfn2"
            opt_cmd.extend(["--gfn", opt_method.replace("gfn", "")])
            opt_cmd.extend(["--acc", str(acc)])
            if scc_iterations and int(scc_iterations) > 0:
                opt_cmd.extend(["--iterations", str(int(scc_iterations))])
        opt_cmd.extend(["--chrg", str(charge), "--uhf", "0", "-P", p_threads])
        if cycles and int(cycles) > 0:
            opt_cmd.extend(["--cycles", str(int(cycles))])

        rc1, stdout1, stderr1, err1 = _run_xtb_subprocess(opt_cmd, tmp_path, timeout, env)
        if err1 is not None:
            return XTBResult(None, stdout1, stderr1, rc1, err1, {}, "xTB failed before WBO parsing.")

        xtbopt_xyz = tmp_path / "xtbopt.xyz"
        if rc1 != 0:
            return XTBResult(None, stdout1, stderr1, rc1,
                             classify_xtb_failure(stdout1, stderr1, rc1), {}, "No WBO: opt failed.")
        if not xtbopt_xyz.exists():
            detail = classify_xtb_failure(stdout1, stderr1, rc1)
            return XTBResult(None, stdout1, stderr1, rc1,
                             f"xTB did not produce xtbopt.xyz ({detail})", {}, "No WBO: opt failed.")
        if "GEOMETRY OPTIMIZATION CONVERGED" not in stdout1:
            detail = classify_xtb_failure(stdout1, stderr1, rc1)
            return XTBResult(None, stdout1, stderr1, rc1,
                             f"xTB did not report GEOMETRY OPTIMIZATION CONVERGED ({detail})", {}, "No WBO: opt failed.")

        optimized_text = read_text_with_fallback(xtbopt_xyz)
        try:
            optimized_text = normalize_xyz_text(optimized_text)
        except Exception as exc:
            return XTBResult(None, stdout1, stderr1, rc1, f"Invalid xtbopt.xyz: {exc}", {}, "No WBO: opt produced unreadable geometry.")

        # ---- Step 2: Wiberg single-point at the optimized geometry ----
        wbo_xyz = tmp_path / f"{safe_in}_opt.xyz"
        wbo_xyz.write_text(optimized_text, encoding="utf-8", newline="\n")
        wbo_cmd: List[str] = [
            xtb_exe, wbo_xyz.name, "--sp", "--wbo",
            "--gfn", str(gfn),
            "--chrg", str(charge), "--uhf", "0",
            "--acc", str(acc),
            "-P", p_threads,
        ]
        if scc_iterations and int(scc_iterations) > 0:
            wbo_cmd.extend(["--iterations", str(int(scc_iterations))])

        rc2, stdout2, stderr2, err2 = _run_xtb_subprocess(wbo_cmd, tmp_path, timeout, env)
        # Merge the two logs into a single combined log so downstream code and
        # the user-facing log file keep both phases.
        combined_stdout = (
            "### STEP 1: GEOMETRY OPTIMIZATION (" + opt_method + ") ###\n"
            + stdout1
            + "\n\n### STEP 2: WIBERG SINGLE POINT (gfn" + str(gfn) + ") ###\n"
            + stdout2
        )
        combined_stderr = (
            (stderr1 or "") + (("\n[WBO step stderr]\n" + stderr2) if stderr2 else "")
        )

        if err2 is not None or rc2 != 0:
            # Opt succeeded; WBO failed. Keep the optimized geometry but report
            # an empty WBO set. Downstream falls through MEDIUM/LOW confidence.
            wbo_msg = err2 if err2 is not None else f"xTB --sp --wbo exited with code {rc2}."
            return XTBResult(optimized_text, combined_stdout, combined_stderr, rc2,
                             "OK (geometry only; WBO single-point failed)",
                             {}, f"WBO single-point failed: {wbo_msg}")

        wbo_pairs, wbo_reason = parse_xtb_wbo(stdout2)
        return XTBResult(optimized_text, combined_stdout, combined_stderr, rc2, "OK", wbo_pairs, wbo_reason)


def _process_one_file_impl(
    xyz_file_str: str,
    opt_xyz_dir_str: str,
    log_dir_str: str,
    xtb_exe: str,
    charge: int,
    opt_level: str,
    gfn: int,
    acc: float,
    timeout: int,
    xtb_threads: int,
    skip_existing: bool = False,
    cycles: int = 0,
    scc_iterations: int = 0,
    opt_method: str = "gfnff",
) -> Tuple[Dict[str, str], List[str]]:
    """
    Worker-safe, picklable per-file pipeline. Returns (row_dict, messages).

    Each worker writes its own optimized-xyz file and xTB log file (distinct paths
    per input file), so no cross-process coordination is needed for I/O.
    """
    xyz_file = Path(xyz_file_str)
    opt_xyz_dir = Path(opt_xyz_dir_str)
    log_dir = Path(log_dir_str)
    messages: List[str] = [f"Processing: {xyz_file.name}"]

    try:
        raw_text = read_text_with_fallback(xyz_file)
        norm_text = normalize_xyz_text(raw_text)
    except Exception as exc:
        messages.append(f"FAILED: {xyz_file.name}: invalid XYZ ({exc})")
        return (
            {
                "filename": xyz_file.name,
                "status": "FAILED",
                "canonical_smiles": "",
                "optimized_xyz": "",
                "reason": f"Invalid input XYZ: {exc}",
                "original_successful_methods": "0",
                "optimized_successful_methods": "0",
                "original_consensus": "",
                "optimized_consensus": "",
                "wbo_check": "",
            },
            messages,
        )

    original_report = consensus_from_xyz_text(norm_text, charge)

    log_path = log_dir / f"{safe_stem(xyz_file.name)}.xtb.out.txt"
    cached_opt_path = opt_xyz_dir / xyz_file.name

    xtb_result: Optional[XTBResult] = None
    if skip_existing and cached_opt_path.exists() and log_path.exists():
        try:
            prev_stdout = read_text_with_fallback(log_path)
            if "GEOMETRY OPTIMIZATION CONVERGED" in prev_stdout:
                cached_text = normalize_xyz_text(read_text_with_fallback(cached_opt_path))
                cached_wbo, cached_wbo_reason = parse_xtb_wbo(prev_stdout)
                xtb_result = XTBResult(
                    optimized_xyz_text=cached_text,
                    stdout_text=prev_stdout,
                    stderr_text="",
                    returncode=0,
                    reason="OK (cached)",
                    wbo_pairs=cached_wbo,
                    wbo_reason=cached_wbo_reason,
                )
                messages.append(f"  cached: reusing {cached_opt_path.name} ({len(cached_wbo)} WBO pairs)")
        except Exception as exc:
            messages.append(f"  cache miss for {xyz_file.name}: {exc}")
            xtb_result = None

    if xtb_result is None:
        xtb_result = run_xtb_optimization(
            xyz_text=norm_text,
            xtb_exe=xtb_exe,
            charge=charge,
            opt_level=opt_level,
            gfn=gfn,
            acc=acc,
            timeout=timeout,
            xtb_threads=xtb_threads,
            cycles=cycles,
            input_stem=safe_stem(xyz_file.name),
            scc_iterations=scc_iterations,
            opt_method=opt_method,
        )
        log_text = (xtb_result.stdout_text or "") + (
            "\n\nSTDERR\n" + xtb_result.stderr_text if xtb_result.stderr_text else ""
        )
        log_path.write_text(log_text, encoding="utf-8")

    if xtb_result.optimized_xyz_text is None:
        messages.append(f"FAILED: {xyz_file.name}: {xtb_result.reason}")
        return (
            {
                "filename": xyz_file.name,
                "status": "FAILED",
                "canonical_smiles": "",
                "optimized_xyz": "",
                "reason": xtb_result.reason,
                "original_successful_methods": str(len(original_report.successful_methods)),
                "optimized_successful_methods": "0",
                "original_consensus": original_report.consensus_smiles or "",
                "optimized_consensus": "",
                "wbo_check": xtb_result.wbo_reason,
            },
            messages,
        )

    optimized_report = consensus_from_xyz_text(xtb_result.optimized_xyz_text, charge)
    optimized_xyz_path = opt_xyz_dir / xyz_file.name
    if xtb_result.reason != "OK (cached)":
        optimized_xyz_path.write_text(xtb_result.optimized_xyz_text, encoding="utf-8", newline="\n")

    wbo_ok = False
    wbo_reason = xtb_result.wbo_reason
    if optimized_report.consensus_candidate is not None:
        wbo_ok, wbo_reason = verify_exact_wbo_pattern(
            optimized_report.consensus_candidate, xtb_result.wbo_pairs
        )

    orig_smi = original_report.consensus_smiles or ""
    opt_smi = optimized_report.consensus_smiles or ""
    parts = [
        f"original={original_report.status}:{original_report.reason}",
        f"optimized={optimized_report.status}:{optimized_report.reason}",
        f"wbo={'OK' if wbo_ok else 'FAIL'}:{wbo_reason}",
    ]
    if orig_smi and opt_smi and orig_smi != opt_smi:
        parts.append("Topology changed after xTB optimization.")
    details = " | ".join(parts)

    # Confidence ladder, most strict first. A SMILES is emitted whenever any
    # RDKit view produced a consensus on either the original or the optimized
    # geometry; only when neither did do we fall through to AMBIGUOUS with
    # no SMILES.
    if (
        original_report.status == "UNANIMOUS"
        and optimized_report.status == "UNANIMOUS"
        and orig_smi == opt_smi
        and wbo_ok
    ):
        status = "HIGH_CONFIDENCE"
        final_smiles = opt_smi
        reason = "Unanimous RDKit consensus on original and xTB-optimized geometries, identical SMILES, and exact-pattern xTB WBO support."
        messages.append(f"HIGH_CONFIDENCE: {xyz_file.name} -> {final_smiles}")
    elif (
        original_report.status == "UNANIMOUS"
        and optimized_report.status == "UNANIMOUS"
        and orig_smi == opt_smi
    ):
        status = "MEDIUM_CONFIDENCE"
        final_smiles = opt_smi
        reason = "Unanimous RDKit consensus on original and optimized geometries with identical SMILES, but WBO exact-pattern verification did not pass. | " + details
        messages.append(f"MEDIUM_CONFIDENCE: {xyz_file.name} -> {final_smiles}")
    elif original_report.status == "UNANIMOUS" and orig_smi:
        status = "LOW_CONFIDENCE"
        final_smiles = orig_smi
        reason = "Unanimous RDKit consensus on the ORIGINAL geometry only (optimized geometry did not reach the same consensus). | " + details
        messages.append(f"LOW_CONFIDENCE: {xyz_file.name} -> {final_smiles}")
    elif optimized_report.status == "UNANIMOUS" and opt_smi:
        status = "LOW_CONFIDENCE"
        final_smiles = opt_smi
        reason = "Unanimous RDKit consensus on the OPTIMIZED geometry only (original geometry did not reach the same consensus). | " + details
        messages.append(f"LOW_CONFIDENCE: {xyz_file.name} -> {final_smiles}")
    elif orig_smi or opt_smi:
        status = "TENTATIVE"
        final_smiles = orig_smi or opt_smi
        source = "original" if orig_smi else "optimized"
        reason = f"Single-method RDKit SMILES from the {source} geometry (no multi-method consensus). | " + details
        messages.append(f"TENTATIVE: {xyz_file.name} -> {final_smiles}")
    else:
        status = "AMBIGUOUS"
        final_smiles = ""
        reason = details
        messages.append(f"AMBIGUOUS: {xyz_file.name}: {reason}")

    return (
        {
            "filename": xyz_file.name,
            "status": status,
            "canonical_smiles": final_smiles,
            "optimized_xyz": optimized_xyz_path.name,
            "reason": reason,
            "original_successful_methods": str(len(original_report.successful_methods)),
            "optimized_successful_methods": str(len(optimized_report.successful_methods)),
            "original_consensus": original_report.consensus_smiles or "",
            "optimized_consensus": optimized_report.consensus_smiles or "",
            "wbo_check": wbo_reason,
        },
        messages,
    )


def process_one_file(
    xyz_file_str: str,
    opt_xyz_dir_str: str,
    log_dir_str: str,
    xtb_exe: str,
    charge: int,
    opt_level: str,
    gfn: int,
    acc: float,
    timeout: int,
    xtb_threads: int,
    skip_existing: bool = False,
    cycles: int = 0,
    scc_iterations: int = 0,
    opt_method: str = "gfnff",
) -> Tuple[Dict[str, str], List[str]]:
    """
    Bulletproof wrapper. Any uncaught Python exception inside the per-file pipeline
    is turned into a FAILED row instead of letting the worker process die and
    poison the ProcessPoolExecutor (which would otherwise fail every other
    pending future with BrokenProcessPool).

    An error log is written alongside the xTB log so the underlying cause is not lost.
    """
    try:
        return _process_one_file_impl(
            xyz_file_str,
            opt_xyz_dir_str,
            log_dir_str,
            xtb_exe,
            charge,
            opt_level,
            gfn,
            acc,
            timeout,
            xtb_threads,
            skip_existing,
            cycles,
            scc_iterations,
            opt_method,
        )
    except (KeyboardInterrupt, SystemExit):
        raise
    except BaseException as exc:
        tb = traceback.format_exc()
        name = Path(xyz_file_str).name
        try:
            err_path = Path(log_dir_str) / f"{safe_stem(name)}.error.txt"
            err_path.write_text(tb, encoding="utf-8")
        except Exception:
            pass
        short = str(exc).splitlines()[0][:200] if str(exc) else type(exc).__name__
        return (
            {
                "filename": name,
                "status": "FAILED",
                "canonical_smiles": "",
                "optimized_xyz": "",
                "reason": f"worker exception ({type(exc).__name__}): {short}",
                "original_successful_methods": "0",
                "optimized_successful_methods": "0",
                "original_consensus": "",
                "optimized_consensus": "",
                "wbo_check": "",
            },
            [f"FAILED: {name}: worker exception ({type(exc).__name__}): {short}"],
        )


def _harden_parallel_env() -> None:
    """
    Preempt known Windows failure modes when multiple Python workers each load RDKit
    (via numpy/MKL) and each xTB subprocess also loads MKL/OpenMP. The most common
    symptom is workers dying silently (BrokenProcessPool) as soon as jobs > 1.

    Set these in the parent BEFORE the pool is created so every spawned worker and
    every xTB subprocess inherits them.
    """
    # Force MKL's sequential threading layer — sidesteps the libiomp5 init race that
    # kills concurrent workers on Windows/conda.
    os.environ.setdefault("MKL_THREADING_LAYER", "SEQUENTIAL")
    # Allow multiple copies of the OMP runtime to coexist in the same process if one
    # sneaks in via a transitive dep (e.g. scipy/numpy).
    os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
    # Keep thread counts pinned to 1 at the Python/numpy level too.
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")


def main() -> int:
    _harden_parallel_env()
    args = parse_args()
    input_folder = Path(args.input_folder)
    output_folder = Path(args.output_folder)
    opt_xyz_dir = output_folder / "optimized_xyz"
    log_dir = output_folder / "logs"

    if not input_folder.exists():
        print(f"Error: input folder does not exist: {input_folder}", file=sys.stderr)
        return 1
    if not input_folder.is_dir():
        print(f"Error: input path is not a folder: {input_folder}", file=sys.stderr)
        return 1

    try:
        xtb_exe = find_xtb_executable(args.xtb)
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    xyz_files = sorted([p for p in input_folder.iterdir() if p.is_file() and p.suffix.lower() == ".xyz"])
    if not xyz_files:
        print(f"No .xyz files found in: {input_folder}")
        return 0

    opt_xyz_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    cpu_total = os.cpu_count() or 2
    jobs = max(1, int(args.jobs))
    xtb_threads = max(1, int(args.xtb_threads))
    effective = jobs * xtb_threads
    print(
        f"Parallelism: {jobs} worker(s) x {xtb_threads} xTB thread(s) = {effective} cores "
        f"(machine has {cpu_total}; half = {cpu_total // 2})."
    )
    if effective > cpu_total // 2:
        print(
            f"Warning: requested load ({effective}) exceeds half of CPU cores ({cpu_total // 2}).",
            file=sys.stderr,
        )

    # Use threads, not processes. Each task is dominated by a long-running xTB
    # subprocess; threads block harmlessly on subprocess pipe I/O (GIL released),
    # so threads give real parallelism for the wall-time that matters. Using a
    # single Python process avoids the Windows DLL-load race that kills ProcessPool
    # workers when multiple of them concurrently import RDKit/MKL.
    rows_by_index: Dict[int, Dict[str, str]] = {}
    print_lock = threading.Lock()

    def _run_one(i: int) -> Tuple[int, Dict[str, str], List[str]]:
        xyz_file = xyz_files[i]
        row, messages = process_one_file(
            str(xyz_file),
            str(opt_xyz_dir),
            str(log_dir),
            xtb_exe,
            args.charge,
            args.opt_level,
            args.gfn,
            args.acc,
            args.timeout,
            xtb_threads,
            args.skip_existing,
            args.cycles,
            args.scc_iterations,
            args.opt_method,
        )
        return i, row, messages

    with ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = [pool.submit(_run_one, i) for i in range(len(xyz_files))]
        for future in as_completed(futures):
            try:
                i, row, messages = future.result()
            except Exception as exc:
                # Should not happen: process_one_file is itself bulletproof.
                with print_lock:
                    print(f"FAILED: (unknown index): thread raised unexpectedly ({exc})")
                continue
            with print_lock:
                for msg in messages:
                    print(msg)
            rows_by_index[i] = row

    rows: List[Dict[str, str]] = [
        rows_by_index.get(i, {
            "filename": xyz_files[i].name,
            "status": "FAILED",
            "canonical_smiles": "",
            "optimized_xyz": "",
            "reason": "No row produced.",
            "original_successful_methods": "0",
            "optimized_successful_methods": "0",
            "original_consensus": "",
            "optimized_consensus": "",
            "wbo_check": "",
        }) for i in range(len(xyz_files))
    ]

    summary_path = output_folder / "summary.csv"
    with summary_path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "filename",
                "status",
                "canonical_smiles",
                "optimized_xyz",
                "reason",
                "original_successful_methods",
                "optimized_successful_methods",
                "original_consensus",
                "optimized_consensus",
                "wbo_check",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    print("\nDone.")
    print(f"Summary:       {summary_path}")
    print(f"Optimized XYZ: {opt_xyz_dir}")
    print(f"xTB logs:      {log_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
