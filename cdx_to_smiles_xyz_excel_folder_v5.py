#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import List, Optional, Tuple

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem


def sanitize_filename(name: str) -> str:
    cleaned = "".join(c if c.isalnum() or c in ("-", "_", ".") else "_" for c in name.strip())
    return cleaned[:150] if cleaned else "molecule"


def get_xyz_writer():
    return getattr(Chem, "MolToXYZFile", None)


def rdkit_can_use_cdxml_params() -> bool:
    return hasattr(Chem, "CDXMLParserParams") and hasattr(Chem, "CDXMLFormat")


def load_with_rdkit(path: Path) -> List[Chem.Mol]:
    reader = getattr(Chem, "MolsFromCDXMLFile", None)
    if reader is None:
        raise RuntimeError("This RDKit build does not expose MolsFromCDXMLFile().")

    RDLogger.DisableLog("rdApp.*")
    try:
        if rdkit_can_use_cdxml_params():
            params = Chem.CDXMLParserParams()
            params.sanitize = True
            params.removeHs = False
            if path.suffix.lower() == ".cdx":
                params.format = Chem.CDXMLFormat.CDX
            elif path.suffix.lower() == ".cdxml":
                params.format = Chem.CDXMLFormat.CDXML
            else:
                params.format = Chem.CDXMLFormat.Auto
            mols = list(reader(str(path), params))
        else:
            if path.suffix.lower() == ".cdx":
                raise RuntimeError(
                    "This RDKit build lacks the CDX/CDXML parser-parameter API needed for binary .cdx files."
                )
            mols = list(reader(str(path), sanitize=True, removeHs=False))
    finally:
        RDLogger.EnableLog("rdApp.*")

    mols = [m for m in mols if m is not None and m.GetNumAtoms() > 0]
    if not mols:
        raise RuntimeError(
            "RDKit returned no molecules. Your RDKit build may lack the optional ChemDraw extensions."
        )
    return mols


def normalize_explicit_obabel_path(explicit: Optional[str]) -> Optional[Path]:
    if not explicit:
        return None
    p = Path(explicit.strip().strip('"'))
    if p.is_dir():
        candidate = p / "obabel.exe"
        if candidate.exists():
            return candidate
        candidate = p / "bin" / "obabel.exe"
        if candidate.exists():
            return candidate
        return p
    if p.name.lower() == "obgui.exe":
        sibling = p.with_name("obabel.exe")
        if sibling.exists():
            return sibling
    return p


def candidate_obabel_paths() -> List[Path]:
    candidates: List[Path] = []

    env = os.environ.get("OBABEL_EXE", "").strip()
    if env:
        normalized = normalize_explicit_obabel_path(env)
        if normalized is not None:
            candidates.append(normalized)

    for which_name in ("obabel", "obabel.exe"):
        which = shutil.which(which_name)
        if which:
            candidates.append(Path(which))

    for p in [
        Path(r"C:\Program Files\OpenBabel-3.1.1\obabel.exe"),
        Path(r"C:\Program Files\OpenBabel-3.1.1\bin\obabel.exe"),
        Path(r"C:\Program Files (x86)\OpenBabel-3.1.1\obabel.exe"),
        Path(r"C:\Program Files (x86)\OpenBabel-3.1.1\bin\obabel.exe"),
        Path(r"C:\Program Files\OpenBabel-2.4.1\obabel.exe"),
        Path(r"C:\Program Files\OpenBabel-2.4.1\bin\obabel.exe"),
        Path(r"C:\Program Files (x86)\OpenBabel-2.4.1\obabel.exe"),
        Path(r"C:\Program Files (x86)\OpenBabel-2.4.1\bin\obabel.exe"),
    ]:
        candidates.append(p)

    seen = set()
    uniq = []
    for p in candidates:
        key = str(p).lower()
        if key not in seen:
            uniq.append(p)
            seen.add(key)
    return uniq


def is_openbabel_executable(exe: Path) -> bool:
    try:
        proc = subprocess.run(
            [str(exe), "-V"],
            capture_output=True,
            text=True,
            timeout=10,
        )
    except Exception:
        return False

    text = f"{proc.stdout}\n{proc.stderr}".lower()
    return "open babel" in text


def find_obabel(explicit: Optional[str] = None) -> Path:
    explicit_path = normalize_explicit_obabel_path(explicit)
    if explicit_path is not None:
        if explicit_path.exists() and is_openbabel_executable(explicit_path):
            return explicit_path
        raise RuntimeError(
            f"The path you provided is not a working Open Babel executable: {explicit_path}"
        )

    for p in candidate_obabel_paths():
        if p.exists() and is_openbabel_executable(p):
            return p
    raise RuntimeError(
        "Open Babel executable not found. Install Open Babel and either add the folder containing obabel.exe to PATH, "
        "set OBABEL_EXE, or pass --obabel with the full path to obabel.exe."
    )


def load_with_openbabel(path: Path, obabel_path: Optional[str]) -> List[Chem.Mol]:
    obabel = find_obabel(obabel_path)

    tmpdir = Path(tempfile.mkdtemp(prefix="obabel_"))
    try:
        sdf_path = tmpdir / "converted.sdf"
        cmd = [str(obabel), str(path.resolve()), "-osdf", "-O", str(sdf_path.resolve())]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(
                f"Open Babel failed.\nSTDERR:\n{proc.stderr.strip()}\nSTDOUT:\n{proc.stdout.strip()}"
            )
        if not sdf_path.exists() or sdf_path.stat().st_size == 0:
            raise RuntimeError(
                f"Open Babel ran but did not create a readable SDF file at {sdf_path}.\n"
                f"STDERR:\n{proc.stderr.strip()}\nSTDOUT:\n{proc.stdout.strip()}"
            )

        supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
        mols = [m for m in supplier if m is not None and m.GetNumAtoms() > 0]
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    if not mols:
        raise RuntimeError("Open Babel conversion produced no readable molecules.")
    return mols


def load_molecules(path: Path, obabel_path: Optional[str]) -> Tuple[List[Chem.Mol], str]:
    rdkit_error = None
    try:
        return load_with_rdkit(path), "rdkit"
    except Exception as e:
        rdkit_error = str(e)

    try:
        return load_with_openbabel(path, obabel_path), "openbabel"
    except Exception as e:
        raise RuntimeError(
            f"Could not parse {path.name}.\n"
            f"RDKit error: {rdkit_error}\n"
            f"Open Babel error: {e}"
        )


def canonical_smiles(mol: Chem.Mol) -> str:
    mol2 = Chem.RemoveHs(Chem.Mol(mol))
    return Chem.MolToSmiles(mol2, canonical=True, isomericSmiles=True)


def to_3d(mol: Chem.Mol) -> Chem.Mol:
    m = Chem.AddHs(Chem.Mol(mol), addCoords=True)
    Chem.SanitizeMol(m)
    params = AllChem.ETKDGv3()
    params.randomSeed = 0xF00D
    status = AllChem.EmbedMolecule(m, params)
    if status != 0:
        params.useRandomCoords = True
        status = AllChem.EmbedMolecule(m, params)
    if status != 0:
        raise RuntimeError("3D embedding failed")
    try:
        if AllChem.MMFFHasAllMoleculeParams(m):
            AllChem.MMFFOptimizeMolecule(m)
        else:
            AllChem.UFFOptimizeMolecule(m)
    except Exception:
        pass
    return m


def write_xyz(mol: Chem.Mol, out_path: Path) -> None:
    writer = get_xyz_writer()
    if writer is None:
        raise RuntimeError("This RDKit build does not expose MolToXYZFile().")
    writer(mol, str(out_path))


def process_folder(input_dir: Path, output_dir: Path, obabel_path: Optional[str] = None) -> int:
    output_dir.mkdir(parents=True, exist_ok=True)
    xyz_dir = output_dir / "xyz"
    xyz_dir.mkdir(exist_ok=True)

    files = sorted([
        *input_dir.glob("*.cdx"),
        *input_dir.glob("*.CDX"),
        *input_dir.glob("*.cdxml"),
        *input_dir.glob("*.CDXML"),
    ])
    if not files:
        print(f"No .cdx or .cdxml files found in {input_dir}", file=sys.stderr)
        return 1

    rows = []
    any_success = False

    for f in files:
        try:
            mols, parser_used = load_molecules(f, obabel_path)
            stem = sanitize_filename(f.stem)
            multi = len(mols) > 1
            for idx, mol in enumerate(mols, start=1):
                smi = canonical_smiles(mol)
                xyz_name = f"{stem}_{idx}.xyz" if multi else f"{stem}.xyz"
                xyz_path = xyz_dir / xyz_name
                mol3d = to_3d(mol)
                write_xyz(mol3d, xyz_path)
                rows.append({"filename": xyz_name, "smiles": smi})
                any_success = True
            print(f"OK: {f.name} ({len(mols)} molecule(s), parser={parser_used})")
        except Exception as e:
            print(f"Failed: {f.name}: {e}", file=sys.stderr)

    excel_path = output_dir / "smiles_results.xlsx"
    pd.DataFrame(rows, columns=["filename", "smiles"]).to_excel(excel_path, index=False)
    print(f"Excel written: {excel_path}")

    if not any_success:
        print(
            "No files were converted successfully. For binary .cdx files, use a working Open Babel obabel.exe or export from ChemDraw to .cdxml/.mol/.sdf first.",
            file=sys.stderr,
        )
        return 2
    return 0


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert all CDX/CDXML files in a folder to XYZ and an Excel file of filename + SMILES.")
    parser.add_argument("input_dir", nargs="?", default="needtoconvertxyzandsmile", help="Folder containing .cdx/.cdxml files")
    parser.add_argument("-o", "--output-dir", default="converted_output", help="Output folder")
    parser.add_argument("--obabel", default=None, help="Full path to obabel.exe if it is not on PATH")
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    if not input_dir.exists() or not input_dir.is_dir():
        print(f"Input folder not found: {input_dir}", file=sys.stderr)
        sys.exit(1)

    code = process_folder(input_dir, output_dir, args.obabel)
    sys.exit(code)


if __name__ == "__main__":
    main()
