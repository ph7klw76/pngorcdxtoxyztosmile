#!/usr/bin/env python3
"""
Convert all .xyz files in a folder into 2D chemical structure .png files.

What it does:
- Reads every .xyz file in the input folder
- Infers bond connectivity and bond orders from 3D coordinates
- Creates a 2D depiction of the molecule
- Saves a .png structure image for each .xyz file
- Writes a CSV summary report

Usage examples:
    python xyz_to_structure_png.py xyz_folder
    python xyz_to_structure_png.py xyz_folder -o png_output
    python xyz_to_structure_png.py xyz_folder -o png_output --charge 0
    python xyz_to_structure_png.py xyz_folder -o png_output --keep-hs
"""

from pathlib import Path
import argparse
import csv
import sys

from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import rdDetermineBonds
from rdkit.Chem import rdDepictor
from rdkit.Chem import Draw


# Silence most RDKit warnings for cleaner output
RDLogger.DisableLog("rdApp.*")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert all .xyz files in a folder into chemical structure .png files."
    )
    parser.add_argument(
        "input_folder",
        help="Folder containing .xyz files"
    )
    parser.add_argument(
        "-o", "--output-folder",
        default="structure_png_output",
        help="Folder to save output PNG files (default: structure_png_output)"
    )
    parser.add_argument(
        "--charge",
        type=int,
        default=0,
        help="Total molecular charge used for bond perception (default: 0)"
    )
    parser.add_argument(
        "--keep-hs",
        action="store_true",
        help="Keep explicit hydrogens in the drawn structure"
    )
    parser.add_argument(
        "--width",
        type=int,
        default=1000,
        help="PNG width in pixels (default: 1000)"
    )
    parser.add_argument(
        "--height",
        type=int,
        default=1000,
        help="PNG height in pixels (default: 1000)"
    )
    return parser.parse_args()


def xyz_to_mol(xyz_file: Path, charge: int = 0) -> Chem.Mol:
    """
    Read an XYZ file and infer bonds.

    Returns an RDKit molecule.
    Raises an exception if parsing or bond determination fails.
    """
    mol = Chem.MolFromXYZFile(str(xyz_file))
    if mol is None:
        raise ValueError("RDKit could not read the XYZ file.")

    # Determine connectivity + bond orders from coordinates
    rdDetermineBonds.DetermineBonds(mol, charge=charge)

    # Sanitize the molecule
    Chem.SanitizeMol(mol)

    return mol


def draw_molecule_to_png(
    mol: Chem.Mol,
    out_png: Path,
    keep_hs: bool = False,
    width: int = 1000,
    height: int = 1000
):
    """
    Convert molecule to a 2D depiction and save as PNG.
    """
    draw_mol = Chem.Mol(mol)

    if not keep_hs:
        draw_mol = Chem.RemoveHs(draw_mol)

    rdDepictor.Compute2DCoords(draw_mol)

    Draw.MolToFile(
        draw_mol,
        str(out_png),
        size=(width, height)
    )


def convert_all_xyz(
    input_folder: Path,
    output_folder: Path,
    charge: int,
    keep_hs: bool,
    width: int,
    height: int
):
    output_folder.mkdir(parents=True, exist_ok=True)

    xyz_files = sorted(input_folder.glob("*.xyz"))

    if not xyz_files:
        print(f"No .xyz files found in: {input_folder}")
        return

    summary_rows = []

    for xyz_file in xyz_files:
        png_file = output_folder / f"{xyz_file.stem}.png"

        try:
            mol = xyz_to_mol(xyz_file, charge=charge)
            smiles = Chem.MolToSmiles(Chem.RemoveHs(Chem.Mol(mol)), canonical=True)

            draw_molecule_to_png(
                mol,
                png_file,
                keep_hs=keep_hs,
                width=width,
                height=height
            )

            summary_rows.append({
                "xyz_filename": xyz_file.name,
                "png_filename": png_file.name,
                "smiles": smiles,
                "status": "SUCCESS",
                "error": ""
            })

            print(f"SUCCESS: {xyz_file.name} -> {png_file.name}")

        except Exception as e:
            summary_rows.append({
                "xyz_filename": xyz_file.name,
                "png_filename": "",
                "smiles": "",
                "status": "FAILED",
                "error": str(e)
            })

            print(f"FAILED: {xyz_file.name}: {e}")

    # Write summary CSV
    summary_csv = output_folder / "conversion_summary.csv"
    with open(summary_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["xyz_filename", "png_filename", "smiles", "status", "error"]
        )
        writer.writeheader()
        writer.writerows(summary_rows)

    print("\nDone.")
    print(f"Output PNG folder: {output_folder}")
    print(f"Summary file:      {summary_csv}")


def main():
    args = parse_args()

    input_folder = Path(args.input_folder)
    output_folder = Path(args.output_folder)

    if not input_folder.exists():
        print(f"Error: input folder does not exist: {input_folder}")
        sys.exit(1)

    if not input_folder.is_dir():
        print(f"Error: input path is not a folder: {input_folder}")
        sys.exit(1)

    convert_all_xyz(
        input_folder=input_folder,
        output_folder=output_folder,
        charge=args.charge,
        keep_hs=args.keep_hs,
        width=args.width,
        height=args.height
    )


if __name__ == "__main__":
    main()