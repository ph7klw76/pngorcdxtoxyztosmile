from pathlib import Path
import csv

from DECIMER import predict_SMILES
from rdkit import Chem
from rdkit.Chem import AllChem, Draw


def image_to_smiles(image_path: str) -> str:
    """Run DECIMER on one image and return the predicted SMILES."""
    smiles = predict_SMILES(image_path)
    if not smiles or not isinstance(smiles, str):
        raise RuntimeError("DECIMER did not return a valid SMILES string.")
    return smiles.strip()


def smiles_to_xyz(
    smiles: str,
    xyz_path: str,
    check_png_path: str,
) -> str:
    """
    Convert SMILES to a 3D RDKit structure, save XYZ and a 2D check image.
    Returns canonical SMILES.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse the predicted SMILES: {smiles}")

    canonical_smiles = Chem.MolToSmiles(mol, canonical=True)

    # Save a 2D structure image for checking against the source PNG
    Draw.MolToFile(mol, check_png_path, size=(1000, 1000))

    # Add hydrogens for 3D embedding
    mol3d = Chem.AddHs(mol)

    # Generate 3D coordinates
    params = AllChem.ETKDGv3()
    params.randomSeed = 0xF00D
    cid = AllChem.EmbedMolecule(mol3d, params)

    # Fallback if embedding fails
    if cid == -1:
        params.useRandomCoords = True
        cid = AllChem.EmbedMolecule(mol3d, params)

    if cid == -1:
        raise RuntimeError("3D embedding failed.")

    # Optimize geometry if possible
    if AllChem.MMFFHasAllMoleculeParams(mol3d):
        AllChem.MMFFOptimizeMolecule(mol3d, maxIters=1000)
    elif AllChem.UFFHasAllMoleculeParams(mol3d):
        AllChem.UFFOptimizeMolecule(mol3d, maxIters=1000)

    # Save XYZ
    Chem.MolToXYZFile(mol3d, xyz_path)

    return canonical_smiles


def process_folder(folder_name: str = "molecule") -> None:
    """
    Process all PNG files in the given folder.
    Outputs:
      - XYZ files in ./molecule/xyz/
      - recognized structure PNGs in ./molecule/recognized/
      - summary CSV in ./molecule/results.csv
    """
    folder = Path(folder_name)

    if not folder.exists() or not folder.is_dir():
        raise FileNotFoundError(f"Folder not found: {folder.resolve()}")

    png_files = sorted(folder.glob("*.png"))
    if not png_files:
        print(f"No .png files found in {folder.resolve()}")
        return

    xyz_dir = folder / "xyz"
    recognized_dir = folder / "recognized"
    xyz_dir.mkdir(exist_ok=True)
    recognized_dir.mkdir(exist_ok=True)

    results = []

    for img_path in png_files:
        base_name = img_path.stem
        xyz_path = xyz_dir / f"{base_name}.xyz"
        check_png_path = recognized_dir / f"{base_name}_recognized.png"

        print(f"Processing: {img_path.name}")

        try:
            raw_smiles = image_to_smiles(str(img_path))
            canonical_smiles = smiles_to_xyz(
                smiles=raw_smiles,
                xyz_path=str(xyz_path),
                check_png_path=str(check_png_path),
            )

            results.append({
                "file": img_path.name,
                "status": "success",
                "raw_smiles": raw_smiles,
                "canonical_smiles": canonical_smiles,
                "xyz_file": str(xyz_path),
                "recognized_png": str(check_png_path),
                "error": "",
            })

            print(f"  Success -> {xyz_path.name}")

        except Exception as e:
            results.append({
                "file": img_path.name,
                "status": "failed",
                "raw_smiles": "",
                "canonical_smiles": "",
                "xyz_file": "",
                "recognized_png": "",
                "error": str(e),
            })

            print(f"  Failed -> {e}")

    # Save summary CSV
    csv_path = folder / "results.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "file",
                "status",
                "raw_smiles",
                "canonical_smiles",
                "xyz_file",
                "recognized_png",
                "error",
            ],
        )
        writer.writeheader()
        writer.writerows(results)

    print("\nDone.")
    print(f"Summary written to: {csv_path.resolve()}")
    print(f"XYZ files saved in: {xyz_dir.resolve()}")
    print(f"Recognized images saved in: {recognized_dir.resolve()}")


if __name__ == "__main__":
    process_folder("molecule")
