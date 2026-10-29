"""Generate clean 2D SVG visualizations of monomer PDB structures using RDKit."""

import glob
import os
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdDetermineBonds


def pdb_to_svg(pdb_path, svg_path, size=(400, 300)):
    # Read PDB without sanitization first
    mol = Chem.MolFromPDBFile(pdb_path, removeHs=False, sanitize=False)
    if mol is None:
        raise ValueError(f"Could not read PDB file: {pdb_path}")

    # Infer bond orders from 3D coordinates
    rdDetermineBonds.DetermineBonds(mol, charge=0)

    # Remove explicit hydrogens for a cleaner 2D depiction
    mol = Chem.RemoveHs(mol)
    Chem.SanitizeMol(mol)

    # Generate 2D coordinates for depiction
    AllChem.Compute2DCoords(mol)

    # Draw to SVG
    drawer = rdMolDraw2D.MolDraw2DSVG(size[0], size[1])
    drawer.drawOptions().addStereoAnnotation = True
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()

    with open(svg_path, "w") as f:
        f.write(drawer.GetDrawingText())

    print(f"  {os.path.basename(pdb_path)} -> {os.path.basename(svg_path)}")


def main():
    monomer_dir = os.path.join(os.path.dirname(__file__), "PAMA-monomers")
    pdb_files = sorted(glob.glob(os.path.join(monomer_dir, "*.pdb")))

    if not pdb_files:
        print("No PDB files found.")
        return

    print(f"Generating SVGs in {monomer_dir}/")
    for pdb_path in pdb_files:
        name = os.path.splitext(os.path.basename(pdb_path))[0]
        svg_path = os.path.join(monomer_dir, f"{name}.svg")
        try:
            pdb_to_svg(pdb_path, svg_path)
        except Exception as e:
            print(f"  ERROR processing {os.path.basename(pdb_path)}: {e}")


if __name__ == "__main__":
    main()
