import warnings

import prolif.utils

# Ignorar advertencias específicas
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", message="DCDReader currently makes independent timesteps")
warnings.filterwarnings("ignore", message="The .* interaction has been superseded by a new class")

import MDAnalysis as mda
from MDAnalysis.topology.guessers import guess_types
from prolif import Fingerprint
from prolif import Molecule
from rdkit import Chem


# Cargar topología y trayectoria

u = mda.Universe("step5_input.psf", "step6.6_equilibration.dcd")

# Inferir elementos para los átomos y añadirlos como atributo
guessed_elements = guess_types(u.atoms.names)
u.add_TopologyAttr('elements', guessed_elements)


# Crear selecciones para el ligando y la proteína
ligand_selection = u.select_atoms("index 36790")
#ligand_selection = u.select_atoms("resname NFX")

protein_selection = u.select_atoms("protein")
#protein_selection = u.select_atoms("protein or resname NFX")

print(protein_selection, "\n",  ligand_selection)

#leer
ligand= Molecule.from_mda(ligand_selection, NoImplicit=False)
protein= Molecule.from_mda(protein_selection)
#protein= Molecule.from_mda(protein_selection, NoImplicit=False)
protein_near=prolif.utils.get_residues_near_ligand(ligand, protein)

# Crear el objeto Fingerprint
fp = Fingerprint()

# Ejecutar Fingerprint directamente con selecciones de MDAnalysis
if __name__ == '__main__':
#    fp.run(u.trajectory[::200], ligand, protein, residues=protein_near)
    fp.run(u.trajectory[:10], ligand, protein, residues=protein_near, metadata=True)