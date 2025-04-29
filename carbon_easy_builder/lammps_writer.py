from typing import List, Union
import numpy as np
from .atom_cluster import AtomCluster
from .graphene import Graphene
from .nanotube import CarbonNanotube
from .box import Box

class LAMMPSWriter:
    def __init__(self, box: Box = None):
        """
        Initialize a LAMMPS writer with a box containing atom clusters.
        
        Args:
            box: Box object containing the atom clusters (optional)
        """
        self.box = box
        
    def write_data_file(self, filename: str) -> None:
        """
        Write a LAMMPS data file.
        
        Args:
            filename: output file name
        """
        if self.box is None:
            raise ValueError("Box is not set. Use write_structure instead or initialize with a Box.")
            
        # Get all positions and atom IDs from the box
        all_positions = self.box.get_all_positions()
        all_atom_ids = self.box.get_all_atom_ids()
        total_atoms = len(all_positions)
        
        # Box boundaries
        xlo, xhi = -self.box.box_size[0]/2, self.box.box_size[0]/2
        ylo, yhi = -self.box.box_size[1]/2, self.box.box_size[1]/2
        zlo, zhi = -self.box.box_size[2]/2, self.box.box_size[2]/2
        

        # Write the data file
        with open(filename, 'w') as f:
            # Header
            f.write("LAMMPS data file\n\n")
            
            # Number of atoms
            f.write(f"{total_atoms} atoms\n\n")
            
            # Atom types (only carbon for now)
            f.write("1 atom types\n\n")
            
            # Box dimensions
            f.write(f"{xlo} {xhi} xlo xhi\n")
            f.write(f"{ylo} {yhi} ylo yhi\n")
            f.write(f"{zlo} {zhi} zlo zhi\n\n")
            
            # Atoms section
            f.write("Atoms\n\n")
            
            # Write atom data
            for atom_id, pos in zip(all_atom_ids, all_positions):
                f.write(f"{atom_id} 1 {pos[0]} {pos[1]} {pos[2]}\n")
                