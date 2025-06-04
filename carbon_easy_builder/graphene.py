import numpy as np
from .atom_cluster import AtomCluster

class Graphene(AtomCluster):
    def __init__(self, nx: int, ny: int, bond_length: float = 1.42):
        """
        Initialize a rectangular graphene sheet with given dimensions.
        
        Args:
            nx: number of unit cells in x direction
            ny: half of the number of unit cells in y direction
            bond_length: C-C bond length in Angstroms (default: 1.42)
        """
        self.nx = nx
        self.ny = ny
        self.bond_length = bond_length
        
        # Unit cell vectors
        a1 = np.array([np.sqrt(3) * bond_length, 0, 0])
        a2 = np.array([np.sqrt(3) * bond_length / 2, 3 * bond_length / 2, 0])
        
        # Basis vectors for the two atoms in the unit cell
        basis = [
            np.array([0, 0, 0]),
            np.array([np.sqrt(3) * bond_length / 2, bond_length / 2, 0])
        ]
        
        # Generate positions
        positions = []
        for i in range(nx):
            for j in range(2*ny):
                # Calculate origin of current unit cell
                origin = i * a1 + j * a2
                
                # Add two atoms in the unit cell
                for b in basis:
                    pos = origin + b
                    positions.append(pos)
        
        positions = np.array(positions)
        atom_ids = list(range(1, len(positions) + 1))
        
        super().__init__(positions, atom_ids)
    
    def dig_hole(self, hole_center, radius):
        """
        Remove all atoms within a circular hole.
        
        Args:
            hole_center: array-like, coordinates of the hole center [x, y, z]
            radius: float, radius of the hole in Angstroms
        """
        hole_center = np.array(hole_center)
        
        # Calculate distances from each atom to the hole center
        distances = np.linalg.norm(self.positions - hole_center, axis=1)
        
        # Create mask for atoms to delete (distance <= radius)
        delete_mask = distances <= radius
        
        # Count atoms to be removed before deletion
        num_removed = np.sum(delete_mask)
        
        # Use parent class method to delete atoms
        self.delete_atoms(delete_mask)
        
        print(f"Removed {num_removed} atoms within radius {radius} Å of center {hole_center}")
        print(f"Remaining atoms: {len(self.atom_ids)}")
