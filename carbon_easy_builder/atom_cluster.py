import numpy as np
from typing import List, Tuple, Optional

class AtomCluster:
    def __init__(self, positions: np.ndarray, atom_ids: List[int]):
        """
        Initialize an atom cluster with positions and atom IDs.
        
        Args:
            positions: numpy array of shape (N, 3) containing atomic positions
            atom_ids: list of atom IDs
        """
        self.positions = positions
        self.atom_ids = atom_ids
        self.num_atoms = len(atom_ids)
        
    def translate(self, vector: np.ndarray) -> None:
        """
        Translate all atoms by a given vector.
        
        Args:
            vector: 3D translation vector
        """
        self.positions += vector
        
    def rotate(self, axis: np.ndarray, angle: float, center: Optional[np.ndarray] = None) -> None:
        """
        Rotate all atoms around an axis by a given angle.
        
        Args:
            axis: 3D rotation axis vector
            angle: rotation angle in radians
            center: rotation center point (default: geometric center)
        """
        if center is None:
            center = np.mean(self.positions, axis=0)
            
        # Normalize axis
        axis = axis / np.linalg.norm(axis)
        
        # Rodrigues rotation formula
        cos_theta = np.cos(angle)
        sin_theta = np.sin(angle)
        
        # Move to origin
        positions_centered = self.positions - center
        
        # Apply rotation
        rotated = (cos_theta * positions_centered +
                  sin_theta * np.cross(axis, positions_centered) +
                  (1 - cos_theta) * np.dot(positions_centered, axis)[:, np.newaxis] * axis)
        
        # Move back
        self.positions = rotated + center
        
    def get_bounding_box(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Get the bounding box of the cluster.
        
        Returns:
            Tuple of (min_corner, max_corner) coordinates
        """
        return np.min(self.positions, axis=0), np.max(self.positions, axis=0)
    
    def delete_atoms(self, mask: np.ndarray) -> None:
        """
        Delete atoms based on a boolean mask.
        
        Args:
            mask: boolean array of length num_atoms indicating which atoms to delete
        """
        self.positions = self.positions[~mask]
        self.atom_ids = [id for id, keep in zip(self.atom_ids, ~mask) if keep]
        self.num_atoms = len(self.atom_ids) 