import numpy as np
from typing import List, Tuple, Union
from .atom_cluster import AtomCluster
from .graphene import Graphene
from .nanotube import CarbonNanotube

class Box:
    """
    A class to manage atom clusters in a periodic box.
    The box is centered at the origin with dimensions (-lx/2, lx/2), (-ly/2, ly/2), (-lz/2, lz/2).
    """
    def __init__(self, box_size: Tuple[float, float, float]):
        """
        Initialize a periodic box with given dimensions.
        
        Args:
            box_size: Tuple of (lx, ly, lz) dimensions of the box
        """
        self.box_size = np.array(box_size)
        self.box_min = -self.box_size / 2
        self.box_max = self.box_size / 2
        self.clusters: List[AtomCluster] = []
        
    def add_cluster(self, cluster: AtomCluster) -> None:
        """
        Add a cluster to the box and wrap its atoms into the box.
        
        Args:
            cluster: The atom cluster to add
        """
        self.clusters.append(cluster)
        self.wrap_cluster(cluster)
        
        # If the cluster is a carbon nanotube, calculate and set its axis
        if isinstance(cluster, CarbonNanotube):
            start_point, end_point = self.get_nanotube_axis(cluster)
            cluster.axis = (start_point, end_point)
        
    def wrap_cluster(self, cluster: AtomCluster) -> None:
        """
        Wrap all atoms of a cluster into the box using periodic boundary conditions.
        
        Args:
            cluster: The cluster whose atoms should be wrapped
        """
        # Wrap positions into [box_min, box_max)
        positions = cluster.positions
        for i in range(3):
            positions[:, i] = ((positions[:, i] - self.box_min[i]) % self.box_size[i]) + self.box_min[i]
            
    def wrap_all(self) -> None:
        """
        Wrap all clusters in the box.
        """
        for cluster in self.clusters:
            self.wrap_cluster(cluster)
            
    def get_all_positions(self) -> np.ndarray:
        """
        Get positions of all atoms from all clusters.
        
        Returns:
            Array of all atom positions
        """
        return np.vstack([cluster.positions for cluster in self.clusters])
    
    def get_all_atom_ids(self) -> List[int]:
        """
        Get atom IDs of all atoms from all clusters.
        
        Returns:
            List of all atom IDs
        """
        atom_ids = []
        for cluster in self.clusters:
            atom_ids.extend(cluster.atom_ids)
        return atom_ids
    
    def delete_atoms_in_region(self, axis: int, min_val: float, max_val: float, cluster: AtomCluster = None) -> None:
        """
        Delete atoms within a specified range along a given axis.
        
        Args:
            axis: 0 for x, 1 for y, 2 for z
            min_val: minimum coordinate value
            max_val: maximum coordinate value
            cluster: specific cluster to operate on, if None then operate on all clusters
        """
        if cluster is not None:
            mask = (cluster.positions[:, axis] >= min_val) & (cluster.positions[:, axis] <= max_val)
            cluster.delete_atoms(mask)
        else:
            for cluster in self.clusters:
                mask = (cluster.positions[:, axis] >= min_val) & (cluster.positions[:, axis] <= max_val)
                cluster.delete_atoms(mask)
            
    @classmethod
    def init_from_graphene(cls, graphene: Graphene, vacuum: float = 10.0) -> 'Box':
        """
        Initialize a box from a graphene sheet with vacuum padding.
        
        Args:
            graphene: The graphene sheet to initialize from
            vacuum: Vacuum padding in Angstroms
            
        Returns:
            A new Box instance
        """
        # Get all positions of the graphene sheet
        graphene_positions = graphene.positions
        
        # Calculate x dimension
        x_max = np.max(graphene_positions[:, 0])
        # Find atoms with the same y value as the atom with maximum x
        max_x_y = graphene_positions[np.argmax(graphene_positions[:, 0]), 1]
        same_y_mask = np.abs(graphene_positions[:, 1] - max_x_y) < 1e-6
        x_min_same_y = np.min(graphene_positions[same_y_mask, 0])
        x_dim = x_max - x_min_same_y + np.sqrt(3) * graphene.bond_length
        
        # Calculate y dimension
        y_min = np.min(graphene_positions[:, 1])
        y_max = np.max(graphene_positions[:, 1])
        y_dim = y_max - y_min + graphene.bond_length
        
        # Calculate box size with vacuum padding
        box_size = np.array([x_dim,y_dim,vacuum])
        
        # Create new box
        box = cls(box_size)
        # Add graphene to box
        box.add_cluster(graphene)

        return box

    def get_nanotube_axis(self, cnt: 'CarbonNanotube') -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate the axis of a carbon nanotube by finding end atoms and their midpoints.
        
        Args:
            cnt: The carbon nanotube to calculate axis for
            
        Returns:
            Tuple of (start_point, end_point) representing the axis line segment
            
        Raises:
            ValueError: If the number of end atoms in the two groups are not equal to 2n
                      where n is the chirality index of the nanotube
        """
        # Calculate distances between all atoms
        positions = cnt.positions
        n_atoms = len(positions)
        distances = np.zeros((n_atoms, n_atoms))
        for i in range(n_atoms):
            for j in range(i+1, n_atoms):
                # Consider periodic boundary conditions
                diff = positions[j] - positions[i]
                diff = np.mod(diff + self.box_size/2, self.box_size) - self.box_size/2
                distances[i,j] = distances[j,i] = np.linalg.norm(diff)
        
        # Find atoms with 2 neighbors (end atoms)
        neighbor_count = np.sum(distances < 1.5, axis=1)  # C-C bond length ~1.42 Å
        end_atoms = positions[neighbor_count == 3]
        # Group end atoms by z-coordinate
        z_coords = end_atoms[:, 2]
        z_min, z_max = np.min(z_coords), np.max(z_coords)
        z_mid = (z_min + z_max) / 2
        
        # Split into two groups based on z-coordinate
        group1 = end_atoms[z_coords < z_mid]
        group2 = end_atoms[z_coords > z_mid]
        
        # Check if both groups have exactly 2n atoms
        expected_count = 2 * cnt.n  # n is the chirality index
        if len(group1) != expected_count or len(group2) != expected_count:
            raise ValueError(f"Invalid number of end atoms: expected {expected_count} atoms per group, "
                           f"but got {len(group1)} and {len(group2)} atoms")
        
        # Calculate midpoints of each group
        midpoint1 = np.mean(group1, axis=0)
        midpoint2 = np.mean(group2, axis=0)
        
        return midpoint1, midpoint2 

    def delete_atoms_in_cnt(self, cluster: AtomCluster, cnt: 'CarbonNanotube') -> None:
        """
        Delete atoms from cluster that are within the cross-section of a carbon nanotube.
        
        Args:
            cluster: The cluster containing atoms that may need to be deleted
            cnt: The carbon nanotube whose cross-section defines the deletion region
        """
        # Get the axis of the nanotube
        start_point, end_point = cnt.axis
        axis_vector = end_point - start_point
        axis_length = np.linalg.norm(axis_vector)
        axis_direction = axis_vector / axis_length
        
        # Calculate distances and projections for all atoms
        positions = cluster.positions
        relative_positions = positions - start_point
        projections = np.dot(relative_positions, axis_direction)
        
        # Check if atoms are within the nanotube's length
        in_length = (projections >= -0.01) & (projections <= axis_length + 0.01)
        
        # Calculate perpendicular distances
        parallel_components = np.outer(projections, axis_direction)
        perpendicular_vectors = relative_positions - parallel_components
        perpendicular_distances = np.linalg.norm(perpendicular_vectors, axis=1)
        
        # Create deletion mask
        delete_mask = in_length & (perpendicular_distances <= cnt.r)
        
        # Delete atoms
        if np.any(delete_mask):
            cluster.delete_atoms(delete_mask) 