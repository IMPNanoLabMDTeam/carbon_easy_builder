import numpy as np
from .atom_cluster import AtomCluster
from typing import List, Union
from collections import OrderedDict           # 用来去重

class CarbonNanotube(AtomCluster):
    """
    Armchair (n,n) SWCNT
      n : 手性指数
      p : 周期单胞数——决定总长  L = p·|T|   (|T| ≈ 2.46 Å)
      d : C–C 键长
    返回   4*n*p  原子 ➜  z 坐标已 wrap 进 [0, L)，首尾不再重复。
    """
    def __init__(self, n: int, p: int, d: float = 1.42):
        # ----------1. 基本几何量----------
        a1 = np.array([np.sqrt(3)*d,       0.0, 0.0])
        a2 = np.array([np.sqrt(3)*d/2, 3*d/2, 0.0])

        Ch   = n * (a1 + a2)                      # 周向矢量
        T    = a1 - a2                            # 轴向平移
        R    = 3*d*n / (2*np.pi)                  # 半径
        Ch_u = Ch / np.linalg.norm(Ch)
        T_u  = T  / np.linalg.norm(T)
        Ch_len = np.linalg.norm(Ch)
        Lz     = np.linalg.norm(T)                # ≈ 2.46 Å
        Ltot   = p * Lz                           # 非周期 CNT 总长

        # ----------2. 枚举 2-D 单胞 (n 列 × 2 行 × 2 基矢 = 4n)----------
        basis = [(0.0, 0.0), (1/3, 1/3)]          # 子晶格 A / B
        atoms2D = []
        for u in range(2*n):                        # n columns
            for v in (0, 1):                      # 2 rows
                P = u*a1 + v*a2
                for fx, fy in basis:
                    atoms2D.append(P + fx*a1 + fy*a2)

        # ----------3. 轴向复制 p 单胞 → 得到周期 CNT----------
        raw_pos = []
        for k in range(p):
            shift = k * T
            for r2 in atoms2D:
                Q   = r2 + shift
                φ   = (2*np.pi * np.dot(Q, Ch_u) / Ch_len) % (2*np.pi)
                z   =  np.dot(Q, T_u)             # 轴向坐标（未 wrap）
                raw_pos.append([R*np.cos(φ), R*np.sin(φ), z])

        # ----------4. Wrap 到 [0, Ltot) 并去重----------
        unique = OrderedDict()
        for x, y, z in raw_pos:
            z_wrapped = z % Ltot                  # 严格限制到管长内部
            # 采用 1e-5 Å 精度做键：去掉 z≈0 与 z≈Ltot 的重叠
            key = (round(x, 5), round(y, 5), round(z_wrapped, 5))
            if key not in unique:
                unique[key] = [x, y, z_wrapped]

        positions = np.asarray(list(unique.values()))
        atom_ids  = list(range(1, len(positions)+1))

        # Center the nanotube at the origin
        center = np.mean(positions, axis=0)
        positions -= center
        self.r = R
        self.n = n
        super().__init__(positions, atom_ids)

    def delete_inner_atoms(self, clusters: List[Union[AtomCluster, 'CarbonNanotube', 'Graphene']], 
                          tolerance: float = 0.1) -> None:
        """
        Delete atoms from input clusters that overlap with this nanotube.
        
        Args:
            clusters: list of other atom clusters to check for overlap
            tolerance: minimum distance between atoms to consider them overlapping (in Angstroms)
        """
        # Check overlap with each cluster
        for cluster in clusters:
            # Get bounding box of the cluster
            min_corner, max_corner = cluster.get_bounding_box()
            
            # Check if any cluster atoms are within the nanotube's bounding box
            in_box = np.all((cluster.positions >= min_corner - tolerance) & 
                          (cluster.positions <= max_corner + tolerance), axis=1)
            
            if np.any(in_box):
                # For atoms within the bounding box, check exact distances
                delete_mask = np.zeros(len(cluster.positions), dtype=bool)
                for i in range(len(cluster.positions)):
                    if in_box[i]:
                        # Calculate distances to all atoms in the nanotube
                        dists = np.linalg.norm(cluster.positions[i] - self.positions, axis=1)
                        if np.any(dists < tolerance):
                            delete_mask[i] = True
                
                # Delete overlapping atoms from the cluster
                if np.any(delete_mask):
                    cluster.delete_atoms(delete_mask)