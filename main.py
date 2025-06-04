#!/usr/bin/env python3

from carbon_easy_builder.graphene import Graphene
import numpy as np

def test_dig_hole():
    """测试 dig_hole 方法"""
    
    # 创建一个小的石墨烯片
    graphene = Graphene(nx=5, ny=3)
    
    print(f"初始原子数量: {len(graphene.atom_ids)}")
    print(f"石墨烯尺寸范围:")
    print(f"  X: {np.min(graphene.positions[:, 0]):.2f} 到 {np.max(graphene.positions[:, 0]):.2f} Å")
    print(f"  Y: {np.min(graphene.positions[:, 1]):.2f} 到 {np.max(graphene.positions[:, 1]):.2f} Å")
    
    # 选择石墨烯中心附近作为挖洞的位置
    center_x = (np.min(graphene.positions[:, 0]) + np.max(graphene.positions[:, 0])) / 2
    center_y = (np.min(graphene.positions[:, 1]) + np.max(graphene.positions[:, 1])) / 2
    hole_center = [center_x, center_y, 0]
    
    print(f"\n在位置 ({center_x:.2f}, {center_y:.2f}, 0) 挖一个半径为 3.0 Å 的洞:")
    
    # 挖洞
    graphene.dig_hole(hole_center, radius=3.0)
    
    return graphene

if __name__ == "__main__":
    test_dig_hole() 