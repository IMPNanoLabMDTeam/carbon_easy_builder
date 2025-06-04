import unittest
import os
import numpy as np
from carbon_easy_builder.graphene import Graphene
from carbon_easy_builder.lammps_writer import LAMMPSWriter
from carbon_easy_builder.box import Box

class TestGraphene(unittest.TestCase):
    def setUp(self):
        # 初始化 Graphene 时提供 nx 和 ny 参数
        self.graphene = Graphene(nx=5, ny=5)
        # 创建测试输出目录
        self.test_output_dir = "test_output"
        if not os.path.exists(self.test_output_dir):
            os.makedirs(self.test_output_dir)

    def test_initialization(self):
        self.assertIsNotNone(self.graphene)
        # 验证初始化参数是否正确设置
        self.assertEqual(self.graphene.nx, 5)
        self.assertEqual(self.graphene.ny, 5)

    def test_create_sheet(self):
        # 测试不同尺寸的石墨烯片
        test_cases = [
            {"nx": 2, "ny": 2, "name": "2"},
            {"nx": 3, "ny": 3, "name": "3"},
            {"nx": 4, "ny": 4, "name": "4"},
            {"nx": 5, "ny": 5, "name": "5"},
            {"nx": 6, "ny": 6, "name": "6"},
            {"nx": 7, "ny": 7, "name": "7"},
            {"nx": 8, "ny": 8, "name": "8"},
            {"nx": 9, "ny": 9, "name": "9"},
            {"nx": 10, "ny": 10, "name": "10"},
        ]
        
        for case in test_cases:
            # 生成石墨烯片
            sheet = Graphene(nx=case["nx"], ny=case["ny"])
            
            # 使用石墨烯创建 Box
            box = Box.init_from_graphene(sheet, vacuum=20.0)
            # box = Box(box_size=(100, 100, 100))
            # box.add_cluster(sheet)
            
            # 使用 Box 创建 LAMMPSWriter
            writer = LAMMPSWriter(box)
            
            # 保存为data文件
            output_file = os.path.join(self.test_output_dir, f"graphene_{case['name']}.data")
            writer.write_data_file(output_file)
            
            # 打印提示信息，让用户检查生成的文件
            print(f"\nGenerated graphene sheet test case: {case['name']}")
            print(f"Output file: {output_file}")
            print(f"Please visually inspect the generated data file for correctness")
            print(f"Expected dimensions: nx={case['nx']}, ny={case['ny']}\n")

    def test_graphene_initialization(self):
        graphene = Graphene(nx=2, ny=2)
        self.assertEqual(len(graphene.positions), 16)  # 2*2*2*2 = 16 atoms
        self.assertEqual(len(graphene.atom_ids), 16)

    def test_graphene_structure(self):
        graphene = Graphene(nx=1, ny=1)
        self.assertEqual(len(graphene.positions), 4)  # 1*1*2*2 = 4 atoms
        
        # Check that all atoms are at z=0 (2D structure)
        z_coords = graphene.positions[:, 2]
        np.testing.assert_array_almost_equal(z_coords, np.zeros(4))

    def test_dig_hole_basic(self):
        """测试基本的dig_hole功能"""
        graphene = Graphene(nx=3, ny=3)
        initial_count = len(graphene.atom_ids)
        
        # 在中心挖一个小洞
        center = np.mean(graphene.positions, axis=0)
        graphene.dig_hole(center, radius=1.0)
        
        # 应该移除了一些原子
        self.assertLess(len(graphene.atom_ids), initial_count)
        self.assertEqual(len(graphene.positions), len(graphene.atom_ids))
    
    def test_dig_hole_no_atoms_removed(self):
        """测试当洞的半径很小时不应该移除原子"""
        graphene = Graphene(nx=2, ny=2)
        initial_count = len(graphene.atom_ids)
        
        # 使用一个远离所有原子的位置和很小的半径
        center = [1000, 1000, 0]  # 远离石墨烯的位置
        graphene.dig_hole(center, radius=0.1)
        
        self.assertEqual(len(graphene.atom_ids), initial_count)
    
    def test_dig_hole_remove_all(self):
        """测试当洞的半径很大时应该移除所有原子"""
        graphene = Graphene(nx=2, ny=2)
        
        # 使用很大的半径，应该移除所有原子
        center = np.mean(graphene.positions, axis=0)
        graphene.dig_hole(center, radius=100.0)
        
        self.assertEqual(len(graphene.atom_ids), 0)
        self.assertEqual(len(graphene.positions), 0)

    def tearDown(self):
        # 测试结束后可以在这里添加清理代码
        pass

if __name__ == '__main__':
    unittest.main() 