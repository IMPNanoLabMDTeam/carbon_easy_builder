import unittest
import os
import numpy as np
from carbon_easy_builder.nanotube import CarbonNanotube
from carbon_easy_builder.lammps_writer import LAMMPSWriter
from carbon_easy_builder.box import Box

class TestNanotube(unittest.TestCase):
    def setUp(self):
        # 使用正确的类名和参数初始化碳纳米管
        self.nanotube = CarbonNanotube(n=5, p=5)
        # 创建测试输出目录
        self.test_output_dir = "test_output"
        if not os.path.exists(self.test_output_dir):
            os.makedirs(self.test_output_dir)

    def test_create_tube(self):
        # 测试不同参数的碳纳米管
        test_cases = [
            {"n": 5, "p": 10, "name": "small_tube"},
            {"n": 8, "p": 15, "name": "medium_tube"},
            {"n": 10, "p": 20, "name": "large_tube"}
        ]
        
        for case in test_cases:
            # 生成碳纳米管
            tube = CarbonNanotube(n=case["n"], p=case["p"])
            
            # 计算合适的 box 尺寸
            min_corner, max_corner = tube.get_bounding_box()
            box_size = max_corner - min_corner
            # 添加真空区
            vacuum = 30.0  # Angstroms
            box_size += vacuum * 2
            
            # 创建 Box
            box = Box(tuple(box_size))
            # 添加纳米管到 Box
            box.add_cluster(tube)
            
            # 使用 Box 创建 LAMMPSWriter
            writer = LAMMPSWriter(box)
            
            # 保存为data文件
            output_file = os.path.join(self.test_output_dir, f"nanotube_{case['name']}.data")
            writer.write_data_file(output_file)
            
            # 打印提示信息，让用户检查生成的文件
            print(f"\nGenerated carbon nanotube test case: {case['name']}")
            print(f"Output file: {output_file}")
            print(f"Please visually inspect the generated data file for correctness")
            print(f"Expected parameters: n={case['n']}, p={case['p']}\n")

    def tearDown(self):
        # 测试结束后可以在这里添加清理代码
        pass

if __name__ == '__main__':
    unittest.main() 