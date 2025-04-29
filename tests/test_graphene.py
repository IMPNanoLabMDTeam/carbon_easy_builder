import unittest
import os
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
            {"nx": 10, "ny": 10, "name": "small_sheet"},
            {"nx": 20, "ny": 20, "name": "medium_sheet"},
            {"nx": 30, "ny": 30, "name": "large_sheet"}
        ]
        
        for case in test_cases:
            # 生成石墨烯片
            sheet = Graphene(nx=case["nx"], ny=case["ny"])
            
            # 使用石墨烯创建 Box
            box = Box.init_from_graphene(sheet, vacuum=20.0)
            
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

    def tearDown(self):
        # 测试结束后可以在这里添加清理代码
        pass

if __name__ == '__main__':
    unittest.main() 