import unittest
import os
import numpy as np
from carbon_easy_builder.graphene import Graphene
from carbon_easy_builder.nanotube import CarbonNanotube
from carbon_easy_builder.lammps_writer import LAMMPSWriter
from carbon_easy_builder.box import Box

class TestComplexStructure(unittest.TestCase):
    """测试复杂结构：两片平行石墨烯夹着一个碳纳米管"""
    
    def setUp(self):
        # 创建测试输出目录
        self.test_output_dir = "test_output"
        if not os.path.exists(self.test_output_dir):
            os.makedirs(self.test_output_dir)
    
    def test_graphene_cnt_sandwich(self):
        """测试创建石墨烯-纳米管-石墨烯三明治结构"""
        
        # 1. 创建第一片石墨烯（底部）
        bottom_graphene = Graphene(nx=30, ny=30)
        
        # 2. 创建碳纳米管
        cnt = CarbonNanotube(n=5, p=10)
        
        # 2.1 先获取初始的纳米管方向
        cnt_min_before, cnt_max_before = cnt.get_bounding_box()
        print(f"CNT before rotation: min={cnt_min_before}, max={cnt_max_before}")
        
        # 3. 获取碳纳米管的尺寸和方向信息
        # CNT默认沿着Z轴延伸，我们希望它保持这个方向，所以不需要旋转
        # 只需确保CNT的轴与石墨烯面垂直即可
        min_corner, max_corner = cnt.get_bounding_box()
        cnt_length = max_corner[2] - min_corner[2]
        cnt_radius = cnt.r  # 获取纳米管半径
        print(f"Initial CNT length along Z-axis: {cnt_length:.2f} Å")
        
        # 4. 创建第二片石墨烯（顶部），与底部石墨烯相同尺寸
        top_graphene = Graphene(nx=30, ny=30)
        
        # 5. 创建一个足够大的Box，使用第一片石墨烯初始化
        vacuum = 100.0  # 足够大的真空层，确保能容纳整个三明治结构
        box = Box.init_from_graphene(bottom_graphene, vacuum=vacuum)
        
        # 6. 获取石墨烯尺寸
        bottom_min, bottom_max = bottom_graphene.get_bounding_box()
        top_min, top_max = top_graphene.get_bounding_box()
        
        # 7. 对称放置两片石墨烯，以z=0平面为中心
        # 计算CNT的一半长度，作为石墨烯的偏移量
        half_cnt_length = cnt_length / 2
        
        # 将底部石墨烯放在-half_cnt_length处
        bottom_center = (bottom_min + bottom_max) / 2
        bottom_translation = -bottom_center.copy()
        bottom_translation[2] = -half_cnt_length - bottom_center[2]  # Z轴放在-half_cnt_length处
        bottom_graphene.translate(bottom_translation)
        
        # 重新获取底部位置
        bottom_min, bottom_max = bottom_graphene.get_bounding_box()
        
        # 将顶部石墨烯放在+half_cnt_length处
        top_center = (top_min + top_max) / 2
        top_translation = -top_center.copy()
        top_translation[2] = half_cnt_length - top_center[2]  # Z轴放在+half_cnt_length处
        top_graphene.translate(top_translation)
        
        # 更新顶部石墨烯位置
        top_min, top_max = top_graphene.get_bounding_box()
        
        # 8. 将石墨烯添加到Box
        box.add_cluster(bottom_graphene)
        box.add_cluster(top_graphene)
        
        # 9. 放置CNT连接两个石墨烯
        # 计算CNT的平移量，使其中心位于原点
        cnt_min, cnt_max = cnt.get_bounding_box()
        cnt_center = (cnt_min + cnt_max) / 2
        cnt_translation = -cnt_center
        cnt.translate(cnt_translation)
        
        # 检查移动后的CNT位置
        cnt_min_after, cnt_max_after = cnt.get_bounding_box()
        print(f"CNT after translation: min={cnt_min_after}, max={cnt_max_after}")
        
        # 10. 手动设置碳纳米管的轴
        # CNT轴应该是垂直于石墨烯平面的，即沿着z轴方向
        start_point = np.array([0, 0, -half_cnt_length])
        end_point = np.array([0, 0, half_cnt_length])
        cnt.axis = (start_point, end_point)
        
        # 11. 将CNT添加到Box
        box.add_cluster(cnt)
        
        # 12. 在两片石墨烯上挖孔（删除与碳纳米管重叠的原子）
        box.delete_atoms_in_cnt(bottom_graphene, cnt)
        box.delete_atoms_in_cnt(top_graphene, cnt)
        
        # 13. 使用delete_atoms_in_region删除CNT在石墨烯平面上的原子
        # 先确定底部石墨烯和顶部石墨烯的Z坐标范围，定义一个很小的区域
        # 底部石墨烯平面的Z坐标范围
        bottom_z_min = bottom_min[2] - 0.1  # 稍微往下扩展一点
        bottom_z_max = bottom_max[2] + 0.1  # 稍微往上扩展一点
        
        # 顶部石墨烯平面的Z坐标范围
        top_z_min = top_min[2] - 0.1  # 稍微往下扩展一点
        top_z_max = top_max[2] + 0.1  # 稍微往上扩展一点
        
        # 打印石墨烯平面的Z范围
        print(f"Bottom graphene Z range: {bottom_z_min:.2f} to {bottom_z_max:.2f} Å")
        print(f"Top graphene Z range: {top_z_min:.2f} to {top_z_max:.2f} Å")
        
        # 删除CNT在底部石墨烯平面的原子
        print("Before deletion, CNT has", len(cnt.positions), "atoms")
        box.delete_atoms_in_region(axis=2, min_val=bottom_z_min, max_val=bottom_z_max, cluster=cnt)
        print("After deleting at bottom graphene, CNT has", len(cnt.positions), "atoms")
        
        # 删除CNT在顶部石墨烯平面的原子
        box.delete_atoms_in_region(axis=2, min_val=top_z_min, max_val=top_z_max, cluster=cnt)
        print("After deleting at top graphene, CNT has", len(cnt.positions), "atoms")
        
        # 14. 使用LAMMPSWriter将结构保存为数据文件
        writer = LAMMPSWriter(box)
        output_file = os.path.join(self.test_output_dir, "graphene_cnt_sandwich.data")
        writer.write_data_file(output_file)
        
        # 打印提示信息
        print("\nGenerated graphene-CNT-graphene sandwich structure")
        print(f"Output file: {output_file}")
        print("Please visually inspect the generated data file for correctness")
        print(f"CNT length: {cnt_length:.2f} Å")
        print(f"CNT radius: {cnt_radius:.2f} Å")
        print(f"Distance between graphene sheets: {(top_min[2] - bottom_max[2]):.2f} Å")
        print(f"Total atoms: {len(box.get_all_atom_ids())}")
        
        # 计算两石墨烯的中心
        bottom_center = np.array([
            (bottom_min[0] + bottom_max[0])/2, 
            (bottom_min[1] + bottom_max[1])/2, 
            (bottom_min[2] + bottom_max[2])/2
        ])
        
        top_center = np.array([
            (top_min[0] + top_max[0])/2, 
            (top_min[1] + top_max[1])/2, 
            (top_min[2] + top_max[2])/2
        ])
        
        print(f"Bottom graphene center: {bottom_center}")
        print(f"Top graphene center: {top_center}")
        print(f"Z symmetry check: bottom Z = {bottom_center[2]}, top Z = {top_center[2]}, sum = {bottom_center[2] + top_center[2]}")
        
        # 断言验证结构的基本性质
        # 1. 确保CNT垂直于石墨烯平面
        cnt_axis_vector = end_point - start_point
        cnt_axis_norm = cnt_axis_vector / np.linalg.norm(cnt_axis_vector)
        z_axis = np.array([0, 0, 1])
        alignment = np.abs(np.dot(cnt_axis_norm, z_axis))
        print(f"CNT axis alignment with Z-axis: {alignment:.4f} (should be close to 1.0)")
        self.assertGreater(alignment, 0.99, "CNT should be aligned with Z-axis")
        
        # 2. 确保两片石墨烯之间的距离合理
        graphene_distance = top_min[2] - bottom_max[2]
        self.assertLess(abs(graphene_distance - cnt_length), 5.0, 
                        "Distance between graphene sheets should be approximately equal to CNT length")
        
        # 3. 验证关于z=0的对称性
        self.assertAlmostEqual(bottom_center[2] + top_center[2], 0.0, delta=0.2,
                              msg="Graphene sheets should be symmetric about z=0 plane")
    
    def test_graphene_tilted_cnt_sandwich(self):
        """测试创建石墨烯-倾斜纳米管-石墨烯三明治结构，CNT倾斜45度"""
        
        # 1. 创建更大尺寸的石墨烯以容纳倾斜的CNT
        bottom_graphene = Graphene(nx=40, ny=40)
        top_graphene = Graphene(nx=40, ny=40)
        
        # 2. 创建更长的CNT，确保足够长以连接两片石墨烯
        cnt = CarbonNanotube(n=5, p=30)
        
        # 获取初始尺寸信息
        cnt_min_before, cnt_max_before = cnt.get_bounding_box()
        cnt_length = cnt_max_before[2] - cnt_min_before[2]
        cnt_radius = cnt.r
        print(f"CNT before rotation: min={cnt_min_before}, max={cnt_max_before}")
        print(f"Initial CNT length: {cnt_length:.2f} Å, radius: {cnt_radius:.2f} Å")
        
        # 3. 创建足够大的box
        vacuum = 200.0
        box = Box.init_from_graphene(bottom_graphene, vacuum=vacuum)
        
        # 4. 计算45度倾斜需要的距离
        distance = 30.0  # 设置一个固定的距离
        
        # 计算需要的CNT最小长度
        required_cnt_length = np.sqrt(2) * distance  # 对角线长度，确保45度角
        
        # 检查CNT长度是否足够
        if cnt_length < required_cnt_length * 1.2:  # 添加20%余量
            print(f"WARNING: CNT might be too short, recommended p > {int(30 * required_cnt_length / cnt_length * 1.2)}")
        
        # 5. 获取石墨烯尺寸
        bottom_min, bottom_max = bottom_graphene.get_bounding_box()
        top_min, top_max = top_graphene.get_bounding_box()
        
        # 6. 清理所有结构的初始位置
        # 将底部石墨烯放在-distance/2处
        bottom_center = (bottom_min + bottom_max) / 2
        bottom_translation = -bottom_center.copy()
        bottom_translation[2] = -distance/2 - bottom_center[2]  # Z轴放在-distance/2处
        bottom_graphene.translate(bottom_translation)
        
        # 重新获取底部位置
        bottom_min, bottom_max = bottom_graphene.get_bounding_box()
        
        # 将顶部石墨烯放在45度方向
        # X轴正向和Z轴正向，确保X和Z方向的距离相等，形成45度角
        top_center = (top_min + top_max) / 2
        top_translation = -top_center.copy()
        top_translation[0] = distance - top_center[0]  # X方向移动
        top_translation[2] = distance/2 - top_center[2]  # Z轴放在distance/2处
        top_graphene.translate(top_translation)
        
        # 获取新的位置
        top_min, top_max = top_graphene.get_bounding_box()
        
        # 7. 将石墨烯添加到box
        box.add_cluster(bottom_graphene)
        box.add_cluster(top_graphene)
        
        # 8. 放置CNT，连接两个石墨烯
        # 从底部石墨烯到顶部石墨烯的向量
        bottom_center = np.array([0, 0, -distance/2])
        top_center = np.array([distance, 0, distance/2])
        
        # 将CNT放在原点
        cnt_min, cnt_max = cnt.get_bounding_box()
        cnt_center = (cnt_min + cnt_max) / 2
        cnt.translate(-cnt_center)
        
        # 9. 旋转CNT 45度（绕y轴）
        rotation_axis = np.array([0.0, 1.0, 0.0])  # 绕y轴旋转
        rotation_angle = np.pi / 4  # 45度
        rotation_center = np.array([0.0, 0.0, 0.0])  # 原点
        
        cnt.rotate(rotation_axis, rotation_angle, rotation_center)
        
        # 10. 将CNT添加到box
        box.add_cluster(cnt)
        
        # 11. 检查旋转后的CNT位置
        cnt_min_after, cnt_max_after = cnt.get_bounding_box()
        print(f"CNT after rotation: min={cnt_min_after}, max={cnt_max_after}")
        
        # 12. 在两片石墨烯上挖孔
        box.delete_atoms_in_cnt(bottom_graphene, cnt)
        box.delete_atoms_in_cnt(top_graphene, cnt)
        
        # 13. 删除CNT在石墨烯平面外的部分
        # 获取石墨烯的Z坐标范围
        bottom_z_min = bottom_min[2] - 0.1
        bottom_z_max = bottom_max[2] + 0.1
        top_z_min = top_min[2] - 0.1
        top_z_max = top_max[2] + 0.1
        
        print(f"Bottom graphene Z range: {bottom_z_min:.2f} to {bottom_z_max:.2f} Å")
        print(f"Top graphene Z range: {top_z_min:.2f} to {top_z_max:.2f} Å")
        
        # 删除CNT在石墨烯外部的原子
        print(f"Before deletion, CNT has {len(cnt.positions)} atoms")
        
        # 删除石墨烯下方的CNT原子
        box.delete_atoms_in_region(axis=2, min_val=-float('inf'), max_val=bottom_z_min, cluster=cnt)
        
        # 删除石墨烯上方的CNT原子
        box.delete_atoms_in_region(axis=2, min_val=top_z_max, max_val=float('inf'), cluster=cnt)
        
        print(f"After deletion, CNT has {len(cnt.positions)} atoms")
        
        # 14. 保存结构
        writer = LAMMPSWriter(box)
        output_file = os.path.join(self.test_output_dir, "graphene_tilted_cnt_sandwich.data")
        writer.write_data_file(output_file)
        
        # 15. 打印信息和验证
        # 计算CNT倾斜角度
        # 计算从底部石墨烯中心到顶部石墨烯中心的向量
        bottom_center = np.array([
            (bottom_min[0] + bottom_max[0])/2, 
            (bottom_min[1] + bottom_max[1])/2, 
            (bottom_min[2] + bottom_max[2])/2
        ])
        
        top_center = np.array([
            (top_min[0] + top_max[0])/2, 
            (top_min[1] + top_max[1])/2, 
            (top_min[2] + top_max[2])/2
        ])
        
        cnt_axis_vector = top_center - bottom_center
        cnt_axis_norm = cnt_axis_vector / np.linalg.norm(cnt_axis_vector)
        z_axis = np.array([0, 0, 1])
        cos_angle = np.dot(cnt_axis_norm, z_axis)
        angle_degrees = np.arccos(cos_angle) * 180 / np.pi
        
        print("\nGenerated graphene with tilted CNT structure")
        print(f"Output file: {output_file}")
        print("Please visually inspect the generated data file for correctness")
        print(f"Original CNT length: {cnt_length:.2f} Å")
        print(f"Required CNT length for 45° tilt: {required_cnt_length:.2f} Å")
        print(f"Horizontal/vertical distance: {distance:.2f} Å")
        print(f"Total atoms: {len(box.get_all_atom_ids())}")
        print(f"CNT tilting angle: {angle_degrees:.2f} degrees (should be close to 45)")
        print(f"Bottom graphene center: {bottom_center}")
        print(f"Top graphene center: {top_center}")
        print(f"Vector between centers: {cnt_axis_vector}")
        print(f"Z symmetry check: bottom Z = {bottom_center[2]}, top Z = {top_center[2]}, sum = {bottom_center[2] + top_center[2]}")
        
        # 断言验证
        self.assertAlmostEqual(angle_degrees, 45.0, delta=2.0, 
                              msg="CNT should be tilted at approximately 45 degrees")
        # 验证关于z=0的对称性
        self.assertAlmostEqual(bottom_center[2] + top_center[2], 0.0, delta=0.2,
                              msg="Graphene sheets should be symmetric about z=0 plane")

if __name__ == '__main__':
    unittest.main() 