#!/usr/bin/env python3

import os
import sys
import numpy as np

# 添加项目根目录到Python路径
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, project_root)

from carbon_easy_builder.graphene import Graphene
from carbon_easy_builder.lammps_writer import LAMMPSWriter
from carbon_easy_builder.box import Box

def generate_holey_graphene_sheet(nx, ny, holes_config, output_name, vacuum, output_dir):
    """
    生成打孔石墨烯片
    
    Args:
        nx: x方向单元胞数量
        ny: y方向单元胞数量（实际是一半）
        holes_config: 洞的配置列表，每个元素是 [center_x, center_y, radius]
        output_name: 输出文件名
        vacuum: 真空层厚度（Å）
        output_dir: 输出目录路径
    """
    
    print(f"\n{'='*50}")
    print(f"生成打孔石墨烯: {output_name}")
    print(f"尺寸: nx={nx}, ny={ny}")
    print(f"洞的配置: {len(holes_config)} 个洞")
    print(f"{'='*50}")
    
    # 创建石墨烯片
    sheet = Graphene(nx=nx, ny=ny)
    initial_atoms = len(sheet.atom_ids)

    # 创建Box并添加石墨烯
    box = Box.init_from_graphene(sheet, vacuum=vacuum)

    # 获取石墨烯尺寸信息
    min_coords, max_coords = sheet.get_bounding_box()
    sheet_length = max_coords[0] - min_coords[0]
    sheet_width = max_coords[1] - min_coords[1]
    
    print(f"\n石墨烯片信息:")
    print(f"  初始原子数: {initial_atoms}")
    print(f"  尺寸: {sheet_length:.2f} × {sheet_width:.2f} Å")
    print(f"  X范围: {min_coords[0]:.2f} 到 {max_coords[0]:.2f} Å")
    print(f"  Y范围: {min_coords[1]:.2f} 到 {max_coords[1]:.2f} Å")

    # 打孔
    total_removed = 0
    for i, (center_x, center_y, radius) in enumerate(holes_config):
        print(f"\n打第 {i+1} 个洞:")
        print(f"  位置: ({center_x:.2f}, {center_y:.2f}, 0.0)")
        print(f"  半径: {radius:.2f} Å")
        
        atoms_before = len(sheet.atom_ids)
        hole_center = [center_x, center_y, 0.0]
        sheet.dig_hole(hole_center, radius)
        atoms_after = len(sheet.atom_ids)
        
        removed_in_this_hole = atoms_before - atoms_after
        total_removed += removed_in_this_hole
        
        print(f"  本次移除原子: {removed_in_this_hole}")
    
    print(f"\n打孔完成:")
    print(f"  总共移除原子: {total_removed}")
    print(f"  剩余原子: {len(sheet.atom_ids)}")
    print(f"  原子移除率: {total_removed/initial_atoms*100:.1f}%")
    
    # 创建LAMMPS writer
    writer = LAMMPSWriter(box)
    
    # 确保输出目录存在
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # 保存文件
    output_file = os.path.join(output_dir, f"{output_name}.data")
    writer.write_data_file(output_file)
    
    print(f"\n输出文件: {output_file}")
    print(f"文件已保存成功!")
    
    return sheet, output_file

def generate_single_hole_series(single_hole_params, vacuum, output_dir):
    """生成单个洞的系列石墨烯"""
    
    print("\n" + "="*60)
    print("生成单孔石墨烯系列")
    print("="*60)
    
    results = []
    for config in single_hole_params['configs']:
        sheet, output_file = generate_holey_graphene_sheet(
            nx=single_hole_params['nx'], 
            ny=single_hole_params['ny'], 
            holes_config=config["holes"],
            output_name=config["name"],
            vacuum=vacuum,
            output_dir=output_dir
        )
        results.append((config["name"], sheet, output_file))
    
    return results

def generate_multiple_holes_series(multiple_holes_params, vacuum, output_dir):
    """生成多孔石墨烯系列"""
    
    print("\n" + "="*60)
    print("生成多孔石墨烯系列")
    print("="*60)
    
    results = []
    for config in multiple_holes_params['configs']:
        sheet, output_file = generate_holey_graphene_sheet(
            nx=multiple_holes_params['nx'], 
            ny=multiple_holes_params['ny'], 
            holes_config=config["holes"],
            output_name=config["name"],
            vacuum=vacuum,
            output_dir=output_dir
        )
        results.append((config["name"], sheet, output_file))
    
    return results

def generate_size_study_series(size_study_params, vacuum, output_dir):
    """生成不同尺寸石墨烯的打孔研究"""
    
    print("\n" + "="*60)
    print("生成尺寸研究系列")
    print("="*60)
    
    results = []
    for config in size_study_params['size_configs']:
        # 计算每个尺寸石墨烯的中心位置
        temp_sheet = Graphene(nx=config["nx"], ny=config["ny"])
        min_coords, max_coords = temp_sheet.get_bounding_box()
        center_x = (min_coords[0] + max_coords[0]) / 2
        center_y = (min_coords[1] + max_coords[1]) / 2
        
        holes_config = [[center_x, center_y, size_study_params['fixed_hole_radius']]]
        
        sheet, output_file = generate_holey_graphene_sheet(
            nx=config["nx"], 
            ny=config["ny"], 
            holes_config=holes_config,
            output_name=f"holey_{config['name']}",
            vacuum=vacuum,
            output_dir=output_dir
        )
        results.append((config["name"], sheet, output_file))
    
    return results

def main():
    """主函数 - 在这里定义所有可调参数"""
    
    print("石墨烯打孔结构生成器")
    print("="*60)
    
    # ==================== 全局参数配置 ====================
    vacuum = 20.0  # 真空层厚度（Å）
    output_dir = "output"  # 输出目录
    
    # ==================== 单孔系列参数配置 ====================
    single_hole_params = {
        'nx': 6,  # x方向单元胞数量
        'ny': 4,   # y方向单元胞数量
        'configs': [
            {
                "name": "single_hole_null",
                "holes": [[0.0, 0.0, 0.1]]  # [x, y, radius]
            },
            {
                "name": "single_hole_medium_2", 
                "holes": [[0, 1.42, 3]]
            },
            {
                "name": "single_hole_medium_1",
                "holes": [[0.0, 0.0, 3.0]]
            }
        ]
    }
    
    # ==================== 多孔系列参数配置 ====================
    multiple_holes_params = {
        'nx': 15,  # x方向单元胞数量
        'ny': 10,  # y方向单元胞数量
        'configs': [
            {
                "name": "three_holes_line",
                "holes": [
                    [8.0, 7.0, 1.5],   # 左 [x, y, radius]
                    [15.0, 7.0, 1.5],  # 中
                    [22.0, 7.0, 1.5]   # 右
                ]
            },
            {
                "name": "four_holes_square",
                "holes": [
                    [10.0, 5.0, 1.2],   # 左下
                    [20.0, 5.0, 1.2],   # 右下
                    [10.0, 10.0, 1.2],  # 左上
                    [20.0, 10.0, 1.2]   # 右上
                ]
            },
            {
                "name": "honeycomb_pattern",
                "holes": [
                    [8.0, 6.0, 1.0],    # 底部左
                    [15.0, 6.0, 1.0],   # 底部中
                    [22.0, 6.0, 1.0],   # 底部右
                    [11.5, 9.0, 1.0],   # 顶部左
                    [18.5, 9.0, 1.0]    # 顶部右
                ]
            }
        ]
    }
    
    # ==================== 尺寸研究系列参数配置 ====================
    size_study_params = {
        'fixed_hole_radius': 3.0,  # 固定洞半径（Å）
        'size_configs': [
            {"nx": 8, "ny": 6, "name": "size_8x6"},
            {"nx": 12, "ny": 8, "name": "size_12x8"},
            {"nx": 16, "ny": 10, "name": "size_16x10"},
            {"nx": 20, "ny": 12, "name": "size_20x12"}
        ]
    }
    
    # ==================== 执行生成任务 ====================
    try:
        print("\n配置参数总览:")
        print(f"  真空层厚度: {vacuum} Å")
        print(f"  输出目录: {output_dir}")
        print(f"  单孔系列: {single_hole_params['nx']}×{single_hole_params['ny']} 单元胞")
        print(f"  多孔系列: {multiple_holes_params['nx']}×{multiple_holes_params['ny']} 单元胞")
        print(f"  尺寸研究: 固定洞半径 {size_study_params['fixed_hole_radius']} Å")
        
        # 生成不同系列的打孔石墨烯
        single_results = generate_single_hole_series(single_hole_params, vacuum, output_dir)
        # multiple_results = generate_multiple_holes_series(multiple_holes_params, vacuum, output_dir)
        # size_results = generate_size_study_series(size_study_params, vacuum, output_dir)
        
        # # 汇总结果
        # all_results = single_results + multiple_results + size_results
        
        # print(f"\n" + "="*60)
        # print(f"所有打孔石墨烯生成完成!")
        # print(f"总共生成了 {len(all_results)} 个结构")
        # print(f"="*60)
        
        # print(f"\n生成的文件列表:")
        # for name, sheet, output_file in all_results:
        #     print(f"  {name}: {len(sheet.atom_ids)} 原子 -> {output_file}")
        
        print(f"\n所有文件保存在: {output_dir}/")
        print(f"可以用LAMMPS或其他可视化软件查看这些结构")
        
    except Exception as e:
        print(f"\n错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main() 