# Carbon Easy Builder

一个用于构建基于碳材料结构的Python库，包括石墨烯、碳纳米管和其他碳基材料。

## 🚨 重要更新

**CNT原子数量bug已修复 (2024-12-19)**
- 修复了 `CarbonNanotube` 类中周期性边界条件处理的关键bug
- 现在 (n,n) Armchair CNT 一个周期生成正确的 4×n 个原子
- 例：(6,6) CNT 现在正确生成24个原子（之前错误生成30个）
- **建议重新生成之前创建的CNT结构**

📋 详细信息请参见 [CHANGELOG.md](./CHANGELOG.md)

## 特性

- 快速创建石墨烯片层（单层或多层）
- 构建各种手性和长度的碳纳米管
- 支持复杂的原子团簇构建和操作
- 导出结构至LAMMPS数据格式，用于分子动力学模拟
- 提供强大的盒子操作工具，支持结构修改、旋转和平移
- 可视化辅助功能，帮助验证生成的结构

## 安装

### 开发安装

如果你想参与开发或进行本地测试：

```bash
# 克隆仓库
git clone https://github.com/yourusername/carbon_easy_builder.git
cd carbon_easy_builder

# 以开发模式安装
pip install -e .
```

### 常规安装

对于一般用户：

```bash
pip install carbon_easy_builder
```

## 快速入门

以下是一些基本用法示例：

```python
from carbon_easy_builder import Graphene, CarbonNanotube, Box, LAMMPSWriter

# 创建10x10单元的石墨烯片
graphene = Graphene(nx=10, ny=10)

# 创建(5,5)手性的碳纳米管，长度为10个单位
nanotube = CarbonNanotube(n=5, p=10)

# 将石墨烯放入模拟盒子
graphene_box = Box.init_from_graphene(graphene, vacuum=20.0)
graphene_writer = LAMMPSWriter(graphene_box)
graphene_writer.write_data_file("graphene.data")

# 为纳米管创建适当大小的盒子
min_corner, max_corner = nanotube.get_bounding_box()
box_size = max_corner - min_corner + 30.0  # 添加30埃的真空层
nanotube_box = Box(tuple(box_size))
nanotube_box.add_cluster(nanotube)
nanotube_writer = LAMMPSWriter(nanotube_box)
nanotube_writer.write_data_file("nanotube.data")
```

## 高级示例

### 石墨烯-碳纳米管复合结构

```python
# 创建石墨烯片和碳纳米管
graphene = Graphene(nx=30, ny=30)
nanotube = CarbonNanotube(n=5, p=20)

# 创建足够大的盒子容纳两个结构
combined_box = Box.init_from_graphene(graphene, vacuum=80.0)

# 平移碳纳米管到合适的位置
# 例如，使纳米管垂直穿过石墨烯中心
min_corner, max_corner = nanotube.get_bounding_box()
cnt_center = (min_corner + max_corner) / 2

# 获取石墨烯的中心位置
graphene_min, graphene_max = graphene.get_bounding_box()
graphene_center = (graphene_min + graphene_max) / 2

# 计算平移向量
translation = graphene_center - cnt_center
translation[2] = 0  # 保持Z坐标不变
nanotube.translate(translation)

# 添加碳纳米管到盒子
combined_box.add_cluster(nanotube)

# 删除石墨烯中与碳纳米管重叠的原子
combined_box.delete_atoms_in_cnt(graphene, nanotube)

# 导出结构
writer = LAMMPSWriter(combined_box)
writer.write_data_file("graphene_nanotube.data")
```

### 倾斜的碳纳米管连接两片石墨烯

```python
# 参见tests/test_complex_structure.py中的test_graphene_tilted_cnt_sandwich方法
# 该示例展示了如何创建两片石墨烯被倾斜45度的碳纳米管连接的结构
```

## 项目结构

```
carbon_easy_builder/
├── __init__.py        # 包初始化
├── graphene.py        # 石墨烯类和相关功能
├── nanotube.py        # 碳纳米管类和相关功能
├── atom_cluster.py    # 原子团簇基类
├── box.py             # 模拟盒子类
└── lammps_writer.py   # LAMMPS数据文件输出工具
```

## 测试

运行测试：

```bash
# 运行所有测试
python -m unittest discover tests

# 运行特定测试
python -m unittest tests.test_graphene
```

## 文档

详细的API文档请参见[项目文档](https://yourusername.github.io/carbon_easy_builder/)。

## 依赖

- NumPy (>=1.19.0)
- SciPy (>=1.5.0)
- Matplotlib (可选，用于可视化)

## 贡献

欢迎贡献！请随时提交Pull Request或创建Issue。

## 许可

本项目采用MIT许可证 - 详情请参见LICENSE文件。 