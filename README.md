# SASA Project 2

这个仓库对应你们项目二当前已经推进完成的 `B + D` 前半段基础工作，重点放在两部分：

1. 基于 `Shrake-Rupley` 打点法的 `SASA` 计算
2. 基于复合物结构的 `ΔSASA` 界面残基标签生成

在这两部分之上，仓库已经进一步整理出：

- `100` 个可直接使用的蛋白复合物数据集
- 批量界面标签结果
- 面向后续 `MLP` 训练的残基层主表

也就是说，这个仓库现在不只是“脚本集合”，而是一个已经能支撑后续建模工作的结构化项目目录。

## 项目流程

```text
单链/复合物 PDB
-> SASA 计算
-> 目标链 apo / holo 残基暴露面积对比
-> ΔSASA 标签生成
-> 100 个复合物批量汇总
-> 训练主表输出
-> 后续可继续接 ESM-2 / MLP
```

## 当前完成情况

目前已经完成并可直接复用的内容：

- 第一部分 `SASA` 基础模块
- 第二部分 `ΔSASA` 标签生成模块
- `100` 个复合物 biological assembly 数据集
- 全部残基标签总表
- 不同阈值下的正负样本统计
- 默认阈值 `ΔSASA > 2.0` 的训练主表

当前最关键的结果文件：

- 复合物清单：[complex_manifest.csv](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/data/processed/complex_manifest.csv:1)
- 全部残基标签：[interface_labels_all.csv](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/data/processed/interface_labels_all.csv:1)
- 总体阈值统计：[threshold_statistics_overall.csv](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/data/processed/threshold_statistics_overall.csv:1)
- 训练主表：[ml_residue_dataset.csv](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/data/processed/ml_residue_dataset.csv:1)

## 目录结构

```text
README.md
requirements.txt
data/
  raw/
    examples/                  # 示例 PDB、Dot 点集、图片
    pdb_complexes/             # 100 个复合物 biological assembly PDB
  processed/
    examples/                  # 示例运行结果
    interface_labels_per_complex/
    threshold_stats_per_complex/
    complex_manifest.csv
    interface_labels_all.csv
    ml_residue_dataset.csv
    threshold_statistics_by_complex.csv
    threshold_statistics_overall.csv
docs/
  README.md                    # 文档导航
  module_1_sasa.md
  module_2_delta_sasa.md
  pipeline_overview.md
  project_spec.md
scripts/                       # 统一运行入口
src/sasa_project/              # 核心代码包
```

## 快速开始

推荐在仓库根目录执行，并通过 `PYTHONPATH=src` 让脚本找到项目包。

1. 运行单结构 `SASA` 示例

```bash
PYTHONPATH=src python3 scripts/run_sasa_example.py
```

2. 运行 `ΔSASA` 界面标签示例

```bash
PYTHONPATH=src python3 scripts/run_delta_sasa_example.py --target-chain C --partner-chains D
```

3. 重新生成 `100` 个复合物的批量标签

```bash
PYTHONPATH=src python3 scripts/run_batch_labeling.py
```

4. 重新整理给后续模型使用的训练主表

```bash
PYTHONPATH=src python3 scripts/run_prepare_mlp_dataset.py --default-threshold 2.0
```

## 代码入口

如果要继续扩展项目，优先从下面这些入口看起：

- [sasa.py](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/src/sasa_project/sasa.py:1)
  负责 PDB 解析、球面打点、原子/残基/链级 SASA 聚合。
- [delta_sasa_label.py](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/src/sasa_project/delta_sasa_label.py:1)
  负责 apo/holo 对比、`ΔSASA` 计算和单复合物标签输出。
- [batch_generate_interface_labels.py](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/src/sasa_project/batch_generate_interface_labels.py:1)
  负责对整个复合物数据集批量生成标签和统计。
- [prepare_mlp_dataset.py](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/src/sasa_project/prepare_mlp_dataset.py:1)
  负责整理残基层训练主表，便于第三位同学继续接 `ESM + MLP`。

## 文档入口

如果你是组员或老师，建议这样阅读：

- 项目整体思路：[pipeline_overview.md](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/docs/pipeline_overview.md:1)
- 第一部分说明：[module_1_sasa.md](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/docs/module_1_sasa.md:1)
- 第二部分说明：[module_2_delta_sasa.md](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/docs/module_2_delta_sasa.md:1)
- 文档导航：[docs/README.md](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/docs/README.md:1)

## 备注

当前仓库运行只依赖 Python 标准库，因此 [requirements.txt](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/requirements.txt:1) 暂时没有额外三方包。后续如果接入 `ESM-2`、`PyTorch` 或绘图模块，再把依赖补进去会更合适。
