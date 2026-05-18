# SASA 项目二

## 项目简介

本仓库用于实现项目二中与蛋白质-蛋白质相互作用界面相关的结构预处理流程，当前重点包括两部分：

- 基于结构的 `SASA` 计算
- 基于残基层 `ΔSASA` 的界面标签生成

在核心算法之外，仓库还整理了可直接复用的复合物数据集、批量标签结果、阈值统计，以及供后续分类模型使用的残基层训练主表。

## 项目目标

本项目的目标是利用蛋白复合物在结合前后表面暴露程度的变化，识别可能参与相互作用的界面残基。

整体流程如下：

```text
PDB 复合物结构
-> SASA 计算
-> 目标链 apo / holo 暴露面积对比
-> ΔSASA 计算
-> 界面残基二分类标签生成
-> 数据集汇总
-> 下游模型训练输入表
```

其中：

- `SASA` 表示溶剂可及表面积
- `ΔSASA = SASA_apo - SASA_holo`
- 当某个残基在结合后被遮挡较多时，其 `ΔSASA` 会变大，更可能是界面残基

## 方法概述

### 1. SASA 计算

`SASA` 模块基于 `Shrake-Rupley` 打点法实现，主要完成：

- 解析 `PDB` 中的原子坐标与残基信息
- 读取球面采样点文件 `Dot.txt`
- 判断每个原子的表面采样点是否被邻近原子遮挡
- 将原子级表面积聚合到残基和链两个层面

### 2. ΔSASA 标签生成

对于复合物中的目标链：

- 仅保留目标链，计算 `apo` 状态残基层 `SASA`
- 保留目标链及其配对链，计算 `holo` 状态残基层 `SASA`
- 逐残基计算 `ΔSASA`
- 在多个阈值下生成界面标签

当前仓库中使用的阈值为：

- `0.5`
- `1.0`
- `2.0`
- `5.0`

## 仓库结构

```text
README.md
requirements.txt
data/
  raw/
    examples/                  # 示例 PDB、Dot 点集、图片等原始文件
    pdb_complexes/             # 复合物 PDB 数据集
  processed/
    examples/                  # 示例输出
    interface_labels_per_complex/
    threshold_stats_per_complex/
    complex_manifest.csv
    interface_labels_all.csv
    ml_residue_dataset.csv
    threshold_statistics_by_complex.csv
    threshold_statistics_overall.csv
docs/
  README.md
  module_1_sasa.md
  module_2_delta_sasa.md
  pipeline_overview.md
  project_spec.md
scripts/
  run_sasa_example.py
  run_delta_sasa_example.py
  run_collect_dataset.py
  run_batch_labeling.py
  run_prepare_mlp_dataset.py
src/
  sasa_project/
```

## 核心模块

- [sasa.py](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/src/sasa_project/sasa.py:1)
  负责 `PDB` 解析、球面打点、原子级 `SASA` 计算，以及残基/链级汇总。

- [delta_sasa_label.py](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/src/sasa_project/delta_sasa_label.py:1)
  负责单个复合物的 `ΔSASA` 计算和标签生成。

- [batch_generate_interface_labels.py](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/src/sasa_project/batch_generate_interface_labels.py:1)
  负责对整个复合物数据集批量生成界面标签和统计结果。

- [prepare_mlp_dataset.py](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/src/sasa_project/prepare_mlp_dataset.py:1)
  负责整理供后续分类模型使用的残基层训练主表。

## 数据资源

当前仓库已包含以下关键数据：

- `100` 个复合物结构，位于 [data/raw/pdb_complexes](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/data/raw/pdb_complexes)
- 复合物清单 [complex_manifest.csv](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/data/processed/complex_manifest.csv:1)
- 全部残基标签总表 [interface_labels_all.csv](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/data/processed/interface_labels_all.csv:1)
- 总体阈值统计表 [threshold_statistics_overall.csv](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/data/processed/threshold_statistics_overall.csv:1)
- 下游训练主表 [ml_residue_dataset.csv](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/data/processed/ml_residue_dataset.csv:1)

## 使用方式

建议在仓库根目录执行以下命令，并通过 `PYTHONPATH=src` 让脚本找到项目包。

### 运行 SASA 示例

```bash
PYTHONPATH=src python3 scripts/run_sasa_example.py
```

### 运行单个复合物的 ΔSASA 标签示例

```bash
PYTHONPATH=src python3 scripts/run_delta_sasa_example.py --target-chain C --partner-chains D
```

### 批量生成复合物界面标签

```bash
PYTHONPATH=src python3 scripts/run_batch_labeling.py
```

### 生成下游模型使用的训练主表

```bash
PYTHONPATH=src python3 scripts/run_prepare_mlp_dataset.py --default-threshold 2.0
```

## 输出结果

仓库中的主要输出包括：

- [data/processed/examples](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/data/processed/examples)
  单结构 `SASA` 与单复合物 `ΔSASA` 的示例结果

- [data/processed/interface_labels_per_complex](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/data/processed/interface_labels_per_complex)
  每个复合物对应的残基标签文件

- [data/processed/threshold_stats_per_complex](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/data/processed/threshold_stats_per_complex)
  每个复合物对应的阈值统计结果

- [interface_labels_all.csv](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/data/processed/interface_labels_all.csv:1)
  全数据集范围内的残基标签总表

- [ml_residue_dataset.csv](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/data/processed/ml_residue_dataset.csv:1)
  面向后续分类模型的训练输入表，默认标签为 `ΔSASA > 2.0`

## 文档说明

- [docs/README.md](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/docs/README.md:1)
  文档导航

- [docs/module_1_sasa.md](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/docs/module_1_sasa.md:1)
  第一部分 `SASA` 模块说明

- [docs/module_2_delta_sasa.md](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/docs/module_2_delta_sasa.md:1)
  第二部分 `ΔSASA` 模块说明

- [docs/pipeline_overview.md](/Users/kzh/Documents/MyWorkspace/02_Projects/SASA/docs/pipeline_overview.md:1)
  项目整体思路说明

## 当前状态

目前仓库已经完成以下结构预处理工作：

- 复合物数据收集与筛选
- 残基层 `apo / holo SASA` 计算
- `ΔSASA` 界面标签生成
- 多阈值统计
- 残基层训练主表整理

后续可以在此基础上继续接入序列特征，例如 `ESM-2`，并训练残基级分类模型。
