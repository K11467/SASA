# SASA 项目二

## 项目简介

本仓库用于实现项目二中与蛋白质-蛋白质相互作用界面相关的结构预处理与下游预测流程。项目从 PDB 复合物结构出发，先实现自研 `SASA` 计算，再通过 `ΔSASA` 生成界面残基弱监督标签，最后结合 ESM-2 残基层表征和结构特征训练界面残基预测模型。

整体流程如下：

```text
PDB 复合物结构
-> 自研 SASA 计算
-> 目标链 apo / holo SASA 对比
-> ΔSASA 弱监督标签生成
-> ESM-2 残基层 embedding 提取
-> SASA / ESM / 坐标多模态数据融合
-> MLP / GCN / EGNN / cross-chain EGNN 界面残基预测
-> 特征消融与阈值敏感性分析
```

其中：

- `SASA` 表示溶剂可及表面积。
- `ΔSASA = SASA_apo - SASA_holo`。
- 当某个残基在结合后被遮挡较多时，其 `ΔSASA` 会变大，更可能是蛋白质相互作用界面残基。
- `ΔSASA` 默认只用于生成弱监督标签，不作为模型输入特征，避免信息泄漏。

## 方法概述

### 1. SASA 计算

`SASA` 模块基于 `Shrake-Rupley` 打点法实现，主要完成：

- 解析 `PDB` 中的原子坐标、链 ID、残基编号和残基名称。
- 读取球面采样点文件 `Dot.txt`。
- 判断每个原子的表面采样点是否被邻近原子遮挡。
- 计算原子级 `SASA`。
- 将原子级结果聚合到残基级和链级。

### 2. ΔSASA 标签生成

对于复合物中的目标链：

- 仅保留目标链，计算 `apo` 状态残基层 `SASA`。
- 保留目标链及其配对链，计算 `holo` 状态残基层 `SASA`。
- 逐残基计算 `ΔSASA`。
- 在多个阈值下生成界面残基二分类标签。

当前常用阈值包括：

- `0.5`
- `1.0`
- `2.0`
- `5.0`

### 3. 多模态界面预测

模块 3 在前两部分基础上加入 ESM-2 蛋白语言模型特征，形成残基层多模态训练数据：

- `sasa_apo`
- `sasa_holo`
- 残基坐标
- ESM-2 embedding
- 由 `ΔSASA` 规则生成的弱监督标签

下游模型支持：

- `MLP`：普通多层感知机，用于快速 baseline 和特征消融。
- `GCN`：基于残基空间距离构图的轻量图卷积网络。
- `EGNN`：使用坐标更新和距离消息传递的 E(n)-equivariant 图网络。
- `cross_egnn`：在 EGNN 基础上增加 partner-chain 距离注意力。

支持的特征组合：

- `sasa`：只使用 `sasa_apo + sasa_holo`
- `esm`：只使用 ESM-2 embedding
- `esm_sasa`：使用 ESM-2 embedding + SASA 特征
- `esm_sasa_struct`：使用 ESM-2 embedding + SASA + 二面角 + HSE + 疏水性特征

## 仓库结构

```text
README.md
requirements.txt
data/
  raw/
    examples/
    pdb_complexes/
  processed/
    examples/
    interface_labels_per_complex/
    threshold_stats_per_complex/
    complex_manifest.csv
    interface_labels_all.csv
    ml_residue_dataset.csv
    esm_residue_embeddings*.csv
    multimodal_residue_dataset*.csv
    threshold_statistics_by_complex.csv
    threshold_statistics_overall.csv
docs/
  README.md
  module_1_sasa.md
  module_2_delta_sasa.md
  module_3.md
  pipeline_overview.md
  project_spec.md
  故事线.md
  项目二.md
  自研SASA与FreeSASA对照及DeltaSASA方案.md
  assets/
scripts/
  run_*.py
src/
  sasa_project/
    *.py
```

## 目录说明

| 路径 | 作用 |
|---|---|
| `.git/` | Git 版本管理目录，记录提交历史、分支和暂存区。 |
| `.gitignore` | 指定不需要纳入 Git 跟踪的文件。 |
| `data/` | 数据目录，包含原始 PDB、示例文件和处理后的 CSV。 |
| `data/raw/` | 原始输入数据。 |
| `data/raw/examples/` | 示例 PDB、`Dot.txt` 球面采样点和示意图片。 |
| `data/raw/pdb_complexes/` | 批量复合物 PDB 数据集。 |
| `data/processed/` | 程序生成的中间结果、轻量结果摘要和最终训练数据。大文件默认不提交。 |
| `data/processed/examples/` | 单结构 `SASA` 和单复合物 `ΔSASA` 的示例输出。 |
| `data/processed/interface_labels_per_complex/` | 每个复合物单独的残基界面标签 CSV。本地可再生缓存，不纳入 Git。 |
| `data/processed/threshold_stats_per_complex/` | 每个复合物在不同阈值下的标签统计。本地可再生缓存，不纳入 Git。 |
| `docs/` | 项目文档、模块说明、实验分析和报告材料。 |
| `docs/assets/` | 文档中使用的图片和实验图。 |
| `scripts/` | 命令行运行入口。大多数脚本只负责调用 `src/sasa_project/` 中的主逻辑。 |
| `src/sasa_project/` | 项目核心源码。 |
| `__pycache__/` | Python 自动生成的字节码缓存，可忽略。 |

## 核心数据文件

| 文件 | 作用 |
|---|---|
| `data/processed/complex_manifest.csv` | 复合物清单，记录 PDB ID、目标链、配对链和文件路径等信息。 |
| `data/processed/interface_labels_all.csv` | 所有复合物合并后的残基层界面标签总表。 |
| `data/processed/esm_residue_embeddings_650m.csv` | ESM-2 650M 模型提取的残基层 embedding。 |
| `data/processed/multimodal_residue_dataset_650m.csv` | 正式多模态训练表，包含 1280 维 ESM 和 9 个 SASA / 结构特征。 |
| `data/processed/benchmark_dset186_*` | Dset_186-local 标签、manifest、预测和指标。 |
| `data/processed/benchmark_pdbtest315_*` | PDBtest_315-local 标签、manifest、预测和指标。 |
| `data/processed/artifact_manifest.csv` | 本地大文件的大小、用途和 Git 交付策略。 |
| `data/processed/threshold_statistics_by_complex.csv` | 每个复合物的阈值统计汇总。 |
| `data/processed/threshold_statistics_overall.csv` | 全数据集整体阈值统计。 |

## 源码文件说明

| 文件 | 作用 |
|---|---|
| `src/sasa_project/__init__.py` | 将 `sasa_project` 标记为 Python 包。 |
| `src/sasa_project/paths.py` | 管理项目根目录、数据目录和示例目录等路径。 |
| `src/sasa_project/sasa.py` | 模块 1 核心。负责 PDB 解析、球面打点、原子级 `SASA` 计算、残基级和链级汇总。 |
| `src/sasa_project/delta_sasa_label.py` | 模块 2 核心。负责单个复合物的 `apo / holo SASA` 对比、`ΔSASA` 计算和界面标签生成。 |
| `src/sasa_project/batch_generate_interface_labels.py` | 批量处理复合物清单，生成所有复合物的界面标签和阈值统计。 |
| `src/sasa_project/collect_complex_dataset.py` | 从 RCSB / PDB 查询、下载或整理复合物结构，并生成 `complex_manifest.csv`。 |
| `src/sasa_project/prepare_mlp_dataset.py` | 将界面标签整理成早期下游分类模型可用的残基层训练表。 |
| `src/sasa_project/residue_features.py` | 残基层公共工具，包括三字母氨基酸转一字母序列、残基提取、`sample_id` 生成、残基坐标提取和空间距离构图。 |
| `src/sasa_project/extract_esm_embeddings.py` | 调用 ESM-2 模型，为目标链每个残基提取语言模型 embedding。 |
| `src/sasa_project/build_multimodal_dataset.py` | 合并 `ΔSASA` 标签、SASA 特征、残基坐标和 ESM embedding，生成多模态训练表。 |
| `src/sasa_project/train_interface_model.py` | 训练或评测 `MLP / GCN / EGNN / cross_egnn`，并输出 Accuracy、Precision、Recall、F1、AUROC、AUPRC。 |

## 脚本入口说明

| 脚本 | 作用 |
|---|---|
| `scripts/run_sasa_example.py` | 运行单个示例 PDB 的 `SASA` 计算。 |
| `scripts/run_delta_sasa_example.py` | 运行单个复合物的 `ΔSASA` 标签生成示例。 |
| `scripts/run_collect_dataset.py` | 收集或下载复合物数据。 |
| `scripts/run_batch_labeling.py` | 批量生成所有复合物的界面标签。 |
| `scripts/run_prepare_mlp_dataset.py` | 生成早期 MLP 训练表。 |
| `scripts/run_extract_esm_embeddings.py` | 提取 ESM-2 残基 embedding。 |
| `scripts/run_build_multimodal_dataset.py` | 构建 SASA + ESM + 标签的多模态训练数据集。 |
| `scripts/run_train_interface_model.py` | 训练或评测 MLP / GCN / EGNN / cross-chain EGNN。 |
| `scripts/run_benchmark_eval.py` | 下载并生成 Dset_186-local 或 PDBtest_315-local 的 ΔSASA 标签与 manifest。 |
| `scripts/run_freesasa_comparison.py` | 在单个示例上对比自研 SASA 与 FreeSASA。 |
| `scripts/run_batch_freesasa_comparison.py` | 批量对比自研 SASA 与 FreeSASA，并生成验证图。 |

## 环境安装

建议在仓库根目录创建虚拟环境后安装依赖：

```bash
pip install -r requirements.txt
```

如果需要运行 ESM-2 或 GPU 训练，还需要安装兼容 CUDA 的 PyTorch，并安装 `transformers`：

```bash
pip install transformers
pip install --upgrade torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu128
```

Windows PowerShell 下运行脚本前，可以设置：

```powershell
$env:PYTHONPATH="src"
```

评测脚本已经兼容 Windows 默认 GBK 控制台，不再要求额外设置
`PYTHONUTF8=1`。设置该变量仍然可用。

Linux / macOS / Git Bash 下可以使用：

```bash
export PYTHONPATH=src
```

## 使用方式

### 1. 运行 SASA 示例

```bash
PYTHONPATH=src python scripts/run_sasa_example.py
```

Windows PowerShell：

```powershell
$env:PYTHONPATH="src"
python scripts/run_sasa_example.py
```

### 2. 运行单个复合物的 ΔSASA 标签示例

```bash
PYTHONPATH=src python scripts/run_delta_sasa_example.py --target-chain C --partner-chains D
```

### 3. 批量生成复合物界面标签

```bash
PYTHONPATH=src python scripts/run_batch_labeling.py
```

### 4. 生成早期 MLP 训练主表

```bash
PYTHONPATH=src python scripts/run_prepare_mlp_dataset.py --default-threshold 2.0
```

### 5. 提取 ESM-2 残基 embedding

```bash
PYTHONPATH=src python scripts/run_extract_esm_embeddings.py \
  --model-name facebook/esm2_t33_650M_UR50D \
  --device cuda \
  --output data/processed/esm_residue_embeddings_650m.csv
```

### 6. 构建多模态训练数据集

```bash
PYTHONPATH=src python scripts/run_build_multimodal_dataset.py \
  --embeddings data/processed/esm_residue_embeddings_650m.csv \
  --label-threshold 2.0 \
  --output data/processed/multimodal_residue_dataset_650m.csv
```

### 7. 训练界面残基预测模型

SASA-only baseline：

```bash
PYTHONPATH=src python scripts/run_train_interface_model.py \
  --input data/processed/multimodal_residue_dataset_650m.csv \
  --model mlp \
  --feature-set sasa \
  --epochs 30 \
  --device cuda
```

ESM-only baseline：

```bash
PYTHONPATH=src python scripts/run_train_interface_model.py \
  --input data/processed/multimodal_residue_dataset_650m.csv \
  --model mlp \
  --feature-set esm \
  --epochs 30 \
  --device cuda
```

ESM + SASA：

```bash
PYTHONPATH=src python scripts/run_train_interface_model.py \
  --input data/processed/multimodal_residue_dataset_650m.csv \
  --model mlp \
  --feature-set esm_sasa \
  --epochs 30 \
  --device cuda
```

GCN + ESM + SASA：

```bash
PYTHONPATH=src python scripts/run_train_interface_model.py \
  --input data/processed/multimodal_residue_dataset_650m.csv \
  --model gcn \
  --feature-set esm_sasa \
  --epochs 30 \
  --device cuda
```

正式 EGNN：

```bash
PYTHONPATH=src python scripts/run_train_interface_model.py \
  --input data/processed/multimodal_residue_dataset_650m.csv \
  --model egnn \
  --feature-set esm_sasa_struct \
  --device cuda
```

### 8. 运行标准 benchmark

完整 PDBtest_315-local 流程：

```bash
python scripts/run_benchmark_eval.py --benchmark pdbtest315

PYTHONPATH=src python -m sasa_project.extract_esm_embeddings \
  --manifest data/processed/benchmark_pdbtest315_manifest.csv \
  --model-name facebook/esm2_t33_650M_UR50D \
  --device cuda \
  --output data/processed/benchmark_pdbtest315_esm_embeddings_650m.csv

PYTHONPATH=src python -m sasa_project.build_multimodal_dataset \
  --labels data/processed/benchmark_pdbtest315_labels.csv \
  --manifest data/processed/benchmark_pdbtest315_manifest.csv \
  --embeddings data/processed/benchmark_pdbtest315_esm_embeddings_650m.csv \
  --output data/processed/benchmark_pdbtest315_multimodal_650m.csv
```

模型推断命令见 [data/processed/README.md](data/processed/README.md)。

## 当前实验结果摘要

正式训练使用 `500` 个复合物、`95,005` 个残基样本，默认标签为：

```text
label = 1 if ΔSASA > 2.0 else 0
```

主模型结果如下：

| Setting | Feature dim | Test F1 | Test AUROC | Test AUPRC |
|---|---:|---:|---:|---:|
| EGNN + ESM-2 650M + SASA + structure | 1289 | 0.8948 | 0.9768 | 0.9497 |
| Cross-chain EGNN + ESM-2 650M + SASA + structure | 1289 | 0.8938 | 0.9755 | 0.9485 |

外部 benchmark 结果：

| Dataset | Complexes | Residues | Model | F1 | AUROC | AUPRC |
|---|---:|---:|---|---:|---:|---:|
| Dset_186-local | 158 | 39,505 | EGNN | 0.6448 | 0.8868 | 0.6770 |
| Dset_186-local | 158 | 39,505 | Cross-chain EGNN | 0.6439 | 0.8891 | 0.5766 |
| PDBtest_315-local | 314 | 65,119 | EGNN | 0.6761 | 0.8946 | 0.7408 |
| PDBtest_315-local | 314 | 65,119 | Cross-chain EGNN | 0.6816 | 0.8983 | 0.6906 |

`-local` 表示按当前仓库的结构可用性和链质量规则处理后的本地子集。

更完整的实验说明见 [docs/module_3.md](docs/module_3.md)。

## 文档说明

| 文档 | 作用 |
|---|---|
| `docs/README.md` | 文档导航。 |
| `docs/module_1_sasa.md` | 模块 1：自研 `SASA` 计算说明。 |
| `docs/module_2_delta_sasa.md` | 模块 2：`ΔSASA` 标签生成说明。 |
| `docs/module_3.md` | 模块 3：ESM-2 + SASA 多模态界面预测实验。 |
| `docs/pipeline_overview.md` | 项目整体流程和思路说明。 |
| `docs/project_spec.md` | 项目规格和设计说明。 |
| `docs/故事线.md` | 项目汇报或论文叙事逻辑。 |
| `docs/自研SASA与FreeSASA对照及DeltaSASA方案.md` | 自研 SASA、FreeSASA 对照与 ΔSASA 方案说明。 |

## 当前状态

目前仓库已经完成：

- 自研 `SASA` 计算器。
- 自研 SASA 与 FreeSASA 的示例和批量对照。
- 复合物数据收集与筛选。
- 残基层 `apo / holo SASA` 计算。
- `ΔSASA` 界面弱监督标签生成。
- 多阈值标签统计。
- ESM-2 残基层 embedding 提取。
- SASA / ESM 多模态训练数据构建。
- MLP / GCN / EGNN / cross-chain EGNN 下游界面残基预测实验。
- 特征消融与 `ΔSASA` 阈值敏感性实验。
- Dset_186-local 和 PDBtest_315-local 全量 GPU 评测。
- 自动化 smoke test。

## 后续改进方向

- 引入独立真实标注测试集，降低弱监督标签带来的评价偏差。
- 尝试更大的 ESM-2 模型或其他蛋白语言模型。
- 扩展至 1,000+ 训练复合物。
- 尝试更强的跨链图结构或 Graph Transformer。
- 与 GraphPPIS 等已有方法的数据划分和 baseline 做更直接对比。
