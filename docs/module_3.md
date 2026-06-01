# 模块 3：ESM-2 + SASA 多模态界面预测实验

本模块在前两部分 `SASA` 和 `ΔSASA` 的基础上，进一步完成 B + D 思路中的下游机器学习闭环：

```text
自研 SASA 计算
-> apo / holo 残基层 SASA
-> ΔSASA 弱监督标签
-> ESM-2 残基层序列表征
-> SASA / ESM 多模态特征融合
-> MLP / GCN 界面残基预测
-> 消融实验与阈值敏感性分析
```

## 1. 新增代码说明

### 1.1 `src/sasa_project/residue_features.py`

该文件提供残基层特征处理的公共工具，主要功能包括：

- 将三字母氨基酸名转换为一字母序列，用于输入 ESM-2。
- 从 PDB 原子列表中提取目标链残基列表。
- 为每个残基生成稳定的 `sample_id`，用于对齐标签、SASA 特征和 ESM 特征。
- 提取残基坐标，优先使用 `CA` 原子坐标；如果缺失 `CA`，则使用该残基全部原子的平均坐标。
- 根据残基空间距离构建图边，用于 GCN。

其中 GCN 默认使用残基中心距离阈值构图：

```text
edge(i, j) = 1 if distance(i, j) <= 8 Å
```

### 1.2 `src/sasa_project/extract_esm_embeddings.py`

该文件负责提取目标链残基层 ESM-2 embedding。

当前使用的大模型为：

```text
facebook/esm2_t33_650M_UR50D
```

该模型每个残基输出 `1280` 维 embedding。

脚本流程：

```text
读取 complex_manifest.csv
-> 解析每个 PDB
-> 提取目标链序列
-> 输入 ESM-2
-> 输出每个残基的 esm_0 ... esm_1279
```

默认输出文件：

```text
data/processed/esm_residue_embeddings.csv
```

本项目同时保留了 8M 小模型的历史结果，但正式实验使用 650M 模型：

```text
data/processed/esm_residue_embeddings_650m.csv
```

### 1.3 `src/sasa_project/build_multimodal_dataset.py`

该文件将不同来源的残基层信息合并成一个多模态训练表：

- `ΔSASA` 标签表：`interface_labels_all.csv`
- ESM-2 embedding：正式实验使用 `esm_residue_embeddings_650m.csv`
- 残基坐标：从 PDB 重新解析得到
- SASA 特征：`sasa_apo`、`sasa_holo`

输出文件：

```text
data/processed/multimodal_residue_dataset_650m.csv
```

需要特别注意：`delta_sasa` 只用于生成弱监督标签，不作为默认模型输入特征，避免数据泄漏。

### 1.4 `src/sasa_project/train_interface_model.py`

该文件负责训练下游界面残基预测模型，支持两类模型：

- `mlp`：普通多层感知机，用于快速 baseline 和特征消融。
- `gcn`：轻量级图卷积网络，用于结构感知节点分类。

支持三种特征组合：

```text
sasa      = sasa_apo + sasa_holo
esm       = ESM-2 embedding
esm_sasa  = ESM-2 embedding + sasa_apo + sasa_holo
```

训练脚本会按复合物划分训练集、验证集和测试集，避免同一复合物的残基同时出现在不同数据划分中。

评价指标包括：

- Accuracy
- Precision
- Recall
- F1
- AUROC
- AUPRC

## 2. CUDA 环境修复

原环境中 GPU 可以被 `nvidia-smi` 识别，但 PyTorch 版本不支持 RTX 5060 的 `sm_120` 架构。

原版本：

```text
torch 2.7.1+cu118
```

问题表现：

```text
NVIDIA GeForce RTX 5060 Laptop GPU with CUDA capability sm_120
is not compatible with the current PyTorch installation.
```

修复后版本：

```text
torch       2.11.0+cu128
torchvision 0.26.0+cu128
torchaudio  2.11.0+cu128
```

验证结果：

```text
CUDA available: True
GPU: NVIDIA GeForce RTX 5060 Laptop GPU
CUDA runtime: 12.8
```

修复后，全量 ESM-2 embedding 提取可以在 GPU 上完成。

## 3. 关键运行命令

### 3.1 安装依赖

```bash
pip install transformers
pip install --upgrade torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu128
```

### 3.2 提取全量 ESM-2 embedding

```bash
PYTHONPATH=src python scripts/run_extract_esm_embeddings.py \
  --model-name facebook/esm2_t33_650M_UR50D \
  --device cuda \
  --output data/processed/esm_residue_embeddings_650m.csv
```

输出：

```text
data/processed/esm_residue_embeddings_650m.csv
```

本次运行结果：

```text
100 complexes
18050 residue embeddings
```

首次运行需要下载 650M 模型。本次在模型缓存后，全量 `100` 个复合物、`18050` 个残基的 embedding 提取约 `30` 秒完成。

### 3.3 构建多模态数据集

默认使用 `ΔSASA > 2.0` 作为标签：

```bash
PYTHONPATH=src python scripts/run_build_multimodal_dataset.py \
  --embeddings data/processed/esm_residue_embeddings_650m.csv \
  --label-threshold 2.0 \
  --output data/processed/multimodal_residue_dataset_650m.csv
```

输出：

```text
data/processed/multimodal_residue_dataset_650m.csv
```

本次运行结果：

```text
18050 multimodal residue rows
```

### 3.4 特征消融实验

SASA only：

```bash
PYTHONPATH=src python scripts/run_train_interface_model.py \
  --input data/processed/multimodal_residue_dataset_650m.csv \
  --model mlp \
  --feature-set sasa \
  --epochs 30 \
  --device cuda
```

ESM only：

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

### 3.5 ΔSASA 阈值敏感性实验

生成不同阈值的数据集：

```bash
PYTHONPATH=src python scripts/run_build_multimodal_dataset.py \
  --embeddings data/processed/esm_residue_embeddings_650m.csv \
  --label-threshold 0.5 \
  --output data/processed/multimodal_residue_dataset_t0_5.csv

PYTHONPATH=src python scripts/run_build_multimodal_dataset.py \
  --embeddings data/processed/esm_residue_embeddings_650m.csv \
  --label-threshold 1.0 \
  --output data/processed/multimodal_residue_dataset_t1_0.csv

PYTHONPATH=src python scripts/run_build_multimodal_dataset.py \
  --embeddings data/processed/esm_residue_embeddings_650m.csv \
  --label-threshold 5.0 \
  --output data/processed/multimodal_residue_dataset_t5_0.csv
```

对应训练命令示例：

```bash
PYTHONPATH=src python scripts/run_train_interface_model.py \
  --input data/processed/multimodal_residue_dataset_t5_0.csv \
  --model mlp \
  --feature-set sasa \
  --epochs 30 \
  --device cuda
```

## 4. 实验结果

### 4.1 特征消融实验

使用全量 `100` 个复合物、`18050` 个残基样本。默认标签为：

```text
label = 1 if ΔSASA > 2.0 else 0
```

实验结果如下：

| Setting | Feature dim | Test F1 | Test AUROC | Test AUPRC |
|---|---:|---:|---:|---:|
| MLP + SASA only | 2 | 0.8754 | 0.9439 | 0.9136 |
| MLP + ESM-2 650M only | 1280 | 0.6480 | 0.8736 | 0.7242 |
| MLP + ESM-2 650M + SASA | 1282 | 0.6642 | 0.8767 | 0.7389 |
| GCN + ESM-2 650M + SASA | 1282 | 0.6025 | 0.8368 | 0.6156 |

完整输出摘要：

```text
MLP + SASA only
test: accuracy=0.9551 precision=1.0000 recall=0.7785 f1=0.8754 auroc=0.9439 auprc=0.9136

MLP + ESM-2 650M only
test: accuracy=0.8430 precision=0.5934 recall=0.7136 f1=0.6480 auroc=0.8736 auprc=0.7242

MLP + ESM-2 650M + SASA
test: accuracy=0.8542 precision=0.6224 recall=0.7120 f1=0.6642 auroc=0.8767 auprc=0.7389

GCN + ESM-2 650M + SASA
test: accuracy=0.8068 precision=0.5164 recall=0.7231 f1=0.6025 auroc=0.8368 auprc=0.6156
```

### 4.2 结果分析

从结果可以看到：

1. `ESM + SASA` 相比 `ESM only` 有提升：

```text
ESM-2 650M only:      F1=0.6480, AUROC=0.8736, AUPRC=0.7242
ESM-2 650M + SASA:    F1=0.6642, AUROC=0.8767, AUPRC=0.7389
```

这说明几何暴露信息可以补充 ESM-2 的序列语义信息，多模态融合是有意义的。

2. `SASA only` 表现最强：

```text
SASA only: F1=0.8754, AUROC=0.9439, AUPRC=0.9136
```

原因是本项目的弱监督标签由 `ΔSASA` 规则生成，而 `sasa_apo`、`sasa_holo` 与标签生成机制天然高度相关。因此，SASA-only baseline 在当前实验中非常强。

这并不否定多模态框架的价值，而是说明当前标签体系下，几何暴露特征本身就是最直接的监督信号来源。更严格的泛化验证需要引入独立真实标注测试集。

3. 当前轻量 GCN 没有超过 MLP：

GCN 使用简单距离图和两层图卷积，尚未引入边距离权重、attention、等变 GNN 或 Graph Transformer。因此它更适合作为结构感知 baseline，而不是最终 SOTA 模型。

## 5. ΔSASA 阈值敏感性实验

本实验考察不同 `ΔSASA` 阈值对伪标签分布和模型性能的影响。

### 5.1 标签分布

| ΔSASA threshold | Total | Positive | Negative | Positive ratio |
|---:|---:|---:|---:|---:|
| 0.5 | 18050 | 4162 | 13888 | 0.2306 |
| 1.0 | 18050 | 4162 | 13888 | 0.2306 |
| 2.0 | 18050 | 3981 | 14069 | 0.2206 |
| 5.0 | 18050 | 3684 | 14366 | 0.2041 |

其中 `0.5` 和 `1.0` 的结果完全相同，说明当前数据集中没有残基的 `ΔSASA` 落在 `(0.5, 1.0]` 这个区间内。

### 5.2 模型结果

使用 `MLP + SASA only`，训练 `30` epoch。

| ΔSASA threshold | Test F1 | Test AUROC | Test AUPRC |
|---:|---:|---:|---:|
| 0.5 | 0.8710 | 0.9384 | 0.9066 |
| 1.0 | 0.8710 | 0.9384 | 0.9066 |
| 2.0 | 0.8754 | 0.9439 | 0.9136 |
| 5.0 | 0.9026 | 0.9643 | 0.9395 |

### 5.3 阈值实验分析

随着阈值从 `0.5 / 1.0` 提高到 `5.0`：

- 正样本比例从 `0.2306` 降低到 `0.2041`。
- 标签变得更严格，只保留更明显的界面残基。
- 在当前数据和 SASA-only 模型下，`5.0` 阈值取得了最高的 F1、AUROC 和 AUPRC。

这说明更严格的 `ΔSASA` 阈值可能降低伪标签噪声，使模型更容易学习到高置信度界面残基。但过高阈值也可能遗漏弱接触位点，因此最终阈值选择需要结合任务目标和真实标注数据进一步验证。

## 6. 当前结论

本模块完成了从几何计算到机器学习预测的工程闭环：

```text
SASA 计算器
-> ΔSASA 弱监督伪标签
-> ESM-2 残基层 embedding
-> SASA / ESM 多模态特征融合
-> MLP / GCN 节点二分类
-> 消融实验与阈值敏感性实验
```

可以在报告中概括为：

> 本项目并未停留在 SASA 数值计算本身，而是进一步将自研 SASA 计算器作为结构生物学特征提取与弱监督标签生成工具，结合 ESM-2 蛋白语言模型和结构图学习方法，构建了一个用于蛋白质相互作用界面残基预测的多模态弱监督框架。

## 7. 局限性与后续改进

当前版本仍有以下局限：

1. ESM-2 已升级为 `650M` 模型，embedding 维度为 `1280`，但仍未使用更大的 `3B / 15B` 模型。
2. 标签来自 `ΔSASA` 规则，是弱监督伪标签，不等同于人工标注或高质量数据库真值。
3. `SASA only` 表现很强，部分原因是输入特征与标签生成规则高度相关。
4. GCN 是轻量实现，尚未加入边距离权重、GAT、Graph Transformer 或 E(3)-equivariant GNN。
5. 目前实验主要验证框架可行性，若要进一步提升含金量，应引入独立真实标签测试集或 GraphPPIS 标准数据划分进行外部评估。

后续可以优先尝试：

- 使用更大的 ESM-2 模型。
- 引入距离权重或 attention-based GNN。
- 做 GCN 邻域半径实验：`6 Å / 8 Å / 10 Å / 12 Å`。
- 与 GraphPPIS 的数据划分和 baseline 进行更直接对比。

## 8. 参考项目与文献

- GraphPPIS: https://github.com/biomed-AI/GraphPPIS
- EquiPPIS: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011435
- DeepProSite: https://pubmed.ncbi.nlm.nih.gov/38015872/
- GTE-PPIS: https://academic.oup.com/bib/article/26/3/bbaf290/8164226
- EDG-PPIS: https://pubmed.ncbi.nlm.nih.gov/41023600/
- ASCE-PPIS: https://academic.oup.com/bioinformatics/article/41/8/btaf423/8211827
- MPBind: https://academic.oup.com/bioinformatics/article/doi/10.1093/bioinformatics/btaf589/8300842

## 问题
- 阈值呈现效果，GCN效果
- 整理更多创新点:多链复合体
- 近几年论文的baseline和数据集:已完成补充500条双链复合物，并且增加Dset_186评测集