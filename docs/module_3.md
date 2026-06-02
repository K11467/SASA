# 模块 3：ESM-2 650M + SASA + EGNN 界面预测

## 1. 当前主流程

```text
PDB 复合物
-> apo / holo SASA
-> Delta-SASA 弱监督标签
-> ESM-2 650M 残基层 embedding
-> SASA + 坐标 + 二面角 + HSE + 疏水性特征
-> EGNN / cross-chain EGNN
-> 残基层界面概率
```

`delta_sasa` 仅用于生成标签，不进入模型特征，避免直接信息泄漏。

## 2. 正式数据

| 数据集 | 复合物 | 残基 | 用途 |
|---|---:|---:|---|
| Main corpus | 500 | 95,005 | 训练、验证和内部测试 |
| Dset_186-local | 158 | 39,505 | 外部 benchmark |
| PDBtest_315-local | 314 | 65,119 | 外部 benchmark |

`-local` 表示结构下载成功后，继续经过当前链质量规则得到的本地可分析子集。

正式 ESM 模型为：

```text
facebook/esm2_t33_650M_UR50D
```

每个残基包含：

| 特征 | 维度 |
|---|---:|
| ESM-2 650M embedding | 1280 |
| `sasa_apo`、`sasa_holo` | 2 |
| `sin_phi`、`cos_phi`、`sin_psi`、`cos_psi` | 4 |
| `hse_up`、`hse_dn` | 2 |
| `hydrophobicity` | 1 |
| 合计 | 1289 |

## 3. 模型

`src/sasa_project/train_interface_model.py` 支持：

- `mlp`：残基层 MLP baseline。
- `gcn`：距离图轻量 GCN baseline。
- `egnn`：E(n)-equivariant residue graph network。
- `cross_egnn`：加入 partner-chain 距离注意力的 EGNN。

EGNN 默认使用 3 层、hidden dim 128、dropout 0.4、8 A 图 cutoff。训练支持梯度累积、梯度裁剪和基于 validation F1 的 early stopping。

## 4. CUDA 环境

本机验证环境：

```text
GPU: NVIDIA GeForce RTX 5060 Laptop GPU
torch: 2.11.0+cu128
CUDA available: True
CUDA runtime: 12.8
```

Windows PowerShell：

```powershell
$env:PYTHONPATH="src"
```

`scripts/run_benchmark_eval.py --help` 已兼容默认 GBK 控制台，不再依赖
`PYTHONUTF8=1`。

## 5. 主训练结果

| 模型 | Cutoff | Accuracy | Precision | Recall | F1 | AUROC | AUPRC |
|---|---:|---:|---:|---:|---:|---:|---:|
| EGNN | 8 A | 0.9590 | 0.9023 | 0.8874 | **0.8948** | **0.9768** | **0.9497** |
| Cross-chain EGNN | 8 A | 0.9584 | 0.8976 | 0.8899 | 0.8938 | 0.9755 | 0.9485 |

EGNN cutoff 消融：

| Cutoff | F1 | AUROC | AUPRC |
|---:|---:|---:|---:|
| 6 A | 0.8902 | 0.9750 | 0.9420 |
| 8 A | **0.8948** | **0.9768** | **0.9497** |
| 10 A | 0.8907 | 0.9758 | 0.9488 |
| 12 A | 0.8908 | 0.9727 | 0.9449 |

## 6. 外部 benchmark

| 数据集 | 模型 | Accuracy | Precision | Recall | F1 | AUROC | AUPRC |
|---|---|---:|---:|---:|---:|---:|---:|
| Dset_186-local | EGNN | 0.9184 | 0.6690 | 0.6223 | **0.6448** | 0.8868 | **0.6770** |
| Dset_186-local | Cross-chain EGNN | 0.9113 | 0.6167 | 0.6735 | 0.6439 | **0.8891** | 0.5766 |
| PDBtest_315-local | EGNN | 0.9136 | 0.7360 | 0.6252 | 0.6761 | 0.8946 | **0.7408** |
| PDBtest_315-local | Cross-chain EGNN | 0.9107 | 0.7014 | 0.6629 | **0.6816** | **0.8983** | 0.6906 |

PDBtest_315-local 已在单卡 RTX 5060 Laptop GPU 上完整运行：

- 315 个官方 chain-aware 条目均取得结构。
- 314 个条目通过当前链质量筛选。
- 65,119 个残基完成 ESM-2 650M embedding。
- EGNN 与 cross-chain EGNN 均完成全量推断。

## 7. 运行命令

完整命令见 [data/processed/README.md](../data/processed/README.md)。

基础验证：

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
python scripts/run_benchmark_eval.py --help
```

## 8. 交付策略

Git 提交轻量的源码、测试、manifest、标签、预测、指标和论文表格。
接近 GB 级的 embedding、多模态表和 checkpoint 保持忽略状态，通过命令重新生成或通过外部制品存储分发。文件清单见：

```text
data/processed/artifact_manifest.csv
```

## 9. 局限性

- 标签来自 Delta-SASA 规则，不等同于人工实验真值。
- 输入包含 SASA 与坐标，因此内部测试指标应视为弱监督场景下的结果。
- 外部 benchmark 是按当前管线筛选后的 local 子集。
- 仍需扩展至 1,000+ 训练复合物，并与更多公开方法统一协议对比。
