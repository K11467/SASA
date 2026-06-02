# 模块 3：ESM-2 650M + Geometry + EGNN 界面预测

## 1. 严格主流程

```text
PDB 复合物
-> apo / holo SASA
-> Delta-SASA 弱监督标签
-> ESM-2 650M 残基层 embedding
-> ESM + 坐标 + 二面角 + HSE + 疏水性
-> EGNN
-> 残基层界面概率
```

`delta_sasa` 只用于生成标签。严格主模型不输入 `sasa_holo`，也不输入任何
SASA 特征，避免模型从标签定义 `sasa_apo - sasa_holo` 中获得近邻泄漏。

## 2. 特征集

| 名称 | 输入 | 用途 |
|---|---|---|
| `apo` | `sasa_apo` | 简单 baseline |
| `sasa` | `sasa_apo + sasa_holo` | 泄漏诊断，不作为正式预测结果 |
| `esm` | ESM-2 650M embedding | MLP baseline |
| `esm_struct` | ESM + geometry | 严格主模型 |
| `esm_apo_struct` | apo SASA + ESM + geometry | apo-only 消融 |
| `esm_sasa_struct` | apo/holo SASA + ESM + geometry | 泄漏诊断 |

严格主模型共 `1287` 个标量输入：`1280` 维 ESM、`4` 个二面角编码、
`2` 个 HSE 和 `1` 个疏水性特征。坐标用于 EGNN 构图和等变更新。

## 3. 内部测试结果

| 角色 | 模型 / 特征 | F1 | AUROC | AUPRC |
|---|---|---:|---:|---:|
| Strict baseline | apo-SASA rank | 0.3317 | 0.5841 | 0.2550 |
| Strict baseline | MLP / apo | 0.4292 | 0.7000 | 0.3098 |
| Strict baseline | MLP / ESM | 0.6290 | 0.8949 | 0.7387 |
| Strict baseline | GCN / ESM + geometry | 0.5961 | 0.8730 | 0.6570 |
| Strict primary | EGNN / ESM + geometry | **0.7745** | **0.9322** | **0.8421** |
| Strict ablation | EGNN / apo + ESM + geometry | 0.7787 | 0.9348 | 0.8514 |
| Holo-aware diagnostic | MLP / apo + holo | 0.8724 | 0.9457 | 0.9117 |
| Holo-aware diagnostic | EGNN / apo + holo + ESM + geometry | 0.8948 | 0.9768 | 0.9497 |

含 `sasa_holo` 的分数明显更高，但不能作为可部署模型的主结果。

## 4. 外部 benchmark

| 数据集 | 角色 | 模型 | F1 | AUROC | AUPRC |
|---|---|---|---:|---:|---:|
| Dset_186-local | Strict primary | EGNN / ESM + geometry | 0.3529 | 0.7277 | 0.3050 |
| PDBtest_315-local | Strict primary | EGNN / ESM + geometry | 0.3040 | 0.6799 | 0.2675 |
| PDBtest_315-local | Holo-aware diagnostic | EGNN | 0.6761 | 0.8946 | 0.7408 |
| PDBtest_315-local | Holo-aware diagnostic | Cross-chain EGNN | 0.6816 | 0.8983 | 0.6906 |

严格外部结果显示仍存在明显泛化差距。cross-chain EGNN 只作为分析模块：
它在部分 holo-aware 外部结果上提高 recall 和 F1，但收益不稳定。

## 5. 复现

完整命令见 [data/processed/README.md](../data/processed/README.md)。

```powershell
$env:PYTHONPATH="src"
python -m unittest discover -s tests -v
python scripts\run_benchmark_eval.py --help
```

Windows 默认 GBK 控制台已兼容，不再要求设置 `PYTHONUTF8=1`。

## 6. 交付文件

- `data/processed/leakage_ablation_summary_650m.csv`
- `data/processed/benchmark_dset186_metrics_650m.csv`
- `data/processed/benchmark_pdbtest315_metrics_650m.csv`
- `data/processed/best_egnn_650m_esm_struct_d8.pt`
- `data/processed/artifact_manifest.csv`

GraphPPIS 等文献方法暂不填入数值横向表。当前 local manifest、标签协议与公开论文
协议不一定完全一致；正式对比需要在同一 manifest 和标签上运行双方方法。
