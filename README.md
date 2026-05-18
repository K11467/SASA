# SASA Project 2

本仓库整理了项目二当前已经完成的核心流程：

1. 单结构 `SASA` 计算
2. 基于复合物的 `ΔSASA` 界面标签生成
3. `100` 个蛋白复合物数据集整理
4. 面向后续 `MLP` 训练的残基层主表输出

## 目录结构

```text
README.md
data/
  raw/
    examples/             # 示例 PDB、Dot 点集、图片
    pdb_complexes/        # 100 个复合物 biological assembly PDB
  processed/
    examples/             # 示例运行结果
    interface_labels_per_complex/
    threshold_stats_per_complex/
    complex_manifest.csv
    interface_labels_all.csv
    ml_residue_dataset.csv
    threshold_statistics_by_complex.csv
    threshold_statistics_overall.csv
docs/
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

## 运行方式

推荐在仓库根目录执行，并通过 `PYTHONPATH=src` 让脚本找到包：

```bash
PYTHONPATH=src python3 scripts/run_sasa_example.py
PYTHONPATH=src python3 scripts/run_delta_sasa_example.py --target-chain C --partner-chains D
PYTHONPATH=src python3 scripts/run_collect_dataset.py --count 100
PYTHONPATH=src python3 scripts/run_batch_labeling.py
PYTHONPATH=src python3 scripts/run_prepare_mlp_dataset.py --default-threshold 2.0
```

## 当前关键产物

- 复合物清单：`data/processed/complex_manifest.csv`
- 全部残基标签总表：`data/processed/interface_labels_all.csv`
- 总体阈值统计：`data/processed/threshold_statistics_overall.csv`
- 给后续模型的训练主表：`data/processed/ml_residue_dataset.csv`

## 说明文档

- 第一部分说明：`docs/module_1_sasa.md`
- 第二部分说明：`docs/module_2_delta_sasa.md`
- 项目整体思路：`docs/pipeline_overview.md`
- 原始任务说明：`docs/project_spec.md`
