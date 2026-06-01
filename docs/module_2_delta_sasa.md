# ΔSASA 标签生成说明

本模块负责基于复合物结构，为目标链残基生成界面标签。

脚本文件：
- `src/sasa_project/delta_sasa_label.py`

依赖基础模块：
- `src/sasa_project/sasa.py`

## 计算逻辑

1. 读取复合物 PDB。
2. 按链拆分原子。
3. 对目标链单独计算 apo 状态残基 SASA。
4. 对目标链和配对链一起计算 holo 状态 SASA。
5. 只统计目标链残基：

```text
ΔSASA = SASA_apo - SASA_holo
```

6. 用多个阈值自动生成界面标签：
- `0.5`
- `1.0`
- `2.0`
- `5.0`

## 输出文件

- `data/processed/examples/interface_labels.csv`
  - 每个目标残基的 `SASA_apo`、`SASA_holo`、`ΔSASA`
  - 每个阈值对应的二分类标签

- `data/processed/examples/threshold_statistics.csv`
  - 每个阈值下的正负样本数量
  - 每个阈值下的正负样本比例

## 运行方式

```bash
PYTHONPATH=src python3 scripts/run_delta_sasa_example.py \
  --pdb your_complex.pdb \
  --target-chain A \
  --partner-chains B
```

如果不写 `--partner-chains`，脚本会自动把除目标链外的其他链当作配对链。

## 当前数据说明

当前仓库中的 `data/raw/examples/2iww_H.pdb` 只有 `A` 链，没有复合物配对链，因此**不能直接用于真实的 ΔSASA 标签生成**。脚本在这种情况下会报错并提示当前 PDB 不适合做 holo 计算。

为完成这一部分的真实标签构造，仓库中已补充官方复合物结构：

- `data/raw/examples/2WWM.pdb`

该结构来自 RCSB PDB，条目 `2WWM`，是一个官方标注的 `Hetero 2-mer` 蛋白复合物，适合做界面残基标签生成。

当前已跑通两组链对：

- `C-D`
- `O-T`

其中默认交付结果使用 `C-D` 这组链对，并生成：

- `data/processed/examples/interface_labels.csv`
- `data/processed/examples/threshold_statistics.csv`

同时保留更明确命名的原始输出：

- `data/processed/examples/2WWM_CD_interface_labels.csv`
- `data/processed/examples/2WWM_CD_threshold_statistics.csv`
- `data/processed/examples/2WWM_OT_interface_labels.csv`
- `data/processed/examples/2WWM_OT_threshold_statistics.csv`

## 批量复合物数据集

为了支持后续 `MLP` 训练，仓库中已补充一批可直接用于界面标签构造的复合物 biological assembly 数据：

- 数据目录：`data/raw/pdb_complexes/`
- 清单文件：`data/processed/complex_manifest.csv`
- 当前规模：`100` 个复合物

筛选原则：

- 来源为 `RCSB PDB`
- biological assembly 中蛋白链数为 `2`
- `X-RAY DIFFRACTION`
- 分辨率不高于 `3.0 Å`
- 本地再次过滤，只保留每条链残基数都不少于 `20` 的样本

清单中包含：

- `pdb_id`
- `assembly_id`
- `target_chain`
- `partner_chain`
- 两条链的残基数
- 两条链的原子数
- 本地 `pdb` 文件路径

## 批量标签结果

基于这 `100` 个复合物，仓库中已进一步生成批量 `ΔSASA` 标签结果：

- 单个复合物标签目录：`data/processed/interface_labels_per_complex/`
- 单个复合物阈值统计目录：`data/processed/threshold_stats_per_complex/`
- 全部残基汇总表：`data/processed/interface_labels_all.csv`
- 按复合物汇总的阈值统计：`data/processed/threshold_statistics_by_complex.csv`
- 全数据集总体阈值统计：`data/processed/threshold_statistics_overall.csv`

当前汇总表规模：

- `100` 个复合物
- `18050` 个目标链残基样本

## 给后续 MLP 的直接输入

为了方便第三位同学直接对接建模，仓库中进一步整理了一个更适合机器学习读取的残基层主表：

- 输出文件：`data/processed/ml_residue_dataset.csv`
- 生成脚本：`src/sasa_project/prepare_mlp_dataset.py`

该表保留了：

- `pdb_id`
- `target_chain`
- `partner_chain`
- `chain_id`
- `residue_id`
- `insertion_code`
- `residue_name`
- `sasa_apo`
- `sasa_holo`
- `delta_sasa`
- 四个阈值下的标签列

同时额外提供：

- `sample_id`：每个残基样本的唯一标识
- `label`：默认用于训练的二分类标签
- `label_threshold`：当前默认标签对应的阈值

默认使用 `ΔSASA > 2.0` 作为训练标签，但如果后续需要切换到 `0.5 / 1.0 / 5.0`，可以直接重新运行：

```bash
PYTHONPATH=src python3 scripts/run_prepare_mlp_dataset.py --default-threshold 1.0
```
