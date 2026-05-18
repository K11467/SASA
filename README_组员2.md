# ΔSASA 标签生成说明

本模块负责基于复合物结构，为目标链残基生成界面标签。

脚本文件：
- `delta_sasa_label.py`

依赖基础模块：
- `sasa.py`

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

- `interface_labels.csv`
  - 每个目标残基的 `SASA_apo`、`SASA_holo`、`ΔSASA`
  - 每个阈值对应的二分类标签

- `threshold_statistics.csv`
  - 每个阈值下的正负样本数量
  - 每个阈值下的正负样本比例

## 运行方式

```bash
python3 delta_sasa_label.py \
  --pdb your_complex.pdb \
  --target-chain A \
  --partner-chains B
```

如果不写 `--partner-chains`，脚本会自动把除目标链外的其他链当作配对链。

## 当前数据说明

当前仓库中的 `2iww_H.pdb` 只有 `A` 链，没有复合物配对链，因此**不能直接用于真实的 ΔSASA 标签生成**。脚本在这种情况下会报错并提示当前 PDB 不适合做 holo 计算。
