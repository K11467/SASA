# 自研 SASA 计算器与 FreeSASA 对照及后续 ΔSASA 标签生成方案

## 一、总体建议

建议采用“双轨策略”：

> **自研 SASA 计算器用于完成 Project 2 的核心算法实现与原理验证；FreeSASA 用作成熟工具 baseline，并在后续批量 ΔSASA 标签生成中作为主计算器或校验器。**

这样做既不会削弱“自研 SASA”的贡献，反而会让整个项目更可信、更专业。

---

## 二、为什么后续可以用 FreeSASA？

后续的 ΔSASA 标签质量会直接影响 MLP / GCN / EGNN 的训练效果。

你们自研 SASA 计算器适合体现算法理解和课程项目实现，但它可能存在一些工程误差来源，例如：

```text
1. 原子半径表是否完整
2. PDB 特殊原子名是否正确解析
3. 氢原子是否处理一致
4. 非标准残基是否跳过
5. 多链之间的遮挡判断是否严谨
6. 球面采样点数量是否足够
7. 残基编号和插入码是否正确对齐
```

FreeSASA 是比较成熟的开源 SASA 工具，支持命令行、C API 和 Python 接口，并实现了 Lee-Richards 和 Shrake-Rupley 两类经典 SASA 近似算法。因此，在批量处理几十个甚至上百个 PDB 复合物时，用 FreeSASA 会更稳。

---

## 三、这样会不会让自研 SASA 显得没用？

不会。

你们要把贡献定位讲清楚：

> **自研 SASA 的作用不是取代 FreeSASA，而是完成算法理解、实现验证和项目二原始要求。**

可以在报告中这样表述：

```text
我们首先基于 Shrake-Rupley 打点法实现了自研 SASA 计算器，用于完成 Project 2 的核心算法要求，并输出原子级、残基级和链级 SASA。随后，我们将自研结果与成熟工具 FreeSASA 进行对照，以验证实现的合理性。在后续大规模 ΔSASA 伪标签生成阶段，为提高标签稳定性和可复现性，我们采用 FreeSASA 作为主计算工具，同时保留自研工具作为算法验证与小规模复核模块。
```

这个逻辑体现了三层价值：

```text
1. 你们确实实现了 SASA 算法；
2. 你们知道用成熟工具做对照验证；
3. 你们为了下游 AI 训练质量，选择了更稳定的工程工具。
```

这比“全程只用自己写的代码”更像一个成熟项目。

---

## 四、推荐的最终方案

### 阶段一：实现并验证自研 SASA 计算器

使用 `2iww_H.pdb` 或少量 PDB 文件测试。

目标：

```text
1. 完成原始 Project 2 要求；
2. 输出 total SASA；
3. 输出 atom-level SASA；
4. 输出 residue-level SASA；
5. 输出 chain-level SASA；
6. 证明自己理解 Shrake-Rupley 算法。
```

同时和 FreeSASA 做对照：

```text
1. 自研 total SASA vs FreeSASA total SASA；
2. 自研 residue-level SASA vs FreeSASA residue-level SASA；
3. 计算相对误差和相关性。
```

可以设计如下表格：

| PDB | 自研 SASA | FreeSASA SASA | 相对误差 |
|---|---:|---:|---:|
| 2iww_H | xxxx | xxxx | x.x% |
| complex_1 | xxxx | xxxx | x.x% |
| complex_2 | xxxx | xxxx | x.x% |

如果误差在可接受范围内，例如几个百分点以内，就说明自研实现基本合理。

---

### 阶段二：后续 ΔSASA 用 FreeSASA 主导

对大量复合物 PDB：

```text
1. FreeSASA 计算 apo 状态下的 residue-level SASA；
2. FreeSASA 计算 holo 状态下的 residue-level SASA；
3. ΔSASA = SASA_apo - SASA_holo；
4. 根据阈值生成 interface / non-interface label。
```

典型流程：

```text
1. 对 Chain A 单独保存为 A_only.pdb；
2. 用 FreeSASA 计算 A_only.pdb 的 residue-level SASA；
3. 对 A+B 复合物保存为 AB_complex.pdb；
4. 用 FreeSASA 计算 AB_complex.pdb 的 residue-level SASA；
5. 提取其中 Chain A 的 residue-level SASA；
6. 做差得到 ΔSASA；
7. 对 Chain B 重复相同流程。
```

---

## 五、holo SASA 计算时的关键注意点

holo SASA 不能简单地算完整个复合物后乱减。

对于 Chain A：

```text
SASA_apo(A) = A 链单独存在时，A 链每个残基的 SASA
SASA_holo(A) = A+B 复合物中，A 链每个残基的 SASA
ΔSASA(A) = SASA_apo(A) - SASA_holo(A)
```

对于 Chain B：

```text
SASA_apo(B) = B 链单独存在时，B 链每个残基的 SASA
SASA_holo(B) = A+B 复合物中，B 链每个残基的 SASA
ΔSASA(B) = SASA_apo(B) - SASA_holo(B)
```

也就是说：

> holo 状态下虽然使用的是复合物结构，但最后只提取目标链的残基 SASA。

这样才能得到每条链每个残基的界面标签。

---

## 六、使用 FreeSASA 作为主计算器的好处

后续如果老师问：

> 你们的弱监督标签可靠吗？会不会是自己写的 SASA 算法误差导致模型训练结果不稳定？

可以回答：

```text
我们先实现了自研 SASA 算法，并与成熟工具 FreeSASA 进行结果对照。自研工具用于完成课程项目的算法要求和小规模验证；在大规模 ΔSASA 标签生成阶段，为了提高伪标签稳定性和实验可复现性，我们采用 FreeSASA 作为主计算器。
```

这样回答会非常稳。

---

## 七、方法部分可以怎么写？

### 英文版

```text
To ensure both algorithmic transparency and computational reliability, we adopted a two-stage SASA computation strategy. First, we implemented an in-house Shrake-Rupley-based SASA calculator according to the requirement of Project 2. The calculator parses atomic coordinates from PDB files, expands each atom by the solvent probe radius, samples surface points on the expanded sphere, and estimates the accessible area according to the proportion of exposed points. This implementation was used to validate our understanding of SASA computation and to produce atom-level, residue-level, and chain-level SASA values.

Second, for large-scale ΔSASA-based pseudo-label generation, we used FreeSASA as the primary SASA engine. FreeSASA is a mature open-source SASA calculation library that supports both Lee-Richards and Shrake-Rupley algorithms and provides residue-level SASA output. We compared our in-house calculator with FreeSASA on representative PDB structures to verify consistency. The subsequent ΔSASA labels were generated using FreeSASA to improve robustness and reproducibility.
```

### 中文版

```text
为了兼顾算法透明性与工程可靠性，本文采用两阶段 SASA 计算策略。首先，根据 Project 2 的要求，我们实现了一个基于 Shrake-Rupley 打点法的自研 SASA 计算器。该计算器从 PDB 文件中解析原子坐标，将每个原子的半径扩展为原子半径与溶剂探针半径之和，并通过球面采样点判断其是否被其他原子遮挡，从而估计原子级、残基级和链级 SASA。该实现主要用于验证 SASA 计算原理并完成课程项目的基础算法要求。

其次，在大规模 ΔSASA 弱监督标签生成阶段，本文采用成熟开源工具 FreeSASA 作为主要 SASA 计算引擎。FreeSASA 支持 Lee-Richards 和 Shrake-Rupley 两类经典算法，并能够输出残基层面的 SASA。我们在代表性 PDB 结构上对自研计算器与 FreeSASA 的结果进行一致性对比；在后续批量处理蛋白复合物时，使用 FreeSASA 生成 apo 与 holo 状态下的残基级 SASA，以提高伪标签的稳定性与可复现性。
```

---

## 八、实验对照怎么设计？

建议做两个对照实验。

---

### 实验一：SASA 工具一致性验证

目的：证明自研工具没有明显错误。

比较对象：

```text
自研 SASA
FreeSASA
```

评价指标：

| 指标 | 含义 |
|---|---|
| Total SASA relative error | 总 SASA 相对误差 |
| Residue-level Pearson correlation | 残基级 SASA 趋势是否一致 |
| Mean Absolute Error | 残基级平均绝对误差 |

可以写成：

```text
We first evaluated the consistency between our in-house SASA calculator and FreeSASA. The comparison was conducted at both protein level and residue level. Total SASA relative error was used to measure the global consistency, while residue-level Pearson correlation and mean absolute error were used to evaluate whether the two tools produced similar residue-wise exposure patterns.
```

---

### 实验二：标签生成稳定性验证

目的：证明两种 SASA 计算器在界面标签判断上基本一致。

比较方法：

```text
1. 用自研 SASA 生成 ΔSASA label；
2. 用 FreeSASA 生成 ΔSASA label；
3. 比较两组 label 的一致性。
```

评价指标：

| 指标 | 含义 |
|---|---|
| Label agreement | 两种工具生成标签的一致率 |
| Positive ratio | 正样本比例 |
| Interface overlap | 预测为界面残基集合的重合度 |

如果一致率较高，可以写：

```text
The high agreement between labels generated by our in-house calculator and FreeSASA indicates that our implementation captures the main geometric principle of SASA-based interface identification. Therefore, FreeSASA was used for large-scale pseudo-label generation due to its robustness and efficiency.
```

中文可写：

```text
自研计算器与 FreeSASA 生成的界面标签具有较高一致性，说明自研实现能够捕捉基于 SASA 差分判断蛋白质界面的主要几何规律。考虑到 FreeSASA 在工程稳定性和批量处理方面更成熟，后续大规模 ΔSASA 弱监督标签生成主要采用 FreeSASA 完成。
```

---

## 九、报告中的最终叙事逻辑

最终报告可以这样组织：

```text
1. 我们首先实现了自研 SASA 计算器，完成原始 Project 2 的核心算法要求；
2. 为了验证自研实现的可靠性，我们将结果与成熟工具 FreeSASA 进行对照；
3. 对照结果表明，自研 SASA 在总面积和残基级暴露趋势上与 FreeSASA 基本一致；
4. 考虑到后续 ΔSASA 标签将直接影响 MLP / GCN / EGNN 训练质量，我们在大规模标签生成阶段使用 FreeSASA 作为主计算器；
5. 自研工具保留为算法解释、结果复核和小规模对照模块。
```

这套逻辑非常稳：

```text
不是因为不会写所以用 FreeSASA；
而是因为已经写了，并且用 FreeSASA 验证；
后续为了提高弱监督标签质量和实验可复现性，采用成熟工具作为大规模数据生产引擎。
```

---

## 十、最终建议

最终建议如下：

```text
1. 自研 SASA：必须做，用于完成课程项目核心要求和算法解释；
2. FreeSASA 对照：建议做，用于验证自研实现的正确性；
3. 后续 ΔSASA 批量标签：建议主要用 FreeSASA；
4. MLP / GCN / EGNN 训练：建议使用 FreeSASA 生成的标签，更稳定、更可信。
```

一句话总结：

> **自研 SASA 证明我们会做，FreeSASA 保证后续标签更稳。**
