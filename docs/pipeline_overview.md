# 项目二 B + D 组合整体思路介绍

> 当前实现状态：正式实验已使用 ESM-2 650M、1289 维多模态特征、EGNN 和
> cross-chain EGNN。轻量 GCN 保留为 baseline。主训练集为 500 个复合物、
> 95,005 个残基；Dset_186-local 和 PDBtest_315-local 外部评测均已完成。

下面内容面向人工智能专业学生，尽量用 AI 任务建模的方式解释这个计算生物学项目。

---

## 1. 这个项目到底在预测什么？

这个项目要预测的是 **PPI interface residue**，也就是：**蛋白质-蛋白质相互作用界面残基**。

一个蛋白质由很多氨基酸残基组成。两个蛋白质结合时，不是整个蛋白都参与结合，而是表面的一部分残基真正接触另一个蛋白。这些真正参与结合的残基，就叫 protein-protein interaction site，或者更具体地说，叫 interface residues。

从 AI 角度看，这个任务就是一个**节点分类问题**：

```text
输入：一个蛋白质结构 / 序列
输出：每个残基是不是结合界面
标签：0 或 1
```

可以类比成图神经网络任务：

```text
每个节点 = 一个氨基酸残基
每条边 = 两个残基在三维空间中距离较近
节点特征 = 序列特征 + 几何特征 + SASA 特征
任务 = 对每个节点做二分类
```

---

## 2. SASA 在这里起什么作用？

SASA 是 **Solvent Accessible Surface Area**，即溶剂可及表面积。

可以把它理解成：**一个残基暴露在外面的程度**。

如果一个残基埋在蛋白内部，它几乎接触不到水，SASA 会比较低。如果一个残基在蛋白表面，它能接触到溶剂，SASA 会比较高。

对于蛋白质相互作用来说，界面残基通常有一个特点：

> 在单独存在时，它可能暴露在溶剂中；但当另一个蛋白靠上来之后，这部分表面被遮住了。

也就是说：

```text
单独蛋白状态：这个残基暴露，SASA 较大
复合物状态：这个残基被另一个蛋白挡住，SASA 变小
```

所以可以通过 SASA 的变化来判断一个残基是不是界面残基。

---

## 3. 什么是 ΔSASA？

ΔSASA 叫做**差分 SASA**。

核心公式是：

```text
ΔSASA = SASA_apo - SASA_holo
```

其中：

```text
apo：蛋白单独存在的状态
holo：蛋白和另一个蛋白结合后的复合物状态
```

通俗解释：

```text
SASA_apo  = 这个残基单独存在时有多少面积暴露在外面
SASA_holo = 这个残基和别的蛋白结合后还有多少面积暴露在外面
ΔSASA     = 结合后被遮住了多少面积
```

如果一个残基在结合前暴露很多，结合后被遮住很多，那么 ΔSASA 很大，说明它很可能在蛋白质结合界面上。

所以方向 B 的核心就是：**不人工标注界面残基，而是用 ΔSASA 自动生成标签。**

例如可以设定：

```text
如果 ΔSASA > 1.0 Å²，则认为该残基是界面残基，label = 1
否则 label = 0
```

阈值可以实验比较，例如 0.5 Å²、1.0 Å²、2.0 Å²、5.0 Å²，然后观察哪个阈值构造出来的伪标签训练效果最好。

---

## 4. 为什么说 B 是“弱监督”？

传统监督学习需要真实标签，例如：

```text
残基 A 是界面
残基 B 不是界面
残基 C 是界面
...
```

但真实标签往往需要生物实验或高质量数据库，获取成本比较高。

而 B 的做法是：**用几何计算自动生成标签**。

这些标签不是人工确认的严格真值，而是根据 ΔSASA 推断出来的，所以叫 **weak supervision，弱监督**。

从 AI 角度看，这很像：

```text
用规则 / 启发式方法自动打伪标签
然后训练神经网络
```

类似 NLP 里用规则生成 pseudo-label，再训练分类器。

---

## 5. ESM-2 是什么？

ESM-2 是一个大型蛋白语言模型，可以把它类比成 NLP 里的 BERT / Transformer。

区别是：

```text
BERT 学的是自然语言序列
ESM-2 学的是蛋白质氨基酸序列
```

蛋白质序列可以看成一个由 20 种氨基酸字母组成的字符串，例如：

```text
MKTFFVLLLCTFT...
```

ESM-2 会把每个氨基酸残基编码成一个高维向量。例如一个残基可以被表示为 1280 维 embedding。这个 embedding 中包含很多进化、功能、结构倾向信息。虽然它只看序列，但由于训练数据巨大，它能学到很多蛋白质规律。

所以方向 D 的意思是：**不只用 SASA 这种人工几何特征，还加入 ESM-2 这种深度蛋白语言模型特征。**

---

## 6. 为什么 B + D 组合比较好？

因为 B 和 D 分别解决两个问题。

### 6.1 B 解决“标签从哪里来”

原本项目二只是算 SASA，算完一个总面积就结束了，偏像一个几何算法作业。

但 B 让 SASA 变成了标签生成工具：

```text
PDB 复合物结构
↓
计算单链 SASA
↓
计算复合物 SASA
↓
得到 ΔSASA
↓
生成界面 / 非界面标签
```

这样 Project 2 就不只是一个孤立算法，而是成为后续机器学习任务的数据构造模块。

### 6.2 D 解决“模型输入太弱”

如果只用 SASA，一个残基只有一个或几个几何数值，信息量太少。

但加入 ESM-2 后，每个残基有：

```text
序列语义信息：ESM embedding
几何暴露信息：SASA / ΔSASA / residue-level SASA
空间邻域信息：GCN / EGNN 图结构
```

这就变成了真正的多模态融合：

```text
序列模态 + 几何模态 + 结构图模态
```

---

## 7. 整体 pipeline 设计

整个项目可以拆成 6 个模块。

### 模块 1：读取 PDB 结构

输入是蛋白质复合物结构文件。例如一个 PDB 里面有两条链：

```text
Chain A：蛋白 A
Chain B：蛋白 B
```

需要解析每个原子的信息：

```text
原子名
残基名
链名
残基编号
x, y, z 坐标
```

然后按残基聚合，最终得到：

```text
Residue 1: 包含若干原子
Residue 2: 包含若干原子
Residue 3: 包含若干原子
...
```

### 模块 2：用 Project 2 计算 residue-level SASA

原始项目二要求输出整个蛋白的总 SASA。但在 B + D 项目里，需要扩展成：不仅计算 protein-level SASA，还要计算 residue-level SASA。

也就是每个残基的 SASA。做法是：

```text
每个原子有自己的 SASA
一个残基的 SASA = 该残基所有原子 SASA 之和
```

例如：

```text
Residue 25 的 SASA = N 原子 SASA + CA 原子 SASA + C 原子 SASA + O 原子 SASA + 侧链原子 SASA
```

### 模块 3：计算 apo SASA 和 holo SASA

这是 B 的核心。假设关注 Chain A，需要计算两种状态下 Chain A 每个残基的 SASA。

第一种：apo 状态，也就是只看 Chain A。

```text
输入结构：Chain A
输出：Chain A 每个残基的 SASA_apo
```

第二种：holo 状态，也就是 Chain A 和 Chain B 一起存在。

```text
输入结构：Chain A + Chain B
输出：在复合物中 Chain A 每个残基的 SASA_holo
```

注意：holo 状态下，仍然只统计 Chain A 的残基 SASA，但遮挡时要考虑 Chain B 的原子。

然后：

```text
ΔSASA = SASA_apo - SASA_holo
```

如果某个残基因为结合被 Chain B 挡住，它的 SASA_holo 会下降，ΔSASA 就会变大。

### 模块 4：用 ΔSASA 生成界面标签

设定阈值：

```text
label = 1 if ΔSASA > threshold
label = 0 otherwise
```

例如：

```text
ΔSASA > 1.0 Å² → interface residue
ΔSASA ≤ 1.0 Å² → non-interface residue
```

这个标签就是训练 GCN / EGNN 的监督信号。

可以在实验里比较不同阈值：

| 阈值 | 可能效果 |
|---:|---|
| 0.5 Å² | 标签更宽松，界面残基更多，但噪声可能更大 |
| 1.0 Å² | 比较常用的中间选择 |
| 2.0 Å² | 标签更严格，界面更可信，但正样本更少 |
| 5.0 Å² | 只保留强界面残基，可能漏掉弱接触位点 |

这一部分很适合写成报告中的实验分析：

> ΔSASA 阈值会影响伪标签质量和正负样本比例，从而影响 PPI 位点预测性能。

### 模块 5：提取 ESM-2 序列嵌入

对于每条蛋白链，把氨基酸序列输入 ESM-2。

例如 Chain A 序列：

```text
MKTAYIAKQRQISFVKSHFSRQ...
```

ESM-2 输出每个残基的向量：

```text
Residue 1  → e1
Residue 2  → e2
Residue 3  → e3
...
```

如果用 ESM-2 650M，常见 embedding 是 1280 维。所以每个残基会有：

```text
ESM_feature_i ∈ R^1280
```

这个向量可以理解为：该残基在蛋白语言模型语义空间中的表示。

### 模块 6：构建图并训练模型

现在每个残基都有特征和标签。

节点：

```text
每个残基是一个节点
```

边：

```text
如果两个残基的 CA 原子距离小于某个阈值，例如 8Å 或 10Å，就连一条边
```

节点特征可以拼接：

```text
x_i = [ESM_embedding_i, SASA_apo_i, SASA_holo_i, 其他几何特征]
```

不过要注意一个问题：如果把 ΔSASA 同时作为标签来源，又作为输入特征，可能会有信息泄漏。

所以更严谨的设计是：训练预测器时输入：

```text
ESM_embedding_i + apo SASA_i + holo SASA_i + residue 几何特征
```

标签来自：

```text
ΔSASA_i > threshold
```

也就是说：**ΔSASA 用来生成 label，不要直接作为模型输入**。否则模型可能学到一个太直接的规则，实验意义下降。

baseline 可以用轻量 GCN；正式实验使用 EGNN：

```text
Input node features
↓
EGNN layer
↓
ReLU
↓
EGNN layer
↓
ReLU
↓
MLP classifier
↓
interface probability
```

输出：

```text
每个残基是界面的概率 p_i
```

损失函数可以使用 Binary Cross Entropy Loss。因为界面残基通常比非界面残基少，所以可以加入 class weight 或 focal loss。

---

## 8. 这个项目和原始 Project 2 的关系

原始 Project 2 是：

```text
输入 PDB + 溶剂半径
输出 总 SASA
```

升级后的项目是：

```text
输入 PDB 复合物 + 溶剂半径 + ESM-2
输出 每个残基是否为 PPI 界面
```

原始 Project 2 在新项目里承担三个作用：

1. 计算每个原子的 SASA；
2. 聚合得到每个残基的 SASA；
3. 计算 apo / holo 的 SASA 差值，用于生成弱监督标签。

所以报告中可以强调：

> 本项目并未停留在 SASA 数值计算本身，而是进一步将自研 SASA 计算器作为结构生物学特征提取与弱监督标签生成工具，构建了一个面向蛋白质相互作用界面预测的机器学习流程。

---

## 9. 推荐实验设计

### 实验一：特征消融实验

比较：

```text
ESM only
SASA only
ESM + SASA
```

对应含义：

| 模型 | 输入特征 | 目的 |
|---|---|---|
| ESM only | 只有序列语言模型嵌入 | 看纯序列信息能做到什么程度 |
| SASA only | 只有几何暴露特征 | 看 Project 2 的 SASA 特征是否有效 |
| ESM + SASA | 序列 + 几何 | 验证多模态融合是否提升性能 |

如果结果是：

```text
ESM + SASA > ESM only > SASA only
```

就很好写分析：

> ESM-2 捕获了进化和序列上下文信息，SASA 补充了残基空间暴露程度，两者具有互补性。

### 实验二：ΔSASA 阈值敏感性实验

比较不同阈值：

```text
0.5, 1.0, 2.0, 5.0 Å²
```

观察：

```text
正样本比例
Precision
Recall
F1
AUC
AUPRC
```

可以分析：

```text
阈值太低 → 正样本多，但标签噪声大
阈值太高 → 标签更干净，但正样本太少，召回率下降
```

这部分可以体现对弱监督标签质量的理解。

### 实验三：和 baseline 对比

最简单 baseline：

```text
MLP on ESM
MLP on SASA
GCN / EGNN on ESM + SASA + structure
```

更完整一点可以参考 GraphPPIS。它显式使用 SASA，并且提供了数据集划分、代码和预训练模型。不一定要完整复现 GraphPPIS，但可以借鉴它的任务设定和数据划分。

---

## 10. 这个项目为什么适合人工智能专业学生？

因为它本质上是一个非常标准的 AI pipeline，只是数据来自生物学。

| 生物学概念 | AI 视角 |
|---|---|
| 蛋白质序列 | token 序列 |
| 氨基酸残基 | token / graph node |
| 蛋白质结构 | 3D graph |
| SASA | 几何特征 |
| ΔSASA | 伪标签生成规则 |
| PPI interface | 节点二分类标签 |
| ESM-2 | 预训练 Transformer encoder |
| GCN / EGNN | 图神经网络 baseline / 正式等变图模型 |
| 消融实验 | 多模态特征贡献分析 |

所以不需要一开始就从生物化学角度理解所有细节。可以把它看成：

> 用一个预训练蛋白语言模型提取序列 embedding，用一个几何算法提取结构特征，再用图神经网络做节点分类。

---

## 11. 报告里的主线写法

可以把项目故事线写成：

> 原始 Project 2 仅要求计算蛋白质的溶剂可及表面积。为了进一步提升项目的计算生物学意义，本文将 SASA 计算扩展到残基层面，并利用蛋白质复合物结合前后的 SASA 变化构造 ΔSASA 弱监督标签。具体而言，若某个残基在复合物形成后 SASA 显著下降，则说明该残基可能被另一条蛋白链遮挡，因此被视为潜在蛋白质相互作用界面残基。在此基础上，本文进一步引入 ESM-2 蛋白语言模型提取残基级序列嵌入，并与 SASA 几何特征进行融合，构建轻量级图神经网络完成 PPI 界面残基预测。该设计使原本的几何表面积计算任务扩展为一个弱监督、多模态、结构感知的蛋白质界面预测任务。

这段可以直接作为报告的“整体方法概述”。

---

## 12. 最小可行版本

为了避免项目太大，建议先做一个最小可行版本：

```text
1. 读取 PDB 复合物
2. 按链拆分结构
3. 用自己的 SASA 程序计算 apo SASA
4. 用自己的 SASA 程序计算 holo SASA
5. 得到每个残基的 ΔSASA
6. 用阈值生成 interface label
7. 用 ESM-2 提取每个残基 embedding
8. 拼接 apo SASA 作为特征
9. 训练 MLP、GCN、EGNN 或 cross-chain EGNN 做残基二分类
10. 做 ESM only / SASA only / ESM + SASA 消融
```

如果只做最小 baseline，可以先不用图模型，先做：

```text
ESM embedding + SASA → MLP → interface / non-interface
```

这已经能讲清楚 B + D 的核心。

如果时间够，再升级成：

```text
ESM embedding + SASA + structure → EGNN → interface / non-interface
```

EGNN 版本是当前正式的“结构感知预测”实现。

---

## 13. 一句话总结

这个 B + D 组合的核心思想是：

> **用 SASA 的变化 ΔSASA 自动判断哪些残基在蛋白结合时被遮挡，从而生成界面伪标签；再把 ESM-2 学到的序列语义信息和 SASA 提供的三维几何暴露信息融合起来，训练一个模型预测蛋白质相互作用界面。**

也就是说，项目从：

```text
计算一个蛋白表面积
```

升级成了：

```text
用自研 SASA 工具 + 蛋白语言模型 + 图神经网络，做蛋白相互作用界面预测
```
