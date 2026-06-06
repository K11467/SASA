# 3.5 等变残基图建模 & 3.6 跨链空间注意力变体 —— 详细说明

> 本文档对论文 §3.5(等变残基图建模)与 §3.6(跨链空间注意力变体)做逐步展开,
> 并把每一条数学公式映射到代码实现 [`src/sasa_project/train_interface_model.py`](../src/sasa_project/train_interface_model.py)。
> 目的:让读者既能理解"为什么这样设计",也能定位到"代码里是哪一行"。

---

## 一、整体动机

蛋白质界面预测的本质是一个**逐残基二分类**问题:给定目标链(target chain)的每个残基,判断它在与伙伴链(partner chain)结合后是否处于结合界面上。

一个好的几何模型应当满足两点物理直觉:

1. **结构等变性**:把整个复合物在空间中平移 / 旋转,界面残基的身份不应改变。因此模型对节点坐标应是 **E(3)-等变** 的(坐标随之变换,但预测不变)。
2. **几何邻近性**:界面残基在几何上紧邻伙伴链。模型需要能利用"残基与伙伴链的空间距离"这一强信号。

§3.5 的 EGNN 满足第 1 点(在单链图上建模等变几何);§3.6 的跨链变体在此基础上显式引入第 2 点(让目标残基直接"看到"伙伴链坐标)。

---

## 二、3.5 等变残基图建模(EGNN)

### 2.1 图的构建

目标链被形式化为几何图 $G=(V,E)$:

- **节点** $V$:每个残基的 Cα 原子,携带
  - 节点特征 $h_i$(ESM 序列嵌入 + 结构特征,经 `input_proj` 线性映射到 `hidden_dim`);
  - 三维坐标 $p_i \in \mathbb{R}^3$(Cα 坐标)。
- **空间边** $E$:连接欧几里得距离在截止半径内的残基对
  $$e_{ij}\in E \iff \lVert p_i - p_j \rVert_2 \le r.$$
  截止半径 $r \in \{6, 8, 10, 12\}\ \text{Å}$ 是被系统评估的超参数(代码中以 `distance_cutoff` / `d6/d8/d10/d12` 体现,默认最优为 8 Å)。

> 代码位置:图构建在 `build_egnn_graphs(...)`,`edge_index` 形如 `(2, E)` 的 `(src, dst)` 张量。

### 2.2 消息传递:三个更新函数

采用 **3 层 EGNN**(`n_layers=3`)。每一层 $l$ 同时更新**不变的节点特征** $h$ 和**等变的空间坐标** $p$。

代码定义在 `EGNNModel._EGNNLayer`([train_interface_model.py:236](../src/sasa_project/train_interface_model.py#L236) 起)。

#### (1) 边消息 $m_{ij}$ —— 公式 (2)

论文形式:
$$m_{ij} = \phi_e\!\left(h_i^l,\ h_j^l,\ \lVert p_i^l - p_j^l\rVert^2 / 10,\ e_{ij}\right)$$

代码实现:
```python
diff = x[src] - x[dst]                                   # (E,3) 坐标差
dist = (diff**2).sum(-1, keepdim=True).clamp(min=1e-8).sqrt()  # (E,1) 欧氏距离
diff_norm = diff / dist                                  # 单位方向向量
dist_norm = dist / 10.0                                  # 归一化标量距离
m_input = torch.cat([h[src], h[dst], dist_norm], dim=-1) # 2*node_dim+1
m = self.phi_e(m_input)                                  # (E, hidden_dim)
```
其中 $\phi_e$ 是两层 `Linear + SiLU` 的 MLP。

> **实现说明(与公式 (2) 的差异)**:论文写的是距离平方 $\lVert\cdot\rVert^2/10$,而代码用的是**归一化后的标量距离** $\lVert\cdot\rVert/10$(而非平方)。这是一个有意的**数值稳定性修正**:原始 Cα 坐标范围可达 10–300 Å,直接喂距离平方会让数值爆炸(代码注释 "Fix 2")。两者携带的几何信息等价(都是 $i,j$ 间距离的单调函数),只是数值范围更友好。

#### (2) 坐标更新 $p^{l+1}$ —— 公式 (3),等变部分

论文形式:
$$p_i^{l+1} = p_i^l + \sum_{j\in\mathcal N(i)} \frac{p_i^l - p_j^l}{\lVert p_i^l - p_j^l\rVert}\,\phi_x(m_{ij})$$

代码实现:
```python
coord_weight = torch.tanh(self.phi_x(m))   # (E,1) 标量权重，tanh 约束到 [-1,1]
coord_agg = torch.zeros_like(x)
coord_agg.scatter_add_(0, dst..., diff_norm * coord_weight)  # 按目标节点聚合
x = x + coord_agg
```
- $\phi_x$ 把每条边的消息映射为一个**标量权重**;
- 权重乘以**单位方向向量** $\frac{p_i-p_j}{\lVert p_i-p_j\rVert}$,再按邻居求和,加回坐标。

> **为什么这一步是等变的**:坐标更新量是一组"方向向量的加权和"。当整个结构旋转 $R$、平移 $t$ 时,方向向量随 $R$ 同步旋转,而权重 $\phi_x(m_{ij})$ 只依赖**距离**(旋转平移不变量),因此更新量也随 $R$ 旋转、随平移不变 —— 这正是 E(3)-等变。
>
> **实现说明**:代码用 `tanh` 把坐标权重约束在 $[-1,1]$("Fix 3"),防止多层堆叠时坐标无限漂移。

#### (3) 特征更新 $h^{l+1}$ —— 公式 (4),不变部分

论文形式:
$$h_i^{l+1} = \phi_h\!\left(h_i^l,\ \sum_{j\in\mathcal N(i)} m_{ij}\right)$$

代码实现:
```python
m_agg = torch.zeros(N, hidden_dim)
m_agg.scatter_add_(0, dst..., m)          # 聚合邻居消息
h = h + self.phi_h(torch.cat([h, m_agg], dim=-1))   # 残差连接
```
- $\phi_h$ 把"自身特征 + 邻居消息和"映射回节点维度;
- 由于消息 $m_{ij}$ 只依赖距离等不变量,$h$ 始终是**旋转平移不变**的标量特征。

> **实现说明**:代码加了**显式残差连接** `h = h + φ_h(...)`("Fix 1"),稳定 3 层堆叠的训练。

其中 $\phi_e,\phi_x,\phi_h$ 均为参数化的多层感知机(MLP)。

### 2.3 分类头

3 层消息传递后,取最终的节点特征 $h_i$,经一个 `Linear → ReLU → Dropout → Linear` 的分类头输出每个残基的 logit:
```python
return self.classifier(h).squeeze(-1)     # (N,) 每个残基一个界面打分
```

### 2.4 小结:EGNN 的不变量 / 等变量

| 量 | 变换性质 | 含义 |
|---|---|---|
| 节点特征 $h_i$ | E(3)-**不变** | 用于最终分类的标量表示 |
| 节点坐标 $p_i$ | E(3)-**等变** | 随结构旋转平移而变换 |
| 边距离 $\lVert p_i-p_j\rVert$ | **不变** | 喂给消息函数的几何信号 |
| 最终预测 | **不变** | 旋转/平移复合物不改变界面判断 |

---

## 三、3.6 跨链空间注意力变体(Cross-Chain Variant)

### 3.1 动机

§3.5 的 EGNN 只在**目标链单链图**上传播,无法直接利用伙伴链的几何布局。一种朴素做法是把目标链 + 伙伴链拼成一张**完整二分复合物图**再做 GNN,但代价高、且需要完整复合物结构。

跨链变体的问题是:**能否在不构建完整二分图的情况下,用一个轻量模块显式引入伙伴链定位信息?** 做法是在 EGNN 输出之后,接一个**距离权重的跨链软注意力**,让每个目标残基直接"感知"伙伴链 Cα 的空间分布。

> 代码:`CrossChainEGNNModel`([train_interface_model.py:325](../src/sasa_project/train_interface_model.py#L325) 起),其 EGNN 主干与 §3.5 完全相同,仅在分类头之前插入 `_CrossChainAttention`。

### 3.2 距离权重软注意力 —— 公式 (5)

对目标残基 $i$ 和伙伴残基 $j$,计算连续的空间注意力权重:
$$\alpha_{ij} = \mathrm{softmax}_j\!\left(-\frac{\lVert p_i - q_j\rVert_2}{\sigma}\right)$$

其中 $q_j$ 是伙伴链的 Cα 坐标,$\sigma$ 是**可学习的缩放(温度)参数**。

代码实现(`_CrossChainAttention.forward`):
```python
sigma = self.log_sigma.exp().clamp(min=0.5, max=20.0)   # 可学习温度，约束范围
diff = pos.unsqueeze(1) - partner_pos.unsqueeze(0)       # (N, P, 3)
dist = diff.norm(dim=-1)                                  # (N, P) 目标-伙伴距离矩阵
attn = F.softmax(-dist / sigma, dim=-1)                   # (N, P) 公式 (5)
partner_h = self.partner_proj(partner_pos)               # (P, H) 伙伴坐标投影
context = attn @ partner_h                                # (N, H) 上下文聚合
return self.combine(torch.cat([h, context], dim=-1))     # 拼接 + MLP
```

逐步解读:

1. **距离即注意力**:用负距离 $-\lVert p_i - q_j\rVert/\sigma$ 做 softmax —— 离目标残基越近的伙伴残基,权重越大。这把"界面残基紧邻伙伴链"的物理直觉直接编码进注意力。
2. **可学习温度 $\sigma$**:以 `log_sigma`(初值 2.0)参数化并 `exp` 还原,约束在 $[0.5, 20]$ Å。$\sigma$ 小 → 注意力尖锐(只看最近的伙伴残基);$\sigma$ 大 → 注意力平滑(综合更多伙伴残基)。让模型自己学出合适的"感受野"。
3. **Key/Value = 伙伴坐标投影**:伙伴链坐标 $q_j$ 经 `partner_proj`(3→H 的 MLP)投影为 $\text{partner\_h}_j$。
4. **上下文向量**:$c_i = \sum_j \alpha_{ij}\,\text{partner\_h}_j$(矩阵乘 `attn @ partner_h`)。
5. **融合**:把上下文 $c_i$ 与目标残基的 EGNN 表示 $h_i$ **拼接**,再经 `combine`(2H→H 的 MLP)得到增强表示,送入最终分类 MLP。

> 对应论文末句:"伙伴坐标被投影后与目标表示拼接,然后进入最终的 MLP 进行分类。"

### 3.3 前向流程对比

```
EGNN (§3.5):      x_feat → input_proj → [3×EGNNLayer] → classifier → logits
Cross-Chain(§3.6): x_feat → input_proj → [3×EGNNLayer] → cross_attn(partner_pos) → classifier → logits
                                                          └── 唯一新增模块 ──┘
```
代码见 `_CrossChainEGNN.forward`([train_interface_model.py:414](../src/sasa_project/train_interface_model.py#L414)):
```python
h = self.input_proj(x_feat)
for layer in self.egnn_layers:
    h, pos = layer(h, pos, edge_index)
h = self.cross_attn(h, pos, partner_pos)   # ← 跨链注意力
return self.classifier(h).squeeze(-1)
```

### 3.4 伙伴链坐标的加载与退化处理

跨链变体在数据侧需要额外加载伙伴链 Cα 坐标(`build_cross_chain_graphs` / `_load_partner_positions`):

- 从 manifest 定位伙伴链,抽取其残基的 Cα 坐标作为 `partner_pos`(形如 `(P,3)`)。
- **退化兜底**:若某复合物缺失伙伴链坐标,用目标链坐标的几何中点作为单个 dummy 伙伴节点,并打印 `[warn] ... missing partner chain positions`,保证训练不中断。

---

## 四、与基线的关系 & 实验结论

- **EGNN(§3.5)** 是几何主干基线;**Cross-Chain(§3.6)** 在其上增加跨链注意力,是消融对比项。
- 在 holo-aware PDBtest_315 诊断设置下观察到的**指标权衡**(来自 [`data/processed/benchmark_pdbtest315_metrics_650m.csv`](../data/processed/benchmark_pdbtest315_metrics_650m.csv)):

  | 指标 | EGNN 基线 | + 跨链注意力 | 变化 |
  |---|---|---|---|
  | Recall | 0.6252 | **0.6629** | ▲ |
  | F1 | 0.6761 | **0.6816** | ▲ |
  | AUPRC | 0.7408 | **0.6906** | ▼ |

- **结论(保守表述)**:跨链上下文能提升**召回**,但对**排序质量(AUPRC)** 的影响是 benchmark 相关、并非一致为正。因此论文中跨链注意力被定位为**分析性变体**而非主推方法,其收益需在更严格的设置(no-SASA 输入、匹配的 manifest)下重新验证。

---

## 五、关键超参数一览

| 超参数 | 取值 | 说明 |
|---|---|---|
| EGNN 层数 `n_layers` | 3 | 消息传递层数 |
| 图截止半径 $r$ | {6, 8, 10, 12} Å | 系统评估,8 Å 通常最优 |
| 隐藏维度 `hidden_dim` | 可配 | 节点表示维度 |
| 激活函数 | SiLU(主干)/ ReLU(分类头) | — |
| 距离归一化 | `/10`(EGNN)、可学习 $\sigma\in[0.5,20]$(跨链) | 数值稳定 + 自适应感受野 |
| 坐标权重约束 | `tanh` ∈ [-1,1] | 防止坐标漂移 |
| 正则 | Dropout + 残差连接 | 稳定 3 层堆叠训练 |

---

## 六、实现与论文的差异速查

| 论文 | 代码实现 | 原因 |
|---|---|---|
| 公式 (2) 用距离平方 $\lVert\cdot\rVert^2/10$ | 用归一化标量距离 $\lVert\cdot\rVert/10$ | 数值稳定(避免大坐标下平方爆炸),信息等价 |
| 公式 (3) 坐标更新 | 额外 `tanh` 约束权重 | 防止多层坐标无限漂移 |
| 公式 (4) 特征更新 | 额外显式残差 `h + φ_h(...)` | 稳定深层训练 |
| 公式 (5) 跨链注意力 | 与论文一致,$\sigma$ 以 `log_sigma` 参数化并裁剪 | 保证 $\sigma>0$ 且范围合理 |

> 这些差异都是**工程化的数值稳定性增强**,不改变模型的等变性与几何语义。若撰写论文方法部分,建议把公式 (2) 的 $\lVert\cdot\rVert^2$ 与代码统一(改为 $\lVert\cdot\rVert$ 或在代码中改回平方),避免审稿人质疑。
