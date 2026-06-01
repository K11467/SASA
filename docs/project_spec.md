# 项目1：序列聚类

## 1. 背景

在蛋白的早期设计中，需要对序列的多个点进行突变，或取大量相似序列进行研究，从而形成一个序列库。此时，需要将这些序列进行相似聚类，以将庞大的序列库压缩到少量的聚类中，以便于后续针对蛋白结构的进一步研究。

## 2. 程序输入

1. 输入Fasta文件路径

Fasta文件的格式如下：

```
>A
CYSS
>B
TVGGA
...
```

即，每条序列由`>...`作为标题行，下一行为序列。

实际的Fasta文件未规定序列只能为1行，而是仅以`>`标记新序列的开始。是否考虑此问题是实现定义的。

**保证每条序列是等长的。**

2. 一个整数，表示单个聚类中，最大允许的序列不相似度

3. 输出文件路径

允许实现定义额外的输入参数。

## 3. 程序输出

一个N行的文件，表示序列被分成了N组，每行为一组。

每一行为当前组的所有序列，以空格隔开。

## 4. 实现说明

### 4.1 序列不相似度的定义

对于两条等长序列，计算其对应位置上的不同的残基数量，即为序列不相似度。

例：

```
// 不相似度为2
AAAAAA
AABABA
  ^ ^
```

### 4.2 聚类

给定N条序列，需要将其尽可能的安排在较少的组中，使得每一组中的所有序列之间的相似度小于给定值。

例如，给定序列：

```
AAAAAA
AABAAA
CCCCCC
CCDCCC
```

若最大允许的序列不相似度为1，则可以得到以下分组：

```
AAAAAA AABAAA
CCCCCC CCDCCC
```

若最大允许的序列不相似度为6（或更高），则可以得到以下分组：

```
AAAAAA AABAAA CCCCCC CCDCCC
```

## 5. 测试文件

```
Seq.fasta
```

# 项目2：溶剂可及表面积

## 1. 背景

蛋白和细胞膜一定存在于溶剂（如水，油）中，因此其表面会与溶剂相互结合。然而，蛋白表面不是光滑的平面，蛋白也不是平板状，而是会形成各种形状，如封闭球形，管状等，溶剂本身也不是无限小的点，而是有体积的。这就导致蛋白-溶剂界面不是完全贴合的，而是会形成空腔，在极端情况下，甚至会形成以下情况：

* 完全封闭的球会使其内部形成空腔
* 管道状蛋白会使管道形成空腔，将比较大的溶剂卡在外面进不来

通过对这些空腔的研究，可以获知蛋白在溶剂中的状态。进一步的，可以研究蛋白的表面电性，溶解度，蛋白结构的变化带来的影响等，这些都是非常重要的问题。

## 2. 程序输入

1. PDB文件路径

PDB文件是蛋白结构文件，每一行表示一个原子，格式如下：

| 列号     | 含义                                         |
| -------- | -------------------------------------------- |
| [0, 4)   | "ATOM"，无需关注不以此开头的行               |
| [6, 11)  | 原子编号                                     |
| [12, 16) | 原子名                                       |
| [16, 17) | 原子变位指示，无需关注                       |
| [17, 20) | 残基名                                       |
| [21, 22) | 链名                                         |
| [22, 26) | 残基编号                                     |
| [26, 27) | 残基插入符，其与残基编号共同构成残基完整编号 |
| [30, 38) | 原子x坐标                                    |
| [38, 46) | 原子y坐标                                    |
| [46, 54) | 原子z坐标                                    |
| [54, 80) | 无需关注                                     |

2. 溶剂半径

## 3. 程序输出

一个浮点数，表示溶剂可及表面积。

## 4. 实现说明

### 4.1 原子半径

蛋白中的所有原子半径如下表所示：

| 原子 | 原子半径(Am) |
| ---- | ------------ |
| C    | 1.77         |
| N    | 1.66         |
| O    | 1.50         |
| S    | 1.89         |
| H    | 1.20         |

注：PDB文件中的原子名并非元素符号，但其第一个字母一定为元素符号。


### 4.2 溶剂可及表面积

将单个溶剂分子简化成一个半径为R的球，则溶剂可及表面积为这个球在整个蛋白表面滚动，所形成的滚动面的面积。如下图所示：

![](./Roll.PNG)

### 4.3 格点坐标

半径为1的正60面体格点坐标见`Dot.txt`文件。

## 5. 测试文件

```
2iww_H.pdb
```

# 项目3：氢键能量最小化

## 1. 背景

氢键是蛋白分子中最重要的非键作用之一，蛋白的高级结构受氢键影响很大。然而，由于氢原子太小，其在x光衍射中解析困难，导致现存的大量蛋白结果缺失氢原子坐标。因此，基于计算化学技术，将氢原子坐标直接计算出来，是一个非常重要的技术。

## 2. 程序输入

1. 输入PDB文件路径
2. 输出PDB文件路径

## 3. 程序输出

能量最小化的PDB文件。

## 4. 实现说明

### 4.1 氢键

氢键是由`氢键供体-氢原子-氢键受体`构成的有极性的非键作用。

氢键供体如下表所示：

| 残基名 | 原子名       |
| ------ | ------------ |
| 任何   | N            |
| ARG    | NE，NH1，NH2 |
| ASN    | ND2          |
| GLN    | NE2          |
| HIS    | NE2，ND1     |
| LYS    | NZ           |
| SER    | OG           |
| THR    | OG1          |
| TRP    | NE1          |
| TYR    | OH           |

氢键受体如下表所示：

| 残基名 | 原子名   |
| ------ | -------- |
| 任何   | O        |
| ASN    | OD1      |
| ASP    | OD1，OD2 |
| GLN    | OE1      |
| GLU    | OE1，OE2 |
| HIS    | ND1      |
| SER    | OG       |
| THR    | OG1      |
| TYR    | OH       |

### 4.2 氨基酸

20种氨基酸的结构见`amino-acids.pdf`。

氨基酸经脱水缩合后剩下的部分被称为氨基酸残基。

### 4.3 化学键的旋转

在各种化学键中，只有一部分单键可以旋转，不可旋转的键如下：

* 双键，三键
* 酰胺键，即肽键。此键为共振式单双键，不是一个纯粹的单键
* 离域键，最常见的即为苯环。离域键处于自稳定的状态，极难被破坏

### 4.4 氢键的判定

形成氢键的原子应同时满足以下所有条件：

1. 供受体原子间距离小于3.5
2. `供体原子-氢原子-受体原子`的夹角大于`2Pi/3`

### 4.5 氢键能量最小化

化学分子具有降低自由能的倾向，而氢键的形成将会降低自由能，因此，分子倾向于形成尽可能多的氢键。又因为单键是可以旋转的，所以，程序的目标是通过旋转各个单键，使得蛋白分子形成尽可能多的氢键。

### 4.6 氢原子初始化

输入PDB是没有氢原子的，因此，程序首先需要为每个残基添加氢原子。实现方法是：取一个含有氢原子的完整残基，将两个残基的`N + CA + C`三个原子叠合在一起，从而，顺带将所有的氢原子安装到目标位置。

`StdResidue.pdb`文件中含有20个完整残基，可用于氢原子初始化。

将两组坐标叠合在一起的算法见下文。

## 5. 测试文件

```
2iww.pdb
```

## 6. 提示

### 6.1 距离

两个三维坐标之间的距离为：

$$
d=\sqrt{\left( x_0-x_1 \right) ^2+\left( y_0-y_1 \right) ^2+\left( z_0-z_1 \right) ^2}
$$

### 6.2 键角

两个三维向量之间的夹角为：

$$
a=\mathrm{arc}\cos \left( \frac{\overrightarrow{v_0}\cdot \overrightarrow{v_1}}{\left| v_0 \right|\left| v_1 \right|} \right)
$$

### 6.3 轴角旋转矩阵

绕**单位长度**轴$(x, y, z)$旋转$\theta$角度的旋转变换矩阵为：
$$
R=\left[ \begin{matrix}
	\cos \theta +x^2\left( 1-\cos \theta \right)&		xy\left( 1-\cos \theta \right) -z\sin \theta&		xz\left( 1-\cos \theta \right) +y\sin \theta\\
	xy\left( 1-\cos \theta \right) +z\sin \theta&		\cos \theta +y^2\left( 1-\cos \theta \right)&		yz\left( 1-\cos \theta \right) -x\sin \theta\\
	xz\left( 1-\cos \theta \right) -y\sin \theta&		yz\left( 1-\cos \theta \right) +x\sin \theta&		\cos \theta +z^2\left( 1-\cos \theta \right)\\
\end{matrix} \right]
$$
此外，有以下结论：
$$
\text{因为：}R^{-1}=R^T\text{，所以：}Rx=xR^T
$$
后式在编程中比较方便。

### 6.4 计算叠合矩阵

两组坐标的叠合矩阵可由以下代码计算得到：

```cpp
#include <Eigen/Dense>

using namespace Eigen;

tuple<RowVector3d, Matrix3d, RowVector3d> calcSuperimposeRotationMatrix(
    const Matrix<double, Dynamic, 3> &tarCoordArray,
    const Matrix<double, Dynamic, 3> &srcCoordArray)
{
    RowVector3d srcCenterCoord = srcCoordArray.colwise().mean();
    RowVector3d tarCenterCoord = tarCoordArray.colwise().mean();

    JacobiSVD<Matrix3d> svd(
        (srcCoordArray.rowwise() - srcCenterCoord).transpose() *
        (tarCoordArray.rowwise() - tarCenterCoord),
        ComputeFullU | ComputeFullV);

    Matrix3d U = svd.matrixU(), V = svd.matrixV().transpose();

    if (U.determinant() * V.determinant() < 0.)
    {
        U.col(2) = -U.col(2);
    }

    auto rotationMatrix = U * V;

    // 对于一个向量v，其叠合后的向量为：(v - srcCenterCoord) * rotationMatrix + tarCenterCoord
    return {srcCenterCoord, rotationMatrix, tarCenterCoord};
}
```
