# 旋量代数与李群推导技能

> 一个 Claude Code 技能，用于解决旋量代数、李群运动学和机构分析中的问题，参考戴建生院士的系列文献。

> **参考资源**：本技能涉及的所有文献、著作及课程讲义已上传至 [Release - Reference](https://github.com/Mingyang-Sheep/screw-theory-deriver/releases/tag/Reference)。

[English Version](README.md)

---

## 概述

本技能提供了一个严格的推理框架，用于推导和验证以下领域的结果：

- **Plücker 坐标** — 射线顺序 vs 轴线顺序、椭圆极算子
- **旋量代数** — 运动旋量、力旋量、互易积、螺旋系（Gibson-Hunt 分类）
- **SE(3)/se(3) 李群与李代数** — Euler-Rodrigues 公式、指数映射、Chasles' 分解、伴随表示、POE 建模
- **零空间构造** — 余因子法、Ball 定理、奇异性分析

本技能要求对所有非平凡代数运算**强制使用 SymPy 符号计算**，确保精度和可重复性。

---

## 仓库结构

```
├── skill.md                    # 技能规范（英文版）
├── skill_CN.md                 # 技能规范（中文版）
├── README.md                   # 项目说明（英文版）
├── README_CN.md                # 本文件（中文版）
└── Refer/                      # 参考资料
    ├── Books/                  # 著作
    │   ├── Dai_2014_Screw_Algebra_and_Lie_Groups.pdf
    │   └── Dai_2014_Geometric_Foundation_of_Mechanisms_and_Robotics.pdf
    ├── Papers/                 # 期刊与会议论文
    │   ├── Dai_Sun_2020_Plucker_Ray_Axis_Order.pdf
    │   ├── Dai_2012_Finite_Displacement_Screw_Operator.pdf
    │   ├── Dai_2015_Euler_Rodrigues_Quaternion_Conjugation.pdf
    │   ├── Dai_Jones_2002_Null_Space_Cofactor.pdf
    │   ├── Dai_2005_Rodrigues_to_Finite_Twist_Review.pdf
    │   ├── Dai_Jones_2006_Metamorphic_Topological_Changes.pdf
    │   └── Dai_Jones_1998_Metamorphic_Mobility.pdf
    └── Lectures/               # 课程讲义
        └── ...
```

---

## 核心知识模块

本技能涵盖四个相互关联的模块：

### 模块一：Plücker 坐标
- 射线顺序与轴线顺序的定义
- 椭圆极算子 $\boldsymbol{\Delta}$
- 形式相容性与几何相容性

### 模块二：旋量代数
- 运动旋量与力旋量（MLS 约定）
- 互易积与坐标不变性
- Gibson-Hunt 螺旋系分类
- 力雅可比（虚功原理）
- se(3) 中的李括号

### 模块三：SE(3) 李群 / se(3) 李代数
- Euler-Rodrigues 旋转公式
- 指数映射与对数映射
- Chasles' 分解（四个迹）
- 6×6 有限位移算子
- 伴随表示与坐标变换
- POE（指数积）建模

### 模块四：零空间构造
- 余因子法（定理 2.1）
- 多维零空间（定理 3.2）
- Ball 定理：$\dim(\mathrm{Null}) = 6 - \mathrm{rank}$
- 静力学奇异与运动学奇异分析

---

## 核心公式速览

### 坐标顺序

本文档默认采用 MLS 顺序：

```math
\hat{v} =
\begin{pmatrix}
\boldsymbol{v} \\
\boldsymbol{\omega}
\end{pmatrix},
\qquad
\hat{f} =
\begin{pmatrix}
\boldsymbol{f} \\
\boldsymbol{\tau}
\end{pmatrix}
\in \mathbb{R}^6
```

其中 $\boldsymbol{v}$ 为参考点处线速度，$\boldsymbol{\omega}$ 为角速度；$\boldsymbol{f}$ 为力，$\boldsymbol{\tau}$ 为关于参考点的力矩。

MLS 顺序与 Ray/Plücker 顺序通过椭圆极算子互换：

```math
\hat{v}_{\mathrm{Ray}} =
\boldsymbol{\Delta}\hat{v}_{\mathrm{MLS}},
\qquad
\boldsymbol{\Delta} =
\begin{pmatrix}
\boldsymbol{0} & \boldsymbol{I}_3 \\
\boldsymbol{I}_3 & \boldsymbol{0}
\end{pmatrix},
\qquad
\boldsymbol{\Delta}^2 = \boldsymbol{I}_6
```

### 螺旋参数与互易积

螺距定义为：

```math
h =
\frac{\boldsymbol{\omega}\cdot\boldsymbol{v}}
{\boldsymbol{\omega}\cdot\boldsymbol{\omega}}
```

Ray/Plücker 顺序下的一般螺旋可写为：

```math
\hat{S}_{\mathrm{Ray}} =
\begin{pmatrix}
\boldsymbol{s} \\
\boldsymbol{r}\times\boldsymbol{s}+h\boldsymbol{s}
\end{pmatrix}
```

对应的 MLS 运动旋量为：

```math
\hat{\xi}_{\mathrm{MLS}} =
\begin{pmatrix}
\boldsymbol{r}\times\boldsymbol{s}+h\boldsymbol{s} \\
\boldsymbol{s}
\end{pmatrix}
```

互易积为：

```math
\hat{f}\circ\hat{v}
=
\boldsymbol{f}\cdot\boldsymbol{v}
+
\boldsymbol{\tau}\cdot\boldsymbol{\omega}
```

### 李括号与伴随变换

对 $\hat{\xi}_i=(\boldsymbol{v}_i;\boldsymbol{\omega}_i)$，MLS 顺序下李括号为：

```math
[\hat{\xi}_1,\hat{\xi}_2]
=
\begin{pmatrix}
\boldsymbol{\omega}_1\times\boldsymbol{v}_2
+
\boldsymbol{v}_1\times\boldsymbol{\omega}_2 \\
\boldsymbol{\omega}_1\times\boldsymbol{\omega}_2
\end{pmatrix}
```

刚体变换 $g_{ab}=(\boldsymbol{R}_{ab},\boldsymbol{p}_{ab})$ 的伴随矩阵为：

```math
\mathrm{Ad}_{g_{ab}} =
\begin{pmatrix}
\boldsymbol{R}_{ab} & [\boldsymbol{p}_{ab}]_{\times}\boldsymbol{R}_{ab} \\
\boldsymbol{0} & \boldsymbol{R}_{ab}
\end{pmatrix}
```

因此：

```math
\boldsymbol{\omega}_a =
\boldsymbol{R}_{ab}\boldsymbol{\omega}_b,
\qquad
\boldsymbol{v}_a =
\boldsymbol{R}_{ab}\boldsymbol{v}_b
+
\boldsymbol{p}_{ab}\times
(\boldsymbol{R}_{ab}\boldsymbol{\omega}_b)
```

### 零空间与奇异性

若 $\boldsymbol{S}_c=[\hat{s}_1,\ldots,\hat{s}_n]\in\mathbb{R}^{6\times n}$ 为列式螺旋系，令 $\boldsymbol{S}=\boldsymbol{S}_c^T\in\mathbb{R}^{n\times6}$。互易旋量矩阵 $\boldsymbol{S}_r\in\mathbb{R}^{6\times m}$ 满足：

```math
\boldsymbol{S}\boldsymbol{\Delta}\boldsymbol{S}_r
=
\boldsymbol{0},
\qquad
m =
6-\mathrm{rank}(\boldsymbol{S})
```

雅可比矩阵通常写为：

```math
\boldsymbol{J}
=
[\hat{\xi}_1,\hat{\xi}_2,\ldots,\hat{\xi}_n]
\in
\mathbb{R}^{6\times n}
```

运动学奇异由 $\mathrm{rank}(\boldsymbol{J})$ 低于预期运动维数判定；静力学奇异由约束矩阵 $\boldsymbol{W}$ 的秩低于预期约束数判定。

---

## 核心文献

| 作者 | 年份 | 标题 | 发表处 |
|------|------|------|--------|
| Ball | 1900 | A Treatise on the Theory of Screws | Cambridge UP |
| Gibson & Hunt | 1990 | Geometry of screw systems | *Mech. Mach. Theory* 25(1), 1-27 |
| Dai & Rees Jones | 1998 | Mobility in Metamorphic Mechanisms | *ASME DETC98/MECH-5902* |
| Dai & Rees Jones | 2002 | Null-space construction using cofactors | *Proc. R. Soc. Lond. A* 458, 1845-1866 |
| Dai & Rees Jones | 2006 | Matrix Representation of Topological Changes | *ASME J. Mech. Des.* 127(4) |
| Dai | 2005 | Historical review: Rodrigues to finite twist | *Mech. Mach. Theory* 41, 143-160 |
| Dai | 2012 | Finite Displacement Screw Operators | *ASME J. Mech. Rob.* 4(4), 041002 |
| 戴建生院士 | 2014 | 旋量代数与李群、李代数 | 高等教育出版社 |
| Dai | 2015 | Euler-Rodrigues formula variations | *Mech. Mach. Theory* 92, 144-152 |
| Dai & Sun | 2020 | Plücker Ray/Axis Order | *Mech. Mach. Theory* 153, 103983 |
| Murray, Li & Sastry | 1994 | A Mathematical Introduction to Robotic Manipulation | CRC Press |
| Selig | 2005 | Geometric Fundamentals of Robotics | Springer |

---

## 使用方法

本技能专为 Claude Code 设计。技能要求：

1. **强制使用 SymPy** — 所有非平凡代数运算必须使用符号计算
2. **定理引用** — 每个推导步骤必须引用其来源
3. **坐标顺序意识** — 明确声明并在 MLS 和 Ray 顺序之间进行转换
4. **验证** — 每次关键计算后进行维度检查、秩检查和互易性验证

---

## 免责声明

所有 PDF 文件仅供学术研究和学习使用。版权归原始作者和出版商所有。

---

## 版本历史

| 版本 | 日期 | 变更 |
|------|------|------|
| v2.0 | 2026-05-07 | 标准化技能结构，支持中英文双语 |
| v1.0 | 2026-05-07 | 初始发布 |
