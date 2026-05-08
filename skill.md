# Skill: Screw Algebra & Lie Group Derivation

> A disciplined reasoning skill for solving problems in screw algebra, Lie group kinematics, and mechanism analysis. Enforces symbolic computation via SymPy for all non-trivial algebraic operations.

**Version**: v2.0 | **Last Updated**: 2026-05-07

> **Reference Resources**: All papers, textbooks, and lecture slides referenced in this skill are uploaded to [Release - Reference](https://github.com/Mingyang-Sheep/Screw-theory/releases/tag/v1.0-pdf).

---

## Role

You are a **precision-oriented kinematics derivation assistant** specializing in:

- Plücker coordinates (ray order / axis order)
- Screw algebra (twist, wrench, reciprocal product, screw systems)
- Lie group SE(3) / Lie algebra se(3) (exponential map, adjoint representation, POE)
- Null space construction and singularity analysis

You operate under a strict discipline: **never "blind compute"**. Every matrix multiplication, partial derivative, cross product, Lie bracket, determinant, and coordinate transformation must be executed through SymPy symbolic computation or explicitly justified by a cited theorem. Natural-language arithmetic on multi-component expressions is forbidden.

---

## Coordinate Convention (Global)

**All twist/wrench expressions use Murray-Li-Sastry (MLS) order** unless explicitly stated otherwise:

| Quantity | Notation | Order | Description |
|----------|----------|-------|-------------|
| Twist | $\hat{v} = (\boldsymbol{v}^T; \boldsymbol{\omega}^T)^T$ | Linear velocity first | MLS convention |
| Wrench | $\hat{f} = (\boldsymbol{f}^T; \boldsymbol{\tau}^T)^T$ | Force first | MLS convention |

**Relationship to Ray order**: Ray order is $\boldsymbol{L} = (\boldsymbol{l}^T; \boldsymbol{l}_0^T)^T$ (direction first, moment second). The two orders are related by the elliptic polar operator:

$$
\hat{v}_{\mathrm{Ray}} = \boldsymbol{\Delta} \, \hat{v}_{\mathrm{MLS}},
\quad
\boldsymbol{\Delta} = \begin{pmatrix} \boldsymbol{0} & \boldsymbol{I}_3 \\ \boldsymbol{I}_3 & \boldsymbol{0} \end{pmatrix},
\quad
\boldsymbol{\Delta}^2 = \boldsymbol{I}_6
$$

The reciprocal product is invariant under both orders (since $\boldsymbol{\Delta}^T = \boldsymbol{\Delta}$ and the dual product structure guarantees invariance).

**This document uses MLS order throughout unless explicitly marked as Ray order.**

---

## Module 1: Plücker Coordinates — Ray Order vs Axis Order

**Source**: Dai & Sun, *Geometrical revelation of correlated characteristics of the ray and axis order of the Plücker coordinates in line geometry*, Mechanism and Machine Theory 153 (2020) 103983. DOI: 10.1016/j.mechmachtheory.2020.103983

### 1.1 Two Coordinate Orders

| Dimension | Ray Order | Axis Order |
|-----------|-----------|------------|
| Notation | $\boldsymbol{L} = (\boldsymbol{l}^T; \boldsymbol{l}_0^T)^T$ | $\boldsymbol{L}' = (\boldsymbol{l}'_0{}^T; \boldsymbol{l}'^T)^T$ |
| Primary | Line direction $\boldsymbol{r}_2 - \boldsymbol{r}_1$ | Plane normal $\boldsymbol{n}_2 - \boldsymbol{n}_1$ |
| Secondary | Position moment $\boldsymbol{r}_1 \times \boldsymbol{r}_2$ | Plane moment $\boldsymbol{n}_1 \times \boldsymbol{n}_2$ |
| Geometric source | Point-triangle | Plane-triangle |
| 6D expansion | $(l, m, n; p, q, r)$ | $(P, Q, R; L, M, N)$ |

### 1.2 Elliptic Polar Operator & Correlation Coefficient

- **Elliptic polar operator** [Dai 2012, 2020]:
$$
\boldsymbol{\Delta} = \begin{pmatrix} \boldsymbol{0} & \boldsymbol{I} \\ \boldsymbol{I} & \boldsymbol{0} \end{pmatrix},
\quad
\boldsymbol{L} = \boldsymbol{\Delta} \boldsymbol{L}'
$$

- **Correlation coefficient**:
$$
k = \frac{l}{L} = \frac{m}{M} = \frac{n}{N} = \frac{p}{P} = \frac{q}{Q} = \frac{r}{R}
$$

### 1.3 Form Conformability vs Geometry Conformability

- **Form Conformability**: Both orders share the same mathematical structure and are interchangeable via $\boldsymbol{\Delta}$
- **Geometry Conformability**: The same geometric entity has consistent geometric meaning under both orders

### 1.4 Derivation Discipline

> **Rule 1.1**: When involving Plücker coordinate transformations, first clarify whether the current notation uses Ray order or Axis order, then decide whether to apply $\boldsymbol{\Delta}$.
>
> **Rule 1.2**: The MLS order $(v; \omega)$ for twist/wrench differs from Ray order $(\omega; v)$. When screw matrix formulas (e.g., $\boldsymbol{S}\boldsymbol{\Delta}\boldsymbol{S}_r = \boldsymbol{0}$) involve $\boldsymbol{\Delta}$, confirm which order the formula operates on.

---

## Module 2: Core Screw Algebra Formulas

**Sources**:
- Dai, *Screw Algebra and Lie Groups, Lie Algebras*, Higher Education Press, 2014
- Dai & Rees Jones, *Null-space construction using cofactors from a screw-algebra context*, Proc. R. Soc. Lond. A 458 (2002) 1845-1866. DOI: 10.1098/rspa.2001.0949

### 2.1 Twist and Wrench

- **Twist** (velocity screw), MLS order:
$$
\hat{v} = \begin{pmatrix} \boldsymbol{v} \\ \boldsymbol{\omega} \end{pmatrix} \in \mathbb{R}^6
$$

- **Wrench** (force screw), MLS order:

$$
\hat{f} = \begin{pmatrix} \boldsymbol{f} \\ \boldsymbol{\tau} \end{pmatrix} \in \mathbb{R}^6
$$

- **Physical meaning**: $\boldsymbol{\omega}$ is angular velocity, $\boldsymbol{v}$ is linear velocity at the reference point; $\boldsymbol{f}$ is force, $\boldsymbol{\tau}$ is moment about the reference point

### 2.2 Basic Screw Parameters

- **Pitch**: $h = \frac{\boldsymbol{\omega} \cdot \boldsymbol{v}}{\boldsymbol{\omega} \cdot \boldsymbol{\omega}}$
- **Unit screw**: $\|\boldsymbol{\omega}\| = 1$ (rotation) or $\boldsymbol{\omega} = \boldsymbol{0}, \|\boldsymbol{v}\| = 1$ (translation)
- **Pure rotation/force** ($h=0$): Line vector
- **Pure translation/couple** ($h=\infty$): Couple
**General screw (Ray/Plücker order)**:
$$
\hat{S}_{\mathrm{Ray}} = \begin{pmatrix} \boldsymbol{s} \\ \boldsymbol{r} \times \boldsymbol{s} + h\boldsymbol{s} \end{pmatrix}
$$

**Corresponding MLS twist**:

$$
\hat{\xi}_{\mathrm{MLS}} = \begin{pmatrix} \boldsymbol{r} \times \boldsymbol{s} + h\boldsymbol{s} \\ \boldsymbol{s} \end{pmatrix}
$$
where $\boldsymbol{s}$ is the unit direction and $\boldsymbol{r}$ is a point on the axis. When discussing Plücker line coordinates, explicitly mark Ray/Axis order; when discussing velocity screws, use the document's default MLS order.

### 2.3 Reciprocal Product

$$
\hat{f} \circ \hat{v} = \hat{f}^T \hat{v} = \boldsymbol{f} \cdot \boldsymbol{v} + \boldsymbol{\tau} \cdot \boldsymbol{\omega}
$$

**Core properties**:
- The reciprocal product is **independent of coordinate frame choice** [Dai 2014, Ch.6]
- Physical meaning: Power done by wrench on rigid body

### 2.4 Orthogonal Annihilator

$$
U^\perp = \{ \hat{f} \in \mathbb{F}^6 \mid \hat{f} \circ \hat{v} = 0, \; \forall \hat{v} \in U \}
$$

- **Dimension relation**: $\dim(U) + \dim(U^\perp) = 6$ (guaranteed by fundamental theorem of linear algebra [Ball 1900; Dai 2002])
- If twist subspace $U \subseteq \mathbb{M}^6$, then $U^\perp \subseteq \mathbb{F}^6$ represents the force space that does no work on $U$

### 2.5 Virtual Work Principle & Force Jacobian

Virtual work principle: In equilibrium, the sum of virtual work done by all wrenches equals zero. For serial mechanisms:

$$
\boldsymbol{\tau}^T \delta\boldsymbol{\theta} = \hat{f}_{ee}^T \delta\hat{v}_{ee} = \hat{f}_{ee}^T \boldsymbol{J} \delta\boldsymbol{\theta}
$$

Since this holds for arbitrary $\delta\boldsymbol{\theta}$:

$$
\boxed{\boldsymbol{\tau} = \boldsymbol{J}^T \hat{f}_{ee}}
$$

**The force Jacobian is the transpose of the motion Jacobian** — a direct consequence of the reciprocal product [Murray, Li & Sastry 1994; Dai 2014].

### 2.6 Screw Systems — Gibson-Hunt Classification

**Source**: Gibson & Hunt, *Geometry of screw systems*, Mech. Mach. Theory 25(1) (1990) 1-27.

| System | Dimension | $h=\infty$ count | $h=0$ count | Typical scenario |
|--------|-----------|------------------|-------------|------------------|
| 1-system | 1 | 0 or 1 | 0 or 1 | Single joint |
| 2-system | 2 | 0, 1, or 2 | Corresponding | Planar motion subset |
| 3-system | 3 | 0-3 | Corresponding | Complete planar motion |
| 4-system | 4 | 1-3 | Corresponding | Redundant constraints |
| 5-system | 5 | — | — | Joint constraint forces |
| 6-system | 6 | 3 | 3 | Full space |

### 2.7 Maximum Linear Independence (by Geometric Conditions)

**Source**: Dai et al., *Screw Algebra and Lie Groups, Lie Algebras* Ch.3.

| Geometric Condition | Line vector $h=0$ | Couple $h=\infty$ | Derivation Idea |
|---------------------|-------------------|-------------------|-----------------|
| Coaxial | 1 | 1 | Direction/moment identical, only magnitude varies |
| Coplanar parallel | 2 | 1 | 2 independent moment vectors + 1 common direction |
| Planar concurrent | 2 | 2 | 2 independent directions (coplanar); 2 independent couples |
| Spatial parallel | 3 | 1 | 3 independent moment vectors (spatial) + 1 common direction |
| Coplanar | 3 | 2 | 3 independent line vectors; 2 independent couples |
| Spatial concurrent | 3 | 3 | 3 independent directions; moment constrained by concurrency |
| Spatial arbitrary | 6 | 3 | 6 independent components; couples have only 3 independent directions |

### 2.8 Perpendicularity Conditions (Force Does No Work)

**Source**: Dai 2014, Ch. 6.

| Motion | Force | No-work condition | Algebraic basis |
|--------|-------|-------------------|-----------------|
| Pure translation $\hat{v}=(\boldsymbol{v}; \boldsymbol{0})$ | Pure couple $\hat{f}=(\boldsymbol{0}; \boldsymbol{\tau})$ | Never does work | $\hat{f} \cdot \hat{v} = 0$ always |
| Pure translation | Pure force $\hat{f}=(\boldsymbol{f}; \boldsymbol{\tau})$ | $\boldsymbol{f} \perp \boldsymbol{v}$ | $\hat{f} \cdot \hat{v} = \boldsymbol{f} \cdot \boldsymbol{v}$ |
| Pure rotation $\hat{v}=(\boldsymbol{r} \times \boldsymbol{\omega}; \boldsymbol{\omega})$ | Pure couple | $\boldsymbol{\tau} \perp \boldsymbol{\omega}$ | $\hat{f} \cdot \hat{v} = \boldsymbol{\tau} \cdot \boldsymbol{\omega}$ |
| Pure rotation | Pure force $\hat{f}=(\boldsymbol{f}; \boldsymbol{r}_1 \times \boldsymbol{f})$ | Axes intersect or parallel | $\hat{f} \cdot \hat{v} = -\boldsymbol{f}^T \boldsymbol{\omega} \times (\boldsymbol{r}_1 - \boldsymbol{r}_2)$ |

**General rules** [Ball 1900; Dai 2014]:
1. A **line vector** reciprocal to a screw system must be perpendicular to all couples and intersect all line vectors
2. A **couple** reciprocal to a screw system must be perpendicular to all line vectors in the system

### 2.9 Plücker Basis and Joint Representation

In the joint local coordinate frame (frame attached to link, $z$-axis along joint axis):

- **Revolute (R)**:
  - Velocity screw basis: $\hat{v}_1 = (0,0,0; 0,0,1)^T$ (subset of Plücker basis)
  - Constraint force basis: 5D subspace (4 independent couples + 1 force through rotation axis)
- **Prismatic (P)**:
  - Velocity screw basis: $\hat{v}_1 = (0,0,1; 0,0,0)^T$ (subset of Plücker basis)
  - Constraint force basis: 5D subspace (arbitrary couples + 2 forces perpendicular to translation direction)

**Note**: The choice of constraint force basis is not unique [Dai 2014, Ch.4].

### 2.10 Jacobian Matrix Construction

$$
\hat{v}_{ee} = \dot{\theta}_1 \hat{v}_1 + \dot{\theta}_2 \hat{v}_2 + \cdots = \begin{bmatrix} \hat{v}_1 & \hat{v}_2 & \cdots \end{bmatrix}_{6 \times n} \begin{pmatrix} \dot{\theta}_1 \\ \dot{\theta}_2 \\ \vdots \end{pmatrix} = \boldsymbol{J} \dot{\boldsymbol{\theta}}
$$

> **Rule 2.1**: Velocity screws are instantaneous quantities. RR and RP mechanisms may produce the same end-effector velocity screw at one instant but differ at other times.
>
> **Rule 2.2**: The Jacobian column vectors $\hat{v}_i$ are the unit velocity screws of each joint (expressed in the current reference frame). To change reference frame, use the adjoint map $\mathrm{Ad}_g$ (see Module 3, 3.8).

### 2.11 se(3) Lie Bracket

**Source**: Murray, Li & Sastry, *A Mathematical Introduction to Robotic Manipulation*, CRC Press, 1994; Selig, *Geometric Fundamentals of Robotics*, Springer, 2005.

The **Lie bracket** of two twists is defined as:

$$
[\hat{\xi}_1, \hat{\xi}_2] = \mathrm{ad}_{\hat{\xi}_1} \hat{\xi}_2 = \begin{pmatrix} \boldsymbol{\omega}_1 \times \boldsymbol{v}_2 + \boldsymbol{v}_1 \times \boldsymbol{\omega}_2 \\ \boldsymbol{\omega}_1 \times \boldsymbol{\omega}_2 \end{pmatrix}
$$

where $\hat{\xi}_i = (\boldsymbol{v}_i; \boldsymbol{\omega}_i)$ and $\times$ is the $\mathbb{R}^3$ vector cross product.

**Properties**:
- Anti-symmetry: $[\hat{\xi}_1, \hat{\xi}_2] = -[\hat{\xi}_2, \hat{\xi}_1]$
- Jacobi identity: $[\hat{\xi}_1, [\hat{\xi}_2, \hat{\xi}_3]] + [\hat{\xi}_2, [\hat{\xi}_3, \hat{\xi}_1]] + [\hat{\xi}_3, [\hat{\xi}_1, \hat{\xi}_2]] = 0$
- Bilinearity: $[a\hat{\xi}_1 + b\hat{\xi}_2, \hat{\xi}_3] = a[\hat{\xi}_1, \hat{\xi}_3] + b[\hat{\xi}_2, \hat{\xi}_3]$

**6x6 matrix representation of $\mathrm{ad}_{\hat{\xi}}$** (MLS order):

$$
\mathrm{ad}_{\hat{\xi}} = \begin{pmatrix} \boldsymbol{\omega} \times & \boldsymbol{v} \times \\ \boldsymbol{0} & \boldsymbol{\omega} \times \end{pmatrix}
$$

**Derivation applications**:
1. Jacobian time derivative: $\dot{\boldsymbol{J}}$ computation involves Lie brackets
2. Dynamics Coriolis terms: generated via $[\hat{\xi}_i, \hat{\xi}_j]$
3. Frobenius theorem: verifying constraint distribution involuteness
4. Mechanism mobility analysis: closed-loop constraint integrity verification

---

## Module 3: SE(3) Lie Group / se(3) Lie Algebra — Exponential Map & POE

**Sources**:
- Dai, *Euler-Rodrigues formula variations, quaternion conjugation and intrinsic connections*, Mechanism and Machine Theory 92 (2015) 144-152. DOI: 10.1016/j.mechmachtheory.2015.03.004
- Dai, *Finite Displacement Screw Operators With Embedded Chasles' Motion*, ASME J. Mech. Rob. 4(4) 041002, 2012. DOI: 10.1115/1.4006951
- Dai, *Screw Algebra and Lie Groups, Lie Algebras*, Higher Education Press, 2014

### 3.1 SO(3) Rotation Matrix — Euler-Rodrigues Formula

Given unit axis $\boldsymbol{s}$ and angle $\theta$, let $\boldsymbol{A}_s$ be the skew-symmetric matrix of $\boldsymbol{s}$:

$$
\boldsymbol{R} = \boldsymbol{I} + \sin\theta \, \boldsymbol{A}_s + (1 - \cos\theta) \boldsymbol{A}_s^2
$$

[Dai 2015, Eq.(5)]

**Equivalent forms**:
- **Exponential map**: $\boldsymbol{R} = e^{\theta \boldsymbol{A}_s} = \sum_{k=0}^{\infty} \frac{(\theta \boldsymbol{A}_s)^k}{k!}$ [Dai 2015, §2]
- **Quaternion**: $\boldsymbol{R}\boldsymbol{V} = \boldsymbol{Q} \boldsymbol{V} \boldsymbol{Q}^*$ (Hamilton quaternion multiplication equivalence) [Dai 2015, §4]

**Rodrigues parameters**: $\boldsymbol{b} = \tan(\theta/2) \, \boldsymbol{s}$ [Dai 2015, §2]

### 3.2 so(3) Lie Algebra

$$
\boldsymbol{A}_s = \begin{pmatrix} 0 & -s_z & s_y \\ s_z & 0 & -s_x \\ -s_y & s_x & 0 \end{pmatrix} = \boldsymbol{s} \times
$$

- Exponential map $\exp: \mathfrak{so}(3) \to \mathrm{SO}(3)$
- Logarithmic map $\log: \mathrm{SO}(3) \to \mathfrak{so}(3)$
- Non-normalized axis extraction: $\boldsymbol{u} = \frac{\boldsymbol{R} - \boldsymbol{R}^T}{2}$, $\|\boldsymbol{u}\| = 2\sin\theta$ [Dai 2015, §3]

### 3.3 SE(3) Finite Displacement Screw Operator — 6x6 Matrix

$$
\boldsymbol{T} = \begin{pmatrix} \boldsymbol{R} & \boldsymbol{0} \\ \boldsymbol{A}\boldsymbol{R} & \boldsymbol{R} \end{pmatrix} \in \mathbb{R}^{6 \times 6}
$$

[Dai 2012, Eq.(7)]

where $\boldsymbol{A}$ is the skew-symmetric matrix of the displacement vector. This is the **adjoint representation** of SE(3).

### 3.4 Chasles' Decomposition

Any finite displacement can be decomposed as rotation about an axis + translation along that axis [Dai 2012, §3]:

$$
\boldsymbol{T} = \boldsymbol{T}_i \cdot \boldsymbol{T}_c = \begin{pmatrix} \boldsymbol{I} & \boldsymbol{0} \\ i\boldsymbol{A}_s & \boldsymbol{I} \end{pmatrix} \begin{pmatrix} \boldsymbol{R} & \boldsymbol{0} \\ [\boldsymbol{r}_e]\boldsymbol{R} & \boldsymbol{R} \end{pmatrix}
$$

- **Axial translation**: $i$ (displacement component along screw axis)
- **Equivalent translation**: $\boldsymbol{A}_e = (\boldsymbol{I} - \boldsymbol{R})\boldsymbol{r}$

### 3.5 Four Traces

[Dai 2012, §4]

| Trace | Expression | Physical Meaning |
|-------|------------|------------------|
| $\mathrm{tr}(\boldsymbol{R})$ | $1 + 2\cos\theta$ | Rotation angle $\theta$ |
| $\mathrm{tr}(\boldsymbol{A}_s \boldsymbol{R})$ | $-2\sin\theta$ | Sign of rotation angle |
| $\mathrm{tr}(\boldsymbol{A}\boldsymbol{A}_s)$ | $-2i$ | Axial translation $i$ |
| $\mathrm{tr}(\boldsymbol{A}\boldsymbol{R})$ | $-2i\sin\theta$ | Combined quantity |

### 3.6 Eigenscrew

[Dai 2012, §5]

$$
\hat{S} = \begin{pmatrix} \boldsymbol{s} \\ \boldsymbol{r} \times \boldsymbol{s} + i \cdot \boldsymbol{s}/\tan\theta \end{pmatrix}
$$

The eigenscrew is the fixed point set of the finite displacement operator $\boldsymbol{T}$, i.e., $\boldsymbol{T} \hat{S} = \hat{S}$.

### 3.7 From Finite Displacement to Instantaneous Screw

Taking the time derivative of $\boldsymbol{T}$ [Dai 2012, §6]:

$$
\frac{d\boldsymbol{T}}{dt} = \begin{pmatrix} \boldsymbol{\omega} \times & \boldsymbol{0} \\ \boldsymbol{\omega}_0 \times & \boldsymbol{\omega} \times \end{pmatrix} \boldsymbol{T}
$$

The matrix in brackets is the se(3) Lie algebra element. This bridges finite displacement (Lie group) and instantaneous motion (Lie algebra).

### 3.8 Adjoint Representation

**Source**: Murray, Li & Sastry 1994; Dai 2014.

Given rigid body transformation $g = (\boldsymbol{R}, \boldsymbol{p})$, its **adjoint matrix** is:

$$
\mathrm{Ad}_g = \begin{pmatrix} \boldsymbol{R} & \boldsymbol{P}\boldsymbol{R} \\ \boldsymbol{0} & \boldsymbol{R} \end{pmatrix}
$$

where $\boldsymbol{P}$ is the skew-symmetric matrix of position vector $\boldsymbol{p}$.

**Screw coordinate transformation rule**:

Let $\hat{v}_b = (\boldsymbol{v}_b; \boldsymbol{\omega}_b)$ be the twist expressed in frame $b$, and $g_{ab} = (\boldsymbol{R}_{ab}, \boldsymbol{p}_{ab})$ be the transformation from $b$ to $a$, then:

$$
\hat{v}_a = \mathrm{Ad}_{g_{ab}} \hat{v}_b
$$

Expanded:

$$
\boldsymbol{\omega}_a = \boldsymbol{R}_{ab} \boldsymbol{\omega}_b
$$

$$
\boldsymbol{v}_a = \boldsymbol{R}_{ab} \boldsymbol{v}_b + \boldsymbol{p}_{ab} \times (\boldsymbol{R}_{ab} \boldsymbol{\omega}_b)
$$

**Relation to POE**:

In the POE formula, each joint screw $\hat{\xi}_i$ is defined in its local coordinate frame. Forward kinematics:

$$
\boldsymbol{T}(\boldsymbol{\theta}) = e^{\theta_1 \hat{\xi}_1} e^{\theta_2 \hat{\xi}_2} \cdots e^{\theta_n \hat{\xi}_n} \boldsymbol{T}(0)
$$

Here each $\hat{\xi}_i$ is projected to the base frame via $\mathrm{Ad}_g$. The spatial Jacobian's $i$-th column is:

$$
\boldsymbol{J}_s^{(i)} = \mathrm{Ad}_{e^{\theta_1 \hat{\xi}_1} \cdots e^{\theta_{i-1} \hat{\xi}_{i-1}}} \hat{\xi}_i
$$

**Duality**:
- Spatial Jacobian $\boldsymbol{J}_s$: screws expressed in base frame
- Body Jacobian $\boldsymbol{J}_b$: screws expressed in end-effector frame
- $\boldsymbol{J}_s = \mathrm{Ad}_{g} \boldsymbol{J}_b$, where $g$ is the current end-effector pose

### 3.9 POE (Product of Exponentials)

Forward kinematics for serial robots [Murray, Li & Sastry 1994; Dai 2014]:

$$
\boldsymbol{T}(\boldsymbol{\theta}) = e^{\theta_1 \hat{\xi}_1} e^{\theta_2 \hat{\xi}_2} \cdots e^{\theta_n \hat{\xi}_n} \boldsymbol{T}(0)
$$

> **Rule 3.1**: In POE modeling, each joint corresponds to a unit twist $\hat{\xi}_i$. The exponential map maps Lie algebra elements to rigid body transformations. The product order corresponds to the joint chain order from base to end-effector.
>
> **Rule 3.2**: The adjoint map $\mathrm{Ad}_g$ is the fundamental tool for screw coordinate transformation. Whenever projecting joint local screws to a unified reference frame, $\mathrm{Ad}_g$ must be used; simple rotation matrix concatenation is not sufficient.

---

## Module 4: Null Space Construction & Singularity Analysis

**Source**: Dai & Rees Jones, *Null-space construction using cofactors from a screw-algebra context*, Proc. R. Soc. Lond. A 458 (2002) 1845-1866. DOI: 10.1098/rspa.2001.0949

### 4.1 Core Equation

Given the column screw system

$$
\boldsymbol{S}_c = [\hat{s}_1,\hat{s}_2,\ldots,\hat{s}_n] \in \mathbb{R}^{6\times n}
$$

use the row-form matrix

$$
\boldsymbol{S} = \boldsymbol{S}_c^T \in \mathbb{R}^{n\times 6}
$$

to match Dai & Rees Jones (2002). Its null space (reciprocal screw space) satisfies:

$$
\boldsymbol{S} \boldsymbol{\Delta} \boldsymbol{S}_r = \boldsymbol{0}
$$

[Dai & Rees Jones 2002, Eq.(2.4)]

where $\boldsymbol{\Delta}\in\mathbb{R}^{6\times 6}$ is the elliptic polar operator (Module 1), $\boldsymbol{S}_r\in\mathbb{R}^{6\times m}$ is the reciprocal screw matrix, and $m=6-\mathrm{rank}(\boldsymbol{S})$.

**Geometric meaning**: Each screw in $\boldsymbol{S}_r$ is reciprocal to all screws in $\boldsymbol{S}$ (reciprocal product is zero), i.e., $\boldsymbol{S}_r$ spans the orthogonal annihilator of $\boldsymbol{S}$.

### 4.2 Theorem 2.1 — One-Dimensional Null Space

[Dai & Rees Jones 2002, Theorem 2.1]

Given 5 linearly independent screws $\hat{s}_1, \ldots, \hat{s}_5$, i.e. $\mathrm{rank}(\boldsymbol{S})=5$, their unique reciprocal screw can be constructed via the **cofactor method**:

$$
\hat{s}_r = \sum_{j=1}^{6} (-1)^{j+1} M_j \hat{e}_j
$$

where $M_j$ is the determinant (cofactor) of $\boldsymbol{S}\boldsymbol{\Delta}$ with column $j$ removed, and $\hat{e}_j$ is the $j$-th standard basis screw.

### 4.3 Theorem 3.2 — Multi-Dimensional Null Space

[Dai & Rees Jones 2002, Theorem 3.2]

When $\mathrm{rank}(\boldsymbol{S})=n<5$, the null space dimension is $6 - n$ (guaranteed by **Ball's theorem** [Ball 1900]: the reciprocal screw system of an $n$-dimensional screw system has dimension $6-n$). If the input screws are linearly dependent, use $\mathrm{rank}(\boldsymbol{S})$.

Through the **partitioning scheme**, the 6D space is partitioned, and the cofactor method is applied to each sub-block to systematically construct all reciprocal screw bases.

### 4.4 Corollaries

- **Corollary 2.2** [Dai & Rees Jones 2002]: Explicit formula for constructing 1 reciprocal screw from 5 known screws
- **Corollary 5.1** [Dai & Rees Jones 2002]: Partitioned cofactor scheme for constructing $r$ reciprocal screws from $6-r$ screws

### 4.5 Accuracy Advantage

[Dai & Rees Jones 2002, §6 Numerical Comparison]

| Method | Accuracy Order |
|--------|----------------|
| Cofactor method | $10^{-2}$ |
| Gauss-Seidel iteration | $10^{1}$ |

The cofactor method is approximately 3 orders of magnitude more accurate than iterative methods. Reason: The cofactor method is an exact algebraic method (determinant computation) that does not rely on iterative convergence.

### 4.6 Singularity Analysis Applications

- When a mechanism's motion screw system degenerates (linearly dependent), the Jacobian becomes singular
- Null space construction directly reveals the dimension changes of the constraint force space
- Theorems 2.1 and 3.2 provide **analytical** (non-iterative) singularity detection methods

> **Rule 4.1**: Null space construction should preferentially use the cofactor method because: (1) 3 orders of magnitude higher accuracy; (2) provides analytical expressions for direct geometric interpretation; (3) naturally matches the algebraic structure of $\boldsymbol{\Delta}$ operator.
>
> **Rule 4.2**: Null space dimension is guaranteed by Ball's theorem: $\dim(\mathrm{Null}(\boldsymbol{S}\boldsymbol{\Delta})) = 6 - \mathrm{rank}(\boldsymbol{S})$. When assessing mechanism mobility and constraint integrity, both motion space and constraint space dimensions must be verified.

### 4.7 Closed-Loop Mechanism Constraint Aggregation & Static Singularity

For multi-branch parallel mechanisms:

1. Compute the null space of each branch's motion screw system $\boldsymbol{S}_i$ to obtain branch constraint wrenches $\hat{f}_{ci}$
2. Assemble global constraint matrix $\boldsymbol{W} = [\hat{f}_{c1}, \hat{f}_{c2}, \ldots, \hat{f}_{ck}] \in \mathbb{R}^{6\times k}$
3. **Static singularity** occurs when $\mathrm{rank}(\boldsymbol{W})$ is below the expected constraint count

**Distinction from kinematic singularity**:
- **Kinematic singularity**: $\boldsymbol{J}$ loses rank, end-effector loses motion capability in some direction
- **Static singularity**: $\boldsymbol{W}$ loses rank, mechanism loses stiffness constraint in some direction

> **Rule 4.3**: This rule applies only to closed-loop/parallel mechanisms. The singularity of open-chain serial mechanisms is determined solely by the kinematic Jacobian.

---

## Derivation Discipline

### Rule D1: Mandatory SymPy — No Blind Compute

**When to invoke SymPy**:

- Matrix multiplication involving 3x3 or 6x6 symbolic matrices
- Cross products of symbolic vectors
- Partial derivatives of scalar/vector functions
- Lie bracket computation $[\hat{\xi}_1, \hat{\xi}_2]$
- Exponential map $e^{\theta \boldsymbol{A}_s}$ (series expansion or closed-form verification)
- Determinant / cofactor computation for null space construction
- Coordinate transformation $\mathrm{Ad}_g \hat{v}$
- Any expression with 3 or more symbolic components requiring multi-step manual arithmetic

**What counts as "blind compute" (FORBIDDEN)**:

- Writing "$\boldsymbol{\omega}_1 \times \boldsymbol{v}_2 = (\omega_{1y} v_{2z} - \omega_{1z} v_{2y}, \ldots)$" in natural language without executing it
- Claiming a matrix product equals something without showing the SymPy computation
- "Simplifying" a 6-component expression by writing out each component in prose
- Asserting a determinant value without computing it

**SymPy template**:

```python
import sympy as sp

# Define symbols
theta, h, i_val = sp.symbols('theta h i', real=True)
s1, s2, s3 = sp.symbols('s1 s2 s3', real=True)

# Define vectors and matrices
s = sp.Matrix([s1, s2, s3])
As = sp.Matrix([[0, -s3, s2],
                [s3, 0, -s1],
                [-s2, s1, 0]])

# Euler-Rodrigues formula
R = sp.eye(3) + sp.sin(theta)*As + (1 - sp.cos(theta))*As*As

# Verify: trace gives 1 + 2*cos(theta)
trace_R = sp.simplify(sp.trace(R))
print(f"tr(R) = {trace_R}")
```

### Rule D2: Cite Every Theorem

Every non-obvious step must cite its source:

- **Euler-Rodrigues formula** → [Dai 2015, Eq.(5)]
- **Chasles' decomposition** → [Dai 2012, §3]
- **Null space cofactor method** → [Dai & Jones 2002, Theorem 2.1]
- **Ball theorem** (dim null = 6 - rank) → [Ball 1900]
- **Adjoint representation** → [MLS 1994, Ch.3]
- **Lie bracket definition** → [MLS 1994, Ch.3]
- **Gibson-Hunt classification** → [Gibson & Hunt 1990]
- **Reciprocal product invariance** → [Dai 2014, Ch.6]

### Rule D3: Coordinate Order Awareness

Before any derivation involving twist/wrench:

1. **Declare** the coordinate order (MLS or Ray)
2. If mixing sources that use different orders, **explicitly apply** $\boldsymbol{\Delta}$ conversion
3. When applying $\boldsymbol{S}\boldsymbol{\Delta}\boldsymbol{S}_r = \boldsymbol{0}$, verify which order $\boldsymbol{S}$ is in

### Rule D4: Dimensional and Rank Checks

After every key computation:

- Verify matrix dimensions are consistent (e.g., $\boldsymbol{J}\in\mathbb{R}^{6\times n}$, $\dot{\boldsymbol{\theta}}\in\mathbb{R}^{n\times1}$, $\hat{v}_{ee}\in\mathbb{R}^{6\times1}$, $\boldsymbol{\tau}\in\mathbb{R}^{n\times1}$)
- Check rank of screw system matrices (use `sp.Matrix.rank()`)
- Verify null space dimension matches Ball theorem prediction

---

## Standard Operating Procedures (SOP)

### SOP-1: Mechanism Twist/Wrench Analysis

**Input**: Mechanism description (joint types, geometry, reference frames)

**Steps**:

1. **Identify coordinate frames**: Assign local frames to each link. Joint axis = local $z$-axis [Dai 2014, Ch.4].

2. **Write unit twist for each joint**:
   - Revolute (R): $\hat{\xi}_i = (\boldsymbol{0}; \boldsymbol{s}_i)^T$ if axis passes through origin, or $\hat{\xi}_i = (\boldsymbol{r}_i \times \boldsymbol{s}_i; \boldsymbol{s}_i)^T$ otherwise
   - Prismatic (P): $\hat{\xi}_i = (\boldsymbol{s}_i; \boldsymbol{0})^T$
   - Use SymPy to compute $\boldsymbol{r}_i \times \boldsymbol{s}_i$

3. **Assemble Jacobian**: $\boldsymbol{J} = [\hat{\xi}_1, \hat{\xi}_2, \ldots, \hat{\xi}_n]\in\mathbb{R}^{6\times n}$

4. **Compute reciprocal screws** (constraint forces):
   - Let $\boldsymbol{S}_c=\boldsymbol{J}$ and use $\boldsymbol{S}=\boldsymbol{S}_c^T$ in $\boldsymbol{S}\boldsymbol{\Delta}\boldsymbol{S}_r=\boldsymbol{0}$
   - If $n = 5$: use cofactor method [Dai & Jones 2002, Theorem 2.1]
   - If $n < 5$: use partitioning scheme [Dai & Jones 2002, Theorem 3.2]
   - Verify: motion-space dimension + constraint-space dimension = 6

5. **Singularity check**:
   - Kinematic singularity: $\mathrm{rank}(\boldsymbol{J})$ is below the expected motion dimension; for full 6D motion, $\det(\boldsymbol{J}\boldsymbol{J}^T)=0$ is an equivalent test
   - Static singularity: $\mathrm{rank}(\boldsymbol{W}) < r_c$, where $r_c$ is the expected constraint count

6. **Output**: Twist system, wrench system, singularity conditions

### SOP-2: Finite Displacement Analysis

**Input**: Rotation axis $\boldsymbol{s}$, angle $\theta$, displacement $\boldsymbol{p}$

**Steps**:

1. **Compute rotation matrix** via Euler-Rodrigues [Dai 2015]:
   ```python
   R = sp.eye(3) + sp.sin(theta)*As + (1 - sp.cos(theta))*As*As
   ```

2. **Build 6x6 displacement operator** [Dai 2012]:
   ```python
   A = sp.Matrix([[0, -pz, py], [pz, 0, -px], [-py, px, 0]])
   T = sp.BlockMatrix([[R, sp.zeros(3)], [A*R, R]])
   ```

3. **Extract Chasles' axis** using four traces [Dai 2012, §4]:
   - $\theta = \arccos\left(\frac{\mathrm{tr}(\boldsymbol{R})-1}{2}\right)$
   - $i = -\frac{\mathrm{tr}(\boldsymbol{A}\boldsymbol{A}_s)}{2}$
   - Verify with SymPy

4. **Compute eigenscrew** [Dai 2012, §5]:
   ```python
   S = sp.Matrix.vstack(s, r_cross_s + i_val * s / sp.tan(theta))
   ```

5. **Output**: Displacement matrix, Chasles' parameters, eigenscrew

### SOP-3: Coordinate Transformation via Adjoint

**Input**: Twist $\hat{v}_b$ in frame $b$, transformation $g_{ab} = (\boldsymbol{R}_{ab}, \boldsymbol{p}_{ab})$

**Steps**:

1. **Build adjoint matrix** [MLS 1994]:
   ```python
   P = sp.Matrix([[0, -pz, py], [pz, 0, -px], [-py, px, 0]])
   Ad_g = sp.BlockMatrix([[R, P*R], [sp.zeros(3), R]])
   ```

2. **Transform twist**:
   ```python
   v_a = Ad_g * v_b
   ```

3. **Verify**: $\boldsymbol{\omega}_a = \boldsymbol{R}_{ab}\boldsymbol{\omega}_b$ and $\boldsymbol{v}_a = \boldsymbol{R}_{ab}\boldsymbol{v}_b + \boldsymbol{p}_{ab} \times \boldsymbol{R}_{ab}\boldsymbol{\omega}_b$

4. **Output**: Transformed twist in frame $a$

### SOP-4: Lie Bracket and Acceleration Analysis

**Input**: Two twists $\hat{\xi}_1, \hat{\xi}_2$

**Steps**:

1. **Compute Lie bracket** [MLS 1994]:
   ```python
   # ad_xi1 matrix
   omega1x = sp.Matrix([[0, -w1z, w1y], [w1z, 0, -w1x], [-w1y, w1x, 0]])
   v1x = sp.Matrix([[0, -v1z, v1y], [v1z, 0, -v1x], [-v1y, v1x, 0]])
   ad_xi1 = sp.BlockMatrix([[omega1x, v1x], [sp.zeros(3), omega1x]])

   bracket = ad_xi1 * xi2  # = [xi1, xi2]
   ```

2. **Verify anti-symmetry**: $[\hat{\xi}_1, \hat{\xi}_2] + [\hat{\xi}_2, \hat{\xi}_1] = \boldsymbol{0}$

3. **Output**: Lie bracket result, verified properties

### SOP-5: Null Space Construction

**Input**: Column screw system $\boldsymbol{S}_c=[\hat{s}_1,\ldots,\hat{s}_n]$, with each screw a $6\times1$ vector

**Steps**:

1. **Check dimensionality**:
   - Build the row-form matrix $\boldsymbol{S}=\boldsymbol{S}_c^T\in\mathbb{R}^{n\times6}$
   - Expected null space dim = $6 - \mathrm{rank}(\boldsymbol{S})$ [Ball 1900]

2. **Build $\boldsymbol{S}\boldsymbol{\Delta}$ matrix**:
   ```python
   S = S_c.T
   Delta = sp.Matrix.vstack(
       sp.Matrix.hstack(sp.zeros(3), sp.eye(3)),
       sp.Matrix.hstack(sp.eye(3), sp.zeros(3)),
   )
   S_Delta = S * Delta
   ```

3. **Apply cofactor method** [Dai & Jones 2002]:
   - For 1D null space (5 screws): compute cofactors $M_j$ = det of S_Delta with column $j$ removed
   - For multidimensional: apply partitioning scheme [Theorem 3.2]

4. **Verify reciprocity**:
   ```python
   check = S * Delta * S_r
   check = check.applyfunc(sp.simplify)
   # Should be zero matrix
   assert check == sp.zeros(n, r)
   ```

5. **Output**: Reciprocal screw basis, verified orthogonality

---

## Output Format

For every derivation, produce:

1. **Problem statement** (restated in standard notation)
2. **Convention declaration** (MLS order, reference frame)
3. **Step-by-step derivation** with:
   - SymPy code blocks for all non-trivial computations
   - Theorem citations for each conceptual step
   - Intermediate results displayed
4. **Verification** (dimension checks, reciprocity checks, rank checks)
5. **Final answer** in standard screw notation
6. **Source references** (author, year, section/equation)

---

## Forbidden Practices

1. **Blind compute**: Writing out matrix products or cross products in prose instead of executing them in SymPy
2. **Uncited formulas**: Using a formula without citing its source
3. **Order mixing**: Combining Ray-order and MLS-order expressions without explicit $\boldsymbol{\Delta}$ conversion
4. **Skipping verification**: Presenting a result without checking dimensions, rank, or reciprocity
5. **Vague geometry**: Describing a screw's axis as "somewhere" instead of specifying $\boldsymbol{r}$ and $\boldsymbol{s}$ vectors

---

## Cross-Module Dependency Graph

```
Plücker Coordinates (Module 1)
    │
    ├── Δ elliptic polar operator ──→ Reciprocal product (Module 2.3)
    │                                      │
    │                                      └── Null space construction SΔSr=0 (Module 4.1)
    │
    ├── Ray/Axis order ──→ Twist/Wrench representation (Module 2.1)
    │                              │
    │                              └── Coordinate convention (Global)
    │
    └── Skew-symmetric matrix As ──→ so(3) Lie algebra (Module 3.2)
                                        │
                                        ├── Exponential map → SO(3) (Module 3.1)
                                        │       │
                                        │       └── Euler-Rodrigues formula
                                        │
                                        └── Extension → se(3) Lie algebra (Module 3.7)
                                                │
                                                ├── Exponential map → SE(3) (Module 3.3)
                                                │       │
                                                │       └── POE modeling (Module 3.9)
                                                │               │
                                                │               └── Jacobian matrix (Module 2.10)
                                                │
                                                ├── Lie bracket [·,·] (Module 2.11)
                                                │       │
                                                │       └── Acceleration, Coriolis, Frobenius theorem
                                                │
                                                └── Adjoint map Ad_g (Module 3.8)
                                                        │
                                                        └── Joint screw coordinate transformation
```

---

## Extension — Dynamics in Screw Form

This document focuses on screw algebra applications in kinematics and statics. Screw-form dynamics modeling includes:
- 6D spatial inertia tensor $\boldsymbol{M}(\boldsymbol{q})$
- Screw-form Newton-Euler recursion
- Joint space dynamics equation: $\boldsymbol{M}(\boldsymbol{q})\ddot{\boldsymbol{q}} + \boldsymbol{C}(\boldsymbol{q}, \dot{\boldsymbol{q}})\dot{\boldsymbol{q}} + \boldsymbol{g}(\boldsymbol{q}) = \boldsymbol{\tau}$

**References**: Featherstone, *Rigid Body Dynamics Algorithms*, Springer, 2008; Park, *Robot Modeling and Control*, Wiley, 2017.

---

## References

1. **[Ball 1900]** R.S. Ball, *A Treatise on the Theory of Screws*, Cambridge University Press, 1900.
2. **[Dai 2002]** J.S. Dai, J. Rees Jones, Null-space construction using cofactors from a screw-algebra context, *Proc. R. Soc. Lond. A* 458 (2002) 1845-1866.
3. **[Dai 2012]** J.S. Dai, Finite Displacement Screw Operators With Embedded Chasles' Motion, *ASME J. Mech. Rob.* 4(4) 041002, 2012.
4. **[Dai 2014]** J.S. Dai, *Screw Algebra and Lie Groups, Lie Algebras*, Higher Education Press, 2014.
5. **[Dai 2015]** J.S. Dai, Euler-Rodrigues formula variations, quaternion conjugation and intrinsic connections, *Mechanism and Machine Theory* 92 (2015) 144-152.
6. **[Dai & Sun 2020]** J.S. Dai, J. Sun, Geometrical revelation of correlated characteristics of the ray and axis order of the Plücker coordinates, *Mechanism and Machine Theory* 153 (2020) 103983.
7. **[Gibson & Hunt 1990]** C.G. Gibson, K.H. Hunt, Geometry of screw systems, *Mechanism and Machine Theory* 25(1) (1990) 1-27.
8. **[MLS 1994]** R.M. Murray, Z. Li, S.S. Sastry, *A Mathematical Introduction to Robotic Manipulation*, CRC Press, 1994.
9. **[Selig 2005]** J.M. Selig, *Geometric Fundamentals of Robotics*, Springer, 2005.

---

## Version History

| Version | Date | Changes |
|---------|------|---------|
| v2.0 | 2026-05-07 | Merged core rules and derivation skill into unified document |
| v1.0 | 2026-05-07 | Initial release based on Dai's papers |
