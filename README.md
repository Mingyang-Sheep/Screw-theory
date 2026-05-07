# Screw Algebra & Lie Group Derivation Skill

> A Claude Code skill for solving problems in screw algebra, Lie group kinematics, and mechanism analysis, based on the works of **Academician Jian S. Dai**.

> **Reference Resources**: All papers, textbooks, and lecture slides referenced in this skill are uploaded to [Release - Reference](https://github.com/Mingyang-Sheep/screw-theory-deriver/releases/tag/Reference).

[中文版](README_CN.md)

---

## Overview

This skill provides a disciplined reasoning framework for deriving and verifying results in:

- **Plücker coordinates** — ray order vs axis order, elliptic polar operator
- **Screw algebra** — twist, wrench, reciprocal product, screw systems (Gibson-Hunt classification)
- **SE(3)/se(3) Lie group & algebra** — Euler-Rodrigues formula, exponential map, Chasles' decomposition, adjoint representation, POE modeling
- **Null space construction** — cofactor method, Ball theorem, singularity analysis

The skill enforces **mandatory SymPy symbolic computation** for all non-trivial algebraic operations, ensuring precision and reproducibility.

---

## Repository Structure

```
├── skill.md                    # Skill specification (English)
├── skill_CN.md                 # Skill specification (中文)
├── README.md                   # This file (English)
├── README_CN.md                # Project description (中文)
└── Refer/                      # Reference materials
    ├── Books/                  # Textbooks
    │   ├── Dai_2014_Screw_Algebra_and_Lie_Groups.pdf
    │   └── Dai_2014_Geometric_Foundation_of_Mechanisms_and_Robotics.pdf
    ├── Papers/                 # Journal & conference papers
    │   ├── Dai_Sun_2020_Plucker_Ray_Axis_Order.pdf
    │   ├── Dai_2012_Finite_Displacement_Screw_Operator.pdf
    │   ├── Dai_2015_Euler_Rodrigues_Quaternion_Conjugation.pdf
    │   ├── Dai_Jones_2002_Null_Space_Cofactor.pdf
    │   ├── Dai_2005_Rodrigues_to_Finite_Twist_Review.pdf
    │   ├── Dai_Jones_2006_Metamorphic_Topological_Changes.pdf
    │   └── Dai_Jones_1998_Metamorphic_Mobility.pdf
    └── Lectures/               # Course lecture slides
        └── ...
```

---

## Core Knowledge Modules

The skill covers four interconnected modules:

### Module 1: Plücker Coordinates
- Ray order vs Axis order definitions
- Elliptic polar operator $\boldsymbol{\Delta}$
- Form conformability vs geometry conformability

### Module 2: Screw Algebra
- Twist and wrench (MLS convention)
- Reciprocal product and coordinate invariance
- Gibson-Hunt screw system classification
- Force Jacobian (virtual work principle)
- Lie bracket in se(3)

### Module 3: SE(3) Lie Group / se(3) Lie Algebra
- Euler-Rodrigues rotation formula
- Exponential map and logarithmic map
- Chasles' decomposition (four traces)
- 6x6 finite displacement operator
- Adjoint representation and coordinate transformation
- POE (Product of Exponentials) modeling

### Module 4: Null Space Construction
- Cofactor method (Theorem 2.1)
- Multi-dimensional null space (Theorem 3.2)
- Ball theorem: $\dim(\mathrm{Null}) = 6 - \mathrm{rank}$
- Static vs kinematic singularity analysis

---

## Key References

| Author(s) | Year | Title | Venue |
|-----------|------|-------|-------|
| Ball | 1900 | A Treatise on the Theory of Screws | Cambridge UP |
| Gibson & Hunt | 1990 | Geometry of screw systems | *Mech. Mach. Theory* 25(1), 1-27 |
| Dai & Rees Jones | 1998 | Mobility in Metamorphic Mechanisms | *ASME DETC98/MECH-5902* |
| Dai & Rees Jones | 2002 | Null-space construction using cofactors | *Proc. R. Soc. Lond. A* 458, 1845-1866 |
| Dai & Rees Jones | 2006 | Matrix Representation of Topological Changes | *ASME J. Mech. Des.* 127(4) |
| Dai | 2005 | Historical review: Rodrigues to finite twist | *Mech. Mach. Theory* 41, 143-160 |
| Dai | 2012 | Finite Displacement Screw Operators | *ASME J. Mech. Rob.* 4(4), 041002 |
| Dai | 2014 | 旋量代数与李群、李代数 | Higher Education Press |
| Dai | 2015 | Euler-Rodrigues formula variations | *Mech. Mach. Theory* 92, 144-152 |
| Dai & Sun | 2020 | Plücker Ray/Axis Order | *Mech. Mach. Theory* 153, 103983 |
| Murray, Li & Sastry | 1994 | A Mathematical Introduction to Robotic Manipulation | CRC Press |
| Selig | 2005 | Geometric Fundamentals of Robotics | Springer |

---

## Usage

This skill is designed for use with Claude Code. The skill enforces:

1. **Mandatory SymPy** — all non-trivial algebraic operations must be computed symbolically
2. **Theorem citation** — every derivation step must reference its source
3. **Coordinate order awareness** — explicit declaration and conversion between MLS and Ray orders
4. **Verification** — dimension checks, rank checks, and reciprocity verification after every key computation

---

## Disclaimer

All PDF files are for academic research and study purposes only. Copyright belongs to the original authors and publishers.

---

## Version History

| Version | Date | Changes |
|---------|------|---------|
| v2.0 | 2026-05-07 | Standardized skill structure with bilingual support |
| v1.0 | 2026-05-07 | Initial release |
