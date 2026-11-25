# Finite-Element-Analysis-of-the-2026-Olympic-Ski-Jump

Finite Element static and dynamic structural analysis of the 2026 Olympic ski jump, implemented in MATLAB using 2D Euler-Bernoulli beam elements.

## Project Overview

This repository contains a custom Finite Element Method (FEM) solver built from scratch to model and analyze the structural behavior of a ski jump truss/frame structure. The model consists of **55 nodes** (constrained at the base supports: nodes `1`, `10`, `18`, and `23`) and **65 beam elements** divided into three distinct cross-sectional and material property groups (main vertical columns, bracing elements, and the inclined ski jump ramp deck), plus a lumped vertical spring ($k_y = 5 \times 10^7\text{ N/m}$) at the take-off tip (node `55`).

### Key Analyses Performed
1. **Modal Analysis:** Computation of natural frequencies and mode shapes in the $0\text{–}6\text{ Hz}$ frequency range (extracting the first 4 dominant modes).
2. **Frequency Response Functions (FRF):**
   * Collocated and cross-FRFs (horizontal excitation at node `9` / Point A, vertical response at node `55` / Point B) using Rayleigh proportional damping ($\alpha = 0.1$, $\beta = 2 \times 10^{-4}$).
   * Comparison between the full FEM direct frequency response and the **Modal Superposition** approach (using the first 2 modes).
3. **Static Snow Load Analysis:** Evaluation of static nodal displacements and maximum vertical deflection under a distributed snow load ($p_0 = -7000\text{ N/m}$) acting along the ramp elements (`55` to `65`).
4. **Moving Load Dynamic Simulation:** Time-domain dynamic response (`ode45` state-space integration over $30\text{ s}$) simulating a $90\text{ kg}$ ski jumper accelerating along the curved and straight sections of the ramp (`elements 55–65`), mapping equivalent nodal forces and bending moments onto the FE mesh and tracking the vertical displacement of the take-off tip (Point B).

## Repository Structure

```text
├── docs/               # Project documentation, reports, and problem statements
├── src/                # MATLAB source code and input mesh data
│   ├── main.m          # Main execution script for static, modal, FRF, and transient analyses
│   ├── loadstructure.m # Parser for the input file; builds nodal coordinates, DoF mapping (idb), and incidences
│   ├── assem.m         # Global mass (M) and stiffness (K) matrix assembler
│   ├── el_tra.m        # Local-to-global 2D Euler-Bernoulli beam element mass and stiffness matrices
│   ├── dis_stru.m      # Plots the undeformed structure with boundary condition markers
│   ├── diseg2.m        # Plots deformed shapes and modal eigenvectors using Hermite shape functions
│   ├── findcard.m      # Utility to locate section headers (*NODES, *BEAMS, *PROPERTIES) in input.inp
│   ├── scom.m          # Utility to read uncommented lines (ignoring '!') from input.inp
│   └── input.inp       # Mesh definition (node coordinates, constraints, beam connectivity, and properties)
├── LICENSE             # MIT License
└── README.md           # Project overview and documentation
