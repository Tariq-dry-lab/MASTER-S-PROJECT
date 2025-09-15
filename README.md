# Effect Of An Electron Donating Group On Overlap-Based Noncovalent Interactions


### A) Project Overview

This repository contains my Master’s project, completed under the supervision of Dr.Jayashree Nagesh.
The project investigates how attaching an electron donating group (EDG) influences overlap-based noncovalent interactions (NCIs) in a model molecular system.

Noncovalent interactions are weak but crucial forces in biology, playing a role in protein folding, molecular recognition, and stability of biomolecules. While hydrogen bonding and electrostatics are well studied, overlap-based interactions such as n→π* and n→σ* are less explored.

The central question of this project:
 **How do electron donating substituents affect not just the familiar n→π interaction, but also other classes of weak interactions in a model propanaldehyde derivative?**

### B) Background

Covalent bonds: strong (~100 kcal/mol)

Noncovalent interactions: weak (~1–10 kcal/mol)

Hydrogen bond: ~5 kcal/mol

Dispersion: ~1 kcal/mol

Overlap-based NCIs arise from delocalisation of electron density between donor and acceptor orbitals.
Examples include:

H-bond → lone pair (O) → σ*(N–H)

n→π* → lone pair (O) → π*(C=O)

n→σ* → lone pair → σ* orbital

These subtle interactions, though individually weak, collectively stabilize biomolecular structures.

### C) Workflow

#### 1) Model system

Propanaldehyde derivative studied previously (IISER Pune).

Four conformers (P1, P2, N1, N2) considered.

#### 2) Objective 1 – Structural analysis

Extract atomic coordinates.

Compute bond lengths, bond angles, dihedrals.

Identify significant structural differences between conformers.

#### 3) Objective 2 – Weak interaction analysis

Perform Natural Bond Orbital (NBO) calculations.

Extract donor→acceptor interactions (second-order perturbation energies).

Apply filter: donor and acceptor separated by ≥3 bonds.

Compare interaction strengths across conformers.

### D) Key Findings

n→π* interaction is indeed enhanced in P1/P2 conformers.

However, other overlap interactions also contribute significantly (e.g., n→σ*).

Structural differences (bond angle/dihedral changes) correlate strongly with changes in interaction energies.

Total weak interaction energy across conformers: ~3.9–4.5 kcal/mol (excluding n→π)* → showing non-negligible contributions from other classes.
