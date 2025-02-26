# Nuclear Lamina Constitutive Model

## Overview
This repository contains the FORTRAN implementation of a microstructure-based constitutive model for the nuclear lamina, as described in the paper "Microstructure-based nuclear lamina constitutive model" by Nima Mostafazadeh and Zhangli Peng (Cytoskeleton, 2024).

The code is implemented as a VUMAT (Vectorized User Material) subroutine for use with Abaqus/Explicit finite element software to simulate the mechanics of the nuclear lamina during nuclear deformation, particularly during transendothelial migration.

## Model Description
The constitutive model is based on the microstructure of the nuclear lamina, accounting for:
1. Lamin dimer force-extension relationship with three distinct phases:
   - Initial coiled-coil stretching
   - Uncoiling of the dimer
   - α-β transition (secondary structure transition)
2. Orientational arrangement of lamin filaments in the network
3. Network density of lamin filaments

The model combines a semiflexible Blundell and Terentjev (BT) model for the coiled state with a modified worm-like chain model for the uncoiled state, connected through a two-state kinetic model for the uncoiling transition.

## Implementation Details
The VUMAT subroutine:

1. Pre-computes a discrete force-extension dataset for the lamin dimer using the `UNFOLD_MODEL` function
2. Creates a 2D lookup table of principal stresses based on principal stretches
3. Uses bilinear interpolation to determine stress values during simulation
4. Incorporates the lipid bilayer mechanics through an area modulus term

## Citation
If you use this code in your research, please cite:
```
Mostafazadeh, N., & Peng, Z. (2024). Microstructure-based nuclear lamina constitutive model. 
Cytoskeleton, 81(8), 297-309. https://doi.org/10.1002/cm.21835
```

## Contact
nmosta5@uic.edu
