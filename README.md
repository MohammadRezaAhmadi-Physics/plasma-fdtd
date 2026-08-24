# Plasma-FDTD: Electromagnetic Wave Propagation in a Dispersive Plasma Slab

## Overview
This repository provides a standalone 2D Finite-Difference Time-Domain (FDTD) solver in Python to simulate the interaction of transverse magnetic ($TM_z$) electromagnetic waves with a dispersive and lossy plasma slab. The plasma medium is modeled using the macroscopic **Drude model**. 

This project is developed from scratch without relying on commercial software, offering a transparent, numerically consistent framework to study wave-plasma interactions across different frequency regimes (below, near, and above the plasma frequency).

## Key Features
*   **Core Algorithm:** 2D staggered Yee grid updating the $E_z$, $H_x$, and $H_y$ field components.
*   **Material Dispersion:** Semi-implicit update method for the Drude model's polarization current ($J_z$), accurately incorporating both plasma frequency ($\omega_p$) and electron collision frequency ($\nu_e$).
*   **Source Implementation:** Total-Field/Scattered-Field (TFSF) formulation for pure plane wave generation, carefully corrected for the Yee grid's numerical dispersion.
*   **Boundary Conditions:** Unsplit Convolutional Perfectly Matched Layer (CPML) for efficient absorption of outgoing waves at the simulation boundaries.
*   **Data Analysis:** Single-frequency Discrete Fourier Transform (DFT) at steady-state to compute Reflection ($R$), Transmission ($T$), and Absorption ($A$) power fractions. The result ($R + T + A = 1$) indicates no major power-balance error and supports the numerical consistency of the solver.
