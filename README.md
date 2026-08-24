<img width="4532" height="1534" alt="2D_contour_Ez_Theta_0 0" src="https://github.com/user-attachments/assets/33723589-daa8-40e8-b515-5d5b434a29e5" />
<img width="4532" height="1534" alt="2D_contour_Ez_Theta_0 0" src="https://github.com/user-attachments/assets/48faf7c0-6cac-4845-9aee-a4e986dd75e3" />
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


## Requirements
The code is written in Python 3. Standard scientific libraries are required:
* `numpy`
* `scipy`
* `matplotlib`

You can install the dependencies using the provided `requirements.txt` file:
```bash
pip install -r requirements.txt

## How to Run
1. **Clone the repository:**
bash
git clone https://github.com/MohammadRezaAhmadi-Physics/plasma-fdtd.git
cd plasma-fdtd

2. **Run the main simulation:**
bash
python main.py
This will run the FDTD solver for three different frequency regimes (below, near, and above the plasma frequency). It will automatically generate and save the steady-state field plots, contours, and R/T/A bar charts in your working directory.


## Project Structure
* `fdtd_solver.py`: The core simulation engine containing the `FDTD_Solver` class. It handles the Yee grid updates, CPML boundaries, TFSF source, and the Drude model semi-implicit integration.
* `main.py`: The entry point of the simulation. It defines the grid, material properties, and simulation parameters, running the solver across three distinct frequency regimes.
* `visualization.py`: A collection of utility functions for generating contour plots, steady-state field profiles, and computing/plotting the final R/T/A power fractions.

## Simulation Results
The simulation investigates wave propagation at three different frequency regimes relative to the plasma frequency ($\omega_p$):
1. **Below $\omega_p$:** The plasma acts as a highly reflective medium (evanescent wave inside the slab).
2. **Near $\omega_p$:** Transition regime with significant absorption and partial transmission.
3. **Above $\omega_p$:** The plasma becomes largely transparent to the incident wave.

### 1. Steady-State Electric Field Contours
<!-- Drag and drop your contour map image below this line -->


### 2. Field Profiles Across the Slab
<!-- Drag and drop your Ez and Hy line graph images below this line -->


### 3. Power Analysis (Reflection, Transmission, Absorption)
The solver computes the R, T, and A fractions using single-frequency DFT. The conservation of energy ($R + T + A = 1$) is consistently maintained.
<!-- Drag and drop your R/T/A bar chart images below this line -->
