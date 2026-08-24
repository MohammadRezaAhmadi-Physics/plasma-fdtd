
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
```
## How to Run
1. **Clone the repository:**
```bash
git clone https://github.com/MohammadRezaAhmadi-Physics/plasma-fdtd.git
cd plasma-fdtd
```
2. **Run the main simulation:**
```bash
python main.py
```
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
<img width="4532" height="1534" alt="2D_contour_Ez_Theta_0 0" src="https://github.com/user-attachments/assets/33723589-daa8-40e8-b515-5d5b434a29e5" />
<img width="4524" height="1220" alt="snapshot_timestep_300" src="https://github.com/user-attachments/assets/b7add745-2d25-4b62-8091-4f4e90475455" />
<img width="4524" height="1227" alt="snapshot_timestep_600" src="https://github.com/user-attachments/assets/6dcebd18-2bd1-4bb7-8683-6ba34c10913e" />
<img width="4524" height="1227" alt="snapshot_timestep_1500" src="https://github.com/user-attachments/assets/4814a9ad-522d-4501-8ad0-2d8aa1585848" />

### 2. Field Profiles Across the Slab
<img width="5434" height="1534" alt="1D_profile_Hy_90pct_Theta_0 0" src="https://github.com/user-attachments/assets/2a83a5d2-534c-4656-bb5f-faf26c778c1c" />
<img width="5434" height="1534" alt="1D_profile_Ez_90pct_Theta_0 0" src="https://github.com/user-attachments/assets/099f1d30-73f1-4c69-b3b6-178f0ac702f2" />


### 3. Power Analysis (Reflection, Transmission, Absorption)
The solver computes the R, T, and A fractions using single-frequency DFT. The conservation of energy ($R + T + A = 1$) is consistently maintained.
<img width="1746" height="1424" alt="R_T_A_omega_0 110" src="https://github.com/user-attachments/assets/5b6a61d9-0400-47a5-8007-2b46afe89ac5" />
<img width="1746" height="1424" alt="R_T_A_omega_0 200" src="https://github.com/user-attachments/assets/c9774d01-e0dc-4c7a-b6ed-02f0ab41581c" />
<img width="1746" height="1424" alt="R_T_A_omega_0 400" src="https://github.com/user-attachments/assets/cb7684ac-91ef-46f8-8402-72bdd564c08a" />


## 🚀 Future Work
We plan to expand the capabilities of this FDTD solver in future updates. Planned features include:
* **Oblique Incidence Analysis:** Upgrading the wave source and TFSF boundaries to study the interaction of EM waves with the plasma slab at various oblique incident angles.
* **Magnetized Plasma:** Extending the Drude model to an anisotropic tensor to simulate wave propagation in magnetized plasma (e.g., observing Faraday rotation).
* **3D FDTD Implementation:** Expanding the spatial grid from 2D to 3D for complete volumetric scattering analysis.
* **Performance Optimization:** Implementing GPU acceleration using `CuPy` or `Numba` for faster computation of high-resolution grids.

## 👥 Contributors
* **MohammadReza Ahmadi** - *Lead Developer / Physics Formulation* - [GitHub](https://github.com/MohammadRezaAhmadi-Physics)
* **[نام همکار اول به انگلیسی]** - *[نقش او، مثلاً Co-Developer یا Data Analysis]* - [لینک گیت‌هاب یا لینکدین او]
* **[نام همکار دوم به انگلیسی]** - *[نقش او، مثلاً Code Review یا Physics Consultant]* - [لینک گیت‌هاب یا لینکدین او]

