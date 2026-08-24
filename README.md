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


با انجام این کار، اطلاعات پیش‌نیازها و دستورات دقیق اجرای شبیه‌سازی به مخزن شما اضافه می‌شود. پس از ذخیره این بخش، در مرحله بعدی ساختار فایل‌ها (Project Structure) و نتایج را قرار خواهیم داد.
