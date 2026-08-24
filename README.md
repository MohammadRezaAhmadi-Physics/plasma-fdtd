# 2D FDTD Simulation of Electromagnetic Wave Interaction with Lossy Dispersive Plasma

## Overview
This repository contains a custom, independent 2D Finite-Difference Time-Domain (FDTD) solver written in Python. It is specifically designed to simulate and analyze the interaction of a normally incident Transverse Magnetic ($TM_z$) plane wave with a finite lossy dispersive plasma slab. The plasma dynamics are modeled using the **Drude model** to investigate wave propagation across different frequency regimes relative to the plasma frequency ($\omega_p$).

This project was developed from scratch to computationally solve Maxwell's equations without relying on commercial software, showcasing numerical methods, theoretical electromagnetics, and scientific programming.

## Key Features
- **Core FDTD Engine:** 2D Yee grid implementation for the $TM_z$ mode ($E_z$, $H_x$, $H_y$).
- **Dispersive Medium Modeling:** Semi-implicit current update scheme for the Drude model, accounting for plasma frequency ($\omega_p$) and collisional losses ($\nu_e$).
- **Advanced Boundaries:** Convolutional Perfectly Matched Layers (CPML) to flawlessly for pure plane wave generation.
- **Data Analysis:** Steady-state single-frequency Discrete Fourier Transform (DFT) to calculate Reflectance for pure plane wave generation.
- **Data Analysis:** Steady-state single-frequency Discrete Fourier Transform (DFT) to calculate Reflectance ($R$), Transmittance ($T$), and Absorptance ($A$), numerically validating energy conservation ($R + T + A = 1$).

## Requirements
The simulation is built using standard scientific Python libraries. You can install the required dependencies using:
```bash
pip install -r requirements.txt
**`requirements.txt` includes:**
- `numpy`
- `matplotlib`
- `scipy`

## How to Run (Usage)
Clone the repository and run the main Python script to start the simulation. The code will generate visualizations for the fields and power fractions.

bash
git clone https://github.com/MohammadRezaAhmadi-Physics/your-repo-name.git
cd your-repo-name
python main.py
*(Note: You can adjust the simulation parameters such as $\omega$, $\omega_p$, and $\nu_e$ directly inside the main script).*

## Project Structure
text
├── main.py              # Main execution script and simulation loop
├── fdtd_solver.py       # Core FDTD grid, CPML, and Drude model classes
├── visualization.py     # Scripts for plotting 1D profiles, 2D contours, and bar charts
├── requirements.txt     # Python dependencies
├── README.md            # Project documentation
└── images/              # Simulation output figures

## Simulation Results
The simulation investigates three distinct frequency regimes. The inclusion of the collision frequency ($\nu_e = 0.01$) results in Joule heating, physically leading to wave absorption.

### 1. 2D Field Contours (Wave Propagation)
The 2D spatial contours of the electric field ($E_z$) demonstrate the behavior of the wave as it hits the plasma slab (green box). 
- **Sub-plasma ($\omega < \omega_p$):** Strong attenuation and evanescent wave formation.
- **Near-plasma ($\omega \approx \omega_p$):** Epsilon-Near-Zero (ENZ) behavior with stretched wavelengths.
- **Above-plasma ($\omega > \omega_p$):** Wave propagation with an effective refractive index $n < 1$.

![2D Contours](images/image-(7).png)  
*(Please upload your 2D contour image to the `images` folder and update this path)*

### 2. 1D Spatial Projections
Steady-state profiles of the $E_z$ and $H_y$ fields across the simulation domain. The shaded region denotes the plasma slab, clearly showing total reflection for low frequencies and transmission for high frequencies.

![1D Profile Ez](images/image-(2).png)  
*(Please upload your 1D Ez profile image to the `images` folder and update this path)*

### 3. Power Fraction Analysis (R/T/A)
Reflectance ($R$), Transmittance ($T$), and Absorptance ($A$) calculated via DFT. The maximum absorption occurs near the resonance frequency due to the extended wave-matter interaction time. Numerical stability is confirmed as $R + T + A = 1.0000$ is maintained.

![RTA Bar Chart](images/image-(5).png)  
*(Please upload one or all of your RTA bar chart images to the `images` folder and update this path)*

## Author
**Mohammad Reza Ahmadi**  
Undergraduate Student in Physics  
- [GitHub Profile](https://github.com/MohammadRezaAhmadi-Physics)
- [LinkedIn Profile](https://www.linkedin.com/in/your-linkedin-profile)  
- Email: your.email@example.com


### چند نکته تکمیلی برای شما:
۱. در بخش **How to Run**، به جای `your-repo-name` نام واقعی مخزنی که ساخته‌اید را قرار دهید.
۲. در بخش **Author** در انتهای فایل، لینک لینکدین و ایمیل خود را جایگزین مقادیر پیش‌فرض (`your-linkedin-profile` و `your.email@example.com`) کنید.
۳. در بخش **Simulation Results**، در خطوطی که با `![Name](images/...)` شروع می‌شوند، مسیر و نام عکس‌ها را بر اساس نامی که هنگام آپلود در پوشه `images` گیت‌هاب انتخاب می‌کنید تنظیم کنید تا عکس‌ها به درستی در صفحه اصلی نمایش داده شوند.
