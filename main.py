THETA_DEG = float(input("please enter an angle in degrees: "))

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from fdtd_solver import FDTD_Solver
from visualization import *

plt.rcParams.update({
    'font.size': 12,
    'axes.titlesize': 14,
    'axes.labelsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10
})

# ---------------- Parameters ----------------
OMEGA_P = 0.2
NU_E = 0.01

NX, NY = 440, 350
SLAB_BOUNDS = (NX // 2 - 50, NX // 2 + 50, 0, NY ) # 100 و280 

TFSF_MARGIN = 30
N_STEPS = 2000
SNAP_EVERY = 50

REGIMES = [
    (r"$\omega < \omega_p$", 0.55 * OMEGA_P),
    (r"$\omega \approx \omega_p$", 1.0 * OMEGA_P),
    (r"$\omega > \omega_p$", 2 * OMEGA_P), #2.2
]

histories = {}

for label, freq in REGIMES:
    print(f"running: {label}, freq={freq:.4f}, theta={THETA_DEG:.1f} deg")

    sim = FDTD_Solver(
        Nx=NX,
        Ny=NY,
        freq=freq,
        theta_deg=THETA_DEG,
        slab_bounds=SLAB_BOUNDS,
        omega_p=OMEGA_P,
        nu_e=NU_E,
        tfsf_margin=TFSF_MARGIN,
        cfl_factor=0.45,
        pml_thickness=30,
        pml_order=3,
        pml_sigma_max=1.2,
        pml_kappa_max=5.0,
        pml_alpha_max=0.08,
        ramp_cycles=7.0,
    )

    lhs, rhs, err = sim.check_numerical_dispersion()
    print(f"dispersion check: lhs={lhs:.12e}, rhs={rhs:.12e}, abs_err={err:.12e}")
    print(f"k={sim.k:.12e}, kx={sim.kx:.12e}, ky={sim.ky:.12e}")

    frames_Ez = []
    frames_Hy = []

    for n in range(N_STEPS):
        sim.step(n)

        if not np.isfinite(sim.Ez).all():
            print(f"NaN/Inf detected at step {n}. stop.")
            break

        if n % SNAP_EVERY == 0:
            frames_Ez.append(sim.Ez.copy())
            frames_Hy.append(sim.Hy.copy())

    histories[label] = {
        "Ez": frames_Ez,
        "Hy": frames_Hy,
        "freq": freq,
        "theta_deg": THETA_DEG,
        "dx": sim.dx,
        "dy": sim.dy,
    }

    print(f"done: {label}, frames={len(frames_Ez)}, max|Ez|={np.max(np.abs(sim.Ez)):.3e}")

#===================================================
for label, freq in REGIMES:
    print(f"\n{'='*50}")
    print(f"Running Regime: {label}, omega={freq:.4f}, theta={THETA_DEG:.1f} deg")
    print(f"{'='*50}")

        # 1. ران مرجع (Reference Run) - پلاسما غیرفعال
    print("-> 1. Running Reference (Vacuum)...")
    sim_ref = FDTD_Solver(
            Nx=NX, Ny=NY, freq=freq , theta_deg=THETA_DEG,
            slab_bounds=SLAB_BOUNDS, omega_p=0.0, nu_e=0.0, # پلاسما صفر
            tfsf_margin=TFSF_MARGIN, ramp_cycles=7.0
        )
    for n in range(N_STEPS):
        sim_ref.step(n)
            
        # میدان فرودی را از مانیتور عبور (که در خلأ همان میدان فرودی است) استخراج می‌کنیم
    E_inc_arr = np.array(sim_ref.E_trans_trace)

    # 2. ران اصلی (Main Run) - پلاسما فعال
    print("-> 2. Running Main Simulation (Plasma)...")
    sim_main = FDTD_Solver(
        Nx=NX, Ny=NY, freq=freq, theta_deg=THETA_DEG,
        slab_bounds=SLAB_BOUNDS, omega_p=OMEGA_P, nu_e=NU_E, # پلاسما فعال
        tfsf_margin=TFSF_MARGIN, ramp_cycles=7.0
    )
    for n in range(N_STEPS):
        sim_main.step(n)

    # میدان‌های بازتابی و عبوری را استخراج می‌کنیم
    E_ref_arr = np.array(sim_main.E_ref_trace)
    E_trans_arr = np.array(sim_main.E_trans_trace)

    #حالا این خروجی هایی ک گرفتیم رو باید تبدیل ب یک بعدی کنیم تا بره تو تابع و داستان دار نشه 

    # 2) انتگرال spatial روی خط: جمع روی محور j و ضرب در Δy
    dy = sim_main.dy    # یا هرجایی که dy تعریف شده
    s_inc_time   = E_inc_arr.sum(axis=1)   * dy   # shape -> (N_steps,)
    s_ref_time   = E_ref_arr.sum(axis=1)   * dy
    s_trans_time = E_trans_arr.sum(axis=1) * dy


    # 3. محاسبه ضرایب و رسم
    print("-> 3. Extracting R, T, A (Steady-State)...")
    dt = sim_main.dt
    time_array = np.arange(N_STEPS) * dt
    
    R, T, A = extract_CW_RTA(
        time_array=time_array,
        E_inc_time=s_inc_time,
        E_ref_time=s_ref_time,
        E_trans_time=s_trans_time,
        omega_0=freq
    )

    print(f"Results for {label}: R={R:.4f}, T={T:.4f}, A={A:.4f} | Sum={R+T+A:.4f}")
    
    # رسم نمودار میله‌ای برای این رژیم
    plot_RTA_barchart(R, T, A, omega_val=freq, theta_val=THETA_DEG, regime_name=label)

#============================================================


# ---------------- Animation ----------------
fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), constrained_layout=True)
ims = []
vmax = max(np.max(np.abs(histories[k]["Ez"])) for k, _ in REGIMES) + 1e-12

for ax, (label, _) in zip(axes, REGIMES):
    ez0 = histories[label]["Ez"][0]
    im = ax.imshow(
        ez0.T, origin="lower", cmap="RdBu_r",
        vmin=-vmax, vmax=vmax, aspect="auto"
    )
    x0, x1, y0, y1 = SLAB_BOUNDS
    ax.plot([x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0], "k--", lw=1.2)

    ax.set_title(label)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ims.append(im)

fig.colorbar(ims[0], ax=axes, shrink=0.8, label="Ez")

nframes = min(len(histories[REGIMES[0][0]]["Ez"]),
              len(histories[REGIMES[1][0]]["Ez"]),
              len(histories[REGIMES[2][0]]["Ez"]))

def animate(it):
    for im, (label, _) in zip(ims, REGIMES):
        im.set_data(histories[label]["Ez"][it].T)
    return ims

ani = animation.FuncAnimation(fig, animate, frames=nframes, interval=45, blit=False, repeat=True)



# ---------------- Run extra plots ----------------
plot_contour(histories, REGIMES, SLAB_BOUNDS, THETA_DEG)
plot_clean_snapshot(histories, REGIMES, SLAB_BOUNDS, THETA_DEG, field_name="Ez")
plot_clean_snapshot(histories, REGIMES, SLAB_BOUNDS, THETA_DEG, field_name="Hy")

    # === اضافه شده برای ذخیره تصاویر در تایم‌استپ‌های دلخواه ===
    # گام‌های زمانی مورد نظر برای ذخیره عکس (مثلاً 400، 700 و 900)
target_timesteps = [300, 600, 1500]  # قبل، حین و بعد از برخورد
save_snapshots_at_timesteps(
    histories,
    REGIMES,
    SLAB_BOUNDS,
    target_timesteps,
    snap_every=SNAP_EVERY
)

plt.show()
