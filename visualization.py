import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


def create_animation(histories, REGIMES, SLAB_BOUNDS):
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
    return fig, ani


def plot_contour(histories, REGIMES, SLAB_BOUNDS, THETA_DEG):
    x0, x1, y0, y1 = SLAB_BOUNDS

    fig, axes = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)

    for ax, (label, _) in zip(axes, REGIMES):
        Ez = histories[label]["Ez"][-1]

        vmax = np.max(np.abs(Ez))
        if vmax < 1e-12:
            vmax = 1e-12

        levels = np.linspace(-vmax, vmax, 50)

        cs = ax.contourf(
            Ez.T,
            levels=levels,
            cmap="RdBu_r",
            origin="lower"
        )

        ax.contour(
            Ez.T,
            levels=levels,
            colors="k",
            linewidths=0.25,
            origin="lower"
        )

        ax.add_patch(
            plt.Rectangle(
                (x0, y0),
                x1 - x0,
                y1 - y0,
                fill=False,
                edgecolor="lime",
                linewidth=2
            )
        )

        ax.set_title(label)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        plt.colorbar(cs, ax=ax, shrink=0.8)

    fig.suptitle("Electric Field Contour (Ez)", fontsize=14)
    fig.savefig(f"contour_Ez_Theta_{THETA_DEG}.png", dpi=300, bbox_inches="tight")
    plt.show()


def bilinear_sample(field, x, y, x_origin=0.0, y_origin=0.0, dx=1.0, dy=1.0):
    xi = (x - x_origin) / dx
    yi = (y - y_origin) / dy

    x0 = np.floor(xi).astype(int)
    y0 = np.floor(yi).astype(int)
    x1 = x0 + 1
    y1 = y0 + 1

    x0 = np.clip(x0, 0, field.shape[0] - 1)
    x1 = np.clip(x1, 0, field.shape[0] - 1)
    y0 = np.clip(y0, 0, field.shape[1] - 1)
    y1 = np.clip(y1, 0, field.shape[1] - 1)

    wx = xi - x0
    wy = yi - y0

    f00 = field[x0, y0]
    f10 = field[x1, y0]
    f01 = field[x0, y1]
    f11 = field[x1, y1]

    return (
        (1.0 - wx) * (1.0 - wy) * f00
        + wx * (1.0 - wy) * f10
        + (1.0 - wx) * wy * f01
        + wx * wy * f11
    )


def _line_box_intersections(xc, yc, ux, uy, xmin, xmax, ymin, ymax, eps=1e-14):
    vals = []

    if abs(ux) > eps:
        s = (xmin - xc) / ux
        y = yc + s * uy
        if ymin <= y <= ymax:
            vals.append(s)

        s = (xmax - xc) / ux
        y = yc + s * uy
        if ymin <= y <= ymax:
            vals.append(s)

    if abs(uy) > eps:
        s = (ymin - yc) / uy
        x = xc + s * ux
        if xmin <= x <= xmax:
            vals.append(s)

        s = (ymax - yc) / uy
        x = xc + s * ux
        if xmin <= x <= xmax:
            vals.append(s)

    vals = np.array(vals, dtype=float)
    if vals.size < 2:
        return None

    vals = np.unique(np.round(vals, 12))
    if vals.size < 2:
        return None

    return float(np.min(vals)), float(np.max(vals))


def make_sampling_line(theta_deg, x_center, y_center, xmin, xmax, ymin, ymax, n_samples=900):
    theta = np.deg2rad(theta_deg)
    ux = np.cos(theta)
    uy = np.sin(theta)

    lims = _line_box_intersections(x_center, y_center, ux, uy, xmin, xmax, ymin, ymax)
    if lims is None:
        raise RuntimeError("Could not build propagation-aligned sampling line.")

    smin, smax = lims
    s = np.linspace(smin, smax, n_samples)
    x = x_center + s * ux
    y = y_center + s * uy
    return s, x, y


def slab_interval_on_sampling_line(theta_deg, x_center, y_center, slab_bounds, dx=1.0, dy=1.0):
    theta = np.deg2rad(theta_deg)
    ux = np.cos(theta)
    uy = np.sin(theta)

    i0, i1, j0, j1 = slab_bounds

    xmin = i0 * dx
    xmax = i1 * dx
    ymin = j0 * dy
    ymax = j1 * dy

    lims = _line_box_intersections(x_center, y_center, ux, uy, xmin, xmax, ymin, ymax)
    return lims


def plot_clean_snapshot(histories, REGIMES, SLAB_BOUNDS, THETA_DEG, field_name="Ez", pct=0.9):
    fig, axes = plt.subplots(1, 3, figsize=(18, 5), constrained_layout=True)

    max_val = 0
    for label in REGIMES:
        max_val = max(max_val, np.max(np.abs(histories[label[0]][field_name])))

    for ax, (label, _) in zip(axes, REGIMES):
        dx = histories[label]["dx"]
        dy = histories[label]["dy"]
        x0, x1, y0, y1 = SLAB_BOUNDS
        xc = 0.5 * (x0 + x1) * dx
        yc = 0.5 * (y0 + y1) * dy

        num_frames = len(histories[label][field_name])
        idx = int(pct * (num_frames - 1))
        data = histories[label][field_name][idx]

        s, xs, ys = make_sampling_line(
            histories[label]["theta_deg"],
            xc, yc,
            0, (data.shape[0]-1)*dx,
            0, (data.shape[1]-1)*dy
        )
        x_origin = 0.0 if field_name == "Ez" else 0.5 * dx
        profile = bilinear_sample(data, xs, ys, x_origin=x_origin, y_origin=0.0, dx=dx, dy=dy)

        ax.plot(s, profile, color="#d62728", lw=1.5, label=f"{int(pct*100)}% progress")

        slab_lims = slab_interval_on_sampling_line(
            histories[label]["theta_deg"], xc, yc, SLAB_BOUNDS, dx=dx, dy=dy
        )
        if slab_lims:
            ax.axvspan(slab_lims[0], slab_lims[1], color="gray", alpha=0.15, label="Slab")

        ax.set_ylim(-max_val*1.1, max_val*1.1)
        ax.set_xlim(-220, 100)
        ax.set_title(f"{label} (t ~ {int(pct*100)}%)")
        ax.set_xlabel(r"$s$ (spatial projection)")
        ax.set_ylabel(field_name)
        ax.legend(fontsize=8, loc="upper right")
        ax.grid(alpha=0.3)

    fig.suptitle(f"Steady-State Field Distribution ({field_name})", fontsize=16)
    fig.savefig(f"profile_{field_name}_{int(pct*100)}pct_Theta_{THETA_DEG}.png", dpi=300, bbox_inches="tight")
    plt.show()


def plot_Ez_profile(histories, REGIMES, SLAB_BOUNDS, THETA_DEG, SNAP_EVERY):
    plot_clean_snapshot(histories, REGIMES, SLAB_BOUNDS, THETA_DEG, field_name="Ez")


def plot_Hy_profile(histories, REGIMES, SLAB_BOUNDS, THETA_DEG, SNAP_EVERY):
    plot_clean_snapshot(histories, REGIMES, SLAB_BOUNDS, THETA_DEG, field_name="Hy")







def extract_CW_RTA(time_array, E_inc_time, E_ref_time, E_trans_time, omega_0):
    """
    محاسبه ضرایب R, T, A برای منبع پیوسته (CW) با استفاده از DFT تک‌فرکانس
    """
    dt = time_array[1] - time_array[0]
    period = 2 * np.pi / omega_0
    
    # برای اطمینان از حالت ماندگار، فقط 5 سیکل آخر را تحلیل می‌کنیم
    steady_state_duration = 5 * period
    start_time = time_array[-1] - steady_state_duration
    
    idx_start = np.argmin(np.abs(time_array - start_time))
    
    t_steady = time_array[idx_start:]
    E_inc_steady = E_inc_time[idx_start:]
    E_ref_steady = E_ref_time[idx_start:]
    E_trans_steady = E_trans_time[idx_start:]
    
    # هسته تبدیل فوریه (Fourier Kernel)
    kernel = np.exp(-1j * omega_0 * t_steady)
    
    # محاسبه فازورها (دامنه مختلط)
    E_inc_phasor = np.sum(E_inc_steady * kernel) * dt
    E_ref_phasor = np.sum(E_ref_steady * kernel) * dt
    E_trans_phasor = np.sum(E_trans_steady * kernel) * dt
    
    # محاسبه توان‌ها (مربع دامنه‌ها)
    P_inc = np.abs(E_inc_phasor)**2
    P_ref = np.abs(E_ref_phasor)**2
    P_trans = np.abs(E_trans_phasor)**2
    
    # محاسبه ضرایب
    R = P_ref / P_inc
    T = P_trans / P_inc
    A = 1.0 - (R + T)
    
    return R, T, A

def plot_RTA_barchart(R, T, A, omega_val, theta_val, regime_name="super_plasma"):
    """
    رسم نمودار میله‌ای برای مقادیر R, T, A در یک فرکانس مشخص.
    """
    labels = ['R', 'T', 'A']
    values = [R, T, A]
    
    fig, ax = plt.subplots(figsize=(6, 5))
    
    colors = ['#1f77b4', '#ff7f0e', '#1b9e77']
    bars = ax.bar(labels, values, color=colors, edgecolor='black', width=0.5)
    
    ax.set_ylim(0, 1.2)
    ax.set_ylabel('Power fraction', fontsize=12)
    
    # خط چین افقی در مقدار 1.0 برای بررسی پایستگی انرژی (Energy Conservation)
    ax.axhline(y=1.0, color='gray', linestyle=':', linewidth=1.5, label='R + T + A = 1')
    
    for bar, val in zip(bars, values):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + 0.02,
                f'{val:.4f}', ha='center', va='bottom', fontweight='bold', fontsize=11)
        
    total = R + T + A
    textstr = f'$R + T + A = {total:.4f}$'
    props = dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray')
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=11,
            verticalalignment='top', bbox=props)
    
    ax.set_title(f'R/T/A | {regime_name}\n($\\omega={omega_val:.3f}$, $\\theta={theta_val}^\\circ$)', fontsize=12)
    ax.legend(loc='upper right')
    ax.grid(axis='y', linestyle='-', alpha=0.3)
    
    plt.tight_layout()
    plt.show()
