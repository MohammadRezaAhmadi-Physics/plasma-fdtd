
import numpy as np
import matplotlib.pyplot as plt

class FDTD_Solver:

    def __init__(
        self,
        Nx,
        Ny,
        freq,
        theta_deg,
        slab_bounds,
        omega_p,
        nu_e,
        tfsf_margin=20,
        cfl_factor=0.45,
        E0=1.0,
        ramp_cycles=5.0,
        dx=1.0,
        dy=1.0,
        pml_thickness=16,
        pml_order=3,
        pml_sigma_max=1.2,
        pml_kappa_max=5.0,
        pml_alpha_max=0.08,
    ):
       
        self.Nx = int(Nx)
        self.Ny = int(Ny)

        self.freq = float(freq)
        self.omega = self.freq
        self.theta_deg = float(theta_deg)
        self.theta = np.deg2rad(self.theta_deg)

        self.dx = float(dx)
        self.dy = float(dy)

        # normalized units
        self.c0 = 1.0
        self.eps0 = 1.0
        self.mu0 = 1.0

        self.cfl_factor = float(cfl_factor)
        self.dt = self.cfl_factor / (self.c0 * np.sqrt((1.0 / self.dx**2) + (1.0 / self.dy**2)))

        self.E0 = float(E0)
        self.ramp_cycles = float(ramp_cycles)
        self.tfsf_margin = int(tfsf_margin)

#------------------------------------

        self.pml_thickness = int(pml_thickness)
        self.pml_order = int(pml_order)
        self.pml_sigma_max = float(pml_sigma_max)
        self.pml_kappa_max = float(pml_kappa_max)
        self.pml_alpha_max = float(pml_alpha_max)
# --- Monitor Y-range (Avoid Top/Bottom CPML) ---
        self.jm0 = self.pml_thickness + 5
        self.jm1 = self.Ny - self.pml_thickness - 5
#--------------------------------------
        # wave vector from numerical dispersion
        self.k = self._solve_numerical_k(self.omega, self.theta)
        self.kx = self.k * np.cos(self.theta)
        self.ky = self.k * np.sin(self.theta)
        if not np.isfinite(self.k) or self.k <= 0:
            raise ValueError(f"Invalid numerical k: {self.k}")

        self.initialize_grid()
        self._build_cpml_profiles()


        # plasma maps
        self.wp = np.zeros((self.Nx, self.Ny), dtype=float)
        self.nu = np.zeros((self.Nx, self.Ny), dtype=float)
        self.Jz = np.zeros((self.Nx, self.Ny), dtype=float)

        self.update_materials(slab_bounds=slab_bounds, omega_p=omega_p, nu_e=nu_e)

        # TFSF box
        self.i1 = self.tfsf_margin
        self.i2 = self.Nx - self.tfsf_margin
        self.j1 = self.tfsf_margin
        self.j2 = self.Ny - self.tfsf_margin

 #-----------------------------------------------
        # تنظیمات مانیتورهای ثبت زمان (Time-History Monitors)
        # مشخص کردن جایگاه مانیتور در i و بازه‌ی خطی در j
        self.mon_i_ref   = self.i1 - 3        # i برای scattered-field
        self.mon_i_trans = self.i2 - 100       # i برای total-field

        # بازه‌ی خطی در جهت y (مثلا از j1 تا j2 از همان TFSF box)
        self.mon_j1 = self.j1 + 60
        self.mon_j2 = self.j2 - 30

        # لیست‌ها برای ذخیره‌ی آرایه‌های خطی در هر گام زمانی
        self.E_ref_trace   = []   # برای scattered-field line
        self.E_trans_trace = []   # برای transmitted-field line
        self.E_inc_trace   = []   # در صورت نیاز برای incident-field line

        #----------------------------------------------
    # (اضافه کردن متد نمونه‌برداری به کلاس)
    def _sample_fields(self, n):
        """
        ثبت یک برش خطی از میدان Ez در جایگاه‌های مشخص برای هر گام زمانی.
        خروجی به صورت آرایه‌ی ۱بعدی طول (mon_j2-mon_j1+1) خواهد بود.
        """
        # slice along y
        j1, j2 = self.mon_j1, self.mon_j2

        # آرایه‌ی خطی scattered field
        Ez_ref_line = self.Ez[self.mon_i_ref, j1:j2+1]
        # آرایه‌ی خطی transmitted field
        Ez_trans_line = self.Ez[self.mon_i_trans, j1:j2+1]
        # آرایه‌ی خطی incident field (در صورت نیاز)
        # Ez_inc_line = self._get_incident_field_line(j1, j2, n)

        # ذخیره (با copy() برای ایمن بودن از تغییر غیرمنتظره)
        self.E_ref_trace.append(Ez_ref_line.copy())
        self.E_trans_trace.append(Ez_trans_line.copy())
        # self.E_inc_trace.append(Ez_inc_line.copy())

#----------------------------------------------------------


    # انتهای متد step را به این شکل تغییر بده تا نمونه‌برداری انجام شود
    def step(self, n):
        t_e = n * self.dt
        t_h = (n + 0.5) * self.dt

        self._update_H()
        self._tfsf_correct_H(t_e)

        self._update_plasma_current()
        self._update_E()
        self._tfsf_correct_E(t_h)
        
        # فراخوانی تابع جدید نمونه‌برداری
        self._sample_fields(n)

        # این متد را برای محاسبه تبدیل فوریه (FFT) و ضرایب R و T اضافه کن
    


    def initialize_grid(self):
        self.Ez = np.zeros((self.Nx, self.Ny), dtype=float)
        self.Hx = np.zeros((self.Nx, self.Ny - 1), dtype=float)
        self.Hy = np.zeros((self.Nx - 1, self.Ny), dtype=float)

        # CPML memory vars (unsplit)
        self.psi_Ezx_y = np.zeros((self.Nx, self.Ny), dtype=float)      # for dHx/dy term
        self.psi_Ezy_x = np.zeros((self.Nx, self.Ny), dtype=float)      # for dHy/dx term
        self.psi_Hx_z_y = np.zeros((self.Nx, self.Ny - 1), dtype=float) # for dEz/dy in Hx
        self.psi_Hy_z_x = np.zeros((self.Nx - 1, self.Ny), dtype=float) # for dEz/dx in Hy

    def update_materials(self, slab_bounds=None, omega_p=0.0, nu_e=0.0):
        self.wp.fill(0.0)
        self.nu.fill(0.0)

        if slab_bounds is None:
            return

        x0, x1, y0, y1 = slab_bounds
        x0 = max(0, int(x0))
        x1 = min(self.Nx, int(x1))
        y0 = max(0, int(y0))
        y1 = min(self.Ny, int(y1))

        self.wp[x0:x1, y0:y1] = float(omega_p)
        self.nu[x0:x1, y0:y1] = float(nu_e)





    def _poly_profile(self, n, thickness, m, sigma_max, kappa_max, alpha_max, half_shift=False):
        sigma = np.zeros(n, dtype=float)
        kappa = np.ones(n, dtype=float)
        alpha = np.zeros(n, dtype=float)

        for i in range(thickness):
            s = (thickness - i - (0.5 if half_shift else 0.0)) / thickness
            s = max(s, 0.0)
            val_sigma = sigma_max * (s ** m)
            val_kappa = 1.0 + (kappa_max - 1.0) * (s ** m)
            val_alpha = alpha_max * (1.0 - s)

            l = i
            r = n - 1 - i
            if 0 <= l < n:
                sigma[l] = max(sigma[l], val_sigma)
                kappa[l] = max(kappa[l], val_kappa)
                alpha[l] = max(alpha[l], val_alpha)
            if 0 <= r < n:
                sigma[r] = max(sigma[r], val_sigma)
                kappa[r] = max(kappa[r], val_kappa)
                alpha[r] = max(alpha[r], val_alpha)

        return sigma, kappa, alpha




    def _build_cpml_profiles(self):
        # Ez grid
        sx_e, kx_e, ax_e = self._poly_profile(
            self.Nx, self.pml_thickness, self.pml_order,
            self.pml_sigma_max, self.pml_kappa_max, self.pml_alpha_max, half_shift=False
        )
        sy_e, ky_e, ay_e = self._poly_profile(
            self.Ny, self.pml_thickness, self.pml_order,
            self.pml_sigma_max, self.pml_kappa_max, self.pml_alpha_max, half_shift=False
        )

        # Hx grid (Ny-1 half shift in y)
        sy_hx, ky_hx, ay_hx = self._poly_profile(
            self.Ny - 1, self.pml_thickness, self.pml_order,
            self.pml_sigma_max, self.pml_kappa_max, self.pml_alpha_max, half_shift=True
        )

        # Hy grid (Nx-1 half shift in x)
        sx_hy, kx_hy, ax_hy = self._poly_profile(
            self.Nx - 1, self.pml_thickness, self.pml_order,
            self.pml_sigma_max, self.pml_kappa_max, self.pml_alpha_max, half_shift=True
        )

        # store as broadcastable
        self.kx_e = kx_e[:, None]
        self.ky_e = ky_e[None, :]
        self.bx_e = np.exp(-(sx_e / kx_e + ax_e) * self.dt)[:, None]
        self.by_e = np.exp(-(sy_e / ky_e + ay_e) * self.dt)[None, :]
        self.cx_e = (sx_e * (self.bx_e[:, 0] - 1.0)) / (sx_e * kx_e + kx_e**2 * ax_e + 1e-15)
        self.cy_e = (sy_e * (self.by_e[0, :] - 1.0)) / (sy_e * ky_e + ky_e**2 * ay_e + 1e-15)
        self.cx_e = self.cx_e[:, None]
        self.cy_e = self.cy_e[None, :]

        self.ky_hx = ky_hx[None, :]
        self.by_hx = np.exp(-(sy_hx / ky_hx + ay_hx) * self.dt)[None, :]
        cy_hx = (sy_hx * (self.by_hx[0, :] - 1.0)) / (sy_hx * ky_hx + ky_hx**2 * ay_hx + 1e-15)
        self.cy_hx = cy_hx[None, :]

        self.kx_hy = kx_hy[:, None]
        self.bx_hy = np.exp(-(sx_hy / kx_hy + ax_hy) * self.dt)[:, None]
        cx_hy = (sx_hy * (self.bx_hy[:, 0] - 1.0)) / (sx_hy * kx_hy + kx_hy**2 * ax_hy + 1e-15)
        self.cx_hy = cx_hy[:, None]




    def _solve_numerical_k(self, omega, theta):
        # Yee numerical dispersion:
        # sin^2(omega*dt/2)/dt^2 = sin^2(kx*dx/2)/dx^2 + sin^2(ky*dy/2)/dy^2
        target = np.sin(0.5 * omega * self.dt) ** 2 / (self.dt ** 2)

        def f(k):
            kx = k * np.cos(theta)
            ky = k * np.sin(theta)
            rhs = (
                np.sin(0.5 * kx * self.dx) ** 2 / (self.dx ** 2)
                + np.sin(0.5 * ky * self.dy) ** 2 / (self.dy ** 2)
            )
            return rhs - target

        k_low = 1e-12
        k_high = max(5.0 * omega / self.c0 + 1e-6, 10.0)

        f_low = f(k_low)
        f_high = f(k_high)

        n_expand = 0
        while f_low * f_high > 0.0 and n_expand < 60:
            k_high *= 1.5
            f_high = f(k_high)
            n_expand += 1

        if f_low * f_high > 0.0:
            return omega / self.c0

        for _ in range(120):
            k_mid = 0.5 * (k_low + k_high)
            f_mid = f(k_mid)

            if abs(f_mid) < 1e-14:
                return k_mid

            if f_low * f_mid <= 0.0:
                k_high = k_mid
                f_high = f_mid
            else:
                k_low = k_mid
                f_low = f_mid

        return 0.5 * (k_low + k_high)


    def check_numerical_dispersion(self):
        lhs = np.sin(0.5 * self.omega * self.dt) ** 2 / (self.dt ** 2)
        rhs = (
            np.sin(0.5 * self.kx * self.dx) ** 2 / (self.dx ** 2)
            + np.sin(0.5 * self.ky * self.dy) ** 2 / (self.dy ** 2)
        )
        return lhs, rhs, abs(rhs - lhs)



    def _source_envelope(self, t):
        t_ramp = self.ramp_cycles * 2.0 * np.pi / max(self.omega, 1e-12)
        if t <= 0.0:
            return 0.0
        if t < t_ramp:
            return np.sin(0.5 * np.pi * t / t_ramp) ** 2
        return 1.0

    def _Ez_inc(self, X, Y, t):
        phase = self.kx * X + self.ky * Y - self.omega * t
        return self.E0 * self._source_envelope(t) * np.sin(phase)

    def _Hx_inc(self, X, Y, t):
        ez = self._Ez_inc(X, Y, t)
        factor = (
            self.dt / (self.mu0 * self.dy)
            * np.sin(0.5 * self.ky * self.dy)
            / np.sin(0.5 * self.omega * self.dt)
        )
        return factor * ez

    def _Hy_inc(self, X, Y, t):
        ez = self._Ez_inc(X, Y, t)
        factor = (
            self.dt / (self.mu0 * self.dx)
            * np.sin(0.5 * self.kx * self.dx)
            / np.sin(0.5 * self.omega * self.dt)
        )
        return -factor * ez



    def _update_H(self):
        # dEz/dy -> Hx
        dEz_dy = (self.Ez[:, 1:] - self.Ez[:, :-1]) / self.dy
        self.psi_Hx_z_y = self.by_hx * self.psi_Hx_z_y + self.cy_hx * dEz_dy
        term_y = dEz_dy / self.ky_hx + self.psi_Hx_z_y
        self.Hx -= (self.dt / self.mu0) * term_y

        # dEz/dx -> Hy
        dEz_dx = (self.Ez[1:, :] - self.Ez[:-1, :]) / self.dx
        self.psi_Hy_z_x = self.bx_hy * self.psi_Hy_z_x + self.cx_hy * dEz_dx
        term_x = dEz_dx / self.kx_hy + self.psi_Hy_z_x
        self.Hy += (self.dt / self.mu0) * term_x

    def _update_plasma_current(self):
        a = (1.0 - 0.5 * self.nu * self.dt) / (1.0 + 0.5 * self.nu * self.dt)
        b = (self.dt * self.wp**2) / (1.0 + 0.5 * self.nu * self.dt)
        self.Jz = a * self.Jz + b * self.Ez

    def _update_E(self):
        dHy_dx = np.zeros_like(self.Ez)
        dHx_dy = np.zeros_like(self.Ez)

        dHy_dx[1:-1, :] = (self.Hy[1:, :] - self.Hy[:-1, :]) / self.dx
        dHx_dy[:, 1:-1] = (self.Hx[:, 1:] - self.Hx[:, :-1]) / self.dy

        self.psi_Ezy_x = self.bx_e * self.psi_Ezy_x + self.cx_e * dHy_dx
        self.psi_Ezx_y = self.by_e * self.psi_Ezx_y + self.cy_e * dHx_dy

        curl_cpml = (dHy_dx / self.kx_e + self.psi_Ezy_x) - (dHx_dy / self.ky_e + self.psi_Ezx_y)

        self.Ez += (self.dt / self.eps0) * (curl_cpml - self.Jz)

    def _tfsf_correct_H(self, t_e):
        i1, i2, j1, j2 = self.i1, self.i2, self.j1, self.j2

        # Hx(i, j) lives at (i*dx, (j+0.5)*dy), but its correction comes from
        # the incident Ez samples on the Ez nodes that enter the Yee difference.
        x_e = np.arange(i1, i2) * self.dx

        y_b_e = j1 * self.dy
        y_t_e = j2 * self.dy

        ez_b = self._Ez_inc(x_e, y_b_e, t_e)
        ez_t = self._Ez_inc(x_e, y_t_e, t_e)

        # Bottom and top TFSF edges for Hx
        self.Hx[i1:i2, j1 - 1] += (self.dt / (self.mu0 * self.dy)) * ez_b
        self.Hx[i1:i2, j2 - 1] -= (self.dt / (self.mu0 * self.dy)) * ez_t

        # Hy(i, j) lives at ((i+0.5)*dx, j*dy), but its correction comes from
        # the incident Ez samples on the Ez nodes that enter the Yee difference.
        y_e = np.arange(j1, j2) * self.dy

        x_l_e = i1 * self.dx
        x_r_e = i2 * self.dx

        ez_l = self._Ez_inc(x_l_e, y_e, t_e)
        ez_r = self._Ez_inc(x_r_e, y_e, t_e)

        # Left and right TFSF edges for Hy
        self.Hy[i1 - 1, j1:j2] -= (self.dt / (self.mu0 * self.dx)) * ez_l
        self.Hy[i2 - 1, j1:j2] += (self.dt / (self.mu0 * self.dx)) * ez_r


    def _tfsf_correct_E(self, t_h):
        i1, i2, j1, j2 = self.i1, self.i2, self.j1, self.j2

        y_e = np.arange(j1, j2) * self.dy
        x_e = np.arange(i1, i2) * self.dx

        # Ez correction uses incident H sampled at the actual H Yee locations.
        # Left / right: Hy at x = (i1-0.5)dx and (i2-0.5)dx, time t_h
        x_l_hy = (i1 - 0.5) * self.dx
        x_r_hy = (i2 - 0.5) * self.dx

        hy_l = self._Hy_inc(x_l_hy, y_e, t_h)
        hy_r = self._Hy_inc(x_r_hy, y_e, t_h)

        self.Ez[i1,     j1:j2] -= (self.dt / (self.eps0 * self.dx)) * hy_l
        self.Ez[i2 - 1, j1:j2] += (self.dt / (self.eps0 * self.dx)) * hy_r

        # Bottom / top: Hx at y = (j1-0.5)dy and (j2-0.5)dy, time t_h
        y_b_hx = (j1 - 0.5) * self.dy
        y_t_hx = (j2 - 0.5) * self.dy

        hx_b = self._Hx_inc(x_e, y_b_hx, t_h)
        hx_t = self._Hx_inc(x_e, y_t_hx, t_h)

        self.Ez[i1:i2, j1]     += (self.dt / (self.eps0 * self.dy)) * hx_b
        self.Ez[i1:i2, j2 - 1] -= (self.dt / (self.eps0 * self.dy)) * hx_t




