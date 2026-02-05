import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def tff_model_main():
    # Membrane / Fluid Properties
    L_p = 5e-8
    mu0 = 0.001
    k_mu = 0.02
    n = 1.5

    # Pressures
    P_feed = 200e3
    P_retentate = 190e3
    P_permeate = 10e3
    TMP = (P_feed + P_retentate)/2 - P_permeate

    iRT = 2.57e6
    S = 0.0
    k_drop = 10

    # Flow / Module
    A_m = 0.1
    Q_f = 1e-4
    A_c = 0.01
    D = 1e-10

    # Vessel
    V_s = 0.5
    V_d = 0.1
    k_d = 0.1
    beta = 0.1
    alpha = V_d / (V_s + V_d)

    # Fouling
    k_f = 1e-5

    # Membrane resistance
    R_m = 1/(L_p * mu0)

    # ================= INITIAL CONDITIONS =================
    y0 = [
        10.0,           # C_s
        10.0,           # C_d
        V_s + V_d,      # V_r
        0.0             # R_f
    ]

    # ================= TIME =================
    t_span = (0, 10000)
    t_eval = np.linspace(t_span[0], t_span[1], 500)

    # ================= SOLVE =================
    sol = solve_ivp(
        lambda t, y: tff_ode(t, y, {
            'mu0': mu0, 'k_mu': k_mu, 'n': n, 'TMP': TMP, 'iRT': iRT,
            'k_drop': k_drop, 'A_m': A_m, 'Q_f': Q_f, 'A_c': A_c,
            'D': D, 'V_s': V_s, 'k_d': k_d, 'beta': beta,
            'alpha': alpha, 'k_f': k_f, 'R_m': R_m
        }),
        t_span, y0, t_eval=t_eval
    )

    t = sol.t
    C_s, C_d, V_r, R_f = sol.y
    C_r = (1 - alpha)*C_s + alpha*C_d

    J_profile = np.zeros_like(t)
    DeltaP_profile = np.zeros_like(t)

    for i in range(len(t)):
        mu = mu0 * (1 + k_mu * (C_r[i]**n))
        u_f = Q_f*(1-beta)/A_c
        k_m = 0.036 * (D**(1/3)) * (u_f**(4/5))

        J = compute_J(mu, R_f[i], TMP, C_r[i], k_drop, iRT, k_m, R_m)
        J_profile[i] = J
        DeltaP_profile[i] = k_drop * mu * J

    # ================= PLOTS =================
    plt.figure()
    plt.plot(t, V_r, linewidth=2)
    plt.xlabel("Time [s]")
    plt.ylabel("Retentate Volume, V_r [m³]")
    plt.grid()

    plt.figure()
    plt.plot(t, C_r, linewidth=2)
    plt.xlabel("Time [s]")
    plt.ylabel("Overall Retentate Concentration, C_r [kg/m³]")
    plt.grid()

    plt.figure()
    plt.plot(t, DeltaP_profile, linewidth=2)
    plt.xlabel("Time [s]")
    plt.ylabel("Pressure Drop, ΔP [Pa]")
    plt.grid()

    plt.show()


# ================= ODE SYSTEM =================
def tff_ode(t, y, p):
    C_s, C_d, V_r, R_f = y
    C_r = (1 - p['alpha'])*C_s + p['alpha']*C_d

    mu = p['mu0'] * (1 + p['k_mu']*(C_r**p['n']))
    u_f = p['Q_f']*(1 - p['beta']) / p['A_c']
    k_m = 0.036 * (p['D']**(1/3)) * (u_f**(4/5))

    J = compute_J(mu, R_f, p['TMP'], C_r, p['k_drop'], p['iRT'], k_m, p['R_m'])

    dC_s = (p['Q_f']*(1-p['beta'])/p['V_s'])*(C_r - C_s) + p['k_d']*(C_d - C_s)
    dC_d = p['k_d']*(C_s - C_d)
    dV_r = -J * p['A_m']
    dR_f = p['k_f'] * J

    return [dC_s, dC_d, dV_r, dR_f]


# ================= FLUX SOLVER =================
def compute_J(mu, R_f, TMP, C_b, k_drop, iRT, k_m, R_m):
    J_guess = 0.0
    for _ in range(20):
        DeltaP = k_drop * mu * J_guess
        TMP_eff = TMP - DeltaP
        C_w = C_b * np.exp(J_guess / k_m)
        J_new = (TMP_eff - iRT * C_w) / ((R_m + R_f) * mu)
        J_new = max(J_new, 0.0)
        if abs(J_new - J_guess) < 1e-6:
            break
        J_guess = J_new
    return J_guess

if __name__ == "__main__":
    tff_model_main()
