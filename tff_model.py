import numpy as np
from scipy.integrate import solve_ivp

# default parameters (not from Excel)
DEFAULT_PARAMS = {
    "L_p": 5e-8,
    "mu0": 0.001,
    "k_mu": 0.02,
    "n": 1.5,
    "iRT": 2.57e6,
    "k_drop": 10,
    "A_m": 0.1,
    "A_c": 0.01,
    "D": 1e-10,
    "k_f": 1e-5,
}

# derived parameter
DEFAULT_PARAMS["R_m"] = 1 / (DEFAULT_PARAMS["L_p"] * DEFAULT_PARAMS["mu0"])


# flux calculation with fixed-point iteration
def compute_J(mu, R_f, TMP, C_b, k_drop, iRT, k_m, R_m):
    J = 0.0
    for _ in range(20):
        DeltaP = k_drop * mu * J
        TMP_eff = TMP - DeltaP
        C_w = C_b * np.exp(J / k_m)
        J_new = (TMP_eff - iRT * C_w) / ((R_m + R_f) * mu)
        J_new = max(J_new, 0.0)
        if abs(J_new - J) < 1e-6:
            break
        J = J_new
    return J


# ODE system definition
def tff_ode(t, y, p):
    C_s, C_d, V_r, R_f = y

    C_r = (1 - p["alpha"]) * C_s + p["alpha"] * C_d

    mu = p["mu0"] * (1 + p["k_mu"] * C_r**p["n"])
    u_f = p["Q_f"] * (1 - p["beta"]) / p["A_c"]
    k_m = 0.036 * p["D"]**(1/3) * u_f**(4/5)

    J = compute_J(mu, R_f, p["TMP"], C_r,
                p["k_drop"], p["iRT"], k_m, p["R_m"])

    dC_s = (p["Q_f"] * (1 - p["beta"]) / p["V_s"]) * (C_r - C_s) \
           + p["k_d"] * (C_d - C_s)

    dC_d = p["k_d"] * (C_s - C_d)
    dV_r = -J * p["A_m"]
    dR_f = p["k_f"] * J

    return [dC_s, dC_d, dV_r, dR_f]


# main function to run the TFF model
def run_tff_model(excel_inputs, t_span, t_eval):
    # merge default parameters with Excel overrides
    params = DEFAULT_PARAMS.copy()
    params.update(excel_inputs)

    # solve ODE system
    sol = solve_ivp(
        lambda t, y: tff_ode(t, y, params),
        t_span,
        excel_inputs["y0"],
        t_eval=t_eval
    )

    return sol
