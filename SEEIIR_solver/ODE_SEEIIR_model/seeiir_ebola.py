"""SEEIIR model for the 2026 Bundibugyo Ebola outbreak (DRC / Uganda).


Compartments
    S   susceptible
    E1  exposed, not infectious
    E2  exposed, can infect others (scaled by r_e2)
    I   symptomatic infectious
    Ia  asymptomatic infectious (scaled by r_ia)
    R   removed (recovered or dead)

Two auxiliary, non-feedback states are appended so the model can be compared
with surveillance data:
    C   cumulative symptomatic cases   dC/dt = lmbda_2 * E2
    D   cumulative deaths              dD/dt = cfr * mu * I


Usage
    python seeiir_ebola.py --data ebola_drc_uganda_template.csv --country DRC \
        --min-date 2026-06-02 --max-beta 0.5

Notes
    * Early confirmed counts are driven partly by lab/testing scale-up, which
      inflates early beta. --min-date drops those points; --max-beta caps beta.
    * Fit is to cumulative confirmed cases and deaths (log scale). Uganda has
      too few cases (20) to fit on its own.
    * --population is an assumption, not data. With ~1e4 cases it barely
      matters until the epidemic is a sizeable fraction of N.
"""
import argparse

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from ODESolver import RungeKutta4


# --------------------------------------------------------------------------
# Model
# --------------------------------------------------------------------------
class SEEIIR:
    def __init__(self, beta=0.30, r_ia=0.1, r_e2=0.0,
                 lmbda_1=0.25, lmbda_2=0.25, p_a=0.05, mu=0.125, cfr=0.47):
        """beta may be a float or a function beta(t).

        Defaults are working assumptions for Bundibugyo, NOT fitted values:
          1/lmbda_1 + 1/lmbda_2 = 8 days mean incubation
          1/mu                  = 8 days mean symptomatic infectious period
          p_a, r_ia             = poorly known; set p_a = 0 to switch off
          cfr                   = fraction of symptomatic cases that die
        """
        self.beta = beta if callable(beta) else (lambda t, b=beta: b)
        self.r_ia = r_ia
        self.r_e2 = r_e2
        self.lmbda_1 = lmbda_1
        self.lmbda_2 = lmbda_2
        self.p_a = p_a
        self.mu = mu
        self.cfr = cfr

    def __call__(self, t, u):
        S, E1, E2, I, Ia, R, C, D = u
        N = S + E1 + E2 + I + Ia + R
        beta = self.beta(t)
        r_ia, r_e2 = self.r_ia, self.r_e2
        l1, l2, p_a, mu = self.lmbda_1, self.lmbda_2, self.p_a, self.mu

        force = beta * (I + r_ia * Ia + r_e2 * E2) / N
        dS = -force * S
        dE1 = force * S - l1 * E1
        dE2 = l1 * (1 - p_a) * E1 - l2 * E2
        dI = l2 * E2 - mu * I
        dIa = l1 * p_a * E1 - mu * Ia
        dR = mu * (I + Ia)
        dC = l2 * E2
        dD = self.cfr * mu * I
        return [dS, dE1, dE2, dI, dIa, dR, dC, dD]

    def R0(self, beta=None):
        """Next-generation R0 (weights the two branches by p_a)."""
        b = self.beta(0.0) if beta is None else beta
        return b * ((1 - self.p_a) * (self.r_e2 / self.lmbda_2 + 1 / self.mu)
                    + self.p_a * self.r_ia / self.mu)

    def R0_book(self, beta=None):
        """Formula as printed in the book (unweighted); kept for comparison."""
        b = self.beta(0.0) if beta is None else beta
        return (self.r_e2 * b / self.lmbda_2 + self.r_ia * b / self.mu + b / self.mu)


def smooth_step_beta(beta_a, beta_b, t_change, width=3.0):
    """beta drops (or rises) from beta_a to beta_b around t_change.

    A smoothed version of the piecewise-constant beta suggested in the book,
    which keeps the least-squares fit well behaved. Use width -> 0 for a hard step.
    """
    return lambda t: beta_b + (beta_a - beta_b) / (1.0 + np.exp((t - t_change) / width))


def piecewise_beta(t_starts, values):
    """Exact piecewise-constant beta: values[i] applies for t >= t_starts[i]."""
    t_starts, values = np.asarray(t_starts, float), np.asarray(values, float)
    return lambda t: values[max(np.searchsorted(t_starts, t, side="right") - 1, 0)]


def simulate(model, t_end, N_pop, e0, dt=0.25):
    """Run the SEEIIR model with RungeKutta4 from t=0 to t_end (days)."""
    solver = RungeKutta4(model)
    solver.set_initial_condition([N_pop - e0, e0, 0, 0, 0, 0, 0, 0])
    n_steps = int(np.ceil(t_end / dt))
    return solver.solve((0.0, t_end), n_steps)


# --------------------------------------------------------------------------
# Data
# --------------------------------------------------------------------------
def load_data(path, country="DRC", lead_days=14.0, min_date=None):
    """Read the CSV template. t = 0 is `lead_days` before the first observation."""
    df = pd.read_csv(path, parse_dates=["date"])
    df = df[df["country"] == country].dropna(subset=["cum_confirmed_cases"])
    if min_date:
        df = df[df["date"] >= pd.Timestamp(min_date)]
    df = df.sort_values("date").reset_index(drop=True)
    day0 = df["date"].iloc[0] - pd.Timedelta(days=lead_days)
    df["t"] = (df["date"] - day0).dt.days.astype(float)
    return df, day0


# --------------------------------------------------------------------------
# Fitting
# --------------------------------------------------------------------------
FIT_NAMES = ["beta_a", "beta_b", "t_change", "log10_e0", "cfr"]


def build_model(theta, fixed):
    beta_a, beta_b, t_change, _, cfr = theta
    return SEEIIR(beta=smooth_step_beta(beta_a, beta_b, t_change),
                  cfr=cfr, **fixed)


def residuals(theta, df, N_pop, fixed, reporting=1.0):
    model = build_model(theta, fixed)
    t, u = simulate(model, df["t"].max(), N_pop, 10 ** theta[3])
    C = reporting * np.interp(df["t"], t, u[:, 6])
    D = np.interp(df["t"], t, u[:, 7])
    res = list(np.log1p(C) - np.log1p(df["cum_confirmed_cases"]))
    dd = df.dropna(subset=["cum_confirmed_deaths"])
    if len(dd):
        Dm = np.interp(dd["t"], t, u[:, 7])
        res += list(np.log1p(Dm) - np.log1p(dd["cum_confirmed_deaths"]))
    return np.array(res)


def fit(df, N_pop, fixed=None, reporting=1.0, theta0=None, max_beta=3.0):
    fixed = fixed or {}
    theta0 = theta0 or [min(0.45, 0.9 * max_beta), 0.15, 0.5 * df["t"].max(), 0.5, 0.45]
    lo = [0.01, 0.01, 5.0, -1.0, 0.2]
    hi = [max_beta, max_beta, df["t"].max(), 3.0, 0.9]
    sol = least_squares(residuals, theta0, bounds=(lo, hi),
                        args=(df, N_pop, fixed, reporting), x_scale=[0.1, 0.1, 10, 0.5, 0.1])
    return sol


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="ebola_drc_uganda_template.csv")
    ap.add_argument("--country", default="DRC")
    ap.add_argument("--population", type=float, default=1.0e7,
                    help="effective population at risk (ASSUMPTION - set this)")
    ap.add_argument("--horizon", type=float, default=60.0, help="days to project past last datum")
    ap.add_argument("--reporting", type=float, default=1.0, help="fraction of symptomatic cases confirmed")
    ap.add_argument("--min-date", default=None, help="ignore observations before this date (YYYY-MM-DD)")
    ap.add_argument("--max-beta", type=float, default=3.0,
                    help="upper bound on beta; ~0.5 keeps early R0 near 2-3 with default durations")
    ap.add_argument("--out", default="ebola_seeiir_fit.png")
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    df, day0 = load_data(args.data, args.country, min_date=args.min_date)
    print(f"{len(df)} observations for {args.country}, t=0 is {day0.date()}")
    if len(df) < 6:
        print("WARNING: fewer than 6 observations - the fit will be poorly determined.")

    sol = fit(df, args.population, reporting=args.reporting, max_beta=args.max_beta)
    th = sol.x
    print("\nFitted parameters")
    for n, v in zip(FIT_NAMES, th):
        print(f"  {n:10s} {v:8.3f}")
    model = build_model(th, {})
    print(f"  R0 early (beta_a) = {model.R0(th[0]):.2f}   R0 late (beta_b) = {model.R0(th[1]):.2f}")
    print(f"  (book's unweighted formula would give {model.R0_book(th[0]):.2f} / {model.R0_book(th[1]):.2f})")
    rmse = np.sqrt(np.mean(sol.fun ** 2))
    print(f"  RMS log-residual = {rmse:.3f}")

    t_end = df["t"].max() + args.horizon
    t, u = simulate(model, t_end, args.population, 10 ** th[3])
    S, E1, E2, I, Ia, R, C, D = u.T

    # ---- plot
    fig, ax = plt.subplots(2, 2, figsize=(11, 7.5))
    tobs = df["t"]
    ax[0, 0].semilogy(t, args.reporting * C, label="model")
    ax[0, 0].semilogy(tobs, df["cum_confirmed_cases"], "ko", label="confirmed (data)")
    ax[0, 0].axvline(tobs.max(), color="grey", ls=":")
    ax[0, 0].set(title="Cumulative confirmed cases", xlabel="days since t=0")
    ax[0, 0].legend()
    dd = df.dropna(subset=["cum_confirmed_deaths"])
    ax[0, 1].plot(t, D, label="model")
    ax[0, 1].plot(dd["t"], dd["cum_confirmed_deaths"], "ko", label="deaths (data)")
    ax[0, 1].axvline(tobs.max(), color="grey", ls=":")
    ax[0, 1].set(title="Cumulative deaths", xlabel="days since t=0")
    ax[0, 1].legend()
    ax[1, 0].plot(t, E1, label="E1"); ax[1, 0].plot(t, E2, label="E2")
    ax[1, 0].plot(t, I, label="I"); ax[1, 0].plot(t, Ia, label="Ia")
    ax[1, 0].set(title="Active compartments", xlabel="days since t=0"); ax[1, 0].legend()
    ax[1, 1].plot(t, [model.beta(x) for x in t])
    ax[1, 1].set(title="beta(t)", xlabel="days since t=0")
    fig.suptitle(f"SEEIIR / RungeKutta4 fit to {args.country} (N={args.population:.0e})")
    fig.tight_layout()
    fig.savefig(args.out, dpi=130)
    print(f"\nSaved {args.out}")

    # ---- conservation check, as recommended in the book
    Ntot = S + E1 + E2 + I + Ia + R
    print(f"Population conservation error: {np.max(np.abs(Ntot - args.population)):.2e}")


if __name__ == "__main__":
    main()
