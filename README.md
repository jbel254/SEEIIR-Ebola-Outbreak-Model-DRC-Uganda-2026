Give all in one readme file, an appropriate title and description of the repo
SEEIIR Ebola Outbreak Model — DRC / Uganda 2026
A from-scratch numerical toolkit for epidemiological modelling, applied to the 2026 Bundibugyo Ebola outbreak in the Democratic Republic of the Congo and Uganda. The project couples a hand-written Runge–Kutta 4 ODE solver with a SEEIIR compartmental model, fits the model to cumulative confirmed cases and deaths, and produces a four-panel diagnostic figure. No SciPy ODE solver is used — the numerics are transparent, inspectable, and written to be read.

https://ebola_seeiir_fit.png/

Table of contents
What this is

The model

The solver

Data

Installation

Usage

Output

Project layout

Fitting procedure

Design notes

Limitations

Roadmap

References

License

What this is
A small, self-contained Python project that:

Implements a SEEIIR compartmental model for Ebola transmission, with two asymptomatic / pre-symptomatic infectious compartments so the model can represent the known epidemiology of Bundibugyo virus.

Solves the model with a hand-written RK4 solver implemented as a subclass of a general ODESolver base class (forward Euler and RK4 are both included for comparison).

Fits the model to surveillance data (cumulative confirmed cases and cumulative deaths) by least squares on a log scale.

Produces a four-panel diagnostic figure: cases, deaths, active compartments, and the fitted beta(t).

The point of the project is not to produce a policy-grade forecast — it is to demonstrate, on a real outbreak, that a small, clearly written ODE solver plus a structured model class is enough to do meaningful epidemic modelling.

The model
Compartments
Symbol	Meaning
S	Susceptible
E1	Exposed, not yet infectious
E2	Exposed, infectious (pre-symptomatic, scaled by r_e2)
I	Symptomatic infectious
Ia	Asymptomatic infectious (scaled by r_ia)
R	Removed (recovered or dead)
C	Cumulative symptomatic cases (auxiliary, non-feedback)
D	Cumulative deaths (auxiliary, non-feedback)
The two auxiliary states C and D are appended so the model's output can be compared directly with surveillance data without post-hoc integration. They do not feed back into the transmission dynamics.

Equations
With N = S + E1 + E2 + I + Ia + R and force of infection

text
force = beta(t) * (I + r_ia * Ia + r_e2 * E2) / N
the system is

text
dS/dt  = -force * S
dE1/dt =  force * S - lambda_1 * E1
dE2/dt =  lambda_1 * (1 - p_a) * E1 - lambda_2 * E2
dI/dt  =  lambda_2 * E2 - mu * I
dIa/dt =  lambda_1 * p_a * E1 - mu * Ia
dR/dt  =  mu * (I + Ia)
dC/dt  =  lambda_2 * E2
dD/dt  =  cfr * mu * I
Parameters
Symbol	Meaning	Default
beta	Transmission rate (may be a function of t)	0.30 /day
r_ia	Relative infectiousness of Ia	0.1
r_e2	Relative infectiousness of E2	0.0
lambda_1	Rate E1 → E2 or E1 → Ia	0.25 /day
lambda_2	Rate E2 → I	0.25 /day
p_a	Fraction of infections that are asymptomatic	0.05
mu	Rate I → R (symptomatic)	0.125 /day
cfr	Case fatality ratio among symptomatic cases	0.47
Default durations implied: mean incubation 1/lambda_1 + 1/lambda_2 = 8 days, mean symptomatic infectious period 1/mu = 8 days.

Basic reproduction number
The next-generation R0, weighting the two branches by p_a, is

text
R0 = beta * [ (1 - p_a) * (r_e2 / lambda_2 + 1 / mu)
            + p_a       * (r_ia / mu) ]
A second method R0_book() reproduces the unweighted formula printed in the book, kept only for side-by-side comparison.

beta(t)
Two time-dependent transmission profiles are provided:

smooth_step_beta(beta_a, beta_b, t_change, width) — a logistic drop (or rise) from beta_a to beta_b, used in the fit.

piecewise_beta(t_starts, values) — exact piecewise-constant beta(t), matching the book's suggestion for modelling interventions.

Fitting uses the smooth step so the least-squares Jacobian stays well-behaved; width → 0 recovers a hard step.

The solver
ODESolver.py implements a small class hierarchy in the style of Section 2.2 of the book:

text
ODESolver            # base class: initial conditions, time loop, advance()
├── ForwardEuler     # first-order, explicit
└── RungeKutta4      # fourth-order, classical RK4
The base class handles initial-condition bookkeeping, array allocation (scalar vs. system), and the time-stepping loop.

Each subclass implements a single advance() method — the only thing that differs between methods.

The right-hand side is always called as f(t, u) and its return value is coerced to a NumPy array, so models may return lists, tuples, or arrays.

RungeKutta4 is used throughout the fitting and plotting in this project. ForwardEuler is included so the two can be compared directly (e.g. on the convergence tests in the book's Chapter 2).

Data
ebola_drc_uganda_template.csv contains cumulative confirmed cases, cumulative deaths, and (where available) recoveries for the 2026 Bundibugyo outbreak in the DRC and Uganda, with source notes for each row.

Column	Meaning
date	Observation date (ISO)
country	DRC or Uganda
cum_confirmed_cases	Cumulative confirmed cases
cum_confirmed_deaths	Cumulative confirmed deaths
cum_recovered	Cumulative recovered
quality	reported or milestone
source_note	Provenance
Caveats:

Early confirmed counts are inflated by testing/lab scale-up, which biases early beta upward. --min-date drops these points; --max-beta caps beta at a defensible upper bound.

Uganda has only 20 confirmed cases — too few to fit on its own. The --country flag lets you fit to DRC only (recommended) or pool.

--population is an assumption, not data. With ~10⁴ cases it barely matters until the epidemic is a sizeable fraction of N.

Installation
Requires Python 3.9+ and the scientific stack.

With conda (recommended)
bash
conda create -n seeiir python=3.10 -y
conda activate seeiir
conda install -y numpy scipy pandas matplotlib
With pip
bash
python -m venv .venv
# Windows:
.venv\Scripts\activate
# macOS / Linux:
source .venv/bin/activate

pip install numpy scipy pandas matplotlib
No ODESolver package is installed from PyPI — ODESolver.py lives in this repository and is imported directly. If a stale odesolver package is present in your environment from earlier experiments, uninstall it with pip uninstall odesolver so it does not shadow the local file on case-insensitive filesystems.

Usage
Fit and plot the DRC outbreak
bash
python seeiir_ebola.py \
    --data ebola_drc_uganda_template.csv \
    --country DRC \
    --min-date 2026-06-02 \
    --max-beta 0.5
This produces ebola_seeiir_fit.png (four panels) and prints the fitted parameters, R0 before and after the fitted change point, the RMS log-residual, and a population-conservation check.

Common options
Flag	Default	Meaning
--data	ebola_drc_uganda_template.csv	CSV with surveillance data
--country	DRC	DRC or Uganda
--population	1.0e7	Assumed effective population at risk
--horizon	60.0	Days to project past the last observation
--reporting	1.0	Fraction of symptomatic cases confirmed
--min-date	(none)	Drop observations before this date
--max-beta	3.0	Upper bound on beta
--out	ebola_seeiir_fit.png	Output figure path
Example: fit through a hard change point
bash
python seeiir_ebola.py \
    --data ebola_drc_uganda_template.csv \
    --country DRC \
    --min-date 2026-05-27 \
    --max-beta 0.5 \
    --out ebola_seeiir_fit_drc.png
Output
The script prints something like:

text
N observations for DRC, t=0 is 2026-05-04

Fitted parameters
  beta_a      0.xxx
  beta_b      0.xxx
  t_change    xx.xxx
  log10_e0     x.xxx
  cfr          0.xxx
  R0 early (beta_a) = x.xx   R0 late (beta_b) = x.xx
  (book's unweighted formula would give x.xx / x.xx)
  RMS log-residual = x.xxx

Saved ebola_seeiir_fit.png
Population conservation error: x.xxe-xx
and writes a four-panel figure:

Panel	Contents
Top-left	Cumulative confirmed cases, model vs. data (log scale)
Top-right	Cumulative deaths, model vs. data
Bottom-left	Active compartments E1, E2, I, Ia
Bottom-right	Fitted beta(t)
The dotted vertical line in the top panels marks the last observation used in the fit; everything to the right is projection.

Project layout
text
SEEIIR_solver/
├── ODESolver.py                      # ODE solver class hierarchy (RK4, Euler)
├── seeiir_ebola.py                   # Model, fitting, plotting, CLI
├── ebola_drc_uganda_template.csv     # Surveillance data
├── ebola_seeiir_fit.png              # Output figure (generated)
└── README.md
Fitting procedure
Fitting is done by scipy.optimize.least_squares on the residual

text
r_i = log1p(model_cases(t_i)) - log1p(data_cases(t_i))     # for each case point
   ⊕ log1p(model_deaths(t_i)) - log1p(data_deaths(t_i))    # for each death point
The model is solved with RungeKutta4 at a fixed internal step of dt = 0.25 days and evaluated at observation times by linear interpolation. Five parameters are fitted:

Parameter	Bounds	Notes
beta_a	[0.01, max_beta]	Early transmission
beta_b	[0.01, max_beta]	Late transmission
t_change	[5, t_max]	Midpoint of the beta drop
log10_e0	[-1, 3]	Initial exposed seed, in log₁₀
cfr	[0.2, 0.9]	Case fatality ratio
x_scale is set per parameter so the optimiser treats beta, t_change, and e0 on comparable footing.

Design notes
No SciPy ODE solver. solve_ivp is deliberately not used. The point is to have a fully explicit RungeKutta4 whose behaviour can be reasoned about from first principles, matching the book's approach.

Model as callable. SEEIIR implements __call__(t, u) so the solver treats it as a plain right-hand side, exactly as in the book's Chapter 5.

beta is a callable. Transmission may be constant, logistic, or piecewise-constant, without changing the model class.

Auxiliary states. C and D are part of the state vector and are integrated by RK4 along with everything else — no post-hoc integration or trapezoidal quadrature.

Conservation check. The script reports the maximum absolute deviation of S + E1 + E2 + I + Ia + R from N, as recommended in the book.

Limitations
The fit is to cumulative counts, so residuals are highly correlated across time. This is fine for parameter estimation but invalidates standard confidence intervals — none are reported.

beta(t) is constrained to a single sigmoid drop. Real outbreaks have multiple interventions; extending to a multi-step beta is a one-line change to smooth_step_beta usage.

r_ia, r_e2, p_a, and the two lambda rates are fixed, not fitted. They should be varied in a sensitivity analysis before drawing conclusions.

The default cfr of 0.47 is a working assumption for Bundibugyo, not a fitted quantity. Treat the fitted cfr with caution: it is identified mainly by the ratio of deaths to cases, and the death series is short.

--population is an assumption. For the current case counts it has almost no effect on the fit; it matters only once the epidemic is a sizeable fraction of N.

Roadmap
Ideas for extending this project:

□ Add a Heun / ExplicitMidpoint subclass to ODESolver.py for convergence-order tests (as in the book's Chapter 2).
□ Add an adaptive-step RungeKuttaFehlberg45 following Chapter 4.
□ Add a backward-Euler / SDIRK solver for the stiff regime (Chapter 3).
□ Fit r_ia, r_e2, and p_a jointly with beta, with priors.
□ Multi-step beta(t) to represent successive interventions.
□ Bootstrap or profile-likelihood confidence intervals on R0.
□ Add an SEIR-with-vaccination variant.
References
The solver hierarchy, the model structure, and the fitting/validation workflow are all adapted from:

Joakim Sundnes, Solving Ordinary Differential Equations in Python, Springer, 2023.
Companion site: https://sundnes.github.io/solving_odes_in_python/

Additional background on the SEEIIR structure and its use during the COVID-19 pandemic:

Kvaerno, A. Singly diagonally implicit Runge–Kutta methods with an explicit first stage. BIT Numerical Mathematics, 44:489–502, 2004.
