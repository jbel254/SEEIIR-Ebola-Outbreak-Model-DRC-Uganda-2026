# SEEIIR Ebola Outbreak Model — DRC / Uganda 2026

A from-scratch numerical toolkit for epidemiological modelling, applied to the
2026 Bundibugyo Ebola outbreak in the Democratic Republic of the Congo and
Uganda. The project couples a hand-written Runge–Kutta 4 ODE solver with a
SEEIIR compartmental model, fits the model to cumulative confirmed cases and
deaths, and produces a four-panel diagnostic figure.

No SciPy ODE solver is used — the numerics are transparent, inspectable, and
written to be read. The full pipeline runs end-to-end from a public data
source: **fetch → clean → fit → plot**.

![Fit](seeiir_new_fit.png)

---

## Table of contents

1. [What this is](#what-this-is)
2. [Pipeline](#pipeline)
3. [The model](#the-model)
4. [The solver](#the-solver)
5. [Data](#data)
6. [Installation](#installation)
7. [Usage](#usage)
8. [Output](#output)
9. [Project layout](#project-layout)
10. [Fitting procedure](#fitting-procedure)
11. [Design notes](#design-notes)
12. [Limitations](#limitations)
13. [Roadmap](#roadmap)
14. [References](#references)
15. [License](#license)

---

## What this is

A small, self-contained Python project that:

1. **Fetches** the latest surveillance data from a public GitHub tracker
   (INSP situation reports compiled by a third party), with an ECDC page
   cross-check.
2. **Cleans** the raw data: fills missing calendar days by forward-fill
   (correct for cumulative counts), normalizes column names, and enforces
   monotonicity.
3. **Models** the outbreak with a SEEIIR compartmental structure that has
   two asymptomatic / pre-symptomatic infectious compartments, so the
   epidemiology of Bundibugyo virus is represented faithfully.
4. **Solves** the model with a hand-written `RungeKutta4` implemented as a
   subclass of a general `ODESolver` base class.
5. **Fits** five parameters to cumulative confirmed cases and deaths by
   least squares on a log scale.
6. **Plots** a four-panel diagnostic: cases, deaths, active compartments,
   and the fitted `beta(t)`.

The point of the project is not to produce a policy-grade forecast — it is
to demonstrate, on a real outbreak, that a small, clearly written ODE solver
plus a structured model class is enough to do meaningful epidemic modelling.

---

## Pipeline
fetch_data.py clean_data.py seeiir_ebola.py
───────────── ───────────── ───────────────
tracker CSV ──▶ ebola_drc_uganda.csv ──▶ ebola_drc_uganda_clean.csv ──▶ seeiir_new_fit.png
(+ECDC) raw dense daily series fitted figure

text

Each stage is a standalone script and can be run on its own. The middle
file (`ebola_drc_uganda.csv`) is intentionally kept as the raw, unedited
download so that the cleaning step is reproducible and auditable.

---

## The model

### Compartments

| Symbol | Meaning |
|--------|---------|
| `S`    | Susceptible |
| `E1`   | Exposed, **not** yet infectious |
| `E2`   | Exposed, infectious (pre-symptomatic, scaled by `r_e2`) |
| `I`    | Symptomatic infectious |
| `Ia`   | Asymptomatic infectious (scaled by `r_ia`) |
| `R`    | Removed (recovered or dead) |
| `C`    | Cumulative symptomatic cases *(auxiliary, non-feedback)* |
| `D`    | Cumulative deaths *(auxiliary, non-feedback)* |

The two auxiliary states `C` and `D` are appended so the model's output can
be compared directly with surveillance data without post-hoc integration.
They do **not** feed back into the transmission dynamics.

### Equations

With `N = S + E1 + E2 + I + Ia + R` and force of infection
force = beta(t) * (I + r_ia * Ia + r_e2 * E2) / N

text

the system is
dS/dt = -force * S
dE1/dt = force * S - lambda_1 * E1
dE2/dt = lambda_1 * (1 - p_a) * E1 - lambda_2 * E2
dI/dt = lambda_2 * E2 - mu * I
dIa/dt = lambda_1 * p_a * E1 - mu * Ia
dR/dt = mu * (I + Ia)
dC/dt = lambda_2 * E2
dD/dt = cfr * mu * I

text

### Parameters

| Symbol | Meaning | Default |
|--------|---------|---------|
| `beta`  | Transmission rate (may be a function of `t`) | 0.30 /day |
| `r_ia`  | Relative infectiousness of `Ia` | 0.1 |
| `r_e2`  | Relative infectiousness of `E2` | 0.0 |
| `lambda_1` | Rate `E1 → E2` or `E1 → Ia` | 0.25 /day |
| `lambda_2` | Rate `E2 → I` | 0.25 /day |
| `p_a`   | Fraction of infections that are asymptomatic | 0.05 |
| `mu`    | Rate `I → R` (symptomatic) | 0.125 /day |
| `cfr`   | Case fatality ratio among symptomatic cases | 0.47 |

Default durations implied: mean incubation `1/lambda_1 + 1/lambda_2 = 8` days,
mean symptomatic infectious period `1/mu = 8` days.

### Basic reproduction number

The next-generation `R0`, weighting the two branches by `p_a`, is
R0 = beta * [ (1 - p_a) * (r_e2 / lambda_2 + 1 / mu)

p_a * (r_ia / mu) ]

text

A second method `R0_book()` reproduces the *unweighted* formula printed in
the book, kept only for side-by-side comparison.

### `beta(t)`

Two time-dependent transmission profiles are provided:

- `smooth_step_beta(beta_a, beta_b, t_change, width)` — a logistic drop
  (or rise) from `beta_a` to `beta_b`, used in the fit.
- `piecewise_beta(t_starts, values)` — exact piecewise-constant `beta(t)`,
  matching the book's suggestion for modelling interventions.

Fitting uses the smooth step so the least-squares Jacobian stays well-behaved;
`width → 0` recovers a hard step.

---

## The solver

`ODESolver.py` implements a small class hierarchy in the style of
Section 2.2 of the book:
ODESolver # base class: initial conditions, time loop, advance()
├── ForwardEuler # first-order, explicit
└── RungeKutta4 # fourth-order, classical RK4

text

- The base class handles initial-condition bookkeeping, array allocation
  (scalar vs. system), and the time-stepping loop.
- Each subclass implements a single `advance()` method — the only thing
  that differs between methods.
- The right-hand side is always called as `f(t, u)` and its return value is
  coerced to a NumPy array, so models may return lists, tuples, or arrays.

`RungeKutta4` is used throughout the fitting and plotting in this project.
`ForwardEuler` is included so the two can be compared directly (e.g. on
the convergence tests in the book's Chapter 2).

---

## Data

### Sources

- **Primary**: a public GitHub tracker of DRC INSP situation reports
  (`webbegole/ebola-bundibugyo-tracker-2026`), downloaded as plain CSV.
- **Cross-check**: the ECDC Ebola outbreak page, parsed from prose or
  fetched with Selenium if the plain request is blocked.

### Files in the repo

| File | What it is |
|------|-----------|
| `fetch_data.py` | Downloader for both sources; writes the raw CSV |
| `ebola_drc_uganda.csv` | Raw download, exactly as fetched |
| `clean_data.py` | Cleans the raw file into a dense daily series |
| `ebola_drc_uganda_clean.csv` | Cleaned output, ready for the fit |
| `ebola_drc_uganda_template.csv` | Earlier hand-curated template, kept for reference |

### Schema (raw and cleaned)

| Column | Meaning |
|--------|---------|
| `date` | Observation date (ISO) |
| `country` | `DRC` or `Uganda` |
| `cum_confirmed_cases` | Cumulative confirmed cases |
| `cum_confirmed_deaths` | Cumulative confirmed deaths |
| `cum_recovered` | Cumulative recovered |
| `quality` | `reported`, `compiled`, or `milestone` *(raw only)* |
| `source_note` | Provenance *(raw only)* |

### Cleaning rules

`clean_data.py` applies four rules, in order:

1. **Normalize column names** — strip whitespace and lowercase, so a BOM
   or trailing space in the header doesn't break parsing.
2. **Insert missing calendar days** — `reindex` to a complete daily range
   within the observed window, per country.
3. **Forward-fill cumulative columns** — a missing day means no new
   report, so the cumulative value is carried forward. Interpolation would
   invent growth that was never reported and can produce non-integer
   cumulative counts.
4. **Enforce monotonicity** — `cummax()` on cumulative columns, to guard
   against any source-to-source dip that would break the log-residual.

Deaths and recoveries are allowed to stay `NaN` where the surveillance
system genuinely had no value yet (as for DRC deaths before 2026-05-24).
Only `cum_confirmed_cases` is required for a row to be kept.

---

## Installation

Requires Python 3.9+.

### With conda (recommended)

```bash
conda create -n seeiir python=3.10 -y
conda activate seeiir
conda install -y numpy scipy pandas matplotlib
Optional, only needed if you use the fetch/cross-check features:

bash
conda install -y requests beautifulsoup4
pip install selenium
With pip
bash
python -m venv .venv
# Windows:
.venv\Scripts\activate
# macOS / Linux:
source .venv/bin/activate

pip install numpy scipy pandas matplotlib requests beautifulsoup4
No ODESolver package is installed from PyPI — ODESolver.py lives in
this repository and is imported directly. If a stale odesolver package is
present in your environment from earlier experiments, uninstall it with
pip uninstall odesolver so it does not shadow the local file on
case-insensitive filesystems.

Usage
End-to-end, three steps
bash
python fetch_data.py --crosscheck
python clean_data.py ebola_drc_uganda.csv ebola_drc_uganda_clean.csv
python seeiir_ebola.py
Each script has sensible defaults that pick up the previous step's output.
If you already have a cleaned CSV (e.g. the one committed in this repo),
you can skip straight to seeiir_ebola.py.

fetch_data.py
bash
python fetch_data.py                       # tracker -> ebola_drc_uganda.csv
python fetch_data.py --crosscheck          # also compare latest row with ECDC
python fetch_data.py --crosscheck --selenium
python fetch_data.py --source ecdc --append  # append an ECDC snapshot
Flag	Default	Meaning
--out	ebola_drc_uganda.csv	Output CSV
--source	tracker	tracker or ecdc
--crosscheck	off	Also fetch ECDC and compare
--selenium	off	Use headless Chrome for the ECDC page
--append	off	Append to / update an existing --out
--enforce-monotone	off	Apply cummax() during validation
--raw-dir	raw	Where to archive the raw download
clean_data.py
bash
python clean_data.py ebola_drc_uganda.csv ebola_drc_uganda_clean.csv
python clean_data.py ebola_drc_uganda.csv ebola_drc_uganda_clean.csv DRC
Args, in order: input CSV, output CSV, optional country filter.

seeiir_ebola.py
bash
python seeiir_ebola.py \
    --data ebola_drc_uganda_clean.csv \
    --country DRC \
    --max-beta 0.5 \
    --out seeiir_new_fit.png
Flag	Default	Meaning
--data	ebola_drc_uganda_clean.csv	Cleaned CSV
--country	DRC	DRC or Uganda
--population	1.0e7	Assumed effective population at risk
--horizon	60.0	Days to project past the last observation
--reporting	1.0	Fraction of symptomatic cases confirmed
--min-date	(none)	Drop observations before this date
--deaths-start	(none)	Drop death observations before this date
--max-beta	0.5	Upper bound on beta
--out	ebola_seeiir_fit.png	Output figure path
Output
The script prints something like:

text
Loaded 145 DRC case observations, 136 with a matching death observation.
  t=0 is 2026-05-01, t_max = 158 days (2026-10-06)

Fitted parameters
  beta_a      0.xxx
  beta_b      0.xxx
  t_change    xx.xxx
  log10_e0     x.xxx
  cfr          0.xxx
  R0 early (beta_a) = x.xx   R0 late (beta_b) = x.xx
  (book's unweighted formula would give x.xx / x.xx)
  RMS log-residual = x.xxx

Saved seeiir_new_fit.png
Population conservation error: x.xxe-xx
and writes a four-panel figure:

Panel	Contents
Top-left	Cumulative confirmed cases, model vs. data (log scale)
Top-right	Cumulative deaths, model vs. data
Bottom-left	Active compartments E1, E2, I, Ia
Bottom-right	Fitted beta(t)
The dotted vertical line in the top panels marks the last observation used
in the fit; everything to the right is projection.

Interpreting the current figure
The figure shows the fit through the full DRC series (2026-05-15 to
2026-10-06, t = 0 on 2026-05-01). Several features are worth noting:

The cases panel shows the model tracking the data closely through
the initial exponential phase and the subsequent slowdown. The data
points crowd the model curve because the cleaned series is daily; the
visible "thickening" of the black dots between days 20 and 150 is the
daily cadence, not scatter.

The deaths panel shows the model slightly under-predicting deaths
in the middle of the fit window (around days 90–150) and then projecting
forward well above the last observed value. This is expected: the CFR
parameter is identified mainly by the ratio of deaths to cases, and once
beta has dropped, the model's cumulative-deaths curve keeps rising
because the existing I compartment is still large.

The beta(t) panel shows the fitted transmission dropping from
≈ 0.50 /day to ≈ 0.15 /day around day 30–50. That change point is
earlier than the visual slowdown in the case curve, which is typical
when the fit is dominated by the early explosive phase.

The active-compartment panel shows I (symptomatic infectious)
dominating the later dynamics, with E1 and E2 following and Ia
staying small because p_a = 0.05 is low.

If the deaths projection looks too aggressive, the two knobs to turn are
--cfr (fix it and refit without it) and --max-beta (lower it to force
an earlier drop).

Project layout
text
SEEIIR_solver/ODE_SEEIIR_model/
├── ODESolver.py                      # ODE solver class hierarchy (RK4, Euler)
├── seeiir_ebola.py                   # Model, fitting, plotting, CLI
├── fetch_data.py                     # Download raw data from public sources
├── clean_data.py                     # Clean raw data into a dense daily series
├── ebola_drc_uganda.csv              # Raw download
├── ebola_drc_uganda_clean.csv        # Cleaned input to the fit
├── ebola_drc_uganda_template.csv     # Earlier hand-curated template
├── seeiir_new_fit.png                # Output figure from the latest run
├── .gitignore
├── LICENSE
└── README.md
Fitting procedure
Fitting is done by scipy.optimize.least_squares on the residual

text
r_i = log1p(model_cases(t_i)) - log1p(data_cases(t_i))     # for each case point
   ⊕ log1p(model_deaths(t_i)) - log1p(data_deaths(t_i))    # for each death point
The model is solved with RungeKutta4 at a fixed internal step of
dt = 0.25 days and evaluated at observation times by linear interpolation.
Five parameters are fitted:

Parameter	Bounds	Notes
beta_a	[0.01, max_beta]	Early transmission
beta_b	[0.01, max_beta]	Late transmission
t_change	[5, t_max]	Midpoint of the beta drop
log10_e0	[-1, 3]	Initial exposed seed, in log₁₀
cfr	[0.2, 0.9]	Case fatality ratio
x_scale is set per parameter so the optimiser treats beta, t_change,
and e0 on comparable footing.

Design notes
No SciPy ODE solver. solve_ivp is deliberately not used. The point
is to have a fully explicit RungeKutta4 whose behaviour can be reasoned
about from first principles, matching the book's approach.

Model as callable. SEEIIR implements __call__(t, u) so the solver
treats it as a plain right-hand side, exactly as in the book's Chapter 5.

beta is a callable. Transmission may be constant, logistic, or
piecewise-constant, without changing the model class.

Auxiliary states. C and D are part of the state vector and are
integrated by RK4 along with everything else — no post-hoc integration
or trapezoidal quadrature.

Conservation check. The script reports the maximum absolute deviation
of S + E1 + E2 + I + Ia + R from N, as recommended in the book.

Raw data is preserved. fetch_data.py writes the raw download to
disk, and clean_data.py reads it without modifying it. The intermediate
file is committed, so the cleaning step is reproducible and auditable.

Limitations
The fit is to cumulative counts, so residuals are highly correlated
across time. This is fine for parameter estimation but invalidates
standard confidence intervals — none are reported.

beta(t) is constrained to a single sigmoid drop. Real outbreaks have
multiple interventions; extending to a multi-step beta is a one-line
change to smooth_step_beta usage.

r_ia, r_e2, p_a, and the two lambda rates are fixed, not fitted.
They should be varied in a sensitivity analysis before drawing conclusions.

The default cfr of 0.47 is a working assumption for Bundibugyo, not a
fitted quantity. Treat the fitted cfr with caution: it is identified
mainly by the ratio of deaths to cases, and the death series is short.

--population is an assumption. For the current case counts it has almost
no effect on the fit; it matters only once the epidemic is a sizeable
fraction of N.

Cumulative-deaths projection is uncertain. Because the model has no
death-specific dynamics beyond cfr * mu * I, the projected deaths curve
can overshoot once beta has dropped and the case curve has plateaued.
This is a structural limitation of the model, not a fitting artefact.

Roadmap
Ideas for extending this project:

□ Add a Heun / ExplicitMidpoint subclass to ODESolver.py for
convergence-order tests (as in the book's Chapter 2).
□ Add an adaptive-step RungeKuttaFehlberg45 following Chapter 4.
□ Add a backward-Euler / SDIRK solver for the stiff regime (Chapter 3).
□ Fit r_ia, r_e2, and p_a jointly with beta, with priors.
□ Multi-step beta(t) to represent successive interventions.
□ Bootstrap or profile-likelihood confidence intervals on R0.
□ Add an SEIR-with-vaccination variant.
□ Add a --fit-cfr-only mode that fixes beta_a/beta_b and fits only
the CFR, for sensitivity analysis.
References
The solver hierarchy, the model structure, and the fitting/validation
workflow are all adapted from:

Joakim Sundnes, Solving Ordinary Differential Equations in Python,
Springer, 2023.
Companion site: https://sundnes.github.io/solving_odes_in_python/

Additional background on the SEEIIR structure and its use during the
COVID-19 pandemic:

Kvaerno, A. Singly diagonally implicit Runge–Kutta methods with an
explicit first stage. BIT Numerical Mathematics, 44:489–502, 2004.

License
Add a license of your choice (e.g. MIT) before publishing. If you use the
code or data in academic work, please cite both this repository and the
Sundnes book above.

text
