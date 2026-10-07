# SEEIIR Ebola Outbreak Model — DRC / Uganda 2026

A deterministic SEIR-type compartmental model of the 2026 Bundibugyo Ebola
outbreak, fitted to cumulative confirmed cases and deaths from the DRC and
Uganda. The solver is a from-scratch `RungeKutta4` implementation built on a
small ODE-solver class hierarchy — no SciPy `solve\_ivp`, no external ODE
libraries — so the numerics are transparent and inspectable.

The model code, the solver, and the fitting pipeline together make up a
self-contained reproduction of the modelling workflow described in
Joakim Sundnes' *Solving Ordinary Differential Equations in Python*
(see **References** below).

!\[Fit](ebola\_seeiir\_fit.png)

\---

## What this is

A small, self-contained Python project that:

1. Implements a **SEEIIR** compartmental model for Ebola transmission,
with two asymptomatic / pre-symptomatic infectious compartments so the
model can represent the known epidemiology of Bundibugyo virus.
2. Solves the model with a **hand-written RK4 solver** implemented as a
subclass of a general `ODESolver` base class (forward Euler and RK4 are
both included for comparison).
3. **Fits** the model to surveillance data (cumulative confirmed cases and
cumulative deaths) by least squares on a log scale.
4. Produces a **four-panel diagnostic figure**: cases, deaths, active
compartments, and the fitted `beta(t)`.

The point of the project is not to produce a policy-grade forecast — it is
to demonstrate, on a real outbreak, that a small, clearly written ODE solver
plus a structured model class is enough to do meaningful epidemic modelling.

\---

## The model

### Compartments

|Symbol|Meaning|
|-|-|
|`S`|Susceptible|
|`E1`|Exposed, **not** yet infectious|
|`E2`|Exposed, infectious (pre-symptomatic, scaled by `r\_e2`)|
|`I`|Symptomatic infectious|
|`Ia`|Asymptomatic infectious (scaled by `r\_ia`)|
|`R`|Removed (recovered or dead)|
|`C`|Cumulative symptomatic cases *(auxiliary, non-feedback)*|
|`D`|Cumulative deaths *(auxiliary, non-feedback)*|

The two auxiliary states `C` and `D` are appended so the model's output can
be compared directly with surveillance data without post-hoc integration.
They do **not** feed back into the transmission dynamics.

### Equations

With `N = S + E1 + E2 + I + Ia + R` and force of infection

