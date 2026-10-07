"""ODESolver class hierarchy, following the structure used in the book
(Section 2.2): the superclass handles initial conditions and the time loop,
subclasses implement a single `advance` step.

The right-hand side is called as f(t, u), as in the book's ODESolver chapter.
"""
import numpy as np


class ODESolver:
    def __init__(self, f):
        # Wrap the user's f so list/tuple output always becomes an array
        self.model = f
        self.f = lambda t, u: np.asarray(f(t, u), float)

    def set_initial_condition(self, u0):
        if np.isscalar(u0):                 # scalar ODE
            self.neq = 1
            u0 = float(u0)
        else:                               # system of ODEs
            u0 = np.asarray(u0, float)
            self.neq = u0.size
        self.u0 = u0

    def solve(self, t_span, N):
        """Compute solution for t_span[0] <= t <= t_span[1] using N steps."""
        t0, T = t_span
        self.dt = (T - t0) / N
        self.t = np.zeros(N + 1)
        self.u = np.zeros(N + 1) if self.neq == 1 else np.zeros((N + 1, self.neq))

        assert hasattr(self, "u0"), "Please set initial condition before calling solve"

        self.t[0] = t0
        self.u[0] = self.u0
        for n in range(N):
            self.n = n
            self.t[n + 1] = self.t[n] + self.dt
            self.u[n + 1] = self.advance()
        return self.t, self.u

    def advance(self):
        raise NotImplementedError("advance is not implemented in the base class")


class ForwardEuler(ODESolver):
    def advance(self):
        u, f, n, t = self.u, self.f, self.n, self.t
        return u[n] + self.dt * f(t[n], u[n])


class RungeKutta4(ODESolver):
    def advance(self):
        u, f, n, t = self.u, self.f, self.n, self.t
        dt = self.dt
        dt2 = dt / 2.0
        k1 = f(t[n], u[n])
        k2 = f(t[n] + dt2, u[n] + dt2 * k1)
        k3 = f(t[n] + dt2, u[n] + dt2 * k2)
        k4 = f(t[n] + dt, u[n] + dt * k3)
        return u[n] + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
