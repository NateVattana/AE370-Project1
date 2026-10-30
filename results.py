# Results for Section IV (A-D).
# Uses the model and RK4 solver in gimli_simulation.py.
# Run:  python results.py
# Saves figures to results/ and prints the numbers.

import math
import os
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import gimli_simulation as gs

DT = 1.0


def dynamics_for(passengers):
    mass = gs.m_fixed + passengers * gs.m_p
    return lambda t, y: gs.aircraft_dynamics(t, y, mass=mass)


def trim_start(h0, passengers):
    """Steady glide at h0 for the fixed C_L (speed and angle)."""
    mass = gs.m_fixed + passengers * gs.m_p
    rho = gs.rho_SL * math.exp(-h0 / gs.H)
    LD = gs.C_L / (gs.C_D0 + gs.K * gs.C_L**2)
    gamma = -math.atan(1 / LD)
    V = math.sqrt(2 * mass * gs.g * math.cos(gamma)
                  / (rho * gs.S * gs.C_L))
    return V, gamma


def run(h0=gs.h0, passengers=gs.N_p, trimmed=False):
    if trimmed:
        V0, g0 = trim_start(h0, passengers)
    else:
        V0, g0 = gs.V0, gs.gamma0
    times, states = gs.simulate_glide(
        dt=DT, initial_altitude=h0, initial_speed=V0,
        passengers=passengers, initial_gamma=g0)
    return times, states, passengers


def arrival(times, states, passengers, x_target=gs.x_G):
    """State when x = x_target, found by bisection on the
    last RK4 step (same approach as the h_G crossing).
    Returns None if the glide ends short of x_target."""
    f = dynamics_for(passengers)
    for k in range(len(states) - 1):
        if states[k + 1][0] >= x_target:
            lo, hi = 0.0, times[k + 1] - times[k]
            for _ in range(50):
                mid = (lo + hi) / 2
                y = gs.rk4_step(times[k], states[k], mid, f=f)
                if y[0] < x_target:
                    lo = mid
                else:
                    hi = mid
            y = gs.rk4_step(times[k], states[k], hi, f=f)
            return times[k] + hi, y
    return None


def summary(label, h0=gs.h0, passengers=gs.N_p, trimmed=False):
    t, s, p = run(h0, passengers, trimmed)
    a = arrival(t, s, p)
    if a is None:
        print(f"{label}: falls short by "
              f"{(gs.x_G - s[-1][0]) / 1000:.2f} km")
    else:
        ta, ya = a
        print(f"{label}: {ya[1] - gs.h_G:.0f} m above field, "
              f"V = {ya[2]:.1f} m/s, t = {ta:.0f} s")
    return t, s, p


def margin(h0, trimmed):
    _, s, _ = run(h0=h0, trimmed=trimmed)
    return s[-1][0] - gs.x_G


def critical_h0(trimmed, lo=1000.0, hi=gs.h0):
    for _ in range(40):
        mid = (lo + hi) / 2
        if margin(mid, trimmed) >= 0:
            hi = mid
        else:
            lo = mid
    return hi

def main():
    os.makedirs("results", exist_ok=True)


    # ---------- A. Baseline ----------
    print("--- A. Baseline ---")
    t, s, p = summary("Gimli")
    print(f"Range: {s[-1][0] / 1000:.1f} km, "
        f"margin: {(s[-1][0] - gs.x_G) / 1000:.1f} km")


    # ---------- C. Initial altitude ----------
    print("\n--- C. Initial altitude ---")
    h0s = [1000 * k for k in range(2, 13)]
    curves = {}
    for trimmed in (False, True):
        curves[trimmed] = [(gs.x_G + margin(h, trimmed)) / 1000
                        for h in h0s]
    hc_fixed = critical_h0(False)
    hc_trim = critical_h0(True)
    for name, hc in (("fixed V0", hc_fixed), ("trimmed", hc_trim)):
        print(f"{name}: critical h0 = {hc:.0f} m "
            f"({hc / 0.3048:,.0f} ft)")

    plt.figure(figsize=(7, 4))
    plt.plot([h / 1000 for h in h0s], curves[False], "o-",
            label=r"Fixed start ($V_0$ = 246.8 m/s, $\gamma_0$ = 0)")
    plt.plot([h / 1000 for h in h0s], curves[True], "s-",
            label="Trimmed start (steady glide)")
    plt.axhline(gs.x_G / 1000, color="red", ls="--", label="Gimli")
    plt.scatter([hc_fixed / 1000, hc_trim / 1000],
                [gs.x_G / 1000] * 2, color="k", zorder=5)
    plt.axvline(gs.h0 / 1000, color="gray", ls=":",
                label="Flight 143 (35,000 ft)")
    plt.xlabel("Altitude at flameout (km)")
    plt.ylabel("Glide range (km)")
    plt.legend(fontsize=8)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("results/C_altitude.png", dpi=300)
    plt.close()


    # ---------- D. Combined boundary ----------
    print("\n--- D. Combined boundary ---")
    grid_pax = list(range(0, 251, 25))
    grid_h = [2500 + 500 * k for k in range(18)]
    M = []
    for h in grid_h:
        row = []
        for n in grid_pax:
            _, s, _ = run(h0=h, passengers=n)
            row.append((s[-1][0] - gs.x_G) / 1000)
        M.append(row)

    plt.figure(figsize=(7, 4.5))
    X, Y = grid_pax, [h / 1000 for h in grid_h]
    norm = TwoSlopeNorm(vcenter=0)
    cf = plt.contourf(X, Y, M, levels=20, cmap="RdBu", norm=norm)
    plt.colorbar(cf, label="Range minus distance to Gimli (km)")
    plt.contour(X, Y, M, levels=[0], colors="k", linewidths=2)
    plt.plot(61, gs.h0 / 1000, "*", color="gold", ms=15,
            mec="k", label="Flight 143")
    plt.xlabel("Passengers")
    plt.ylabel("Altitude at flameout (km)")
    plt.title("Reaches Gimli (blue) vs. falls short (red)")
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig("results/D_boundary.png", dpi=300)
    plt.close()
    print("Figures saved to results/")

if __name__ == "__main__":
    main()