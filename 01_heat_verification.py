"""
01_heat_verification.py
=======================
Verifies the MODFLOW 6 groundwater energy transport model (GWE) for a streambed under a daily stream
temperature cycle, against the analytical solution of Stallman (1965).

Setup
-----
Column      2 m of saturated sand, 1 cm cells, uniform vertical Darcy flux q
Flux        q from -0.5 m/d (gaining stream, upwelling) to +0.5 m/d (losing stream, downwelling)
Stream      12 +/- 4 deg C, one cycle per day; groundwater at 12 deg C at the bottom
Sediment    porosity 0.30, bulk thermal conductivity 1.8 W/m/K, quartz grains (Stonestrom and Constantz, 2003)
Time        6 days, 3.75 min steps; the last 3 days are analysed
Sensors     0.10, 0.20 and 0.30 m below the streambed surface

Checks
------
- amplitude ratio and phase lag of the daily wave at each sensor, against Stallman (1965)
- time-step and mesh convergence: MODFLOW 6 solves the transient term with implicit Euler, which damps
  the daily wave; the error should fall in proportion to the time step

Output: results/01_heat.json, figures/01_heat_verification.png
"""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import plotstyle
from column import Sediment, amplitude_ratio, fit_sinusoid, phase_lag_days, run_heat

ROOT = Path(__file__).resolve().parent
SED = Sediment()
FLUXES = [-0.5, -0.25, 0.0, 0.25, 0.5]          # m/d, positive downward
SENSORS = [0.10, 0.20, 0.30]                   # m below the streambed surface
DZ, DT, DAYS, SPIN = 0.01, 1 / 384, 6.0, 3.0


def wave(t, depth, T, d):
    """Amplitude ratio and phase lag (h) of the daily wave at depth d, relative to the surface cell."""
    sel = t >= SPIN
    a0, p0 = fit_sinusoid(t[sel], T[sel, 0])
    i = int(np.argmin(np.abs(depth - d)))
    a, p = fit_sinusoid(t[sel], T[sel, i])
    return a / a0, ((p0 - p) % (2 * np.pi)) / (2 * np.pi) * 24, depth[i]


if __name__ == "__main__":
    plotstyle.apply()
    out = {"sediment": {"porosity": SED.porosity, "lambda_bulk_W_mK": SED.lam_bulk,
                        "kappa_m2_per_d": SED.kappa, "rhoc_bulk_J_m3K": SED.rhoc_b},
           "runs": [], "convergence": []}
    series = {}
    for q in FLUXES:
        t, depth, T = run_heat(SED, q, dz=DZ, dt=DT, days=DAYS, name=f"heat_{q:+.2f}")
        v = SED.thermal_velocity(q)
        for d in SENSORS:
            ar, lag, z = wave(t, depth, T, d)
            ar_a, lag_a = amplitude_ratio(v, SED.kappa, z), phase_lag_days(v, SED.kappa, z) * 24
            out["runs"].append({"q_m_per_d": q, "depth_m": d, "Ar_mf6": ar, "Ar_analytical": ar_a,
                                "lag_h_mf6": lag, "lag_h_analytical": lag_a,
                                "Ar_rel_error": ar / ar_a - 1, "lag_error_min": (lag - lag_a) * 60})
        series[q] = (t, depth, T)
        worst = max(abs(r["Ar_rel_error"]) for r in out["runs"] if r["q_m_per_d"] == q)
        print(f"q = {q:+.2f} m/d: largest amplitude-ratio error {100 * worst:.2f} %")

    # Convergence at q = +0.25 m/d, sensor at 0.30 m
    q = 0.25
    v = SED.thermal_velocity(q)
    for dz in (0.01, 0.005):
        for dt in (1 / 96, 1 / 192, 1 / 384, 1 / 768):
            t, depth, T = run_heat(SED, q, dz=dz, dt=dt, days=5.0, name="heat_conv")
            ar, lag, z = wave(t, depth, T, 0.30)
            e_ar = ar / amplitude_ratio(v, SED.kappa, z) - 1
            e_lag = (lag - phase_lag_days(v, SED.kappa, z) * 24) * 60
            out["convergence"].append({"dz_m": dz, "dt_min": dt * 1440, "Ar_rel_error": e_ar, "lag_error_min": e_lag})
            print(f"dz {100 * dz:.1f} cm, dt {dt * 1440:5.2f} min: amplitude ratio {100 * e_ar:+.2f} %, lag {e_lag:+.1f} min")
    json.dump(out, open(ROOT / "results" / "01_heat.json", "w"), indent=2)

    # ------------------------------------------------------------------ figure
    fig, ax = plt.subplots(2, 2, figsize=(11, 8.2))
    cols = [plotstyle.C1, plotstyle.C2, plotstyle.C3]
    for a, q, title in ((ax[0, 0], 0.25, "Losing stream (downwelling, q = +0.25 m/d)"),
                        (ax[0, 1], -0.25, "Gaining stream (upwelling, q = -0.25 m/d)")):
        t, depth, T = series[q]
        sel = t >= SPIN
        a.plot(t[sel], T[sel, 0], color=plotstyle.INK_2, lw=1.2, label="streambed surface")
        for c, d in zip(cols, SENSORS):
            i = int(np.argmin(np.abs(depth - d)))
            a.plot(t[sel], T[sel, i], color=c, lw=1.8, label=f"{100 * d:.0f} cm")
        a.set_title(title, fontsize=11)
        a.set_xlabel("Day")
        a.set_ylabel("Temperature (°C)")
        a.set_ylim(7.5, 16.5)
    ax[0, 0].legend(loc="lower left", ncol=4, fontsize=8.5)

    qq = np.linspace(-0.6, 0.6, 121)
    a = ax[1, 0]
    for c, d in zip(cols, SENSORS):
        a.plot(qq, amplitude_ratio(SED.thermal_velocity(qq), SED.kappa, d), color=c, lw=1.6, label=f"Stallman, {100 * d:.0f} cm")
        pts = [r for r in out["runs"] if r["depth_m"] == d]
        a.plot([r["q_m_per_d"] for r in pts], [r["Ar_mf6"] for r in pts], "o", color=c, mec=plotstyle.SURFACE, mew=1.5)
    a.set_xlabel("Darcy flux q (m/d), positive downward")
    a.set_ylabel("Amplitude ratio (depth / surface)")
    a.set_title("Daily wave: MODFLOW 6 (dots) against Stallman (lines)", fontsize=11)
    a.legend(loc="upper left", fontsize=8.5)

    a = ax[1, 1]
    for dz, mk in ((0.01, "o"), (0.005, "s")):
        pts = [r for r in out["convergence"] if r["dz_m"] == dz]
        a.loglog([r["dt_min"] for r in pts], [abs(r["Ar_rel_error"]) * 100 for r in pts], mk + "-",
                 color=plotstyle.C1 if dz == 0.01 else plotstyle.C2, label=f"cells of {100 * dz:.1f} cm")
    xs = np.array([1.875, 15])
    a.loglog(xs, 0.45 * xs / 1.875, ":", color=plotstyle.INK_2, lw=1.2, label="slope 1")
    a.set_xticks([2, 4, 8, 15], ["2", "4", "8", "15"])
    a.set_yticks([0.5, 1, 2, 4], ["0.5", "1", "2", "4"])
    a.minorticks_off()
    a.set_xlabel("Time step (min)")
    a.set_ylabel("Amplitude-ratio error at 30 cm (%)")
    a.set_title("The time step sets the error, not the mesh", fontsize=11)
    a.legend(loc="upper left", fontsize=8.5)
    fig.tight_layout()
    fig.savefig(ROOT / "figures" / "01_heat_verification.png")
    print("Figure saved: figures/01_heat_verification.png")
