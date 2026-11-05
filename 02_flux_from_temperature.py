"""
02_flux_from_temperature.py
===========================
Two temperature loggers buried in a streambed: how well do they give the water flux between the stream
and the aquifer?

MODFLOW 6 (GWF + GWE, verified in script 01) plays the field: it produces the "observed" temperatures at
two sensors, which are then sampled every 15 minutes with noise, as a logger would record them. The flux is
recovered with the two classic heat-tracer methods, the amplitude ratio and the phase lag of the daily wave
between the sensors (Hatch et al., 2006; Keery et al., 2007), and compared with the true flux.

Setup
-----
Stream      16 +/- 4 deg C per day; groundwater at 11 deg C at 2 m (different means, as in the field)
Sensors     0.10 and 0.30 m below the streambed surface, 15-minute logging, Gaussian noise of 0.05 deg C
Flux        seven cases from -0.6 m/d (upwelling) to +0.6 m/d (downwelling)
Analysis    days 3 to 7; 200 noise realisations per case
Uncertainty bulk thermal conductivity misjudged by -20 % and +20 % in the inversion

Output: results/02_flux.json, figures/02_flux_from_temperature.png
"""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

import plotstyle
from column import Sediment, amplitude_ratio, fit_sinusoid, phase_lag_days, run_heat

ROOT = Path(__file__).resolve().parent
SED = Sediment()
FLUXES = [-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6]
Z1, Z2 = 0.10, 0.30
NOISE, N_REAL, SPIN, DAYS = 0.05, 200, 3.0, 7.0
RNG = np.random.default_rng(7)


def kappa_with(lam_bulk):
    return lam_bulk * 86400 / SED.rhoc_b


def q_from_amplitude(ar, kappa, dz):
    """Darcy flux (m/d) from the amplitude ratio between two sensors dz apart (Hatch et al., 2006)."""
    f = lambda v: amplitude_ratio(v, kappa, dz) - ar
    v = brentq(f, -50, 50)
    return v * SED.rhoc_b / SED.rhoc_w


def q_from_phase(lag_days, kappa, dz):
    """Magnitude of the Darcy flux (m/d) from the phase lag; the phase does not tell the direction."""
    f = lambda v: phase_lag_days(v, kappa, dz) - lag_days
    if f(0.0) <= 0:          # lag at least as long as for pure conduction: no detectable flow
        return 0.0
    v = brentq(f, 0.0, 50)
    return v * SED.rhoc_b / SED.rhoc_w


def analyse(t, T1, T2):
    a1, p1 = fit_sinusoid(t, T1)
    a2, p2 = fit_sinusoid(t, T2)
    return a2 / a1, ((p1 - p2) % (2 * np.pi)) / (2 * np.pi)


if __name__ == "__main__":
    plotstyle.apply()
    dz = Z2 - Z1
    k_true = SED.kappa
    k_lo, k_hi = kappa_with(0.8 * SED.lam_bulk), kappa_with(1.2 * SED.lam_bulk)
    out = {"sensors_m": [Z1, Z2], "noise_degC": NOISE, "cases": []}
    example = None
    for q in FLUXES:
        t, depth, T = run_heat(SED, q, dt=1 / 384, days=DAYS, t_mean=16.0, t_amp=4.0, t_gw=11.0,
                               name=f"flux_{q:+.2f}")
        i1, i2 = int(np.argmin(np.abs(depth - Z1))), int(np.argmin(np.abs(depth - Z2)))
        log = (t >= SPIN) & (np.isclose((t * 96) % 1, 0, atol=1e-6) | np.isclose((t * 96) % 1, 1, atol=1e-6))
        tl, T1, T2 = t[log], T[log, i1], T[log, i2]
        ar0, lag0 = analyse(tl, T1, T2)
        q_amp, q_ph = [], []
        for _ in range(N_REAL):
            ar, lag = analyse(tl, T1 + RNG.normal(0, NOISE, tl.size), T2 + RNG.normal(0, NOISE, tl.size))
            q_amp.append(q_from_amplitude(ar, k_true, dz))
            q_ph.append(q_from_phase(lag, k_true, dz))
        q_amp, q_ph = np.array(q_amp), np.array(q_ph)
        case = {"q_true": q,
                "q_amplitude_median": float(np.median(q_amp)),
                "q_amplitude_p05": float(np.percentile(q_amp, 5)), "q_amplitude_p95": float(np.percentile(q_amp, 95)),
                "q_amplitude_lambda_minus20": q_from_amplitude(ar0, k_lo, dz),
                "q_amplitude_lambda_plus20": q_from_amplitude(ar0, k_hi, dz),
                "abs_q_phase_median": float(np.median(q_ph)),
                "abs_q_phase_p05": float(np.percentile(q_ph, 5)), "abs_q_phase_p95": float(np.percentile(q_ph, 95)),
                "amplitude_ratio": ar0, "phase_lag_h": lag0 * 24}
        out["cases"].append(case)
        print(f"q true {q:+.2f} m/d | amplitude method {case['q_amplitude_median']:+.3f} "
              f"[{case['q_amplitude_p05']:+.3f}, {case['q_amplitude_p95']:+.3f}], conductivity -20/+20 %: "
              f"{case['q_amplitude_lambda_minus20']:+.3f} / {case['q_amplitude_lambda_plus20']:+.3f} | "
              f"phase method |q| {case['abs_q_phase_median']:.3f} [{case['abs_q_phase_p05']:.3f}, {case['abs_q_phase_p95']:.3f}]")
        if q == 0.2:
            example = (tl, T[log, 0], T1, T2)
    json.dump(out, open(ROOT / "results" / "02_flux.json", "w"), indent=2)

    # ------------------------------------------------------------------ figure
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.6), gridspec_kw={"width_ratios": [1, 1.1]})
    tl, T0, T1, T2 = example
    a = ax[0]
    a.plot(tl, T0, color=plotstyle.INK_2, lw=1.2, label="stream (streambed surface)")
    a.plot(tl, T1 + RNG.normal(0, NOISE, tl.size), ".", ms=2.5, color=plotstyle.C1, label=f"logger at {100 * Z1:.0f} cm")
    a.plot(tl, T2 + RNG.normal(0, NOISE, tl.size), ".", ms=2.5, color=plotstyle.C2, label=f"logger at {100 * Z2:.0f} cm")
    a.set_xlabel("Day")
    a.set_ylabel("Temperature (°C)")
    a.set_title("What the loggers record (losing reach, q = +0.2 m/d)", fontsize=11)
    a.legend(loc="lower left", fontsize=8.5, ncol=1)

    a = ax[1]
    cs = out["cases"]
    qt = np.array([c["q_true"] for c in cs])
    a.plot([-0.7, 0.7], [-0.7, 0.7], color=plotstyle.GRID, lw=1.5, zorder=0)
    a.fill_between(qt, [c["q_amplitude_lambda_minus20"] for c in cs], [c["q_amplitude_lambda_plus20"] for c in cs],
                   color=plotstyle.C1, alpha=0.15, lw=0, label="amplitude, conductivity ±20 %")
    med = np.array([c["q_amplitude_median"] for c in cs])
    lo, hi = np.array([c["q_amplitude_p05"] for c in cs]), np.array([c["q_amplitude_p95"] for c in cs])
    a.errorbar(qt, med, yerr=[med - lo, hi - med], fmt="o", color=plotstyle.C1, mec=plotstyle.SURFACE, mew=1.5,
               capsize=3, label="amplitude ratio (sign and size)")
    ph = np.array([c["abs_q_phase_median"] for c in cs])
    a.plot(qt, ph, "^", color=plotstyle.C2, mec=plotstyle.SURFACE, mew=1.2, ms=7, label="phase lag (size only, |q|)")
    a.set_xlim(-0.7, 0.7)
    a.set_ylim(-0.75, 0.75)
    a.set_xlabel("True Darcy flux (m/d), positive downward")
    a.set_ylabel("Flux recovered from the two loggers (m/d)")
    a.set_title("Two loggers recover the flux", fontsize=11)
    a.legend(loc="lower right", fontsize=8.5)
    fig.tight_layout()
    fig.savefig(ROOT / "figures" / "02_flux_from_temperature.png")
    print("Figure saved: figures/02_flux_from_temperature.png")
