"""
03_oxygen_in_the_streambed.py
=============================
Where does a streambed run out of oxygen, and can the two temperature loggers of script 02 tell us?

Stream water carries dissolved oxygen into the bed under downwelling; microbes consume it as it travels.
Under upwelling, oxygen-free groundwater pushes against the stream and oxygen can only enter by diffusion.
The depth at which oxygen falls below 2 mg/L (hypoxia) marks where denitrification and other anaerobic
processes can start, the oxic-anoxic switch that Zarnetske et al. (2012) framed with a Damkohler number.

MODFLOW 6 (GWF + GWT) is verified against the steady analytical profile C = C0 exp(lam z), then used to map
the hypoxic depth against the flux, and fed with the fluxes recovered from the loggers in script 02.

Setup
-----
Oxygen      9 mg/L in the stream (about saturation at 20 deg C); 0 in upwelling groundwater
Consumption first-order, k = 0.5 and 1.4 /d (flume dunes with leaf litter, Reeder et al., 2018:
            0.6e-5 to 1.6e-5 /s) and 5 /d (field value at 20 deg C used by Tonina et al., 2015)
Diffusion   molecular diffusion of oxygen in water, 2.0e-9 m2/s (1.73e-4 m2/d), divided by the tortuosity
            factor 1 - ln(n^2) of Boudreau (1996): 5.1e-5 m2/d in the pore water; dispersivity left out
Grid        geometric, 10 micrometres at the surface, 2 m deep; 30 days to steady state
Hypoxia     2 mg/L

Output: results/03_oxygen.json, figures/03_oxygen_in_the_streambed.png
"""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import plotstyle
from column import Sediment, oxygen_profile, run_oxygen

ROOT = Path(__file__).resolve().parent
SED = Sediment()
C0, C_HYP = 9.0, 2.0
D = 1.73e-4 / (1 - np.log(SED.porosity ** 2))
RATES = [0.5, 1.4, 5.0]
FLUXES = [-0.4, -0.2, -0.1, 0.0, 0.1, 0.2, 0.4, 0.6]
BED = 2.0


def hypoxic_depth(z, C):
    """Depth (m) where the profile first falls below C_HYP; inf if it never does within the bed."""
    below = np.flatnonzero(C < C_HYP)
    if below.size == 0:
        return np.inf
    i = below[0]
    return float(np.interp(C_HYP, [C[i], C[i - 1]], [z[i], z[i - 1]]))


def hypoxic_depth_analytical(q, k):
    v = q / SED.porosity
    lam = (v - np.sqrt(v * v + 4 * D * k)) / (2 * D)
    return np.log(C_HYP / C0) / lam


if __name__ == "__main__":
    plotstyle.apply()
    out = {"C0_mg_L": C0, "hypoxia_mg_L": C_HYP, "D_m2_per_d": D, "runs": [], "from_loggers": []}
    profiles = {}
    for k in RATES:
        for q in FLUXES:
            z, C = run_oxygen(SED, q, k=k, D=D, name=f"o2_{k}_{q:+.2f}")
            Ca = oxygen_profile(z, C0, q / SED.porosity, D, k)
            zh, za = hypoxic_depth(z, C), hypoxic_depth_analytical(q, k)
            inside = z <= min(za, BED) * 1.0
            err = float(np.max(np.abs(C[inside] - Ca[inside]))) if inside.any() else 0.0
            out["runs"].append({"k_per_d": k, "q_m_per_d": q, "hypoxic_depth_mf6_m": None if np.isinf(zh) else zh,
                                "hypoxic_depth_analytical_m": za, "max_abs_error_mg_L": err})
            profiles[(k, q)] = (z, C)
            zs = "beyond 2 m" if np.isinf(zh) else f"{1000 * zh:8.2f} mm"
            print(f"k {k:3.1f} /d, q {q:+.2f} m/d: hypoxic depth MODFLOW {zs}, analytical {1000 * za:9.2f} mm, "
                  f"largest profile difference {err:.3f} mg/L")

    # From two temperature loggers to the hypoxic depth (fluxes recovered in script 02)
    flux = json.load(open(ROOT / "results" / "02_flux.json"))
    for c in flux["cases"]:
        row = {"q_true": c["q_true"]}
        for k in RATES:
            est = [hypoxic_depth_analytical(c[key], k) for key in
                   ("q_amplitude_median", "q_amplitude_lambda_minus20", "q_amplitude_lambda_plus20")]
            row[f"k{k}"] = {"true_m": hypoxic_depth_analytical(c["q_true"], k), "from_loggers_m": est[0],
                            "range_conductivity_pm20_m": [min(est[1:]), max(est[1:])]}
        out["from_loggers"].append(row)
    json.dump(out, open(ROOT / "results" / "03_oxygen.json", "w"), indent=2)
    print("\nHypoxic depth for k = 1.4 /d, true flux against flux from the loggers (conductivity +/-20 %):")
    for r in out["from_loggers"]:
        x = r["k1.4"]
        print(f"  q {r['q_true']:+.2f} m/d: {100 * x['true_m']:7.2f} cm true, {100 * x['from_loggers_m']:7.2f} cm from loggers "
              f"[{100 * x['range_conductivity_pm20_m'][0]:.2f}, {100 * x['range_conductivity_pm20_m'][1]:.2f}]")

    # ------------------------------------------------------------------ figure
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.9))
    a = ax[0]
    cols = {-0.2: plotstyle.C2, 0.0: plotstyle.INK_2, 0.2: plotstyle.C1, 0.6: plotstyle.C3}
    for q, c in cols.items():
        z, C = profiles[(1.4, q)]
        zz = np.geomspace(1e-5, BED, 400)
        a.plot(oxygen_profile(zz, C0, q / SED.porosity, D, 1.4), zz, color=c, lw=1.6,
               label=f"q = {q:+.1f} m/d" + (" (upwelling)" if q < 0 else " (downwelling)" if q > 0 else " (diffusion only)"))
        sel = np.unique(np.geomspace(1, len(z) - 1, 28).astype(int))
        a.plot(C[sel], z[sel], "o", ms=4, color=c, mec=plotstyle.SURFACE, mew=1)
    a.axvline(C_HYP, color=plotstyle.INK_2, lw=1, ls=":")
    a.text(C_HYP + 0.15, 1.3e-5, "hypoxia, 2 mg/L", fontsize=8.5, color=plotstyle.INK_2)
    a.set_yscale("log")
    a.set_ylim(BED, 1e-5)
    a.set_xlim(0, 9.5)
    a.set_xlabel("Dissolved oxygen (mg/L)")
    a.set_ylabel("Depth below the streambed surface (m)")
    a.set_title("Oxygen profiles, k = 1.4 /d: MODFLOW 6 (dots), exact (lines)", fontsize=11)
    a.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2, fontsize=8.5)

    a = ax[1]
    qq = np.linspace(-0.45, 0.65, 221)
    for k, c in zip(RATES, (plotstyle.C1, plotstyle.C2, plotstyle.C3)):
        a.plot(qq, [hypoxic_depth_analytical(x, k) for x in qq], color=c, lw=1.6, label=f"k = {k} /d")
        pts = [r for r in out["runs"] if r["k_per_d"] == k and r["hypoxic_depth_mf6_m"] is not None]
        a.plot([r["q_m_per_d"] for r in pts], [r["hypoxic_depth_mf6_m"] for r in pts], "o", color=c,
               mec=plotstyle.SURFACE, mew=1.2)
    a.axhline(BED, color=plotstyle.GRID, lw=1.2)
    a.axvline(0, color=plotstyle.GRID, lw=1.2)
    a.set_yscale("log")
    a.set_ylim(5e-5, 6)
    a.set_xlabel("Darcy flux q (m/d), positive downward")
    a.set_ylabel("Depth where oxygen falls below 2 mg/L (m)")
    a.set_title("The flux decides where the streambed turns hypoxic", fontsize=11)
    a.text(-0.43, 2.4e-4, "upwelling:\nhypoxic within\na fraction of a mm", fontsize=8.5, color=plotstyle.INK_2)
    a.text(0.08, 0.012, "downwelling:\noxic over decimetres\nto metres", fontsize=8.5, color=plotstyle.INK_2)
    a.legend(loc="upper left", fontsize=8.5)
    fig.tight_layout()
    fig.savefig(ROOT / "figures" / "03_oxygen_in_the_streambed.png")
    print("Figure saved: figures/03_oxygen_in_the_streambed.png")
