"""
column.py
=========
A vertical streambed column in MODFLOW 6, built with FloPy, plus the analytical solutions used to check it.

Units: metres, days, degrees Celsius, joules, kilograms, mg/L.

The column runs from the streambed surface (cell 1, whose temperature or concentration is imposed by the
stream) down to a depth L. A constant head difference between the top and bottom cells drives a uniform
vertical Darcy flux q: q > 0 is downward (losing stream, downwelling), q < 0 is upward (gaining stream,
upwelling). Depths are measured from the centre of the top cell, where the stream boundary acts.
"""
import os
import shutil
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import flopy

ROOT = Path(__file__).resolve().parent
WORK = ROOT / "_mf6_runs"

# --- MODFLOW 6 executable: environment variable, then PATH, then a local ./bin folder
MF6 = os.environ.get("MF6_EXE") or shutil.which("mf6") or str(ROOT / "bin" / "mf6")

SECONDS_PER_DAY = 86400.0


@dataclass
class Sediment:
    """Saturated sandy streambed.

    Typical values for quartz sand and water, within the ranges compiled by Stonestrom and Constantz
    (2003, USGS Circular 1260): bulk thermal conductivity about 1.8 W/m/K, volumetric heat capacity of the
    bulk sediment about 2.8 MJ/m3/K.
    """
    porosity: float = 0.30
    K: float = 10.0                 # hydraulic conductivity (m/d), medium sand
    rho_w: float = 1000.0           # kg/m3
    c_w: float = 4180.0             # J/kg/K
    rho_s: float = 2650.0           # quartz, kg/m3
    c_s: float = 840.0              # quartz, J/kg/K
    lam_w: float = 0.58             # thermal conductivity of water, W/m/K
    lam_bulk: float = 1.8           # bulk thermal conductivity of the saturated sediment, W/m/K

    @property
    def lam_s(self):
        """Solid conductivity such that MODFLOW's arithmetic mean n*k_w + (1-n)*k_s equals lam_bulk."""
        return (self.lam_bulk - self.porosity * self.lam_w) / (1.0 - self.porosity)

    @property
    def rhoc_w(self):
        return self.rho_w * self.c_w

    @property
    def rhoc_b(self):
        n = self.porosity
        return n * self.rho_w * self.c_w + (1 - n) * self.rho_s * self.c_s

    @property
    def kappa(self):
        """Effective thermal diffusivity (m2/d)."""
        return self.lam_bulk * SECONDS_PER_DAY / self.rhoc_b

    def thermal_velocity(self, q):
        """Velocity of the thermal front (m/d) for a Darcy flux q (m/d)."""
        return q * self.rhoc_w / self.rhoc_b


# ----------------------------------------------------------------------------- analytical solutions
def heat_lambda(v, kappa, period=1.0):
    """Complex decay rate of a periodic temperature signal in steady 1D vertical flow.

    Stallman (1965) solution of  dT/dt = kappa d2T/dz2 - v dT/dz  (z downward) for a sinusoidal surface
    temperature: T - Tmean = A exp(Re(lam) z) cos(w t + Im(lam) z), with lam the root of
    kappa lam^2 - v lam - i w = 0 that decays with depth. Hatch et al. (2006) and Keery et al. (2007) use
    the same relation to infer the flux from the amplitude ratio and the phase lag between two sensors.
    """
    w = 2 * np.pi / period
    return (v - np.sqrt(v * v + 4j * kappa * w)) / (2 * kappa)


def amplitude_ratio(v, kappa, dz, period=1.0):
    return np.exp(np.real(heat_lambda(v, kappa, period)) * dz)


def phase_lag_days(v, kappa, dz, period=1.0):
    w = 2 * np.pi / period
    return -np.imag(heat_lambda(v, kappa, period)) * dz / w


def oxygen_profile(z, C0, v, D, k):
    """Steady dissolved-oxygen profile with first-order consumption k (1/d), pore velocity v (m/d, downward
    positive) and diffusion-dispersion D (m2/d): C = C0 exp(lam z), lam = (v - sqrt(v^2 + 4 D k)) / (2 D)."""
    lam = (v - np.sqrt(v * v + 4 * D * k)) / (2 * D)
    return C0 * np.exp(lam * z)


# ----------------------------------------------------------------------------- MODFLOW 6 models
def _grid(L, dz, growth=1.0):
    """Layers from the streambed surface to depth L: uniform (growth = 1) or geometric, starting at dz.

    Returns the number of layers, their bottoms, and the depth of each cell centre below the centre of
    the top cell (where the stream boundary acts).
    """
    if growth == 1.0:
        thick = np.full(int(round(L / dz)), dz)
    else:
        thick = [dz]
        while sum(thick) < L:
            thick.append(thick[-1] * growth)
        thick = np.array(thick) * L / sum(thick)
    botm = -np.cumsum(thick)
    centre = botm + thick / 2
    return len(thick), botm, centre[0] - centre


def _flow(sim, sed, q, n, botm, depth):
    """Steady GWF model with a uniform vertical Darcy flux q (m/d)."""
    gwf = flopy.mf6.ModflowGwf(sim, modelname="gwf", save_flows=True)
    ims = flopy.mf6.ModflowIms(sim, complexity="simple", inner_dvclose=1e-12, outer_dvclose=1e-11,
                               filename="gwf.ims")
    sim.register_ims_package(ims, [gwf.name])
    flopy.mf6.ModflowGwfdis(gwf, nlay=n, nrow=1, ncol=1, delr=1.0, delc=1.0, top=0.0, botm=botm)
    flopy.mf6.ModflowGwfnpf(gwf, k=sed.K, icelltype=0, save_specific_discharge=True)
    # Head difference between the centres of the top and bottom cells that gives the flux q
    dist = depth[-1]
    h_top, h_bot = 10.0, 10.0 - q * dist / sed.K
    flopy.mf6.ModflowGwfic(gwf, strt=h_top)
    flopy.mf6.ModflowGwfchd(gwf, auxiliary=["TEMPERATURE", "CONCENTRATION"],
                            stress_period_data=[[(0, 0, 0), h_top, 0.0, 0.0], [(n - 1, 0, 0), h_bot, 0.0, 0.0]],
                            pname="CHD-1")
    flopy.mf6.ModflowGwfoc(gwf, budget_filerecord="gwf.cbc", saverecord=[("BUDGET", "LAST")])
    return gwf


def run_heat(sed, q, *, L=2.0, dz=0.01, days=8.0, dt=1 / 96, t_mean=12.0, t_amp=4.0, t_gw=None,
             period=1.0, scheme="TVD", name="heat"):
    """Streambed temperatures under a daily stream temperature cycle (MODFLOW 6 GWF + GWE).

    Returns times (d), cell-centre depths below the top cell centre (m), and temperatures [time, depth].
    """
    t_gw = t_mean if t_gw is None else t_gw
    ws = WORK / name
    if ws.exists():
        shutil.rmtree(ws)
    n, botm, depth = _grid(L, dz)
    nstp = int(round(days / dt))
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=str(ws), exe_name=MF6, verbosity_level=0)
    flopy.mf6.ModflowTdis(sim, time_units="days", nper=1, perioddata=[(days, nstp, 1.0)])
    _flow(sim, sed, q, n, botm, depth)

    gwe = flopy.mf6.ModflowGwe(sim, modelname="gwe")
    ims = flopy.mf6.ModflowIms(sim, complexity="simple", linear_acceleration="BICGSTAB",
                               inner_dvclose=1e-10, outer_dvclose=1e-9, filename="gwe.ims")
    sim.register_ims_package(ims, [gwe.name])
    flopy.mf6.ModflowGwedis(gwe, nlay=n, nrow=1, ncol=1, delr=1.0, delc=1.0, top=0.0, botm=botm)
    flopy.mf6.ModflowGweic(gwe, strt=t_gw)
    flopy.mf6.ModflowGweadv(gwe, scheme=scheme)
    k_w = sed.lam_w * SECONDS_PER_DAY            # W/m/K to J/d/m/K
    k_s = sed.lam_s * SECONDS_PER_DAY
    flopy.mf6.ModflowGwecnd(gwe, alh=0.0, alv=0.0, ath1=0.0, ath2=0.0, atv=0.0, ktw=k_w, kts=k_s)
    flopy.mf6.ModflowGweest(gwe, porosity=sed.porosity, heat_capacity_water=sed.c_w, density_water=sed.rho_w,
                            heat_capacity_solid=sed.c_s, density_solid=sed.rho_s)
    # Stream temperature: daily sinusoid, given as a time series at every time step (15 min at most)
    ts_t = np.arange(0.0, days + 1e-9, min(dt, 1 / 96))
    ts = [(t, t_mean + t_amp * np.sin(2 * np.pi * t / period)) for t in ts_t]
    flopy.mf6.ModflowGwectp(gwe, stress_period_data=[[(0, 0, 0), "stream"], [(n - 1, 0, 0), t_gw]],
                            timeseries={"filename": "stream.ts", "time_series_namerecord": "stream",
                                        "timeseries": ts, "interpolation_methodrecord": "linear"})
    flopy.mf6.ModflowGwessm(gwe, sources=[("CHD-1", "AUX", "TEMPERATURE")])
    flopy.mf6.ModflowGweoc(gwe, temperature_filerecord="gwe.ucn", saverecord=[("TEMPERATURE", "ALL")])
    flopy.mf6.ModflowGwfgwe(sim, exgtype="GWF6-GWE6", exgmnamea="gwf", exgmnameb="gwe")

    sim.write_simulation(silent=True)
    ok, buff = sim.run_simulation(silent=True)
    if not ok:
        raise RuntimeError("MODFLOW 6 failed:\n" + "\n".join(buff[-20:]))
    hds = flopy.utils.HeadFile(str(ws / "gwe.ucn"), text="TEMPERATURE", precision="double")
    times = np.array(hds.get_times())
    T = np.array([hds.get_data(totim=t)[:, 0, 0] for t in times])
    return times, depth, T


def run_oxygen(sed, q, *, k, D, C0=9.0, L=2.0, dz=1e-5, growth=1.05, days=30.0, dt=0.02, name="oxygen"):
    """Dissolved oxygen in the streambed (MODFLOW 6 GWF + GWT), first-order consumption k (1/d).

    The stream supplies C0 at the top cell; groundwater arriving from below carries no oxygen.
    The grid is geometric, 10 micrometres at the surface, because oxygen can vanish within a fraction of a
    millimetre under upwelling. Returns depths below the top cell centre (m) and the final (steady) profile.
    """
    ws = WORK / name
    if ws.exists():
        shutil.rmtree(ws)
    n, botm, depth = _grid(L, dz, growth)
    nstp = int(round(days / dt))
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=str(ws), exe_name=MF6, verbosity_level=0)
    flopy.mf6.ModflowTdis(sim, time_units="days", nper=1, perioddata=[(days, nstp, 1.0)])
    _flow(sim, sed, q, n, botm, depth)

    gwt = flopy.mf6.ModflowGwt(sim, modelname="gwt")
    ims = flopy.mf6.ModflowIms(sim, complexity="simple", linear_acceleration="BICGSTAB",
                               inner_dvclose=1e-12, outer_dvclose=1e-11, filename="gwt.ims")
    sim.register_ims_package(ims, [gwt.name])
    flopy.mf6.ModflowGwtdis(gwt, nlay=n, nrow=1, ncol=1, delr=1.0, delc=1.0, top=0.0, botm=botm)
    flopy.mf6.ModflowGwtic(gwt, strt=0.0)
    flopy.mf6.ModflowGwtadv(gwt, scheme="TVD")
    flopy.mf6.ModflowGwtdsp(gwt, diffc=D, alh=0.0, alv=0.0, ath1=0.0, ath2=0.0, atv=0.0)
    flopy.mf6.ModflowGwtmst(gwt, porosity=sed.porosity, first_order_decay=True, decay=k)
    # Stream water at the top; under upwelling, oxygen-free groundwater enters at the bottom. Under
    # downwelling the bottom is a free outflow, so the column behaves as a semi-infinite bed.
    cells = [[(0, 0, 0), C0]] + ([[(n - 1, 0, 0), 0.0]] if q <= 0 else [])
    flopy.mf6.ModflowGwtcnc(gwt, stress_period_data=cells)
    flopy.mf6.ModflowGwtssm(gwt, sources=[("CHD-1", "AUX", "CONCENTRATION")])
    flopy.mf6.ModflowGwtoc(gwt, concentration_filerecord="gwt.ucn", saverecord=[("CONCENTRATION", "LAST")])
    flopy.mf6.ModflowGwfgwt(sim, exgtype="GWF6-GWT6", exgmnamea="gwf", exgmnameb="gwt")

    sim.write_simulation(silent=True)
    ok, buff = sim.run_simulation(silent=True)
    if not ok:
        raise RuntimeError("MODFLOW 6 failed:\n" + "\n".join(buff[-20:]))
    ucn = flopy.utils.HeadFile(str(ws / "gwt.ucn"), text="CONCENTRATION", precision="double")
    C = ucn.get_data()[:, 0, 0]
    return depth, C


def fit_sinusoid(t, y, period=1.0):
    """Least-squares fit y = m + s t + a cos(w t) + b sin(w t); returns amplitude and phase (rad).

    The linear term absorbs a slow drift, as in field records.
    """
    w = 2 * np.pi / period
    X = np.column_stack([np.ones_like(t), t, np.cos(w * t), np.sin(w * t)])
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    a, b = coef[2], coef[3]
    return np.hypot(a, b), np.arctan2(a, b)   # y ~ A sin(w t + phi)
