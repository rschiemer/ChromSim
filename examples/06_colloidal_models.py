# ---
# jupyter:
#   jupytext:
#     formats: py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Tutorial 06 — Colloidal isotherms side by side: CPA vs Colloidal (Xu & Lenhoff)
#
# ChromSim has two isotherms built on colloid theory:
#
# | | `CPA` (Briskot et al. 2021) | `Colloidal` (Xu & Lenhoff 2009) |
# |---|---|---|
# | Solver | reference backend (Python) | CADET (native), also reference |
# | Henry coefficient | double-layer physics (Z, radius, ligands, I) | empirical ln K_e(salt) |
# | Saturation | hard-disc exclusion **and** lateral repulsion | lattice repulsion only |
# | Screening length | Debye length from ionic strength | κ(salt) is a fit function |
#
# This notebook uses **synthetic data only**. One mAb on a cation exchanger is described with
# CPA, and the Colloidal model is then *calibrated to reproduce CPA*. Whatever still differs
# afterwards comes from the structure of the two models, not from arbitrary parameter choices.
#
# 1. Inside CPA: from salt to Henry coefficient, step by step.
# 2. Inside Colloidal: calibrating ln K_e to CPA's Henry coefficient.
# 3. The nonlinear part: what makes each isotherm saturate, and a global calibration.
# 4. Columns at **low load** (1 g/L): both models agree.
# 5. Columns at **high load** (40 g/L, breakthrough, load sweep): where they part ways.
# 6. Handling: parameters, kinetics, solvers, speed.
#
# Theory pages: `docs/theory/isotherms/cpa.md`, `docs/theory/isotherms/colloidal.md` and
# `docs/theory/colloidal_models.md`.

# %%
import time

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares

from chromsim import api
from chromsim.isotherms.custom.colloidal_particle_adsorption import (
    AVOGADRO,
    ELEMENTARY_CHARGE,
    available_surface_function,
    debye_kappa,
    grahame_potential,
    henry_factor,
    lateral_energy,
    well_depth,
)
from chromsim.isotherms.python_isotherm import binding_parameters
from chromsim.isotherms.registry import get_isotherm

plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})


def log_salt_axis(ax):
    """Logarithmic salt axis with readable tick labels."""
    ax.set_xscale("log")
    ax.set_xticks([20, 50, 100, 200, 500])
    ax.xaxis.set_major_formatter(plt.ScalarFormatter())
    ax.xaxis.set_minor_formatter(plt.NullFormatter())
    ax.set_xlabel("salt / mM")


CPA_COLOR, COLLOIDAL_COLOR = "#1f77b4", "#d62728"

NAMES = ["salt", "mAb"]
MW_KG_PER_MOL = 150.0  # mAb, 150 kDa: c in mol/m³ × 150 = g/L
LOAD_SALT = 50.0  # mM

cpa = get_isotherm("CPA")
colloidal = get_isotherm("Colloidal")

# %% [markdown]
# ## The synthetic "true" protein: CPA parameters
#
# A mAb (radius 5 nm, net charge +40 at the operating pH) on a strong cation exchanger with
# 2 µmol/m² ligands (1.2 charges/nm²) and 2·10⁸ m² of surface per m³ of solid. A lateral charge
# Z_lat = 10 sets the repulsion between adsorbed proteins. These are plausible magnitudes,
# not a specific resin.

# %%
CPA_PARAMETERS = {
    "mAb": {"radius": 5e-9, "charge": 40.0, "Delta": 0.3, "lateral_charge": 10.0, "k_kin": 10.0}
}
CPA_GLOBAL = {"ligand_density": 2e-6, "ligand_charge": -1.0, "specific_surface": 2e8}
cpa_p = binding_parameters(NAMES, CPA_PARAMETERS, CPA_GLOBAL, None)

RADIUS = CPA_PARAMETERS["mAb"]["radius"]
SURFACE = CPA_GLOBAL["specific_surface"]
MONOLAYER = SURFACE / (np.pi * AVOGADRO * RADIUS**2)  # q at full coverage Θ = 1, mol/m³ solid
print(
    f"full monolayer: q = {MONOLAYER:.2f} mol/m³ solid = {MONOLAYER * MW_KG_PER_MOL:.0f} g/L solid"
)


def cpa_factor(salt_mM, q):
    """q*/c of CPA for bound concentration(s) q (mol/m³ solid) at a salt concentration (mM)."""
    q = np.atleast_1d(q).astype(float)
    states = np.vstack([np.zeros_like(q), q])
    return cpa.equilibrium_factor(np.full_like(q, salt_mM), states, cpa_p)[1]


JAMMED = 1e3  # mol/m³ (150 kg/L): liquid concentrations above this mean "surface full"


def cpa_liquid(salt_mM, q, p=cpa_p):
    """CPA's liquid concentration in equilibrium with q, mol/m³ (NaN where the surface is full)."""
    q = np.atleast_1d(q).astype(float)
    factor = cpa.equilibrium_factor(np.full_like(q, salt_mM), np.vstack([np.zeros_like(q), q]), p)[
        1
    ]
    with np.errstate(divide="ignore", over="ignore"):
        liquid = q / factor
    return np.where(liquid < JAMMED, liquid, np.nan)


# %% [markdown]
# ## 1. Inside CPA: from salt to Henry coefficient
#
# At low loading CPA is linear, q = K_H·c. K_H follows from four steps, each a public function:
#
# 1. **Screening.** `debye_kappa(I)`: the Debye length 1/κ shrinks as 1/√I.
# 2. **Surface potentials.** `grahame_potential(σ, κ)`: the same surface charge gives a lower
#    potential when it is better screened.
# 3. **Energy well.** `well_depth(ψ_A, ψ_i, a)`: depth y of the sphere–plane double-layer
#    minimum, in k_BT. Attraction from 2ψ_Aψ_i, repulsion from the "image" term ψ_A² + ψ_i².
# 4. **Henry coefficient.** K_H = Δ·(e^y − 1)/y.

# %%
salt = np.geomspace(20, 600, 200)  # mM
kappa = debye_kappa(salt)
sigma_adsorber = -ELEMENTARY_CHARGE * AVOGADRO * CPA_GLOBAL["ligand_density"]  # C/m²
sigma_protein = 40.0 * ELEMENTARY_CHARGE / (4 * np.pi * RADIUS**2)  # C/m²
psi_adsorber = grahame_potential(sigma_adsorber, kappa)
psi_protein = grahame_potential(sigma_protein, kappa)
depth = well_depth(psi_adsorber, psi_protein, RADIUS)
henry = CPA_PARAMETERS["mAb"]["Delta"] * henry_factor(depth)
# Same numbers as the isotherm itself:
assert np.allclose(henry, [cpa_factor(s, 0.0)[0] for s in salt])

fig, axes = plt.subplots(1, 4, figsize=(15, 3.4))
axes[0].plot(salt, 1e9 / kappa)
axes[0].set(ylabel="Debye length / nm", title="1 · screening")
axes[1].plot(salt, 1e3 * psi_adsorber, label="adsorber ψ_A")
axes[1].plot(salt, 1e3 * psi_protein, label="protein ψ_i")
axes[1].set(ylabel="surface potential / mV", title="2 · potentials")
axes[1].legend()
axes[2].plot(salt, depth)
axes[2].set(ylabel="well depth y / k_BT", title="3 · energy well")
axes[3].semilogy(salt, henry, color=CPA_COLOR)
axes[3].set(ylabel="K_H (q per solid / c)", title="4 · Henry coefficient")
for ax in axes:
    log_salt_axis(ax)
fig.tight_layout()

# %% [markdown]
# Between 50 and 150 mM the well depth drops from ~14 to ~5 k_BT, and K_H falls by more than
# three decades. At high salt the well disappears and K_H levels off at Δ: the protein then
# only partitions into the adsorption layer.
#
# On a log–log plot SMA would give a straight line of slope −ν. CPA's line is curved: its
# local slope (an "apparent ν") changes with salt. This is one of its main differences from
# stoichiometric models.

# %%
apparent_nu = -np.gradient(np.log(henry), np.log(salt))
for s in (50, 100, 150, 250):
    print(f"apparent ν at {s:3d} mM: {np.interp(s, salt, apparent_nu):5.1f}")

# %% [markdown]
# ## 2. Inside Colloidal: an empirical ln K_e(salt)
#
# Colloidal's low-load limit is also linear, with Henry constant K_e:
#
# ln K_e = k_a·c_s^(−k_b) + k_c·exp(k_d·c_s)
#
# Here k_d = 0, so the exponential term is a constant offset k_c. We fit k_a, k_b and k_c so
# that K_e reproduces CPA's K_H between 30 and 500 mM.

# %%
fit_salt = np.geomspace(30, 500, 40)
ln_k_cpa = np.log([cpa_factor(s, 0.0)[0] for s in fit_salt])


def ln_ke(x, s):
    k_a, k_b, k_c = x
    return k_a * s ** (-k_b) + k_c


henry_fit = least_squares(lambda x: ln_ke(x, fit_salt) - ln_k_cpa, x0=[100.0, 1.0, 0.0])
k_a, k_b, k_c = henry_fit.x
print(f"k_a = {k_a:.1f}, k_b = {k_b:.3f}, k_c = {k_c:.3f}")
print(
    f"largest error in ln K: {np.abs(henry_fit.fun).max():.2f} "
    f"(factor {np.exp(np.abs(henry_fit.fun).max()):.2f})"
)

fig, axes = plt.subplots(1, 2, figsize=(10, 3.4))
axes[0].semilogy(salt, henry, color=CPA_COLOR, label="CPA K_H (physics)")
axes[0].semilogy(
    salt, np.exp(ln_ke(henry_fit.x, salt)), "--", color=COLLOIDAL_COLOR, label="Colloidal K_e (fit)"
)
axes[0].set(ylabel="Henry coefficient", title="Henry coefficients")
axes[0].legend()
axes[1].plot(fit_salt, henry_fit.fun, "o-", color=COLLOIDAL_COLOR)
axes[1].axhline(0, color="k", lw=0.8)
axes[1].set(ylabel="ln K_e − ln K_H", title="fit residual (fit range 30–500 mM)")
for ax in axes:
    log_salt_axis(ax)
fig.tight_layout()

# %% [markdown]
# The power law plus offset follows CPA within about 20 %. It can't follow CPA's curvature
# exactly, because CPA's K_H is not a power law. The two Henry coefficients play the same role:
# both are q (per solid volume) / c at vanishing load. **At low load the two models should
# therefore give nearly the same chromatograms** (section 4).
#
# ## 3. The nonlinear part: what makes each isotherm saturate
#
# With q > 0 both models add a "crowding" factor that lowers q*/c:
#
# * **CPA**: q*/c = K_H · B(Θ) · e^(−u_lat/k_BT).
#   * B(Θ) is hard-disc exclusion (scaled particle theory). It does not depend on salt and goes
#     to 0 as the surface jams.
#   * e^(−u_lat) is screened repulsion from neighbours, which weakens as salt rises.
# * **Colloidal**: q*/c = K_e · e^(−S(q)). S is the Yukawa repulsion on a hexagonal lattice,
#   with strength b_pp(salt) and range 1/κ(salt). Colloidal has **no separate excluded-area
#   term**, so S has to do all the work.

# %%
q = np.linspace(1e-4, 0.8 * MONOLAYER, 400)
surface_density = np.vstack([np.zeros_like(q), q]) / SURFACE  # Γ, mol/m²
radii = np.array([[0.0], [RADIUS]])
crowding_spt = available_surface_function(surface_density, radii)[1]

Z_LAT = np.array([[0.0], [CPA_PARAMETERS["mAb"]["lateral_charge"]]])
fig, axes = plt.subplots(1, 3, figsize=(15, 3.6), sharey=True)
for ax, s in zip(axes, (10, 50, 150), strict=True):
    lateral = lateral_energy(surface_density, radii, Z_LAT, np.full(q.size, debye_kappa(s)))[1]
    ax.semilogy(q / MONOLAYER, crowding_spt, label="B(Θ): hard discs")
    ax.semilogy(q / MONOLAYER, np.exp(-lateral), label="e^(−u_lat): lateral repulsion")
    ax.semilogy(q / MONOLAYER, crowding_spt * np.exp(-lateral), "k", lw=2, label="CPA total")
    ax.set(
        xlabel="coverage Θ = q / q_monolayer",
        title=f"CPA crowding factors at {s} mM",
        ylim=(1e-8, 2),
    )
axes[0].set_ylabel("factor multiplying K_H")
axes[0].legend(loc="lower left")
fig.tight_layout()

# %% [markdown]
# **Hard-disc exclusion does almost all the work.** Lateral repulsion only matters at low ionic
# strength (10 mM: up to ~3 k_BT at Θ ≈ 0.65 with Z_lat = 10). By 50 mM it is screened away,
# because at Θ = 0.5 the gap between neighbouring proteins (~3.5 nm) is already 2.5 Debye
# lengths. So at typical CEX load conditions, CPA's capacity is set by **geometry**: radius
# and specific surface. Salt enters only through K_H. The q*/c of CPA is then nearly
# K_H(salt) · g(Θ) with a salt-independent g.
#
# ### Calibrating Colloidal's crowding term to CPA
#
# The equilibrium liquid concentration for a given q is c = q/(q*/c) in CPA. In Colloidal it is
# c*(q), the quantity inside its driving force k_kin(c − c*). Setting c = 0 in the flux returns
# −k_kin·c*, which lets us read c*(q) straight from the model code.
#
# We fit Colloidal's lattice parameters (b_a, b_b, κ_f, κ_ef, κ_c) to CPA's isotherms at 50,
# 100 and 150 mM. The fit covers the liquid range that matters in practice (up to 15 g/L
# mAb). φ = A_s and the radius are taken from CPA.


# %%
def colloidal_parameters(lattice):
    """Colloidal parameters: Henry part from section 2, ``lattice`` = (b_a, b_b, κ_f, κ_ef, κ_c)."""
    b_a, b_b, kappa_f, kappa_ef, kappa_c = lattice
    parameters = {
        "mAb": {
            "radius": RADIUS,
            "ln_ke_salt_factor": float(k_a),
            "ln_ke_salt_exponent": float(k_b),
            "ln_ke_exp_factor": float(k_c),
            "ln_ke_exp_salt_coefficient": 0.0,
            "b_pp_salt_factor": float(b_a),
            "b_pp_salt_exponent": float(b_b),
        }
    }
    global_parameters = {
        "phase_ratio": SURFACE,
        "kappa_factor": float(kappa_f),
        "kappa_exponent": float(kappa_ef),
        "kappa_constant": float(kappa_c),
    }
    return parameters, global_parameters


def colloidal_c_star(salt_mM, q, lattice):
    """Colloidal's equilibrium liquid concentration c*(q), mol/m³.

    With c = 0 and k_kin = 1 the flux k_kin (c − c*) is −c*.

    NaN where the surface is full.
    """
    parameters, global_parameters = colloidal_parameters(lattice)
    parameters["mAb"]["k_kin"] = 1.0
    p = binding_parameters(NAMES, parameters, global_parameters, None)
    q = np.atleast_1d(q).astype(float)
    liquid = np.vstack([np.full_like(q, salt_mM), np.zeros_like(q)])
    with np.errstate(over="ignore"):
        c_star = -colloidal.flux(liquid, np.vstack([np.zeros_like(q), q]), p)[1]
    return np.where(c_star < JAMMED, c_star, np.nan)


FIT_SALTS = (50.0, 100.0, 150.0)
C_MAX = 15.0 / MW_KG_PER_MOL  # 15 g/L in mol/m³
fit_grids = {}
for s in FIT_SALTS:
    q_grid = np.linspace(0.01, 0.95 * MONOLAYER, 600)
    c_grid = cpa_liquid(s, q_grid)
    fit_grids[s] = q_grid[(c_grid > 1e-7) & (c_grid < C_MAX)]


def isotherm_residuals(x):
    b_a, b_b, log_kappa_f, kappa_ef, log_kappa_c = x
    lattice = (b_a, b_b, 10**log_kappa_f, kappa_ef, 10**log_kappa_c)
    return np.concatenate(
        [
            np.log(colloidal_c_star(s, q_s, lattice)) - np.log(cpa_liquid(s, q_s))
            for s, q_s in fit_grids.items()
        ]
    )


# κ_ef ≥ 0: CADET-Process 0.12 refuses negative exponents (ADR 0011).
lattice_fit = least_squares(
    isotherm_residuals,
    x0=[3.0, 0.0, -11.0, 1.0, -9.5],
    bounds=([0, -3, -16, 0, -11], [1e3, 3, -6, 3, -7]),
)
b_a, b_b, log_kappa_f, kappa_ef, log_kappa_c = lattice_fit.x
LATTICE = (b_a, b_b, 10**log_kappa_f, kappa_ef, 10**log_kappa_c)
print(
    f"b_pp = {b_a:.2f}·c_s^{b_b:.3f};  "
    f"1/κ = {LATTICE[2]:.2e}·c_s^{kappa_ef:.3f} + {LATTICE[4]:.1e} m"
)
print(f"rms error in ln c: {np.sqrt(np.mean(lattice_fit.fun**2)):.2f}")
for s in FIT_SALTS:
    fitted = (LATTICE[2] * s**kappa_ef + LATTICE[4]) * 1e9
    print(
        f"  {s:4.0f} mM: Colloidal screening length {fitted:4.1f} nm   "
        f"vs Debye length {1e9 / debye_kappa(s):4.2f} nm"
    )

# %% [markdown]
# **The calibrated screening length grows with salt, against Debye physics.** Colloidal has no
# excluded-area term, so to mimic the salt-independent hard-disc exclusion that dominates CPA
# at high salt, its repulsion must reach *further* as salt rises. CADET-Process only allows
# κ_ef ≥ 0, i.e. a screening length that grows with salt, which happens to be the direction
# this fit needs. The parameters still work, but Colloidal's κ should be read as a fit
# function, not as a Debye length.

# %%
fig, axes = plt.subplots(1, 3, figsize=(15, 3.6), sharey=True)
c_plot = np.geomspace(1e-5, C_MAX, 200)
for ax, s in zip(axes, FIT_SALTS, strict=True):
    q_dense = np.linspace(1e-4, 0.97 * MONOLAYER, 2000)
    c_cpa = cpa_liquid(s, q_dense)
    c_col = colloidal_c_star(s, q_dense, LATTICE)
    ax.plot(c_cpa * MW_KG_PER_MOL, q_dense * MW_KG_PER_MOL, color=CPA_COLOR, lw=2, label="CPA")
    ax.plot(
        c_col * MW_KG_PER_MOL,
        q_dense * MW_KG_PER_MOL,
        "--",
        color=COLLOIDAL_COLOR,
        lw=2,
        label="Colloidal (calibrated)",
    )
    ax.axhline(MONOLAYER * MW_KG_PER_MOL, color="grey", ls=":", label="monolayer")
    ax.set(xlim=(0, 15), xlabel="c / (g/L)", title=f"{s:.0f} mM salt")
axes[0].set_ylabel("q / (g/L solid)")
axes[0].legend(loc="lower right")
fig.tight_layout()

# %% [markdown]
# Both isotherms are rectangular at load conditions and get flatter with salt. They agree well
# at low loading. Close to the plateau CPA bends over more sharply, because exclusion diverges
# as the surface jams; Colloidal's lattice repulsion approaches it more gradually.
#
# **Sensitivity of each model's crowding parameter.** In CPA, Z_lat changes the isotherm only at
# low salt (20 mM, left), and the hard-disc limit stays in place even with Z_lat = 0. In
# Colloidal, b_pp is the *only* source of saturation: without it the isotherm is linear, with
# no capacity limit at all. That is why b_pp has a strong effect at any salt (right, 50 mM).

# %%
fig, axes = plt.subplots(1, 2, figsize=(11, 3.6), sharey=True)
q_dense = np.linspace(1e-4, 0.97 * MONOLAYER, 2000)
for z_lat in (0.0, 10.0, 20.0):
    p = binding_parameters(
        NAMES, {"mAb": CPA_PARAMETERS["mAb"] | {"lateral_charge": z_lat}}, CPA_GLOBAL, None
    )
    axes[0].plot(
        cpa_liquid(20.0, q_dense, p) * MW_KG_PER_MOL,
        q_dense * MW_KG_PER_MOL,
        label=f"Z_lat = {z_lat:g}",
    )
for scale in (0.25, 1.0, 4.0):
    lattice = (b_a * scale,) + LATTICE[1:]
    axes[1].plot(
        colloidal_c_star(50.0, q_dense, lattice) * MW_KG_PER_MOL,
        q_dense * MW_KG_PER_MOL,
        label=f"b_pp × {scale:g}",
    )
axes[0].set(
    xlim=(0, 15), xlabel="c / (g/L)", ylabel="q / (g/L solid)", title="CPA: lateral charge (20 mM)"
)
axes[1].set(xlim=(0, 15), xlabel="c / (g/L)", title="Colloidal: b_pp (50 mM)")
for ax in axes:
    ax.axhline(MONOLAYER * MW_KG_PER_MOL, color="grey", ls=":")
    ax.legend(loc="lower right")
fig.tight_layout()

# %% [markdown]
# ## 4. Column simulations
#
# Same column for both: LRMP, 10 cm × 0.5 cm, 150 cm/h. The method is load at 50 mM,
# 3 CV wash, then a 20 CV gradient from 50 to 500 mM.
#
# * **CPA** runs on the reference backend. It is kinetic: k_kin = 10 s⁻¹ is the desorption
#   rate, fast compared with mass transfer.
# * **Colloidal** runs on CADET in rapid equilibrium (no k_kin). Its k_kin would have a
#   different meaning: a rate on the driving force c − c*.

# %%
COLLOIDAL_PARAMETERS, COLLOIDAL_GLOBAL = colloidal_parameters(LATTICE)
ISOTHERMS = {
    "CPA": (
        {"model": "CPA", "parameters": CPA_PARAMETERS, "global_parameters": CPA_GLOBAL},
        "reference",
    ),
    "Colloidal": (
        {
            "model": "Colloidal",
            "parameters": COLLOIDAL_PARAMETERS,
            "global_parameters": COLLOIDAL_GLOBAL,
        },
        "cadet",
    ),
}
COLORS = {"CPA": CPA_COLOR, "Colloidal": COLLOIDAL_COLOR}


def make_project(
    model, load_g_per_L, feed_g_per_L=5.0, load_salt_mM=LOAD_SALT, gradient=True, n_col=60
):
    """LRMP column; load ``load_g_per_L`` mAb per litre of column, optionally a 20 CV gradient."""
    isotherm, backend = ISOTHERMS[model]
    load_cv = load_g_per_L / feed_g_per_L
    buffer = {"salt": f"{load_salt_mM} mM"}
    method = [
        {"name": "equilibrate", "duration": "1 CV", "flow": "150 cm/h", "inlet": buffer},
        {
            "name": "load",
            "duration": f"{load_cv} CV",
            "flow": "150 cm/h",
            "inlet": buffer | {"mAb": f"{feed_g_per_L} g/L"},
        },
    ]
    if gradient:
        method += [
            {"name": "wash", "duration": "3 CV", "flow": "150 cm/h", "inlet": buffer},
            {
                "name": "elute",
                "duration": "20 CV",
                "flow": "150 cm/h",
                "inlet": buffer,
                "inlet_end": {"salt": "500 mM"},
            },
        ]
    return api.Project.model_validate(
        {
            "schema_version": 1,
            "name": f"{model} {load_g_per_L} g/L",
            "components": [{"name": "salt"}, {"name": "mAb", "mw_kDa": 150}],
            "column": {
                "model": "LRMP",
                "length": "10 cm",
                "diameter": "0.5 cm",
                "axial_dispersion": "1e-7 m^2/s",
                "film_diffusion": {"salt": "1e-5 m/s", "mAb": "5e-6 m/s"},
            },
            "resin": {
                "name": "synthetic CEX",
                "particle_radius": "30 um",
                "ionic_capacity": "500 mM",
                "col_porosity": 0.35,
                "par_porosity": 0.6,
            },
            "isotherm": isotherm,
            "method": method,
            "solver": {"backend": backend, "n_col": n_col},
        }
    )


def run(model, *args, **kwargs):
    project = make_project(model, *args, **kwargs)
    start = time.perf_counter()
    result = api.simulate(project)
    return project, result, time.perf_counter() - start


def cv_from_gradient_start(project, result):
    """Volume in CV, counted from the start of the elution gradient."""
    gradient_start_s = project.resolved_method()[-1].start_s
    return result.column_volumes - np.interp(gradient_start_s, result.time_s, result.column_volumes)


def peak_summary(project, result):
    """Peak position (CV into gradient), salt at the peak (mM) and peak height (g/L)."""
    x = cv_from_gradient_start(project, result)
    protein = result.outlet("mAb")
    top = np.argmax(protein)
    return x[top], result.outlet("salt")[top], protein[top] * MW_KG_PER_MOL


# %% [markdown]
# ### Low load: 1 g/L column

# %%
low = {model: run(model, 1.0) for model in ISOTHERMS}

fig, ax = plt.subplots(figsize=(8, 3.8))
for model, (project, result, seconds) in low.items():
    x = cv_from_gradient_start(project, result)
    ax.plot(
        x,
        result.outlet("mAb") * MW_KG_PER_MOL,
        color=COLORS[model],
        label=f"{model} ({seconds:.1f} s)",
    )
    position, salt_at_peak, height = peak_summary(project, result)
    print(
        f"{model:<10} peak at {position:5.2f} CV into gradient, "
        f"{salt_at_peak:5.1f} mM salt, {height:.3f} g/L"
    )
salt_ax = ax.twinx()
salt_ax.plot(x, result.outlet("salt"), "k:", lw=1)
salt_ax.set_ylabel("salt / mM (dotted)")
ax.set(
    xlim=(0, 14),
    xlabel="CV from gradient start",
    ylabel="mAb / (g/L)",
    title="Low load (1 g/L): linear range",
)
ax.legend()
fig.tight_layout()

# %% [markdown]
# At low load both models are essentially in their Henry range, and the peaks almost coincide.
# The small offset is the ~20 % Henry-fit error from section 2: a gradient peak elutes where K
# drops to a value set by the gradient slope, so a 20 % error in K moves it by only a few mM.
#
# As a cross-check of the *solvers* rather than the models, Colloidal also runs on the
# reference backend (kinetic form, with a fast k_kin). This shows the difference between the
# two models above is not a CADET-vs-Python artefact.

# %%
project_ref = make_project("Colloidal", 1.0)
project_ref.solver.backend = "reference"
project_ref.isotherm.parameters["mAb"]["k_kin"] = 1e3  # fast driving-force kinetics ≈ equilibrium
result_ref = api.simulate(project_ref)
_, result_cadet, _ = low["Colloidal"]
print(
    "Colloidal, reference vs CADET peak:",
    peak_summary(project_ref, result_ref)[0].round(3),
    "vs",
    peak_summary(low["Colloidal"][0], result_cadet)[0].round(3),
    "CV",
)

# %% [markdown]
# ## 5. High load
#
# ### 5a. Gradient elution after loading 40 g/L

# %%
high = {model: run(model, 40.0) for model in ISOTHERMS}

fig, axes = plt.subplots(1, 2, figsize=(14, 3.8))
for model in ISOTHERMS:
    for ax, runs, title in ((axes[0], low, "1 g/L"), (axes[1], high, "40 g/L")):
        project, result, _ = runs[model]
        x = cv_from_gradient_start(project, result)
        ax.plot(x, result.outlet("mAb") * MW_KG_PER_MOL, color=COLORS[model], label=model)
        ax.set(
            xlim=(0, 14),
            xlabel="CV from gradient start",
            ylabel="mAb / (g/L)",
            title=f"Load {title}",
        )
    position, salt_at_peak, height = peak_summary(*high[model][:2])
    print(
        f"40 g/L {model:<10} peak at {position:5.2f} CV, {salt_at_peak:5.1f} mM, "
        f"{height:5.2f} g/L ({high[model][2]:.1f} s)"
    )
for ax in axes:
    ax.legend()
fig.tight_layout()

# %% [markdown]
# At high load both peaks move forward and sharpen at the front. This is the classic
# self-sharpening "triangular" elution of a convex isotherm: the crowded front sees a lower
# effective K and runs ahead. The two models are calibrated to the same equilibrium data,
# but the peaks still differ in height and shape. The difference comes from how each model's
# nonlinearity *changes along the gradient*: CPA's exclusion term doesn't depend on salt,
# while Colloidal's fitted repulsion grows with salt (section 3).
#
# ### 5b. Breakthrough: dynamic binding capacity at two load salts

# %%
fig, axes = plt.subplots(1, 2, figsize=(14, 3.8), sharey=True)
total_porosity = 0.35 + (1 - 0.35) * 0.6  # hold-up volume in CV
for ax, load_salt in zip(axes, (50.0, 100.0), strict=True):
    for model in ISOTHERMS:
        project, result, seconds = run(
            model, 150.0, load_salt_mM=load_salt, gradient=False, n_col=150
        )  # 30 CV
        x = result.column_volumes - 1.0  # CV loaded (equilibration was 1 CV)
        ratio = result.outlet("mAb") * MW_KG_PER_MOL / 5.0
        v10 = np.interp(0.1, ratio, x) if ratio.max() > 0.1 else np.nan
        dbc = 5.0 * (v10 - total_porosity)
        ax.plot(
            x, ratio, color=COLORS[model], label=f"{model}: DBC₁₀ ≈ {dbc:.0f} g/L ({seconds:.1f} s)"
        )
    ax.set(
        xlim=(0, 30),
        xlabel="CV loaded (5 g/L feed)",
        ylabel="c / c₀",
        title=f"Breakthrough at {load_salt:.0f} mM",
    )
    ax.legend(loc="upper left")
fig.tight_layout()

# %% [markdown]
# Breakthrough is the classic high-load experiment, yet here both models give **the same dynamic
# binding capacity**, and both drop by a factor of ~3 between 50 and 100 mM. That is expected:
# a breakthrough curve probes mainly one point of the isotherm (q at the feed concentration),
# and that point was part of the calibration. The fronts are steep in both models because
# the isotherms are strongly favourable, and the front width is set by mass transfer, which is
# identical. (DBC₁₀ here: loaded volume at 10 % breakthrough minus the column hold-up, times
# the feed concentration.)
#
# **Lesson: breakthrough data alone can't tell the two models apart.**
#
# ### 5c. Load sweep: how retention moves with load

# %%
loads = [1, 10, 20, 30, 40, 60]
sweep = {model: [peak_summary(*run(model, load)[:2]) for load in loads] for model in ISOTHERMS}

fig, axes = plt.subplots(1, 2, figsize=(12, 3.6))
for model, rows in sweep.items():
    rows = np.array(rows)
    axes[0].plot(loads, rows[:, 1], "o-", color=COLORS[model], label=model)
    axes[1].plot(loads, rows[:, 2], "o-", color=COLORS[model], label=model)
axes[0].set(
    xlabel="load / (g/L column)", ylabel="salt at peak maximum / mM", title="Elution salt vs load"
)
axes[1].set(xlabel="load / (g/L column)", ylabel="peak height / (g/L)", title="Peak height vs load")
for ax in axes:
    ax.legend()
fig.tight_layout()

# %% [markdown]
# The elution salt falls as load rises in both models, and up to ~30 g/L the two agree within
# a few mM. Above that they separate. At 60 g/L CPA elutes ~10 mM earlier, because its
# isotherm bends over more sharply near the plateau (section 3). **High-load gradient elution
# discriminates between the models better than breakthrough does.**
#
# ### 5d. Outside the calibration window: loading at 25 mM
#
# So far all conditions lay inside the 50–150 mM window used to calibrate Colloidal. Many
# processes load at lower conductivity, though. There CPA uses the same physics (Debye
# screening, lateral repulsion that now matters, see section 3). Colloidal *extrapolates*
# its fitted power laws for ln K_e, b_pp and κ.

# %%
fig, axes = plt.subplots(1, 2, figsize=(13, 3.8))
q_dense = np.linspace(1e-4, 0.97 * MONOLAYER, 2000)
for s, style in ((25.0, "-"), (50.0, ":")):
    axes[0].plot(
        cpa_liquid(s, q_dense) * MW_KG_PER_MOL,
        q_dense * MW_KG_PER_MOL,
        style,
        color=CPA_COLOR,
        label=f"CPA {s:.0f} mM",
    )
    axes[0].plot(
        colloidal_c_star(s, q_dense, LATTICE) * MW_KG_PER_MOL,
        q_dense * MW_KG_PER_MOL,
        style,
        color=COLLOIDAL_COLOR,
        label=f"Colloidal {s:.0f} mM",
    )
axes[0].set(
    xlim=(0, 15),
    xlabel="c / (g/L)",
    ylabel="q / (g/L solid)",
    title="Isotherms: inside (50 mM) and outside (25 mM) the fit",
)
axes[0].legend(loc="lower right")
for model in ISOTHERMS:
    project, result, seconds = run(model, 40.0, load_salt_mM=25.0)
    x = cv_from_gradient_start(project, result)
    axes[1].plot(x, result.outlet("mAb") * MW_KG_PER_MOL, color=COLORS[model], label=model)
    position, salt_at_peak, height = peak_summary(project, result)
    print(
        f"load at 25 mM, 40 g/L, {model:<10}: peak {position:5.2f} CV, "
        f"{salt_at_peak:5.1f} mM, {height:5.2f} g/L"
    )
axes[1].set(
    xlim=(0, 14),
    xlabel="CV from gradient start (25 → 500 mM)",
    ylabel="mAb / (g/L)",
    title="40 g/L loaded at 25 mM",
)
axes[1].legend()
fig.tight_layout()

# %% [markdown]
# Here the extrapolation goes surprisingly well. At 25 mM Colloidal's capacity is ~5 % below
# CPA's, and the elution peaks differ by ~3 mM, about as much as inside the window. Halving
# the salt is a mild extrapolation, and at 25 mM CPA's lateral term is still small. The
# structural point stands, though: Colloidal's ln K_e, b_pp and κ are power laws with no
# physics behind them outside the fitted range. CPA carries its physics into the new
# condition: Debye screening, and lateral repulsion that grows at low salt. In practice:
#
# * **calibrate Colloidal over the full salt range** you will simulate (load *and* elution);
# * use CPA when you need to **extrapolate** across conductivity, or to reason about how a
#   different resin (ligand density, surface) or molecule (charge, size) would behave.
#
# ## 6. Handling: what to know when you use them
#
# | | `CPA` | `Colloidal` |
# |---|---|---|
# | Backend | `reference` (LRM/LRMP) | `cadet` (any column) or `reference` |
# | Per protein | radius, charge, Delta, lateral_charge, k_kin | radius, 4 ln K_e + 4 b_pp, k_kin |
# | Resin-wide | ligand_density, ligand_charge, specific_surface | phase_ratio, κ_f, κ_ef, κ_c |
# | Physical meaning | yes: check against titration, DLS, resin data | fit functions |
# | k_kin | desorption rate (ChromSim) | rate on c − c* (CADET); omit = equilibrium |
# | Screening vs salt | Debye, from I | fitted κ(salt), κ_ef ≥ 0 (CADET-Process) |
# | Not (yet) included | pH dependence, charge regulation | pH (CADET's `use_ph`) |
#
# **Practical advice**
#
# * **Start low-load.** Both models reduce to a Henry coefficient. Calibrate it from LGE runs
#   (CPA: charge and Delta; Colloidal: k_a, k_b, k_c), then add the nonlinear parameters from
#   breakthrough or high-load runs. This is the same order as for SMA (Tutorial 03).
# * **Use log scales.** Delta, ligand_density, specific_surface and phase_ratio are set with
#   `log_scale`, so fitting (`FitParameter`) works on log10.
# * **Correlations.** In CPA, charge and Delta both shift retention. Charge sets the *slope* of
#   log K_H vs log salt, Delta the offset. In Colloidal, k_a and k_b correlate strongly,
#   and so do b_pp and κ. Check the correlation matrix after a fit.
# * **Speed.** Colloidal on CADET is the fast route for big studies (optimisation, GRM). CPA
#   on the reference backend runs seconds per simulation at these grid sizes, which is fine for
#   calibration and comparisons.
# * **Stiffness.** On the reference backend, Colloidal with a large k_kin and very high
#   loading (q beyond hexagonal close packing) becomes too stiff for the finite-difference
#   Jacobian. Use CADET for those cases.

# %%
timings = {model: (low[model][2], high[model][2]) for model in ISOTHERMS}
for model, (t_low, t_high) in timings.items():
    print(f"{model:<10} simulation time: {t_low:4.1f} s (1 g/L), {t_high:4.1f} s (40 g/L)")


# %% [markdown]
# ## Key takeaways
#
# 1. **Same low-load physics, different origins.** Both reduce to a Henry coefficient. CPA
#    derives it from charge, radius, ligand density and Debye screening, and its log K–log salt
#    curve is not a straight line ("apparent ν" from 11 to 2). Colloidal fits it with
#    ln K_e = k_a c^(−k_b) + k_c.
# 2. **Different sources of saturation.** In CPA at typical load salts, hard-disc exclusion
#    (geometry) sets the capacity; lateral charge matters only below ~20 mM. Colloidal has
#    only its lattice repulsion. To mimic CPA it needs a screening length that *grows* with
#    salt (4 → 35 nm), so its κ is a fit function, not a Debye length.
# 3. **Calibrated to the same isotherms, the models agree** at low load, in breakthrough
#    (DBC within ~2 %) and in gradient elution up to ~30 g/L. They separate at the highest
#    loads, where the shape of the isotherm near its plateau matters (≈ 10 mM in elution
#    salt at 60 g/L).
# 4. **Discriminating experiments**: high-load gradient elutions at several loads, not
#    breakthrough alone.
# 5. **Handling**: Colloidal is the fast, CADET-native choice for large studies within a
#    calibrated window. CPA has fewer and physically meaningful parameters and is the better
#    tool for extrapolation and mechanistic reasoning, but runs on the reference backend
#    (LRM/LRMP only, seconds per run).
