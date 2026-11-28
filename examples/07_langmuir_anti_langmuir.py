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
# # Tutorial 07 — Langmuirian vs anti-Langmuirian elution: SMA, SAS, Colloidal and CPA
#
# When more protein is loaded onto an ion exchanger, the elution peak moves.
#
# * **Langmuirian** behaviour: the peak elutes **earlier** (at lower salt) with increasing
#   load. It has a sharp front and a tailing rear ("shark fin" leaning backwards). The isotherm is
#   concave: q/c *falls* as c rises, because the binding sites fill up.
# * **Anti-Langmuirian** behaviour: the peak elutes **later** with increasing load. It
#   fronts and ends with a sharp rear. The isotherm is convex: q/c *rises* with c. Something
#   makes binding easier when more protein is already around, e.g. proteins binding to
#   already bound proteins (self-association).
#
# Real proteins often show **anti-Langmuir-to-Langmuir (AL-L)** behaviour: anti-Langmuirian at
# moderate load, Langmuirian once the column gets crowded.
#
# ## What the literature reports
#
# * **Bispecific mAb on POROS XS** (Seelinger et al. 2022, *J. Chromatogr. A* 1689):
#   anti-Langmuirian elution under high load, attributed to self-association. Modeled with a
#   pH-dependent SAS-SMA; plain SMA could only produce Langmuirian elution.
# * **Polypeptide on a cation exchanger** (Koch et al. 2022, *J. Sep. Sci.* 45, 2008):
#   fronting (anti-Langmuirian) at low load, Langmuirian at high load. Modeled with SMA plus
#   activity coefficients.
# * **Fc-fusion charge variants on CEX** (*J. Chromatogr. A*, 2025): AL-L transition with
#   load, modeled with the Mollerup activity model in a generalised SMA (nGIEX).
# * **mAb, industrial CEX** (Huuk et al. 2017, *Biotechnol. J.* 12, 1600336): unusual
#   high-load peak shapes, modeled with GIEX and asymmetric activity coefficients.
# * **mAb on Fractogel SO3⁻** (Diedrich et al. 2017, *J. Chromatogr. A* 1525, 60): shoulders
#   under overload, modeled with a multi-state SMA.
#
# The common thread is protein–protein interactions on the surface that SMA ignores. Here we
# compare four ChromSim isotherms on one synthetic mAb:
#
# * **SMA**: Langmuirian by construction.
# * **SAS** (self-association SMA, Mollerup 2008): SMA × (1 + k_sa·c). Can be anti-Langmuirian.
# * **Colloidal** (Xu & Lenhoff) and **CPA** (Briskot et al.): colloid-theory models from
#   Tutorial 06.
#
# The question for each: which way does the peak move with load, and why?

# %%
import time

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

from chromsim import api
from chromsim.isotherms.python_isotherm import binding_parameters
from chromsim.isotherms.registry import get_isotherm

plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})

NAMES = ["salt", "mAb"]
MW_KG_PER_MOL = 150.0  # c in mol/m³ × 150 = g/L
CAPACITY = 500.0  # ionic capacity Λ, mol/m³ solid (SMA, SAS)
LOAD_SALT, ELUTION_SALT = 50.0, 170.0  # mM: load condition; typical salt at the low-load peak

cpa = get_isotherm("CPA")
colloidal = get_isotherm("Colloidal")

# %% [markdown]
# ## 1. One protein, four isotherms
#
# The "true" protein is the CPA mAb from Tutorial 06 (radius 5 nm, charge +40). The other
# models are matched to it, so all four agree at low load:
#
# * **Colloidal**: parameters calibrated to CPA in Tutorial 06 (copied here).
# * **SMA**: ν and k_eq from CPA's Henry coefficient over the elution range (120–260 mM, a
#   straight line in log K vs log salt). σ is chosen so that SMA binds the same amount as CPA
#   at the load condition (5 g/L feed, 50 mM).
# * **SAS**: SMA plus self-association k_sa: weak (k_sa = 30 m³/mol) and strong (300 m³/mol).

# %%
CPA_PARAMETERS = {
    "mAb": {"radius": 5e-9, "charge": 40.0, "Delta": 0.3, "lateral_charge": 10.0, "k_kin": 10.0}
}
CPA_GLOBAL = {"ligand_density": 2e-6, "ligand_charge": -1.0, "specific_surface": 2e8}
cpa_p = binding_parameters(NAMES, CPA_PARAMETERS, CPA_GLOBAL, None)

# Calibrated to CPA in Tutorial 06 (section 3).
COLLOIDAL_PARAMETERS = {
    "mAb": {
        "radius": 5e-9,
        "ln_ke_salt_factor": 442.594,
        "ln_ke_salt_exponent": 0.92572,
        "ln_ke_exp_factor": -1.82151,
        "ln_ke_exp_salt_coefficient": 0.0,
        "b_pp_salt_factor": 153.310,
        "b_pp_salt_exponent": -0.99158,
    }
}
COLLOIDAL_GLOBAL = {
    "phase_ratio": 2e8,
    "kappa_factor": 1.92414e-13,
    "kappa_exponent": 2.40298,
    "kappa_constant": 1.87561e-9,
}


def cpa_henry(salt_mM):
    """CPA Henry coefficient K_H (q per solid / c) at salt concentration(s) in mM."""
    salt_mM = np.atleast_1d(salt_mM).astype(float)
    return cpa.equilibrium_factor(salt_mM, np.zeros((2, salt_mM.size)), cpa_p)[1]


def cpa_c(salt_mM, q):
    """CPA's equilibrium liquid concentration for bound q, mol/m³."""
    q = np.atleast_1d(q).astype(float)
    factor = cpa.equilibrium_factor(
        np.full_like(q, salt_mM), np.vstack([np.zeros_like(q), q]), cpa_p
    )[1]
    with np.errstate(divide="ignore", over="ignore"):
        c = q / factor
    return np.where(c < 1e3, c, np.nan)


# SMA Henry coefficient K_H = k_eq (Λ/c_s)^ν: fit a straight line in log–log space.
fit_salt = np.geomspace(120, 260, 30)
slope, intercept = np.polyfit(np.log(fit_salt), np.log(cpa_henry(fit_salt)), 1)
NU = -slope
K_EQ = np.exp(intercept) / CAPACITY**NU

# σ from the load condition: CPA binds q_load at 5 g/L and 50 mM. SMA's equilibrium
# q = K_H(50) c (1 − (ν+σ) q/Λ)^ν, solved for σ.
C_FEED = 5.0 / MW_KG_PER_MOL  # mol/m³
q_load = brentq(lambda q: cpa_c(LOAD_SALT, q)[0] - C_FEED, 1e-3, 3.0)
henry_load = K_EQ * (CAPACITY / LOAD_SALT) ** NU
SIGMA = CAPACITY / q_load * (1 - (q_load / (henry_load * C_FEED)) ** (1 / NU)) - NU
print(f"CPA binds {q_load * MW_KG_PER_MOL:.0f} g/L solid at the load condition")
print(f"SMA: ν = {NU:.2f}, k_eq = {K_EQ:.3g}, σ = {SIGMA:.1f}")

SMA_VALUES = {"k_eq": K_EQ, "nu": NU, "sigma": SIGMA, "k_kin": 100.0}
SAS_WEAK, SAS_STRONG = 30.0, 300.0  # k_sa, m³/mol (= 1/mM)

MODELS = {
    "SMA": ({"model": "SMA", "parameters": {"mAb": SMA_VALUES}}, "cadet"),
    "SAS weak": (
        {"model": "SAS", "parameters": {"mAb": SMA_VALUES | {"self_association": SAS_WEAK}}},
        "cadet",
    ),
    "SAS strong": (
        {"model": "SAS", "parameters": {"mAb": SMA_VALUES | {"self_association": SAS_STRONG}}},
        "cadet",
    ),
    "Colloidal": (
        {
            "model": "Colloidal",
            "parameters": COLLOIDAL_PARAMETERS,
            "global_parameters": COLLOIDAL_GLOBAL,
        },
        "cadet",
    ),
    "CPA": (
        {"model": "CPA", "parameters": CPA_PARAMETERS, "global_parameters": CPA_GLOBAL},
        "reference",
    ),
}
COLORS = {
    "SMA": "#7f7f7f",
    "SAS weak": "#ff7f0e",
    "SAS strong": "#9467bd",
    "Colloidal": "#d62728",
    "CPA": "#1f77b4",
}

# %% [markdown]
# ## 2. The diagnostic: does q/c rise or fall with c?
#
# The direction of the shift is decided by the isotherm. Plot the **affinity q/c against the
# liquid concentration c**: if it falls, the isotherm is Langmuirian; if it rises,
# anti-Langmuirian. Two salt levels matter:
#
# * at the **load salt** (50 mM), where the protein is taken up;
# * at the **elution salt** (≈ 170 mM), where the peak forms. This one sets the peak shape and
#   the direction of the shift.
#
# Each model's equilibrium is computed from its own code: from the flux functions for SMA/SAS
# and Colloidal (c*(q) where the flux vanishes), and from `equilibrium_factor` for CPA.


# %%
def sma_sas_q(c, salt_mM, k_sa):
    """SMA/SAS equilibrium q (mol/m³ solid): root of K_H (1 + k_sa c)(q̄0/Λ)^ν c − q = 0."""
    henry = K_EQ * (CAPACITY / salt_mM) ** NU
    q_max = CAPACITY / (NU + SIGMA)

    def residual(q):
        return henry * (1 + k_sa * c) * (1 - q / q_max) ** NU * c - q

    return brentq(residual, 0.0, q_max)


def colloidal_c_star(salt_mM, q):
    """Colloidal's equilibrium liquid concentration c*(q), mol/m³.

    With c = 0 and k_kin = 1 the flux k_kin (c − c*) is −c*.
    """
    parameters = {"mAb": COLLOIDAL_PARAMETERS["mAb"] | {"k_kin": 1.0}}
    p = binding_parameters(NAMES, parameters, COLLOIDAL_GLOBAL, None)
    q = np.atleast_1d(q).astype(float)
    liquid = np.vstack([np.full_like(q, salt_mM), np.zeros_like(q)])
    with np.errstate(over="ignore"):
        c_star = -colloidal.flux(liquid, np.vstack([np.zeros_like(q), q]), p)[1]
    return np.where(c_star < 1e3, c_star, np.nan)


def affinity_curves(salt_mM):
    """(c, q/c) per model over the liquid range 0.01–20 g/L."""
    c = np.geomspace(0.01, 20.0, 120) / MW_KG_PER_MOL
    curves = {}
    for label, k_sa in (("SMA", 0.0), ("SAS weak", SAS_WEAK), ("SAS strong", SAS_STRONG)):
        curves[label] = (c, np.array([sma_sas_q(ci, salt_mM, k_sa) for ci in c]) / c)
    q_grid = np.linspace(1e-5, 4.0, 4000)
    for label, c_of_q in (("Colloidal", colloidal_c_star), ("CPA", cpa_c)):
        c_model = c_of_q(salt_mM, q_grid)
        keep = np.isfinite(c_model) & (c_model > c[0]) & (c_model < c[-1])
        curves[label] = (c_model[keep], q_grid[keep] / c_model[keep])
    return curves


fig, axes = plt.subplots(1, 2, figsize=(14, 4))
for ax, salt_mM in zip(axes, (LOAD_SALT, ELUTION_SALT), strict=True):
    for label, (c, affinity) in affinity_curves(salt_mM).items():
        ax.loglog(c * MW_KG_PER_MOL, affinity, color=COLORS[label], lw=2, label=label)
    ax.set(xlabel="c / (g/L)", ylabel="q/c (affinity)", title=f"{salt_mM:.0f} mM salt")
axes[0].legend()
fig.tight_layout()

# %% [markdown]
# * **SMA, Colloidal, CPA**: q/c only falls, so all three are **Langmuirian at every salt**.
# * **SAS** is Langmuirian (weak) or nearly flat (strong) at the load salt, but **rises with c
#   at the elution salt**. It is anti-Langmuirian over a concentration window, then falls once
#   the surface gets crowded: the equilibrium origin of AL-L behaviour.
# * At 50 mM, SMA's Henry coefficient is ~25× below CPA's. SMA's straight log K–log salt line
#   was fitted over the *elution* range and can't follow CPA's curvature down to load salt
#   (Tutorial 06). For gradient elution this barely matters: everything binds at 50 mM, and
#   the peak forms near 170 mM, where the models agree.
#
# ### When does SAS turn anti-Langmuirian?
#
# Linearising the SAS isotherm at low loading gives
#
# d ln(q/c)/dc |₍c→0₎ = k_sa − ν(ν+σ)·K_H(c_s)/Λ
#
# Self-association wins when **k_sa > ν(ν+σ)K_H/Λ**. K_H drops steeply with salt, so above
# some salt every SAS protein starts anti-Langmuirian. That salt is lower the stronger the
# self-association. (Derivation on `docs/theory/isotherms/sas.md`.)

# %%
salt = np.geomspace(40, 400, 200)
threshold = NU * (NU + SIGMA) * K_EQ * (CAPACITY / salt) ** NU / CAPACITY
fig, ax = plt.subplots(figsize=(7.5, 3.8))
ax.loglog(salt, threshold, "k", lw=2, label="threshold ν(ν+σ)K_H/Λ")
for label, k_sa in (("SAS weak", SAS_WEAK), ("SAS strong", SAS_STRONG)):
    ax.axhline(k_sa, color=COLORS[label], ls="--", label=f"{label}: k_sa = {k_sa:g}")
    crossing = salt[np.argmax(threshold < k_sa)]
    ax.axvline(crossing, color=COLORS[label], ls=":", lw=1)
    print(f"{label}: anti-Langmuirian start above ≈ {crossing:.0f} mM salt")
ax.axvspan(LOAD_SALT * 0.95, LOAD_SALT * 1.05, color="grey", alpha=0.2)
ax.text(LOAD_SALT * 1.07, threshold.max() * 2, "load salt", color="grey")
ax.set_xticks([40, 50, 100, 200, 400])
ax.xaxis.set_major_formatter(plt.ScalarFormatter())
ax.xaxis.set_minor_formatter(plt.NullFormatter())
ax.set(
    xlabel="salt / mM",
    ylabel="k_sa / (m³/mol)",
    title="SAS: anti-Langmuirian where the dashed line is above the black",
)
ax.legend(loc="lower left")
fig.tight_layout()

# %% [markdown]
# ## 3. Load series: linear gradient elution at 0.5 to 40 g/L
#
# Same column and method as Tutorial 06: LRMP, 10 cm × 0.5 cm, 150 cm/h. Load 5 g/L mAb at
# 50 mM, wash 3 CV, then a 20 CV gradient to 500 mM. SMA, SAS and Colloidal run on CADET, CPA
# on the reference backend.


# %%
def make_project(model, load_g_per_L):
    isotherm, backend = MODELS[model]
    buffer = {"salt": f"{LOAD_SALT} mM"}
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
                "ionic_capacity": f"{CAPACITY} mM",
                "col_porosity": 0.35,
                "par_porosity": 0.6,
            },
            "isotherm": isotherm,
            "method": [
                {"name": "equilibrate", "duration": "1 CV", "flow": "150 cm/h", "inlet": buffer},
                {
                    "name": "load",
                    "duration": f"{load_g_per_L / 5.0} CV",
                    "flow": "150 cm/h",
                    "inlet": buffer | {"mAb": "5 g/L"},
                },  # fmt: skip
                {"name": "wash", "duration": "3 CV", "flow": "150 cm/h", "inlet": buffer},
                {
                    "name": "elute",
                    "duration": "20 CV",
                    "flow": "150 cm/h",
                    "inlet": buffer,
                    "inlet_end": {"salt": "500 mM"},
                },  # fmt: skip
            ],
            "solver": {"backend": backend, "n_col": 60},
        }
    )


def gradient_cv(project, result):
    """Volume in CV counted from the start of the gradient."""
    start_s = project.resolved_method()[-1].start_s
    return result.column_volumes - np.interp(start_s, result.time_s, result.column_volumes)


LOADS = [0.5, 2, 5, 10, 20, 40]
runs, seconds = {}, {}
for model in MODELS:
    start = time.perf_counter()
    for load in LOADS:
        project = make_project(model, load)
        runs[model, load] = (project, api.simulate(project))
    seconds[model] = time.perf_counter() - start
    print(f"{model:<11} {len(LOADS)} runs in {seconds[model]:5.1f} s")

# %%
fig, axes = plt.subplots(1, len(MODELS), figsize=(19, 3.8), sharey=True)
shades = plt.cm.viridis(np.linspace(0.0, 0.85, len(LOADS)))
for ax, model in zip(axes, MODELS, strict=True):
    for shade, load in zip(shades, LOADS, strict=True):
        project, result = runs[model, load]
        protein = result.outlet("mAb")
        ax.plot(
            gradient_cv(project, result),
            protein / protein.max(),
            color=shade,
            label=f"{load:g} g/L",
        )
    ax.set(xlim=(1, 13), xlabel="CV from gradient start", title=model)
axes[0].set_ylabel("mAb / peak maximum")
axes[-1].legend(title="load", fontsize=8)
fig.tight_layout()

# %% [markdown]
# Normalised to their maxima, the peaks show the direction of the shift and the peak shape at
# a glance:
#
# * **SMA, Colloidal, CPA**: every increase in load moves the peak to the left, with a sharp
#   front and a long tail. Textbook Langmuirian.
# * **SAS weak**: the peak first moves *right* at low loads, then left: a mild AL-L pattern.
# * **SAS strong**: a pronounced anti-Langmuirian regime. Fronting peaks with sharp rears
#   move to later elution up to ~10 g/L (by ~60 mM in elution salt). Then the peak moves back
#   and turns Langmuirian (tailing) at high load, though it still elutes later than at 0.5
#   g/L. This is the AL-L picture described by Seelinger et al. and Koch et al.
#
# ### Elution salt and peak asymmetry against load


# %%
def peak_statistics(project, result):
    """Salt at the peak maximum (mM) and asymmetry factor A_s at 10 % height.

    A_s > 1: tailing (Langmuirian); A_s < 1: fronting (anti-Langmuirian).
    """
    x = gradient_cv(project, result)
    protein = result.outlet("mAb")
    top = np.argmax(protein)
    return result.outlet("salt")[top], api.asymmetry_factor(x, protein, window=(0.0, 20.0))


fig, axes = plt.subplots(1, 2, figsize=(13, 3.9))
for model in MODELS:
    statistics = np.array([peak_statistics(*runs[model, load]) for load in LOADS])
    axes[0].plot(LOADS, statistics[:, 0], "o-", color=COLORS[model], label=model)
    axes[1].plot(LOADS, statistics[:, 1], "o-", color=COLORS[model], label=model)
axes[0].set(
    xscale="log",
    xlabel="load / (g/L column)",
    ylabel="salt at peak maximum / mM",
    title="Retention shift with load",
)
axes[1].axhline(1.0, color="k", lw=0.8)
axes[1].set(
    xscale="log",
    xlabel="load / (g/L column)",
    ylabel="asymmetry factor A_s (10 %)",
    title="Peak shape: >1 tailing, <1 fronting",
)
axes[0].legend()
fig.tight_layout()

for model in MODELS:
    salts = [peak_statistics(*runs[model, load])[0] for load in LOADS]
    turn = LOADS[int(np.argmax(salts))]
    direction = (
        "always earlier (Langmuirian)"
        if turn == LOADS[0]
        else f"later up to {turn:g} g/L, then earlier (AL-L)"
    )
    print(
        f"{model:<11} elution salt {salts[0]:5.1f} → {max(salts):5.1f} → {salts[-1]:5.1f} mM: "
        f"{direction}"
    )

# %% [markdown]
# The asymmetry factor makes the classification quantitative. Langmuirian peaks tail
# (A_s > 1); the anti-Langmuirian SAS peaks front (A_s < 1) in the load range where they move
# to later elution. Measured on real data, **elution salt and A_s against load** are the
# fingerprint that tells you which class of isotherm you need.
#
# ## 4. Why the colloidal models can't be anti-Langmuirian (as implemented)
#
# Both colloidal isotherms in ChromSim are Langmuirian for any parameter set, for structural
# reasons:
#
# * **CPA**: q*/c = K_H · B(Θ) · e^(−u_lat/k_BT). B(Θ) ≤ 1 (excluded area) and u_lat ≥ 0
#   (screened *repulsion* between like-charged proteins: β_ij ∝ Z_lat,i·Z_lat,j > 0). Both
#   factors fall as q rises, so q/c can only fall.
# * **Colloidal**: c* = q·exp(S − ln K_e) with S ∝ Σ q_j √(b_pp,i·b_pp,j)·… ≥ 0, so S rises
#   with q. Even a *negative* b_pp doesn't help: CADET uses √(b_pp,i·b_pp,j), and for a single
#   protein b_pp² is positive, so the sign is lost.

# %%
negative = {"mAb": COLLOIDAL_PARAMETERS["mAb"] | {"b_pp_salt_factor": -153.310, "k_kin": 1.0}}
p_negative = binding_parameters(NAMES, negative, COLLOIDAL_GLOBAL, None)
q_test = np.array([0.5, 1.5, 2.5])
liquid = np.vstack([np.full(3, ELUTION_SALT), np.zeros(3)])
c_negative = -colloidal.flux(liquid, np.vstack([np.zeros(3), q_test]), p_negative)[1]
print("Colloidal c*(q), b_pp > 0:", colloidal_c_star(ELUTION_SALT, q_test))
print("Colloidal c*(q), b_pp < 0:", c_negative, "→ identical: the sign of b_pp is lost")

# %% [markdown]
# Colloid theory itself *can* produce anti-Langmuirian adsorption. It needs a net
# **attractive** protein–protein interaction on the surface: short-range van der Waals or
# hydrophobic contacts, attraction between oppositely charged patches, or a negative second
# virial coefficient B₂₂ in solution. That is the same physics behind the self-association
# that SAS describes empirically. Adding an attractive lateral term (u_lat < 0 at short range)
# to CPA would be the colloidal route to AL-L behaviour. It is not part of the published CPA
# formulation that ChromSim implements.
#
# ## 5. Practical notes
#
# * **Diagnose before you model.** Run LGEs at 3–5 loads. Peaks moving later and fronting mean
#   an anti-Langmuirian isotherm, and SMA, Colloidal and CPA can't fit them, whatever the
#   parameters.
# * **Calibrating SAS.** Low-load LGEs give k_eq and ν, exactly as for SMA (Yamamoto, Tutorial
#   03). k_sa and σ then come from the load series: k_sa sets *where* AL behaviour starts, σ
#   *where it ends*. They correlate, so fit them together on several loads, not one.
# * **Salt matters.** The criterion k_sa > ν(ν+σ)K_H/Λ means AL behaviour shows up at elution
#   salt even when the isotherm measured at load salt (batch uptake) looks perfectly
#   Langmuirian. Batch isotherms at load conditions alone can miss it.
# * **pH.** Seelinger et al. needed a *pH-dependent* SAS. ChromSim's SAS is salt-only for
#   now; at constant pH you can calibrate one parameter set per pH.
# * **Kinetics.** Here SMA/SAS use k_kin = 100 s⁻¹ (near equilibrium). Slow binding kinetics
#   also make peaks front and broaden, so don't mistake them for anti-Langmuirian equilibrium.
#   Check by changing the flow rate: kinetic fronting changes with residence time, equilibrium
#   AL behaviour does not.

# %% [markdown]
# ## Key takeaways
#
# 1. **Langmuirian** (concave isotherm, q/c falls): peaks move to **earlier** elution and tail
#    with load. SMA, Colloidal and CPA are always Langmuirian, for structural reasons.
# 2. **Anti-Langmuirian** (convex isotherm, q/c rises): peaks move **later** and front. In
#    ChromSim, SAS reproduces it, and with crowding at high load the AL-L transition seen for
#    bispecifics, Fc-fusions and polypeptides.
# 3. SAS is anti-Langmuirian where **k_sa > ν(ν+σ)K_H/Λ**: typically at elution salt, not at
#    load salt.
# 4. The fingerprint in data: **elution salt and asymmetry factor against load**.
# 5. A colloidal model with attractive lateral interactions would be the mechanistic
#    counterpart of SAS. ChromSim doesn't have one yet.
