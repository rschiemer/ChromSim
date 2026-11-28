# 0011 — Colloidal isotherms: Colloidal (CADET) and CPA (Python)

- Status: accepted
- Date: 2026-11-27

## Context
Colloidal adsorption models are a published, non-stoichiometric alternative to SMA for ion
exchange (and Protein A). Two families exist: the colloidal isotherm of Xu & Lenhoff (2009),
which CADET-Core 5 implements (`MULTI_COMPONENT_COLLOIDAL`), and the colloidal particle
adsorption (CPA) model of Briskot et al. (2021), which CADET lacks.

## Decision
1. **`Colloidal`** is a CADET-native descriptor (`MultiComponentColloidal`) that also carries a
   Python `flux` copied from CADET-Core's `ColloidalBinding.cpp`, so it runs on both backends.
   The reference backend matches CADET to ≤ 0.3 % NRMSE on a two-protein nonlinear case.
2. **Parameter mapping follows the C++ source, not CADET's interface docs.** The source computes
   ln K_e = c^(−POWEXP)·POWFACT + exp(c·EXPARGMULT)·EXPFACT (same for b_pp with +POWEXP);
   the docs describe POWEXP/POWFACT the other way round. ChromSim uses descriptive names
   (`ln_ke_salt_factor`, `ln_ke_salt_exponent`, …) and the backend cross-check pins the mapping.
3. **k_kin keeps CADET's meaning for `Colloidal`**: dq/dt = k_kin (c − c*), i.e. a rate on the
   liquid-phase driving force, not the desorption rate of ChromSim's `k_eq`/`k_kin` convention
   (ADR 0005). The colloidal model has no separable k_eq, so a translation would hide CADET's
   behaviour and make published parameters harder to reuse. Documented on the theory page.
4. **κ exponent ≥ 0.** CADET-Process 0.12 declares `kappa_exponential` unsigned, although a
   Debye length ∝ c^(−1/2) needs κ_ef = −0.5. ChromSim keeps the same bound (CLAUDE.md rule 4)
   and documents the workaround. Revisit when CADET-Process lifts it (worth an upstream issue).
5. **`CPA`** is a Python isotherm (reference backend). It implements Briskot et al. (2021) with
   documented simplifications: strong ligand without charge regulation, no pH dependence,
   constant k_kin, one specific surface A_s. Each physical step is a small public function with
   its own textbook-limit test.

## Consequences
- Colloidal models are available today for IEX studies; CPA is limited to LRM/LRMP and to the
  speed of the reference backend. Its route to production is ADR 0008 (CADET-Core 6
  data-driven binding or upstreaming to C++).
- Extending CPA with pH (Z(pH), Δ(pH)) needs a second modifier component (`n_modifiers = 2`),
  which is a new isotherm name or a breaking change; decide when needed.
