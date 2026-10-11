# The Yamamoto method for ChromSim's model

The Yamamoto method estimates the SMA parameters ν and k_eq from linear-gradient elution (LGE)
runs at low load. This page derives the relation for the column model ChromSim simulates, so
that the definitions match the simulator exactly.

## Derivation

At low load, the SMA isotherm is linear in the protein concentration, with a distribution
coefficient that depends on salt:

$$
K(c_s) = \frac{q}{c_p} = k_{eq}\left(\frac{\Lambda}{c_s}\right)^{\nu}
$$

Here q is bound protein per solid volume, c_p is protein in the pore liquid, c_s is salt and Λ is
the ionic capacity. With solid volume fraction φ = (1 − ε_c)(1 − ε_p) = 1 − ε_t, an equilibrium
protein zone at relative position x = z/L moves as

$$
\frac{dx}{dV} = \frac{1}{V_c(\varepsilon_t + \varphi K)}
$$

where V is the pumped volume. The salt at the zone was fed when V − ε_t V_c x had been pumped.
In a linear gradient of slope g (salt per pumped volume) the zone therefore sees
c_s = c_0 + g W, with W = V − ε_t V_c x − V_{start}. Differentiating,
dW = dV − ε_t V_c dx = φ K V_c dx, and so

$$
\frac{dc_s}{K(c_s)} = g\,\varphi\,V_c\,dx
\quad\Rightarrow\quad
\int_{c_0}^{I_R}\frac{dc_s}{K(c_s)} = g\,(1-\varepsilon_t)\,V_c \equiv GH .
$$

I_R is the salt concentration at the column outlet when the zone leaves the column, which is
when the peak maximum elutes. Inserting K(c_s) and assuming I_R ≫ c_0:

$$
GH = \frac{I_R^{\nu+1}}{(\nu+1)\,k_{eq}\,\Lambda^{\nu}}
\quad\Longrightarrow\quad
\log_{10} GH = (\nu+1)\log_{10} I_R - \log_{10}\!\big((\nu+1)\,k_{eq}\,\Lambda^{\nu}\big).
$$

This is the classical result of Yamamoto et al. (1983, 1988), with G·H = g·(V_t − V_0) in their
notation and V_t − V_0 = (1 − ε_t)V_c.

## Validation

`tests/workflows/test_calibration_recovery.py` simulates four low-load LGE runs with CADET (LRM,
rapid equilibrium) and recovers ν to 0.06 % and k_eq to 0.9 %. With film and pore diffusion
(GRM), peaks elute slightly late: ν comes out a few percent high and k_eq lower. That is why the
Yamamoto result is a starting point for an inverse fit, not the final answer.

## Assumptions

- low load, so steric shielding (σ) does not matter;
- protein and salt access the same pore volume;
- I_R ≫ c_0, the starting salt;
- mass transfer does not shift the peak maximum much.

References: Yamamoto, Nakanishi, Matsuno & Kamijubo, *Biotechnol. Bioeng.* 25 (1983) 1465;
Yamamoto, Nakanishi & Matsuno, *Ion-Exchange Chromatography of Proteins*, Dekker 1988;
Pedersen, Mollerup, Hansen & Jungbauer, *J. Chromatogr. B* 790 (2003) 161.
