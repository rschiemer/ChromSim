# Self-association steric mass action (SAS)

`isotherm: {model: SAS}` · **CADET-native** (`SelfAssociation`) — also runs on the reference
backend (LRM, LRMP).

SMA in which a protein can also adsorb onto an already adsorbed molecule of its own kind
(surface dimerisation, "self-association"; Mollerup 2008). The extra adsorption term grows with
the liquid concentration, so the isotherm can turn **convex** (anti-Langmuirian). Peaks then
elute *later* as the load increases, with fronting and a sharp rear, until the steric term takes
over and the behaviour turns Langmuirian again: the **anti-Langmuir-to-Langmuir (AL-L)**
pattern reported for bispecific antibodies, Fc-fusion proteins and polypeptides on cation
exchangers. Plain SMA can only produce Langmuirian behaviour.

$$
\frac{dq_i}{dt} = k_{a,i}\,(1 + k_{sa,i}\,c_i)\,c_i\left(\frac{\bar q_0}{\Lambda}\right)^{\nu_i}
- k_{d,i}\,q_i\left(\frac{c_0}{\Lambda}\right)^{\nu_i},\qquad
\bar q_0 = \Lambda - \sum_{j\ge1}(\nu_j + \sigma_j)\,q_j
$$

At equilibrium

$$
q_i = k_{eq,i}\,(1 + k_{sa,i}\,c_i)\left(\frac{\bar q_0}{c_0}\right)^{\nu_i} c_i .
$$

| Parameter | Unit | Meaning |
|---|---|---|
| `k_eq` | – | k_a/k_d, as in [SMA](sma.md) |
| `nu` | – | characteristic charge ν |
| `sigma` | – | steric factor σ |
| `self_association` | m³/mol (= 1/mM) | k_sa = k_a2/k_a; 0 gives SMA exactly |
| `k_kin` | 1/s | desorption rate k_d (omit for rapid equilibrium) |

ChromSim passes k_a2 = k_sa·k_a to CADET (`SAI_KA2`) and, as for SMA, sets both reference
concentrations to Λ.

## When is it anti-Langmuirian?

Expanding q/c for a single protein at low loading (q̄_0 ≈ Λ − (ν+σ)q, q ≈ K_H c):

$$
\frac{d\ln(q/c)}{dc}\Big|_{c\to0} = k_{sa} - \frac{\nu(\nu+\sigma)}{\Lambda}\,K_H(c_0),
\qquad K_H = k_{eq}\left(\frac{\Lambda}{c_0}\right)^{\nu}
$$

The isotherm starts convex (anti-Langmuirian) when **k_sa > ν(ν+σ)K_H/Λ**. K_H falls steeply
with salt, so the same protein is Langmuirian at load conditions and anti-Langmuirian at the
salt where it elutes. At higher concentration the steric term always wins, so every SAS isotherm
eventually turns Langmuirian: hence the AL-L transition with increasing load.

Tutorial 07 (`examples/07_langmuir_anti_langmuir.py`) compares SAS with SMA and the colloidal
isotherms over a load sweep.

## Checks in the test suite

- k_sa = 0 gives exactly SMA on CADET (`tests/physics/test_isotherm_equivalences.py`).
- Reference backend vs CADET with strong self-association, LRM and LRMP: ≤ 1 % NRMSE.
- q/c rises, then falls with c when k_sa exceeds the criterion; with k_sa = 0 it only falls.
- The flux vanishes at the equilibrium isotherm.

## References

- Mollerup, J. M. A review of the thermodynamics of protein association to ligands, protein
  adsorption, and adsorption isotherms. *Chem. Eng. Technol.* 31 (2008) 864–874.
- Westerberg, K., Broberg Hansen, E., Degerman, M., Budde Hansen, T. & Nilsson, B. Model-based
  process challenge of an industrial ion-exchange chromatography step. *Chem. Eng. Technol.* 35
  (2012) 183–190.
- Seelinger, F., Wittkopp, F., von Hirschheydt, T. & Frech, C. Anti-Langmuir elution behavior of
  a bispecific monoclonal antibody in cation exchange chromatography: Mechanistic modeling using
  a pH-dependent Self-Association Steric Mass Action isotherm. *J. Chromatogr. A* 1689 (2022).
- Koch, J., Scheps, D., Gunne, M., Boscheinen, O., Hafner, M. & Frech, C. Mechanistic modeling
  and simulation of a complex low and high loading elution behavior of a polypeptide in cation
  exchange chromatography. *J. Sep. Sci.* 45 (2022) 2008–2023. doi:10.1002/jssc.202200098
- CADET documentation, *Self Association*; CADET-Core `SelfAssociationBinding.cpp` (v5.1).
