# Colloidal adsorption models: background and literature

Stoichiometric models such as SMA treat protein binding as an exchange of a fixed number ν of
counter-ions. Colloidal models instead treat the protein as a charged particle and the resin
surface as a charged wall, and derive adsorption from colloid physics: electrostatic
double-layer forces screened by salt, the area a protein excludes on the surface, and repulsion
between adsorbed neighbours. The appeal for process modeling is twofold: parameters with a
physical meaning (radius, charge), and a more realistic shape of the isotherm at high load,
where SMA tends to overestimate nonlinearity (Briskot et al. 2021a).

ChromSim ships two models of this family. Tutorial 06 (`examples/06_colloidal_models.py`)
compares them side by side on synthetic data at low and high load.


| | [Colloidal](isotherms/colloidal.md) (Xu & Lenhoff) | [CPA](isotherms/cpa.md) (Briskot et al.) |
|---|---|---|
| Backend | CADET (native) + reference | reference (Python) |
| Protein–surface term | empirical ln K_e(salt) | double-layer well depth from Z, Γ_L, I |
| Crowding | hexagonal lattice spacing | hard-disc SPT available surface + lattice |
| Protein–protein term | Yukawa b_pp(salt) | Yukawa from Z_lat with DLVO prefactor |
| Salt dependence of screening | κ(salt), restricted (κ_ef ≥ 0) | κ ∝ √I, physical |
| Typical parameters per protein | 5–9 | 4–5 |
| Use when | you need CADET speed, GRM, or the published Lenhoff parameterisations (IEX, Protein A) | you want physically interpretable IEX parameters; LRM/LRMP is enough |

## How the theory developed

1. **Electrostatic retention models (1990s).** Ståhlberg, Jönsson & Horváth (1991) and Roth,
   Unger & Lenhoff (1996) replaced the stoichiometric displacement picture by double-layer
   interaction energies between a charged protein and a charged surface, explaining the
   log k′–log I relationship without a stoichiometric ν.
2. **Colloidal energetics for isotherms.** Oberholzer & Lenhoff (1999) computed full isotherms
   from protein–surface and protein–protein energies, using Brownian dynamics to describe the
   structure of the adsorbed layer. Xu & Lenhoff (2008, 2009) condensed this into an algebraic
   isotherm with a hexagonal adsorbed layer, extended to binary mixtures — the model CADET
   implements as *Multi Component Colloidal*.
3. **Excluded area.** Scaled particle theory for hard discs (Helfand, Frisch & Lebowitz 1961)
   gives the chemical potential of a crowded 2D layer; Chatelier & Minton (1996) and Talbot et
   al. (2000) brought it to protein adsorption as the *available surface function*.
4. **CPA for preparative IEX (Hubbuch group / GoSilico).** Briskot et al. (2020) built a linear
   isotherm from the linearised Poisson–Boltzmann equation with charge regulation of protein and
   ligand. Briskot et al. (2021a) extended it to the nonlinear range with SPT crowding and
   Yukawa lateral repulsion, and compared it with SMA on several mAbs. Briskot et al. (2021b)
   applied the kinetic form to complex elution behaviour (peak splitting, pH effects) of
   industrial mAbs. The derivations are collected in Briskot's thesis (2022, open access).
5. **Applications since 2022.** Parameter estimation workflows for CPA (Lorenz-Cristea et al.
   2024); Protein A capture with the colloidal isotherm (Bhoyar et al. 2024) and with CPA on
   fiber adsorbents (Hahn et al. 2024); mixed-mode chromatography using the CPA surface coverage
   function (Biotechnol. Prog. 2022, doi:10.1002/btpr.3316); CPA vs SMA for parvovirus clearance
   on AEX (Döring et al., KITopen 1000184789); robust IEX modeling with the colloidal isotherm in
   CADET (Kumar et al. 2021).

## Langmuirian and anti-Langmuirian behaviour

Both colloidal isotherms are **Langmuirian for every parameter set**: peaks move to earlier
elution with increasing load. In CPA, the excluded-area factor B(Θ) ≤ 1 and the lateral
term is repulsive (β_ij ∝ Z_lat,i·Z_lat,j > 0). In Colloidal, the lattice term
√(b_pp,i·b_pp,j)·… is non-negative, and even a negative b_pp loses its sign for a single
protein. Anti-Langmuirian elution (later elution and fronting with load), as reported for
bispecific antibodies, Fc-fusion proteins and polypeptides on cation exchangers, needs
attractive protein–protein interactions. In ChromSim the [SAS](isotherms/sas.md) isotherm
describes it empirically. Tutorial 07 (`examples/07_langmuir_anti_langmuir.py`) compares SMA,
SAS, Colloidal and CPA over a load series.

Literature on anti-Langmuirian protein IEX

- Seelinger, F., Wittkopp, F., von Hirschheydt, T. & Frech, C. Anti-Langmuir elution behavior of
  a bispecific monoclonal antibody in cation exchange chromatography: Mechanistic modeling using
  a pH-dependent Self-Association Steric Mass Action isotherm. *J. Chromatogr. A* 1689 (2022).
- Koch, J., Scheps, D., Gunne, M., Boscheinen, O., Hafner, M. & Frech, C. Mechanistic modeling
  and simulation of a complex low and high loading elution behavior of a polypeptide in cation
  exchange chromatography. *J. Sep. Sci.* 45 (2022) 2008–2023.
- Mechanistic modeling of anti-Langmuirian to Langmuirian behavior of Fc-fusion proteins in
  cation exchange chromatography. *J. Chromatogr. A* (2025) (nGIEX model).
- Huuk, T. C., Hahn, T., Doninger, K., Griesbach, J., Hepbildikler, S. & Hubbuch, J. Modeling of
  complex antibody elution behavior under high protein load densities in ion exchange
  chromatography using an asymmetric activity coefficient. *Biotechnol. J.* 12 (2017) 1600336.
- Diedrich, J. et al. Multi-state steric mass action model and case study on complex high loading
  behavior of mAb on ion exchange tentacle resin. *J. Chromatogr. A* 1525 (2017) 60–70.
- Mollerup, J. M. A review of the thermodynamics of protein association to ligands, protein
  adsorption, and adsorption isotherms. *Chem. Eng. Technol.* 31 (2008) 864–874.

## Literature

Core model papers

- Xu, X. & Lenhoff, A. M. Binary adsorption of globular proteins on ion-exchange media.
  *J. Chromatogr. A* 1216 (2009) 6177–6195. doi:10.1016/j.chroma.2009.06.082
- Oberholzer, M. R. & Lenhoff, A. M. Protein adsorption isotherms through colloidal energetics.
  *Langmuir* 15 (1999) 3905–3914. doi:10.1021/la981199k
- Briskot, T., Hahn, T., Huuk, T. & Hubbuch, J. Adsorption of colloidal proteins in
  ion-exchange chromatography under consideration of charge regulation. *J. Chromatogr. A* 1611
  (2020) 460608. doi:10.1016/j.chroma.2019.460608
- Briskot, T., Hahn, T., Huuk, T. & Hubbuch, J. Protein adsorption on ion exchange adsorbers:
  A comparison of a stoichiometric and non-stoichiometric modeling approach. *J. Chromatogr. A*
  1653 (2021) 462397. doi:10.1016/j.chroma.2021.462397
- Briskot, T. et al. Analysis of complex protein elution behavior in preparative ion exchange
  processes using a colloidal particle adsorption model. *J. Chromatogr. A* 1654 (2021) 462439.
- Briskot, T. *Science and Risk-Based Development of Ion Exchange Chromatography.* PhD thesis,
  Karlsruhe Institute of Technology (2022). Open access:
  publikationen.bibliothek.kit.edu/1000144198

Foundations

- Roth, C. M., Unger, K. K. & Lenhoff, A. M. Mechanistic model of retention in protein
  ion-exchange chromatography. *J. Chromatogr. A* 726 (1996) 45–56.
- Ståhlberg, J., Jönsson, B. & Horváth, Cs. Theory for electrostatic interaction chromatography
  of proteins. *Anal. Chem.* 63 (1991) 1867–1874.
- Helfand, E., Frisch, H. L. & Lebowitz, J. L. Theory of the two- and one-dimensional rigid
  sphere fluids. *J. Chem. Phys.* 34 (1961) 1037.
- Chatelier, R. C. & Minton, A. P. Adsorption of globular proteins on locally planar surfaces:
  models for the effect of excluded surface area and aggregation of adsorbed protein on
  adsorption equilibria. *Biophys. J.* 71 (1996) 2367–2374.
- Talbot, J., Tarjus, G., Van Tassel, P. R. & Viot, P. From car parking to protein adsorption:
  an overview of sequential adsorption processes. *Colloids Surf. A* 165 (2000) 287–324.
- Usui, S. Interaction of electrical double layers at constant surface charge. *J. Colloid
  Interface Sci.* 44 (1973) 107–113.

Applications and calibration

- Kumar, V. et al. Robust mechanistic modeling of protein ion-exchange chromatography.
  *J. Chromatogr. A* 1660 (2021) 462669.
- Lorenz-Cristea, O. et al. A systematic approach for estimating colloidal particle adsorption
  model parameters. *J. Chromatogr. A* (2024) 465512. doi:10.1016/j.chroma.2024.465512
- Bhoyar, S., Kumar, V., Foster, M., Xu, X., Traylor, S. J., Guo, J. & Lenhoff, A. M.
  Predictive mechanistic modeling of loading and elution in protein A chromatography.
  *J. Chromatogr. A* (2024) 464558. doi:10.1016/j.chroma.2023.464558
- Hahn, T. et al. Predictive scaling of fiber-based protein A capture chromatography using
  mechanistic modeling. *Biotechnol. Bioeng.* 121 (2024) 2388–2399.
- Mechanistic modeling, simulation, and optimization of mixed-mode chromatography for an
  antibody polishing step. *Biotechnol. Prog.* (2022). doi:10.1002/btpr.3316
- Döring, L., Winderl, J., Müller, M., Kron, M. W. & Hubbuch, J. Comparison of stoichiometric
  and colloidal adsorption models for mechanistic modeling of parvovirus clearance by anion
  exchange chromatography. Open access: publikationen.bibliothek.kit.edu/1000184789
