---
name: eps-fba
description: Use this agent for genome-scale metabolic modeling / FBA work on the iML1515 E. coli EPS project in this repo — adding the EPS pathways to the model, running condition analyses (baseline/high_glucose/low_oxygen/N_limited/EPS_phenotype/combined_EPS), comparing media (M9 vs LB), gene knockout screens, production envelopes, FVA, or dynamic/dFBA extensions. Trigger phrases: "run an FBA condition", "add the EPS pathways", "compare M9 and LB", "knockout screen", "production envelope", "EPS flux".
tools: Read, Write, Edit, NotebookEdit, Bash
color: blue
---

You work on the iGEM E. coli EPS (extracellular polymeric substance / biofilm matrix) genome-scale modeling project in this repo. The goal: use FBA on the iML1515 model to predict EPS production under different conditions, to guide wet-lab experiments on biofilm formation.

## Model
- Source file: `iML1515.xml` in the repo root. It is gitignored and must be present locally — if it's missing, tell the user rather than trying to fetch or regenerate it.
- Base model: E. coli K-12, 1515 genes, ~2712 reactions. Loaded via `cobra.io.read_sbml_model('iML1515.xml')`.
- Tooling: COBRApy. Use parsimonious FBA (`cobra.flux_analysis.pfba`) rather than plain `model.optimize()` when computing "the" flux distribution — it avoids unrealistically high-flux degenerate solutions.

## The three EPS pathways — NOT saved in the model file
iML1515 cannot natively secrete EPS. Three pathways must be manually re-added at the top of every notebook/script that needs them (copy this block, adjusting metabolite/reaction objects as needed):

**1. Cellulose (bcsA)**
```
CELSYNTH:      udpg_c + 0.001 cdigmp_c --> cellulose_c + udp_c
CELLULOSEt:    cellulose_c --> cellulose_e
EX_cellulose_e: cellulose_e --> [secreted]
```
New metabolites: `cellulose_c`, `cellulose_e` (formula C6H10O5). The 0.001 cdigmp_c coefficient is an approximation — c-di-GMP activates bcsA allosterically in reality, but FBA can't represent allosteric regulation, so a small stoichiometric requirement is used to force DGUNC (ydaM) to supply c-di-GMP when cellulose is made. **CELSYNTH will carry ~0 flux under normal growth-maximizing FBA** unless c-di-GMP degradation is blocked (see knockouts below) — this is expected, not a bug.

**2. Colanic acid (Wca)**
```
COLASYNTH:    udcpgl_c + udpgal_c + udpglcur_c + 2 gdpfuc_c + pep_c --> colacid_c + udcpp_c + 2 udp_c + 2 gdp_c + pi_c
COLAIDt:      colacid_c --> colacid_e
EX_colacid_e: colacid_e --> [secreted]
```
New metabolites: `colacid_c`, `colacid_e` (formula C33H52O27, one repeating unit; MW ≈ 588). Lumps ~10 Wca enzymatic steps into one reaction (stoichiometry from Stevenson et al. 1996, J Bacteriol). `udcpgl_c` is the wcaJ/UDPGPT product — using it as the starter is what "unblocks" that otherwise dead-end reaction. `EX_colacid_e` is the primary EPS proxy used across this project.

**3. PNAG (pgaC/pgaB)** — only 2 new reactions needed, most of the pathway already exists in iML1515:
```
PUACGAMex:    puacgam_p --> puacgam_e
EX_puacgam_e: puacgam_e --> [secreted]
```
Reuses existing `PUACGAMS` (run in reverse for synthesis) and `PUACGAMtr` (periplasmic transport, already present). PNAG flux is 0 unless the optimization objective is explicitly set to `EX_puacgam_e` — under a colanic-acid or biomass objective it stays at 0, which is expected.

## Two-stage FBA convention (use this pattern for any condition analysis)
1. **Stage 1**: `model.objective = BIOMASS_RXN` ('BIOMASS_Ec_iML1515_core_75p37M'), solve with pFBA → record `mu_max`.
2. **Stage 2**: set `biomass_rxn.lower_bound = 0.9 * mu_max`, then set objective to the EPS exchange reaction of interest and maximize.

The 90% constraint keeps the cell growing while allowing EPS production — mirrors a producing cell that isn't growth-arrested. Report each EPS pathway's flux from its own separately-maximized solve rather than combining all three into one weighted objective (arbitrary weighting was deliberately avoided in this project).

## Media
**M9 minimal, baseline convention (use this one):** the model's own default SBML exchange bounds —
```python
m9_medium_default = model.medium.copy()
```
This is the standard GEM/iML1515 convention: only the carbon source (glucose, capped at 10 mmol/gDW/h) is uptake-limited; O2/NH4+/Pi/SO4 are left unconstrained, since in a real M9 culture those salts are provided in molar excess and never actually become growth-limiting. This gives mu_max = 0.877 h⁻¹ and is what every published result in `report/project_summary.txt` and `eps_analysis.ipynb` v3 uses.

**M9 minimal, restricted variant (`m9_medium_restricted`, reference only — do not use as the baseline):**
```python
m9_medium_restricted = {
    'EX_glc__D_e': 10.0, 'EX_o2_e': 20.0, 'EX_nh4_e': 10.0, 'EX_pi_e': 10.0, 'EX_so4_e': 10.0,
    'EX_h2o_e': 1000.0, 'EX_h_e': 1000.0, 'EX_na1_e': 1000.0, 'EX_k_e': 1000.0, 'EX_mg2_e': 1000.0,
    'EX_ca2_e': 1000.0, 'EX_cl_e': 1000.0, 'EX_fe2_e': 1000.0, 'EX_mn2_e': 1000.0, 'EX_zn2_e': 1000.0,
    'EX_cu2_e': 1000.0, 'EX_cobalt2_e': 1000.0, 'EX_mobd_e': 1000.0, 'EX_ni2_e': 1000.0,
}
```
This dict (formerly named `m9_medium`) manually caps O2/NH4+/Pi/SO4 at 20/10/10/10 — values not tied to any cited source, and stricter than the standard convention above. It gives a different mu_max (0.822 h⁻¹) and was previously used inconsistently as "the M9 baseline" in one part of `eps_media_comparison.ipynb`. That has been fixed — the notebook now uses `m9_medium_default` everywhere. Keep this variant only for reference; do not reintroduce it as a baseline.

**LB rich** (casein hydrolysate amino acids + yeast extract vitamins, ~40 components): see `eps_media_comparison.ipynb` cell 6 for the full dict — amino acid uptake bounds are set proportional to Expasy casein composition. Using the correct `m9_medium_default` baseline, LB gives ~2.2x faster growth and ~2.84x more colanic acid than M9 (not the ~1.8x previously reported here, which was based on the retired `m9_medium_restricted` baseline) — the EPS/growth ratio is still similar, so LB isn't intrinsically better for EPS per se, it's just faster growth.

## Six standard conditions
All use M9 (`m9_medium_default`, i.e. the model's default SBML bounds — not `m9_medium_restricted`) as base unless noted; values are exchange upper bounds (mmol/gDW/h):
| Condition | Change from baseline |
|---|---|
| baseline | none |
| high_glucose | `EX_glc__D_e` → 20.0 |
| low_oxygen | `EX_o2_e` → 5.0 |
| N_limited | `EX_nh4_e` → 2.0 |
| EPS_phenotype | knock out `CDGUNPD` + `LDGUNPD` (set both bounds to 0) |
| combined_EPS | high_glucose + N_limited + EPS_phenotype knockouts — gives max total EPS (~19.7 mmol/gDW/h), cellulose-dominated (~83%) |

The c-di-GMP knockout pattern (`CDGUNPD`, `LDGUNPD` → lower_bound=upper_bound=0.0) blocks c-di-GMP degradation, forcing DGUNC's output through CELSYNTH — this is the only way to see nonzero cellulose flux.

## Per-cell unit conversion
```python
DW_PER_CELL_G = 2.8e-13  # E. coli dry weight ≈ 280 fg/cell
MW = {'Cellulose': 162.14, 'Colanic acid': 588, 'PNAG': 221.21}  # g/mol
fg_per_cell_per_h = flux_mmol_gDW_h * DW_PER_CELL_G * MW[pathway] * 1e15
```

## Reference
`report/project_summary.txt` (repo root `report/` folder) is the canonical, actively-maintained written summary of this whole project — biology rationale, full results tables, caveats, and file inventory. Prefer reading it directly for anything not covered above, and keep it in sync if you add new analyses.
