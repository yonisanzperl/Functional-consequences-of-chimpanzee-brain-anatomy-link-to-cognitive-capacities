# Code: whole-brain modelling of chimpanzee brain dynamics and cognition

MATLAB code to reproduce the simulations, sweeps, statistics and figures of
*"Functional-consequences-of-chimpanzee-brain-anatomy-link-to-cognitive-capacities"* (Sanz Perl et al.).


## Folder structure

```
functions/                    helper functions used by several scripts
01_chimp_model/               chimpanzee model: coupling sweep, measures at G_opt, surrogate connectomes
02_human_data_and_model/      HCP 7T: empirical measures, empirical FC, human model working point
03_cross_species_model/       comparative human-chimpanzee dataset: SC harmonisation (V0/V1/V2), sweep + measures
figures/Fig1 ... Fig6         statistics and figure panels, one folder per figure

```

## Pipeline (run in this order)

| Step | Script | Input | Output | Used by |
|---|---|---|---|---|
| 1. Chimp coupling sweep (metastability vs G, 0-0.01, step 1e-4) | `01_chimp_model/CHIMP_CWAS_metastability_empirica.m` | `CHIMP_CWAS_DATA_FORDECO.mat` | `meta_G_CHIMP_CWAS_Cmeannorm_narrowfilt2.mat` | 2, 3, Fig 2-4 |
| 2. Chimp measures at G_opt (200 simulations; synchrony, metastability, Ignition, edge turbulence, entropy eFCD, FC) | `01_chimp_model/functional_measures_CHIMP_MVH.m` | step 1 + data | `functionalCHIM_MVH_filtnarr_LZ_granger4_200_1200_REAL_0926_001.mat` | Fig 2-4 |
| 3. Surrogate connectomes (degree-preserving and degree+strength-preserving, 50 sets each, simulated at the real G_opt) | `01_chimp_model/Fig2_surrogate_connectomes.m` | step 1 + data | folder `surrogates_sweepfreq/` | Fig 2c |
| 4. Human empirical FC and metastability in the model band (0.02-0.03 Hz) | `02_human_data_and_model/HCP7T_empirical_FC_chimpband.m` | `hcp7t_rfMRI_REST1_PA_schaefer1000.mat`, `schaefer1000to100.mat` | `HCP7T_empirical_FC_chimpband.mat` | 5 |
| 5. Human model working point (group SC, chimp settings) | `02_human_data_and_model/model_Humans_HCP7T_asChimp.m` | `SC_schaefer100_17Networks_32fold_groupconnectome_2mm_symm.mat`, step 4 | `results_model_Humans_HCP7T_asCHIMP_bandFC.mat` | Fig 1c |
| 6. Human dynamical measures from empirical 7T fMRI (181 participants) | `02_human_data_and_model/HCP7T_measures_for_chimp_comparison.m` | 7T time series, `schaefer1000to100.mat`, `hcpbehaviouraldata.mat` | `HCP7T_measures_for_chimp_comparison.mat` | Fig 5 |
| 7. Cross-species SC harmonisation (V0 original, V1 density-matched, V2 density + weight-matched) | `03_cross_species_model/HumChimp_SC_harmonize.m` | `connectivity_HUMAN_BB50_v7.mat`, `connectivity_CHIMP_BB50_v7.mat` | `SC_HumChimp_harmonized.mat` | 8 |
| 8. Cross-species sweep and measures (identical model in both species) | `03_cross_species_model/HumChimp_sweep_and_measures.m` | step 7 | `HumChimp_V0/V1/V2_sweep_measures.mat` | Fig 6 |

## Figures

| Figure | Script(s) | Needs |
|---|---|---|
| Fig 1 (a, b, d, e, f: schematic, assembled in Illustrator) | panel c: `figures/Fig1/Fig5a_working_point.m`; panel d from step 1; panels e-f are Fig 5d and Fig 4b | step 5 |
| Fig 2 (dynamics vs cognition, robustness, surrogates) | `figures/Fig2/Fig2_REAL0926_full.m` (figure + stats); `Fig2_partial_correlations_controls.m` (partial r, FDR); `Fig2c_surrogate_rank_test.m` (rank test) | steps 1-3 |
| Fig 3 (SC vs inferred FC network topology) | `figures/Fig3/Fig3_network_metrics_SC_FC.m` (sign test, Williams tests, G_opt control) | steps 1-2 |
| Fig 4 (held-out explanation) | `figures/Fig4/Fig4_heldout_explanation.m` | step 2 + `results_Fig3_network_metrics.mat` (Fig 3) |
| Fig 5 (humans, PMAT24) | `figures/Fig5/Fig5_PMAT_figure.m` (figure); `HCP7T_PMAT_vs_chimp_measures.m` (full stats, species comparison) | step 6, `HCP7T_sex.mat` |
| Fig 6 (cross-species, model level) | `figures/Fig6/Fig6_alt_model_level.m` | step 8 |


## Dependencies

- MATLAB (Parallel Computing Toolbox for the `parfor` loops in steps 3, 8 and Fig 3).
- Brain Connectivity Toolbox (BCT): https://sites.google.com/site/bctnet/ (copies of `modularity_und`, `null_model_und_sign`, `randmio_und_signed` are in `functions/`).
- `cmocean.m` (colour maps) and `plot_mean_std.m` are in `functions/`.
- Python 3 with numpy, scipy and matplotlib for the supplementary figures.
