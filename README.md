# rdoc-hcp-code

Analysis code for the preprint:

> Quah SKL, Pirzada S, Madsen S, Jo B, Uddin LQ, Mumford JA, Barch DM, Gotlib IH, Fair DA, Poldrack RA, Saggar M.
> **Refining RDoC Using Individual-Level Task fMRI Factor Models Reveals Reproducible and Clinically Relevant Brain-Wide Motifs.**
> bioRxiv (2025). https://www.biorxiv.org/content/10.1101/2025.10.13.682124v2
> doi: [10.1101/2025.10.13.682124](https://doi.org/10.1101/2025.10.13.682124)

From the [Brain Dynamics Lab](https://github.com/braindynamicslab), Stanford University.

## Contents

| File | Language | Purpose |
|---|---|---|
| `rdochcp_analysis.Rmd` | R | Main analysis |
| `rdochcp_othermodels.Rmd` | R | Comparison with alternative factor models |
| `rdochcp_individualfa.R` | R | Individual-level factor analysis of task-fMRI contrasts |
| `rdochcp_variablepca.R` | R | Variable-level PCA |
| `rdochcp_plotmultipleviolins.Rmd` | R | Violin plots |
| `rdochcp_mapper_modularity.mlx` | MATLAB | Mapper graph modularity |
| `rdochcp_mapper_modularitybs1.m`, `rdochcp_mapper_modularitybs2.m` | MATLAB | Jackknife modularity for cohorts 1 and 2 |
| `rdochcp_community_corrheatmap.mlx` | MATLAB | Community correlation heatmaps |
| `rdochcp_behavcorr.mlx` | MATLAB | Brain-behavior correlations |

## Data

No data are included in this repository. The scripts expect data files (e.g. cohort lists, HCP task-contrast
maps, behavioral data) at local paths that you will need to set at the top of each script.
Human Connectome Project data, including restricted behavioral variables, can be obtained from the
[HCP](https://www.humanconnectome.org/) subject to its data use terms.

## Mapper functions

The MATLAB scripts depend on Mapper and plotting functions that are not included in this repository
(e.g. `runBDLMapper`, `calMod`, `drawd3graph`, `computeProportionCategorical`, `computeHubness`,
`testModularitySignificance`). `rdochcp_mapper_modularity.mlx` loads precomputed Mapper outputs
(`rdochcp_contrastmapper_subs1.mat` / `subs2.mat`) rather than running Mapper itself.

For Mapper implementations from our lab, see:

- [neuMapper](https://github.com/braindynamicslab/neumapper): a scalable Mapper algorithm for neuroimaging data analysis
- [deMapper](https://github.com/braindynamicslab/demapper): the deconstructed Mapper toolbox for MATLAB

## Citation

If you use this code, please cite the preprint above.
