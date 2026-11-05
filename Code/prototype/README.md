# QTLshuffle prototype (R only)

The full design on one tissue and a small gene set. Use it to settle the experimental design before the transcriptome-wide pipeline in the parent folder.

## Scope

| Item | Choice |
|---|---|
| Tissue | `TISSUE` = Muscle_Skeletal |
| Genes | the `N_GENES` = 50 genes with the smallest 24 h harmonic regression p (cd-eQTL step 01), autosomes |
| Variants | ±1 Mb of the TSS; MAF ≥ 1 % among the tissue's donors; minor allele count ≥ 10 in each donor set tested |
| Windows | 8 non-overlapping 3 h DIP windows; windows with < 60 donors skipped |
| Model | `y ~ g + sex + platform + pcr + PC1-5 + 10 PEER factors`; PEER recomputed on all genes in every donor set |
| eGene | per gene, 10 000 permutations, beta-approximated p with tensorQTL's effective degrees of freedom, then Storey q-value (R `qvalue`, λ = 0.85) across the 50 genes, q < 0.05 |
| 24 h reference | per window, as many donors as in the window drawn from all donors (about 1/8 shared); PEER recomputed; window lead variant re-tested |
| Vanishing hit | p_24h > 0.05 and z = (β_win − β_24h) / √(SE_win² + SE_24h²) with p < 0.05 |
| Null | DIP permuted across donors 1000 times; everything re-run |

## Figure: effect against inferred phase (`04_phasePlot.R`)

One panel per gene. One point per 3 h DIP window: the variant effect re-estimated on that window's donors, with PEER recomputed. The thick bar is ±1 SE (the SD of the estimate) and the thin bar is the 95% CI. The dashed line is the 24 h effect. The blue curve is the cosinor G × time fit on all donors, attenuated for the window width. The script makes two versions, one for each way of choosing the variant:

- `betaByPhase_24hLead.pdf`: the lead variant of a standard 24 h scan;
- `betaByPhase_GxTLead.pdf`: the variant with the strongest G × time interaction.

A pure phase or amplitude QTL has no 24 h marginal effect, so the 24 h scan usually does not pick its variant. The variant is chosen using all donors, the window donors included, so the figure is descriptive and not a test.

## Differences from tensorQTL

- **Permutations.** tensorQTL permutes the covariate-residualised expression and correlates it with the residualised genotypes. Permuted residuals are no longer orthogonal to the covariates, so their correlations have variance ≈ 1/(n − 1) instead of 1/(n − k − 1). With about 19 covariates and about 100 window donors, that null is too narrow. Here, the permuted residuals are residualised again (Freedman–Lane), which restores trueDof ≈ n − k − 1. `tests/null_calibration.R` simulates 2000 genes with no genetic effect (n = 80, 19 covariates, 1000 permutations each):

  | Scheme | p < 0.05 | p < 0.01 | p < 0.001 | KS test vs uniform |
  |---|---|---|---|---|
  | Freedman–Lane (this prototype) | 4.8 % | 0.85 % | 0 % | p = 0.13 |
  | Plain residual permutation (tensorQTL) | 22 % | 5.6 % | 1.0 % | p < 10⁻¹⁵ |

  The 2000 genes share 50 cis variant sets, so these rates are not fully independent. This matters at window sample sizes; at full GTEx sizes the difference is small.
- **q-value with 50 genes.** If no gene has p ≥ 0.85, π0 is estimated as 0 and `qvalue` fails. The window then falls back to BH and records `qMethod = 'BH'`.

## Run order (Bunya)

```bash
cd Code && mkdir -p logs && cd prototype   # logs go to Code/logs
sbatch 01_prepareData.sh                               # donors, 50 genes, cis genotypes read from plink
sbatch --dependency=afterok:<01> 02_windowScan.sh      # array 1-51: blocks of 20 DIP permutations
Rscript 03_summarise.R                                 # empirical p for nHits and nVanish
sbatch --dependency=afterok:<01> 04_phasePlot.sh       # beta of one variant per gene in each window
```

Outputs go to `Results/prototype/`:
- `genes.csv`: the selected genes;
- `windowCounts/block<b>.csv`;
- `windowHits/observed.csv`;
- `summary/`.

## Test

`bash tests/run_toy_test.sh` uses the toy data from `../tests/make_toy_data.R`, with 39 DIP permutations, 1000 permutations per gene, and PCA in place of PEER. It took 2 min on 4 cores.

Result:
- 116 observed hits against a null mean of 0.36 (p = 0.025).
- 100 vanishing hits against a null mean of 0.36.
- 53/63 phase-QTL hits and 43/47 amplitude-QTL hits vanish.
- 6 hits among the rhythmic genes without a genetic effect. Most likely the latent factors (PCA here) absorb the very strong toy cis effects of neighbouring genes, and adjusting for them induces an association.

## Requirements

R 4.4 (`module load r/4.4.2-heavy`) with tidyverse, data.table, R.utils and qvalue (Bioconductor). PEER does not build under R 4.4. It runs in its own R 3.4 conda environment (`PEER_RSCRIPT`, installed by the rhyh2 repo's `Code/peer/install_peer.sh`) through `peer_fit.R`, one process per fit, with R 4.4's environment variables removed. Smoke test: `module purge; /home/uqschauq/software/peer_env/bin/Rscript -e 'library(peer); cat("ok\n")'`. No plink reader package is needed: `readPlinkVariants` reads the cis variants from the .bed directly.
