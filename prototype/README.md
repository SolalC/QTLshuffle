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

## Differences from tensorQTL

- **Permutations.** tensorQTL permutes the covariate-residualised expression and correlates it with the residualised genotypes. Permuted residuals are no longer orthogonal to the covariates, so their correlations have variance ≈ 1/(n − 1) instead of 1/(n − k − 1). With about 19 covariates and about 100 window donors, that null is too narrow. In a simulation (n = 80, 19 covariates, no genetic effect), this scheme gave beta-approximated p < 0.05 for 15 % of genes. Here, the permuted residuals are residualised again (Freedman–Lane), which gave 8 % in 200 null genes, with trueDof ≈ n − k − 1. This matters at window sample sizes; at full GTEx sizes the difference is small.
- **q-value with 50 genes.** If no gene has p ≥ 0.85, π0 is estimated as 0 and `qvalue` fails. The window then falls back to BH and records `qMethod = 'BH'`.

## Run order (Bunya)

```bash
cd prototype && mkdir -p logs
sbatch 01_prepareData.sh                               # donors, 50 genes, cis genotypes read from plink
sbatch --dependency=afterok:<01> 02_windowScan.sh      # array 1-51: blocks of 20 DIP permutations
Rscript 03_summarise.R                                 # empirical p for nHits and nVanish
```

Outputs go to `Results_prototype/`:
- `genes.csv`: the selected genes;
- `windowCounts/block<b>.csv`;
- `windowHits/observed.csv`;
- `summary/`.

## Test

`bash tests/run_toy_test.sh` uses the toy data from `../tests/make_toy_data.R`, with 39 DIP permutations, 1000 permutations per gene, and PCA in place of PEER. It took 2 min on 4 cores.

Result:
- 41 observed hits against a null mean of 0.36 (p = 0.025).
- 36 vanishing hits against a null mean of 0.33.
- 21/25 phase-QTL hits and 15/16 amplitude-QTL hits vanish.
- No hits among the rhythmic genes without a genetic effect.

## Requirements

R ≥ 4.4 with tidyverse, data.table, R.utils, peer, and qvalue (Bioconductor). No plink reader package is needed: `readPlinkVariants` reads the cis variants from the .bed directly.
