# QTLshuffle

Do cis-eQTL found in a narrow window of circadian time disappear when the rest of the day is added?

## Hypothesis

A genotype × time-of-day (G×time) effect can look like an ordinary eQTL when all samples come from a narrow window of internal time. Under a cosinor model, the genotype effect measured in a window centred on t* is approximately

β(t*) ≈ b + λ (γ1 cos ωt* + γ2 sin ωt*),  λ = sinc(ωw/2) ≈ 0.97 for w = 3 h,

where b is the 24 h mean effect and γ1, γ2 the genotype effect on the rhythm. If b ≈ 0 (pure phase or amplitude effect), the hit is real in the window and absent over 24 h. It is a circadian-dependent eQTL (cd-eQTL), not a false positive.

## Design

| Item | Choice |
|---|---|
| Data | GTEx v10 expression and covariates, GTEx WGS genotypes, donor internal phase (DIP, CHIRAL) |
| Tissues | the `N_TISSUES` = 10 tissues with the most donors that have a DIP |
| Windows | 8 non-overlapping 3 h DIP windows, fixed in advance; windows with < 60 donors are skipped |
| Genes | all expressed genes, minus the genes CHIRAL uses to infer DIP |
| Variants | ±1 Mb of the TSS; MAF ≥ 1 % in the tissue; minor allele count ≥ 10 in each donor subset tested |
| Model | `y ~ g + sex + platform + pcr + PC1-5 + 10 PEER factors`; PEER recomputed in every donor subset |
| eGene | lead p × Meff (eigenMT-style), then BH across genes, q < 0.05 |
| Re-test | the window lead variant in the out-of-window donors (PEER recomputed there) |
| Vanishing hit | p_out > 0.05 and z = (β_win − β_out) / √(SE_win² + SE_out²) with p < 0.05 (Cochran's Q, K = 2) |
| Null | DIP permuted across donors within tissue, 1000 times; every step re-run |
| Test | empirical p = (1 + #{null ≥ observed}) / (1 + B), for nHits, nVanish, and nVanish in rhythmic vs non-rhythmic genes |

Expression is the GTEx inverse-normal transformed matrix, computed once across all donors.

## Run order (Bunya)

```bash
mkdir -p logs
sbatch 01_prepareGenotypes.sh                           # per-chromosome dosages for donors with a DIP
sbatch --dependency=afterok:<01> 02_prepareTissues.sh   # array 1-10: tissue objects, cis map, Meff
sbatch --dependency=afterok:<02> 03_windowScan.sh       # array 1-510: tissue x block of 20 permutations
Rscript 04_summarise.R                                  # empirical tests and null histograms
```

All paths and constants are in `00_config.R`. Each input path can be overridden by an environment variable prefixed `QTLSHUFFLE_` (see `tests/run_toy_test.sh`).

## Outputs (`Results/`)

| Folder | Content |
|---|---|
| `genotypes/` | `<chr>.rds`: donor × variant dosages and variant info |
| `tissues/` | `selectedTissues.csv`; `<tissue>.rds`: expression, fixed covariates, DIP, cis map, Meff, rhythmic labels |
| `windowCounts/` | `<tissue>_block<b>.csv`: one row per permutation × window (nW, nO, nGenes, nHits, nVanish, rhythmic split, nOutUntestable) |
| `windowHits/` | `<tissue>_observed.csv`: hit-level table for the observed DIP |
| `summary/` | `overall.csv`, `perTissue.csv`, null histograms, `missingPermutations.csv` if any task is incomplete |

## Test

`bash tests/run_toy_test.sh` simulates two toy tissues (phase, amplitude and mean eQTL, rhythmic-only and null genes), runs 01 to 04 with 39 permutations, and checks that

- no permutation is missing, the CHIRAL stand-in gene is dropped, and 8 fixed covariates are kept;
- observed vanishing hits exceed the 95th percentile of the null;
- mean eQTL rarely vanish (≤ 20 %).

Toy result: 81 hits (null mean 31), 26 vanishing (null mean 2.8, empirical p = 0.025, the minimum with 39 permutations); 17/29 phase and 8/17 amplitude hits vanish, 0/34 mean eQTL hits do. The test uses PCA in place of PEER (`PEER_BACKEND = 'pca_test_only'`) and an `.rds` genotype file in place of plink.

## Not yet tested

- The `peer` R package on Bunya and its run time. PEER is the critical path: 10 tissues × 1001 DIP sets × 8 windows × 2 subsets ≈ 160 000 fits. Time one `03` task before submitting the array.
- `genio::read_plink` on the real genotype files, and the real input paths.
- `CHIRAL_GENES`: the file listing the CHIRAL input genes does not exist yet. Until it does, no gene is excluded and a warning is printed.

## Requirements

R ≥ 4.4 with tidyverse, data.table, R.utils, genio, peer.
