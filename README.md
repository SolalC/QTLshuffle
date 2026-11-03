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
| Genes | all expressed genes |
| Variants | ±1 Mb of the TSS; MAF ≥ 1 % in the tissue; minor allele count ≥ 10 in each donor subset tested |
| Model | `y ~ g + sex + platform + pcr + PC1-5 + 10 PEER factors`; PEER recomputed in every donor subset |
| eGene | tensorQTL `map_cis` as in GTEx: 10 000 genotype permutations per gene, beta-approximated p, Storey q-value (λ = 0.85; computed in Python with the `qvalue` formula, so rpy2 is not needed), q < 0.05 |
| 24 h reference | per window, n_W donors drawn at random from all donors over 24 h (about 1/8 are window donors); PEER recomputed there; the window lead variant re-tested |
| Vanishing hit | p_24h > 0.05 and z = (β_win − β_24h) / √(SE_win² + SE_24h²) with p < 0.05 (Cochran's Q, K = 2; conservative because of the shared donors) |
| Null | DIP permuted across donors within tissue, 1000 times; every step re-run |
| Test | empirical p = (1 + #{null ≥ observed}) / (1 + B), for nHits, nVanish, and nVanish in rhythmic vs non-rhythmic genes |

Expression is the GTEx inverse-normal transformed matrix, computed once across all donors.

## Run order (Bunya)

```bash
mkdir -p logs
sbatch 01_prepareTissues.sh                             # array 1-10 (R): tissues, rhythmic labels, config.json
sbatch --dependency=afterok:<01> 02_peerCovariates.sh   # array 1-510 (R, CPU): PEER per window and 24 h reference
sbatch --dependency=afterok:<02> 03_tensorQTLscan.sh    # array 1-510 (Python, GPU): tensorQTL scans and re-tests
Rscript 04_summarise.R                                  # empirical tests and null histograms
```

Arrays 02 and 03 use the same task numbering (tissue × block of 20 permutations), so task k of 03 can start once task k of 02 has finished. All paths and constants are in `00_config.R`; 01 copies what the Python step needs to `Results/config.json`. Each input path can be overridden by an environment variable prefixed `QTLSHUFFLE_` (see `tests/run_toy_test.sh`).

## Outputs (`Results/`)

| Folder | Content |
|---|---|
| `config.json` | paths and constants for 03 |
| `tissues/` | `selectedTissues.csv`; `<tissue>.rds`: expression, fixed covariates, DIP; `<tissue>_genes.csv`: rhythmic labels |
| `peerCovariates/` | `<tissue>_block<b>.csv.gz`: donors and covariates (fixed + PEER) of every window and 24 h reference sample |
| `windowCounts/` | `<tissue>_block<b>.csv`: one row per permutation × window (nW, n24h, nShared, nGenes, nHits, nVanish, rhythmic split, nUntestable24h) |
| `windowHits/` | `<tissue>_observed.csv`: hit-level table for the observed DIP |
| `summary/` | `overall.csv`, `perTissue.csv`, null histograms, `missingPermutations.csv` if any task is incomplete |

## Test

`bash tests/run_toy_test.sh` simulates two toy tissues (phase, amplitude and mean eQTL, rhythmic-only and null genes), runs 01 to 04 with 39 DIP permutations and 1000 tensorQTL permutations per gene (on CPU), and checks that

- no permutation is missing and 8 fixed covariates are kept;
- observed vanishing hits exceed the 95th percentile of the null;
- mean eQTL rarely vanish (≤ 20 %).

Toy result: 87 hits (null mean 32), 46 vanishing (null mean 3.2, empirical p = 0.025, the minimum with 39 permutations); 25/28 phase and 16/17 amplitude hits vanish, 2/38 mean eQTL hits do.

The test uses PCA in place of PEER (`PEER_BACKEND = 'pca_test_only'`). Run it with `PYTHON=/path/to/python bash tests/run_toy_test.sh`.

## Not yet tested

- The `peer` R package on Bunya and its run time: 10 tissues × 1001 DIP sets × 8 windows × 2 subsets ≈ 160 000 fits. Time one `02` task before submitting the array.
- tensorQTL on a Bunya GPU and its run time: ≈ 80 000 `map_cis` calls with 10 000 permutations each. Check the GPU partition and `--gres` syntax in `03_tensorQTLscan.sh`, then time one task.
- Loading the real plink genotypes in tensorQTL (memory), and the real input paths.

## Requirements

R ≥ 4.4 with tidyverse, data.table, R.utils, jsonlite, peer. Python ≥ 3.10 with tensorQTL (PyTorch with CUDA), pandas, numpy, scipy.
