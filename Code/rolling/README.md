# QTLshuffle: rolling-window eQTL (R only)

Are cis-eQTL found in a narrow window of inferred circadian phase (DIP) larger there than in random samples of the same size taken over 24 h? One tissue (`TISSUE` = Muscle_Skeletal), all expressed genes, 3 h windows first (6 h and 12 h later, by `WINDOW_WIDTH_H`).

## Design

| Step | Choice |
|---|---|
| Windows | 24 centres, 0 h to 23 h, 1 h apart; width `WINDOW_WIDTH_H` = 3 h; a donor is in the window if its DIP is in [centre − 1.5 h, centre + 1.5 h), wrapping around midnight (the 0 h window covers 22:30 to 01:30). Centres with < 60 donors are skipped |
| Covariates | sex, platform, pcr, PC1–5 (GTEx v10 covariate file) and PEER factors recomputed on the window donors (all genes). Number of PEER factors as GTEx: 15 if n < 150, 30 if n < 250, 45 if n < 350, 60 otherwise |
| Genotypes | GTEx v9 WGS VCF (953 donors, all variants), converted once with plink2: tissue donors, FILTER = PASS, biallelic, chr1–22, MAF ≥ 1 % among the tissue's donors. Dosage counts the REF allele |
| Scan | every gene, variants within ±1 Mb of the TSS, minor allele count ≥ 10 in the window; `y ~ g + covariates` |
| eGenes (GTEx) | 10 000 expression permutations per gene, beta-approximated p (tensorQTL's effective-dof fit), Storey q-value with λ = 0.85 across all genes of the window, q ≤ 0.05 |
| Significant pairs (GTEx) | for each eGene, variants with nominal p ≤ qbeta(p_threshold, shape1, shape2), where p_threshold is midway between the largest pBeta with q ≤ 0.05 and the smallest with q > 0.05 (tensorQTL `post.calculate_qvalues`) |
| Independent signals | greedy clumping of each eGene's significant variants (as plink `--clump`): the most significant remaining variant is an index variant, and variants with r² > 0.01 with it join its clump. No distance limit within the cis window. LD from all the tissue's donors (see below) |
| 24 h reference | per centre, 1000 samples of n_W donors drawn from all the tissue's donors, whatever their DIP (window donors can be drawn); PEER recomputed in each with the same rule; every independent signal of the centre re-tested |

### Choices to check

- **Permutations.** Freedman–Lane in place of tensorQTL's plain residual permutation, which is anti-conservative at window sizes (about 24 covariates for about 90 donors; see `../prototype/README.md`, 22 % against 4.8 % of null genes at p < 0.05).
- **LD for clumping.** With about 90 window donors, r² between unlinked variants averages about 1/90, above the 0.01 cut-off. Window-donor LD would therefore merge independent signals at random. LD is computed on all the tissue's donors instead.
- **Stored pairs.** 03 keeps the pairs with nominal p < 0.01, and 04 warns if an eGene's nominal threshold is above that.

## Run order (Bunya)

```bash
cd /scratch/user/uqschauq/QTLshuffle/Code && mkdir -p logs && cd rolling   # logs go to Code/logs
j1=$(sbatch --parsable 01_prepareTissue.sh)          # tissue data; PEER for the 24 windows
j2a=$(sbatch --parsable 02a_vcfToPlink.sh)          # full WGS VCF to plink, one pass (set PLINK2 in the script)
j2=$(sbatch --parsable --dependency=afterok:$j2a 02b_prepareGenotypes.sh)  # array 1-22: genotypes per chromosome
sbatch --array 1 --dependency=afterok:$j1:$j2 03_windowScan.sh   # time one task first (0 h, chr1)
j3=$(sbatch --parsable --dependency=afterok:$j1:$j2 03_windowScan.sh)  # array 1-528: 24 centres x 22 chromosomes
j4=$(sbatch --parsable --dependency=afterok:$j3 04_windowCall.sh)      # eGenes, significant pairs, clumping
sbatch --array 1 --dependency=afterok:$j4 05_referenceRetest.sh        # time one task first (24 PEER fits)
j5=$(sbatch --parsable --dependency=afterok:$j4 05_referenceRetest.sh) # array 1-1000: one task per reference sample
sbatch --dependency=afterok:$j5 --wrap "module load r/4.4.2-heavy; Rscript 06_collect.R" -A a_imb_cpdg -p general -t 2:00:00 --mem 60G
```

Cost: one gene scan (90 donors, 7000 variants, 10 000 permutations) took 7 s on one core here, so about 40 core-hours per window and 1000 core-hours for 24 windows if muscle has 20 000 genes. PEER time is not known: 24 window fits plus 24 000 reference fits (1000 × 24). Time one 05 task (24 fits, 8 at a time) before submitting the array.

## Genotype files (`GENO_DIR` = /QRISdata/Q5776/QTLshuffle/)

All individual-level genotypes written by the pipeline go to `GENO_DIR`, not to the results folder:

- `plink/GTEx_v9_WGS_<tissue>.{bed,bim,fam,keep,log}`: plink file of the tissue's donors (02a);
- `<tissue>/genotypes_chr<c>.rds`: per-chromosome genotypes, raw codes (02b);
- `<tissue>/w<W>h_indexGenotypes.rds`: genotypes of the index variants (04, read by 05).

## Outputs (`Results/rolling/`)

Everything that is expensive is written to disk and re-used.

| Path | Content |
|---|---|
| `data/tissue.rds` | donors, DIP, expression (all genes), fixed covariates |

| `w3h/windows.csv` | per centre: number of donors, PEER factors, used |
| `w3h/windowCovariates.rds` | per centre: window donors and covariates (fixed + PEER), donor IDs as row names. 01 re-uses it if the window donors are unchanged |
| `w3h/scan/c<centre>_chr<c>_{genes,pairs}.csv.gz` | gene-level results (lead variant, pBeta, beta shapes, effective dof) and pairs with p < 0.01 |
| `w3h/calls/genes.csv.gz` | all genes with q-value and nominal threshold, per centre |
| `w3h/calls/significantPairs.csv.gz` | significant pairs with their index variant |
| `w3h/calls/independent.csv.gz` | index variants (one row per independent signal and centre) |
| `w3h/referenceCovariates/ref<r>.rds` | per reference sample and centre: donors and covariates (fixed + PEER). 05 re-uses them on a rerun, for example after changing a threshold in 04 |
| `w3h/reference/ref<r>.csv.gz` | re-test of every signal of every centre in reference sample r (β, SE, p, nShared) |
| `w3h/results/rollingWindow.csv.gz` | **main table**: one row per signal × centre, with the window β, SE and p; nRef, refMean, refSD; zRef = (β − refMean)/refSD; two-sided pEmp; one-sided pEmpLarger (reference β at least as large in the window's direction); then `beta_ref_1` to `beta_ref_1000` |
| `w3h/results/referenceLong.csv.gz` | β, SE and p of every reference association |

## Test

`bash tests/run_toy_test.sh` uses the toy data from `../tests/make_toy_data.R`: 560 donors, 300 genes, 2 chromosomes. It runs with 1000 permutations per gene, 40 reference samples and PCA in place of PEER, and took 4 min on 4 cores. 21 of 24 centres had ≥ 60 donors.

| Toy gene type | Signals (centre × gene) | pEmp < 0.05 | pEmpLarger < 0.05 |
|---|---|---|---|
| phase QTL | 81 | 100 % | 100 % |
| amplitude QTL | 52 | 100 % | 100 % |
| mean eQTL | 74 | 12 % | 20 % |
| null gene | 13 | 100 % | 100 % |

The null-gene row is the winner's curse: false-positive eGenes are selected for a large window β, so they always look larger than over 24 h. The 20 % of mean eQTL are the same effect on true eQTL. This is the comparison deferred until after the first run.
