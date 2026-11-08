# QTLshuffle: rolling-window eQTL (R only)

Are cis-eQTL found in a narrow window of inferred circadian phase (DIP) larger there than in random samples of the same size taken over 24 h? One tissue (`TISSUE` = Muscle_Skeletal), all expressed genes, 3 h windows first (6 h and 12 h later, by `WINDOW_WIDTH_H`).

## Design

| Step | Choice |
|---|---|
| Windows | 24 centres, 0 h to 23 h, 1 h apart; width `WINDOW_WIDTH_H` = 3 h; a donor is in the window if its DIP is in [centre − 1.5 h, centre + 1.5 h), wrapping around midnight (the 0 h window covers 22:30 to 01:30). Centres with < 60 donors are skipped |
| Covariates | sex, platform, pcr, PC1–5 (GTEx v10 covariate file) and PEER factors recomputed on the window donors (all genes). Number of PEER factors as GTEx: 15 if n < 150, 30 if n < 250, 45 if n < 350, 60 otherwise |
| Genotypes | GTEx v9 WGS VCF (953 donors, all variants), converted once with plink2: tissue donors, FILTER = PASS, biallelic, chr1–22, MAF ≥ 1 % among the tissue's donors. Dosage counts the ALT allele (plink A1) |
| Scan | MatrixEQTL (Shabalin 2012, Bioinformatics 28:1353), linear model `y ~ g + covariates`, every gene, variants within ±1 Mb of the TSS with minor allele count ≥ 10 in the window; missing genotypes mean-imputed. Nominal p-values; **no permutations** |
| Significant pairs and eGenes | Benjamini–Hochberg q < 0.05 over all cis tests of the window (pooled across chromosomes); an eGene has ≥ 1 significant pair. This is a pair-level FDR, not the GTEx gene-level permutation threshold |
| Independent signals | `plink2 --clump` (Purcell et al. 2007; Chang et al. 2015) of each eGene's significant variants: r² > 0.01 within 2000 kb (the whole cis window), p1 = p2 = 1. LD from all the tissue's donors (see below) |
| 24 h reference | per centre, 1000 samples of n_W donors drawn from all the tissue's donors, whatever their DIP (window donors can be drawn); PEER recomputed in each with the same rule; every independent signal of the centre re-tested with `lm(y ~ g + covariates)` |
| Selection-bias control | the window β is the β that selected the signal (threshold, then lead variant within the clump), so it is inflated relative to the unselected reference β even without a time effect (winner's curse). References 1–100 of each centre are therefore treated as windows (same scan, BH and clumping: 07); their signals are tested against the other 99 references and the window signals against all 100, with the same MatrixEQTL model (08). A window signal's calibrated p is the share of these pseudo-window signals, in the same window \|z\| bin, with \|zRef\| at least as large; BH within the window gives qCal; time dependent: qCal < 0.05. Window level: % of signals with pEmp < 0.05 against the same % in each pseudo-window |

### Choices to check

- **No permutations.** The threshold is nominal (BH over all cis pairs), so eGene counts are not directly comparable with GTEx's gene-level permutation threshold; state this in any write-up.
- **LD for clumping.** With about 90 window donors, r² between unlinked variants averages about 1/90, above the 0.01 cut-off. Window-donor LD would therefore merge independent signals at random. LD is computed on all the tissue's donors instead.
- **Stored pairs.** 03 keeps pairs with p < 0.01; 04 warns if significant pairs come close to that.
- **Selection bias.** Uncorrected, 66 % of the Muscle signals had pEmp < 0.05 (62 % after BH), falling from 80 % at window |z| ≤ 5 to 10 % at |z| > 15: the pattern of selection on the window β, not of a time effect. Only the calibrated calls (08) correct for it. On synthetic data (`no time effect` vs `+2 SD in the window`, lead variant chosen among 20 correlated proxies), pEmp < 0.05 called 21 % of time-invariant signals; calibrated qCal < 0.05 called 0.2 % of them and 21 % of time-varying ones.
- **Pseudo-window resolution.** Pseudo-window signals are compared with 99 references, window signals with 100 (pEmp ≥ 1/101); zRef, not pEmp, is calibrated, so that the ordering within the tail is kept.

## Run order (Bunya)

```bash
cd /scratch/user/uqschauq/QTLshuffle/Code && mkdir -p logs && cd rolling   # logs go to Code/logs
j1=$(sbatch --parsable 01_prepareTissue.sh)          # tissue data; PEER for the 24 windows
j2=$(sbatch --parsable 02a_vcfToPlink.sh)           # full WGS VCF to plink (skipped if done), split by chromosome
sbatch --array 1 --dependency=afterok:$j1:$j2 03_windowScan.sh   # time one task first (0 h, chr1)
j3=$(sbatch --parsable --dependency=afterok:$j1:$j2 03_windowScan.sh)  # array 1-528: 24 centres x 22 chromosomes
j4=$(sbatch --parsable --dependency=afterok:$j3 04_windowCall.sh)      # FDR, eGenes, plink2 clumping
sbatch --array 1 --dependency=afterok:$j4 05_referenceRetest.sh        # time one task first (24 PEER fits)
j5=$(sbatch --parsable --dependency=afterok:$j4 05_referenceRetest.sh) # array 1-1000: one task per reference sample
j6=$(sbatch --parsable --dependency=afterok:$j5 06_collect.sh)
j7=$(sbatch --parsable --dependency=afterok:$j5 07_nullScan.sh)       # array 1-100: reference k scanned and clumped as a window
j8=$(sbatch --parsable --dependency=afterok:$j7 08_nullCalibrate.sh)  # array 1-24: re-tests and calibration per centre
sbatch --dependency=afterok:$j6:$j8 09_plotWindows.sh   # figures; args: [p=0.05] [test=pEmp|pEmpLarger] [independent=c,gene,var] [dependent=c,gene,var]
```

Cost: 03 is one MatrixEQTL run per centre × chromosome on one core (time one task first). 01 took 2 h 57 min for the 24 window PEER fits in parallel, so a PEER fit takes about 3 h; 05 is 1000 tasks of 24 fits in parallel (24 cores, ~3 h each), about 72 000 core-hours in all, capped at 1536 cores at a time by the `normal` QOS. 07 and 08 fit no PEER (05 saved the covariates): 07 is 100 tasks of 24 × 22 MatrixEQTL scans plus clumping, 08 is 24 tasks of 100 × 22 MatrixEQTL runs on the signal pairs only; time one task of each first.

## Genotype files (`GENO_DIR` = /QRISdata/Q5776/QTLshuffle/)

The only genotype files the pipeline writes are the plink files of the tissue's donors (02a), in `GENO_DIR`, not in the results folder: `plink/GTEx_v9_WGS_<tissue>.{bed,bim,fam,keep,log}` and its split by chromosome, `plink/<tissue>_chr/` (for `plink2 --clump`). Steps 03 to 05 read it directly with BEDMatrix (memory-mapped); variant and sample tables are read with genio. Dosage counts the A1 (ALT) allele.

## Outputs (`Results/rolling/`)

Everything that is expensive is written to disk and re-used.

| Path | Content |
|---|---|
| `data/tissue.rds` | donors, DIP, expression (all genes), fixed covariates |

| `w3h/windows.csv` | per centre: number of donors, PEER factors, used |
| `w3h/windowCovariates.rds` | per centre: window donors and covariates (fixed + PEER), donor IDs as row names. 01 re-uses it if the window donors are unchanged |
| `w3h/scan/c<centre>_chr<c>_{genes,pairs}.csv.gz` | MatrixEQTL output: per gene, number of cis variants and smallest p; pairs with p < 0.01 (β, SE, t, p); number of cis tests |
| `w3h/calls/genes.csv.gz` | per centre: every gene tested, its number of cis variants, smallest p and number of significant pairs |
| `w3h/calls/significantPairs.csv.gz` | significant pairs with their index variant |
| `w3h/calls/independent.csv.gz` | index variants (one row per independent signal and centre) |
| `w3h/referenceCovariates/ref<r>.rds` | per reference sample and centre: donors and covariates (fixed + PEER). 05 re-uses them on a rerun, for example after changing a threshold in 04 |
| `w3h/reference/ref<r>.csv.gz` | re-test of every signal of every centre in reference sample r (β, SE, p, nShared) |
| `w3h/results/rollingWindow.csv.gz` | **main table**: one row per signal × centre, with the window β, SE and p; nRef, refMean, refSD; zRef = (β − refMean)/refSD; two-sided pEmp; one-sided pEmpLarger (reference β at least as large in the window's direction); then `beta_ref_1` to `beta_ref_1000` |
| `w3h/results/referenceLong.csv.gz` | β, SE and p of every reference association |
| `w3h/null/scan/ref<k>.csv.gz`, `ref<k>_summary.csv` | 07: pseudo-window signals of reference k at every centre (β, SE, p, q, clump), and per centre the cis tests, significant pairs, eGenes and signals |
| `w3h/null/calibration/c<centre>_betas.rds` | 08: β of every window and pseudo-window signal in references 1–100 (re-used on a rerun) |
| `w3h/null/calibration/c<centre>_window.csv.gz` | 08: window signals with zW, \|z\| bin, zRefK and pEmpK (against the 100 references), pCal, qCal, timeDependent |
| `w3h/null/calibration/c<centre>_{null.csv.gz,nullByRef.csv,summary.csv}` | 08: pseudo-window signals (zRef, pEmp against the other 99), % pEmp < 0.05 per pseudo-window, and per centre the window % against the null % (mean, SD, range, pWindow) |
| `w3h/figures/rollingWindow_{nominal,BH,calibrated}.{pdf,png}` | 09: (a) time-independent and (c) time-dependent example, histogram of the 1000 reference betas with the window beta; (b) number of independent and of time-dependent signals per centre; (d) stacked bars to 100%: time dependent vs time independent. Time dependent: pEmp < 0.05 (nominal), BH-adjusted pEmp < 0.05 within each window (BH), or qCal < 0.05 (calibrated, 08) |
| `w3h/figures/nullCalibration.{pdf,png,csv}` | 09: per centre, % of window signals with pEmp < 0.05 (red) against the same % in its 100 pseudo-windows (box) |
| `w3h/figures/associationsPerWindow_{...}.csv`, `exampleSignals_{...}.csv` | 09: the counts behind b and d, and the two example rows (with gene symbol and rsID) |

## Requirements

R ≥ 4.4 (`module load r/4.4.2-heavy`) with tidyverse, data.table, MatrixEQTL, BEDMatrix and genio; plink2 (`PLINK2` in `00_config.R`) for 02a and 04; PEER in its own R 3.4 environment (`PEER_RSCRIPT`).

## Test

`bash tests/run_toy_test.sh` runs 01 to 08 on the toy data from `../tests/make_toy_data.R` (560 donors, 300 genes, 2 chromosomes; no real data), with 40 reference samples (10 of them as pseudo-windows) and PCA in place of PEER, then checks the output (`tests/check_toy_results.R`). It needs MatrixEQTL, BEDMatrix, genio and plink2, so run it on Bunya:

```bash
cd /scratch/user/uqschauq/QTLshuffle/Code/rolling
sbatch tests/run_toy_test_slurm.sh      # logs: ../logs/rollToyTest.{out,err}
```

The checks: 24 centres with the 0 h window wrapping around midnight; no missing reference samples; phase and amplitude QTL fall outside the 24 h distribution more often than mean eQTL.
