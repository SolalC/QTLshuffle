# =============================================================================
# 03  tensorQTL scan of every window and 24 h reference re-test (one tissue x
#     permutation block per SLURM array task; same task numbering as 02)
#
# For each permutation and window with covariates from 02:
#   1. cis-eQTL scan of the window donors with tensorQTL map_cis (GTEx
#      protocol: +/- CIS_WINDOW, TQ_NPERM genotype permutations per gene,
#      beta-approximated p-values, Storey q-value with lambda QVALUE_LAMBDA);
#      variants need MAF >= MAF_THRESHOLD across the tissue's donors and a
#      minor allele count >= MIN_MAC in the window;
#   2. eGenes: q < GENE_FDR;
#   3. the lead variant of each eGene re-tested in the 24 h reference sample,
#      same model (covariates of the reference sample);
#   4. vanishing hit: reference p > VANISH_P and window vs reference
#      difference z test p < HET_P (conservative: the two samples share
#      about 1/8 of their donors).
#
# All paths and constants come from Results/main/config.json, written by
# 01_prepareTissues.R from 00_config.R.
#
# Run:     sbatch 03_tensorQTLscan.sh       (after 02)
# Usage:   python 03_tensorQTLscan.py <Results/main/config.json> <task index>
# Outputs: Results/main/windowCounts/<tissue>_block<b>.csv  (one row per perm x window)
#          Results/main/windowHits/<tissue>_observed.csv    (block 1 only: hits, perm 0)
# =============================================================================

import json
import sys
import time

import numpy as np
import pandas as pd
import tensorqtl
import torch
from scipy import stats
from tensorqtl import cis, genotypeio, post

cfg = json.load(open(sys.argv[1]))
task = int(sys.argv[2])
device = 'cuda' if torch.cuda.is_available() else 'cpu'
print(f'tensorQTL device: {device}', flush=True)

tissues = pd.read_csv(cfg['TISSUE_DIR'] + 'selectedTissues.csv')
nBlock = int(np.ceil((cfg['N_PERM'] + 1) / cfg['PERM_BLOCK']))
tissueIdx = (task - 1) // nBlock + 1
block = (task - 1) % nBlock + 1
if tissueIdx > len(tissues):
    sys.exit(f'Task {task} is beyond the last tissue')
tissueName = tissues['tissue'][tissueIdx - 1]
print(f'Tissue: {tissueName}; block {block}', flush=True)

# --- Inputs --------------------------------------------------------------------

covAll = pd.read_csv(cfg['PEER_DIR'] + f'{tissueName}_block{block}.csv.gz')
genes = pd.read_csv(cfg['TISSUE_DIR'] + f'{tissueName}_genes.csv').set_index('gene_id')
donors = sorted(covAll['donor'].unique())

exprFile = cfg['GTEX_EXPR_DIR'] + f'{tissueName}.v10.normalized_expression.bed.gz'
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(exprFile)
phenotype_df = phenotype_df[donors]

# Genotypes: WGS sample IDs (GTEX-XXXX-...) renamed to subject IDs (GTEX-XXXX),
# restricted to chr1-22, the tissue's donors and MAF >= MAF_THRESHOLD among them.
genotype_df, variant_df = genotypeio.load_genotypes(cfg['GENOTYPE_PLINK'], dosages=True)
genotype_df.columns = ['-'.join(c.split('-')[:2]) for c in genotype_df.columns]
genotype_df = genotype_df[donors]
# Missing genotypes are -9 (int8), as tensorQTL expects; it mean-imputes them.
autosomes = [f'chr{i}' for i in range(1, 23)]
af = genotype_df.where(genotype_df >= 0).mean(axis=1, skipna=True) / 2
keep = variant_df['chrom'].isin(autosomes).values & (np.minimum(af, 1 - af) >= cfg['MAF_THRESHOLD']).values
genotype_df, variant_df = genotype_df[keep], variant_df[keep]
phenotype_df = phenotype_df[phenotype_pos_df['chr'].isin(autosomes)]
phenotype_pos_df = phenotype_pos_df.loc[phenotype_df.index]
print(f'{genotype_df.shape[0]} variants, {phenotype_df.shape[0]} genes, {len(donors)} donors', flush=True)

# --- Re-test of lead variants in another donor set -----------------------------

def leadVariantTest(hits, rows, cov):
    """Lead variant of each hit regressed on expression in donors `rows`, adjusted
    for `cov` (donors x covariates) and an intercept; dof as in tensorQTL."""
    C = np.column_stack([np.ones(len(rows)), cov.loc[rows].values])
    Q, _ = np.linalg.qr(C)
    dof = len(rows) - 2 - cov.shape[1]
    out = []
    for gene, variant in zip(hits.index, hits['variant_id']):
        g = genotype_df.loc[variant, rows].values.astype(float)
        g[g < 0] = np.nan
        g[np.isnan(g)] = np.nanmean(g)
        ac = g.sum()
        if min(ac, 2 * len(g) - ac) < cfg['MIN_MAC']:
            out.append((np.nan, np.nan, np.nan))
            continue
        y = phenotype_df.loc[gene, rows].values.astype(float)
        gr = g - Q @ (Q.T @ g)
        yr = y - Q @ (Q.T @ y)
        sxx = np.sum(gr ** 2)
        b = np.sum(gr * yr) / sxx
        se = np.sqrt(max(np.sum(yr ** 2) - b ** 2 * sxx, 0) / dof / sxx)
        out.append((b, se, 2 * stats.t.sf(abs(b / se), dof)))
    return pd.DataFrame(out, index=hits.index, columns=['beta_24h', 'se_24h', 'p_24h'])

# --- One window ----------------------------------------------------------------

def covariates(perm, window, setName):
    d = covAll[(covAll['perm'] == perm) & (covAll['window'] == window) & (covAll['set'] == setName)]
    c = d.drop(columns=['perm', 'window', 'set']).set_index('donor')
    return c.dropna(axis=1, how='all')

def scanWindow(perm, window):
    covW = covariates(perm, window, 'window')
    covR = covariates(perm, window, 'reference')
    rowsW, rowsR = list(covW.index), list(covR.index)
    res = cis.map_cis(genotype_df[rowsW], variant_df, phenotype_df[rowsW], phenotype_pos_df,
                      covariates_df=covW, nperm=int(cfg['TQ_NPERM']), window=int(cfg['CIS_WINDOW']),
                      maf_threshold=cfg['MIN_MAC'] / (2 * len(rowsW)),
                      seed=int(cfg['SEED']) + 100 * perm + window, verbose=False)
    post.calculate_qvalues(res, fdr=cfg['GENE_FDR'], qvalue_lambda=cfg['QVALUE_LAMBDA'])
    hits = res[res['qval'] < cfg['GENE_FDR']].copy()
    hits['rhythmic'] = genes.reindex(hits.index)['rhythmic'].fillna(False).astype(bool).values
    hits = hits.join(leadVariantTest(hits, rowsR, covR)) if len(hits) else \
        hits.assign(beta_24h=np.nan, se_24h=np.nan, p_24h=np.nan)
    hits['z'] = (hits['slope'] - hits['beta_24h']) / np.sqrt(hits['slope_se'] ** 2 + hits['se_24h'] ** 2)
    hits['pHet'] = 2 * stats.norm.sf(np.abs(hits['z']))
    hits['vanish'] = hits['p_24h'].notna() & (hits['p_24h'] > cfg['VANISH_P']) & (hits['pHet'] < cfg['HET_P'])
    counts = dict(nW=len(rowsW), n24h=len(rowsR), nShared=len(set(rowsW) & set(rowsR)),
                  nGenes=len(res), nHits=len(hits), nVanish=int(hits['vanish'].sum()),
                  nHitsRhythmic=int(hits['rhythmic'].sum()),
                  nVanishRhythmic=int((hits['vanish'] & hits['rhythmic']).sum()),
                  nUntestable24h=int(hits['p_24h'].isna().sum()))
    return counts, hits

# --- All permutations and windows of the block ---------------------------------

perms = range((block - 1) * cfg['PERM_BLOCK'], min(block * cfg['PERM_BLOCK'] - 1, cfg['N_PERM']) + 1)
countRows, hitTables = [], []
for perm in perms:
    t0 = time.time()
    for window in range(1, cfg['N_WINDOWS'] + 1):
        sub = covAll[(covAll['perm'] == perm) & (covAll['window'] == window)]
        if sub.empty:   # fewer than MIN_WINDOW_N donors in the window (see 02)
            countRows.append(dict(tissue=tissueName, perm=perm, window=window))
            continue
        counts, hits = scanWindow(perm, window)
        countRows.append(dict(tissue=tissueName, perm=perm, window=window, **counts))
        if perm == 0:
            hitTables.append(hits.assign(tissue=tissueName, window=window))
    print(f'perm {perm}: {(time.time() - t0) / 60:.2f} min', flush=True)

pd.DataFrame(countRows).to_csv(cfg['COUNT_DIR'] + f'{tissueName}_block{block}.csv', index=False)
if 0 in perms:
    hitCols = ['tissue', 'window', 'variant_id', 'num_var', 'slope', 'slope_se', 'pval_nominal',
               'pval_beta', 'qval', 'rhythmic', 'beta_24h', 'se_24h', 'p_24h', 'z', 'pHet', 'vanish']
    obs = pd.concat(hitTables) if hitTables else pd.DataFrame(columns=hitCols)
    obs.rename_axis('gene').reset_index()[['gene'] + hitCols] \
       .to_csv(cfg['HIT_DIR'] + f'{tissueName}_observed.csv', index=False)
