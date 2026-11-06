# One-off check on Bunya: BEDMatrix genotypes (new 03-05) against the
# per-chromosome files written by the former 02b reader, on chr1.
# Prints counts only (no genotypes). Run from Code/rolling, before deleting
# GENO_DIR/<tissue>/genotypes_chr*.rds:
#   Rscript tests/check_bedmatrix.R
.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
Sys.setenv(QTLSHUFFLE_CODE_DIR = normalizePath(file.path(dirname(normalizePath(.f)), '..')))
source(file.path(Sys.getenv('QTLSHUFFLE_CODE_DIR'), '00_config.R'))
old <- readRDS(paste0(GENO_DIR, TISSUE, '/genotypes_chr1.rds'))
set.seed(1)
j <- sort(sample(ncol(old$G), min(5000, ncol(old$G))))
geno <- openGenotypes(rownames(old$G))
G <- readGenotypes(geno, seq_along(geno$rows), match(colnames(old$G)[j], geno$bim$id))
O <- matrix(as.integer(old$G[, j]), nrow(old$G)); O[O == 3] <- NA   # old: copies of A2
message('Variants compared: ', length(j), '; donors: ', nrow(G))
message('Missing in the same places: ', identical(is.na(G), is.na(O)))
message('BEDMatrix = old (A2 count): ', isTRUE(all(G == O, na.rm = TRUE)))
message('BEDMatrix = 2 - old (A1 count): ', isTRUE(all(G == 2 - O, na.rm = TRUE)))
