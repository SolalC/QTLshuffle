# PEER fit run under PEER's own R (3.4, conda environment; see the rhyh2 repo,
# Code/peer/install_peer.sh). Base R only. Called by runPEER() in 00_functions.R (copy of prototype/peer_fit.R).
# Settings as GTEx: default priors, mean added, PEER_ITER iterations; returns the
# last k columns of PEER_getX (the hidden factors).
# Usage: <PEER_RSCRIPT> peer_fit.R <input.rds> <output.rds> <k> <iterations> <seed>
suppressMessages(library(peer))
a <- commandArgs(trailingOnly = TRUE)
Y <- readRDS(a[1])                     # donors x genes
k <- as.integer(a[3])
set.seed(as.integer(a[5]))
model <- PEER()
invisible(PEER_setPhenoMean(model, as.matrix(Y)))   # setters print NULL otherwise
invisible(PEER_setAdd_mean(model, TRUE))
invisible(PEER_setNk(model, k))
invisible(PEER_setNmax_iterations(model, as.integer(a[4])))
invisible(PEER_update(model))
X <- PEER_getX(model)
X <- X[, (ncol(X) - k + 1):ncol(X), drop = FALSE]
saveRDS(X, a[2], version = 2)
