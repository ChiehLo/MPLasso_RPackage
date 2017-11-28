#####################################################################
# Log-ratio transform
#
# Adopted from Zachary D. Kurtz et. al., Sparse and Compositionally Robust Inference of Microbial Ecological Networks, Plos One, 2015, http://dx.doi.org/10.1371/journal.pcbi.1004226
#####################################################################

#' Centered log-ratio functions
#' @export
clr <- function(x, ...) {
    UseMethod('clr', x)
}

#' @method clr default
#' @export
clr.default <- function(x.f, base=exp(1), tol=.Machine$double.eps) {
    nzero <- (x.f >= tol)
    LOG <- log(ifelse(nzero, x.f, 1), base)
    ifelse(nzero, LOG - mean(LOG)/mean(nzero), 0.0)
}

#' @method clr matrix
#' @export
clr.matrix <- function(x.f, mar=2, ...) {
    apply(x.f, mar, clr, ...)
}
