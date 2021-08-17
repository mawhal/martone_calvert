
#                                 Appendix 4, R function
#
# An R function to compute the box.cox.chord transformation.


#' Compute the box.cox.chord transformation on quantitative community composition 
#' data for any exponent. Usual exponents are larger than or equal to 0.
#'
#' Arguments --
#' @param mat : matrix or data.frame of quantitative non-negative community 
#'    composition data (frequencies, biomasses, energy measures, etc.)
#' @param bc.exp : Box-Cox exponent to the data before chord transformation. 
#'    Usual exponent values are {1, 0.5, 0.25, 0}, where 
#'    bc.exp=1: no transformation; 
#'    bc.exp=0.5: square-root transformation; 
#'    bc.exp=0.25: fourth-root (or double square-root) transformation; 
#'    bc.exp=0: log(y+1) transformation (default value). 
#'    Default value: bc.exp=0 (log(y+1) transformation).
#'
#' Value --
#' A Box-Cox+chord transformed matrix of the same size as the original data matrix.
#'
#' Author:: Pierre Legendre
#' License: GPL (>=2)

box.cox.chord <- 
	function(mat, 
             bc.exp=0) 
{ 
# Internal function
vec.norm <- function(vec)  sqrt(sum(vec^2))
#
chck <- apply(mat, 1, sum)
if(any(chck == 0)) stop("Rows",which(chck==0)," of the data matrix sum to 0")
#
# Apply the user-selected Box-Cox exponent (bc.exp) to the frequency data
if(bc.exp==0) {
	tmp <- log(mat+1) 
	} else { 
	tmp <- mat^bc.exp 
	}
row.norms <- apply(tmp, 1, vec.norm)
#
# Apply the chord transformation to matrix "tmp" before returning it
res <- sweep(tmp, 1, row.norms, "/")
}
