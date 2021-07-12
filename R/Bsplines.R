#' Computation of a cubic B-spline basis associated to a given vector of knots
#'
#' @param x vector of values where the B-spline basis must be evaluated.
#' @param knots vector of knots spanning the desired B-spline basis.
#'
#' @return A matrix of dimension \code{length(x)} by \code{(length(knots)+2)}.
#' @return Each column of the matrix corresponds to one cubic B-spline in the basis.
#' @export
#' @useDynLib cubicBsplines, .registration = TRUE
#' @examples Bsplines(x=runif(20),knots=seq(0,1,length=11))
Bsplines = function(x, knots){
  nx = length(x)
  nknots = length(knots)
  B = matrix(0,nrow=nx,ncol=nknots+2)
  out = .Fortran("cubicBsplines_general",
                 nx= as.integer(nx),
                 x= as.double(x),
                 nknots = as.integer(nknots),
                 knots = as.double(knots),
                 B = as.double(c(B)),
                 PACKAGE="cubicBsplines")
  return(matrix(out$B,nrow=nx,ncol=nknots+2))
}
