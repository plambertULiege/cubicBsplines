#' Computation of the 2nd derivative of a cubic B-spline basis associated to a given vector of knots
#'
#' @param x vector of values where the 2nd derivative of the B-spline basis must be evaluated.
#' @param knots vector of knots spanning the desired B-spline basis.
#'
#' @return A matrix of dimension \code{length(x)} by \code{(length(knots)+2)}.
#' @return Each column of the matrix corresponds to (the 2nd derivative of) one cubic B-spline in the basis.
#' @export
#' @useDynLib cubicBsplines, .registration = TRUE
#' @examples D2Bsplines(x=runif(20),knots=seq(0,1,length=11))
D2Bsplines = function(x, knots){
  nx = length(x)
  nknots = length(knots)
  d2B = matrix(0,nrow=nx,ncol=nknots+2)
  out = .Fortran("D2_cubicBsplines_general",
                 nx= as.integer(nx),
                 x= as.double(x),
                 nknots = as.integer(nknots),
                 knots = as.double(knots),
                 d2B = as.double(c(d2B)))
  return(matrix(out$d2B,nrow=nx,ncol=nknots+2))
}
