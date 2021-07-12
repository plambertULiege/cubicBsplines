#' Computation of the 1st derivative of a cubic B-spline basis associated to a given vector of knots
#'
#' @param x vector of values where the 1st derivative of the B-spline basis must be evaluated.
#' @param knots vector of knots spanning the desired B-spline basis.
#'
#' @return A matrix of dimension \code{length(x)} by \code{(length(knots)+2)}.
#' @return Each column corresponds to (the 1st derivative of) one cubic B-spline in the basis.
#' @export
#' @useDynLib cubicBsplines, .registration = TRUE
#' @examples D1Bsplines(x=runif(20),knots=seq(0,1,length=11))
D1Bsplines = function(x, knots){
  nx = length(x)
  nknots = length(knots)
  dB = matrix(0,nrow=nx,ncol=nknots+2)
  out = .Fortran("D1_cubicBsplines_general",
                 nx= as.integer(nx),
                 x= as.double(x),
                 nknots = as.integer(nknots),
                 knots = as.double(knots),
                 dB = as.double(c(dB)))
  return(matrix(out$dB,nrow=nx,ncol=nknots+2))
}
