#' Computation of the integral of a cubic B-spline basis over (t0,x) for a given vector of knots
#'
#' @param t0 scalar giving lower value of the integration interval.
#' @param x vector giving the upper values of the integration interval.
#' @param knots vector of knots spanning the desired B-spline basis.
#'
#' @return A matrix of dimension \code{length(x)} by \code{(length(knots)+2)}.
#' @return Each integrated cubic B-spline is within a given column.
#' @export
#' @useDynLib cubicBsplines, .registration = TRUE
#' @examples IBsplines(t0=0,x=runif(20),knots=seq(0,1,length=11))
IBsplines = function(t0, x, knots){
  nx = length(x)
  nknots = length(knots)
  IB = matrix(0,nrow=nx,ncol=nknots+2)
  out = .Fortran("integrated_cubicBsplines_general",
                 t0= as.double(t0),
                 nx= as.integer(nx),
                 x= as.double(x),
                 nknots = as.integer(nknots),
                 knots = as.double(knots),
                 IB = as.double(c(IB)))
  return(matrix(out$IB,nrow=nx,ncol=nknots+2))
}
