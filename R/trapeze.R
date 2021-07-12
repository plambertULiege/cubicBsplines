#' Trapeze integration from a vector of function values evaluated at quadrature points
#'
#' @param x grid of values for the quadrature (vector).
#' @param fx values of the function on the grid (vector).
#'
#' @return vector with a numerical approximation of \eqn{\int_{min(x)}^{max(x)} f(t) dt} on the grid using the trapeze method.
#' @export
#' @examples
#' x = seq(-4,2,length=100) ; fx = dnorm(x) ; res = trapeze(x,fx)
#' cbind(true=pnorm(x),trapeze=res)
#'
trapeze = function(x,fx){
  dx = diff(x)
  nx = length(x)
  temp = .5 * fx[-nx] * dx
  temp = temp + .5 * fx[-1] * dx
  return(c(0,cumsum(temp)))
}
