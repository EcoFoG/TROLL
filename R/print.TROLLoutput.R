#' @import methods
NULL

#' Function to print TROLL outputs.
#'
#' @param TROLLoutput
#'
#' @return Print in console
#'
#' @export
#'
#' @examples
#'
setMethod('print', 'TROLLoutput', function(x, ...) {

  cat('Object of class :', class(x)[1],'\n\n')

  cat('Name :', x@name, '\n\n')

  cat('2D discrete network: horizontal step = ', x@info$step, 'm, one tree per 1 m^2 \n')
  cat('Number of sites      : ', x@info$SitesNb[1], 'x', x@info$SitesNb[2], '\n')
  cat('Number of iterations : ', x@info$IterationsNb, '\n')
  cat('Duration of timestep : ', x@info$timestep, 'years \n')
  cat('Number of Species    : ', x@info$SpeciesNb, '\n\n')

  cat('Average computation time : ', x@info$ComputationTime, 'seconds \n\n')
})
