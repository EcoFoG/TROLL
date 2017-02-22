#' Build TROLL
#'
#' Function to build TROLL code
#'
#' @param src char. path to src file
#' @param app char. path to troll app (e.g. TROLL.out)
#' @param path char. working directory
#' 
#' @return build the TROLL app
#'
#' @export
#'
#' @examples
#'
build <- function(
  src = getOption("TROLL.src"),
  app = getOption("TROLL.app"),
  path = getOption("TROLL.path"),
  verbose = TRUE
){
  output <- file.path(path, app)
  command <- paste('g++', src, '-o', output)
  if(verbose)
    cat(command, '\n')
  system(command)
}
