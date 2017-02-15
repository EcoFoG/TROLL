#' Build TROLL
#'
#' Function to build TROLL code
#'
#' @return build the TROLL.out file in ./src
#'
#' @export
#'
#' @examples
#'
build <- function(){
  system('g++  src/main.cpp  -o src/TROLL.out')
}
