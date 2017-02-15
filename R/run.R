#' Run TROLL
#'
#' Function to run TROLL code
#'
#' @param name char. name of the model
#' @param input char. input file path
#' @param path char. apth to save the ouptuts
#' @param overwrite logical. allow to overwrite existing outputs files
#'
#' @return model output files in the path folder
#'
#' @export
#'
#' @examples
#'
run <- function(name,
                input = './src/input.txt',
                path = './src/OUTPUT',
                overwrite = FALSE){

  if(name %in% list.dirs(path, full.names = FALSE)[-1]){
    if(!overwrite)
      stop('Outputs already exist, use overwrite = T.')
    path <- file.path(path, name)
    unlink(path, recursive = TRUE)
  } else {
    path <- file.path(path, name)
  }
  dir.create(path)

  input <- substring(input, 3)
  input <- paste0("'", input, "'")
  path <- paste0(path, '/', name)
  command <- paste0('./src/TROLL.out -i', input, ' -o', path)
  system(command)
}
