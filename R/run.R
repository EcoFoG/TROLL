#' Run TROLL
#'
#' Function to run TROLL code
#'
#' @param name char. name of the model
#' @param path char. working directory
#' @param app char. path to troll app (e.g. TROLL.out)
#' @param input char. input file
#' @param output char. name of the folder to save outputs
#' @param overwrite logical. allow to overwrite existing outputs files
#' @param verbose logical. allow output in console
#'
#' @return model output files in the path folder
#'
#' @export
#'
#' @examples
#'
run <- function(name,
                path = getOption("TROLL.path"),
                app = getOption("TROLL.app"),
                input = getOption("TROLL.init"),
                output = getOption("TROLL.output"),
                overwrite = TRUE,
                verbose = TRUE){

  if(name %in% list.dirs(file.path(path, output), full.names = FALSE)[-1]){
    if(!overwrite)
      stop('Outputs already exist, use overwrite = T.')
    path_o <- file.path(path, output, name)
    unlink(path_o, recursive = TRUE)
  } else {
    path_o <- file.path(path, output, name)
  }
  dir.create(path_o)

  app_c <- paste0('./', app)
  input_c <- paste0("'", input, "'")
  output_c <- file.path(paste0('./', output), name, name)
  command <- paste0('cd ', path, ' ; ',
                    app_c,
                    ' -i', input_c, 
                    ' -o', output_c)
  if(verbose)
    cat(command, '\n')
  system(command)
}
