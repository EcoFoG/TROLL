#' @export
updateHTML <- function(origin = './inst/doc',
                       destination = './docs'){
  docs <- list.files(path = origin, pattern = '.html')
  for(doc in docs){
    file.copy(from = file.path(origin, doc),
              to = file.path(destination, doc),
              overwrite = TRUE)
    cat(doc, 'updated\n')
  }
}