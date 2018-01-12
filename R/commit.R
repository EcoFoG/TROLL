#' @export
commit <- function(m,
                  HTML = TRUE,
                  check = TRUE,
                  push = FALSE){
  if(HTML)
    updateHTML()
  system('git add *')
  if(check)
    system('git status') ; invisible(readline(prompt="Press [enter] to continue"))
  command <- paste0("git commit -m '", m, "'")
  if(check)
    cat(command) ; invisible(readline(prompt="Press [enter] to continue"))
  system(command)
  if(push)
    system('git push origin master')
}