#' @export
extractData <- function(stack){
  
  time <- seq(1,stack@nbiter,1)/stack@iter
  data <- data.frame(
    time = rep(time, length(names(stack@layers))),
    sim = rep(names(stack@layers), each = length(time))
  )
  data$agb <- unlist(lapply(stack@layers, function(layer) layer@agb$Total))
  data$ba <- unlist(lapply(stack@layers, function(layer) layer@ba$ba$Total))
  data$n10 <- unlist(lapply(stack@layers, function(layer) layer@abundances$abu10$Total))
  data$n30 <- unlist(lapply(stack@layers, function(layer) layer@abundances$abu30$Total))
  
  species <- lapply(stack@layers, function(layer){
    abd <- table(row.names(layer@sp_par)[layer@final_pattern$sp_lab])
    abd <- abd / sum(abd)
    abd <- data.frame(row.names = names(abd), abd = as.vector(abd))
  })
  spdata <- data.frame(row.names = unique(unlist(lapply(species, row.names))))
  species <- lapply(species, function(abd){
    spdata$abd <- 0
    spdata[row.names(abd),] <- abd$abd
    return(spdata)
  })
  species <- do.call('cbind', species)
  names(species) <- names(stack@layers)
  
  return(list(
    time = data,
    species = species
  ))
}