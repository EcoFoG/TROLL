#' @export
extractData <- function(stack,
                        sp = 'Total',
                        vars = c('agb', 'ba', 'n', 'n10', 'n30',
                                 'gpp', 'npp', 'Rday', 'Rnight'),
                        reduce = TRUE){
  
  time <- seq(1,stack@nbiter,1)/stack@iter
  data <- data.frame(time = rep(time, length(names(stack@layers))),
                     sim = rep(names(stack@layers), each = length(time)))
  
  if('agb' %in% vars)
    data$agb <- unlist(lapply(stack@layers, function(layer) layer@agb[[sp]]))
  if('ba' %in% vars)
    data$ba <- unlist(lapply(stack@layers, function(layer) layer@ba$ba[[sp]]))
  if('n' %in% vars)
    data$n <- unlist(lapply(stack@layers, function(layer) layer@abundances$abund[[sp]]))
  if('n10' %in% vars)
    data$n10 <- unlist(lapply(stack@layers, function(layer) layer@abundances$abu10[[sp]]))
  if('n30' %in% vars)
    data$n30 <- unlist(lapply(stack@layers, function(layer) layer@abundances$abu30[[sp]]))
  if('gpp' %in% vars)
    data$gpp <- unlist(lapply(stack@layers, function(layer) layer@gpp[[sp]]))
  if('npp' %in% vars)
    data$npp <- unlist(lapply(stack@layers, function(layer) layer@npp[[sp]]))
  if('Rday' %in% vars)
    data$Rday <- unlist(lapply(stack@layers, function(layer) layer@R$Rday[[sp]]))
  if('Rnight' %in% vars)
    data$Rnight <- unlist(lapply(stack@layers, function(layer) layer@R$Rnight[[sp]]))
  
  if(reduce)
    data <- data[data$time %in% round(data$time),]
  
  return(data)
}