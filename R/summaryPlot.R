#' @export
summaryPlot <- function(stack){
  agb <- plot(stack, what = 'agb', ggplot2 = T)
  ba <- plot(stack, what = 'ba', ggplot2 = T)
  abu30 <- plot(stack, what = 'abu30', ggplot2 = T)
  div <- plot(stack, what = 'diversity', ggplot2 = T)
  cowplot::plot_grid(agb, ba, abu30, div, 
                     labels = c('AGB','BA','N30','Div'),
                     nrow = 2, ncol = 2)
}