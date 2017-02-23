#' @importFrom ggplot2 ggplot aes geom_point xlab ylab scale_colour_hue theme_bw theme scale_x_continuous
NULL

#' Compare plot
#' 
#' function to compare
#'
#' @param sim list of TROLL.output. TROLL model outputs to compare
#' @param output char. Output to compare (abu10, abu30, agb and gpp available)
#' @param legend logical. add legend
#'
#' @return ggplot2 object
#' 
#' @export
#'
#' @examples
#' 
compareplot <- function(sim, output = 'agb', legend = TRUE){
  time <- seq(1,sim[[1]]@par$general$nbiter,1)/sim[[1]]@par$general$iter
  if(grepl('abu', output)){
    table <- data.frame(time = rep(time, length(names(sim))),
                        output = unlist(lapply(sim, function(x) x@abundances[[which(names(x@abundances) == output)]]$Total)),
                        trait = rep(names(sim), each = length(time)))
  } else {
    table <- data.frame(time = rep(time, length(names(sim))),
                        output = unlist(lapply(sim, function(x) slot(x, output)$Total)),
                        trait = rep(names(sim), each = length(time)))
  }
  lab <- switch(output,
                'agb' = 'Aboveground biomass (tonnes/ha)',
                'abu10' = "Number of trees with dbh > 10 cm (in stems/ha)",
                'abu30' = "Number of trees with dbh > 30 cm (in stems/ha)",
                'gpp' = "Total GPPLeaf (in MgC/ha)"
                )
  g <- ggplot(table, aes(x = time, y = output, colour = trait)) +
    geom_point(size=2, fill = "white") +
    xlab('Times (years)') +
    scale_x_continuous(breaks = seq(0, max(time), length.out = sim[[1]]@par$general$nbiter/sim[[1]]@par$general$iter+1)) +
    ylab(lab) +
    theme_bw()
  if(legend)
    g <- g + scale_colour_hue(name="Trait fixed to mean", l=40)
  else
    g <- g + theme(legend.position="none")

  return(g)
}