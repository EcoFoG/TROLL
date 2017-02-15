#' @import methods
NULL

#' Function to plot TROLL outputs.
#'
#' @param x TROLLoutput
#' @param y char. Output to plot, see details.
#'
#' @return Plot the output
#'
#' @details Available plots
#' \describe{
#' \item{defaults}{Not implemented}
#' \item{abund}{To precise}
#' \item{agb}{To precise}
#' \item{ba}{To precise}
#' \item{death}{To precise}
#' \item{final pattern}{To precise}
#' \item{flux}{To precise}
#' \item{gpp}{To precise}
#' \item{LAI}{To precise}
#' \item{litterfall}{To precise}
#' \item{npp}{To precise}
#' \item{R}{To precise}
#' \item{Ra}{To precise}
#' \item{relabund}{To precise}
#' \item{vertd}{To precise}
#' \item{wd}{To precise}
#' }
#'
#' @export
#'
#' @examples
#'
setMethod('plot', 'TROLLoutput', function(x, y = NULL, ...) {

  ##### data ####
  nbiter <- x@par$general$nbiter
  iter <- x@par$general$iter
  age <- round(nbiter / iter, 1) # in years
  surf <- prod(x@info$step * x@info$SitesNb) / 10000 # in ha

  switch (y,

    ##### default ####
    NULL = {
      warning('Defalut not implemented !')
    },

    ##### abund ####
    'abund' = {
      plot(seq(1,nbiter,1)/iter,
           x@abundances$abund$Total, pch=20,
           main="Total abundance (in stems/ha)",
           xlab="Times (in years)",
           ylab="Number of stems/ha",
           xlim=c(0,nbiter/iter))
    },

    'abu10' = {
      plot(seq(1,nbiter,1)/iter,
           x@abundances$abu10$Total, pch=20,
           main="Number of trees with dbh > 10 cm (in stems/ha)",
           xlab="Times (in years)",
           ylab="Number of stems/ha",
           xlim=c(0,nbiter/iter))
    },

    'abu30' = {
      plot(seq(1,nbiter,1)/iter,
           x@abundances$abu30$Total, pch=20,
           main="Number of trees with dbh > 30 cm (in stems/ha)",
           xlab="Times (in years)",
           ylab="Number of stems/ha",
           xlim=c(0,nbiter/iter))
    },

    'allabund' = {
      par(mfrow=c(3,1))
      plot(x, 'abund')
      plot(x, 'abu10')
      plot(x, 'abu30')
      par(mfrow=c(1,1))
    },

    ##### agb ####
    'agb' = {
      plot(seq(1,nbiter,1)/iter,
           x@agb$Total, pch=20,
           main="Total basal area (in m2/ha)",
           xlab="Aboveground biomass (in tonnes/ha)",
           ylab="Aboveground biomass (in tonnes/ha)",
           xlim=c(0,nbiter/iter))
    },

    ##### ba ####
    'ba' = {
      plot(seq(1,nbiter,1)/iter,
           x@ba$ba$Total, pch=20,
           main="Total basal area (in m2/ha)",
           xlab="Times (in years)",
           ylab="Basal area (in m2/ha)",
           xlim=c(0,nbiter/iter))
    },

    'ba10' = {
      plot(seq(1,nbiter,1)/iter,
           x@ba$ba10$Total, pch=20,
           main="Basal area of trees with dbh > 10 cm (in stems/ha)",
           xlab="Times (in years)",
           ylab="Basal area (in m2/ha)",
           xlim=c(0,nbiter/iter))
    },

    'allba' = {
      par(mfrow=c(2,1))
      plot(x, 'ba')
      plot(x, 'ba10')
      par(mfrow=c(1,1))
    },

    ##### death ####
    'death' = {
      warning('death plots not implemented yet !')
    },

    ##### final pattern ####

    'height hist' = {
      hist <- hist(x@final_pattern$height[x@final_pattern$height != 0], breaks=seq(0,50,by=1), plot = FALSE)
      plot(hist,
          main = paste("Tree height histogram after", age, "years (", surf, "ha)"),
          col="green",
          xlim=c(0,30),
          ylim=c(0,40000),
          xlab = "Height class")
    },

    'species count' = {
      barplot(table(x@final_pattern$sp_lab), horiz = T)
    },

    'age' = {
      plot(x@final_pattern['age'],
           col=rev(heat.colors(10)),
           xlim=c(0,x_max),
           ylim=c(0,y_max),
           main="TROLL age distribution",
           xlab="x (m)",
           ylab="y (m)"
      )
    },

    'species' = {
      plot(x@final_pattern['sp_lab'],
           # col=rev(heat.colors(10)),
           col = rainbow(length(unique(final_pattern$sp_lab))),
           xlim=c(0,x_max),
           ylim=c(0,y_max),
           main="TROLL species distribution",
           xlab="x (m)",
           ylab="y (m)"
      )
    },

    ##### flux ####
    'flux' = {
      plot(seq(1,nbiter,1)/iter,
           x@gpp$Total, pch=20,
           xlab="Times (in years)",
           ylab="Total flux (in MgC/ha)",
           main="Flux (in MgC/ha)")
      points(seq(1,nbiter,1)/iter,
           x@npp$Total, pch=20, col = 'red')
      points(seq(1,nbiter,1)/iter,
           x@gpp$Total - x@npp$Total,
           pch=20, col = 'green')
      legend("topright", col=c("black", "green",  "red"), pch=20,
             legend=c("GPP", "Autotrophic respiration", "NPP"))
    },

    ##### gpp ####
    'gpp' = {
      plot(seq(1,nbiter,1)/iter,
           x@gpp$Total, pch=20,
           xlab="Times (in years)",
           ylab="Total GPPLeaf (in MgC/ha)",
           main="Gross primary productivity (in MgC/ha)")
    },

    ##### LAI ####
    'LAI' = {
      warning('LAI plots not implemented yet !')
    },

    ##### litterfall ####
    'litterfall' = {
      plot(seq(1,nbiter,1)/iter,
           x@litterfall$Total * 12, pch=20,
           xlab="Times (in years)",
           ylab="Leaf litterfall per month (in Mg/ha/year dry mass)",
           main="Leaf litterfall per month (in Mg dry mass/ha/year)")
    },

    ##### npp ####
    'npp' = {
      plot(seq(1,nbiter,1)/iter,
           x@npp$Total, pch=20,
           xlab="Times (in years)",
           ylab="Total NPPLeaf (in MgC/ha)",
           main="Net primary productivity (in MgC/ha)")
    },

    ##### R ####
    'R' = {
      warning('R plots not implemented yet !')
    },


    ##### Ra ####
    'Ra' = {
      plot(seq(1,nbiter,1)/iter,
           x@gpp$Total - x@npp$Total, pch=20,
           xlab="Times (in years)",
           ylab="Autotrophic respiration (in MgC/ha)",
           main="Autotrophic respiration (in MgC/ha)")
    },

    ##### relabund ####
    'relabund' = {
      warning('relabund plots not implemented yet !')
    },

    ##### vertd ####

    'vertd' = {
      max_height <- max(x@vertd$height)
      plot(x =c(-diff(tail(x@vertd$vertd, max_height)),0),
           y = seq(1, max_height),
           type = "l",
           main = paste("PAD distribution after", age, "years"),
           xlab = "PAD (m2 per m3)",
           ylab = "Height",
           xlim = c(0,1),
           ylim = c(0,60))
    },

    ##### wd ####
    'wd' = {
      mean_wood_dens <- rep(0,nbiter)
      for (i in 1:nbiter) {
        mean_wood_dens[i] <- sum(x@abundances$relabdund[i,] / 100 * x@sp_par$wsg)
      }
      plot(seq(1,nbiter,1)/iter,
           mean_wood_dens,
           ylim=c(0, 0.7), pch=20,
           main="Total wood density (in g/cm3)",
           xlab="Time (in year)",
           ylab="Average plot wood density",
           xlim=c(0,nbiter/iter))
    },

    'wd10' = {
      mean_wood_dens <- rep(0,nbiter)
      for (i in 1:nbiter) {
        mean_wood_dens[i] <- sum(x@abundances$relabu10[i,] / 100 * x@sp_par$wsg)
      }
      plot(seq(1,nbiter,1)/iter,
           mean_wood_dens,
           ylim=c(0, 0.7), pch=20,
           main="Wood density of trees with dbh > 10 cm (in g/cm3)",
           xlab="Time (in year)",
           ylab="Average plot wood density",
           xlim=c(0,nbiter/iter))
    },

    'wd30' = {
      mean_wood_dens <- rep(0,nbiter)
      for (i in 1:nbiter) {
        mean_wood_dens[i] <- sum(x@abundances$relabu30[i,] / 100 * x@sp_par$wsg)
      }
      plot(seq(1,nbiter,1)/iter,
           mean_wood_dens,
           ylim=c(0, 0.7), pch=20,
           main="Wood density of trees with dbh > 30 cm (in g/cm3)",
           xlab="Time (in year)",
           ylab="Average plot wood density",
           xlim=c(0,nbiter/iter))
    },

    'allwd' = {
      par(mfrow=c(3,1))
      plot(x, 'wd')
      plot(x, 'wd10')
      plot(x, 'wd30')
      par(mfrow=c(1,1))
    },

    ##### else ####
    warning('Not implemented yet.')
  )
})
