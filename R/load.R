#' @include TROLLoutput.R
NULL


#' Function to load TROLL output
#'
#' @param name char. Name given to the model output
#' @param path char. Path where the model is saved
#'
#' @return an S4 \linkS4class{TROLLoutput} class object
#'
#' @export
#'
#' @examples
#' name = 'test2'
#' path = './src/OUTPUT/test2/'
#' load(name, path)
#'
load <- function(name, path = getwd()){

  vertd <- read.table(paste0(path, name, '_0_vertd.txt'))[,2]
  names(vertd) <- read.table(paste0(path, name, '_0_vertd.txt'))[,1]

  TROLLoutput(
    name = name,
    abundances = list(
      adund = read.table(paste0(path, name, '_0_abund.txt'), row.names = 1),
      abu10 = read.table(paste0(path, name, '_0_abund.txt'), row.names = 1),
      abu30 = read.table(paste0(path, name, '_0_abund.txt'), row.names = 1)
    ),
    agb = read.table(paste0(path, name, '_0_agb.txt')),
    ba = list(
      ba = read.table(paste0(path, name, '_0_ba.txt'), row.names = 1),
      ba10 = read.table(paste0(path, name, '_0_ba10.txt'), row.names = 1)
    ),
    dbh = data.frame(), # Don't get it must ask
    death = list(
      death = read.table(paste0(path, name, '_0_death.txt'), row.names = 1),
      death1 = read.table(paste0(path, name, '_0_death1.txt')),
      death2 = read.table(paste0(path, name, '_0_death2.txt')),
      death3 = read.table(paste0(path, name, '_0_death3.txt')),
      deathrate = read.table(paste0(path, name, '_0_deathrate.txt'))
    ),
    gpp = read.table(paste0(path, name, '_0_gpp.txt'), row.names = 1),
    info = list(
      step = scan(paste0(path, name, '_0_info.txt'), character(), skip = 4, n = 7, quiet = T)[7],
      SitesNb = unlist((strsplit(scan(paste0(path, name, '_0_info.txt'),
                                      character(), skip = 12, n = 5, quiet = T)[5], 'x'))),
      IterationsNb = scan(paste0(path, name, '_0_info.txt'), character(), skip = 13, n = 5, quiet = T)[5],
      timestep = scan(paste0(path, name, '_0_info.txt'), character(), skip = 14, n = 5, quiet = T)[5],
      SpeciesNb = scan(paste0(path, name, '_0_info.txt'), character(), skip = 15, n = 5, quiet = T)[5],
      ComputationTime = scan(paste0(path, name, '_0_info.txt'), character(), skip = 17, n = 5, quiet = T)[5]
    ),
    litterfall = read.table(paste0(path, name, '_0_litterfall.txt')),
    npp = read.table(paste0(path, name, '_0_npp.txt'), row.names = 1),
    par = list(
      general = list(),
      species = data.frame(),
      climate = data.frame()
    ), # All parameters !
    paramspace = list(
      proc	= as.integer(scan(paste0(path, name, '_0_paramspace.txt'), character(), 2, quiet = TRUE)[2]),
      phi = as.numeric(scan(paste0(path, name, '_0_paramspace.txt'), character(), 2, skip = 1, quiet = TRUE)[2]),
      k = as.numeric(scan(paste0(path, name, '_0_paramspace.txt'), character(), 2, skip = 2, quiet = TRUE)[2]),
      fallocwood	= as.numeric(scan(paste0(path, name, '_0_paramspace.txt'), character(), 2, skip = 3, quiet = TRUE)[2]),
      falloccanopy	= as.numeric(scan(paste0(path, name, '_0_paramspace.txt'), character(), 2, skip = 4, quiet = TRUE)[2]),
      m	= as.numeric(scan(paste0(path, name, '_0_paramspace.txt'), character(), 2, skip = 5, quiet = TRUE)[2]),
      m1	= as.numeric(scan(paste0(path, name, '_0_paramspace.txt'), character(), 2, skip = 6, quiet = TRUE)[2])
    ),
    ppfd0 = read.table(paste0(path, name, '_0_ppfd0.txt'), row.names = 1),
    R = list(
      Rday = read.table(paste0(path, name, '_0_Rday.txt'), row.names = 1),
      Rnight = read.table(paste0(path, name, '_0_Rnight.txt'), row.names = 1)
    ),
    site = list(
      site1 = read.table(paste0(path, name, '_0_site1.txt'), row.names = 1),
      site2 = read.table(paste0(path, name, '_0_site2.txt'), row.names = 1),
      site3 = read.table(paste0(path, name, '_0_site3.txt'), row.names = 1),
      site4 = read.table(paste0(path, name, '_0_site4.txt'), row.names = 1),
      site5 = read.table(paste0(path, name, '_0_site5.txt'), row.names = 1),
      site6 = read.table(paste0(path, name, '_0_site6.txt'), row.names = 1)
    ),
    sp_par = read.table(paste0(path, name, '_0_sp_par.txt')),
    vertd = vertd
  )
}
