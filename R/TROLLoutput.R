#'An S4 class to represent TROLL outputs
#'
#'This is an S4 class to represent TROLL outputs.
#'
#'@slot name char. model name
#'@slot abundances list. abundances data frames
#'@slot agb df. agb data frame
#'@slot ba list. ba data frames
#'@slot dbh df. dbh data frame
#'@slot death list. death data frames
#'@slot gpp df. gpp data frame
#'@slot info list. model info
#'@slot litterfall df. litterfall data frame
#'@slot npp df. npp data frame
#'@slot par df. par data frame
#'@slot paramspace list. model space parameters
#'@slot ppfd0 df. ground level ppfd data frames
#'@slot R list. respiration data frames
#'@slot site list. site data frames
#'@slot sp_par df. species par data frame
#'@slot vertd int. vertd through iterations
#'
#'@export
setClass('TROLLoutput',
         representation(
           name = 'character',
           abundances = 'list',
           agb = 'data.frame',
           ba = 'list',
           # cica = 'data.frame', To massive
           dbh = 'data.frame',
           death = 'list',
           # final_pattern, To massive
           gpp = 'data.frame',
           info = 'list',
           # leafdens = 'list', To massive
           litterfall = 'data.frame',
           npp = 'data.frame',
           par = 'list',
           paramspace = 'list',
           ppfd0 = 'data.frame',
           R = 'list',
           site = 'list',
           sp_par = 'data.frame',
           vertd = 'numeric'
         ),
         prototype(
           name = character(),
           abundances = list(),
           agb = data.frame(),
           ba = list(),
           dbh = data.frame(),
           death = list(),
           gpp = data.frame(),
           info = list(),
           litterfall = data.frame(),
           npp = data.frame(),
           par = list(),
           paramspace = list(),
           ppfd0 = data.frame(),
           R = list(),
           site = list(),
           sp_par = data.frame(),
           vertd =  numeric()
         )
)

TROLLoutput <- function(
  name = character(),
  abundances = list(
    adund = data.frame(),
    abu10 = data.frame(),
    abu30 = data.frame()
  ),
  agb = data.frame(),
  ba = list(
    ba = data.frame(),
    ba10 = data.frame()
  ),
  dbh = data.frame(),
  death = list(
    death = data.frame(),
    death1 = data.frame(),
    death2 = data.frame(),
    death3 = data.frame(),
    deathrate = data.frame()
  ),
  gpp = data.frame(),
  info = list(
    step = NA,
    SitesNb = c(NA,NA),
    IterationsNb = NA,
    timestep = NA,
    SpeciesNb = NA,
    ComputationTime = NA
  ),
  litterfall = data.frame(),
  npp = data.frame(),
  par = list(),
  paramspace = list(
    proc	= integer(),
    phi = numeric(),
    k = numeric(),
    fallocwood	= numeric(),
    falloccanopy	= numeric(),
    m	= numeric(),
    m1	= numeric()
  ),
  ppfd0 = data.frame(),
  R = list(
    Rday = data.frame(),
    Rnight = data.frame()
  ),
  site = list(
    site1 = data.frame(),
    site2 = data.frame(),
    site3 = data.frame(),
    site4 = data.frame(),
    site5 = data.frame(),
    site6 = data.frame()
  ),
  sp_par = data.frame(),
  vertd = numeric()
){
  return(new('TROLLoutput',
             name = name,
             abundances = abundances,
             agb = agb,
             ba = ba,
             dbh = dbh,
             death = death,
             gpp = gpp,
             info = info,
             litterfall = litterfall,
             npp = npp,
             par = par,
             paramspace = paramspace,
             ppfd0 = ppfd0,
             R = R,
             site = site,
             sp_par = sp_par,
             vertd = vertd
  )
  )
}
