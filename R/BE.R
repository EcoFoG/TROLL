#' @export
getBE <- function(sim,
                  mature,
                  monoculture,
                  vars = c('agb', 'ba', 'n', 'n10', 'n30',
                           'gpp', 'npp', 'Rday', 'Rnight'),
                  reduce = TRUE){
  # Abundances
  mature <- mature@layers[[
    which(names(mature@layers) == 
            paste0(strsplit(sim@name, '.', fixed = T)[[1]][1:2], collapse = '.'))]]
  species <- row.names(sim@sp_par)
  abund <- table(species[mature@final_pattern$sp_lab])
  abund <- sort(abund/sum(abund))
  
  # Species DeltaRY
  d <- strsplit(sim@name, '.', fixed = T)[[1]][3]
  DeltaRY <- array(dim = c(600,length(vars),length(species))) # time*DelatRvars*sp
  M <- array(dim = c(600,length(vars),length(species))) # time*DelatRvars*sp
  for(i in seq_len(length(species))){
    sp <- species[i]
    monosp <- monoculture@layers[[paste0(sp, '.', d)]]
    datamonosp <- extractData(stack(monosp), vars = vars, reduce = reduce)[-c(1:2)]
    M[,,i] <- as.matrix(datamonosp)
    datamixturesp <- extractData(stack(sim), sp = sp, vars = vars, reduce = reduce)[-c(1:2)]
    DeltaRY[,,i] <- as.matrix(cbind(datamixturesp/datamonosp - abund[sp]))
  } ; rm(i, reduce, sp, datamixturesp, datamonosp)
  
  # BNE
  CE <- length(species)*apply(DeltaRY, c(1,2), mean, na.rm = T)*apply(M, c(1,2), mean, na.rm = T)
  SE <- t(length(species)* mapply(function(x,y){
    mapply(function(xp, yp){
      cov(xp,yp)
    }, xp <- plyr::alply(x,1), yp <- plyr::alply(y,1))
  }, plyr::alply(DeltaRY,1), plyr::alply(M,1)))
  NE <- CE + SE
  BE <- list(NE = NE, CE = CE, SE = SE) ; rm(NE, CE, SE)
  BE <- lapply(BE, function(l){
    l <- data.frame(l)
    names(l) <- vars
    sim <- sim@name
    time <- 1:600
    l <- cbind(time, sim, l)
    return(l)})
  BE$NE <- cbind(BE = rep('NE',600), BE$NE)
  BE$CE <- cbind(BE = rep('CE',600), BE$CE)
  BE$SE <- cbind(BE = rep('SE',600), BE$SE)
  BE <- do.call('rbind', BE)
  
  return(BE)
}
