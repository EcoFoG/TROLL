#' @export
sylvicultureInit <- function(sp_par,
                             input,
                             path,
                             harvest = NULL,
                             disturb_iter = 1,
                             disturb_intensity = 0.25,
                             designated_volume = 30,
                             harvested_volume = 25,
                             overwrite = TRUE){
  # Creating and opening file
  if(!overwrite)
    if(input %in% list.files(path))
      stop('The file already exist, use overwrite = T.')
  fileConn <- file(file.path(path, input))
  
  # Writing file
  writeLines(c(
    '##########################################################',
    '###  Sylviculture parameter file for the TROLL program ###',
    '##########################################################',
    '###	GENERAL PARAMETERS',
    paste(disturb_iter,	'/* disturb_iter # iteration step where the disturbation occure */', sep = '\t'),
    '###	disturbance',
    paste(disturb_intensity,	'/* disturb_intensity # intensity of disturbance in percent of BA */', sep = '\t'),
    '###	logging',
    paste(designated_volume,	'/* designated_volume # volume designated for harvesting in m3/ha */', sep = '\t'),
    paste(harvested_volume,	'/* harvested_volume # volume harvested in m3/ha */', sep = '\t'),
    paste(dim(sp_par)[1],	'/* numespharvestable # number of harvestable species */', sep = '\t'),
    '					\n### Species'
  ), fileConn)
  close(fileConn)
  
  # Harvest table
  if(is.null(harvest)){
    harvest <- data.frame(
      species = c('Brosimum_guianense', 'Brosimum_rubescens', 'Caryocar_glabrum',
                  'Dicorynia_guianensis', 'Goupia_glabra', 'Manilkara_bidentata',
                  'Manilkara_huberi', 'Ocotea_argyrophylla', 'Qualea_rosea', 
                  'Ruizterania_albiflora', 'Vouacapoua_americana',
                  'Couma_guianensis', 'Eperua_grandiflora', 'Lecythis_zabucajo',
                  'Moronobea_coccinea', 'Rhodostemonodaphne_grandis',
                  'Sterculia_pruriens', 'Sterculia_speciosa',
                  'Sterculia_villifera', 'Virola_michelii',
                  'Vochysia_guianensis', 'Bocoa_prouacensis',
                  'Couratari_multiflora', 'Eperua_falcata',
                  'Eperua_rubiginosa'),
      dbh_min = 0.55,
      dbh_max = 2,
      interest = c(rep(1,11), rep(2,10), rep(3,4))
    )
    row.names(harvest) <- harvest$species
    harvest['Vouacapoua_americana','dbh_max'] <- 0.70
  }
  harvest$sp_lab <- sp_par$sp_lab[match(harvest$species, row.names(sp_par))]
  harvest <- harvest[-which(is.na(harvest$sp_lab)),]
  harvest <- harvest[c('species', 'sp_lab', 'dbh_min', 'dbh_max', 'interest')]
  names(harvest)[1:2] <- c('****', 'spnum')
  suppressWarnings(
    write.table(harvest, file.path(path, input), sep = "\t", append = TRUE,
                row.names = FALSE, quote = FALSE)
  )
    
}