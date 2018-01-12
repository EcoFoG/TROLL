# # Monoculture 1 
# load('~/Documents/ECOFOG/Results/disturbance/monoculture/mature/mmData.Rdata')
# mmdata <- datal$time
# load('~/Documents/ECOFOG/Results/disturbance/monoculture/disturbance/mdData.Rdata')
# mddata <- datal$time
# mddata$time <- mddata$time + max(mmdata$time)
# data <- rbind(mmdata, mddata)
# ggplot(data[grep('Vouacapoua', data$sim),], aes(x = time, y = agb, color = sim)) + 
#   geom_point()
# 
# # Monoculture 2
# # load('~/Documents/ECOFOG/Results/disturbance/monoculture/mature/mm.Rdata')
# # mature <- mm@layers$Vouacapoua_americana ; rm(mm)
# # load('~/Documents/ECOFOG/Results/disturbance/monoculture/disturbance/md.Rdata')
# # t0 <- md@layers$Vouacapoua_americana.0 ; rm(md)
# # time <- plot(mature, what = 'agb', ggplot = T)$data$x
# # mdata <- data.frame(time = time, agb = mature@agb$Total)
# # ddata <- data.frame(time = time+max(time), agb = t0@agb$Total)
# # data <- rbind(mdata, ddata)
# # ggplot(data, aes(x = time, y = agb)) + 
# #   geom_point() + 
# #   geom_vline(xintercept = 600)
# 
# # Monoculture disturbance is perfectly working !
# 
# # Mixture 1
# load('~/Documents/ECOFOG/Results/disturbance/maturity/matureData.Rdata')
# mdata <- datal$time
# load('~/Documents/ECOFOG/Results/disturbance/disturbance/disturbanceData.Rdata')
# ddata <- datal$time
# ddata$time <- ddata$time + max(mdata$time)
# data <- rbind(mdata, ddata)
# # # All plots
# # for(x in c(5,25,125)){
# #   for(y in 1:20){
# #     g <- ggplot(data[which(data$sim %in% paste0(paste0(x,'.',y), 
# #                                                 c('', '.0', '.25', '.50', '.75'))),], 
# #                 aes(x = time, y = agb, color = sim)) + 
# #       geom_point() + geom_vline(xintercept = 600)
# #     plot(g)
# #     readline(prompt="Press [enter] to continue")
# #   }
# # }
# g1 <- ggplot(data[which(data$sim %in% paste0('25.5', c('', '.0', '.25', '.50', '.75'))),], 
#              aes(x = time, y = agb, color = sim)) + 
#   geom_point() + geom_vline(xintercept = 600)
# g2 <- ggplot(data[which(data$sim %in% paste0('125.14', c('', '.0', '.25', '.50', '.75'))),], 
#              aes(x = time, y = agb, color = sim)) + 
#   geom_point() + geom_vline(xintercept = 600)
# cowplot::plot_grid(g1, g2)
# 
# # Mixture 2
# load('~/Documents/ECOFOG/Results/disturbance/maturity/mature.Rdata')
# mature.5 <- stack(mature@layers[which(names(mature@layers) %in% paste0('5.',1:20))])
# plot(mature.5, what = 'diversity', ggplot = T)
# mature.25 <- stack(mature@layers[which(names(mature@layers) %in% paste0('25.',1:20))])
# plot(mature.25, what = 'diversity', ggplot = T)
