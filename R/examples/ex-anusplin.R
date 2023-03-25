# data(TempBrazil)
# colnames(TempBrazil) <- c("lon", "lat", "temp")
# df = TempBrazil
# 
# # range <- c(70, 140, 15, 55)
# range <- c(-78, -34, -36, 5)
# p = anusplin_params(df, "TempBrazil", "dem.txt", range, alt = NULL)
# str(p)

library(data.table)
library(magrittr)

dat <- fread('inst/extdata/anusplin_demo/RH_01-01.dat') %>% 
  `colnames<-`(c('site', 'lon', 'lat', 'alt', 'RH'))

range <- c(69.625, 140.375, 14.625, 55.375)

params <- anusplin_params(dat, 'RH', range, 'inst/extdata/anusplin_demo/china_dem_025deg.txt')
anusplin_write('R/examples/anusplin_output', params, is.run = T)