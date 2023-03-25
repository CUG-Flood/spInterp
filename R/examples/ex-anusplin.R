library(data.table)
library(magrittr)

dat <- fread('inst/extdata/anusplin_demo/RH_01-01.dat') %>% 
  `colnames<-`(c('site', 'lon', 'lat', 'alt', 'RH'))

range <- c(69.625, 140.375, 14.625, 55.375)

params <- anusplin_params(dat, 'RH', range, 'inst/extdata/anusplin_demo/china_dem_025deg.txt')
anusplin_write('R/examples/anusplin_output', params, is.run = T)