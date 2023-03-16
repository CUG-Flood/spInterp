library(data.table)
library(terra)
library(glue)
library(plyr)

#' Create configures for ANUSPLIN
#'
#' Format the input data and generate the configuration file required for ANUSPLIN interpolation.
#'
#' @param dat A data.frame or data.table to interpolate, for the colnames, lon,
#' lat must be included, site and alt are optional, the others are variable names.
#' @param basename Basename for all output files, extensions should not be included.
#' @param range Range of interpolation grid (xmin, xmax, ymin, ymax).
#' @param file.alt When `type.alt` is not 0, this is the input grid file name; otherwise
#' a manual constant.
#' @param unit The unit of `dat`, a non-negative integer, possible values are:
#' \describe{
#' \item{`0`}{undefined (**default**)}
#' \item{`1`}{meteres}
#' \item{`2`}{feet}
#' \item{`3`}{kilometers}
#' \item{`4`}{miles}
#' \item{`5`}{degrees}
#' \item{`6`}{radians}
#' \item{`7`}{millimetres}
#' \item{`8`}{megajoules}
#' }
#' @param res The grid resolution (degree).
#' @param width The fixed width of numbers in formatted data.
#' @param alt Type of elevation was treated, possible values are:
#' \describe{
#' \item{`NULL`}{no use of elevation}
#' \item{`cov`}{considered as independent covariates (**default**)}
#' \item{`spl`}{considered as independent spline variables}
#' }
#' @param lim.lon A vector containing lower and upper limits, `auto` (default)
#' meant to use the minimum and maximum values in the data, or set manually.
#' Data points outside these limits, augmented by margins, are ignored.
#' @param lim.lat Same as `lim.lon`, but for longitude.
#' @param lim.alt Same as `lim.lon`, but for altitude.
#' @param cvt.lon Transformation and scale factor (**default** is `1`) of longitude.\cr
#' Real Value = (Table Value) * (Scale Factor)\cr
#' The possible transformations are:
#' \describe{
#' \item{`0`}{no transformation (**default**)}
#' \item{`1`}{x/a}
#' \item{`2`}{ax}
#' \item{`3`}{a·log(x+b)}
#' \item{`4`}{(x/b)^a}
#' \item{`5`}{a·exp(x/b)}
#' \item{`6`}{a·tanh(x/b)}
#' \item{`7`}{anisotropy angle in degrees}
#' \item{`8`}{anisotropy factor - in the direction specified by the anisotropy angle}
#' }
#' @param cvt.lat Same as `cvt.lon`, but for latitude.
#' @param cvt.alt Same as `cvt.lon`, but for altitude.
#' @param cvt.coef Parameters used for transformation, one or two real numbers.
#' @param trans.dep Dependent variable transformation, possible values are:
#' \describe{
#' \item{`0`}{no transformation (**default**)}
#' \item{`1`}{fit surface to natural logarithm of the data values}
#' \item{`2`}{fit surface to the square root of the data values}
#' \item{`5`}{occurrence – transform data values by setting all positive value to 1.0 and ignoring all negative values}
#' }
#' @param order Order of spline, a positive integer.
#' @param err.wgt Number of relative error variances, a non-negative integer,
#' possible values are:
#' \describe{
#' \item{`0`}{data points uniformly weighted for each surface (**default**)}
#' \item{`1`}{the same weighting is applied to each surface}
#' \item{Number of surfaces}{a different weighting is applied to each surface}
#' }
#' @param optimize Optimization directive, a non-negative integer, possible values are:
#' \describe{
#' \item{`0`}{common smoothing parameter for all surfaces}
#' \item{`1`}{common smoothing directive for all surfaces (**default**)}
#' \item{`2`}{different smoothing directive for each surface}
#' }
#' @param smooth Smoothing directive for each surface, a non-negative integer,
#' possible values are:
#' \describe{
#' \item{`0`}{fixed smoothing parameter - supply value}
#' \item{`1`}{ minimise GCV (**default**)}
#' \item{`2`}{minimise true mean square error using supplied error standard deviation estimate}
#' \item{`3`}{fixed signal - supply value}
#' \item{`4`}{minimise GML}
#' }
#' @param type.mask Mode of mask grid, a non-negative integer, possible values are:
#' \describe{
#' \item{`0`}{mask grid not supplied (**default**)}
#' \item{`1`}{generic mask grid}
#' \item{`2`}{Arc/Info mask grid}
#' \item{`3`}{Idrisi mask grid}
#' }
#' @param file.mask Filename of mask grid, only valid if `type.mask` set to positive
#' integer.
#' @param type.alt Mode of the independent variable, possible values are:
#' \describe{
#' \item{0}{user supplied constant}
#' \item{1}{user supplied grid in generic row format with the same size as the grid being calculated}
#' \item{2}{user supplied Arc/Info grid with same size as the grid being calculated (**default**)}
#' \item{3}{user supplied Idrisi image with the same size as the grid being calculated}
#' }
#' @param type.grd Same as `type.mask`, but for interpolated grid.
#' @param missing Filling of missing values.
#' @param err.cov test
#' @param grid.pos Grid position option, a non-negative integer, possible values are:
#' \describe{
#' \item{`0`}{grid points at cell corners (**default**)}
#' \item{`1`}{grid points at cell centres}
#' }
#' @param essential If `True`, only export essential process files, large residual
#' file, optimisation parameters file, data list file and validation data file are
#' ignored.
#' @return a list with three components:
#' \describe{
#' \item{data}{fixed length of data.frame (data.table)}
#' \item{splina}{splina parameters}
#' \item{lapgrd}{lapgrd parameters}
#' }
#' @importFrom  glue glue
#'
#' @examples opt_anusplin(dat, 'test', c(70, 140, 15, 55))
opt_anusplin <- function(dat,
                         basename,
                         range,
                         file.alt,
                         unit = 0,
                         res = 0.25,
                         width = 7, 
                         alt = 'cov',
                         lim.lon = 'auto',
                         lim.lat = 'auto',
                         lim.alt = 'auto',
                         cvt.lon = c(0, 1),
                         cvt.lat = c(0, 1),
                         cvt.alt = c(1, 1),
                         cvt.coef = 1000,
                         trans.dep = 0,
                         order = 3,
                         err.wgt = 0,
                         optimize = 1,
                         smooth = 1,
                         type.mask = 0,
                         file.mask = NULL,
                         type.alt = 2,
                         type.grd = 2,
                         missing = -9999,
                         err.cov = 2,
                         grid.pos = 1,
                         essential = T) {
  # check arguments validation
  if (type.alt != 0) {
    if (!is.character(file.alt))
      stop('1')
  } else {
    if (!is.numeric(file.alt))
      stop('2')
  }
  
  names <- colnames(dat)
  
  if (is.null(alt)) {
    if ('alt' %in% names) {
      dat$alt <- NULL
      names <- names[names != 'alt']
    }
    cvt.alt <- NULL
    
  } else {
    cvt.alt <- if (lim.alt == 'auto') {
      dat$alt %>% range %>% c(cvt.alt) %>% paste(collapse = ' ')
    } else {
      paste(lim.alt, cvt.alt, collapse = ' ')
    }
  }
  
  # n.sur <- length(names) - 2
  
  ind.spl <- 2
  ind.cov <- sur.spl <- sur.cov <- 0
  
  cvt.lon <- if (lim.lon == 'auto') {
    dat$lon %>% range %>% c(cvt.lon) %>% paste(collapse = ' ')
  } else {
    paste(lim.lon, cvt.lon, collapse = ' ')
  }
  
  cvt.lat <- if (lim.lat == 'auto') {
    dat$lat %>% range %>% c(cvt.lat) %>% paste(collapse = ' ')
  } else {
    paste(lim.lat, cvt.lat, collapse = ' ')
  }
  
  if ('site' %in% names) {
    pos.site <- which('site' == names)
    chr.site <- dat$site %>% nchar %>% max
    fmt.site <- glue('%{chr.site}i')
    dat$site %<>% sprintf(fmt = fmt.site)
    # n.sur <- n.sur - 1
  } else {
    chr.site <- 0
    pos.site <- -1
    fmt.site <- NULL
  }
  
  pos.coord = which(names %in% c('lon', 'lat'))
  if (length(pos.coord) != 2) {
    stop('Check whether lon lat in colnames of dat')
  }
  dat[, pos.coord] <-
    dat[, lapply(.SD, sprintf, fmt = glue('%{width}.2f')), .SDcols = pos.coord]
  
  if (!is.null(alt) & 'alt' %in% names) {
    dat$alt %<>% sprintf(fmt = glue('%{width}.1f'))
    # n.sur <- n.sur - 1
    if (alt == 'cov') {
      ind.cov <- ind.cov + 1
    } else if (alt == 'spl') {
      ind.spl <- ind.spl + 1
    }
  }
  
  
  pos.var <- which(!names %in% c('site', 'lon', 'lat', 'alt'))
  n.sur <- length(pos.var)
  dat[, pos.var] <-
    dat[, lapply(.SD, sprintf, fmt = glue('%{width}.2f')), .SDcols = pos.var]
  
  file.sur <- glue('{basename}.sur')
  file.cov <- glue('{basename}.cov')
  if (essential) {
    output <- c('', '', file.sur, '', file.cov)
  } else {
    output <- c(
      glue('{basename}.res'),
      glue('{basename}.opt'),
      file.sur,
      glue('{basename}.lis'),
      file.cov
    )
  }
  
  rle <- str_extract(dat[1], "(?<=\\.)[0-9]+") %>% nchar %>% rle
  rle$lengths[which(rle$lengths == 1)] <- ''
  rle$values[which(!is.na(rle$value))] %<>% paste0('f', width, '.', .)
  rle$values[which(is.na(rle$values))] <-
    glue('a{chr.site}')
  
  fmt <- paste0(rle[['lengths']], rle[['values']]) %>%
    paste(collapse = ',') %>%
    paste0('(', ., ')')
  
  splina <- c(
    basename,
    unit,
    ind.spl,
    ind.cov,
    sur.spl,
    sur.cov,
    cvt.lon,
    cvt.lat,
    cvt.alt,
    cvt.coef,
    paste(trans.dep, collapse = ' '),
    order,
    n.sur,
    err.wgt,
    optimize,
    smooth,
    glue('{basename}.dat'),
    round_any(nrow(dat), 5, f = ceiling),
    chr.site,
    fmt,
    output,
    '',
    ''
  )
  
  lapgrd <- c(
    file.sur,
    0,
    1,
    file.cov,
    err.cov,
    '',
    grid.pos,
    pos.coord[1],
    glue('{range[1]} {range[2]} {res}'),
    pos.coord[2],
    glue('{range[3]} {range[4]} {res}'),
    ifelse(type.mask, c(type.mask, file.mask), type.mask),
    type.alt,
    file.alt,
    type.grd,
    missing,
    glue('{names[pos.var]}.grd'),
    glue('({n.sur + 2}f{width}.2)'),
    type.grd,
    missing,
    glue('cov_{names[pos.var]}.grd'),
    glue('({n.sur + 2}f{width}.2)'),
    '',
    ''
  )
  
  return(list(
    data = dat,
    splina = splina,
    lapgrd = lapgrd
  ))
}


write_anusplin <-
  function(dat,
           opt_splina,
           opt_lapgrd,
           file_path,
           na.width = 'auto',
           names = c('splina.txt', 'lapgrd.txt'),
           exe = c('splina.exe', 'lapgrd.exe'),
           cmd = F) {
    
    width <- data[, lapply(.SD, nchar)] %>% unique %>% nrow
    if (width != 1) stop('Data column width is not fixed.')
    
    if (na.width == 'auto') {
      width <- opt_splina[length(opt_splina) - 7] %>%
        str_extract_all('(?<=f)[0-9]+', simplify = T)
      if (uniqueN(fmt) != 1) {
        stop(glue('Conflicted parse, the detected width is not unique, there are {width}'))
      }
    }
    
    file <- opt_splina[length(opt_splina) - 10]
    write.table(
      dat,
      glue("{file_path}/{file}"),
      sep = "",
      col.names = F,
      row.names = F, 
      quote = F,
      na = strrep(" ", width)
    )
    
    splina <- unlist(opt_splina, use.names = F)
    fileConn <- glue("{file_path}/{names[1]}") %>% file()
    writeLines(splina, fileConn)
    close(fileConn)
    
    lapgrd <- unlist(opt_lapgrd, use.names = F)
    fileConn <- glue("{file_path}/{names[2]}") %>% file()
    writeLines(lapgrd, fileConn)
    close(fileConn)
    
    if ((length(exe) == 2) && cmd) {
      cmd <- c(
        glue('{exe[1]}<{names[1]}>splina.log'),
        glue('{exe[2]}<{names[2]}>lapgrd.log')
      )
      fileConn <- glue("{file_path}/{cmd}") %>% file()
      writeLines(lapgrd, fileConn)
      close(fileConn)
    }
  }
