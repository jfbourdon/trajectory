Rcpp::sourceCpp("find_range.cpp")

lasrangecorrection = function(las, flightlines, Rs = 1000, f = 2)
{
  coords <- flightlines@coords
  coords <- data.table::data.table(coords)
  data.table::setnames(coords, c("X", "Y"))
  data   <- data.table::copy(flightlines@data)
  fl     <- cbind(data, coords)
  data.table::setDT(fl)
  data.table::setorder(fl, gpstime)
  
  R <- find_range(las@data, fl) 

  las@data$RawIntensity <- las@data$Intensity
  las@data$Intensity <- las@data$Intensity * (R/Rs)^f
  return(las)
}
