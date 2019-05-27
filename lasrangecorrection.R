lasrangecorrection = function(las, flightlines, Rs = 1000)
{
  fl <- as.data.frame(flightlines)
  data.table::setDT(fl)
  data.table::setorder(fl, gpstime)
  data.table::setkey(fl, gpstime)
  
  ids <- fl[J(las$gpstime), roll = "nearest", which = TRUE]
  R   <- (las$X - fl$X[ids])^2 + (las$Y - fl$Y[ids])^2 + (las$Z - fl$Z[ids])^2
  Rs  <- Rs^2
  las@data$RawIntensity <- las@data$Intensity
  las@data$Intensity <- las@data$Intensity * R/Rs
  return(las)
}
