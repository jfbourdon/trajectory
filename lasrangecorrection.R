Rcpp::sourceCpp("find_range.cpp")

#' Normalize intensity with a range correction
#' 
#' Normalize intensity with a range correction according to the formula (see references): 
#' \deqn{I_{norm} = I_{obs} (\frac{R}{Rs})^f)}
#' To achieve the range correction the position of the sensor must be known at different descrete times. 
#' Using the 'gpstime' of each point, the position of the sensor is interpolated from the reference 
#' and a range correction is applied.
#' 
#' 
#' @template param-las
#' @param flightlines SpatialPointsDataDrame containing the coordinates of the sensor at different
#' instant t. The time is stored in the 'gpstime' attribute and the elevation is stored into 'Z' 
#' attribute. Can be computed with \link{sensor_tracking}.
#' @param Rs numeric. Range of reference.
#' @param f numeric. Exponent. Usually between 2 and 3 in vegetation context.
#' 
#' @return An object of class LAS. The attribute 'Intensity' records the normalised intensity. A extra 
#' attribute 'RawIntensity' records the original intensities (see references).
#' 
#' @references 
#' Gatziolis, D. (2013). Dynamic Range-based Intensity Normalization for Airborne, Discrete Return 
#' Lidar Data of Forest Canopies. Photogrammetric Engineering & Remote Sensing, 77(3), 251–259. 
#' https://doi.org/10.14358/pers.77.3.251
lasrangecorrection = function(las, flightlines, Rs = 1000, f = 2)
{
  if (!"gpstime" %in% names(las@data))         stop("No 'gpstime' attribute found", call. = FALSE)
  if (!"gpstime" %in% names(flightlines@data)) stop("No 'gpstime' attribute found", call. = FALSE)
  if (!"Z" %in% names(flightlines@data))       stop("No 'Z' attribute found", call. = FALSE)
  if (!"Intensity" %in% names(las@data))       stop("No 'Intensity' attribute found", call. = FALSE)
  
  coords <- flightlines@coords
  coords <- data.table::data.table(coords)
  data.table::setnames(coords, c("X", "Y"))
  data   <- data.table::copy(flightlines@data)
  fl     <- cbind(data, coords)
  data.table::setDT(fl)
  data.table::setorder(fl, gpstime)
  
  R <- find_range(las@data, fl) 

  las@data[["RawIntensity"]] <- las@data[["Intensity"]]
  las@data[["Intensity"]]    <- as.integer(las@data[["Intensity"]] * (R/Rs)^f)
  
  #rlas::is_valid_Intensity(las@data, "warning")
  
  return(las)
}
