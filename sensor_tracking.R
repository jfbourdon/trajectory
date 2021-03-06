#' Reconstruct the trajectory of the LiDAR sensor using multiple returns
#'
#' Use multiple returns to estimate the positionning of the sensor by computing the intersection in 
#' space of the line passing througt the first and last returns. To work this function requieres a 
#' dataset where the 'gpstime', 'ReturnNumber', 'NumberOfReturns' and 'PointSourceID' attributes are 
#' properly populated otherwise the output may be incorrect or weird.
#' 
#' When multiple returns from a single pulse are detected, the sensor compute their positions as being 
#' in the center of the footprint and thus being all aligned. Because of that beheaviour, a line
#' drawn between and beyond those returns must cross the sensor. Thus, several consecutive pulses
#' emitted in a tight interval (e.g. 0.1 second) can be used to approximate an intersection
#' point in the sky that correspond to the sensor position given that the sensor carrier hasn't
#' moved much during this interval. A weighed least squares method gives an  approximation of the 
#' intersection by minimising the squared sum of the distances between the intersection point and all 
#' the lines.
#' 
#' @section Test of data integrity:
#' In theory the sensor tracking is a simple problem to solve as long as each pulse is properly 
#' identified from a well populated dataset. In practice many problems may arise from wrongly populated 
#' datasets. Here a list of problem that may happens. Those with a * denote already encountered problem and 
#' internally checked:
#' \itemize{
#' \item 'gpstime' do not record the time at with pulses were emitted and thus pulse are not idenfiable
#' \item *A pulse (two or more points that share the same gpstime) is made of points from different 
#' flightlines (different PointSourceID). This is impossible and denote wrongly populated PointSourceID
#' attribute.
#' }
#' 
#' @template param-las
#' @param interval numeric. Tight interval used to bin the gps times and group the pulses. 
#' @param pmin interger. Minimum number of pulses needed to estimate a sensor position.
#' For a given interval, the sensor position is not computed if the number of pulse is lower than 
#' \code{pmin}.
#' @param extra_check boolean. Dataset are rarely perfectly populated leading to unexpected errors. 
#' Time consuming checks of data integrity are performed. These checks can be skipped as they account 
#' for most of the computation time. See also section 'Tests of data integrity'
#' 
#' @return A SpatialPointDataFrame with the Z elevation stored in the table of attribute. Informations 
#' about the time interval and the number of pulses used to find the points are also in the table of 
#' attributes.
#' 
#' @author Jean-Francois Bourdon & Jean-Romain Roussel
sensor_tracking <- function(las, interval = 0.5, pmin = 200, extra_check = TRUE)
{
  UseMethod("sensor_tracking", las)
}

sensor_tracking.LAS <- function(las, interval = 0.5, pmin = 200, extra_check = TRUE)
{
  if (!"PointSourceID" %in% names(las@data))     stop("No 'PointSourceID' attribute found", call. = FALSE)
  if (!"gpstime" %in% names(las@data))           stop("No 'gpstime' attribute found", call. = FALSE)
  if (!"ReturnNumber" %in% names(las@data))      stop("No 'ReturnNumber' attribute found", call. = FALSE)
  if (!"NumberOfReturns" %in% names(las@data))   stop("No 'NumberOfReturns' attribute found", call. = FALSE)
  if (!any(las@data[["gpstime"]] != 0))          stop("'gpstime' attribute is not populated.", call. = FALSE)
  if (!any(las@data[["ReturnNumber"]] != 0))     stop("'ReturnNumber' attribute is not populated.", call. = FALSE)
  if (!any(las@data[["NumberOfReturns"]] != 0))  stop("'NumberOfReturns' attribute is not populated.", call. = FALSE)
  if (!any(las@data[["PointSourceID"]] != 0))    stop("'PointSourceID' attribute is not populated.", call. = FALSE)

  data <- las@data
  
  # Reordering of input data by gpstime and ReturnNumber
  data.table::setorder(data, gpstime, ReturnNumber)
  
  # Compute an ID for each pulse
  data$pulseID <- lidR:::.lagisdiff(data[["gpstime"]])
  
  # Get only the first and last returns of multiple returns
  data <- data[(ReturnNumber == NumberOfReturns | ReturnNumber == 1) & NumberOfReturns > 1]
  
  # Filter some edge points that may not be paired in a pulse
  count <- lidR:::fast_table(data$pulseID,  max(data$pulseID))
  ii    <- which(count == 1L)
  data  <- data[!pulseID %in% ii]
  
  if (extra_check)
  {
    tests <- data[, list(test1 = .N != 2L, test2 = PointSourceID[1] != PointSourceID[2]), by = .(pulseID)]
    
    multiple_source <- tests$test1
    unpaired_pulse  <- tests$test2

    # This happens if two points share the same 'gpstime' but different 'PointSourceID'
    if (any(multiple_source))
      warning(glue::glue("{sum(multiple_source)} pulses (points with same gpstime) come from different flightlines. The point cloud is likely to be wrongly populated. These points were removed"), call. = FALSE)
    
    # Does this may really happens?
    if (any(unpaired_pulse))
      warning(glue::glue("{sum(unpaired_pulse)} pulses with multiple returns were not actually paired. The point cloud is likely to be wrongly populated. These points were removed"), call. = FALSE)
    
    ii    <- tests$pulseID[unpaired_pulse | multiple_source]
    data  <- data[!pulseID %in% ii]
  }
  
  # Generate the bins
  bins <- lidR:::round_any(data$gpstime, interval)
  
  # Find the position P of the sensor in each interval
  P <- data[, if (.N > 2*pmin) sensor_positions(X,Y,Z, ReturnNumber), by = .(gpstime = bins, PointSourceID = PointSourceID)]
  P <- P[!is.na(X)]
  P <- sp::SpatialPointsDataFrame(P[,3:4], P[,c(5,1,2,6)])
  
  return(P)
}

sensor_tracking.LAScluster <- function(las, interval = 0.5, pmin = 200, extra_check = TRUE)
{
  x <- readLAS(las)
  if (is.empty(x)) return(NULL)
  pos  <- sensor_tracking(x, interval, pmin, extra_check)
  #bbox <- raster::extent(las)
  #pos  <- raster::crop(pos, bbox)
  return(pos)
}

sensor_tracking.LAScatalog <- function(las, interval = 0.5, pmin = 200, extra_check = TRUE)
{
  opt_select(las) <- "xyzrntp"
  opt_filter(las) <- paste("-drop_single", opt_filter(las))
  
  options <- list(need_buffer = TRUE, drop_null = TRUE, need_output_file = FALSE)
  output  <- catalog_apply(las, sensor_tracking, interval = interval, pmin = pmin, extra_check = extra_check, .options = options)
  output  <- do.call(rbind, output)
  output@proj4string <- las@proj4string
  
  data <- as.data.frame(output)
  data.table::setDT(data)
  i <- data[, .I[which.max(npulses)], by = gpstime]$V1
  return(output[i,])
}

# Translated and adapted from MATLAB (Anders Eikenes, 2012)
# http://www.mathworks.com/matlabcentral/fileexchange/37192-intersection-point-of-lines-in-3d-space
sensor_positions <- function(x, y, z, rn)
{
  first <- rn == 1L
  last  <- rn > 1L
  
  # Matrix n x 3 containing XYZ coordinates for each first/last return
  A <- matrix(c(x[first], y[first], z[first]), ncol = 3L)
  B <- matrix(c(x[last], y[last], z[last]), ncol = 3L)
  
  # That should never happen but it might happen in wronly populated dataset
  if (nrow(A) != nrow(B))
  {
    warning("Something went wrong. The point cloud is likely to be wrongly populated in a way not tested internally.", call. = FALSE)
    return(list(X = NA_real_, Y = NA_real_, Z = NA_real_, npulses = NA_integer_))
  }
  
  # Lines described as vectors V
  V <- B - A
  
  # Weights lines by distances 
  D <- matrix(sqrt(.rowSums(V^2, nrow(V), ncol(V))))
  W <- D / sum(D)
  
  # Normalized vectors (|V| = 1)
  ones <- matrix(1, ncol = 3)
  N  <- V / (D %*% ones)
  Nx <- N[, 1L, drop = FALSE]
  Ny <- N[, 2L, drop = FALSE]
  Nz <- N[, 3L, drop = FALSE]
  
  Wxx <- W * (Nx^2 - 1)
  Wyy <- W * (Ny^2 - 1)
  Wzz <- W * (Nz^2 - 1)
  Wxy <- W * Nx * Ny
  Wxz <- W * Nx * Nz
  Wyz <- W * Ny * Nz
  
  Sxx <- sum(Wxx)
  Syy <- sum(Wyy)
  Szz <- sum(Wzz)
  Sxy <- sum(Wxy)
  Sxz <- sum(Wxz)
  Syz <- sum(Wyz)
  S   <- matrix(c(Sxx, Sxy, Sxz, Sxy, Syy, Syz, Sxz, Syz, Szz), nrow = 3L)
  
  Cx <- sum(A[,1, drop = FALSE] * Wxx + A[,2, drop = FALSE] * Wxy + A[,3, drop = FALSE] * Wxz)
  Cy <- sum(A[,1, drop = FALSE] * Wxy + A[,2, drop = FALSE] * Wyy + A[,3, drop = FALSE] * Wyz)
  Cz <- sum(A[,1, drop = FALSE] * Wxz + A[,2, drop = FALSE] * Wyz + A[,3, drop = FALSE] * Wzz)
  C  <- matrix(c(Cx,Cy,Cz))
  
  M  <- t(solve(S,C))
  
  return(list(X = M[1], Y = M[2], Z = M[3], npulses = nrow(A)))
}