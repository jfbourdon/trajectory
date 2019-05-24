#' Reconstruct the trajectory of the LiDAR sensor using multiple returns
#'
#' Description longue ici
#' 
#' @details
#' When returns from a single pulse are detected, the sensor compute their positions as being in
#' the center of the footprint and thus being all aligned. Because of that beheaviour, a line
#' drawn between and beyond those returns must cross the sensor. Thus, several consecutive pulses
#' emitted in a tight interval (e.g. 0.001 second) can be used to approximate an intersection
#' point in the sky that would correspond to the sensor position given that the sensor carrier hasn't
#' moved much during this interval. A weighed least squares method using pseudoinverse gives a "close
#' enough" approximation be minimising the squared sum of the distances between the intersection
#' point and all the lines.
#' @param las An object of class LAS containing the LiDAR data as read with \link[lidR:readLAS]{readLAS}
#' @param bin numeric interval (in second) in which a unique sensor position will be find
#' @param min_length numeric minimum length that vectors from the combination of Start and
#' End must have to be kept in the analysis
#' @param nbpairs minimum number of multiple return pairs needed to estimate a sensor position
#' @author Jean-Francois Bourdon

### Evaluate sensor positions from LiDAR points. Points must all have the same PointSourceID.
sensor_tracking <- function(las, bin = 0.5, min_length = 2, nbpairs = 500)
{
  if (!"PointSourceID" %in% names(las@data))
    stop("No 'PointSourceID' attribute found", call. = FALSE)
  
  if (!"gpstime" %in% names(las@data))
    stop("No 'gpstime' attribute found", call. = FALSE)
  
  if (!"ReturnNumber" %in% names(las@data))
    stop("No 'ReturnNumber' attribute found", call. = FALSE)
  
  if (!"NumberOfReturns" %in% names(las@data))
    stop("No 'NumberOfReturns' attribute found", call. = FALSE)
  
  if (!is.numeric(min_length) || min_length < 0)
    stop("Invalid min_length argument. Must be a positive numeric")

  data <- las@data
  
  # Reordering of input data by gpstime and ReturnNumber
  data.table::setorder(data, gpstime, ReturnNumber)
  
  # Get only the first and last returns of multiple returns
  data$pulseID <- lidR:::.lagisdiff(data[["gpstime"]])
  data <- data[(ReturnNumber == NumberOfReturns | ReturnNumber == 1) & NumberOfReturns > 1]
  
  # Filter some edge points that may not be paired
  count <- lidR:::fast_table(data$pulseID,  max(data$pulseID))
  ii    <- which(count == 1)
  data  <- data[!pulseID %in% ii]
  
  # Generate the bins
  bins <- lidR:::round_any(data$gpstime, bin)
  
  # Find the position P of the sensor in each bin
  P <- data[, if(.N > 2*nbpairs) sensor_positions(X,Y,Z, ReturnNumber, min_length), by = .(bins, PointSourceID)]
  P <- sp::SpatialPointsDataFrame(P[,3:5], P[,c(1,2,6)])
  
  return(P)
}

### Estimate for a specific bin the sensor position
sensor_positions <- function(x, y, z, rn, min_length, weights = NULL)
{
  first <- rn == 1
  last  <- rn > 1
  Start <- matrix(c(x[first], y[first], z[first]), ncol = 3L)
  End   <- matrix(c(x[last], y[last], z[last]), ncol = 3L)
  
  # Start:  matrix n x 3 containing XYZ coordinates for each starting returns
  # End:    matrix n x 3 containing XYZ coordinates for each ending returns
  # min_length: minimum length that vectors from the combination of Start and
  #             and End must have to be kept in the analysis
  # weights:    determine if a weights matrix will be computed
  #              -> NULL         : compute based on vectors length
  #              -> 1            : equal weights
  #              -> matrix object: user-defined weights matrix
  
  # Translated and adapted from MATLAB (Anders Eikenes, 2012)
  # http://www.mathworks.com/matlabcentral/fileexchange/37192-intersection-point-of-lines-in-3d-space
  
  Direction <- End - Start
  length_vectors <- matrix(sqrt(.rowSums(Direction^2, dim(Direction)[1], 3)))
  
  # Filtering out vectors shorter than the specified minimum length
  indexes_valid <- length_vectors >= min_length
  Start <- Start[indexes_valid,, drop = FALSE]
  Direction <- Direction[indexes_valid,, drop = FALSE]
  length_vectors <- length_vectors[indexes_valid,, drop = FALSE]
  
  # Validation of weights argument
  if (is.null(weights))
    weights <- length_vectors / sum(length_vectors)
  else if (length(weights) == 1L)
    weights <- matrix(weights, nrow = nrow(Start))
  
  normalized_vectors <- Direction / (length_vectors %*% matrix(1, ncol = 3))
  
  Nx <- normalized_vectors[,1, drop = FALSE]
  Ny <- normalized_vectors[,2, drop = FALSE]
  Nz <- normalized_vectors[,3, drop = FALSE]
  
  Wxx <- weights * (Nx^2 - 1)
  Wyy <- weights * (Ny^2 - 1)
  Wzz <- weights * (Nz^2 - 1)
  Wxy <- weights * Nx * Ny
  Wxz <- weights * Nx * Nz
  Wyz <- weights * Ny * Nz
  
  Sxx <- sum(Wxx)
  Syy <- sum(Wyy)
  Szz <- sum(Wzz)
  Sxy <- sum(Wxy)
  Sxz <- sum(Wxz)
  Syz <- sum(Wyz)
  
  S <- matrix(c(Sxx,Sxy,Sxz,Sxy,Syy,Syz,Sxz,Syz,Szz), nrow = 3)
  
  Cx <- sum(Start[,1, drop = FALSE] * Wxx + Start[,2, drop = FALSE] * Wxy + Start[,3, drop = FALSE] * Wxz)
  Cy <- sum(Start[,1, drop = FALSE] * Wxy + Start[,2, drop = FALSE] * Wyy + Start[,3, drop = FALSE] * Wyz)
  Cz <- sum(Start[,1, drop = FALSE] * Wxz + Start[,2, drop = FALSE] * Wyz + Start[,3, drop = FALSE] * Wzz)
  
  C <- matrix(c(Cx,Cy,Cz))
  M <- t(solve(S,C))
  
  return(list(X = M[1], Y = M[2], Z = M[3], N = nrow(Start)))
}