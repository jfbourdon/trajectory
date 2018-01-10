#' Reconstruct the trajectory of the LiDAR sensor using multiple returns
#'
#' Description longue ici
#' 
#' @details
#' When returns from a single pulse are detected, the sensor compute their positions as being in
#' the center of the footprint and thus being all aligned. Because of that beheaviour, a line
#' drawn between and beyond those returns must cross the sensor. Thus, several consecutive pulses
#' emitted in a tight interval (e.g. 0.001 millisecond) can be used to approximate an intersection
#' point in the sky that would correspond to the sensor position given that the sensor carrier hasn't
#' move much during this interval. A least squares means method using pseudoinverse gives a "close
#' enough" approximation be minimising the squared sum of the distances between the intersection
#' point and all the lines.
#' @param pts.LiDAR: data.table containing the LiDAR data as read with rlas::readlasdata
#' @param PtSourceID: vector containing flight line numbers to use (if NULL, all flight lines 
#' found will be use)
#' @param bin: interval (in second) in which a unqiue sensor position will be find
#' @param step: interval (in second) at which a new sensor position will be find
#' @param nbpairs: minimal number of multiple return pairs needed to estimate a sensor position
#' @author Jean-Francois Bourdon

### Main function to evaluate sensor positions for all flight lines of a LAS file
# Does not currently work as intented because PtSourceID is not used inside sensor_positions().
# Parallel processing should be add inside trajectory() to process all flight lines and return
# df_sensor_positions containing a PtSourceID field.
trajectory <- function(pts.LiDAR, PtSourceID = NULL, bin = 0.001, step = 2, nbpairs = 20) 
{
  if (!"PointSourceID" %in% names(pts.LiDAR))
    stop("No 'PointSourceID' field found", call. = FALSE)
  
  if (is.null(PtSourceID))
    PtSourceID <- fast_unique(pts.LiDAR$PointSourceID)

  data.table::setorder(pts.LiDAR, gpstime, ReturnNumber)
  
  df_sensor_positions <- sensor_positions(pts.LiDAR, bin, step, nbpairs)
  
  return(df_sensor_positions)
}

Rcpp::sourceCpp("C_fn_interval.cpp")

### Evaluate sensor positions for a specific flight line
sensor_positions <- function(pts.LiDAR, bin, step, nbpairs) 
{
  # Time elapsed from first point to last point
  last <- nrow(pts.LiDAR)
  tf <- pts.LiDAR[[last,"gpstime"]] 
  ti <- pts.LiDAR[[1,"gpstime"]]
  elapsed <- (tf - ti)/(bin + step)
  
  # Number of bins
  if (elapsed - floor(elapsed) > bin) 
    nb_bins <- floor(elapsed) + 1
  else 
    nb_bins <- floor(elapsed)
  
  # Generate a list of the ends (line numbers) of each bin
  cumulative_sum <- fast_cumsum_diff(pts.LiDAR[["gpstime"]])
  ls_bins_ends <- bins_ends(nb_bins, bin, step, nbpairs, cumulative_sum)
  
  # Generate a list of XYZ coordinates corresponding to the estimated sensor position
  ls_sensor_positions <- lapply(ls_bins_ends, XYZ_sensor_positions, pts.LiDAR, nbpairs) 
  mat_sensor_positions <- do.call("rbind", ls_sensor_positions)
  
  if(is.null(mat_sensor_positions)) 
    stop("The algorithm failed computing sensor position")
  
  # Generate a data.frame containing sensor position (XYZ), gpstime and nbpairs
  df_sensor_positions <- data.frame(mat_sensor_positions)
  names(df_sensor_positions) <- c("X", "Y", "Z", "gpstime", "NbPairs")
  
  return(df_sensor_positions)
}

### Find indexes (line numbers) of start/end of each bin
bins_ends <- function(nb_bins, bin, step, nbpairs, cumulative_sum)
{
  indexes <- C_fn_interval(nb_bins, bin, step, nbpairs, cumulative_sum)
  indexes_start <- indexes$debut
  indexes_end <- indexes$fin
  
  output <- list()
  
  for (i in 1:nb_bins)
  {
    npulse <- length(unique(cumulative_sum[indexes_start[i]:indexes_end[i]])) + 1
    
    if (npulse > nbpairs) 
      output[[length(output) + 1]] <- c(indexes_start[i], indexes_end[i])
  }
  
  return(output)
}

### Estimate for a specific bin the sensor position
XYZ_sensor_positions <- function(ends_indexes, pts.LiDAR.sourceID, nbpairs)
{
  # Subset containing only returns within a specific bin
  pts.LiDAR_subset <- pts.LiDAR.sourceID[ends_indexes[1]:ends_indexes[2]]
  
  # Lagged difference of gpstime of returns within a specific bin
  diff_gpstime <- diff(pts.LiDAR_subset[["gpstime"]])
  
  # Generate a vector containing the start position (line) of each pulse
  indexes_start_pulses <- c(1, which(diff_gpstime != 0) + 1)
  
  # Generate a vector containing the start position (line) of each pulse with two or more returns
  #  Can't use the NumberOfReturns and ReturnNumber fields in case the dataset has been cleaned
  #  and that some returns are missing
  diff_nb_returns <- diff(indexes_start_pulses)
  indexes_start_pulses_multiple <- which(diff_nb_returns >= 2)
  
  if (length(indexes_start_pulses_multiple) >= nbpairs) 
  {
    indexes_multiple_start <- indexes_start_pulses[indexes_start_pulses_multiple]
    indexes_multiple_end <- indexes_multiple_start + (diff_nb_returns[indexes_start_pulses_multiple] - 1)
    
    # data.table containing coordinates XYZ of each start/end of each multiple returns pulse
    dt_start <- pts.LiDAR_subset[indexes_multiple_start, c("X","Y","Z")]
    dt_end <- pts.LiDAR_subset[indexes_multiple_end, c("X","Y","Z")]
    
    XYZ_start <- as.matrix(dt_end)
    XYZ_end   <- as.matrix(dt_start)
    time <- pts.LiDAR_subset[[1, "gpstime"]]
    
    # Matrix containing sensor position (XYZ), gpstime and number of pulses used
    mat_XYZ_positions <- XYZ_intersect(XYZ_start, XYZ_end)
    output <- cbind(mat_XYZ_positions, time, nrow(XYZ_start))
    
    return(output)
  }
}

### Calculate intersection point from several returns pairs using least squares means
XYZ_intersect <- function(XYZ_start, XYZ_end)
{
  # XYZ_start: matrix n x 3 containing XYZ coordinates for each starting returns
  # XYZ_end: matrix n x 3 containing XYZ coordinates for each ending returns
  # Translated from MATLAB (Anders Eikenes, 2012)
  # http://www.mathworks.com/matlabcentral/fileexchange/37192-intersection-point-of-lines-in-3d-space
  
  direction_vectors <- XYZ_end - XYZ_start
  length_vectors <- matrix(sqrt(.rowSums(direction_vectors^2, dim(direction_vectors)[1], 3)))
  normalized_vectors <- direction_vectors/(length_vectors %*% matrix(1, ncol = 3))
  
  Nx <- normalized_vectors[,1, drop = FALSE]
  Ny <- normalized_vectors[,2, drop = FALSE]
  Nz <- normalized_vectors[,3, drop = FALSE]
  
  Sxx <- sum(Nx^2 - 1)
  Syy <- sum(Ny^2 - 1)
  Szz <- sum(Nz^2 - 1)
  Sxy <- sum(Nx*Ny)
  Sxz <- sum(Nx*Nz)
  Syz <- sum(Ny*Nz)
  
  S <- matrix(c(Sxx,Sxy,Sxz,Sxy,Syy,Syz,Sxz,Syz,Szz), nrow = 3)
  
  Cx <- sum(XYZ_start[,1, drop = FALSE]*(Nx^2 - 1) + XYZ_start[,2, drop = FALSE]*(Nx*Ny) + XYZ_start[,3, drop = FALSE]*(Nx*Nz))
  Cy <- sum(XYZ_start[,1, drop = FALSE]*(Nx*Ny) + XYZ_start[,2, drop = FALSE]*(Ny^2 - 1) + XYZ_start[,3, drop = FALSE]*(Ny*Nz))
  Cz <- sum(XYZ_start[,1, drop = FALSE]*(Nx*Nz) + XYZ_start[,2, drop = FALSE]*(Ny*Nz) + XYZ_start[,3, drop = FALSE]*(Nz^2 - 1))
  
  C <- matrix(c(Cx,Cy,Cz))
  mat_XYZ_intersect <- t(solve(S,C))
  
  return(mat_XYZ_intersect)
}

# # ========================================================== 
# # A RETIRER SI NE SERT PAS ? UNE QUELCONQUE VALIDATION
# # distances perpendiculaires entre le point et les lignes
# 
# dist3d <- function(ii, a, b, c) 
# {
#   v1 <- b[ii,] - c[ii,]
#   v2sp: <- a[1,] - b[ii,]      
#   v3 <- cross3d_prod(v1,v2)
#   area <- sqrt(sum(v3*v3))/2
#   d <- 2*area/sqrt(sum(v1*v1))
#   return(d)
# }
# 
# cross3d_prod <- function(v1, v2)
# {
#   v3 <- vector()
#   v3[1] <- v1[2]*v2[3]-v1[3]*v2[2]
#   v3[2] <- v1[3]*v2[1]-v1[1]*v2[3]
#   v3[3] <- v1[1]*v2[2]-v1[2]*v2[1]
#   return(v3)
# }
# 
# N <- dim(PA)[1]
# 
# #Distances from intersection point to the input lines
# distances<-sapply(1:N, dist3d, P_intersect, PA, PB) 
# 
# # func_ECDF<-ecdf(distances)
# # 
# # eval.max<-max(distances)
# # eval.min<-min(distances)
# # limite<-0.8
# # 
# # while (abs(mean(c(eval.max, eval.min))-eval)>0.001) {
# #   eval<-mean(c(eval.max, eval.min))
# #   if(func_ECDF(eval)<limite) {
# #     eval.min<-eval
# #   } else {
# #     eval.max<-eval
# #   }
# # }
# # 
# # plot(func_ECDF, main="Distribution cumulative des distances au point", ylab="Proportion", xlab="Distance (m)")
# # abline(h=limite, lty=2, col="red")
# # abline(v=eval, lwd=2, col="red")
# # text(1, paste0(sprintf("%.1f", eval)," m @ ", limite), pos=1, col="red", font=2)
# ## FIN de ? RETIRER SI NE SERT PAS ? UNE QUELCONQUE VALIDATION
# # ======================================================================================
