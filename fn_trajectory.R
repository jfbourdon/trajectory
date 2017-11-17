#' Reconstitue une trajectoire de vol estimée à partir des retours multiples
#'
#' Description longue ici
#' 
#' @details
#' Explication de l'algo ici
#' @param pts.LiDAR: data.table containing the data
#' @param PtSourceID: vecteur des numéros de lignes de vol à traiter (si 0, l'algo le fera 
#' pour toutes les lignes de vols trouv?es)
#' @param bin: intervalle maximal en seconde dans lequel les retours doivent tous se situ?s 
#' pour qu'une position XYZt unique leur soit attribu?e
#' @param saut: intervalle en seconde auquel une nouvelle position doit ?tre trouv?e
#' @param nbpairs: nombre minimal de paires de retours multiples devant ?tre trouv? pour 
#' qu'une position soit estim?e
#' @author Jean-Francois Bourdon
fn_trajectory <- function(pts.LiDAR, PtSourceID = NULL, bin = 0.001, step = 2, nbpairs = 20) 
{
  if (!"PointSourceID" %in% names(pts.LiDAR))
    stop("No 'PointSourceID' field found", call. = FALSE)
  
  if (is.null(PtSourceID))
    PtSourceID <- fast_unique(pts.LiDAR$PointSourceID)

  data.table::setorder(pts.LiDAR, gpstime, ReturnNumber)
  
  SPDF.XYZtp <- fn_SPDF.XYZtp(pts.LiDAR, bin, step, nbpairs)
  
  return(SPDF.XYZtp)
}

Rcpp::sourceCpp("C_fn_interval.cpp")

### ?valuation de la trajectoire de vol pour un PointSource sp?cifique
fn_SPDF.XYZtp<-function(pts.LiDAR, bin, step, nbpairs) 
{
  fin <- nrow(pts.LiDAR)
  
  # Time ellapsed from first point to last point
  tf   <- pts.LiDAR[[fin,"gpstime"]] 
  ti   <- pts.LiDAR[[1,"gpstime"]]
  ellapsed <- (tf - ti)/(bin+step)
  
  # Number of bins
  if (ellapsed - floor(ellapsed) > bin) 
    nb.ellapsed <- floor(ellapsed)+1
  else 
    nb.ellapsed <- floor(ellapsed)
  
  # g?n?ne une somme cumulative des diff?rences de temps afin d'obtenir exactement l'index 
  # d'une impulsion X seconde apr?s une autre
  somme.cumulative <- fast_cumsum_diff(pts.LiDAR[["gpstime"]])
  
  # g?n?re une liste des index pour chaque bin
  #ls.index <- lapply(1:nb.ellapsed, fn_index, bin, step, nbpairs, somme.cumulative) 
  ls.index <- fn_index2(nb.ellapsed, bin, step, nbpairs, somme.cumulative)
  
  # pour chaque liste des index, trouve les retours multiples appartenant ? la m?me impulsion, 
  # trace un prolongement dans le ciel et trouve ultimement les coordonn?es XYZ du point de 
  # rencontre de toutes ces lignes par la technique des moindres carr?s
  ls.XYZ<-lapply(ls.index, fn_XYZ.l2m.complet, pts.LiDAR, nbpairs) 
  
  mat.XYZt<-do.call("rbind", ls.XYZ)
  
  if(is.null(mat.XYZt)) 
    stop("The algorithm failed computing aircraft position")
  
  df.XYZtp <- data.frame(cbind(mat.XYZt, id))
  names(df.XYZtp) <- c("X", "Y", "Z", "gpstime", "NbPaires", "PointSourceID")
  XYZtp.SPDF <- sp::SpatialPointsDataFrame(coords=df.XYZtp[,c("X","Y")], data = df.XYZtp)
  
  return(XYZtp.SPDF) #SpatialPointsDataFrame en sortie avec XYZ + gpstime + PointSourceID
}

fn_index2 = function(nb.ellapsed, bin, step, nbpairs, somme.cumulative)
{
  indexes = C_fn_interval(nb.ellapsed, bin, step, nbpairs, somme.cumulative)
  indexes.debut = indexes$debut
  indexes.fin = indexes$fin
  
  output = list()
  
  for (i in 1:nb.ellapsed)
  {
    npulse = length(unique(somme.cumulative[indexes.debut[i]:indexes.fin[i]]))+1
    
    if (npulse > nbpairs)
      output[[length(output) + 1]] <- c(indexes.debut[i], indexes.fin[i])
  }
  
  return(output)
}

### Trouve les index (num?ro de ligne) du d?but et de la fin de chaque bin
fn_index <- function(ii, bin, step, nbpairs, somme.cumulative) 
{
  intervalle.debut <- (bin+step)*(ii-1)
  index.debut <- match(TRUE, somme.cumulative > intervalle.debut)
  intervalle.fin <- intervalle.debut+bin
  index.fin <- match(TRUE, somme.cumulative>intervalle.fin)
  
  nb.pulse <- length(unique(somme.cumulative[index.debut:index.fin]))+1
  
  # imposition d'un seuil minimal de nombre de retours pour poursuivre l'?valuation
  if (nb.pulse > nbpairs) 
    return(c(index.debut, index.fin))
  else
    return(NULL)
}

fn_XYZ.l2m.complet<-function(index, pts.LiDAR.sourceID, nbpairs) 
{
  index.debut <- index[1]
  index.fin   <- index[2]
  
  # sous-ensemble comprenant uniquement les points de la bo?te
  pts.LiDAR.subset<-pts.LiDAR.sourceID[index.debut:index.fin]
  
  #calcule la diff?rence entre chaque paire d'?l?ments cons?cutifs
  diff.pts.LiDAR<-diff(pts.LiDAR.subset[["gpstime"]])
  
  # trouve la position du d?but de chaque impulsion
  pos.impulsion<-c(1, which(diff.pts.LiDAR != 0) + 1) 
  
  diff.impulsion<-diff(pos.impulsion)
  index.diff.multiple<-which(diff.impulsion >= 2)
  
  if (length(index.diff.multiple) >= nbpairs) 
  {
    index.multiple.debut<-pos.impulsion[index.diff.multiple]
    index.multiple.fin<-index.multiple.debut+(diff.impulsion[index.diff.multiple]-1)
    
    dt.debut <- pts.LiDAR.subset[index.multiple.debut, c("X","Y","Z")]
    dt.fin <- pts.LiDAR.subset[index.multiple.fin, c("X","Y","Z")]
    
    PA   <- as.matrix(dt.fin)   #XYZ des points de d?part
    PB   <- as.matrix(dt.debut) #XYZ des points d'arriv?e
    time <- pts.LiDAR.subset[[1, "gpstime"]]
    
    positions = fn_XYZ.l2m(PA, PB)
    
    return(cbind(positions, time, nrow(PA)))
  }
}

fn_XYZ.l2m <- function(PA, PB) 
{
  ### Calcul du point d'intersection XYZ par la technique des moindres carr?s.
  ### PA: matrice n x 3 contenant les coordonn?es XYZ de chaque point de d?part
  ### PB: matrice n x 3 contenant les coordonn?es XYZ de chaque point de d'arriv?e
  ### Code traduit de MATLAB en R
  ### http://www.mathworks.com/matlabcentral/fileexchange/37192-intersection-point-of-lines-in-3d-space
  
  Si <- PB-PA #N lines described as vectors, direction vector
  ni <- Si/(matrix(sqrt(.rowSums(Si^2, dim(Si)[1], dim(Si)[2]))) %*% matrix(1, nrow=1, ncol=3)) #Normalize vectors
  
  nx <- ni[,1, drop=FALSE]
  ny <- ni[,2, drop=FALSE]
  nz <- ni[,3, drop=FALSE]
  
  SXX <- sum(nx^2-1)
  SYY <- sum(ny^2-1)
  SZZ <- sum(nz^2-1)
  SXY <- sum(nx*ny)
  SXZ <- sum(nx*nz)
  SYZ <- sum(ny*nz)
  
  S <- matrix(c(SXX,SXY,SXZ,SXY,SYY,SYZ,SXZ,SYZ,SZZ), nrow = 3, byrow = TRUE)
  
  CX <- sum(PA[,1, drop=FALSE]*(nx^2-1) + PA[,2, drop=FALSE]*(nx*ny) + PA[,3, drop=FALSE]*(nx*nz))
  CY <- sum(PA[,1, drop=FALSE]*(nx*ny) + PA[,2, drop=FALSE]*(ny^2-1) + PA[,3, drop=FALSE]*(ny*nz))
  CZ <- sum(PA[,1, drop=FALSE]*(nx*nz) + PA[,2, drop=FALSE]*(ny*nz)  + PA[,3, drop=FALSE]*(nz^2-1))
  
  C <- matrix(c(CX,CY,CZ))
  P_intersect <- Conj(t(solve(S,C))) #Best intersection point of the N lines, in least squares sense
  
  return(P_intersect)
  #return(cbind(P_intersect, mean(distances), sd(distances)))
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