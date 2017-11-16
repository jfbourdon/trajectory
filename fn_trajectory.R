fn_trajectory<-function(path, PtSourceID=0, bin=0.001, step=2, nbpairs=20) {
  
  #### Reconstitue une trajectoire de vol estimée à partir des retours multiples
  #### Dernière modification 2017-11-10
  #### Jean-François Bourdon
  #### MFFP - DIF - Division image et photointerprétation 
  
  # path: chemin d'accès au fichier LAS/LAZ à traiter
  # PtSourceID: vecteur des numéros de lignes de vol à traiter (si 0, l'algo le fera pour toutes les lignes de vols trouvées)
  # bin: intervalle maximal en seconde dans lequel les retours doivent tous se situés pour qu'une position XYZt unique leur soit attribuée
  # step: intervalle en seconde auquel une nouvelle position doit être trouvée
  # nbpairs: nombre minimal de paires de retours multiples devant être trouvé pour qu'une position soit estimée
  
  ### AMÉLIORATIONS À APPORTER ###
  # Ajouter une option pour ignorer les paires dont les points sont distants de moins de x.x mètre
  # Utiliser des moindres carrés biaisés en donnant plus de poids aux paires dont la longueur de leur
  #    segement est plus grande. Plus les deux points sont éloignés, moins l'imprécision de la position
  #    des points eux mêmes aura un impact sur l'angle de l'angle de la projection dans le ciel.
  # Ajouter une option pour XX % des paires dont la distance au point milieu sont les plus éloignées
  #    - nécessite d'ajouter une option déterminer le nombre de passes
  #    - recalculer à chaque passe un nouveau point milieu
  
  require(rlas)
  
  pts.LiDAR<-readlasdata(path, filter="-drop_single") #importation du fichier LAZ en excluant les retours uniques
  pts.LiDAR<-pts.LiDAR[order(gpstime, ReturnNumber)] #data.table des points LiDAR en ordre croissant du gpstime puis du ReturnNumber
  
  if (PtSourceID==0) {PtSourceID<-unique(pts.LiDAR[["PointSourceID"]])}
  
  ls.SPDF.XYZtp<-lapply(seq_along(PtSourceID), fn_SPDF.XYZtp, pts.LiDAR, PtSourceID, bin, step, nbpairs)
  
  ls.SPDF.XYZtp[sapply(ls.SPDF.XYZtp, is.null)]<-NULL
  SPDF.XYZtp<-do.call(rbind, ls.SPDF.XYZtp)
  
  return(SPDF.XYZtp)
}


fn_SPDF.XYZtp<-function(jj, pts.LiDAR, PtSourceID, bin, step, nbpairs) {
  ### Évaluation de la trajectoire de vol pour un PointSource spécifique
  
  require(sp)
  
  pts.LiDAR.sourceID<-pts.LiDAR[PointSourceID==PtSourceID[jj],] #data.table des points LiDAR d'un PointSourceID spécifique
  fin<-nrow(pts.LiDAR.sourceID) #nombre total de retours pour la ligne de vol
  eval<-(pts.LiDAR.sourceID[[fin,"gpstime"]]-pts.LiDAR.sourceID[[1,"gpstime"]])/(bin+step) #calcul du nombre de secondes entre les premier et dernier point de la ligne de vol 
  if (eval-floor(eval)>bin) {nb.eval<-floor(eval)+1} else {nb.eval<-floor(eval)} #calcule le nombre de boîtes qu'il sera possible d'évaluer selon les paramètres
  somme.cumulative<-cumsum(diff(pts.LiDAR.sourceID[["gpstime"]])) #génène une somme cumulative des différences de temps afin d'obtenir exactement l'index d'une impulsion X seconde après une autre
  
  ls.index<-lapply(1:nb.eval, fn_index, bin, step, nbpairs, somme.cumulative) #génère une liste des index pour chaque bin
  
  ls.XYZ<-lapply(1:nb.eval, fn_XYZ.l2m.complet, pts.LiDAR.sourceID, ls.index, nbpairs) #pour chaque liste des index, trouve les retours multiples appartenant à la même impulsion, trace un prolongement dans le ciel et trouve ultimement les coordonnées XYZ du point de rencontre de toutes ces lignes par la technique des moindres carrés
  
  mat.XYZt<-do.call("rbind", ls.XYZ)
  
  if(!is.null(mat.XYZt)) {
    df.XYZtp<-data.frame(cbind(mat.XYZt, PtSourceID[jj]))
    names(df.XYZtp)<-c("X", "Y", "Z", "meanDist", "sdDist", "gpstime", "NbPaires", "PointSourceID")
    XYZtp.SPDF<-SpatialPointsDataFrame(coords=df.XYZtp[,c("X","Y")], data=df.XYZtp)
    
    return(XYZtp.SPDF) #SpatialPointsDataFrame en sortie avec XYZ + gpstime + PointSourceID
  }
}


fn_index<-function(ii, bin, step, nbpairs, somme.cumulative) {
  ### Trouve les index (numéro de ligne) du début et de la fin de chaque bin
  
  intervalle.debut<-(bin+step)*(ii-1)
  index.debut<-match(TRUE, somme.cumulative>intervalle.debut)
  intervalle.fin<-intervalle.debut+bin
  index.fin<-match(TRUE, somme.cumulative>intervalle.fin)
  
  nb.pulse<-length(unique(somme.cumulative[index.debut:index.fin]))+1
  
  if (nb.pulse>nbpairs) {return(c(index.debut, index.fin))} #imposition d'un seuil minimal de nombre de retours pour poursuivre l'évaluation
}


fn_XYZ.l2m.complet<-function(kk, pts.LiDAR.sourceID, ls.index, nbpairs) {
  ### l2m = least sqares mean
  if (!is.null(ls.index[[kk]])) {
    index.debut<-ls.index[[kk]][1] #à définir selon les bin
    index.fin<-ls.index[[kk]][2] #à définir selon les bin. Ne pas soustraire 1 ici, le faire directement dans l'extraction de diff.pts.LiDAR[a:b-1]
    
    pts.LiDAR.subset<-pts.LiDAR.sourceID[index.debut:index.fin] #sous-ensemble comprenant uniquement les points de la boîte
    diff.pts.LiDAR<-diff(pts.LiDAR.subset[["gpstime"]]) #calcule la différence entre chaque paire d'éléments consécutifs
    pos.impulsion<-c(1, which(diff.pts.LiDAR!=0)+1) #trouve la position du début de chaque impulsion
    
    diff.impulsion<-diff(pos.impulsion)
    index.diff.multiple<-which(diff.impulsion>=2)
    
    if (length(index.diff.multiple)>=nbpairs) {
      index.multiple.debut<-pos.impulsion[index.diff.multiple]
      index.multiple.fin<-index.multiple.debut+(diff.impulsion[index.diff.multiple]-1)
      
      dt.debut<-pts.LiDAR.subset[index.multiple.debut, c("X","Y","Z")]
      dt.fin<-pts.LiDAR.subset[index.multiple.fin, c("X","Y","Z")]
      
      mat.segment.XYZt<-as.matrix(cbind(dt.fin, dt.debut, pts.LiDAR.subset[index.multiple.debut, "gpstime"])) #temporairement, uniquement le gpstime de la première impulsion est copié partout, mais cela pourrait être perfectionné un peu. Le défi est de trouver le temps qui sera le plus représentatif de l'endroit où se trouvera l'avion. Un temps médian de tous les points pourrait être utilisé. Je doute toutefois que cela ait un gros impact peu importe l'utilisation qui sera faite ensuite de la ligne de vol étant donné qu'elle sera déjà approximative en XYZ
      
      PA<-mat.segment.XYZt[,1:3] #XYZ des points de départ
      PB<-mat.segment.XYZt[,4:6] #XYZ des points d'arrivée
      return(cbind(fn_XYZ.l2m(PA, PB), mat.segment.XYZt[1,7], nrow(mat.segment.XYZt)))
    }
  }
}


fn_XYZ.l2m<-function(PA, PB) {
  ### Calcul du point d'intersection XYZ par la technique des moindres carrés.
  ### PA: matrice n x 3 contenant les coordonnées XYZ de chaque point de départ
  ### PB: matrice n x 3 contenant les coordonnées XYZ de chaque point de d'arrivée
  ### Code traduit de MATLAB en R
  ### http://www.mathworks.com/matlabcentral/fileexchange/37192-intersection-point-of-lines-in-3d-space
  
  Si<-PB-PA #N lines described as vectors, direction vector
  ni<-Si/(matrix(sqrt(.rowSums(Si^2, dim(Si)[1], dim(Si)[2]))) %*% matrix(1, nrow=1, ncol=3)) #Normalize vectors
  
  nx<-ni[,1, drop=FALSE]
  ny<-ni[,2, drop=FALSE]
  nz<-ni[,3, drop=FALSE]
  
  SXX<-sum(nx^2-1)
  SYY<-sum(ny^2-1)
  SZZ<-sum(nz^2-1)
  SXY<-sum(nx*ny)
  SXZ<-sum(nx*nz)
  SYZ<-sum(ny*nz)
  
  S<-matrix(c(SXX,SXY,SXZ,SXY,SYY,SYZ,SXZ,SYZ,SZZ), nrow=3, byrow=TRUE)
  
  CX<-sum(PA[,1, drop=FALSE]*(nx^2-1) + PA[,2, drop=FALSE]*(nx*ny)  + PA[,3, drop=FALSE]*(nx*nz))
  CY<-sum(PA[,1, drop=FALSE]*(nx*ny) + PA[,2, drop=FALSE]*(ny^2-1)  + PA[,3, drop=FALSE]*(ny*nz))
  CZ<-sum(PA[,1, drop=FALSE]*(nx*nz) + PA[,2, drop=FALSE]*(ny*nz)  + PA[,3, drop=FALSE]*(nz^2-1))
  
  C<-matrix(c(CX,CY,CZ))
  P_intersect<-Conj(t(solve(S,C))) #Best intersection point of the N lines, in least squares sense
  
  
  ## À RETIRER SI NE SERT PAS À UNE QUELCONQUE VALIDATION
  # distances perpendiculaires entre le point et les lignes
  dist3d <- function(ii, a, b, c) {
    v1 <- b[ii,] - c[ii,]
    v2 <- a[1,] - b[ii,]      
    v3 <- cross3d_prod(v1,v2)
    area <- sqrt(sum(v3*v3))/2
    d <- 2*area/sqrt(sum(v1*v1))
    return(d)
  }
  
  cross3d_prod <- function(v1, v2){
    v3 <- vector()
    v3[1] <- v1[2]*v2[3]-v1[3]*v2[2]
    v3[2] <- v1[3]*v2[1]-v1[1]*v2[3]
    v3[3] <- v1[1]*v2[2]-v1[2]*v2[1]
    return(v3)
  }
  
  N<-dim(PA)[1]
  distances<-sapply(1:N, dist3d, P_intersect, PA, PB) #Distances from intersection point to the input lines
  # func_ECDF<-ecdf(distances)
  # 
  # eval.max<-max(distances)
  # eval.min<-min(distances)
  # limite<-0.8
  # 
  # while (abs(mean(c(eval.max, eval.min))-eval)>0.001) {
  #   eval<-mean(c(eval.max, eval.min))
  #   if(func_ECDF(eval)<limite) {
  #     eval.min<-eval
  #   } else {
  #     eval.max<-eval
  #   }
  # }
  # 
  # plot(func_ECDF, main="Distribution cumulative des distances au point", ylab="Proportion", xlab="Distance (m)")
  # abline(h=limite, lty=2, col="red")
  # abline(v=eval, lwd=2, col="red")
  # text(1, paste0(sprintf("%.1f", eval)," m @ ", limite), pos=1, col="red", font=2)
  ## FIN de À RETIRER SI NE SERT PAS À UNE QUELCONQUE VALIDATION
  
  return(cbind(P_intersect, mean(distances), sd(distances)))
}