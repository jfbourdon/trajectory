library(lidR)
library(magrittr)

source("sensor_tracking.R")
source("lasrangecorrection.R")


#ctg = catalog("~/Documents/ALS data/Montmorency dataset/las/")
#ctg = ctg[c(380,381,382,406,407, 408, 432, 433, 434, 459,460),]
#plot(ctg)

ctg = catalog("~/Documents/ALS data/BCTS/")
ctg = ctg[c(14,16,25,26),]
plot(ctg, chunk = TRUE)

# Test on a point cloud

las <- readLAS(ctg$filename[4], filter = "-drop_single -thin_pulses_with_time 0.0001")
flightlines = sensor_tracking(las, interval = 0.25, pmin = 200, extra_check = FALSE)
plot(extent(las))
plot(flightlines, add = TRUE)

las <- readLAS(ctg$filename[4], select = "ti")

las <- lasrangecorrection(las, flightlines, Rs = 800)
M   <- grid_metrics(las, list(Imax = max(RawIntensity), cImax = max(Intensity)), 5)

I = M[[1]]
cI = M[[2]]

plot(I, col = heat.colors(50))
plot(cI, col = heat.colors(50))

I[I > 170] <- NA
cI[cI > 100] <- NA

plot(I, col = heat.colors(50))
plot(cI, col = heat.colors(50))

# Apply on a catalog

opt_filter(ctg) <- "-drop_single -thin_pulses_with_time 0.0001"

flightlines = sensor_tracking(ctg, interval = 0.25, pmin = 200, extra_check = FALSE)
shapefile(flightlines, filename = "~/Téléchargements/flightlines.shp")

flightlines <- shapefile("~/Téléchargements/flightlines.shp")

plot(ctg)
plot(flightlines, add = TRUE)

#plot(las) %>% add_treetops3d(flightline, radius = 10)

corrected_metrics = function(chunk, flightlines)
{
  las = readLAS(chunk)
  if (is.empty(las)) return(NULL)

  las <- lasrangecorrection(las, flightlines, Rs = 800)
  M   <- grid_metrics(las, list(Imax = max(RawIntensity), cImax = max(Intensity)), 5)
  M   <- raster::crop(M, raster::extent(chunk))
  return(M)
}

opt_filter(ctg) <- ""
opt_chunk_size(ctg) <- 800
opt_chunk_buffer(ctg) <- 1
opt_select(ctg) <- "ti"
options <- list(raster_alignment = 5)
M <- catalog_apply(ctg, corrected_metrics, flightlines = flightlines, .options = options)
M <- lidR:::merge_rasters(M)

I = M[[1]]
cI = M[[2]]

plot(I)
plot(cI)

I[I > 200] <- NA
cI[cI > 300] <- 300

plot(I, col = heat.colors(50))
plot(cI, col = heat.colors(50))

