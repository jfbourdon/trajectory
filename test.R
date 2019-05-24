library(lidR)
library(magrittr)

ctg = catalog("~/Documents/ALS data/Montmorency dataset/las/")
ctg = ctg[c(380,381,382,406,407, 408, 432, 433, 434, 459,460),]
plot(ctg)

source("trajectory.R")

las = readLAS(ctg$filename, select = "xyzrntp", filter = "-keep_point_source 8 -drop_single")
plot(las)

output = sensor_tracking(las, bin = 0.5)
Z <- output@coords[,3]
output@coords = output@coords[,-3]
output$Z <- Z

plot(las) %>% add_treetops3d(output, radius = 10)

report = profvis::profvis({fn_trajectory(las@data, PtSourceID = NULL, bin = 0.001, step = 2, nbpairs = 20)})

microbenchmark::microbenchmark(fn_trajectory(las@data, PtSourceID = NULL, bin = 0.001, step = 2, nbpairs = 20), times = 5)

# v0.1 : time ~  18 s | memory used ~ 3.4 GB
# v0.2 : time ~   6 s | memory used ~ 1.3 GB
# v0.3 : time ~ 1.5 s | memory used ~ 1.1 GB
# v0.4 : time ~   1 s | memory used ~ 1.0 GB
# v0.5 : time ~ 0.5 s | memory used ~ 150 MB
# v0.6 : time ~ 0.4 s | memory used ~ 150 MB