library(lidR)

file = "~/Téléchargements/Aircraft position JF/L006-1_2017_HARV_4_V01_2017081715_P01_r.laz"

source("fn_trajectory.R")

ctg = catalog(file)
plot(ctg)

las = readLAS(file, select = "xyzrtp", filter= "-drop_single -drop_y_above 4716800")

output = fn_trajectory(las@data, PtSourceID = NULL, bin = 0.001, step = 2, nbpairs = 20)

sp::plot(output)

report = profvis::profvis({fn_trajectory(las@data, PtSourceID = NULL, bin = 0.001, step = 2, nbpairs = 20)})

# v0.1 : time ~ 18 s | memory used ~ 3.4 GB
# v0.2 : time ~  6 s | memory used ~ 1.3 GB
