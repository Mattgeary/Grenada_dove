setwd("C:/Users/Mgeary/Desktop/Future climate data/max_80_50_bio")
season <- raster("bc85bi5015.tif")
setwd("C:/Users/Mgeary/Desktop/Future climate data/max_80_50_ppt")
ppt.aug <- raster("bc85pr508.tif")
setwd("C:/Users/Mgeary/Desktop/Future climate data/max_80_50_tmin")
tmin.jan <- raster("bc85tn501.tif")
tmin.aug <- raster("bc85tn508.tif")

setwd("C:/Users/Mgeary/Dropbox/Research/Grenada Dove/Dummy model/Climate")
max_80 <- stack(tmin.jan, season, ppt.aug, tmin.aug)
max_80 <- crop(max_80, extent(GRD))
fire <- (max_80[[4]]/1000) * 4
max_80 <- stack(max_80, fire)
names(max_80) <- var.names
writeRaster(max_80, "max_80.grd")