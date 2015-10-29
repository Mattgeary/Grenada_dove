# list.of.packages <- c("raster", "sp", "dismo", "MuMIn", "rgdal", "rms")
# new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# if(length(new.packages)) install.packages(new.packages)

library(rgdal)
library(raster)
library(sp)
library(dismo)
# library(MuMIn)
# library(rms)

options("na.action" = na.fail) 

work <- "C:/Users/mgeary/Dropbox/Research/Grenada Dove"


### Load environmental and species location data

# Base map

	GRD <- getData('GADM', country='GRD', level=0)

# Restrict to island of Grenada

	GRD <- SpatialPolygons(list(Polygons(list(GRD@polygons[[1]]@Polygons[[8]]), ID="1")))

# Altitude

	#alt <- getData('alt', country='GRD')

	#alt <- crop(alt, GRD)
	
	setwd(paste(work, "/Climate_risk_model/dem", sep = ""))
		alt <- raster('GRD_alt.asc')
	dem <- raster("DEM.img")
	dem <- projectRaster(dem, crs = CRS("+init=epsg:4326"))
	# dem.r <- resample(dem, alt)
	# rm(dem)
	# writeRaster(dem.r, "dem_resample.asc")
	dem.r <- raster("dem_resample.asc")
  slope <- raster("slope.asc")
  aspect <- raster("aspect.asc")
	setwd(paste(work, "/Climate_risk_model/", sep = ""))



## Precipitation

#	#ppt <- getData('worldclim', var='prec', res=0.5, lat = 12.05, lon = -61.75)

#	#ppt <- crop(ppt, extent(GRD))

#	#ppt.jan <- subset(ppt, 1)

#	#ppt.aug <- subset(ppt, 2)

#	ppt.jan <- raster('GRD_ppt_jan.asc')

#	ppt.aug <- raster('GRD_ppt_aug.asc')

## Temperature

#	#tmin <- getData('worldclim', var='tmin', res=0.5, lat = 12.05, lon = -61.75)

#	#tmin <- crop(tmin, extent(GRD))

#	tmin.jan <- raster('GRD_tmin_jan.asc')

#	tmin.aug <- raster('GRD_tmin_aug.asc')

## BIOCLIM variables

#	#BIO <- getData('worldclim', var='bio', res=0.5, lat = 12.05, lon = -61.75)

#	#BIO <- crop(BIO, extent(GRD))

# ppt.season <- subset(BIO, 15)
# 
# ppt.dry.month <- subset(BIO, 14)
# 
# ppt.dry.qrt <- subset(BIO, 17)
# 
# m.temp.dry.qrt <- subset(BIO, 9)

setwd(paste(work, "/Climate_risk_model/Current_clim", sep = ""))

ppt.season <- raster("ppt.season.asc")

ppt.dry.month <- raster("ppt.dry.month.asc")

ppt.dry.qrt <- raster("ppt.dry.qrt.asc")

m.temp.dry.qrt <- raster("m.temp.dry.qrt.asc")

PET <- raster("perrbfsimf1km.img")

PET <- crop(PET, extent(GRD))

clip<-function(raster,shape) {
  a1_crop<-crop(raster,shape)
  step1<-rasterize(shape,a1_crop)
  a1_crop*step1}

PET <- clip(PET, GRD)

layer.names <- c("firerisk1", "pdqrbfsimf1km", "perrbfsimf1km", "taprerbfsim1k", "tbirbfsimf1km", "temrbfsimf1km")
fire.risk <- raster("firerisk1.img")
total.ppt.dr.qrt <- raster("pdqrbfsimf1km.img")
total.ppt.annual <- raster("taprerbfsim1k.img")
biotemp <- raster("tbirbfsimf1km.img")
temp <- raster("temrbfsimf1km.img")

### Future climate

# Select year - 2050 or 2070

setwd(paste(work, "/Climate_risk_model/Future/Fifty", sep = ""))
future.ppt.dry.month.50 <- stack(raster("ppt.dry.month.26.asc"), raster("ppt.dry.month.45.asc"), raster("ppt.dry.month.60.asc"), raster("ppt.dry.month.85.asc"))

setwd(paste(work, "/Climate_risk_model/Future/Seventy", sep = ""))
future.ppt.dry.month.70 <- stack(raster("ppt.dry.month.26.asc"), raster("ppt.dry.month.45.asc"), raster("ppt.dry.month.60.asc"), raster("ppt.dry.month.85.asc"))

### Select point in future
future.ppt.dry.month <- future.ppt.dry.month.50


#	#temp.range <- subset(BIO, 7)

#	ppt.season <- raster('GRD_ppt_season.asc')

#	temp.range <- raster('GRD_temp_range.asc')

### Soil and Geology
# ext <-  extent (alt)
# xy <- abs(apply(as.matrix(bbox(ext)), 1, diff))
# n <- 5
# r <- raster(ext, ncol=xy[1]*5, nrow=xy[2]*5)
# 
# # Soil 
 setwd(paste(work, "/Climate_risk_model/Land_use", sep = ""))
# 
# soil <-  readOGR(dsn = paste(work, "/Climate_risk_model/Land_use", layer  = "Soils")
# soil <- spTransform(soil, CRS = CRS("+init=epsg:4326"))
# 
#                 ## Rasterize the shapefile
# soil.raster.code <-rasterize(soil, alt, field = "CODE_DESC", background = NA)
# writeRaster(soil.raster.code, "soil.raster.code.asc")
# soil.raster.desc <-rasterize(soil, alt, field = "S_MAP_DESC", background = NA)
# writeRaster(soil.raster.desc, "soil.raster.desc.asc")
soil.raster.code <- raster("soil.raster.code.asc")
soil.raster.desc <- raster("soil.raster.desc.asc")

#  Geology

# geol <-  readOGR(dsn = paste(work, "/Climate_risk_model/Land_use", layer  = "Geology")
# geol <- spTransform(geol, CRS = CRS("+init=epsg:4326"))
# 
#                 ## Rasterize the shapefile
# geol.raster.code <-rasterize(geol, alt, field = "GEO_CODE", background = NA)
# writeRaster(geol.raster.code, "geol.raster.code.asc")
# geol.raster.desc <-rasterize(geol, alt, field = "GEO_DESC", background = NA)
# writeRaster(geol.raster.desc, "geol.raster.desc.asc")
geol.raster.code <- raster("geol.raster.code.asc")
geol.raster.desc <- raster("geol.raster.desc.asc")
setwd(paste(work, "/Climate_risk_model/", sep = ""))


#### Fire risk raster

setwd(paste(work, "/Climate_risk_model/Current_clim", sep = ""))

#fire <- raster("weightedall12.tif")
#fire <- projectRaster(fire, crs = CRS("+init=epsg:4326"))
#fire <- resample(fire, alt)
#writeRaster(fire, "fire.asc")
fire <- raster("fire.asc")
setwd(paste(work, "/Climate_risk_model/", sep = ""))

# Load current climate stack

	# current <- stack("Climate/current.grd")

# Create landcover raster mask

setwd(paste(work, "/Climate_risk_model/Land_use", sep = ""))
# land.use <-  readOGR(dsn = paste(work, "/Climate_risk_model/Land_use", layer  = "Landuse2009")
# land.use <- spTransform(land.use, CRS = CRS("+init=epsg:4326"))
# land.use.levels <- data.frame(land.use = c(levels(land.use$LUDESC)), value = c(1:8))
#                 ## Rasterize the shapefile
# land.use <-rasterize(land.use, alt, field = "LUDESC", background = NA)
# land.rcl <- matrix(c(c(1:8), c(0,1,1,1,0,1,1,0)), ncol = 2)
# land.use.mask <- reclassify(land.use, land.rcl, right = NA)
# writeRaster(land.use, "land.use.asc")
# writeRaster(land.use.mask, "land.use.mask.asc")
land.use <- raster("land.use.asc")
land.use.mask <- raster("land.use.mask.asc")
setwd(paste(work, "/Climate_risk_model/", sep = ""))
# Sea level rise

#	sea.low <- reclassify(alt, matrix(c(-354, 0.6, NA), nrow=3, ncol=3, byrow=T))
#	sea.high <- reclassify(alt, matrix(c(-354, 1.2, NA), nrow=3, ncol=3, byrow=T))

sea.lvl <- reclassify(dem, matrix(c(-354, 6, NA, 6, 3000, 0), nrow=2, ncol=3, byrow=T)) # Change to hi res DEM
sea.lvl <- resample(sea.lvl, alt)

# Dove locations

	species.pa <- read.csv("indicatorPresAbs.csv")
	forest.pa <- data.frame(long = species.pa$long, lat = species.pa$lat, pa = species.pa$ALL)
	
	train.rows <- sample(c(1:length(forest.pa$pa)), 37)
	train.forest <- forest.pa[train.rows,]
	test.forest <- forest.pa[-train.rows,]

### Create environmental data frame for modelling

#	forest.env <- data.frame(pres = dove$pres, alt = extract(alt, dove[,1:2]), tmin.jan = extract(tmin.jan, dove[,1:2]), tmin.aug = extract(tmin.aug, dove[,1:2]), temp.range = extract(temp.range, dove[,1:2]), ppt.jan = extract(ppt.jan, dove[,1:2]), ppt.aug = extract(ppt.aug, dove[,1:2]), ppt.season = extract(ppt.season, dove[,1:2]), fire = extract(fire, dove[,1:2]))

forest.env <- data.frame(lat = train.forest$lat, long = train.forest$long, pres = train.forest$pa, dem.r = extract(dem.r, train.forest[,1:2]), ppt.dry.month  = extract(ppt.dry.month, train.forest[,1:2]), ppt.dry.qrt = extract(ppt.dry.qrt, train.forest[,1:2]), ppt.season = extract(ppt.season, train.forest[,1:2]), m.temp.dry.qrt = extract(m.temp.dry.qrt, train.forest[,1:2]), soil.raster.desc = extract(soil.raster.desc, train.forest[,1:2]), geol.raster.code = extract(soil.raster.code, train.forest[,1:2]), fire = extract(fire, train.forest[,1:2]), slope = extract(slope, train.forest[,1:2]), aspect = extract(aspect, train.forest[,1:2]))
# forest.env <- data.frame(lat = train.forest$lat, long = train.forest$long, pres = train.forest$pa, PET = extract(PET, train.forest[,1:2]), fire.risk = extract(fire.risk, train.forest[,1:2]), total.ppt.dr.qrt = extract(total.ppt.dr.qrt, train.forest[,1:2]), total.ppt.annual = extract(total.ppt.annual, train.forest[,1:2]), biotemp = extract(biotemp, train.forest[,1:2]), temp = extract(temp, train.forest[,1:2]))

test.env <- data.frame(lat = test.forest$lat, long = test.forest$long, pres = test.forest$pa, dem.r = extract(dem.r, test.forest[,1:2]), ppt.dry.month  = extract(ppt.dry.month, test.forest[,1:2]), ppt.dry.qrt = extract(ppt.dry.qrt, test.forest[,1:2]), ppt.season = extract(ppt.season, test.forest[,1:2]), m.temp.dry.qrt = extract(m.temp.dry.qrt, test.forest[,1:2]), soil.raster.desc = extract(soil.raster.desc, test.forest[,1:2]), geol.raster.code = extract(soil.raster.code, test.forest[,1:2]), fire = extract(fire, test.forest[,1:2]), slope = extract(slope, test.forest[,1:2]), aspect = extract(aspect, test.forest[,1:2]))
# test.env <- data.frame(lat = test.forest$lat, long = test.forest$long, pres = test.forest$pa, PET = extract(PET, test.forest[,1:2]), fire.risk = extract(fire.risk, test.forest[,1:2]), total.ppt.dr.qrt = extract(total.ppt.dr.qrt, test.forest[,1:2]), total.ppt.annual = extract(total.ppt.annual, test.forest[,1:2]), biotemp = extract(biotemp, test.forest[,1:2]), temp = extract(temp, test.forest[,1:2]))

future.env <- 	data.frame(ppt.dry.month = as.data.frame(future.ppt.dry.month), fire = as.data.frame(fire), aspect = as.data.frame(aspect))

# Test for correlated variables

	cor(forest.env) 

panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, use = "pairwise.complete.obs"))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
  
  text(0.5, 0.5, txt, cex = cex * r)   
}

pairs(forest.env, lower.panel = panel.smooth, upper.panel = panel.cor)

# Define global model with all candidate predictors

#	 forest.mod <- glm(pres ~ dem.r + ppt.dry.month + ppt.dry.qrt + ppt.season + m.temp.dry.qrt + soil.raster.code + geol.raster.code + fire, family = binomial(link='logit'), data=forest.env)

forest.mod <- glm(pres ~ ppt.dry.month + fire + aspect, family = binomial(link='logit'), data=forest.env)
forest.mod.2 <- glm(pres ~ ppt.dry.month + fire, family = binomial(link='logit'), data=forest.env)
forest.mod.3 <- glm(pres ~ ppt.dry.month + aspect, family = binomial(link='logit'), data=forest.env)
# forest.soil <- glm(pres ~ ppt.dry.month + fire + aspect + soil, family = binomial(link='logit'), data=forest.env)

forest.AICc <- data.frame(mod = c("Full", "Fire", "Aspect"), AICc = c(AICc(forest.mod), AICc(forest.mod.2), AICc(forest.mod.3)))
forest.AICc$delta <- forest.AICc$AICc - min(forest.AICc$AICc)
#forest.lrm <- lrm(pres ~ ppt.dry.month + fire + aspect, data=forest.env)
# Fit all candidate models


	#models.dove <- dredge(forest.mod)

# Average models with delta < 4

#	dove.avg <- model.avg(subset(models.dove, delta < 4), fit=TRUE)

	forest.eval <- evaluate(test.env[test.env$pres == 1,], test.env[test.env$pres == 0,], forest.mod)
  forest.eval.2 <- evaluate(test.env[test.env$pres == 1,], test.env[test.env$pres == 0,], forest.mod.2)

# Predictions

#	env.vars <- data.frame(dem.r = as.data.frame(as(dem.r, "SpatialGridDataFrame"))[,1], ppt.season = as.data.frame(as(ppt.season, "SpatialGridDataFrame"))[,1], ppt.dry.month = as.data.frame(as(ppt.dry.month, "SpatialGridDataFrame"))[,1], ppt.dry.qrt = as.data.frame(as(ppt.dry.qrt, "SpatialGridDataFrame"))[,1], m.temp.dry.qrt = as.data.frame(as(m.temp.dry.qrt, "SpatialGridDataFrame"))[,1], soil.raster.code = as.data.frame(as(soil.raster.code, "SpatialGridDataFrame"))[,1], geol.raster.code = as.data.frame(as(geol.raster.code, "SpatialGridDataFrame"))[,1], fire = as.data.frame(as(fire, "SpatialGridDataFrame"))[,1])
env.vars <- data.frame(ppt.dry.month = as.data.frame(ppt.dry.month), fire = as.data.frame(fire), aspect = as.data.frame(aspect))

	pred <- predict(forest.mod, newdata = env.vars, type="response")
 pred.2 <- predict(forest.mod.2, newdata = env.vars, type="response")

#	pred.df <- data.frame(data = pred, x = as.data.frame(as(fire, "SpatialGridDataFrame"))[,2], y = as.data.frame(as(fire, "SpatialGridDataFrame"))[,3])
 pred.df <- data.frame(data = pred, x = coordinates(fire)[,1], y = coordinates(fire)[,2])
 pred.df.2 <- data.frame(data = pred.2, x = coordinates(fire)[,1], y = coordinates(fire)[,2])

	coordinates(pred.df) <- ~x+y
  coordinates(pred.df.2) <- ~x+y
	
  gridded(pred.df) <- TRUE
  gridded(pred.df.2) <- TRUE

	pred.sp <- raster(pred.df, layer='data')
	pred.sp <- crop(pred.sp, extent(GRD))

  pred.sp.2 <- raster(pred.df.2, layer='data')
  pred.sp.2 <- crop(pred.sp.2, extent(GRD))

# Mask by land cover

	pred.sp.mask <- pred.sp - land.use.mask
	rcl <- matrix(c(-1.5,0,0), nrow=1, ncol=3, byrow=T)
	pred.sp.mask <- reclassify(pred.sp.mask, rcl)

#####################################

## Model simulation functions ##

## Functions to 'grow' patches

# Probability of transition
f.prob <- function(x, p.trans){
          x*p.trans
 }

 # Selects from a binary distribution to determine transition
f.bin <- function(x){
    if(!is.na(x)) x <- rbinom(1,1,x)
    else x <- NA
 }

 # Corrects for probabilities > 1
f.correct <- function(x){
  if(!is.na(x) && x > 1) x <- 1
  else x
}

f.Exp <- function(map, base.map, p.trans){
              mat.0 <- as.matrix(map)
               map.t <- focal(map, w=matrix(1, nrow=3, ncol=3), sum, na.rm=T, pad=T) 
              map.t <- map.t + base.map
                mat.t <- as.matrix(map.t) 
                mat.prob <- apply(mat.t, c(1,2), f.prob, p.trans=p.trans)
               mat.cor <- mat.0 + mat.prob
              mat.cor <- apply(mat.cor, c(1,2), f.correct)
               mat.new <- apply(mat.cor, c(1,2), f.bin)
               map.new <- raster(mat.new, xmn=xmin(map), xmx=xmax(map), ymn=ymin(map), ymx=ymax(map)) # need to add CRS to thins function if we want one
}


## Fire risk increase function

# fire.risk <- function(fire, p.trans=0.05, iterations=50,  fires = 10, threshold = 0.7, total = 40){
	# rcl <- matrix(c(-1,threshold,NA,(threshold + 0.01),1.2,0), nrow=2, ncol=3, byrow=T)
	# potential <- reclassify(fire, rcl)

	# rnd.pts <- randomPoints(potential, fires)
	# rnd.pts <- as.data.frame(rnd.pts)
	# potential.pts <- rasterize(rnd.pts, potential, field = 1, background=0)
  # potential.pts
	# potential.pts <- potential + potential.pts

	# wildfire <- f.Exp(potential.pts, potential.pts, p.trans)
	# for(i in 1:iterations){
		# if(cellStats(wildfire, sum) <= total) wildfire <- f.Exp(wildfire, potential.pts, p.trans)
		# else wildfire <- wildfire
	    # }

	# return(wildfire)
# }

### Hurricane function

hurricane <- function(template, p.trans=0.05, iterations=50, years = 35,  strike = 0.073, threshold = 0.7, total = rpois(300, c(2,4,10))){
	total[total > 19] <- 19
	canes <- strike * years
	h.map <- template
	projection(h.map) <- CRS("+init=epsg:4326")
  h.map <- reclassify(h.map, matrix(c(0, 700, 0), nrow = 1))
	potential <- h.map
	for(j in 1:canes){
		rnd.pts <- randomPoints(potential, 1)
		rnd.pts <- as.data.frame(rnd.pts)
		potential.pts <- rasterize(rnd.pts, potential, field = 1, background=0)
		potential.pts
		potential.pts <- potential + potential.pts

		extreme <- f.Exp(potential.pts, potential.pts, p.trans)
		for(i in 1:iterations){
			if(cellStats(extreme, sum) <= sample(total, 1)) extreme <- f.Exp(extreme, potential.pts, p.trans)
			else extreme <- extreme
	    }
		h.map <- h.map + extreme
	}
	h.map[h.map > 1] <- 1
	return(h.map)
}


## Development increase function

# development <- function(cov, p.trans=0.1, iterations=50,  hotels = 10, total = 40){
	# rnd.pts <- randomPoints(cov, hotels)
	# rnd.pts <- as.data.frame(rnd.pts)
  # potential.pts <- rasterize(rnd.pts, cov, field = 1, background=0)
	# potential.pts <- cov + potential.pts

	# develop <- f.Exp(potential.pts, cov, p.trans)
	# for(i in 1:iterations){
		# if(cellStats(develop, sum) <= total) develop <- f.Exp(develop, potential.pts, p.trans)
		# else develop <- develop
	    # }

	# return(develop)
# }

###### Future scenarios function #####

climate.fun <- function(model, future.ppt){

	env.vars <- data.frame("ppt.dry.month" = as.data.frame(future.ppt), "fire" = as.data.frame(fire), "aspect" = as.data.frame(aspect))
  names(env.vars) <- c("ppt.dry.month", "fire", "aspect")
	
  pred <- predict(model, newdata = env.vars, type="response")

	pred.df <- data.frame(data = pred, x = coordinates(future.ppt)[,1], y = coordinates(future.ppt)[,2])

	coordinates(pred.df) <- ~x+y

	gridded(pred.df) <- TRUE

	pred.sp <- raster(pred.df, layer='data')
	
	return(pred.sp)
}


### Need to supply name of model from above - default predicts under all scenarios

# Threshold is maximum sensitivity and specificity set from threshold()

future <- function(model, years = 35, sea = TRUE, storm = TRUE, develop = TRUE, climate = TRUE, all = TRUE, reps=10 , threshold = 0.52){

	if(years == 35) future.ppt.dry.month <- future.ppt.dry.month.50
	
	if(years == 55) future.ppt.dry.month <- future.ppt.dry.month.70

# Initial area predicted present

	pred.current <- climate.fun(model, ppt.dry.month)

	thresh.rcl <- matrix(c(-1, threshold, 0, threshold, 1.5, 1), nrow=2, ncol=3, byrow=TRUE)

	current.area <- cellStats(reclassify(pred.current, thresh.rcl), mean)

# Create list to contain results

	pred.area <- list(c(reps = reps, threshold = threshold, pred.area = current.area))

# Conditional statements to add elements to the list

	if(climate == TRUE) pred.area$climate <- data.frame(low_26 = numeric(1), mid_45 = numeric(1), high_60 = numeric(1), max_80 = numeric(1))

	if(sea == TRUE) pred.area$sea <- data.frame(sea = numeric(1), sea.low_26 = numeric(1), sea.mid_45 = numeric(1), sea.high_60 = numeric(1), sea.max_80 = numeric(1))

	if(develop == TRUE) pred.area$develop <- data.frame(develop = numeric(reps), dev.low_26 = numeric(reps), dev.mid_45 = numeric(reps), dev.high_60 = numeric(reps), dev.max_80 = numeric(reps))

	if(storm == TRUE) pred.area$hurricane <- data.frame(hurricane = numeric(reps), hurricane.low_26 = numeric(reps), hurricane.mid_45 = numeric(reps), hurricane.high_60 = numeric(reps), hurricane.max_80 = numeric(reps))

	if(all == TRUE) pred.area$all <- data.frame(hurricane = numeric(reps), sea.dev = numeric(reps), dev.hurricane = numeric(reps), sea.dev.hurricane = numeric(reps), sea.hurricane = numeric(reps), all.low_26 = numeric(reps), all.mid_45 = numeric(reps), all.high_60 = numeric(reps), all.max_80 = numeric(reps))

#Calculate predicted areas under scenarios

#Climate - single values
	if(climate == TRUE){

	# Predict suitability across area
		low_26 <- climate.fun(model, future.ppt.dry.month[[1]])
		mid_45 <- climate.fun(model, future.ppt.dry.month[[2]])
		high_60 <- climate.fun(model, future.ppt.dry.month[[3]])
		max_80 <- climate.fun(model, future.ppt.dry.month[[4]])

	# Calculate proportion of island predicted
		pred.area$climate$low_26 <- cellStats(reclassify(low_26, thresh.rcl), mean)
		pred.area$climate$mid_45 <- cellStats(reclassify(mid_45, thresh.rcl), mean)
		pred.area$climate$high_60 <- cellStats(reclassify(high_60, thresh.rcl), mean)
		pred.area$climate$max_80 <- cellStats(reclassify(max_80, thresh.rcl), mean)

		rm(low_26, mid_45, high_60, max_80)
	}

# Sea level change
	if(sea == TRUE){
		
	# Mask initial predictions by sea level change
		sea.mask <- pred.current + sea.lvl
		pred.area$sea$sea <- cellStats(reclassify(sea.mask, thresh.rcl), mean)

	# Create future scenarios and mask by sea level
		sea.map <- climate.fun(model, future.ppt.dry.month[[1]])
			sea.mask <- sea.map + sea.lvl
			pred.area$sea$sea.low_26 <- cellStats(reclassify(sea.mask, thresh.rcl), mean)
		sea.map <- climate.fun(model, future.ppt.dry.month[[2]])
			sea.mask <- sea.map + sea.lvl
			pred.area$sea$sea.mid_45 <- cellStats(reclassify(sea.mask, thresh.rcl), mean)
		sea.map <- climate.fun(model, future.ppt.dry.month[[3]])
			sea.mask <- sea.map + sea.lvl
			pred.area$sea$sea.high_60 <- cellStats(reclassify(sea.mask, thresh.rcl), mean)
		sea.map <- climate.fun(model, future.ppt.dry.month[[4]])
			sea.mask <- sea.map + sea.lvl
			pred.area$sea$sea.max_80 <- cellStats(reclassify(sea.mask, thresh.rcl), mean)
	
		rm(sea.map, sea.mask)
	}

# Development
	if(develop == TRUE){
		low_26 <- climate.fun(model, future.ppt.dry.month[[1]])
		mid_45 <- climate.fun(model, future.ppt.dry.month[[2]])
		high_60 <- climate.fun(model, future.ppt.dry.month[[3]])
		max_80 <- climate.fun(model, future.ppt.dry.month[[4]])

		for(i in 1:reps){
		# Simulate development with default values
		setwd(paste(work, "/Climate_risk_model/Land_use", sep = ""))
			dev.mask <-  raster("land.use.mask.asc")
		setwd(paste(work, "/Climate_risk_model/", sep = ""))
		
		# Use simulated development to mask climate predictions
			pred.area$develop$develop[i] <- cellStats(reclassify(reclassify(pred.current, thresh.rcl) - dev.mask, matrix(c(-1.5,0,0), nrow=1, ncol=3, byrow=T) ), mean)
			pred.area$develop$dev.low_26[i] <- cellStats(reclassify(reclassify(low_26, thresh.rcl) - dev.mask, matrix(c(-1.5,0,0), nrow=1, ncol=3, byrow=T)), mean)
			pred.area$develop$dev.mid_45[i] <- cellStats(reclassify(reclassify(mid_45, thresh.rcl) - dev.mask, matrix(c(-1.5,0,0), nrow=1, ncol=3, byrow=T)), mean)
			pred.area$develop$dev.high_60[i] <- cellStats(reclassify(reclassify(high_60, thresh.rcl) - dev.mask, matrix(c(-1.5,0,0), nrow=1, ncol=3, byrow=T)), mean)
			pred.area$develop$dev.max_80[i] <- cellStats(reclassify(reclassify(max_80, thresh.rcl) - dev.mask, matrix(c(-1.5,0,0), nrow=1, ncol=3, byrow=T)), mean)
		}
		rm(low_26, mid_45, high_60, max_80)
	}

# Hurricane
	if(storm == TRUE){
		low_26 <- climate.fun(model, future.ppt.dry.month[[1]])
		mid_45 <- climate.fun(model, future.ppt.dry.month[[2]])
		high_60 <- climate.fun(model, future.ppt.dry.month[[3]])
		max_80 <- climate.fun(model, future.ppt.dry.month[[4]])

		for(i in 1:reps){
		# Simulate hurricane with default values
			hurricane.dis.current <- hurricane(alt)			
			hurricane.dis.low_26 <- hurricane(alt)	
			hurricane.dis.mid_45 <-hurricane(alt)	
			hurricane.dis.high_60 <- hurricane(alt)	
			hurricane.dis.max_80 <- hurricane(alt)	
			hurricane.dis.current <- reclassify(hurricane.dis.current, matrix(c(-1,0.99999,0, 0.999999, 1, 1), nrow=2, ncol=3, byrow=T))
			hurricane.dis.low_26 <- reclassify(hurricane.dis.low_26, matrix(c(-1,0.99999,0, 0.999999, 1, 1), nrow=2, ncol=3, byrow=T))	
			hurricane.dis.mid_45 <- reclassify(hurricane.dis.mid_45, matrix(c(-1,0.99999,0, 0.999999, 1, 1), nrow=2, ncol=3, byrow=T))
			hurricane.dis.high_60 <- reclassify(hurricane.dis.high_60, matrix(c(-1,0.99999,0, 0.999999, 1, 1), nrow=2, ncol=3, byrow=T))
			hurricane.dis.max_80 <- reclassify(hurricane.dis.max_80, matrix(c(-1,0.99999,0, 0.999999, 1, 1), nrow=2, ncol=3, byrow=T))
		
		# Use simulated development to mask climate predictions
			pred.area$hurricane$hurricane[i] <- cellStats(reclassify(reclassify(pred.current, thresh.rcl) - hurricane.dis.current, matrix(c(-1.5,0,0), nrow=1, ncol=3, byrow=T)), mean)
			pred.area$hurricane$hurricane.low_26[i] <- cellStats(reclassify(reclassify(low_26, thresh.rcl) - hurricane.dis.low_26, matrix(c(-1.5,0,0), nrow=1, ncol=3, byrow=T)), mean)
			pred.area$hurricane$hurricane.mid_45[i] <- cellStats(reclassify(reclassify(mid_45, thresh.rcl) - hurricane.dis.mid_45, matrix(c(-1.5,0,0), nrow=1, ncol=3, byrow=T)), mean)
			pred.area$hurricane$hurricane.high_60[i] <- cellStats(reclassify(reclassify(high_60, thresh.rcl) - hurricane.dis.high_60, matrix(c(-1.5,0,0), nrow=1, ncol=3, byrow=T)), mean)
			pred.area$hurricane$hurricane.max_80[i] <- cellStats(reclassify(reclassify(max_80, thresh.rcl) - hurricane.dis.max_80, matrix(c(-1.5,0,0), nrow=1, ncol=3, byrow=T)), mean)
		}
		rm(low_26, mid_45, high_60, max_80)
	}

# All
	if(all == TRUE){
		low_26 <- climate.fun(model, future.ppt.dry.month[[1]])
		mid_45 <- climate.fun(model, future.ppt.dry.month[[2]])
		high_60 <- climate.fun(model, future.ppt.dry.month[[3]])
		max_80 <- climate.fun(model, future.ppt.dry.month[[4]])
		
		for(i in 1:reps){	
			# Reclassify prediction
			current.bin <- reclassify(pred.current, thresh.rcl)			
			low_26 <- reclassify(low_26, thresh.rcl)
			mid_45 <- reclassify(mid_45, thresh.rcl)			
			high_60 <- reclassify(high_60, thresh.rcl)
			max_80 <- reclassify(max_80, thresh.rcl)

			# Development mask
			setwd(paste(work, "/Climate_risk_model/Land_use", sep = ""))
			dev.mask <-  raster("land.use.mask.asc")
			setwd(paste(work, "/Climate_risk_model/", sep = ""))
		
			# Simulate hurricane with default values
			hurricane.dis.current <- hurricane(alt)			
			hurricane.dis.low_26 <- hurricane(alt)	
			hurricane.dis.mid_45 <-hurricane(alt)	
			hurricane.dis.high_60 <- hurricane(alt)	
			hurricane.dis.max_80 <- hurricane(alt)	
			hurricane.dis.current <- reclassify(hurricane.dis.current, matrix(c(-1,0.99999,0, 0.999999, 1, 1), nrow=2, ncol=3, byrow=T))
			hurricane.dis.low_26 <- reclassify(hurricane.dis.low_26, matrix(c(-1,0.99999,0, 0.999999, 1, 1), nrow=2, ncol=3, byrow=T))	
			hurricane.dis.mid_45 <- reclassify(hurricane.dis.mid_45, matrix(c(-1,0.99999,0, 0.999999, 1, 1), nrow=2, ncol=3, byrow=T))
			hurricane.dis.high_60 <- reclassify(hurricane.dis.high_60, matrix(c(-1,0.99999,0, 0.999999, 1, 1), nrow=2, ncol=3, byrow=T))
			hurricane.dis.max_80 <- reclassify(hurricane.dis.max_80, matrix(c(-1,0.99999,0, 0.999999, 1, 1), nrow=2, ncol=3, byrow=T))
			
			# Use simulated hurricanes to mask climate predictions
			pred.area$all$hurricane[i] <- cellStats(reclassify(reclassify(pred.current, thresh.rcl) - hurricane.dis.current, matrix(c(-1.5,0,0), nrow=1, ncol=3, byrow=T)), mean)
			pred.area$all$hurricane.low_26[i] <- cellStats(reclassify(reclassify(low_26, thresh.rcl) - hurricane.dis.low_26, matrix(c(-1.5,0,0), nrow=1, ncol=3, byrow=T)), mean)
			pred.area$all$hurricane.mid_45[i] <- cellStats(reclassify(reclassify(mid_45, thresh.rcl) - hurricane.dis.mid_45, matrix(c(-1.5,0,0), nrow=1, ncol=3, byrow=T)), mean)
			pred.area$all$hurricane.high_60[i] <- cellStats(reclassify(reclassify(high_60, thresh.rcl) - hurricane.dis.high_60, matrix(c(-1.5,0,0), nrow=1, ncol=3, byrow=T)), mean)
			pred.area$all$hurricane.max_80[i] <- cellStats(reclassify(reclassify(max_80, thresh.rcl) - hurricane.dis.max_80, matrix(c(-1.5,0,0), nrow=1, ncol=3, byrow=T)), mean)

			# Sea level and development
			mask <- current.bin + sea.lvl
			mask <- reclassify(mask - dev.mask, matrix(c(-1.5,0,0)))
			pred.area$all$sea.dev[i] <- cellStats(mask, mean)

			# Sea level and hurricane
			mask <- current.bin + sea.lvl
			mask <- reclassify(mask - hurricane.dis.current, matrix(c(-1.5,0,0)))
			pred.area$all$sea.hurricane[i] <- cellStats(mask, mean)

			# Fire and development
			mask <- reclassify(current.bin - dev.mask, matrix(c(-1.5,0,0)))
			mask <- reclassify(mask - hurricane.dis.current, matrix(c(-1.5,0,0)))
			pred.area$all$dev.hurricane[i] <- cellStats(mask, mean)

			# Sea level and hurricane
			# Current values			
			mask <- current.bin + sea.lvl
			mask <- reclassify(mask - dev.mask, matrix(c(-1.5,0,0)))
			mask <- reclassify(mask - hurricane.dis.current, matrix(c(-1.5,0,0)))
			pred.area$all$sea.dev.hurricane[i] <- cellStats(mask, mean)
				
			# Low_26
			mask <- low_26 + sea.lvl
			mask <- reclassify(mask - dev.mask, matrix(c(-1.5,0,0)))
			mask <- reclassify(mask - hurricane.dis.low_26, matrix(c(-1.5,0,0)))
			pred.area$all$all.low_26[i] <- cellStats(mask, mean)

			# Mid_45
			mask <- mid_45 + sea.lvl
			mask <- reclassify(mask - dev.mask, matrix(c(-1.5,0,0)))
			mask <- reclassify(mask - hurricane.dis.mid_45, matrix(c(-1.5,0,0)))
			pred.area$all$all.mid_45[i] <- cellStats(mask, mean)

			# High_60
			mask <- high_60 + sea.lvl
			mask <- reclassify(mask - dev.mask, matrix(c(-1.5,0,0)))
			mask <- reclassify(mask - hurricane.dis.high_60, matrix(c(-1.5,0,0)))
			pred.area$all$all.high_60[i] <- cellStats(mask, mean)

			# Max_80
			mask <- max_80 + sea.lvl
			mask <- reclassify(mask - dev.mask, matrix(c(-1.5,0,0)))
			mask <- reclassify(mask - hurricane.dis.max_80, matrix(c(-1.5,0,0)))
			pred.area$all$all.max_80[i] <- cellStats(mask, mean)

		}
		rm(low_26, mid_45, high_60, max_80, mask, dev.mask)
	}

	return(simulation.predictions = pred.area)
}

##########
 
