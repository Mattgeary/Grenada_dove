list.of.packages <- c("raster", "sp", "dismo", "MuMIn")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(raster)
library(sp)
library(dismo)
library(MuMIn)

options("na.action" = na.fail) 


### Load environmental and species location data

# Base map

	GRD <- getData('GADM', country='GRD', level=0)

# Restrict to island of Grenada

	GRD <- SpatialPolygons(list(Polygons(list(GRD@polygons[[1]]@Polygons[[8]]), ID="1")))

# Altitude

	#alt <- getData('alt', country='GRD')

	#alt <- crop(alt, GRD)

	alt <- raster('GRD_alt.asc')

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

#	#ppt.season <- subset(BIO, 15)

#	#temp.range <- subset(BIO, 7)

#	ppt.season <- raster('GRD_ppt_season.asc')

#	temp.range <- raster('GRD_temp_range.asc')

#### Create FAKE fire risk raster

#	fire <- (tmin.aug/1000) * 4

# Load current climate stack

	current <- stack("Climate/current.grd")

# Create FAKE landcover raster

	rcl <- matrix(c(-354, 29.999, 1, 29.9999,  3884,0), nrow=3, ncol=3, byrow=T)
	cov <- reclassify(alt, rcl)

# Sea level rise

#	sea.low <- reclassify(alt, matrix(c(-354, 0.6, NA), nrow=3, ncol=3, byrow=T))
#	sea.high <- reclassify(alt, matrix(c(-354, 1.2, NA), nrow=3, ncol=3, byrow=T))
sea.lvl <- reclassify(alt, matrix(c(-354, 1.2, NA, 1.2, 3000, 0), nrow=2, ncol=3, byrow=T))

# Dove locations

	dove <- randomPoints(alt, 100, ext=extent(GRD))

	dove <- data.frame(x = dove[,1], y=dove[,2], pres = 0)
	
	dove$pres[sample(c(1:length(dove$pres)), 80)] <-1

	train.rows <- sample(c(1:length(dove$pres)), 20)
	train.dove <- dove[train.rows,]
	test.dove <- dove[-train.rows,]

### Create environmental data frame for modelling

#	dove.env <- data.frame(pres = dove$pres, alt = extract(alt, dove[,1:2]), tmin.jan = extract(tmin.jan, dove[,1:2]), tmin.aug = extract(tmin.aug, dove[,1:2]), temp.range = extract(temp.range, dove[,1:2]), ppt.jan = extract(ppt.jan, dove[,1:2]), ppt.aug = extract(ppt.aug, dove[,1:2]), ppt.season = extract(ppt.season, dove[,1:2]), fire = extract(fire, dove[,1:2]))

dove.env <- data.frame(pres = train.dove$pres, alt = extract(alt, train.dove[,1:2]), tmin.jan = extract(current$tmin.jan, train.dove[,1:2]), tmin.aug = extract(current$tmin.aug, train.dove[,1:2]), ppt.aug = extract(current$ppt.aug, train.dove[,1:2]), ppt.season = extract(current$ppt.season, train.dove[,1:2]), fire = extract(current$fire, train.dove[,1:2]))

test.env <- data.frame(pres = test.dove$pres, alt = extract(alt, test.dove[,1:2]), tmin.jan = extract(current$tmin.jan, test.dove[,1:2]), tmin.aug = extract(current$tmin.aug, test.dove[,1:2]), ppt.aug = extract(current$ppt.aug, test.dove[,1:2]), ppt.season = extract(current$ppt.season, test.dove[,1:2]), fire = extract(current$fire, test.dove[,1:2]))

# Test for correlated variables

	cor(dove.env) # They are all far too correlated - not to worry for now

# Define global model with all candidate predictors
	
	dove.mod <- glm(pres ~ tmin.jan + ppt.season + fire + ppt.aug * tmin.aug, family = binomial(link='logit'), data=dove.env)

# Fit all candidate models

	models.dove <- dredge(dove.mod)

# Average models with delta < 4

	dove.avg <- model.avg(subset(models.dove, delta < 4), fit=TRUE)

	dove.eval <- evaluate(test.env[test.env$pres == 1,], test.env[test.env$pres == 0,], dove.avg)


# Predictions

	env.vars <- data.frame(tmin.jan = as.data.frame(as(current$tmin.jan, "SpatialGridDataFrame"))[,1], ppt.season = as.data.frame(as(current$ppt.season, "SpatialGridDataFrame"))[,1], ppt.aug = as.data.frame(as(current$ppt.aug, "SpatialGridDataFrame"))[,1], tmin.aug = as.data.frame(as(current$tmin.aug, "SpatialGridDataFrame"))[,1], fire = as.data.frame(as(current$fire, "SpatialGridDataFrame"))[,1])

	pred <- predict(dove.avg, newdata = env.vars)

	pred.df <- data.frame(data = pred, x = as.data.frame(as(current$ppt.aug, "SpatialGridDataFrame"))[,2], y = as.data.frame(as(current$ppt.aug, "SpatialGridDataFrame"))[,3])

	coordinates(pred.df) <- ~x+y

	gridded(pred.df) <- TRUE

	pred.sp <- raster(pred.df, layer='data')
	pred.sp <- crop(pred.sp, extent(GRD))

# Mask by land cover

	pred.sp.mask <- pred.sp - cov
	rcl <- matrix(c(-1.5,0,0), nrow=1, ncol=3, byrow=T)
	pred.sp.mask <- reclassify(pred.sp.mask, rcl)

#####################################

## Model simulation functions ##

## Functions to 'grow' patches

f.prob <- function(x, p.trans){
          x*p.trans
 }

f.bin <- function(x){
    if(!is.na(x)) x <- rbinom(1,1,x)
    else x <- NA
 }

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


### Fire risk increase function

fire.risk <- function(fire, p.trans=0.05, iterations=50,  fires = 10, threshold = 0.7, total = 40){
	rcl <- matrix(c(-1,threshold,NA,(threshold + 0.01),1.2,0), nrow=2, ncol=3, byrow=T)
	potential <- reclassify(fire, rcl)

	rnd.pts <- randomPoints(potential, fires)
	rnd.pts <- as.data.frame(rnd.pts)
	potential.pts <- rasterize(rnd.pts, potential, field = 1, background=0)
  potential.pts
	potential.pts <- potential + potential.pts

	wildfire <- f.Exp(potential.pts, potential.pts, p.trans)
	for(i in 1:iterations){
		if(cellStats(wildfire, sum) <= total) wildfire <- f.Exp(wildfire, potential.pts, p.trans)
		else wildfire <- wildfire
	    }

	return(wildfire)
}

### Development increase function

development <- function(cov, p.trans=0.1, iterations=50,  hotels = 10, total = 40){
	rnd.pts <- randomPoints(cov, hotels)
	rnd.pts <- as.data.frame(rnd.pts)
  potential.pts <- rasterize(rnd.pts, cov, field = 1, background=0)
	potential.pts <- cov + potential.pts

	develop <- f.Exp(potential.pts, cov, p.trans)
	for(i in 1:iterations){
		if(cellStats(develop, sum) <= total) develop <- f.Exp(develop, potential.pts, p.trans)
		else develop <- develop
	    }

	return(develop)
}

###### Future scenarios function #####

climate.fun <- function(model, future){

	env.vars <- data.frame(tmin.jan = as.data.frame(as(future$tmin.jan, "SpatialGridDataFrame"))[,1], ppt.season = as.data.frame(as(future$ppt.season, "SpatialGridDataFrame"))[,1], ppt.aug = as.data.frame(as(future$ppt.aug, "SpatialGridDataFrame"))[,1], tmin.aug = as.data.frame(as(future$tmin.aug, "SpatialGridDataFrame"))[,1], fire = as.data.frame(as(future$fire, "SpatialGridDataFrame"))[,1])

	pred <- predict(model, newdata = env.vars)

	pred.df <- data.frame(data = pred, x = as.data.frame(as(future$ppt.aug, "SpatialGridDataFrame"))[,2], y = as.data.frame(as(future$ppt.aug, "SpatialGridDataFrame"))[,3])

	coordinates(pred.df) <- ~x+y

	gridded(pred.df) <- TRUE

	pred.sp <- raster(pred.df, layer='data')
	
	return(pred.sp)
}


### Need to supply name of model from above - default predicts under all scenarios

future <- function(model, sea = TRUE, develop = TRUE, fire = TRUE, climate = TRUE, all = TRUE, reps=10 , threshold = 0.7){

# Initial area predicted present

	pred.current <- climate.fun(model, stack("Climate/current.grd"))

	thresh.rcl <- matrix(c(-1, threshold, 0, threshold, 1.5, 1), nrow=2, ncol=3, byrow=TRUE)

	current.area <- cellStats(reclassify(pred.current, thresh.rcl), mean)

# Create list to contain results

	pred.area <- list(c(reps = reps, threshold = threshold, pred.area = current.area))

# Conditional statements to add elements to the list

	if(climate == TRUE) pred.area$climate <- data.frame(low_26 = numeric(1), mid_45 = numeric(1), high_60 = numeric(1), max_80 = numeric(1))

	if(sea == TRUE) pred.area$sea <- data.frame(sea = numeric(1), sea.low_26 = numeric(1), sea.mid_45 = numeric(1), sea.high_60 = numeric(1), sea.max_80 = numeric(1))

	if(develop == TRUE) pred.area$develop <- data.frame(develop = numeric(reps), dev.low_26 = numeric(reps), dev.mid_45 = numeric(reps), dev.high_60 = numeric(reps), dev.max_80 = numeric(reps))

	if(fire == TRUE) pred.area$fire <- data.frame(fire = numeric(reps), fire.low_26 = numeric(reps), fire.mid_45 = numeric(reps), fire.high_60 = numeric(reps), fire.max_80 = numeric(reps))

	if(all == TRUE) pred.area$all <- data.frame(sea.dev = numeric(reps), sea.fire = numeric(reps), dev.fire = numeric(reps), sea.dev.fire = numeric(reps), all.low_26 = numeric(reps), all.mid_45 = numeric(reps), all.high_60 = numeric(reps), all.max_80 = numeric(reps))

#Calculate predicted areas under scenarios

#Climate - single values
	if(climate == TRUE){

	# Predict suitability across area
		low_26 <- climate.fun(model, stack("Climate/low_26.grd"))
		mid_45 <- climate.fun(model, stack("Climate/mid_45.grd"))
		high_60 <- climate.fun(model, stack("Climate/high_60.grd"))
		max_80 <- climate.fun(model, stack("Climate/max_80.grd"))

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
		sea.map <- climate.fun(model, stack("Climate/low_26.grd"))
			sea.mask <- sea.map + sea.lvl
			pred.area$sea$sea.low_26 <- cellStats(reclassify(sea.mask, thresh.rcl), mean)
		sea.map <- climate.fun(model, stack("Climate/mid_45.grd"))
			sea.mask <- sea.map + sea.lvl
			pred.area$sea$sea.mid_45 <- cellStats(reclassify(sea.mask, thresh.rcl), mean)
		sea.map <- climate.fun(model, stack("Climate/high_60.grd"))
			sea.mask <- sea.map + sea.lvl
			pred.area$sea$sea.high_60 <- cellStats(reclassify(sea.mask, thresh.rcl), mean)
		sea.map <- climate.fun(model, stack("Climate/max_80.grd"))
			sea.mask <- sea.map + sea.lvl
			pred.area$sea$sea.max_80 <- cellStats(reclassify(sea.mask, thresh.rcl), mean)
	
		rm(sea.map, sea.mask)
	}

# Development
	if(develop == TRUE){
		low_26 <- climate.fun(model, stack("Climate/low_26.grd"))
		mid_45 <- climate.fun(model, stack("Climate/mid_45.grd"))
		high_60 <- climate.fun(model, stack("Climate/high_60.grd"))
		max_80 <- climate.fun(model, stack("Climate/max_80.grd"))

		for(i in 1:reps){
		# Simulate development with default values
			dev.mask <- development(cov)
		
		# Use simulated development to mask climate predictions
			pred.area$develop$develop[i] <- cellStats(reclassify(reclassify(pred.current, thresh.rcl) - dev.mask, matrix(c(-1.5,0,0), nrow=1, ncol=3, byrow=T) ), mean)
			pred.area$develop$dev.low_26[i] <- cellStats(reclassify(reclassify(low_26, thresh.rcl) - dev.mask, matrix(c(-1.5,0,0), nrow=1, ncol=3, byrow=T)), mean)
			pred.area$develop$dev.mid_45[i] <- cellStats(reclassify(reclassify(mid_45, thresh.rcl) - dev.mask, matrix(c(-1.5,0,0), nrow=1, ncol=3, byrow=T)), mean)
			pred.area$develop$dev.high_60[i] <- cellStats(reclassify(reclassify(high_60, thresh.rcl) - dev.mask, matrix(c(-1.5,0,0), nrow=1, ncol=3, byrow=T)), mean)
			pred.area$develop$dev.max_80[i] <- cellStats(reclassify(reclassify(max_80, thresh.rcl) - dev.mask, matrix(c(-1.5,0,0), nrow=1, ncol=3, byrow=T)), mean)
		}
		rm(low_26, mid_45, high_60, max_80)
	}

# Wildfire
	if(fire == TRUE){
		low_26 <- climate.fun(model, stack("Climate/low_26.grd"))
		mid_45 <- climate.fun(model, stack("Climate/mid_45.grd"))
		high_60 <- climate.fun(model, stack("Climate/high_60.grd"))
		max_80 <- climate.fun(model, stack("Climate/max_80.grd"))

		for(i in 1:reps){
		# Simulate fire with default values
			fire.dis.current <- fire.risk(stack("Climate/current.grd")$fire)			
			fire.dis.low_26 <- fire.risk(stack("Climate/low_26.grd")$fire)
			fire.dis.mid_45 <- fire.risk(stack("Climate/mid_45.grd")$fire)
			fire.dis.high_60 <- fire.risk(stack("Climate/high_60.grd")$fire)
			fire.dis.max_80 <- fire.risk(stack("Climate/max_80.grd")$fire)
			fire.dis.current <- reclassify(fire.dis.current, matrix(c(-1,0.99999,0, 0.999999, 1, 1), nrow=2, ncol=3, byrow=T))
			fire.dis.low_26 <- reclassify(fire.dis.low_26, matrix(c(-1,0.99999,0, 0.999999, 1, 1), nrow=2, ncol=3, byrow=T))	
			fire.dis.mid_45 <- reclassify(fire.dis.mid_45, matrix(c(-1,0.99999,0, 0.999999, 1, 1), nrow=2, ncol=3, byrow=T))
			fire.dis.high_60 <- reclassify(fire.dis.high_60, matrix(c(-1,0.99999,0, 0.999999, 1, 1), nrow=2, ncol=3, byrow=T))
			fire.dis.max_80 <- reclassify(fire.dis.max_80, matrix(c(-1,0.99999,0, 0.999999, 1, 1), nrow=2, ncol=3, byrow=T))
		
		# Use simulated development to mask climate predictions
			pred.area$fire$fire[i] <- cellStats(reclassify(reclassify(pred.current, thresh.rcl) - fire.dis.current, matrix(c(-1.5,0,0), nrow=1, ncol=3, byrow=T)), mean)
			pred.area$fire$fire.low_26[i] <- cellStats(reclassify(reclassify(low_26, thresh.rcl) - fire.dis.low_26, matrix(c(-1.5,0,0), nrow=1, ncol=3, byrow=T)), mean)
			pred.area$fire$fire.mid_45[i] <- cellStats(reclassify(reclassify(mid_45, thresh.rcl) - fire.dis.mid_45, matrix(c(-1.5,0,0), nrow=1, ncol=3, byrow=T)), mean)
			pred.area$fire$fire.high_60[i] <- cellStats(reclassify(reclassify(high_60, thresh.rcl) - fire.dis.high_60, matrix(c(-1.5,0,0), nrow=1, ncol=3, byrow=T)), mean)
			pred.area$fire$fire.max_80[i] <- cellStats(reclassify(reclassify(max_80, thresh.rcl) - fire.dis.max_80, matrix(c(-1.5,0,0), nrow=1, ncol=3, byrow=T)), mean)
		}
		rm(low_26, mid_45, high_60, max_80)
	}

# All
	if(all == TRUE){
		low_26 <- climate.fun(model, stack("Climate/low_26.grd"))
		mid_45 <- climate.fun(model, stack("Climate/mid_45.grd"))
		high_60 <- climate.fun(model, stack("Climate/high_60.grd"))
		max_80 <- climate.fun(model, stack("Climate/max_80.grd"))
		
		for(i in 1:reps){	
			# Reclassify prediction
			current.bin <- reclassify(pred.current, thresh.rcl)			
			low_26 <- reclassify(low_26, thresh.rcl)
			mid_45 <- reclassify(mid_45, thresh.rcl)			
			high_60 <- reclassify(high_60, thresh.rcl)
			max_80 <- reclassify(max_80, thresh.rcl)

			# Simulate development with default values
			dev.mask <- development(cov)
		
			# Simulate fire with default values
			fire.dis.current <- fire.risk(stack("Climate/current.grd")$fire)			
			fire.dis.low_26 <- fire.risk(stack("Climate/low_26.grd")$fire)
			fire.dis.mid_45 <- fire.risk(stack("Climate/mid_45.grd")$fire)
			fire.dis.high_60 <- fire.risk(stack("Climate/high_60.grd")$fire)
			fire.dis.max_80 <- fire.risk(stack("Climate/max_80.grd")$fire)
			fire.dis.current <- reclassify(fire.dis.current, matrix(c(-1,0.99999,0, 0.999999, 1, 1), nrow=2, ncol=3, byrow=T))
			fire.dis.low_26 <- reclassify(fire.dis.low_26, matrix(c(-1,0.99999,0, 0.999999, 1, 1), nrow=2, ncol=3, byrow=T))	
			fire.dis.mid_45 <- reclassify(fire.dis.mid_45, matrix(c(-1,0.99999,0, 0.999999, 1, 1), nrow=2, ncol=3, byrow=T))
			fire.dis.high_60 <- reclassify(fire.dis.high_60, matrix(c(-1,0.99999,0, 0.999999, 1, 1), nrow=2, ncol=3, byrow=T))
			fire.dis.max_80 <- reclassify(fire.dis.max_80, matrix(c(-1,0.99999,0, 0.999999, 1, 1), nrow=2, ncol=3, byrow=T))

			# Sea level and development
			mask <- current.bin + sea.lvl
			mask <- reclassify(mask - dev.mask, matrix(c(-1.5,0,0)))
			pred.area$all$sea.dev[i] <- cellStats(mask, mean)

			# Sea level and fire
			mask <- current.bin + sea.lvl
			mask <- reclassify(mask - fire.dis.current, matrix(c(-1.5,0,0)))
			pred.area$all$sea.fire[i] <- cellStats(mask, mean)

			# Fire and development
			mask <- reclassify(current.bin - dev.mask, matrix(c(-1.5,0,0)))
			mask <- reclassify(mask - fire.dis.current, matrix(c(-1.5,0,0)))
			pred.area$all$dev.fire[i] <- cellStats(mask, mean)

			# Sea level, development and fire
			# Current values			
			mask <- current.bin + sea.lvl
			mask <- reclassify(mask - dev.mask, matrix(c(-1.5,0,0)))
			mask <- reclassify(mask - fire.dis.current, matrix(c(-1.5,0,0)))
			pred.area$all$sea.dev.fire[i] <- cellStats(mask, mean)
				
			# Low_26
			mask <- low_26 + sea.lvl
			mask <- reclassify(mask - dev.mask, matrix(c(-1.5,0,0)))
			mask <- reclassify(mask - fire.dis.low_26, matrix(c(-1.5,0,0)))
			pred.area$all$all.low_26[i] <- cellStats(mask, mean)

			# Mid_45
			mask <- mid_45 + sea.lvl
			mask <- reclassify(mask - dev.mask, matrix(c(-1.5,0,0)))
			mask <- reclassify(mask - fire.dis.mid_45, matrix(c(-1.5,0,0)))
			pred.area$all$all.mid_45[i] <- cellStats(mask, mean)

			# High_60
			mask <- high_60 + sea.lvl
			mask <- reclassify(mask - dev.mask, matrix(c(-1.5,0,0)))
			mask <- reclassify(mask - fire.dis.high_60, matrix(c(-1.5,0,0)))
			pred.area$all$all.high_60[i] <- cellStats(mask, mean)

			# Max_80
			mask <- max_80 + sea.lvl
			mask <- reclassify(mask - dev.mask, matrix(c(-1.5,0,0)))
			mask <- reclassify(mask - fire.dis.max_80, matrix(c(-1.5,0,0)))
			pred.area$all$all.max_80[i] <- cellStats(mask, mean)

		}
		rm(low_26, mid_45, high_60, max_80, mask, dev.mask, fire.dis)
	}

	return(simulation.predictions = pred.area)
}

##########
 
