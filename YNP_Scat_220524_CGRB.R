
#############################
# Yosemite Scat Analysis    #
# Scent-detection dog teams #
# David Seth Green          #
# Started 11 August 2021    #
#############################

library(reshape2)
library(DBI)
library(rgdal)
library(raster)
library(rgeos)
library(sp)
library(dismo)
library(dplyr)
library(jagsUI)
library(plotly)
library(tinytex)
library(rjags)
library(dclone)
library(snow)
library(spdep)
library(fields)
library(concaveman)
library(exactextractr)
library(sf)
library(rdist)
library(lubridate)
library(readxl)

# Read in camera datasets
setwd("C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/YosemiteCougar/Data")
photos <- read_xlsx("CPWFlatFileVERIFIEDONLY20210920.xlsx")
visits <- read_xlsx("VisitsFlatFileYosemite20210920.xlsx")

# Add Day/Month/Year to photos
photos <- photos %>% dplyr::mutate(year = lubridate::year(ImageDate),
                                   month = lubridate::month(ImageDate),
                                   day = lubridate::day(ImageDate))
# Add hour to photos
photos <- photos %>% dplyr::mutate(hour = lubridate::hour(ImageDate))


# Add Julian date to look at seasonality
photos$JulianDate <- julian(photos$ImageDate) 

photos$uniquedet <- as.numeric(paste0(photos$LocationID,photos$year,photos$month,photos$day))

# Create a matrix of cougar detections at cameras per year
photos2 <- photos[photos$CommonName == "Mountain Lion" & photos$year > 2018,] # Subset to only cougars and year >2018
photos2$Detection <- ifelse(is.na(photos2$SpeciesID), NA, ifelse(photos2$SpeciesID > 0, 1, 0)) # Add detection to photo database for those other than nothing, and NA for those that have not been tagged
photos2 <- photos2[!is.na(photos2$Detection),] # Remove photos without tagging
unique(photos2$uniquedet)
write.csv(photos2, "cougarphotos_220727.csv", row.names = FALSE)

dets <- melt(photos2, id.vars = c("LocationName","day","month","year"), measure.vars = "Detection", na.rm = F, value.name = "Presence")
dets2 <- acast(dets, LocationName ~ year ~ month ~ day, value.var = "Presence", na.rm = F, fun.aggregate = max)
class(dets2) <- "numeric"
dets2[is.infinite(dets2)] <- 0 # replace "Inf" with 0's
YNPOcc <- dets2

# Determine the total number of days that cougars were detected each year
YNPOcc.Year <- array(NA, dim = c(dim(YNPOcc)[1], dim(YNPOcc)[2]))
for(j in 1:dim(YNPOcc.Year)[1]){
  for(t in 1:dim(YNPOcc.Year)[2]){
    YNPOcc.Year[j,t] <- sum(YNPOcc[j,t,,])
  }
}

# Create deployed matrix of cameras (determining how long they were deployed)
dep <- melt(photos2, id.vars = c("LocationName","day","month","year"), measure.vars = "Detection", na.rm = F, value.name = "Deployed")
dep2 <- acast(dep, LocationName ~ year ~ month ~ day, value.var = "Deployed", na.rm = F, fun.aggregate = max)
class(dep2) <- "numeric"
deployed <- dep2
deployed <- ifelse(is.infinite(deployed), 0, 1)

# Create trap matrix of cameras
trapmat <- photos2[!duplicated(photos2$LocationName),c("LocationName","UTM_E","UTM_N")]
trapmat <- arrange(trapmat, by = LocationName)
trapmat$LocationNum <- c(1:nrow(trapmat))

# Test by plotting
plot(trapmat$UTM_E, trapmat$UTM_N)

# Bring in YNP Boundary
#setwd("E:/Pekania/Projects/Yosemite/YNP_MountainLions/Data/GIS")
YNP_boundary <- readOGR(dsn = ".", layer = "yose_boundaryNAD83Z11")

# Bring in scat and track datasets
#setwd("E:/Pekania/Projects/Yosemite/YNP_MountainLions/Data/DetectionDogsFinal")
scat <- read.csv("ScatsProjected_NAD83Z11.csv", header = T)

# Bring in merged track set to allow for division by year
alltracks <- readOGR(dsn = ".", layer = "Tracks20210712")

# # Double check that they look right
# plot(YNP_boundary)
# plot(alltracks, add = T)
# points(scat$POINT_X, scat$POINT_Y, col = "red", cex = 0.5)
# # They look good and aligned

# Drop the scats that were collected opportunistically but in areas without effort or those without a location
scat2 <- scat[!(scat$SampleIDFi == "20191010_1" | scat$SampleIDFi == "20191012_1" | scat$POINT_X == 0),]

# # Double check again
# plot(scat2$POINT_X, scat2$POINT_Y, col = "red", cex = 0.5)
# plot(alltracks, add = T)
# # Much better


# Add year to track data
alltracks@data$Year <- year(alltracks$TrackDate)

# Reorganize all tracks per season
s19 <- alltracks[alltracks$Year == 2019,]
s20 <- alltracks[alltracks$Year == 2020,]

# # Test plot
# plot(s19, col = "red")
# plot(s20, col = "blue", add = T)
# # Looks good

# Buffer YNP boundary to create a study area
trackpoints <- as(alltracks, "SpatialPointsDataFrame")
t <- concaveman(trackpoints@coords, concavity = 2)
p <- Polygon(t)
ps <- Polygons(list(p),1)
sps <- SpatialPolygons(list(ps))
SA <- gBuffer(sps, width = 22000)  # Buffer by 2-3*sigma (11000 for males)
crs(SA) <- crs(YNP_boundary) # Assign a CRS to the boundary

area <- rgeos::gArea(SA,byid=TRUE)
area <- area/1000000
# Write to shapefile
shapefile(x = SA, file = "YNP_CougarStudyarea.shp")

# Create the grids to use within the defined state-space
cellsize = 9000 # Size of grid cell width/length; in meters; close to size of average home ranges of cougars in YNP (6-11km sigma)
cellsize2 = 3000 # Size of smaller, detection grid
fishnet <- spsample(SA, type = "regular", pretty = T, n = 15000, cellsize = c(cellsize,cellsize), offset = 0)
fishnet2 <- spsample(SA, type = "regular", pretty = T, n = 15000, cellsize = c(cellsize2,cellsize2), offset = 0)

# # Double check the look of these
# plot(SA)
# plot(fishnet, add = T, pch = 1, cex = 0.5)
# plot(fishnet2, add = T, pch = 2, cex = 0.5)
# plot(alltracks, add = T, col = "blue")
# plot(YNP_boundary, add = T)
# # Look good!

# Create gridID for each grid centroid (both large and small)
grid <- data.frame(cbind(fishnet@coords[,1], fishnet@coords[,2]))
grid <- arrange(grid, X1, X2)
grid$GridID <- as.numeric(rownames(grid))
grid2 <- data.frame(cbind(fishnet2@coords[,1], fishnet2@coords[,2]))
grid2 <- arrange(grid2, X1, X2)
grid2$GridID <- as.numeric(rownames(grid2))

# Create square buffers around grid centers (to create the polygon for adding any kind of landscape features within these grids)
radius <- cellsize/2 # Radius is the width/2
xMinus <- grid$X1-radius
xPlus <- grid$X1+radius
yMinus <- grid$X2-radius
yPlus <- grid$X2+radius
radius2 <- cellsize2/2 # Radius is the width/2
xMinus2 <- grid2$X1-radius
xPlus2 <- grid2$X1+radius
yMinus2 <- grid2$X2-radius
yPlus2 <- grid2$X2+radius

# Calculate polygon coordinates for each plot centroid to allow for assigning information to the grid cells
square=cbind(xMinus,yPlus,  # NW corner
             xPlus, yPlus,  # NE corner
             xPlus, yMinus,  # SE corner
             xMinus, yMinus, # SW corner
             xMinus, yPlus)  # NW corner again - close ploygon
ID <- grid$GridID

square2=cbind(xMinus2,yPlus2,  # NW corner
              xPlus2, yPlus2,  # NE corner
              xPlus2, yMinus2,  # SE corner
              xMinus2, yMinus2, # SW corner
              xMinus2, yPlus2)  # NW corner again - close ploygon
ID2 <- grid2$GridID

# Now make these into polygons
grid_poly <- SpatialPolygons(mapply(function(poly, id)
{
  xy <- matrix(poly, ncol=2, byrow=TRUE)
  Polygons(list(Polygon(xy)), ID=id)
},
split(square, row(square)), ID),
proj4string=CRS(as.character(proj4string(alltracks))))

grid_poly2 <- SpatialPolygons(mapply(function(poly, id)
{
  xy <- matrix(poly, ncol=2, byrow=TRUE)
  Polygons(list(Polygon(xy)), ID=id)
},
split(square2, row(square2)), ID2),
proj4string=CRS(as.character(proj4string(alltracks))))

# # Double check the look of these
# plot(grid_poly2) # Smaller detection grid cell
# plot(s19, col = "red", add = T)
# plot(s20, col = "blue", add = T)
# # Look good!

# Determine which small cell each scat was located in
scat_sp <- SpatialPoints(cbind(scat2$POINT_X, scat2$POINT_Y), proj4string = CRS("+proj=utm +zone=11 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))
scatpix <- over(scat_sp, grid_poly2)
scat2$Pix <- scatpix # Add this to the dataframe

# Determine which small cell each camera was located in
camera_sp <- SpatialPoints(cbind(trapmat$UTM_E, trapmat$UTM_N), proj4string = CRS("+proj=utm +zone=11 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))
campix <- over(camera_sp, grid_poly2)
trapmat$Grid <- campix # Add this to the trapmatrix

# Determine which of the larger cells the smaller ones fall in
grid2_sp <- SpatialPoints(cbind(grid2$X1, grid2$X2), proj4string = CRS("+proj=utm +zone=11 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))
gridingridpix <- over(grid2_sp, grid_poly)

# Add year to the scat
scat2$Date <- as.Date(scat2$Date, format = "%m/%d/%Y")
scat2$Year <- year(scat2$Date)

# Subset scat locations by species
scat_cougar <- scat2[scat2$SpeciesIDL == "mountain lion",]
scat_marten <- scat2[scat2$SpeciesIDL == "marten",]
scat_bobcat <- scat2[scat2$SpeciesIDL == "bobcat",]
scat_coyote <- scat2[scat2$SpeciesIDL == "coyote",]
scat_grayfox <- scat2[scat2$SpeciesIDL == "gray fox",]
scatsp_cougar <- SpatialPoints(cbind(scat_cougar$POINT_X, scat_cougar$POINT_Y), proj4string = CRS("+proj=utm +zone=11 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))
scatsp_marten <- SpatialPoints(cbind(scat_marten$POINT_X, scat_marten$POINT_Y), proj4string = CRS("+proj=utm +zone=11 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))
scatsp_bobcat <- SpatialPoints(cbind(scat_bobcat$POINT_X, scat_bobcat$POINT_Y), proj4string = CRS("+proj=utm +zone=11 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))
scatsp_coyote <- SpatialPoints(cbind(scat_coyote$POINT_X, scat_coyote$POINT_Y), proj4string = CRS("+proj=utm +zone=11 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))
scatsp_grayfox <- SpatialPoints(cbind(scat_grayfox$POINT_X, scat_grayfox$POINT_Y), proj4string = CRS("+proj=utm +zone=11 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))

#### summary
scat_cougar2 <- scat_cougar[!(scat_cougar$CougarID == "fail"),]
scat_cougar$count <- rep(1, length(scat_cougar$SampleIDFi))

options(digits=2)
scatsum <- scat_cougar2 %>%
  group_by(CougarID) %>% 
  summarize(ave = sum(count),
            sd = sd(sum(count)))

meancoug <- mean(scatsum$ave) #2.80
sdcoug <- sd(scatsum$ave) #2.65

# # Test plot
# plot(grid_poly2)
# points(scatsp_cougar, col = "purple", cex = 0.5)
# points(scatsp_marten, col = "blue", cex = 0.5)
# points(scatsp_bobcat, col = "orange", cex = 0.5)
# points(scatsp_coyote, col = "brown", cex = 0.5)
# points(scatsp_grayfox, col = "red", cex = 0.5)
# Looks good

# Determine for each year if a non-cougar scat was found in each small cell
carscat_year <- array(NA, dim = c(length(grid_poly2),4,2)) # 4 other species, 2 years

# Determine the cells with detections
for(i in 1:dim(carscat_year)[1]){
  carscat_year[i,1,1] <- ifelse(any(scat_marten$Pix==i & scat_marten$Year == 2019),1,0)
  carscat_year[i,1,2] <- ifelse(any(scat_marten$Pix==i & scat_marten$Year == 2020),1,0)
  carscat_year[i,2,1] <- ifelse(any(scat_bobcat$Pix==i & scat_bobcat$Year == 2019),1,0)
  carscat_year[i,2,2] <- ifelse(any(scat_bobcat$Pix==i & scat_bobcat$Year == 2020),1,0)
  carscat_year[i,3,1] <- ifelse(any(scat_coyote$Pix==i & scat_coyote$Year == 2019),1,0)
  carscat_year[i,3,2] <- ifelse(any(scat_coyote$Pix==i & scat_coyote$Year == 2020),1,0)
  carscat_year[i,4,1] <- ifelse(any(scat_grayfox$Pix==i & scat_grayfox$Year == 2019),1,0)
  carscat_year[i,4,2] <- ifelse(any(scat_grayfox$Pix==i & scat_grayfox$Year == 2020),1,0)
}

# Determine length of lines in each small grid cell per year to account for effort
grid_sf <- st_as_sf(grid_poly2) # Turn grid into sf
y19_sf <- st_set_precision(st_as_sf(s19), 1e5) # Set precision to help with really long UTMs
y20_sf <- st_set_precision(st_as_sf(s20), 1e5) # Set precision to help with really long UTMs

ints19 = st_intersection(y19_sf,grid_sf)
ints19$len = st_length(ints19)
grid_sf$Ids19 = 1:nrow(grid_sf)
joins19 = st_join(grid_sf, ints19)
outs19 = group_by(joins19, Ids19) %>%
  summarize(length = sum(len))

ints20 = st_intersection(y20_sf,grid_sf)
ints20$len = st_length(ints20)
grid_sf$Ids20 = 1:nrow(grid_sf)
joins20 = st_join(grid_sf, ints20)
outs20 = group_by(joins20, Ids20) %>%
  summarize(length = sum(len))

# Merge effort together
library(abind)
effort_year <- abind(outs19$length, outs20$length, along = 2)
effort_year[is.na(effort_year)] <- 0 # Replace NA's with 0

# standardize track effort per small grid cell
effort_yearstd <- array(NA, dim = c(dim(effort_year)[1],2)) 
for(i in 1:dim(effort_yearstd)[1]){
  for(t in 1:dim(effort_yearstd)[2]){
    effort_yearstd[i,t] <- (effort_year[i,t]-mean(effort_year[,]))/sd(effort_year[,])
  }
}

# determine if sampled at all in a season; needs to have 0 probability of detection for sites with no effort
year_sampled <- ifelse(effort_yearstd > min(effort_yearstd), 1, 0)

# For each identified cougar, determine if they were captured in each small cell each year
# Step 1. Select only cougar locations with ID's
scat_cougar2 <- scat_cougar[!(scat_cougar$CougarID == "fail"),]

# Step 2. Merge the scat dataset with the grid ID's so that all are represented
scat_cougar3 <- merge(scat_cougar2, grid2, by.x = "Pix", by.y = "GridID", all.y = T)

# Step 3. Melt and cast this to be detection of each individaul x grid x year
scat_cougar3$present <- 1 # Add 1 to say they were detected
toot <- melt(scat_cougar3, id.vars = c("CougarID", "Year", "Pix"),  measure.var = "present", na.rm = F)
toot2 <- acast(toot, Pix ~ CougarID ~ Year, sum)
y_cougar <- toot2[,1:(dim(toot2)[2])-1,1:(dim(toot2)[3])-1] # Remove the place holder for NA's

# Extract sex for cougars
cougID <- data.frame(dimnames(y_cougar)[[2]])
k <- scat_cougar[,c("CougarID","CougarSex")]
k2 <- unique(k)
k3 <- merge(cougID, k2, by.x = "dimnames.y_cougar...2..", by.y = "CougarID")
Sex <- ifelse(k3$CougarSex == "M", 1, 0)

# Step 4. Add ghost critters (augmentation)
n <- dim(y_cougar)[2] # Number of unique cougars detected
M <- 175              # Total number with augmentation
nz <- M-n             # Ghost critters to add
Y_cougar <- abind(y_cougar, array(0, dim = c(dim(y_cougar)[1], nz, 2)), along = 2)
Sex2 <- c(Sex, rep(NA, nz))

# Habitat covariate stuffs; need to reorgnize the structure of these
p <- st_as_sf(grid_poly) # Pull out grid cells
p2 <- st_as_sf(grid_poly2) # smaller, detection grid

# Covariates previously coded, but determined not worth pursuing: 
# elevation, aspect, canopy, surface fuels, canopy layers, development
# Had tried looking at binary variables for rivers, trails, roads, but distance was actually better

# Slope/aspect (From USGS SRTM)
#setwd("E:/Pekania/Projects/Yosemite/YNP_MountainLions/Data/GIS")
slope <- raster("USGS_SRTM_Slope_NAD83Z11.tif")

# Canopy cover and canopy layers (from California Forest Observatory; 2020)
#setwd("E:/Pekania/Projects/Yosemite/YNP_MountainLions/Data/GIS/")
ndvi <- raster("meanNDVI_NAD83Z11.tif")

# Bring in the shapefiles for development, trails, roads, rivers, and lakes/marshes
trails <- readOGR("Trails_NAD83Z11.shp")
crs(trails) <- crs(grid_poly) # Need to rename to get rid of vertical
roads <- readOGR("YosemiteRoads_NAD83Z11.shp")
usgsrivers <- readOGR("USGS_NHD_Perennial_NAD83Z11_50kbuff.shp") # Perennial streams from USGS
yoserivers <- readOGR("Rivers_NAD83Z11.shp") # Rivers from YOSE
rivers <- bind(usgsrivers, yoserivers) # Merge the 2 river shapefiles

# # Determine distance to nearest trail, road, river
grid2.pts <- SpatialPoints(cbind(grid2$X1, grid2$X2))
trail.dist <- as.matrix(apply(gDistance(grid2.pts, trails, byid=TRUE),2,min))
road.dist <- as.matrix(apply(gDistance(grid2.pts, roads, byid=TRUE),2,min))
river.dist <- as.matrix(apply(gDistance(grid2.pts, rivers, byid=TRUE),2,min))

# Extract averages within grid cells
slope.mean <- exact_extract(slope, p2, 'mean')
ndvi.mean <- exact_extract(ndvi, p, 'mean')

# Standardize covs among grid cells
slope.std <- array(NA, dim = c(nrow(grid2))) # Grid rows and 1 year(doesn't change)
for(i in 1:nrow(grid2)){
  slope.std[i] <- (slope.mean[i]-mean(slope.mean[], na.rm = T))/(sd(slope.mean[], na.rm = T))
}
ndvi.std <- array(NA, dim = c(nrow(grid))) # Grid rows and 1 year(doesn't change)
for(i in 1:nrow(grid)){
  ndvi.std[i] <- (ndvi.mean[i]-mean(ndvi.mean[], na.rm = T))/(sd(ndvi.mean[], na.rm = T))
}
trail.std <- array(NA, dim = c(nrow(grid2))) # Grid rows and 1 year(doesn't change)
for(i in 1:nrow(grid2)){
  trail.std[i] <- (trail.dist[i]-mean(trail.dist[], na.rm = T))/(sd(trail.dist[], na.rm = T))
}
road.std <- array(NA, dim = c(nrow(grid2))) # Grid rows and 1 year(doesn't change)
for(i in 1:nrow(grid2)){
  road.std[i] <- (road.dist[i]-mean(road.dist[], na.rm = T))/(sd(road.dist[], na.rm = T))
}
river.std <- array(NA, dim = c(nrow(grid2))) # Grid rows and 1 year(doesn't change)
for(i in 1:nrow(grid2)){
  river.std[i] <- (river.dist[i]-mean(river.dist[], na.rm = T))/(sd(river.dist[], na.rm = T))
}

# Check for collinearity among habitat covariates
covs <- cbind(river_bin, trail_bin, road_bin, slope.std, river.std, trail.std, road.std, effort_yearstd)
covscor <- cor(covs, use = "complete.obs", method = "pearson")
covscor
# Looks good!

####### CBRB #########
# Define the directory where things will go
library(runjags)
library(jagsUI)

setwd("E:/Pekania/Projects/Yosemite/YNP_MountainLions/Analysis/CGRB/Scat_24May22_175M_22000K_9_3")

# Data
jags_data <- list(y = Y_cougar,
                  effort = effort_yearstd,
                  sex = Sex2,
                  y.camY = YNPOcc.Year,
                  deployed = deployed,
                  camGrid = c(as.matrix(trapmat$Grid)),
                  pixArea = (cellsize*cellsize)/1000000,
                  GG = dim(grid)[1],
                  XX = dim(grid2)[1],
                  M = dim(Y_cougar)[2],
                  TT = dim(Y_cougar)[3],
                  JJ = dim(YNPOcc)[1],
                  MM = dim(YNPOcc)[3],
                  DD = dim(YNPOcc)[4],
                  S = as.matrix(grid[,1:2]),
                  SS = as.matrix(grid2[,1:2]),
                  slope = slope.std,
                  ndvi = ndvi.std,
                  rivers = river.std,
                  trails = trail.std,
                  roads = road.std,
                  gridpix = gridingridpix)

# Initial values
jags_inits1 <- list(".RNG.seed" = 1, ".RNG.name" = "base::Wichmann-Hill", z = array(1, dim = c(dim(Y_cougar)[2])), alpha0 = -5.2, sigma = c(10000,10000), beta0 = -3)#, ro = 0.05)
jags_inits2 <- list(".RNG.seed" = 2, ".RNG.name" = "base::Wichmann-Hill", z = array(1, dim = c(dim(Y_cougar)[2])), alpha0 = -5.2, sigma = c(10000,10000), beta0 = -3)#, ro = 0.05)
jags_inits3 <- list(".RNG.seed" = 3, ".RNG.name" = "base::Wichmann-Hill", z = array(1, dim = c(dim(Y_cougar)[2])), alpha0 = -5.2, sigma = c(10000,10000), beta0 = -3)#, ro = 0.05)

# Write the data to the CGRB
data.dump <- dump.format(jags_data)
write(data.dump, file = "data.txt")
inits.dump <- dump.format(jags_inits1)
write(inits.dump, file = "inits1.txt")
inits.dump <- dump.format(jags_inits2)
write(inits.dump, file = "inits2.txt")
inits.dump <- dump.format(jags_inits3)
write(inits.dump, file = "inits3.txt")

# Write model code
sink("Scat_SCR_DSS_CamY.txt")
cat("
    model {

    # Derived parameters
    N <- sum(z[1:M]) 				      # population size (N)
  
    # SCR encounter process
    for (x in 1:XX){              # For each grid cell
    for (i in 1:M){ 					    # For each individual
    Dist2[x,i] <- sqrt((x0g[i] - SS[x,1])^2 + (y0g[i] - SS[x,2])^2)
    for (t in 1:TT){					    # For each year
    y[x,i,t] ~ dnegbin(p.1[x,i,t],ro)
    p.1[x,i,t] <- ro/(ro + lambda[x,i,t])
    lambda[x,i,t] <- p0[x,i,t]*exp(-gitj[sex2[i]]*Dist2[x,i]*Dist2[x,i])*z[i]
    p0[x,i,t] <- exp(beta0 + beta1*effort[x,t] + beta2*sex[i] + beta3*(t-1) # Males are 1
                + beta4*slope[x] + beta5*trails[x] + beta6*roads[x] + beta7*rivers[x])
    } #t
    } #i
    } #x
    ro ~ dgamma(0.1,0.1)

    # Intensity function for DSS habitat covariates
    EN <- sum(mu[])
    psi <- EN/M
    for(g in 1:GG){
    probs[g] <- mu[g]/EN  # Probability of an activity center being in a grid cell in a year
    mu[g] <- exp(alpha0 + alpha1*ndvi[g])*pixArea
    } #g

    # Make activity centers
    for(i in 1:M){
    sex[i] ~ dbern(psi.sex)
    sex2[i] <- sex[i]+1
    x0g[i] <- S[s[i],1]
    y0g[i] <- S[s[i],2]
    s[i] ~ dcat(probs[])
    z[i] ~ dbern(psi)
    } #i
    
    # Add in the camera count data of days with cougars
    for(j in 1:JJ){
    for(t in 1:TT){
    bigLambda[j,t] <- sum(lambda[camGrid[j],,t])
    y.camY[j,t] ~ dnegbin(p[j,t],r)
    p[j,t] <- r/(r + bigLambda[j,t])
    } #t
    } #j
    r ~ dgamma(0.1,0.1)

    # Movement parameter varies by sex
    for(s in 1:2){
    sigma[s] ~ dunif(0,20000)
    gitj[s] <- 1/(2*sigma[s]*sigma[s])
    } #s

    # Priors
    beta0 ~ dnorm(0, 0.000001)
    alpha0 ~ dnorm(0,0.000001)
    alpha1 ~ dnorm(0,0.000001)
    beta1 ~ dnorm(0, 0.000001)
    beta2 ~ dnorm(0, 0.000001)
    beta3 ~ dnorm(0, 0.000001)
    beta4 ~ dnorm(0, 0.000001)
    beta5 ~ dnorm(0, 0.000001)
    beta6 ~ dnorm(0, 0.000001)
    beta7 ~ dnorm(0, 0.000001)
    psi.sex ~ dunif(0,1)
    } ",fill = TRUE)
sink()

# Write the JAGS Command line
sink("1_Scat.cmd")
cat("
    /* OR Carn
    */
    model in Scat_SCR_DSS_CamY.txt
    data in data.txt
    compile, nchains(1)
    inits in inits1.txt
    initialize
    update 1000
    monitor N, thin(10)
    monitor alpha0, thin(10)
    monitor alpha1, thin(10)
    monitor alpha2, thin(10)
    monitor alpha3, thin(10)
    monitor alpha4, thin(10)
    monitor alpha5, thin(10)
    monitor alpha6, thin(10)
    monitor alpha7, thin(10)
    monitor alpha8, thin(10)
    monitor beta0, thin(10)
    monitor beta1, thin(10)
    monitor beta2, thin(10)
    monitor beta3, thin(10)
    monitor beta4, thin(10)
    monitor beta5, thin(10)
    monitor beta6, thin(10)
    monitor beta7, thin(10)
    monitor sigma, thin(10)
    monitor psi.sex, thin(10)
    monitor r, thin(10)
    monitor ro, thin(10)
    monitor s, thin(10)
    monitor z, thin(10)
    monitor sex, thin(10)
    update 3000
    coda *, stem(Chain1)
    ", fill = T)
sink()

sink("2_Scat.cmd")
cat("
    /* OR Carn
    */
    model in Scat_SCR_DSS_CamY.txt
    data in data.txt
    compile, nchains(1)
    inits in inits2.txt
    initialize
    update 1000
    monitor N, thin(10)
    monitor alpha0, thin(10)
    monitor alpha1, thin(10)
    monitor alpha2, thin(10)
    monitor alpha3, thin(10)
    monitor alpha4, thin(10)
    monitor alpha5, thin(10)
    monitor alpha6, thin(10)
    monitor alpha7, thin(10)
    monitor alpha8, thin(10)
    monitor beta0, thin(10)
    monitor beta1, thin(10)
    monitor beta2, thin(10)
    monitor beta3, thin(10)
    monitor beta4, thin(10)
    monitor beta5, thin(10)
    monitor beta6, thin(10)
    monitor beta7, thin(10)
    monitor sigma, thin(10)
    monitor psi.sex, thin(10)
    monitor r, thin(10)
    monitor ro, thin(10)
    monitor s, thin(10)
    monitor z, thin(10)
    monitor sex, thin(10)
    update 3000
    coda *, stem(Chain2)
    ", fill = T)
sink()

sink("3_Scat.cmd")
cat("
    /* OR Carn
    */
    model in Scat_SCR_DSS_CamY.txt
    data in data.txt
    compile, nchains(1)
    inits in inits3.txt
    initialize
    update 1000
    monitor N, thin(10)
    monitor alpha0, thin(10)
    monitor alpha1, thin(10)
    monitor alpha2, thin(10)
    monitor alpha3, thin(10)
    monitor alpha4, thin(10)
    monitor alpha5, thin(10)
    monitor alpha6, thin(10)
    monitor alpha7, thin(10)
    monitor alpha8, thin(10)
    monitor beta0, thin(10)
    monitor beta1, thin(10)
    monitor beta2, thin(10)
    monitor beta3, thin(10)
    monitor beta4, thin(10)
    monitor beta5, thin(10)
    monitor beta6, thin(10)
    monitor beta7, thin(10)
    monitor sigma, thin(10)
    monitor psi.sex, thin(10)
    monitor r, thin(10)
    monitor ro, thin(10)
    monitor s, thin(10)
    monitor z, thin(10)
    monitor sex, thin(10)
    update 3000
    coda *, stem(Chain3)
    ", fill = T)
sink()

cd Green_Lab/JAGS/KRFP_Scat/Scat_22Mar
SGE_Batch -c 'jags 1_Scat.cmd' -r KRFPScat.1 -q inr
SGE_Batch -c 'jags 2_Scat.cmd' -r KRFPScat.2 -q inr
SGE_Batch -c 'jags 3_Scat.cmd' -r KRFPScat.3 -q inr

# Read in
library(runjags)
library(jagsUI)
library(rjags)
library(R2WinBUGS)
library(MCMCglmm)
library(coda)
library(bayesplot)
library(MCMCvis)
library(rstan)
library(raster)
library(reshape2)
library(dplyr)
library(rgdal)
library(rgeos)

# Bring in Yosemite park footprint
setwd("C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/YosemiteCougar/Data")
yoseNP <- readOGR(".", "yose_boundaryNAD83Z11")

# Set wd for where the output is and bring it in
setwd("C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/YosemiteCougar/Scat_24May22_175M_22000K_9_3_Long")
c1 <- read.coda("Chain1chain1.txt", "Chain1index.txt")
c2 <- read.coda("Chain2chain1.txt", "Chain2index.txt")
c3 <- read.coda("Chain3chain1.txt", "Chain3index.txt")
fit <- mcmc.list(c1,c2,c3)

# Create a traceplot to examine the chains
MCMCtrace(fit, pdf = T)
fit1 <- fit
fit2 <- as.mcmc.list(list(c1,c2,c3))
fit3 <- as.matrix(fit2)
fit4 <- as.data.frame(fit3)
sum1 <- summary(fit1)
data <- read_rdump("data.txt")

# Output for the paper
statssum <- as.data.frame(sum1$statistics[1:30,])
quantsum <- as.data.frame(sum1$quantiles[1:30,])
totalsum <- cbind(statssum,quantsum)
totalsum <- round(totalsum, digits = 3)

write.csv(totalsum, "YosemiteCougar_ParamSummary_220927.csv")

# Significance probabilities
probs <- array(NA, dim = c(8, 2))
for(i in 1:nrow(fit4)){
  probs[1,1] <- sum(ifelse(fit4$beta1 < 0, 1, 0))/(nrow(fit4))
  probs[2,1] <- sum(ifelse(fit4$beta2 < 0, 1, 0))/(nrow(fit4))
  probs[3,1] <- sum(ifelse(fit4$beta3 < 0, 1, 0))/(nrow(fit4))
  probs[4,1] <- sum(ifelse(fit4$beta4 < 0, 1, 0))/(nrow(fit4))
  probs[5,1] <- sum(ifelse(fit4$beta5 < 0, 1, 0))/(nrow(fit4))
  probs[6,1] <- sum(ifelse(fit4$beta6 < 0, 1, 0))/(nrow(fit4))
  probs[7,1] <- sum(ifelse(fit4$beta7 < 0, 1, 0))/(nrow(fit4))
  probs[8,1] <- sum(ifelse(fit4$alpha1 < 0, 1, 0))/(nrow(fit4))
  probs[1,2] <- sum(ifelse(fit4$beta1 > 0, 1, 0))/(nrow(fit4))
  probs[2,2] <- sum(ifelse(fit4$beta2 > 0, 1, 0))/(nrow(fit4))
  probs[3,2] <- sum(ifelse(fit4$beta3 > 0, 1, 0))/(nrow(fit4))
  probs[4,2] <- sum(ifelse(fit4$beta4 > 0, 1, 0))/(nrow(fit4))
  probs[5,2] <- sum(ifelse(fit4$beta5 > 0, 1, 0))/(nrow(fit4))
  probs[6,2] <- sum(ifelse(fit4$beta6 > 0, 1, 0))/(nrow(fit4))
  probs[7,2] <- sum(ifelse(fit4$beta7 > 0, 1, 0))/(nrow(fit4))
  probs[8,2] <- sum(ifelse(fit4$alpha1 > 0, 1, 0))/(nrow(fit4))
}
probs <- data.frame(probs)
row.names(probs) <- c("beta1","beta2","beta3","beta4","beta5","beta6","beta7","alpha1")
colnames(probs) <- c("P < 0", "P > 0")
probs

write.csv(probs, "YosemiteCOugar_CovariatePercentProbs_220728.csv")
################################################
# Constructing a density map for discrete SS #
################################################

# Pull in the predicted sexes
sex1 <- cbind(fit4$`sex[1]`,fit4$`sex[2]`,fit4$`sex[3]`,fit4$`sex[4]`,fit4$`sex[5]`,fit4$`sex[6]`,fit4$`sex[7]`,fit4$`sex[8]`,fit4$`sex[9]`,fit4$`sex[10]`,fit4$`sex[11]`,fit4$`sex[12]`,fit4$`sex[13]`,fit4$`sex[14]`,fit4$`sex[15]`,fit4$`sex[16]`,fit4$`sex[17]`,fit4$`sex[18]`,fit4$`sex[19]`,fit4$`sex[20]`,fit4$`sex[21]`,fit4$`sex[22]`,fit4$`sex[23]`,fit4$`sex[24]`,fit4$`sex[25]`,fit4$`sex[26]`,fit4$`sex[27]`,fit4$`sex[28]`,fit4$`sex[29]`,fit4$`sex[30]`,fit4$`sex[31]`,fit4$`sex[32]`,fit4$`sex[33]`,fit4$`sex[34]`,fit4$`sex[35]`,fit4$`sex[36]`,fit4$`sex[37]`,fit4$`sex[38]`,fit4$`sex[39]`,fit4$`sex[40]`,fit4$`sex[41]`,fit4$`sex[42]`,fit4$`sex[43]`,fit4$`sex[44]`,fit4$`sex[45]`,fit4$`sex[46]`,fit4$`sex[47]`,fit4$`sex[48]`,fit4$`sex[49]`,fit4$`sex[50]`,fit4$`sex[51]`,fit4$`sex[52]`,fit4$`sex[53]`,fit4$`sex[54]`,fit4$`sex[55]`,fit4$`sex[56]`,fit4$`sex[57]`,fit4$`sex[58]`,fit4$`sex[59]`,fit4$`sex[60]`,fit4$`sex[61]`,fit4$`sex[62]`,fit4$`sex[63]`,fit4$`sex[64]`,fit4$`sex[65]`,fit4$`sex[66]`,fit4$`sex[67]`,fit4$`sex[68]`,fit4$`sex[69]`,fit4$`sex[70]`,fit4$`sex[71]`,fit4$`sex[72]`,fit4$`sex[73]`,fit4$`sex[74]`,fit4$`sex[75]`,fit4$`sex[76]`,fit4$`sex[77]`,fit4$`sex[78]`,fit4$`sex[79]`,fit4$`sex[80]`,
              fit4$`sex[81]`,fit4$`sex[82]`,fit4$`sex[83]`,fit4$`sex[84]`,fit4$`sex[85]`,fit4$`sex[86]`,fit4$`sex[87]`,fit4$`sex[88]`,fit4$`sex[89]`,fit4$`sex[90]`,
              fit4$`sex[91]`,fit4$`sex[92]`,fit4$`sex[93]`,fit4$`sex[94]`,fit4$`sex[95]`,fit4$`sex[96]`,fit4$`sex[97]`,fit4$`sex[98]`,fit4$`sex[99]`,
              fit4$`sex[100]`,fit4$`sex[101]`,fit4$`sex[102]`,fit4$`sex[103]`,fit4$`sex[104]`,fit4$`sex[105]`,fit4$`sex[106]`,fit4$`sex[107]`,fit4$`sex[108]`,fit4$`sex[109]`,
              fit4$`sex[110]`,fit4$`sex[111]`,fit4$`sex[112]`,fit4$`sex[113]`,fit4$`sex[114]`,fit4$`sex[115]`,fit4$`sex[116]`,fit4$`sex[117]`,fit4$`sex[118]`,fit4$`sex[119]`,
              fit4$`sex[120]`,fit4$`sex[121]`,fit4$`sex[122]`,fit4$`sex[123]`,fit4$`sex[124]`,fit4$`sex[125]`,fit4$`sex[126]`,fit4$`sex[127]`,fit4$`sex[128]`,fit4$`sex[129]`,
              fit4$`sex[130]`,fit4$`sex[131]`,fit4$`sex[132]`,fit4$`sex[133]`,fit4$`sex[134]`,fit4$`sex[135]`,fit4$`sex[136]`,fit4$`sex[137]`,fit4$`sex[138]`,fit4$`sex[139]`,
              fit4$`sex[140]`,fit4$`sex[141]`,fit4$`sex[142]`,fit4$`sex[143]`,fit4$`sex[144]`,fit4$`sex[145]`,fit4$`sex[146]`,fit4$`sex[147]`,fit4$`sex[148]`,fit4$`sex[149]`,
              fit4$`sex[150]`,fit4$`sex[151]`,fit4$`sex[152]`,fit4$`sex[153]`,fit4$`sex[154]`,fit4$`sex[155]`,fit4$`sex[156]`,fit4$`sex[157]`,fit4$`sex[158]`,fit4$`sex[159]`,
              fit4$`sex[160]`,fit4$`sex[161]`,fit4$`sex[162]`,fit4$`sex[163]`,fit4$`sex[164]`,fit4$`sex[165]`,fit4$`sex[166]`,fit4$`sex[167]`,fit4$`sex[168]`,fit4$`sex[169]`,
              fit4$`sex[170]`,fit4$`sex[171]`,fit4$`sex[172]`,fit4$`sex[173]`,fit4$`sex[174]`,fit4$`sex[175]`)
sex2 <- ifelse(sex1 == 1, 0, 1) 
# Sex1: 1 = males
# Sex2: 1 = females

# Pull in activity centers (s) and alive (z) for 1 year of estimate
t1.1<- cbind(fit4$`s[1]`,fit4$`s[2]`,fit4$`s[3]`,fit4$`s[4]`,fit4$`s[5]`,fit4$`s[6]`,fit4$`s[7]`,fit4$`s[8]`,fit4$`s[9]`,fit4$`s[10]`,fit4$`s[11]`,fit4$`s[12]`,fit4$`s[13]`,fit4$`s[14]`,fit4$`s[15]`,fit4$`s[16]`,fit4$`s[17]`,fit4$`s[18]`,fit4$`s[19]`,fit4$`s[20]`,fit4$`s[21]`,fit4$`s[22]`,fit4$`s[23]`,fit4$`s[24]`,fit4$`s[25]`,fit4$`s[26]`,fit4$`s[27]`,fit4$`s[28]`,fit4$`s[29]`,fit4$`s[30]`,fit4$`s[31]`,fit4$`s[32]`,fit4$`s[33]`,fit4$`s[34]`,fit4$`s[35]`,fit4$`s[36]`,fit4$`s[37]`,fit4$`s[38]`,fit4$`s[39]`,fit4$`s[40]`,fit4$`s[41]`,fit4$`s[42]`,fit4$`s[43]`,fit4$`s[44]`,fit4$`s[45]`,fit4$`s[46]`,fit4$`s[47]`,fit4$`s[48]`,fit4$`s[49]`,fit4$`s[50]`,fit4$`s[51]`,fit4$`s[52]`,fit4$`s[53]`,fit4$`s[54]`,fit4$`s[55]`,fit4$`s[56]`,fit4$`s[57]`,fit4$`s[58]`,fit4$`s[59]`,fit4$`s[60]`,fit4$`s[61]`,fit4$`s[62]`,fit4$`s[63]`,fit4$`s[64]`,fit4$`s[65]`,fit4$`s[66]`,fit4$`s[67]`,fit4$`s[68]`,fit4$`s[69]`,fit4$`s[70]`,fit4$`s[71]`,fit4$`s[72]`,fit4$`s[73]`,fit4$`s[74]`,fit4$`s[75]`,fit4$`s[76]`,fit4$`s[77]`,fit4$`s[78]`,fit4$`s[79]`,fit4$`s[80]`,
             fit4$`s[81]`,fit4$`s[82]`,fit4$`s[83]`,fit4$`s[84]`,fit4$`s[85]`,fit4$`s[86]`,fit4$`s[87]`,fit4$`s[88]`,fit4$`s[89]`,fit4$`s[90]`,
             fit4$`s[91]`,fit4$`s[92]`,fit4$`s[93]`,fit4$`s[94]`,fit4$`s[95]`,fit4$`s[96]`,fit4$`s[97]`,fit4$`s[98]`,fit4$`s[99]`,
             fit4$`s[100]`,fit4$`s[101]`,fit4$`s[102]`,fit4$`s[103]`,fit4$`s[104]`,fit4$`s[105]`,fit4$`s[106]`,fit4$`s[107]`,fit4$`s[108]`,fit4$`s[109]`,
             fit4$`s[110]`,fit4$`s[111]`,fit4$`s[112]`,fit4$`s[113]`,fit4$`s[114]`,fit4$`s[115]`,fit4$`s[116]`,fit4$`s[117]`,fit4$`s[118]`,fit4$`s[119]`,
             fit4$`s[120]`,fit4$`s[121]`,fit4$`s[122]`,fit4$`s[123]`,fit4$`s[124]`,fit4$`s[125]`,fit4$`s[126]`,fit4$`s[127]`,fit4$`s[128]`,fit4$`s[129]`,
             fit4$`s[130]`,fit4$`s[131]`,fit4$`s[132]`,fit4$`s[133]`,fit4$`s[134]`,fit4$`s[135]`,fit4$`s[136]`,fit4$`s[137]`,fit4$`s[138]`,fit4$`s[139]`,
             fit4$`s[140]`,fit4$`s[141]`,fit4$`s[142]`,fit4$`s[143]`,fit4$`s[144]`,fit4$`s[145]`,fit4$`s[146]`,fit4$`s[147]`,fit4$`s[148]`,fit4$`s[149]`,
             fit4$`s[150]`,fit4$`s[151]`,fit4$`s[152]`,fit4$`s[153]`,fit4$`s[154]`,fit4$`s[155]`,fit4$`s[156]`,fit4$`s[157]`,fit4$`s[158]`,fit4$`s[159]`,
             fit4$`s[160]`,fit4$`s[161]`,fit4$`s[162]`,fit4$`s[163]`,fit4$`s[164]`,fit4$`s[165]`,fit4$`s[166]`,fit4$`s[167]`,fit4$`s[168]`,fit4$`s[169]`,
             fit4$`s[170]`,fit4$`s[171]`,fit4$`s[172]`,fit4$`s[173]`,fit4$`s[174]`,fit4$`s[175]`)
Ssim1 <- t1.1
t1.z <- cbind(fit4$`z[1]`,fit4$`z[2]`,fit4$`z[3]`,fit4$`z[4]`,fit4$`z[5]`,fit4$`z[6]`,fit4$`z[7]`,fit4$`z[8]`,fit4$`z[9]`,fit4$`z[10]`,fit4$`z[11]`,fit4$`z[12]`,fit4$`z[13]`,fit4$`z[14]`,fit4$`z[15]`,fit4$`z[16]`,fit4$`z[17]`,fit4$`z[18]`,fit4$`z[19]`,fit4$`z[20]`,fit4$`z[21]`,fit4$`z[22]`,fit4$`z[23]`,fit4$`z[24]`,fit4$`z[25]`,fit4$`z[26]`,fit4$`z[27]`,fit4$`z[28]`,fit4$`z[29]`,fit4$`z[30]`,fit4$`z[31]`,fit4$`z[32]`,fit4$`z[33]`,fit4$`z[34]`,fit4$`z[35]`,fit4$`z[36]`,fit4$`z[37]`,fit4$`z[38]`,fit4$`z[39]`,fit4$`z[40]`,fit4$`z[41]`,fit4$`z[42]`,fit4$`z[43]`,fit4$`z[44]`,fit4$`z[45]`,fit4$`z[46]`,fit4$`z[47]`,fit4$`z[48]`,fit4$`z[49]`,fit4$`z[50]`,fit4$`z[51]`,fit4$`z[52]`,fit4$`z[53]`,fit4$`z[54]`,fit4$`z[55]`,fit4$`z[56]`,fit4$`z[57]`,fit4$`z[58]`,fit4$`z[59]`,fit4$`z[60]`,fit4$`z[61]`,fit4$`z[62]`,fit4$`z[63]`,fit4$`z[64]`,fit4$`z[65]`,fit4$`z[66]`,fit4$`z[67]`,fit4$`z[68]`,fit4$`z[69]`,fit4$`z[70]`,fit4$`z[71]`,fit4$`z[72]`,fit4$`z[73]`,fit4$`z[74]`,fit4$`z[75]`,fit4$`z[76]`,fit4$`z[77]`,fit4$`z[78]`,fit4$`z[79]`,fit4$`z[80]`,
              fit4$`z[81]`,fit4$`z[82]`,fit4$`z[83]`,fit4$`z[84]`,fit4$`z[85]`,fit4$`z[86]`,fit4$`z[87]`,fit4$`z[88]`,fit4$`z[89]`,fit4$`z[90]`,
              fit4$`z[91]`,fit4$`z[92]`,fit4$`z[93]`,fit4$`z[94]`,fit4$`z[95]`,fit4$`z[96]`,fit4$`z[97]`,fit4$`z[98]`,fit4$`z[99]`,
              fit4$`z[100]`,fit4$`z[101]`,fit4$`z[102]`,fit4$`z[103]`,fit4$`z[104]`,fit4$`z[105]`,fit4$`z[106]`,fit4$`z[107]`,fit4$`z[108]`,fit4$`z[109]`,
              fit4$`z[110]`,fit4$`z[111]`,fit4$`z[112]`,fit4$`z[113]`,fit4$`z[114]`,fit4$`z[115]`,fit4$`z[116]`,fit4$`z[117]`,fit4$`z[118]`,fit4$`z[119]`,
              fit4$`z[120]`,fit4$`z[121]`,fit4$`z[122]`,fit4$`z[123]`,fit4$`z[124]`,fit4$`z[125]`,fit4$`z[126]`,fit4$`z[127]`,fit4$`z[128]`,fit4$`z[129]`,
              fit4$`z[130]`,fit4$`z[131]`,fit4$`z[132]`,fit4$`z[133]`,fit4$`z[134]`,fit4$`z[135]`,fit4$`z[136]`,fit4$`z[137]`,fit4$`z[138]`,fit4$`z[139]`,
              fit4$`z[140]`,fit4$`z[141]`,fit4$`z[142]`,fit4$`z[143]`,fit4$`z[144]`,fit4$`z[145]`,fit4$`z[146]`,fit4$`z[147]`,fit4$`z[148]`,fit4$`z[149]`,
              fit4$`z[150]`,fit4$`z[151]`,fit4$`z[152]`,fit4$`z[153]`,fit4$`z[154]`,fit4$`z[155]`,fit4$`z[156]`,fit4$`z[157]`,fit4$`z[158]`,fit4$`z[159]`,
              fit4$`z[160]`,fit4$`z[161]`,fit4$`z[162]`,fit4$`z[163]`,fit4$`z[164]`,fit4$`z[165]`,fit4$`z[166]`,fit4$`z[167]`,fit4$`z[168]`,fit4$`z[169]`,
              fit4$`z[170]`,fit4$`z[171]`,fit4$`z[172]`,fit4$`z[173]`,fit4$`z[174]`,fit4$`z[175]`)

# Reorganize step 1
Sout1 <- Ssim1*t1.z
Sout1Sex <- Ssim1*t1.z*sex2    # Females
Sout1Sex.2 <- Ssim1*t1.z*sex1  # Males

# Reorganize step 2
post1 <- matrix(0, nrow=ncol(Sout1), data$GG+1)
dimnames(post1) <- list(1:nrow(post1), c(0,(1:data$GG)))
for(i in 1:nrow(post1)){
  xx <- table(Sout1[,i])
  idx <- names(xx)
  post1[i,idx] <- xx
}
post1 <- post1[,-1]

post1.males <- matrix(0, nrow=ncol(Sout1Sex.2), data$GG+1)
dimnames(post1.males) <- list(1:nrow(post1.males), c(0,(1:data$GG)))
for(i in 1:nrow(post1.males)){
  xx <- table(Sout1Sex.2[,i])
  idx <- names(xx)
  post1.males[i,idx] <- xx
}
post1.males <- post1.males[,-1]

post1.females <- matrix(0, nrow=ncol(Sout1Sex), data$GG+1)
dimnames(post1.females) <- list(1:nrow(post1.females), c(0,(1:data$GG)))
for(i in 1:nrow(post1.females)){
  xx <- table(Sout1Sex[,i])
  idx <- names(xx)
  post1.females[i,idx] <- xx
}
post1.females <- post1.females[,-1]

# Reorganize step 3
post.mean1 <- data.frame(apply(post1,2,mean))
post.mean1$GridID <- seq.int(nrow(post.mean1))
year1 <- cbind(post.mean1[,1], data$S)

year1dat <- as.data.frame(year1)
emp.hpd(year1[,1], conf = 0.95)
median(year1[,1])

post.meanmales <- data.frame(apply(post1.males,2,mean))
post.meanmales$GridID <- seq.int(nrow(post.meanmales))
year1.males <- cbind(post.meanmales[,1], data$S)

post.meanfemales <- data.frame(apply(post1.females,2,mean))
post.meanfemales$GridID <- seq.int(nrow(post.meanfemales))
year1.females <- cbind(post.meanfemales[,1], data$S)

# Double check
plot(data$S, cex = 0.5)
plot(yoseNP, add = T)

# Create the raster
StSp1_All <- rasterFromXYZ(cbind(year1[,2], year1[,3], year1[,1]))
StSp1_Females <- rasterFromXYZ(cbind(year1.females[,2], year1.females[,3], year1.females[,1]))
StSp1_Males <- rasterFromXYZ(cbind(year1.males[,2], year1.males[,3], year1.males[,1]))

# Plot the entire state space
jpeg("CougarDensityPlot_220927.jpg", width = 1000, height = 1000)
par(mar = c(5,5,2,1))
par(mfrow = c(2,2))
image(StSp1_Females, axes = 2, bty = "n", cex.axis = 1.5, cex.lab = 2, cex.main = 2, main = "Female cougars", ylab ="N-S UTM coordinates", xlab = "E-W UTM coordinates", ylim = c(StSp1_All@extent@ymin-1500,StSp1_All@extent@ymax+1500), xlim = c(StSp1_All@extent@xmin-1500,StSp1_All@extent@xmax+1500))
plot(yoseNP, add = T)
image(StSp1_Males, axes = 2, bty = "n", cex.axis = 1.5, cex.lab = 2, cex.main = 2, main = "Male cougars", ylab ="N-S UTM coordinates", xlab = "E-W UTM coordinates", ylim = c(StSp1_All@extent@ymin-1500,StSp1_All@extent@ymax+1500), xlim = c(StSp1_All@extent@xmin-1500,StSp1_All@extent@xmax+1500))
plot(yoseNP, add = T)
image(StSp1_All, axes = 2, bty = "n", cex.axis = 1.5, cex.lab = 2, cex.main = 2, main = "All cougars", ylab ="N-S UTM coordinates", xlab = "E-W UTM coordinates", ylim = c(StSp1_All@extent@ymin-1500,StSp1_All@extent@ymax+1500), xlim = c(StSp1_All@extent@xmin-1500,StSp1_All@extent@xmax+1500))
plot(yoseNP, add = T)
plot(0.5,0.5, xlim = c(0,1), ylim = c(0,1), type = "n", axes = F, xlab = "", ylab = "")
legend(-0.05,0.45, title = "Posterior density", cex = 2, bty = "n", fill = c('#abdda4','#e6f598','#fee08b','#fdae61','#f46d43','#d53e4f','#9e0142'), legend = c("<4", "4 - 5", "5 - 10", "10 - 15", "15 - 20",  ">20"))
# legend(-0.05,0.32, lty = c(1), lwd = c(2), legend = c("YNP boundary"), cex = 2, bty = "n")
dev.off()


###############################################
# Determine number of males/females in YNP #

# Create the Study Area (inside YNP)
gridSP <- SpatialPoints(data$S)
proj4string(gridSP) <- proj4string(yoseNP)
t <- over(gridSP, yoseNP)
t$GridID <- rownames(t)
insideSA <- t[!is.na(t$PARKNAME),]

# Create the empty arrays for individuals
count1 <- array(NA, dim = c(nrow(Sout1), ncol(Sout1))) 
count1.Sex <- array(NA, dim = c(nrow(Sout1), ncol(Sout1))) 
count1.Sex2 <- array(NA, dim = c(nrow(Sout1), ncol(Sout1))) 
test2 <- matrix(data = NA, nrow = nrow(Sout1), ncol = 1) # ncol is the number of years
test2.Sex <- matrix(data = NA, nrow = nrow(Sout1), ncol = 1) # ncol is the number of years
test2.Sex2 <- matrix(data = NA, nrow = nrow(Sout1), ncol = 1) # ncol is the number of years
both <- matrix(data = NA, nrow = 3, ncol = 5)
colnames(both) <- c("Mean", "SD", "2.5 CI", "Median", "97.5 CI")
rownames(both) <- c("NTot","NMales","NFemales")

library(TeachingDemos)
for(i in 1:nrow(Sout1)){ # number of iterations
  count1[i,] <- Sout1[i,] %in% insideSA$GridID
  test2[i,1] <- length(count1[i,][count1[i,] == TRUE])
  count1.Sex[i,] <- Sout1Sex.2[i,] %in% insideSA$GridID
  test2.Sex[i,1] <- length(count1.Sex[i,][count1.Sex[i,] == TRUE])
  count1.Sex2[i,] <- Sout1Sex[i,] %in% insideSA$GridID
  test2.Sex2[i,1] <- length(count1.Sex2[i,][count1.Sex2[i,] == TRUE])
}

for(t in 1:1){
  both[t,1] <- mean(test2[,t])
  both[t,2] <- sd(test2[,t])
  both[t,3] <- emp.hpd(test2[,t], conf = 0.95)[1]
  both[t,4] <- median(test2[,t])
  both[t,5] <- emp.hpd(test2[,t], conf = 0.95)[2]
  both[t+1,1] <- mean(test2.Sex[,1])
  both[t+1,2] <- sd(test2.Sex[,1])
  both[t+1,3] <- emp.hpd(test2.Sex[,1], conf = 0.95)[1]
  both[t+1,4] <- median(test2.Sex[,1])
  both[t+1,5] <- emp.hpd(test2.Sex[,1], conf = 0.95)[2]
  both[t+2,1] <- mean(test2.Sex2[,1])
  both[t+2,2] <- sd(test2.Sex2[,1])
  both[t+2,3] <- emp.hpd(test2.Sex2[,1], conf = 0.95)[1]
  both[t+2,4] <- median(test2.Sex2[,1])
  both[t+2,5] <- emp.hpd(test2.Sex2[,1], conf = 0.95)[2]
}

write.csv(both, file = "CougarAbundance.csv")


#### create average and variance of activity centers ####
gridcells <- as.data.frame(data$S)
names(gridcells) <- c("UTM_E", "UTM_N")
gridcells$GridCellCenter <- 1:nrow(gridcells)

allcougs <- as.data.frame(Sout1)
allcougs[allcougs==0] <- NA
allcougsmean <- summarize_if(allcougs, is.numeric, mean, na.rm = TRUE)
allcougsmean <- round(allcougsmean, digits=0)
allcougsvar <- summarize_if(allcougs, is.numeric, var, na.rm = TRUE)
allcougssd <- summarize_if(allcougs, is.numeric, sd, na.rm = TRUE)

allcougssum <- rbind(allcougsmean,
                     allcougssd,
                     allcougsvar)
allcougssum <- t(allcougssum)
allcougssum <- as.data.frame(allcougssum)
names(allcougssum) <- c("GridCellCenter", "Var", "SD")
allcougssum <- merge(allcougssum, gridcells[, c("GridCellCenter", "UTM_E", "UTM_N")], by="GridCellCenter")
allcougssum$ID <- 1:nrow(allcougssum)
allcougssum$Sex <- data$sex
write.csv(allcougssum, "allcougs_activitycentersum_220714.csv", row.names = FALSE)


##### females ####
femcougs <- as.data.frame(Sout1Sex)
femcougs[femcougs==0] <- NA
femcougsmean <- summarize_if(femcougs, is.numeric, mean, na.rm = TRUE)
femcougsmean <- round(femcougsmean, digits=0)
femcougsvar <- summarize_if(femcougs, is.numeric, var, na.rm = TRUE)
femcougssd <- summarize_if(femcougs, is.numeric, sd, na.rm = TRUE)

femcougssum <- rbind(femcougsmean,
                     femcougssd,
                     femcougsvar)
femcougssum <- t(femcougssum)
femcougssum <- as.data.frame(femcougssum)
names(femcougssum) <- c("GridCellCenter", "Var", "SD")
femcougssum <- merge(femcougssum, gridcells[, c("GridCellCenter", "UTM_E", "UTM_N")], by="GridCellCenter")
femcougssum$ID <- 1:nrow(femcougssum)
femcougssum$Sex <- data$sex
write.csv(femcougssum, "femcougs_activitycentersum_220720.csv", row.names = FALSE)


#### covariate relationships figure ####
#### occupancy figure ####
library(MCMCvis)
library(ggplot2)
library(viridis)
library(ggridges)
library(forcats)
library(cowplot)
library(ggpubr)
cougarsummary <- as.data.frame(MCMCsummary(fit,
                                           params = "all",
                                           Rhat = TRUE,
                                           n.eff = TRUE,
                                           probs = c(0.025, 0.1, 0.5, 0.9, 0.975),
                                           round = 2))

write.csv(cougarsummary, "YosemiteCougar_ParameterSummary_220728.csv")
cougarcovsummary <- MCMCsummary(fit,
                             params = c("beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta7"),
                             Rhat = TRUE,
                             n.eff = TRUE,
                             probs = c(0.025, 0.1, 0.5, 0.9, 0.975),
                             round = 2)

cougarcovsummary$Covariate <- c("Survey effort","Sex", "Year", "Slope", "Distance to trails", "Distance to roads", "Distance to rivers")
cougarcovsummary$Covariate = forcats::fct_rev(factor(cougarcovsummary$Covariate))
cougarcovsummary$overlap0 <- c("no", "yes", "no", "yes", "no", "yes", "yes")
cougarcovgraph <- ggplot(cougarcovsummary, aes(x = `50%`, y = factor(Covariate), fill = overlap0, color = overlap0)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(shape = 15, size = 3)+
  geom_errorbarh(aes(xmax = `97.5%`, xmin = `2.5%`, height = 0.0), size = 1.0)+
  geom_errorbarh(aes(xmax = `90%`, xmin = `10%`, height = 0.0), size = 2.0)+
  xlim(-3.5,3.5)+
  scale_color_manual(name = "", values = c("yes" = "#cccccc", 
                                           "no" = "#000000" 
                                           )) +
  theme(axis.text = element_text(size=18),
        axis.title= element_text(size = 22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x=element_text(size=18, margin = margin(2,2,2,2)),
        legend.position = "none",
        strip.background = element_rect(colour = "black", fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.line = element_line(colour = "black"))+
  labs(x="Effect on cougar detection probability", y = "")
tiff(filename = "YosemiteCOugar_Detection_220927.tiff", width = 8.5, height = 4.5, units = "in", res = 600, compression = 'lzw')
cougarcovgraph
dev.off()


