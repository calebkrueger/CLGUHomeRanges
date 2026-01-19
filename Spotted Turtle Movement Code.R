# Code to calculate home ranges, displacement, and relevant covariates
# Written by C. J. Krueger
# Last edited: 19-Jan-26

### Check that your telemetry data contain the following named columns:
### ID  
### Date
### Latitude
### Longitude

# Install and load packages

packages <- c(dplyr = "1.1.4",
              units = "1.0-0",
              glue = "1.8.0",
              curl = "7.0.0",
              sf = "1.0-23",
              terra = "1.8-86",
              adehabitatHR = "0.4.22",
              ctmm = "1.3.0",
              landscapemetrics = "2.2.1",
              daymetr = "1.7.1",
              spData = "2.3.4")

install.packages('remotes')
require(remotes)

# When the following code pauses and asks:
# "Enter one or more numbers, or an empty line to skip updates:"
# Type 1 then press enter to install all dependencies
# This can take a while since some packages have a lot of dependencies

mapply(remotes::install_version,
       package = names(packages),
       version = packages,
       MoreArgs = list(dependencies = T))

lapply(names(packages), require, character.only = T)

######################################################################

## STEP ONE: CALCULATE HOME RANGE SIZE ##

# Read in your spotted turtle telemetry data
# Make sure IDs are being treated as factors

data <- read.csv("./Telemetry Data.csv")
data$ID <- as.factor(data$ID)

# Specify date formatting
# Assumes dates are encoded mm/dd/yyyy (see 'format' argument)

data$Date <- as.POSIXct(data$Date, format = "%m/%d/%Y")

# Filter so there is only one observation per date
# Multiple observations per day with irregular sampling makes AKDE crash out

data %>%
  group_by(ID, Date) %>%
  slice(1) -> data

# Create ID-by-Year term for grouping observations to calculate annual home ranges

data$IDY <- interaction(data$ID, 
                        format(data$Date, "%Y"),
                        sep = "_",
                        drop = T)

# Convert lat-long data to an sf object
# Note that the 'crs' argument specifies that data coordinates are using
# the WGS84 datum (default for most handheld GPS units) and is NOT projected

sf.df <- st_as_sf(data,
                  coords = c("Longitude", "Latitude"),
                  crs = 4326)

# Calculate 100% MCP home range size
# Also calculates first and last days of tracking, tracking period, and # of points

# The following code calculates the appropriate UTM zone to reproject your data
# If data are from multiple sites across UTM zones, run them separately!

utm <- 32600 + (floor((mean(data$Longitude) + 180)/6)) + 1
sf.df <- st_transform(sf.df, crs = utm)

sf.df %>%
  group_by(ID, IDY) %>%
  summarise(IDY = unique(IDY),
            Start = format(min(Date), "%j"),
            Finish = format(max(Date), "%j"),
            Days = as.numeric(max(Date) - min(Date)),
            Points = length(geometry),
            geometry = st_combine(geometry)) %>%
  mutate(geometry = st_convex_hull(geometry)) %>%
  reframe(ID,
          IDY,
          Start,
          Finish,
          Days,
          Points,
          MCP100 = as.numeric(set_units(st_area(geometry), ha))) -> out

# Convert object to a SpatialPointsDataFrame for adehabitatHR
# Remove individuals with < 5 relocations in a tracking year
# adehabitatHR doesn't like them

sf.df[!sf.df$IDY %in% out[out$Points < 5,]$IDY,] -> ade.df
ade.df$IDY <- droplevels(ade.df$IDY)
ade.df <- as(ade.df[,"IDY"], "Spatial")

# Calculate 95% MCPs

as.data.frame(mcp(ade.df,
                  percent = 95,
                  unin = "m",
                  unout = "ha")) %>%
  reframe(IDY = id,
          MCP95 = area) %>%
  left_join(x = out, by = "IDY") -> out

# Calculate 50% MCPs

as.data.frame(mcp(ade.df,
                  percent = 50,
                  unin = "m",
                  unout = "ha")) %>%
  reframe(IDY = id,
          MCP50 = area) %>%
  left_join(x = out, by = "IDY") -> out

# Calculate 95% and 50% wAKDEc

# First, convert the data to a format that ctmm likes
# Remove individuals with < 3 points
# This is the lower technical limit for ctmm to estimate movement models
# Individuals will be filtered later on as well
# Project to same UTM Zone from above

data %>%
  filter(!IDY %in% out[out$Points < 3,]$IDY) %>%
  droplevels() %>%
  reframe(individual.local.identifier = IDY,
          timestamp = Date,
          location.long = Longitude,
          location.lat = Latitude) %>%
  as.telemetry(projection = utm) -> tel.df

# Calculate most appropriate movement model for each individual
# Also outputs empirical variograms for assessing home ranging behavior

mods <- list()
first.fits <- list()
fits <- list()
variograms <- list()

out$ctmm.mod <- NA

for(i in 1:length(tel.df)){
  variograms[[i]] <- variogram(tel.df[[i]],
                               dt = sort(unique(diff(tel.df[[i]]$timestamp)) %#% "day"),
                               fast = F,
                               trace = F,
                               CI = "Gauss")
  mods[[i]] <- ctmm.guess(tel.df[[i]],
                          variogram = variograms[[i]],
                          interactive = F)
  first.fits[[i]] <- ctmm.select(tel.df[[i]],
                                 mods[[i]],
                                 method = "pHREML",
                                 IC = "AICc",
                                 verbose = T,
                                 cores = 4)
  x <- summary(first.fits[[i]])
  # Perform parametric bootstrapping for individuals with small Neff
  if(x[1,3] < 5 & x[1,3] > 2.7){
    fits[[i]] <- ctmm.boot(tel.df[[i]],
                           first.fits[[i]][[1]],
                           error = 0.05,
                           iterate = T)
  } else {
    fits[[i]] <- first.fits[[i]][[1]]
  }
  out[out$IDY == tel.df[[i]]@info$identity,]$ctmm.mod <- summary(fits[[i]])$name
  # Print and visualize outputs
  print(variograms[[i]]@info$identity)
  print(summary(fits[[i]]))
  plot(variograms[[i]],
       CTMM = fits[[i]],
       level = c(0.5, 0.95),
       fraction = 1,
       main = variograms[[i]]@info$identity)
}

# Calculate wAKDEc home ranges
# Also output effective sample sizes for downstream filtering

akdehr <- list()
out$AKDE95 <- NA
out$AKDE50 <- NA
out$neff <- NA

for(i in 1:length(tel.df)){
  tmp.fit <- fits[[i]]
  akdehr[[i]] <- akde(tel.df[[i]],
                      tmp.fit,
                      debias = T,
                      weights = T)
  plot(tel.df[[i]],
       UD = akdehr[[i]],
       main = akdehr[[i]]@info$identity,
       error = F)
  out[out$IDY == akdehr[[i]]@info$identity,]$AKDE95 <- as.numeric(summary(akdehr[[i]], units = F)$CI[[2]]) * 0.0001
  out[out$IDY == akdehr[[i]]@info$identity,]$AKDE50 <- as.numeric(summary(akdehr[[i]], level.UD = 0.5, units = F)$CI[[2]]) * 0.0001
  out[out$IDY == akdehr[[i]]@info$identity,]$neff <- summary(akdehr[[i]])$DOF[[1]]
}

# Add info on site, sex, SCL, and mass
# Pulled from 'Extra Info.csv' file with columns:
# ID, Site, Year, Sex, SCL, and Mass

extra_info <- read.csv("./Extra Info.csv")
extra_info$ID <- as.factor(extra_info$ID)
extra_info$IDY <- interaction(extra_info$ID, 
                              extra_info$Year,
                              sep = "_",
                              drop = T)

extra_info %>%
  select(-ID) %>%
  left_join(x = out, y = ., by = "IDY") -> out

######################################################################

## STEP TWO: CALCULATE LANDSCAPE METRICS ##

# Create buffered capture points to crop NLCD raster
# 5 km is well beyond known movement distances of spotted turtles
# So no relevant habitat is being removed

data %>%
  arrange(Date) %>% 
  group_by(IDY) %>%
  st_as_sf(coords = c("Longitude", "Latitude"),
           crs = 4326) %>%
  st_transform(crs = utm) %>%
  st_buffer(dist = set_units(5, "km")) -> buffer.polys

# Download the Annual NLCD raster for each year in the dataset

# First, create list of years, download locations, and download urls
# Currently, 2025 needs to be dropped because data are not available yet

year <- unique(format(data$Date, "%Y"))
year <- as.list(year[year < 2025])

temp_zip <- lapply(seq_along(year), function(x){
  temp <- tempfile(fileext = ".zip")})

zip_url <- lapply(seq_along(year), function(i){
  temp <- glue("https://www.mrlc.gov/downloads/sciweb1/shared/mrlc/data-bundles/Annual_NLCD_LndCov_{year[[i]]}_CU_C1V1.zip")})

# Download the zip file containing the NLCD raster .tif for each year
# This will take at least a few minutes
# NOTE: you will need a couple Gb of space per raster

lapply(seq_along(year), function(i){
  curl_download(url = zip_url[[i]],
                destfile = temp_zip[[i]],
                handle = new_handle(timeout = 1e6))
  })

# Extract just the raster .tifs

lapply(seq_along(year), function(i){
  unzip(temp_zip[[i]],
        files = glue("Annual_NLCD_LndCov_{year[[i]]}_CU_C1V1.tif"),
        exdir = tempdir())
  })

# Read each annual raster into R

r <- lapply(seq_along(year), function(i){
  temp <- rast(file.path(tempdir(), glue("Annual_NLCD_LndCov_{year[[i]]}_CU_C1V1.tif")))
  })

# Make each raster name its year

names(r) <- year

# Finally, crop the raster to your 5km-buffered points

mask.polys <- st_transform(buffer.polys, crs(r[[1]]))
lc <- lapply(r, terra::crop, mask.polys, mask = T)

# Reclassify to end up with:
# 1 = vegetated wetlands (woody and herbaceous)
# 2 = open water (including perennial ice/snow)
# 3 = woody vegetation (including shrubs)
# 4 = open canopy vegetation (excluding crops and pastures)
# 5 = bare and cultivated land (including rock / clay / sand)
# 6 = developed land

lc <- lapply(lc, function(x){
  x[x %in% c(90,95)] <- 1
  return(x)})
lc <- lapply(lc, function(x){
  x[x %in% c(11,12)] <- 2
  return(x)})
lc <- lapply(lc, function(x){
  x[x %in% c(41,42,43,51,52)] <- 3
  return(x)})
lc <- lapply(lc, function(x){
  x[x %in% c(71,72,73,74)] <- 4
  return(x)})
lc <- lapply(lc, function(x){
  x[x %in% c(31,81,82)] <- 5
  return(x)})
lc <- lapply(lc, function(x){
  x[x %in% c(21,22,23,24)] <- 6
  return(x)})

lc <- lapply(lc, function(x){
  coltab(x) <- data.frame(value = c(1,2,3,4,5,6),
                          color = c('deepskyblue',
                                    'deepskyblue4',
                                    'darkseagreen4',
                                    'darkseagreen',
                                    'khaki',
                                    'gray90'))
  return(x)
})

# Plot the first annual raster to make sure it's looking right

terra::plot(lc[[1]])

# Download NWI data for your state
# Need to start by using coordinates to determine state code

states <- st_transform(spData::us_states, crs = crs(sf.df))
state <- state.abb[match(states[["NAME"]][as.integer(st_intersects(sf.df[1,], states))],state.name)]

temp_zip <- tempfile(fileext = ".zip")
zip_url <- glue("https://documentst.ecosphere.fws.gov/wetlands/data/State-Downloads/{state}_geodatabase_wetlands.zip")

curl_download(url = zip_url,
              destfile = temp_zip,
              handle = new_handle(timeout = 1e6))

unzip(temp_zip,
      exdir = tempdir())

# Read in the NWI vector data and convert to raster

wetlands <- terra::vect(paste(tempdir(), glue("/{state}_geodatabase_wetlands.gdb"), sep =""),
                        layer = glue("{state}_Wetlands"),
                        proxy = T)

wet.crop <- terra::project(query(wetlands, extent = terra::project(lc[[1]], crs(wetlands))), crs(lc[[1]]))
wet.rast <- rasterize(wet.crop,
                      disagg(lc[[1]], fact = 2),
                      field = "WETLAND_TYPE")

# Remove lakes, rivers/streams, and marine wetlands and deepwater from NWI wetland raster

nwi.all <- mask(wet.rast,
                wet.rast %in% c("Lake", "Riverine", "Estuarine and Marine Wetland", "Estuarine and Marine Deepwater"),
                maskvalue = 1)

nwi.open <- mask(nwi.all,
                 nwi.all %in% c("Freshwater Forested/Shrub Wetland"),
                 maskvalue = 1)

nwi.closed <- mask(nwi.all,
                   nwi.all %in% c("Freshwater Emergent Wetland", "Freshwater Pond", "Other"),
                   maskvalue = 1)

terra::plot(nwi.all)
terra::plot(nwi.open)
terra::plot(nwi.closed)

binary.nwi.all <- nwi.all / nwi.all
binary.nwi.open <- nwi.open / nwi.open
binary.nwi.closed <- nwi.closed / nwi.closed

# Create buffers within which we'll calculate landscape variables

# 250 meters
# Average home range length from published literature
# Also predicted wetland occupancy in Joyal et al. 2001
# Also ~ average interwetland movement observed in Milam & Melvin 2001, Beaudry et al. 2009

data %>%
  arrange(Date) %>% 
  group_by(IDY) %>%
  slice_head(n = 1) %>%
  st_as_sf(coords = c("Longitude", "Latitude"),
           crs = 4326) %>%
  st_transform(crs = crs(lc[[1]])) %>%
  st_buffer(dist = set_units(250, "m")) -> low.buffer

# 1 kilometer
# Maximum inter-wetland movement of Joyal et al. 2001, Beaudry et al. 2009

data %>%
  arrange(Date) %>% 
  group_by(IDY) %>%
  slice_head(n = 1) %>%
  st_as_sf(coords = c("Longitude", "Latitude"),
           crs = 4326) %>%
  st_transform(crs = crs(lc[[1]])) %>%
  st_buffer(dist = set_units(1, "km")) -> med.buffer

# 3 kilometers
# Maximum movement from Lassiter et al. 2024

data %>%
  arrange(Date) %>% 
  group_by(IDY) %>%
  slice_head(n = 1) %>%
  st_as_sf(coords = c("Longitude", "Latitude"),
           crs = 4326) %>%
  st_transform(crs = crs(lc[[1]])) %>%
  st_buffer(dist = set_units(3, "km")) -> hi.buffer

# Create object of initial capture points for each individual

data %>%
  arrange(Date) %>% 
  group_by(IDY) %>%
  slice_head(n = 1) %>%
  st_as_sf(coords = c("Longitude", "Latitude"),
           crs = 4326) %>%
  st_transform(crs = crs(lc[[1]])) -> pts

sf.df <- st_transform(sf.df, crs = crs(lc[[1]]))

# Create object with 100% MCP polygons

sf.df %>%
  group_by(IDY) %>%
  summarise(geometry = st_combine(geometry)) %>%
  mutate(geometry = st_convex_hull(geometry)) -> mcp.polys

# Create object with 95% AKDE polygons

akde.polys <- lapply(seq_along(akdehr), function(i){
  temp <- st_transform(as.sf(akdehr[[i]], level = 0.95), crs = crs(lc[[1]]))[2,]
})

akde.names <- lapply(seq_along(akde.polys), function(i){
  temp <- sub(" .*$", "", akde.polys[[i]]$name)
})

names(akde.polys) <- akde.names

# Create raster of disjunct wetland patches and count cell size to get wetland areas

lc.water <- lapply(seq_along(lc), function(i){
  temp <- classify(lc[[i]], cbind(c(2,3,4,5,6), c(NA,NA,NA,NA,NA)))
})
wetland.patches <- lapply(seq_along(lc.water), function(i){
  temp <- patches(lc.water[[i]], directions = 8)
})
wetland.cells <- lapply(seq_along(wetland.patches), function(i){
  temp <- freq(wetland.patches[[i]])
})
wetland.cells <- lapply(seq_along(wetland.cells), function(i){
  wetland.cells[[i]]$count <- wetland.cells[[i]]$count * ((30*30) / 1e4)
  wetland.cells[[i]]
})

names(lc.water) <- names(wetland.patches) <- names(wetland.cells) <- year

# Repeat for NWI rasters

nwi.all.wetland.patches <- patches(nwi.all, directions = 8)
nwi.all.wetland.cells <- freq(nwi.all.wetland.patches)
nwi.all.wetland.cells$count <- nwi.all.wetland.cells$count * ((15*15) / 1e4)

nwi.open.wetland.patches <- patches(nwi.open, directions = 8)
nwi.open.wetland.cells <- freq(nwi.open.wetland.patches)
nwi.open.wetland.cells$count <- nwi.open.wetland.cells$count * ((15*15) / 1e4)

nwi.closed.wetland.patches <- patches(nwi.closed, directions = 8)
nwi.closed.wetland.cells <- freq(nwi.closed.wetland.patches)
nwi.closed.wetland.cells$count <- nwi.closed.wetland.cells$count * ((15*15) / 1e4)

# Repeat for developed land

lc.dev <- lapply(seq_along(lc), function(i){
  temp <- classify(lc[[i]], cbind(c(1,2,3,4,5), c(NA,NA,NA,NA,NA)))
})
dev.patches <- lapply(seq_along(lc.dev), function(i){
  temp <- patches(lc.dev[[i]], directions = 8)
})
dev.cells <- lapply(seq_along(dev.patches), function(i){
  temp <- freq(dev.patches[[i]])
})
dev.cells <- lapply(seq_along(dev.cells), function(i){
  dev.cells[[i]]$count <- dev.cells[[i]]$count * ((30*30) / 1e4)
  dev.cells[[i]]
})

names(lc.dev) <- names(dev.patches) <- names(dev.cells) <- year

split_pts <- split(pts, pts$IDY)
split_sf.df <- split(sf.df, sf.df$IDY)
split_low.buffer <- split(low.buffer, low.buffer$IDY)
split_med.buffer <- split(med.buffer, med.buffer$IDY)
split_hi.buffer <- split(hi.buffer, hi.buffer$IDY)
split_mcp.polys <- split(mcp.polys, mcp.polys$IDY)

# Extract landscape variables of interest using each buffer level

for(i in as.character(unique(low.buffer$IDY))) {
  yr <- substr(i, nchar(i)-3, nchar(i))
  if(yr == 2025){
  } else {
    # Calculate distance from capture point to nearest wetland
    out[out$IDY==i, "wetland.dist"] <- min(terra::distance(x = vect(split_pts[[i]]),
                                                           y = as.polygons(wetland.patches[[yr]])))
    # Calculate size of wetland closest to capture point
    nearby(x = vect(split_pts[[i]]),
           y = as.polygons(wetland.patches[[yr]]),
           centroids = F)[2] -> idx
    tmp.wetland.cells <- wetland.cells[[yr]]
    out[out$IDY==i, "wetland.size"] <- tmp.wetland.cells[idx, 3]
    # Calculate maximum and net displacement
    dists <- st_distance(split_pts[[i]], split_sf.df[[i]])
    out[out$IDY==i, "dmax"] <- as.numeric(max(dists))
    out[out$IDY==i, "dnet"] <- as.numeric(dists[length(dists)])
    # Calculate proportion and number of wetland features around capture point, size of largest wetland
    tmp <- terra::extract(wetland.patches[[yr]], split_low.buffer[[i]], exact = T)
    out[out$IDY==i, "pwet250m"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    out[out$IDY==i, "nwet250m"] <- length(unique(na.omit(tmp)$patches))
    tmp.sizes <- tmp.wetland.cells[tmp.wetland.cells$value %in% unique(tmp$patches),] %>% arrange(-count)
    out[out$IDY==i, "wetland.size250m"] <- tmp.sizes[1,3]
    tmp <- terra::extract(wetland.patches[[yr]], split_med.buffer[[i]], exact = T)
    out[out$IDY==i, "pwet1km"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    out[out$IDY==i, "nwet1km"] <- length(unique(na.omit(tmp)$patches))
    tmp.sizes <- tmp.wetland.cells[tmp.wetland.cells$value %in% unique(tmp$patches),] %>% arrange(-count)
    out[out$IDY==i, "wetland.size1km"] <- tmp.sizes[1,3]
    tmp <- terra::extract(wetland.patches[[yr]], split_hi.buffer[[i]], exact = T)
    out[out$IDY==i, "pwet3km"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    out[out$IDY==i, "nwet3km"] <- length(unique(na.omit(tmp)$patches))
    tmp.sizes <- tmp.wetland.cells[tmp.wetland.cells$value %in% unique(tmp$patches),] %>% arrange(-count)
    out[out$IDY==i, "wetland.size3km"] <- tmp.sizes[1,3]
    # Calculate proportion and number of wetland features within home range areas
    if(out[out$IDY==i, "Points"] < 3) { } else {
      tmp <- terra::extract(wetland.patches[[yr]], split_mcp.polys[[i]], exact = T)
      out[out$IDY==i, "pwet.mcp"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
      out[out$IDY==i, "nwet.mcp"] <- length(unique(na.omit(tmp)$patches))
      tmp.sizes <- tmp.wetland.cells[tmp.wetland.cells$value %in% unique(tmp$patches),] %>% arrange(-count)
      out[out$IDY==i, "wetland.size.mcp"] <- tmp.sizes[1,3]
      tmp <- terra::extract(wetland.patches[[yr]], akde.polys[[i]], exact = T)
      out[out$IDY==i, "pwet.akde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
      out[out$IDY==i, "nwet.akde"] <- length(unique(na.omit(tmp)$patches))
      tmp.sizes <- tmp.wetland.cells[tmp.wetland.cells$value %in% unique(tmp$patches),] %>% arrange(-count)
      out[out$IDY==i, "wetland.size.akde"] <- tmp.sizes[1,3]
    }
    # Calculate proportion developed features around capture point
    tmp <- terra::extract(dev.patches[[yr]], split_low.buffer[[i]], exact = T)
    out[out$IDY==i, "pdev250m"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(dev.patches[[yr]], split_med.buffer[[i]], exact = T)
    out[out$IDY==i, "pdev1km"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(dev.patches[[yr]], split_hi.buffer[[i]], exact = T)
    out[out$IDY==i, "pdev3km"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    # Calculate proportion developed features within home range areas
    if(out[out$IDY==i, "Points"] < 3) { } else {
      tmp <- terra::extract(dev.patches[[yr]], split_mcp.polys[[i]], exact = T)
      out[out$IDY==i, "pdev.mcp"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
      tmp <- terra::extract(dev.patches[[yr]], akde.polys[[i]], exact = T)
      out[out$IDY==i, "pdev.akde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    }
    # Calculate wetland clumpiness index and cohesion within each buffer
    tmp <- terra::crop(lc[[yr]], split_low.buffer[[i]], touches = T, mask = T, snap = "out")
    tmp.clump <- lsm_c_clumpy(tmp)
    out[out$IDY==i, "clumpy250m"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
    tmp.coh <- lsm_c_cohesion(tmp)
    out[out$IDY==i, "cohes250m"] <- as.numeric(tmp.coh[tmp.coh$class == 1, 6])
    tmp <- terra::crop(lc[[yr]], split_med.buffer[[i]], touches = T, mask = T, snap = "out")
    tmp.clump <- lsm_c_clumpy(tmp)
    out[out$IDY==i, "clumpy1km"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
    tmp.coh <- lsm_c_cohesion(tmp)
    out[out$IDY==i, "cohes1km"] <- as.numeric(tmp.coh[tmp.coh$class == 1, 6])
    tmp <- terra::crop(lc[[yr]], split_hi.buffer[[i]], touches = T, mask = T, snap = "out")
    tmp.clump <- lsm_c_clumpy(tmp)
    out[out$IDY==i, "clumpy3km"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
    tmp.coh <- lsm_c_cohesion(tmp)
    out[out$IDY==i, "cohes3km"] <- as.numeric(tmp.coh[tmp.coh$class == 1, 6])
    # Calculate wetland clumpiness index and cohesion within each home range polygon
    if(out[out$IDY==i, "Points"] < 3) { } else {
      tmp <- terra::crop(lc[[yr]], split_mcp.polys[[i]], touches = T, mask = T, snap = "out")
      tmp.clump <- lsm_c_clumpy(tmp)
      out[out$IDY==i, "clumpy.mcp"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
      tmp.coh <- lsm_c_cohesion(tmp)
      out[out$IDY==i, "cohes.mcp"] <- as.numeric(tmp.coh[tmp.coh$class == 1, 6])
      tmp <- terra::crop(lc[[yr]], akde.polys[[i]], touches = T, mask = T, snap = "out")
      tmp.clump <- lsm_c_clumpy(tmp)
      out[out$IDY==i, "clumpy.akde"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
      tmp.coh <- lsm_c_cohesion(tmp)
      out[out$IDY==i, "cohes.akde"] <- as.numeric(tmp.coh[tmp.coh$class == 1, 6])
    }
  }
}

# Calculate landscape variables using NWI rasters instead

for(i in as.character(unique(low.buffer$IDY))) {
  # Calculate distance from capture point to nearest wetland
  out[out$IDY==i, "NWI.wetland.dist"] <- min(terra::distance(x = vect(split_pts[[i]]),
                                                         y = as.polygons(nwi.all.wetland.patches)))
  out[out$IDY==i, "NWI.open.dist"] <- min(terra::distance(x = vect(split_pts[[i]]),
                                                          y = as.polygons(nwi.open.wetland.patches)))
  out[out$IDY==i, "NWI.closed.dist"] <- min(terra::distance(x = vect(split_pts[[i]]),
                                                            y = as.polygons(nwi.closed.wetland.patches)))
  # Calculate size of wetland closest to capture point
  nearby(x = vect(split_pts[[i]]),
         y = as.polygons(nwi.all.wetland.patches),
         centroids = F)[2] -> idx
  out[out$IDY==i, "NWI.wetland.size"] <- nwi.all.wetland.cells[idx, 3]
  nearby(x = vect(split_pts[[i]]),
         y = as.polygons(nwi.open.wetland.patches),
         centroids = F)[2] -> idx
  out[out$IDY==i, "NWI.open.size"] <- nwi.open.wetland.cells[idx, 3]
  nearby(x = vect(split_pts[[i]]),
         y = as.polygons(nwi.closed.wetland.patches),
         centroids = F)[2] -> idx
  out[out$IDY==i, "NWI.closed.size"] <- nwi.closed.wetland.cells[idx, 3]
  # Calculate proportion and number of wetland features around capture point, size of largest wetland
  tmp <- terra::extract(nwi.all.wetland.patches, split_low.buffer[[i]], exact = T)
  out[out$IDY==i, "NWI.pwet250m"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
  out[out$IDY==i, "NWI.nwet250m"] <- length(unique(na.omit(tmp)$patches))
  tmp.sizes <- nwi.all.wetland.cells[nwi.all.wetland.cells$value %in% unique(tmp$patches),] %>% arrange(-count)
  out[out$IDY==i, "NWI.wetland.size250m"] <- tmp.sizes[1,3]
  tmp <- terra::extract(nwi.all.wetland.patches, split_med.buffer[[i]], exact = T)
  out[out$IDY==i, "NWI.pwet1km"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
  out[out$IDY==i, "NWI.nwet1km"] <- length(unique(na.omit(tmp)$patches))
  tmp.sizes <- nwi.all.wetland.cells[nwi.all.wetland.cells$value %in% unique(tmp$patches),] %>% arrange(-count)
  out[out$IDY==i, "NWI.wetland.size1km"] <- tmp.sizes[1,3]
  tmp <- terra::extract(nwi.all.wetland.patches, split_hi.buffer[[i]], exact = T)
  out[out$IDY==i, "NWI.pwet3km"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
  out[out$IDY==i, "NWI.nwet3km"] <- length(unique(na.omit(tmp)$patches))
  tmp.sizes <- nwi.all.wetland.cells[nwi.all.wetland.cells$value %in% unique(tmp$patches),] %>% arrange(-count)
  out[out$IDY==i, "NWI.wetland.size3km"] <- tmp.sizes[1,3]
  # Repeat with only open canopy wetlands
  tmp <- terra::extract(nwi.open.wetland.patches, split_low.buffer[[i]], exact = T)
  out[out$IDY==i, "NWI.popen250m"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
  out[out$IDY==i, "NWI.nopen250m"] <- length(unique(na.omit(tmp)$patches))
  tmp.sizes <- nwi.open.wetland.cells[nwi.open.wetland.cells$value %in% unique(tmp$patches),] %>% arrange(-count)
  out[out$IDY==i, "NWI.open.size250m"] <- tmp.sizes[1,3]
  tmp <- terra::extract(nwi.open.wetland.patches, split_med.buffer[[i]], exact = T)
  out[out$IDY==i, "NWI.popen1km"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
  out[out$IDY==i, "NWI.nopen1km"] <- length(unique(na.omit(tmp)$patches))
  tmp.sizes <- nwi.open.wetland.cells[nwi.open.wetland.cells$value %in% unique(tmp$patches),] %>% arrange(-count)
  out[out$IDY==i, "NWI.open.size1km"] <- tmp.sizes[1,3]
  tmp <- terra::extract(nwi.open.wetland.patches, split_hi.buffer[[i]], exact = T)
  out[out$IDY==i, "NWI.popen3km"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
  out[out$IDY==i, "NWI.nopen3km"] <- length(unique(na.omit(tmp)$patches))
  tmp.sizes <- nwi.open.wetland.cells[nwi.open.wetland.cells$value %in% unique(tmp$patches),] %>% arrange(-count)
  out[out$IDY==i, "NWI.open.size3km"] <- tmp.sizes[1,3]
  # Repeat again for closed canopy wetlands
  tmp <- terra::extract(nwi.closed.wetland.patches, split_low.buffer[[i]], exact = T)
  out[out$IDY==i, "NWI.pclosed250m"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
  out[out$IDY==i, "NWI.nclosed250m"] <- length(unique(na.omit(tmp)$patches))
  tmp.sizes <- nwi.closed.wetland.cells[nwi.closed.wetland.cells$value %in% unique(tmp$patches),] %>% arrange(-count)
  out[out$IDY==i, "NWI.closed.size250m"] <- tmp.sizes[1,3]
  tmp <- terra::extract(nwi.closed.wetland.patches, split_med.buffer[[i]], exact = T)
  out[out$IDY==i, "NWI.pclosed1km"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
  out[out$IDY==i, "NWI.nclosed1km"] <- length(unique(na.omit(tmp)$patches))
  tmp.sizes <- nwi.closed.wetland.cells[nwi.closed.wetland.cells$value %in% unique(tmp$patches),] %>% arrange(-count)
  out[out$IDY==i, "NWI.closed.size1km"] <- tmp.sizes[1,3]
  tmp <- terra::extract(nwi.closed.wetland.patches, split_hi.buffer[[i]], exact = T)
  out[out$IDY==i, "NWI.pclosed3km"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
  out[out$IDY==i, "NWI.nclosed3km"] <- length(unique(na.omit(tmp)$patches))
  tmp.sizes <- nwi.closed.wetland.cells[nwi.closed.wetland.cells$value %in% unique(tmp$patches),] %>% arrange(-count)
  out[out$IDY==i, "NWI.closed.size3km"] <- tmp.sizes[1,3]
  # Calculate proportion and number of wetland features within home range areas
  if(out[out$IDY==i, "Points"] < 3) { } else {
    tmp <- terra::extract(nwi.all.wetland.patches, split_mcp.polys[[i]], exact = T)
    out[out$IDY==i, "NWI.pwet.mcp"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    out[out$IDY==i, "NWI.nwet.mcp"] <- length(unique(na.omit(tmp)$patches))
    tmp.sizes <- nwi.all.wetland.cells[nwi.all.wetland.cells$value %in% unique(tmp$patches),] %>% arrange(-count)
    out[out$IDY==i, "NWI.wetland.size.mcp"] <- tmp.sizes[1,3]
    tmp <- terra::extract(nwi.all.wetland.patches, akde.polys[[i]], exact = T)
    out[out$IDY==i, "NWI.pwet.akde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    out[out$IDY==i, "NWI.nwet.akde"] <- length(unique(na.omit(tmp)$patches))
    tmp.sizes <- nwi.all.wetland.cells[nwi.all.wetland.cells$value %in% unique(tmp$patches),] %>% arrange(-count)
    out[out$IDY==i, "NWI.wetland.size.akde"] <- tmp.sizes[1,3]
  }
  # Repeat with open canopy wetlands
  if(out[out$IDY==i, "Points"] < 3) { } else {
    tmp <- terra::extract(nwi.open.wetland.patches, split_mcp.polys[[i]], exact = T)
    out[out$IDY==i, "NWI.popen.mcp"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    out[out$IDY==i, "NWI.nopen.mcp"] <- length(unique(na.omit(tmp)$patches))
    tmp.sizes <- nwi.open.wetland.cells[nwi.open.wetland.cells$value %in% unique(tmp$patches),] %>% arrange(-count)
    out[out$IDY==i, "NWI.open.size.mcp"] <- tmp.sizes[1,3]
    tmp <- terra::extract(nwi.open.wetland.patches, akde.polys[[i]], exact = T)
    out[out$IDY==i, "NWI.popen.akde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    out[out$IDY==i, "NWI.nopen.akde"] <- length(unique(na.omit(tmp)$patches))
    tmp.sizes <- nwi.open.wetland.cells[nwi.open.wetland.cells$value %in% unique(tmp$patches),] %>% arrange(-count)
    out[out$IDY==i, "NWI.open.size.akde"] <- tmp.sizes[1,3]
  }
  # Repeat with closed canopy wetlands
  if(out[out$IDY==i, "Points"] < 3) { } else {
    tmp <- terra::extract(nwi.closed.wetland.patches, split_mcp.polys[[i]], exact = T)
    out[out$IDY==i, "NWI.pclosed.mcp"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    out[out$IDY==i, "NWI.nclosed.mcp"] <- length(unique(na.omit(tmp)$patches))
    tmp.sizes <- nwi.closed.wetland.cells[nwi.closed.wetland.cells$value %in% unique(tmp$patches),] %>% arrange(-count)
    out[out$IDY==i, "NWI.closed.size.mcp"] <- tmp.sizes[1,3]
    tmp <- terra::extract(nwi.closed.wetland.patches, akde.polys[[i]], exact = T)
    out[out$IDY==i, "NWI.pclosed.akde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    out[out$IDY==i, "NWI.nclosed.akde"] <- length(unique(na.omit(tmp)$patches))
    tmp.sizes <- nwi.closed.wetland.cells[nwi.closed.wetland.cells$value %in% unique(tmp$patches),] %>% arrange(-count)
    out[out$IDY==i, "NWI.closed.size.akde"] <- tmp.sizes[1,3]
  }
  # Calculate wetland clumpiness index and cohesion within each buffer
  tmp <- terra::crop(binary.nwi.all, split_low.buffer[[i]], touches = T, mask = T, snap = "out")
  tmp.clump <- lsm_c_clumpy(tmp)
  out[out$IDY==i, "NWI.clumpy250m"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
  tmp.coh <- lsm_c_cohesion(tmp)
  out[out$IDY==i, "NWI.cohes250m"] <- as.numeric(tmp.coh[tmp.coh$class == 1, 6])
  tmp <- terra::crop(binary.nwi.all, split_med.buffer[[i]], touches = T, mask = T, snap = "out")
  tmp.clump <- lsm_c_clumpy(tmp)
  out[out$IDY==i, "NWI.clumpy1km"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
  tmp.coh <- lsm_c_cohesion(tmp)
  out[out$IDY==i, "NWI.cohes1km"] <- as.numeric(tmp.coh[tmp.coh$class == 1, 6])
  tmp <- terra::crop(binary.nwi.all, split_hi.buffer[[i]], touches = T, mask = T, snap = "out")
  tmp.clump <- lsm_c_clumpy(tmp)
  out[out$IDY==i, "NWI.clumpy3km"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
  tmp.coh <- lsm_c_cohesion(tmp)
  out[out$IDY==i, "NWI.cohes3km"] <- as.numeric(tmp.coh[tmp.coh$class == 1, 6])
  # Repeat with open canopy wetlands
  tmp <- terra::crop(binary.nwi.open, split_low.buffer[[i]], touches = T, mask = T, snap = "out")
  tmp.clump <- lsm_c_clumpy(tmp)
  out[out$IDY==i, "NWI.openclumpy250m"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
  tmp.coh <- lsm_c_cohesion(tmp)
  out[out$IDY==i, "NWI.opencohes250m"] <- as.numeric(tmp.coh[tmp.coh$class == 1, 6])
  tmp <- terra::crop(binary.nwi.open, split_med.buffer[[i]], touches = T, mask = T, snap = "out")
  tmp.clump <- lsm_c_clumpy(tmp)
  out[out$IDY==i, "NWI.openclumpy1km"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
  tmp.coh <- lsm_c_cohesion(tmp)
  out[out$IDY==i, "NWI.opencohes1km"] <- as.numeric(tmp.coh[tmp.coh$class == 1, 6])
  tmp <- terra::crop(binary.nwi.open, split_hi.buffer[[i]], touches = T, mask = T, snap = "out")
  tmp.clump <- lsm_c_clumpy(tmp)
  out[out$IDY==i, "NWI.openclumpy3km"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
  tmp.coh <- lsm_c_cohesion(tmp)
  out[out$IDY==i, "NWI.opencohes3km"] <- as.numeric(tmp.coh[tmp.coh$class == 1, 6])
  # Repeat with closed canopy wetlands
  tmp <- terra::crop(binary.nwi.closed, split_low.buffer[[i]], touches = T, mask = T, snap = "out")
  tmp.clump <- lsm_c_clumpy(tmp)
  out[out$IDY==i, "NWI.closedclumpy250m"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
  tmp.coh <- lsm_c_cohesion(tmp)
  out[out$IDY==i, "NWI.closedcohes250m"] <- as.numeric(tmp.coh[tmp.coh$class == 1, 6])
  tmp <- terra::crop(binary.nwi.closed, split_med.buffer[[i]], touches = T, mask = T, snap = "out")
  tmp.clump <- lsm_c_clumpy(tmp)
  out[out$IDY==i, "NWI.closedclumpy1km"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
  tmp.coh <- lsm_c_cohesion(tmp)
  out[out$IDY==i, "NWI.closedcohes1km"] <- as.numeric(tmp.coh[tmp.coh$class == 1, 6])
  tmp <- terra::crop(binary.nwi.closed, split_hi.buffer[[i]], touches = T, mask = T, snap = "out")
  tmp.clump <- lsm_c_clumpy(tmp)
  out[out$IDY==i, "NWI.closedclumpy3km"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
  tmp.coh <- lsm_c_cohesion(tmp)
  out[out$IDY==i, "NWI.closedcohes3km"] <- as.numeric(tmp.coh[tmp.coh$class == 1, 6])
  # Calculate wetland clumpiness index and cohesion within each home range polygon
  if(out[out$IDY==i, "Points"] < 3) { } else {
    tmp <- terra::crop(binary.nwi.all, split_mcp.polys[[i]], touches = T, mask = T, snap = "out")
    tmp.clump <- lsm_c_clumpy(tmp)
    out[out$IDY==i, "NWI.clumpy.mcp"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
    tmp.coh <- lsm_c_cohesion(tmp)
    out[out$IDY==i, "NWI.cohes.mcp"] <- as.numeric(tmp.coh[tmp.coh$class == 1, 6])
    tmp <- terra::crop(binary.nwi.all, akde.polys[[i]], touches = T, mask = T, snap = "out")
    tmp.clump <- lsm_c_clumpy(tmp)
    out[out$IDY==i, "NWI.clumpy.akde"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
    tmp.coh <- lsm_c_cohesion(tmp)
    out[out$IDY==i, "NWI.cohes.akde"] <- as.numeric(tmp.coh[tmp.coh$class == 1, 6])
  }
  # Repeat with open canopy wetlands
  if(out[out$IDY==i, "Points"] < 3) { } else {
    tmp <- terra::crop(binary.nwi.open, split_mcp.polys[[i]], touches = T, mask = T, snap = "out")
    tmp.clump <- lsm_c_clumpy(tmp)
    out[out$IDY==i, "NWI.openclumpy.mcp"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
    tmp.coh <- lsm_c_cohesion(tmp)
    out[out$IDY==i, "NWI.opencohes.mcp"] <- as.numeric(tmp.coh[tmp.coh$class == 1, 6])
    tmp <- terra::crop(binary.nwi.open, akde.polys[[i]], touches = T, mask = T, snap = "out")
    tmp.clump <- lsm_c_clumpy(tmp)
    out[out$IDY==i, "NWI.openclumpy.akde"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
    tmp.coh <- lsm_c_cohesion(tmp)
    out[out$IDY==i, "NWI.opencohes.akde"] <- as.numeric(tmp.coh[tmp.coh$class == 1, 6])
  }
  # Repeat with closed canopy wetlands
  if(out[out$IDY==i, "Points"] < 3) { } else {
    tmp <- terra::crop(binary.nwi.closed, split_mcp.polys[[i]], touches = T, mask = T, snap = "out")
    tmp.clump <- lsm_c_clumpy(tmp)
    out[out$IDY==i, "NWI.closedclumpy.mcp"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
    tmp.coh <- lsm_c_cohesion(tmp)
    out[out$IDY==i, "NWI.closedcohes.mcp"] <- as.numeric(tmp.coh[tmp.coh$class == 1, 6])
    tmp <- terra::crop(binary.nwi.closed, akde.polys[[i]], touches = T, mask = T, snap = "out")
    tmp.clump <- lsm_c_clumpy(tmp)
    out[out$IDY==i, "NWI.closedclumpy.akde"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
    tmp.coh <- lsm_c_cohesion(tmp)
    out[out$IDY==i, "NWI.closedcohes.akde"] <- as.numeric(tmp.coh[tmp.coh$class == 1, 6])
  }
}

######################################################################

## STEP THREE: CALCULATE PRECIPITATION AND TEMPERATURE METRICS ##

data %>%
  group_by(ID) %>%
  slice_head(n = 1) %>%
  left_join(y = extra_info, by = "ID") %>%
  group_by(Site) %>%
  summarise(latitude = mean(Latitude),
            longitude = mean(Longitude)) %>%
  write.table('coords.csv',
              sep = ",",
              row.names = F)

clim.df <- download_daymet_batch(file_location = 'coords.csv',
                                 start = 1980)

lapply(clim.df,
       function(x) {
         mutate(x$data,
                Site = x$site,
                Date = as.Date(paste(year, yday, sep = "-"), "%Y-%j")) %>%
           mutate(Month = format(Date, "%m"),
                  Year = as.numeric(format(Date, "%Y")),
                  Precip = prcp..mm.day.,
                  Temp = (tmax..deg.c.+tmin..deg.c.)/2) %>%
           select(Site, Date, Year, Month, Precip, Temp)
       }) %>%
  bind_rows() -> climate.df

# Add precipitation during water year

climate.df$WY <- ifelse(climate.df$Month %in% c("10","11","12"),
                        as.numeric(climate.df$Year) + 1,
                        as.numeric(climate.df$Year))

climate.df %>%
  group_by(Site, WY) %>%
  summarize(Precip = sum(Precip, na.rm = T)) %>%
  group_by(Site) %>%
  summarize(Av.Precip = mean(Precip, na.rm = T)) -> WY.av

climate.df %>%
  group_by(Site, WY) %>%
  summarize(Year = min(WY),
            WY.Precip = sum(Precip, na.rm = T)) %>%
  left_join(y = WY.av, by = "Site") %>%
  mutate(Rel.WY.Precip = WY.Precip / Av.Precip,
         Abs.WY.Precip = WY.Precip - Av.Precip) %>%
  select(Site, Year, WY.Precip, Rel.WY.Precip, Abs.WY.Precip) -> WY.precip

out %>%
  left_join(y = WY.precip, by = c("Site", "Year")) -> out

# Add mean temperature during calendar year

climate.df %>%
  group_by(Site, Year) %>%
  summarize(Temp = mean(Temp, na.rm = T)) %>%
  group_by(Site) %>%
  summarize(Av.Temp = mean(Temp, na.rm = T)) -> temp.av

climate.df %>%
  group_by(Site, Year) %>%
  summarize(Year = max(Year),
            Temp = mean(Temp, na.rm = T)) %>%
  left_join(y = temp.av, by = "Site") %>%
  mutate(Rel.Temp = Temp / Av.Temp,
         Abs.Temp = Temp - Av.Temp) %>%
  select(Site, Year, Temp, Rel.Temp, Abs.Temp) -> temp

out %>%
  left_join(y = temp, by = c("Site", "Year")) -> out

# Add new ID column that combines ID and Site name to ensure the IDs are unique

out %>%
  mutate(UniqueID = interaction(as.character(Site), 
                                ID,
                                sep = "_",
                                drop = T)) %>%
  select(UniqueID, IDY, ID, Site, Year, Sex, SCL, Mass, everything()) -> out

# Save variogram plots to pdf and out object to "output.csv"

write.csv(out, 
          'output.csv',
          row.names = F)

pdf("variograms.pdf")

par(mfrow = c(3,2))
lapply(seq_along(variograms), function(i){
  plot(variograms[[i]],
       CTMM = fits[[i]],
       level = c(0.5, 0.95),
       fraction = 1,
       main = variograms[[i]]@info$identity)
})

dev.off()

