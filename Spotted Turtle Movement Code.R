# Code to calculate home ranges, displacement distances, and environmental covariates
# Done using all data ("full"), data separated by year ("ann"), and by season ("seas")
# Written by C. J. Krueger
# Last edited: 17-Aug-26

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

install.packages("remotes")
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

## STEP ONE: CALCULATE HOME RANGE SIZES ##

# Read in your spotted turtle telemetry data
# Make sure IDs are being treated as factors

data <- read.csv("./Telemetry Data.csv")
data$ID <- as.factor(data$ID)

# Specify date formatting
# Assumes dates are encoded mm/dd/yyyy

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

# Create ID-by-Year-by-Season term for calculating seasonal home ranges
# First, make sure that December is properly grouped with the following January and February for Winter home ranges

data %>%
  mutate(Day = as.numeric(format(Date, "%d")),
         Month = format(Date, "%m"),
         Year = format(Date, "%Y")) -> data
  
data$Year <- ifelse(data$Month == "12",
                    as.numeric(data$Year) + 1,
                    as.numeric(data$Year))

data$Season <- ifelse(data$Month %in% c("12", "01", "02"),
                      "Winter",
                      ifelse(data$Month %in% c("03", "04") | data$Month == "05" & data$Day < 16,
                             "Pre-nesting",
                             ifelse(data$Month == "05" & data$Day > 15 | data$Month == "06" | data$Month == "07" & data$Day < 16,
                                    "Nesting",
                                    ifelse(data$Month == "07" & data$Day > 15 | data$Month == "08",
                                           "Summer",
                                           "Fall"))))

data$IDYS <- interaction(data$ID,
                         data$Year,
                         data$Season,
                         sep = "_",
                         drop = T)

# Convert telemetry data to an sf object
# Note that the 'crs' argument specifies that data coordinates are using the WGS84 datum (default for most handheld GPS units) and are NOT projected

sf.df <- st_as_sf(data,
                  coords = c("Longitude", "Latitude"),
                  crs = 4326)

# The following code calculates the appropriate UTM zone to project your data
# If data are from multiple sites across UTM zones, run them separately!

utm <- 32600 + (floor((mean(data$Longitude) + 180)/6)) + 1
sf.df <- st_transform(sf.df, crs = utm)

# Create object of (daily) movement distances
# Will be combined across datasets to potentially filter home range estimates later on

sf.df %>%
  arrange(ID, Date) %>%
  group_by(ID, Date) %>%
  slice(1) -> move.df

move.df %>%
  group_by(IDY) %>%
  mutate(movement = as.numeric(st_distance(geometry, lag(geometry), by_element = T)),
         days = as.numeric(difftime(Date, lag(Date), units = "days")),
         daily.movement = movement / days) %>%
  mutate(Date = format(Date, "%j"),
         date = (as.numeric(Date) + as.numeric(lag(Date)))/2) %>%
  as_tibble() %>%
  select(-geometry) %>%
  ungroup() -> move.df

write.csv(move.df, 
          'movement_output.csv',
          row.names = F)

# Calculate 100% MCP home range size at three different levels (full, annual, seasonal)
# Also calculates some other variables

sf.df %>%
  group_by(ID) %>%
  summarise(Start = format(min(Date), "%j"),
            Finish = format(max(Date), "%j"),
            Days = as.numeric(max(Date) - min(Date)),
            Points = length(geometry),
            UniquePoints = length(unique(geometry)),
            geometry = st_combine(geometry)) %>%
  mutate(geometry = st_convex_hull(geometry)) %>%
  reframe(ID,
          Start,
          Finish,
          Days,
          Points,
          UniquePoints,
          MCP100 = as.numeric(set_units(st_area(geometry), ha))) -> full.out

sf.df %>%
  group_by(ID, IDY) %>%
  summarise(IDY = unique(IDY),
            Start = format(min(Date), "%j"),
            Finish = format(max(Date), "%j"),
            Days = as.numeric(max(Date) - min(Date)),
            Points = length(geometry),
            UniquePoints = length(unique(geometry)),
            geometry = st_combine(geometry)) %>%
  mutate(geometry = st_convex_hull(geometry)) %>%
  reframe(ID,
          IDY,
          Start,
          Finish,
          Days,
          Points,
          UniquePoints,
          MCP100 = as.numeric(set_units(st_area(geometry), ha))) -> ann.out

sf.df %>%
  group_by(ID, IDYS) %>%
  summarise(IDY = tail(IDY, n = 1),
            IDYS = unique(IDYS),
            Season = unique(Season),
            Start = format(min(Date), "%j"),
            Finish = format(max(Date), "%j"),
            Days = as.numeric(max(Date) - min(Date)),
            Points = length(geometry),
            UniquePoints = length(unique(geometry)),
            geometry = st_combine(geometry)) %>%
  mutate(geometry = st_convex_hull(geometry)) %>%
  reframe(ID,
          IDY,
          IDYS,
          Season,
          Start,
          Finish,
          Days,
          Points,
          UniquePoints,
          MCP100 = as.numeric(set_units(st_area(geometry), ha))) -> seas.out

# Calculate full home ranges first
# Convert object to a SpatialPointsDataFrame for adehabitatHR
# Remove individuals with < 5 points
# adehabitatHR doesn't like them

sf.df[!sf.df$ID %in% full.out[full.out$Points < 5,]$ID,] -> ade.df
ade.df$ID <- droplevels(ade.df$ID)
ade.df <- as(ade.df[,"ID"], "Spatial")

# Calculate 95% MCPs

as.data.frame(mcp(ade.df,
                  percent = 95,
                  unin = "m",
                  unout = "ha")) %>%
  reframe(ID = id,
          MCP95 = area) %>%
  left_join(x = full.out, by = "ID") -> full.out

# Calculate 50% MCPs

as.data.frame(mcp(ade.df,
                  percent = 50,
                  unin = "m",
                  unout = "ha")) %>%
  reframe(ID = id,
          MCP50 = area) %>%
  left_join(x = full.out, by = "ID") -> full.out

# Calculate 95% and 50% KDE and AKDE

# First, convert the data to a format that ctmm likes
# Remove individuals with < 3 unique points
# This is the lower technical limit for ctmm to estimate movement models
# Individuals will be filtered later on as well
# Project to same UTM Zone from above

data %>%
  filter(!ID %in% full.out[full.out$UniquePoints < 3,]$ID) %>%
  droplevels() %>%
  reframe(individual.local.identifier = ID,
          timestamp = Date,
          location.long = Longitude,
          location.lat = Latitude) %>%
  as.telemetry(projection = utm) -> tel.df

# Calculate most appropriate movement model for each individual
# Also outputs empirical variograms for assessing home ranging behavior

mods <- list()
first.fits <- list()
akde.fits <- list()
kde.fits <- list()
variograms <- list()

full.out$ctmm.mod <- NA

for(i in 1:length(tel.df)){
  # Fit empirical variogram
  variograms[[i]] <- variogram(tel.df[[i]],
                               fast = F,
                               trace = F,
                               CI = "Gauss")
  # Force ctmm to parameterize an OU model
  mods[[i]] <- ctmm.guess(tel.df[[i]],
                          variogram = variograms[[i]],
                          CTMM = ctmm(tau = c(1 %#% "day")),
                          interactive = F)
  # Select best-fit model
  first.fits[[i]] <- ctmm.select(tel.df[[i]],
                                 mods[[i]],
                                 method = "pHREML",
                                 IC = "AICc",
                                 verbose = T,
                                 cores = 4)
  # Force ctmm to fit IID model to estimate KDE
  iid <- ctmm.fit(tel.df[[i]],
                  ctmm())
  iid.fits <- ctmm.select(tel.df[[i]],
                          iid,
                          method = "pHREML",
                          IC = "AICc",
                          verbose = T,
                          cores = 4)
  # Save best-fit IID and non-IID models to estimate KDE and AKDE home ranges later
  x <- summary(first.fits[[i]])
  kde.fits[[i]] <- iid.fits[[1]]
  akde.fits[[i]] <- first.fits[[i]][[min(which(substr(rownames(x),1,2) != "II"))]]
  # Bootstrap non-IID model if Neff is between 2.7 and 5
  if(between(x[min(which(substr(rownames(x),1,2) != "II")),3], 2.7, 5)) {
    akde.fits[[i]] <- ctmm.boot(tel.df[[i]],
                                first.fits[[i]][[min(which(substr(rownames(x),1,2) != "II"))]],
                                error = 0.05,
                                iterate = T)
  } else {}
  full.out[full.out$ID == tel.df[[i]]@info$identity,]$ctmm.mod <- summary(akde.fits[[i]])$name
  # Print and visualize outputs
  print(variograms[[i]]@info$identity)
  print(summary(kde.fits[[i]]))
  print(summary(akde.fits[[i]]))
  # plot(variograms[[i]],
  #      CTMM = first.fits[[i]][1:3],
  #      level = c(0.5, 0.95),
  #      fraction = 1,
  #      main = variograms[[i]]@info$identity,
  #      col.CTMM = c("red","blue","black"))
}

# Calculate KDE and AKDE home ranges
# Also output effective sample sizes for downstream filtering

akdehr <- list()
kdehr <- list()
full.out$AKDE95 <- NA
full.out$AKDE50 <- NA
full.out$neff <- NA
full.out$KDE95 <- NA
full.out$KDE50 <- NA
full.out$kdeneff <- NA

names(akde.fits) <- names(tel.df) <- sapply(akde.fits, function(x){x@info$identity})

akde.fits <- akde.fits[sapply(akde.fits, function(x) {summary(x)$DOF[2]}) > 0.1]
tmp.tel.df <- tel.df[names(tel.df) %in% names(akde.fits)]

akdehr <- akde(tmp.tel.df,
               akde.fits,
               debias = T,
               weights = T,
               grid = list(dr = c(1,1)))

kdehr <- akde(tel.df,
              kde.fits,
              debias = T,
              weights = T,
              grid = list(dr = c(1,1)))

for(i in 1:length(akdehr)){
  # plot(tel.df[[i]],
  #      UD = akdehr[[i]],
  #      main = akdehr[[i]]@info$identity,
  #      error = F)
  full.out[full.out$ID == akdehr[[i]]@info$identity,]$AKDE95 <- as.numeric(summary(akdehr[[i]], units = F)$CI[[2]]) * 0.0001
  full.out[full.out$ID == akdehr[[i]]@info$identity,]$AKDE50 <- as.numeric(summary(akdehr[[i]], level.UD = 0.5, units = F)$CI[[2]]) * 0.0001
  full.out[full.out$ID == akdehr[[i]]@info$identity,]$neff <- summary(akdehr[[i]])$DOF[[1]]
}

for(i in 1:length(kdehr)){
  # plot(tel.df[[i]],
  #      UD = kdehr[[i]],
  #      main = kdehr[[i]]@info$identity,
  #      error = F)
  full.out[full.out$ID == kdehr[[i]]@info$identity,]$KDE95 <- as.numeric(summary(kdehr[[i]], units = F)$CI[[2]]) * 0.0001
  full.out[full.out$ID == kdehr[[i]]@info$identity,]$KDE50 <- as.numeric(summary(kdehr[[i]], level.UD = 0.5, units = F)$CI[[2]]) * 0.0001
  full.out[full.out$ID == kdehr[[i]]@info$identity,]$kdeneff <- summary(kdehr[[i]])$DOF[[1]]
}

# Add info on site, sex, SCL, and mass
# Pulled from 'Extra Info.csv' file with columns:
# ID, Site, Year, Sex, SCL, and Mass

extra_info <- read.csv("./Extra Info.csv")
extra_info$ID <- as.factor(extra_info$ID)

extra_info %>%
  group_by(ID) %>%
  slice_head(n = 1) %>%
  left_join(x = full.out, y = ., by = "ID") -> full.out

rm(extra_info, first.fits, iid, iid.fits, mods)
gc()

### Repeat for annual home ranges ###

sf.df[!sf.df$IDY %in% ann.out[ann.out$Points < 5,]$IDY,] -> ade.df
ade.df$IDY <- droplevels(ade.df$IDY)
ade.df <- as(ade.df[,"IDY"], "Spatial")

as.data.frame(mcp(ade.df,
                  percent = 95,
                  unin = "m",
                  unout = "ha")) %>%
  reframe(IDY = id,
          MCP95 = area) %>%
  left_join(x = ann.out, by = "IDY") -> ann.out

as.data.frame(mcp(ade.df,
                  percent = 50,
                  unin = "m",
                  unout = "ha")) %>%
  reframe(IDY = id,
          MCP50 = area) %>%
  left_join(x = ann.out, by = "IDY") -> ann.out

data %>%
  filter(!IDY %in% ann.out[ann.out$UniquePoints < 3,]$IDY) %>%
  droplevels() %>%
  reframe(individual.local.identifier = IDY,
          timestamp = Date,
          location.long = Longitude,
          location.lat = Latitude) %>%
  as.telemetry(projection = utm) -> ann.tel.df

ann.mods <- list()
ann.first.fits <- list()
ann.akde.fits <- list()
ann.kde.fits <- list()
ann.variograms <- list()

ann.out$ctmm.mod <- NA

for(i in 1:length(ann.tel.df)){
  ann.variograms[[i]] <- variogram(ann.tel.df[[i]],
                                   fast = F,
                                   trace = F,
                                   CI = "Gauss")
  ann.mods[[i]] <- ctmm.guess(ann.tel.df[[i]],
                              variogram = ann.variograms[[i]],
                              CTMM = ctmm(tau = c(1 %#% "day")),
                              interactive = F)
  ann.first.fits[[i]] <- ctmm.select(ann.tel.df[[i]],
                                     ann.mods[[i]],
                                     method = "pHREML",
                                     IC = "AICc",
                                     verbose = T,
                                     cores = 4)
  ann.iid <- ctmm.fit(ann.tel.df[[i]],
                      ctmm())
  ann.iid.fits <- ctmm.select(ann.tel.df[[i]],
                              ann.iid,
                              method = "pHREML",
                              IC = "AICc",
                              verbose = T,
                              cores = 4)
  x <- summary(ann.first.fits[[i]])
  ann.kde.fits[[i]] <- ann.iid.fits[[1]]
  ann.akde.fits[[i]] <- ann.first.fits[[i]][[min(which(substr(rownames(x),1,2) != "II"))]]
  if(between(x[min(which(substr(rownames(x),1,2) != "II")),3], 2.7, 5)) {
    ann.akde.fits[[i]] <- ctmm.boot(ann.tel.df[[i]],
                                    ann.first.fits[[i]][[min(which(substr(rownames(x),1,2) != "II"))]],
                                    error = 0.05,
                                    iterate = T)
  } else {}
  ann.out[ann.out$IDY == ann.tel.df[[i]]@info$identity,]$ctmm.mod <- summary(ann.akde.fits[[i]])$name
  print(ann.variograms[[i]]@info$identity)
  print(summary(ann.kde.fits[[i]]))
  print(summary(ann.akde.fits[[i]]))
  # plot(ann.variograms[[i]],
  #      CTMM = ann.first.fits[[i]][1:3],
  #      level = c(0.5, 0.95),
  #      fraction = 1,
  #      main = ann.variograms[[i]]@info$identity,
  #      col.CTMM = c("red","blue","black"))
}

# Calculate wAKDEc home ranges
# Also output effective sample sizes for downstream filtering

ann.akdehr <- list()
ann.kdehr <- list()
ann.out$AKDE95 <- NA
ann.out$AKDE50 <- NA
ann.out$neff <- NA
ann.out$KDE95 <- NA
ann.out$KDE50 <- NA
ann.out$kdeneff <- NA

names(ann.akde.fits) <- names(ann.tel.df) <- sapply(ann.akde.fits, function(x){x@info$identity})

ann.akde.fits <- ann.akde.fits[sapply(ann.akde.fits, function(x) {summary(x)$DOF[2]}) > 0.1]
tmp.tel.df <- ann.tel.df[names(ann.tel.df) %in% names(ann.akde.fits)]

ids <- unique(sub("_[0-9]{4}", "", names(ann.akde.fits)))

for(i in 1:length(ids)){
  ann.akdehr[[i]] <- akde(tmp.tel.df[sub("_[0-9]{4}", "", names(ann.akde.fits)) == ids[i]],
                          ann.akde.fits[sub("_[0-9]{4}", "", names(ann.akde.fits)) == ids[i]],
                          debias = T,
                          weights = T,
                          grid = list(dr = c(1,1)))
}

ann.akdehr <- do.call(c, ann.akdehr)

names(ann.kde.fits) <- names(ann.tel.df) <- sapply(ann.kde.fits, function(x){x@info$identity})

ids <- unique(sub("_[0-9]{4}", "", names(ann.kde.fits)))

for(i in 1:length(ids)){
  ann.kdehr[[i]] <- akde(ann.tel.df[sub("_[0-9]{4}", "", names(ann.kde.fits)) == ids[i]],
                         ann.kde.fits[sub("_[0-9]{4}", "", names(ann.kde.fits)) == ids[i]],
                         debias = T,
                         weights = T,
                         grid = list(dr = c(1,1)))
}

ann.kdehr <- do.call(c, ann.kdehr)

for(i in 1:length(ann.akdehr)){
  # plot(ann.tel.df[names(ann.akdehr[i])],
  #      UD = ann.akdehr[[i]],
  #      main = ann.akdehr[[i]]@info$identity,
  #      error = F)
  ann.out[ann.out$IDY == ann.akdehr[[i]]@info$identity,]$AKDE95 <- as.numeric(summary(ann.akdehr[[i]], units = F)$CI[[2]]) * 0.0001
  ann.out[ann.out$IDY == ann.akdehr[[i]]@info$identity,]$AKDE50 <- as.numeric(summary(ann.akdehr[[i]], level.UD = 0.5, units = F)$CI[[2]]) * 0.0001
  ann.out[ann.out$IDY == ann.akdehr[[i]]@info$identity,]$neff <- summary(ann.akdehr[[i]])$DOF[[1]]
}

for(i in 1:length(ann.kdehr)){
  # plot(ann.tel.df[names(ann.kdehr[i])],
  #      UD = ann.kdehr[[i]],
  #      main = ann.kdehr[[i]]@info$identity,
  #      error = F)
  ann.out[ann.out$IDY == ann.kdehr[[i]]@info$identity,]$KDE95 <- as.numeric(summary(ann.kdehr[[i]], units = F)$CI[[2]]) * 0.0001
  ann.out[ann.out$IDY == ann.kdehr[[i]]@info$identity,]$KDE50 <- as.numeric(summary(ann.kdehr[[i]], level.UD = 0.5, units = F)$CI[[2]]) * 0.0001
  ann.out[ann.out$IDY == ann.kdehr[[i]]@info$identity,]$kdeneff <- summary(ann.kdehr[[i]])$DOF[[1]]
}

extra_info <- read.csv("./Extra Info.csv")
extra_info$ID <- as.factor(extra_info$ID)
extra_info$IDY <- interaction(extra_info$ID, 
                              extra_info$Year,
                              sep = "_",
                              drop = T)

extra_info %>%
  select(-ID) %>%
  left_join(x = ann.out, y = ., by = "IDY") -> ann.out

rm(extra_info, ann.first.fits, ann.iid, ann.iid.fits, ann.mods)
gc()

### Repeat for seasonal home ranges ###

sf.df[!sf.df$IDYS %in% seas.out[seas.out$Points < 5,]$IDYS,] -> ade.df
ade.df$IDYS <- droplevels(ade.df$IDYS)
ade.df <- as(ade.df[,"IDYS"], "Spatial")

as.data.frame(mcp(ade.df,
                  percent = 95,
                  unin = "m",
                  unout = "ha")) %>%
  reframe(IDYS = id,
          MCP95 = area) %>%
  left_join(x = seas.out, by = "IDYS") -> seas.out

as.data.frame(mcp(ade.df,
                  percent = 50,
                  unin = "m",
                  unout = "ha")) %>%
  reframe(IDYS = id,
          MCP50 = area) %>%
  left_join(x = seas.out, by = "IDYS") -> seas.out

data %>%
  filter(!IDYS %in% seas.out[seas.out$UniquePoints < 3,]$IDYS) %>%
  droplevels() %>%
  reframe(individual.local.identifier = IDYS,
          timestamp = Date,
          location.long = Longitude,
          location.lat = Latitude) %>%
  as.telemetry(projection = utm) -> seas.tel.df

seas.mods <- list()
seas.first.fits <- list()
seas.akde.fits <- list()
seas.kde.fits <- list()
seas.variograms <- list()

seas.out$ctmm.mod <- NA

for(i in 1:length(seas.tel.df)){
  seas.variograms[[i]] <- variogram(seas.tel.df[[i]],
                                   fast = F,
                                   trace = F,
                                   CI = "Gauss")
  seas.mods[[i]] <- ctmm.guess(seas.tel.df[[i]],
                              variogram = seas.variograms[[i]],
                              CTMM = ctmm(tau = c(1 %#% "day")),
                              interactive = F)
  seas.first.fits[[i]] <- ctmm.select(seas.tel.df[[i]],
                                     seas.mods[[i]],
                                     method = "pHREML",
                                     IC = "AICc",
                                     verbose = T,
                                     cores = 4)
  seas.iid <- ctmm.fit(seas.tel.df[[i]],
                      ctmm())
  seas.iid.fits <- ctmm.select(seas.tel.df[[i]],
                              seas.iid,
                              method = "pHREML",
                              IC = "AICc",
                              verbose = T,
                              cores = 4)
  x <- summary(seas.first.fits[[i]])
  seas.kde.fits[[i]] <- seas.iid.fits[[1]]
  seas.akde.fits[[i]] <- seas.first.fits[[i]][[min(which(substr(rownames(x),1,2) != "II"))]]
  if(between(x[min(which(substr(rownames(x),1,2) != "II")),3], 2.7, 5)) {
    seas.akde.fits[[i]] <- ctmm.boot(seas.tel.df[[i]],
                                    seas.first.fits[[i]][[min(which(substr(rownames(x),1,2) != "II"))]],
                                    error = 0.05,
                                    iterate = T)
  } else {}
  seas.out[seas.out$IDYS == seas.tel.df[[i]]@info$identity,]$ctmm.mod <- summary(seas.akde.fits[[i]])$name
  print(seas.variograms[[i]]@info$identity)
  print(summary(seas.kde.fits[[i]]))
  print(summary(seas.akde.fits[[i]]))
  # plot(seas.variograms[[i]],
  #      CTMM = seas.first.fits[[i]][1:3],
  #      level = c(0.5, 0.95),
  #      fraction = 1,
  #      main = seas.variograms[[i]]@info$identity,
  #      col.CTMM = c("red","blue","black"))
}

seas.akdehr <- list()
seas.kdehr <- list()
seas.out$AKDE95 <- NA
seas.out$AKDE50 <- NA
seas.out$neff <- NA
seas.out$KDE95 <- NA
seas.out$KDE50 <- NA
seas.out$kdeneff <- NA

names(seas.akde.fits) <- names(seas.tel.df) <- sapply(seas.akde.fits, function(x){x@info$identity})

seas.akde.fits <- seas.akde.fits[sapply(seas.akde.fits, function(x) {summary(x)$DOF[2]}) > 0.1]
tmp.tel.df <- seas.tel.df[names(seas.tel.df) %in% names(seas.akde.fits)]

ids <- unique(sub("_[0-9]{4}_[^_]+$", "", names(seas.akde.fits)))

for(i in 1:length(ids)){
  seas.akdehr[[i]] <- akde(tmp.tel.df[sub("_[0-9]{4}_[^_]+$", "", names(seas.akde.fits)) == ids[i]],
                           seas.akde.fits[sub("_[0-9]{4}_[^_]+$", "", names(seas.akde.fits)) == ids[i]],
                           debias = T,
                           weights = T,
                           grid = list(dr = c(1,1)))
}

seas.akdehr <- do.call(c, seas.akdehr)

names(seas.kde.fits) <- names(seas.tel.df) <- sapply(seas.kde.fits, function(x){x@info$identity})

seas.kde.fits <- seas.kde.fits[sapply(seas.kde.fits, function(x) {summary(x)$DOF[2]}) > 0.1]
tmp.tel.df <- seas.tel.df[names(seas.tel.df) %in% names(seas.kde.fits)]

ids <- unique(sub("_[0-9]{4}_[^_]+$", "", names(seas.kde.fits)))

for(i in 1:length(ids)){
  seas.kdehr[[i]] <- akde(tmp.tel.df[sub("_[0-9]{4}_[^_]+$", "", names(seas.kde.fits)) == ids[i]],
                          seas.kde.fits[sub("_[0-9]{4}_[^_]+$", "", names(seas.kde.fits)) == ids[i]],
                          debias = T,
                          weights = T,
                          grid = list(dr = c(1,1)))
}

seas.kdehr <- do.call(c, seas.kdehr)

for(i in 1:length(seas.akdehr)){
  # plot(seas.tel.df[names(seas.akdehr[i])],
  #      UD = seas.akdehr[[i]],
  #      main = seas.akdehr[[i]]@info$identity,
  #      error = F)
  seas.out[seas.out$IDYS == seas.akdehr[[i]]@info$identity,]$AKDE95 <- as.numeric(summary(seas.akdehr[[i]], units = F)$CI[[2]]) * 0.0001
  seas.out[seas.out$IDYS == seas.akdehr[[i]]@info$identity,]$AKDE50 <- as.numeric(summary(seas.akdehr[[i]], level.UD = 0.5, units = F)$CI[[2]]) * 0.0001
  seas.out[seas.out$IDYS == seas.akdehr[[i]]@info$identity,]$neff <- summary(seas.akdehr[[i]])$DOF[[1]]
}

for(i in 1:length(seas.kdehr)){
  # plot(seas.tel.df[names(seas.kdehr[i])],
  #      UD = seas.kdehr[[i]],
  #      main = seas.kdehr[[i]]@info$identity,
  #      error = F)
  seas.out[seas.out$IDYS == seas.kdehr[[i]]@info$identity,]$KDE95 <- as.numeric(summary(seas.kdehr[[i]], units = F)$CI[[2]]) * 0.0001
  seas.out[seas.out$IDYS == seas.kdehr[[i]]@info$identity,]$KDE50 <- as.numeric(summary(seas.kdehr[[i]], level.UD = 0.5, units = F)$CI[[2]]) * 0.0001
  seas.out[seas.out$IDYS == seas.kdehr[[i]]@info$identity,]$kdeneff <- summary(seas.kdehr[[i]])$DOF[[1]]
}

extra_info <- read.csv("./Extra Info.csv")
extra_info$ID <- as.factor(extra_info$ID)
extra_info$IDY <- interaction(extra_info$ID, 
                              extra_info$Year,
                              sep = "_",
                              drop = T)

extra_info %>%
  select(-ID) %>%
  left_join(x = seas.out, y = ., by = "IDY") -> seas.out

rm(extra_info, seas.first.fits, seas.iid, seas.iid.fits, seas.mods)
gc()

######################################################################

## STEP TWO: CALCULATE LANDSCAPE METRICS ##

# Create buffered points to crop NLCD and NWI rasters
# 5 km is well beyond known movement distances of spotted turtles
# So no relevant habitat is being removed

data %>%
  st_as_sf(coords = c("Longitude", "Latitude"),
           crs = 4326) %>%
  st_transform(crs = utm) %>%
  st_buffer(dist = set_units(5, "km")) -> buffer.polys

# Download the Annual NLCD raster for each year in the dataset

# First, create list of years, download locations, and download urls
# Only use first year for each individual represented in the data

data %>%
  group_by(IDY) %>%
  arrange(Date) %>%
  slice_head(n = 1) -> tmp.df

year <- unique(format(tmp.df$Date, "%Y"))
year <- as.list(unique(year))

temp_zip <- lapply(seq_along(year), function(x){
  temp <- tempfile(fileext = ".zip")})

zip_url <- lapply(seq_along(year), function(i){
  temp <- glue("https://www.mrlc.gov/downloads/sciweb1/shared/mrlc/data-bundles/Annual_NLCD_LndCov_{year[[i]]}_CU_C1V2.zip")})

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
        files = glue("Annual_NLCD_LndCov_{year[[i]]}_CU_C1V2.tif"),
        exdir = tempdir())
  })

# Read each annual raster into R

r <- lapply(seq_along(year), function(i){
  temp <- rast(file.path(tempdir(), glue("Annual_NLCD_LndCov_{year[[i]]}_CU_C1V2.tif")))
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
# 5 = bare land (including rock / clay / sand)
# 6 = developed land
# 7 = cultivated land

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
  x[x %in% c(31)] <- 5
  return(x)})
lc <- lapply(lc, function(x){
  x[x %in% c(21,22,23,24)] <- 6
  return(x)})
lc <- lapply(lc, function(x){
  x[x %in% c(81,82)] <- 7
  return(x)})

lc <- lapply(lc, function(x){
  coltab(x) <- data.frame(value = c(1,2,3,4,5,6,7),
                          color = c('deepskyblue',
                                    'deepskyblue4',
                                    'darkseagreen4',
                                    'darkseagreen',
                                    'khaki',
                                    'gray40',
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
                      disagg(lc[[1]], fact = 3),
                      field = "WETLAND_TYPE")

wet.rast <- terra::crop(wet.rast, mask.polys, mask = T)

# Remove lakes, rivers/streams, and marine wetlands and deepwater from NWI wetland raster

nwi.all <- mask(wet.rast,
                wet.rast %in% c("Lake", "Riverine", "Estuarine and Marine Wetland", "Estuarine and Marine Deepwater"),
                maskvalue = 1)

terra::plot(nwi.all)

nwi.all <- droplevels(nwi.all)
binary.nwi.all <- classify(nwi.all, cbind(levels(nwi.all)[[1]]$ID, 1))
binary.nwi.all <- classify(binary.nwi.all, cbind(NA, 0))

# Create raster of disjunct wetland patches and count cell size to get wetland areas

lc.water <- lapply(seq_along(lc), function(i){
  temp <- classify(lc[[i]], cbind(c(2,3,4,5,6,7), c(NA,NA,NA,NA,NA,NA)))
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
nwi.all.wetland.cells$count <- nwi.all.wetland.cells$count * ((10*10) / 1e4)

# Repeat for developed land

lc.dev <- lapply(seq_along(lc), function(i){
  temp <- classify(lc[[i]], cbind(c(1,2,3,4,5,7), c(NA,NA,NA,NA,NA,NA)))
})
dev.patches <- lapply(seq_along(lc.dev), function(i){
  temp <- patches(lc.dev[[i]], directions = 8)
})

names(lc.dev) <- names(dev.patches) <- year

# Repeat for cultivated land

lc.ag <- lapply(seq_along(lc), function(i){
  temp <- classify(lc[[i]], cbind(c(1,2,3,4,5,6), c(NA,NA,NA,NA,NA,NA)))
})
ag.patches <- lapply(seq_along(lc.ag), function(i){
  temp <- patches(lc.ag[[i]], directions = 8)
})

names(lc.ag) <- names(ag.patches) <- year

# Create buffers within which we'll calculate landscape variables

# Start by creating buffered activity centers

# First, full home ranges:

lapply(seq_along(akde.fits), function(i){
  st_as_sf(as.data.frame(akde.fits[[i]]$mu),
           coords = c("x","y"),
           crs = utm) %>%
    st_transform(crs = crs(lc[[1]])) %>%
    st_buffer(dist = set_units(200, "m")) %>%
    mutate(ID = akde.fits[[i]]@info$identity)
}) %>%
  bind_rows() -> low.center.buffer

lapply(seq_along(akde.fits), function(i){
  st_as_sf(as.data.frame(akde.fits[[i]]$mu),
           coords = c("x","y"),
           crs = utm) %>%
    st_transform(crs = crs(lc[[1]])) %>%
    st_buffer(dist = set_units(1000, "m")) %>%
    mutate(ID = akde.fits[[i]]@info$identity)
}) %>%
  bind_rows() -> hi.center.buffer

lapply(seq_along(kde.fits), function(i){
  st_as_sf(as.data.frame(kde.fits[[i]]$mu),
           coords = c("x","y"),
           crs = utm) %>%
    st_transform(crs = crs(lc[[1]])) %>%
    st_buffer(dist = set_units(200, "m")) %>%
    mutate(ID = kde.fits[[i]]@info$identity)
}) %>%
  bind_rows() -> low.center.buffer.kde

lapply(seq_along(kde.fits), function(i){
  st_as_sf(as.data.frame(kde.fits[[i]]$mu),
           coords = c("x","y"),
           crs = utm) %>%
    st_transform(crs = crs(lc[[1]])) %>%
    st_buffer(dist = set_units(1000, "m")) %>%
    mutate(ID = kde.fits[[i]]@info$identity)
}) %>%
  bind_rows() -> hi.center.buffer.kde

# Next, annual home ranges:

lapply(seq_along(ann.akde.fits), function(i){
  st_as_sf(as.data.frame(ann.akde.fits[[i]]$mu),
           coords = c("x","y"),
           crs = utm) %>%
    st_transform(crs = crs(lc[[1]])) %>%
    st_buffer(dist = set_units(200, "m")) %>%
    mutate(IDY = ann.akde.fits[[i]]@info$identity)
}) %>%
  bind_rows() -> ann.low.center.buffer

lapply(seq_along(ann.akde.fits), function(i){
  st_as_sf(as.data.frame(ann.akde.fits[[i]]$mu),
           coords = c("x","y"),
           crs = utm) %>%
    st_transform(crs = crs(lc[[1]])) %>%
    st_buffer(dist = set_units(1000, "m")) %>%
    mutate(IDY = ann.akde.fits[[i]]@info$identity)
}) %>%
  bind_rows() -> ann.hi.center.buffer

lapply(seq_along(ann.kde.fits), function(i){
  st_as_sf(as.data.frame(ann.kde.fits[[i]]$mu),
           coords = c("x","y"),
           crs = utm) %>%
    st_transform(crs = crs(lc[[1]])) %>%
    st_buffer(dist = set_units(200, "m")) %>%
    mutate(IDY = ann.kde.fits[[i]]@info$identity)
}) %>%
  bind_rows() -> ann.low.center.buffer.kde

lapply(seq_along(ann.kde.fits), function(i){
  st_as_sf(as.data.frame(ann.kde.fits[[i]]$mu),
           coords = c("x","y"),
           crs = utm) %>%
    st_transform(crs = crs(lc[[1]])) %>%
    st_buffer(dist = set_units(1000, "m")) %>%
    mutate(IDY = ann.kde.fits[[i]]@info$identity)
}) %>%
  bind_rows() -> ann.hi.center.buffer.kde

# Last, seasonal home ranges:

lapply(seq_along(seas.akde.fits), function(i){
  st_as_sf(as.data.frame(seas.akde.fits[[i]]$mu),
           coords = c("x","y"),
           crs = utm) %>%
    st_transform(crs = crs(lc[[1]])) %>%
    st_buffer(dist = set_units(200, "m")) %>%
    mutate(IDYS = seas.akde.fits[[i]]@info$identity)
}) %>%
  bind_rows() -> seas.low.center.buffer

lapply(seq_along(seas.akde.fits), function(i){
  st_as_sf(as.data.frame(seas.akde.fits[[i]]$mu),
           coords = c("x","y"),
           crs = utm) %>%
    st_transform(crs = crs(lc[[1]])) %>%
    st_buffer(dist = set_units(1000, "m")) %>%
    mutate(IDYS = seas.akde.fits[[i]]@info$identity)
}) %>%
  bind_rows() -> seas.hi.center.buffer

lapply(seq_along(seas.kde.fits), function(i){
  st_as_sf(as.data.frame(seas.kde.fits[[i]]$mu),
           coords = c("x","y"),
           crs = utm) %>%
    st_transform(crs = crs(lc[[1]])) %>%
    st_buffer(dist = set_units(200, "m")) %>%
    mutate(IDYS = seas.kde.fits[[i]]@info$identity)
}) %>%
  bind_rows() -> seas.low.center.buffer.kde

lapply(seq_along(seas.kde.fits), function(i){
  st_as_sf(as.data.frame(seas.kde.fits[[i]]$mu),
           coords = c("x","y"),
           crs = utm) %>%
    st_transform(crs = crs(lc[[1]])) %>%
    st_buffer(dist = set_units(1000, "m")) %>%
    mutate(IDYS = seas.kde.fits[[i]]@info$identity)
}) %>%
  bind_rows() -> seas.hi.center.buffer.kde

# Create object of initial points for each individual

data %>%
  arrange(Date) %>% 
  group_by(ID) %>%
  slice_head(n = 1) %>%
  st_as_sf(coords = c("Longitude", "Latitude"),
           crs = 4326) %>%
  st_transform(crs = crs(lc[[1]])) -> pts

data %>%
  arrange(Date) %>% 
  group_by(IDY) %>%
  slice_head(n = 1) %>%
  st_as_sf(coords = c("Longitude", "Latitude"),
           crs = 4326) %>%
  st_transform(crs = crs(lc[[1]])) -> ann.pts

data %>%
  arrange(Date) %>% 
  group_by(IDYS) %>%
  slice_head(n = 1) %>%
  st_as_sf(coords = c("Longitude", "Latitude"),
           crs = 4326) %>%
  st_transform(crs = crs(lc[[1]])) -> seas.pts

sf.df <- st_transform(sf.df, crs = crs(lc[[1]]))

# Create object with 95% AKDE polygons

akde.polys <- lapply(seq_along(akdehr), function(i){
  try(temp <- st_transform(as.sf(akdehr[[i]], level = 0.95), crs = crs(lc[[1]]))[2,], silent = T)
})

akde.polys <- Filter(is.list, akde.polys)

akde.names <- lapply(seq_along(akde.polys), function(i){
  temp <- sub(" .*$", "", akde.polys[[i]]$name)
})

names(akde.polys) <- akde.names

ann.akde.polys <- lapply(seq_along(ann.akdehr), function(i){
  try(temp <- st_transform(as.sf(ann.akdehr[[i]], level = 0.95), crs = crs(lc[[1]]))[2,], silent = T)
})

ann.akde.polys <- Filter(is.list, ann.akde.polys)

ann.akde.names <- lapply(seq_along(ann.akde.polys), function(i){
  temp <- sub(" .*$", "", ann.akde.polys[[i]]$name)
})

names(ann.akde.polys) <- ann.akde.names

seas.akde.polys <- lapply(seq_along(seas.akdehr), function(i){
  try(temp <- st_transform(as.sf(seas.akdehr[[i]], level = 0.95), crs = crs(lc[[1]]))[2,], silent = T)
})

seas.akde.polys <- Filter(is.list, seas.akde.polys)

seas.akde.names <- lapply(seq_along(seas.akde.polys), function(i){
  temp <- sub(" .*$", "", seas.akde.polys[[i]]$name)
})

names(seas.akde.polys) <- seas.akde.names

# Create object with 95% KDE polygons

kde.polys <- lapply(seq_along(kdehr), function(i){
  try(temp <- st_transform(as.sf(kdehr[[i]], level = 0.95), crs = crs(lc[[1]]))[2,], silent = T)
})

kde.polys <- Filter(is.list, kde.polys)

kde.names <- lapply(seq_along(kde.polys), function(i){
  temp <- sub(" .*$", "", kde.polys[[i]]$name)
})

names(kde.polys) <- kde.names

ann.kde.polys <- lapply(seq_along(ann.kdehr), function(i){
  try(temp <- st_transform(as.sf(ann.kdehr[[i]], level = 0.95), crs = crs(lc[[1]]))[2,], silent = T)
})

ann.kde.polys <- Filter(is.list, ann.kde.polys)

ann.kde.names <- lapply(seq_along(ann.kde.polys), function(i){
  temp <- sub(" .*$", "", ann.kde.polys[[i]]$name)
})

names(ann.kde.polys) <- ann.kde.names

seas.kde.polys <- lapply(seq_along(seas.kdehr), function(i){
  try(temp <- st_transform(as.sf(seas.kdehr[[i]], level = 0.95), crs = crs(lc[[1]]))[2,], silent = T)
})

seas.kde.polys <- Filter(is.list, seas.kde.polys)

seas.kde.names <- lapply(seq_along(seas.kde.polys), function(i){
  temp <- sub(" .*$", "", seas.kde.polys[[i]]$name)
})

names(seas.kde.polys) <- seas.kde.names

split_pts <- split(pts, pts$ID)
split_sf.df <- split(sf.df, sf.df$ID)
split_low.center.buffer <- split(low.center.buffer, low.center.buffer$ID)
split_hi.center.buffer <- split(hi.center.buffer, hi.center.buffer$ID)
split_low.center.buffer.kde <- split(low.center.buffer.kde, low.center.buffer.kde$ID)
split_hi.center.buffer.kde <- split(hi.center.buffer.kde, hi.center.buffer.kde$ID)

split_ann.pts <- split(ann.pts, ann.pts$IDY)
split_ann.sf.df <- split(sf.df, sf.df$IDY)
split_ann.low.center.buffer <- split(ann.low.center.buffer, ann.low.center.buffer$IDY)
split_ann.hi.center.buffer <- split(ann.hi.center.buffer, ann.hi.center.buffer$IDY)
split_ann.low.center.buffer.kde <- split(ann.low.center.buffer.kde, ann.low.center.buffer.kde$IDY)
split_ann.hi.center.buffer.kde <- split(ann.hi.center.buffer.kde, ann.hi.center.buffer.kde$IDY)

split_seas.pts <- split(seas.pts, seas.pts$IDYS)
split_seas.sf.df <- split(sf.df, sf.df$IDYS)
split_seas.low.center.buffer <- split(seas.low.center.buffer, seas.low.center.buffer$IDYS)
split_seas.hi.center.buffer <- split(seas.hi.center.buffer, seas.hi.center.buffer$IDYS)
split_seas.low.center.buffer.kde <- split(seas.low.center.buffer.kde, seas.low.center.buffer.kde$IDYS)
split_seas.hi.center.buffer.kde <- split(seas.hi.center.buffer.kde, seas.hi.center.buffer.kde$IDYS)

# Extract landscape variables of interest from annual NLCD rasters

for(i in as.character(unique(full.out$ID))) {
  yr <- as.character(full.out[full.out$ID == i,]$Year)
  # Calculate maximum and net displacement
  dists <- st_distance(split_pts[[i]], split_sf.df[[i]])
  full.out[full.out$ID==i, "dmax"] <- as.numeric(max(dists))
  full.out[full.out$ID==i, "dnet"] <- as.numeric(dists[length(dists)])
  # Calculate the proportion of wetland within each home range and buffered activity center
  if(i %in% names(akde.polys) == TRUE) {
    tmp <- terra::extract(wetland.patches[[yr]], split_low.center.buffer[[i]], exact = T)
    full.out[full.out$ID==i, "pwet.low.center"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(wetland.patches[[yr]], split_hi.center.buffer[[i]], exact = T)
    full.out[full.out$ID==i, "pwet.hi.center"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(wetland.patches[[yr]], akde.polys[[i]], exact = T)
    full.out[full.out$ID==i, "pwet.akde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
  } else {}
  if(i %in% names(kde.polys) == TRUE) {
    tmp <- terra::extract(wetland.patches[[yr]], split_low.center.buffer.kde[[i]], exact = T)
    full.out[full.out$ID==i, "pwet.low.center.kde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(wetland.patches[[yr]], split_hi.center.buffer.kde[[i]], exact = T)
    full.out[full.out$ID==i, "pwet.hi.center.kde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(wetland.patches[[yr]], kde.polys[[i]], exact = T)
    full.out[full.out$ID==i, "pwet.kde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
  } else {}
  # Calculate proportion of developed and agricultural land
  if(i %in% names(akde.polys) == TRUE) {
    tmp <- terra::extract(dev.patches[[yr]], split_low.center.buffer[[i]], exact = T)
    full.out[full.out$ID==i, "pdev.low.center"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(dev.patches[[yr]], split_hi.center.buffer[[i]], exact = T)
    full.out[full.out$ID==i, "pdev.hi.center"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(dev.patches[[yr]], akde.polys[[i]], exact = T)
    full.out[full.out$ID==i, "pdev.akde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(ag.patches[[yr]], split_low.center.buffer[[i]], exact = T)
    full.out[full.out$ID==i, "pag.low.center"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(ag.patches[[yr]], split_hi.center.buffer[[i]], exact = T)
    full.out[full.out$ID==i, "pag.hi.center"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(ag.patches[[yr]], akde.polys[[i]], exact = T)
    full.out[full.out$ID==i, "pag.akde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
  } else {}
  if(i %in% names(kde.polys) == TRUE) {
    tmp <- terra::extract(dev.patches[[yr]], split_low.center.buffer.kde[[i]], exact = T)
    full.out[full.out$ID==i, "pdev.low.center.kde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(dev.patches[[yr]], split_hi.center.buffer.kde[[i]], exact = T)
    full.out[full.out$ID==i, "pdev.hi.center.kde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(dev.patches[[yr]], kde.polys[[i]], exact = T)
    full.out[full.out$ID==i, "pdev.kde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(ag.patches[[yr]], split_low.center.buffer.kde[[i]], exact = T)
    full.out[full.out$ID==i, "pag.low.center.kde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(ag.patches[[yr]], split_hi.center.buffer.kde[[i]], exact = T)
    full.out[full.out$ID==i, "pag.hi.center.kde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(ag.patches[[yr]], kde.polys[[i]], exact = T)
    full.out[full.out$ID==i, "pag.kde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
  } else {}
  # Calculate clumpiness index for the wetland class
  if(i %in% names(akde.polys) == TRUE) {
    tmp <- terra::crop(lc[[yr]], split_low.center.buffer[[i]], touches = T, mask = T, snap = "full.out")
    tmp.clump <- lsm_c_clumpy(tmp)
    full.out[full.out$ID==i, "clumpy.low.center"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
    tmp <- terra::crop(lc[[yr]], split_hi.center.buffer[[i]], touches = T, mask = T, snap = "full.out")
    tmp.clump <- lsm_c_clumpy(tmp)
    full.out[full.out$ID==i, "clumpy.hi.center"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
    tmp <- terra::crop(lc[[yr]], akde.polys[[i]], touches = T, mask = T, snap = "full.out")
    tmp.clump <- lsm_c_clumpy(tmp)
    full.out[full.out$ID==i, "clumpy.akde"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
  } else {}
  if(i %in% names(kde.polys) == TRUE) {
    tmp <- terra::crop(lc[[yr]], split_low.center.buffer.kde[[i]], touches = T, mask = T, snap = "full.out")
    tmp.clump <- lsm_c_clumpy(tmp)
    full.out[full.out$ID==i, "clumpy.low.center.kde"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
    tmp <- terra::crop(lc[[yr]], split_hi.center.buffer.kde[[i]], touches = T, mask = T, snap = "full.out")
    tmp.clump <- lsm_c_clumpy(tmp)
    full.out[full.out$ID==i, "clumpy.hi.center.kde"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
    tmp <- terra::crop(lc[[yr]], kde.polys[[i]], touches = T, mask = T, snap = "full.out")
    tmp.clump <- lsm_c_clumpy(tmp)
    full.out[full.out$ID==i, "clumpy.kde"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
  } else {}
}

# Repeat for annual home ranges

for(i in as.character(unique(ann.out$IDY))){
  yr <- as.character(ann.out[ann.out$IDY == i,]$Year)
  dists <- st_distance(split_ann.pts[[i]], split_ann.sf.df[[i]])
  ann.out[ann.out$IDY==i, "dmax"] <- as.numeric(max(dists))
  ann.out[ann.out$IDY==i, "dnet"] <- as.numeric(dists[length(dists)])
  if(i %in% names(ann.akde.polys) == TRUE) {
    tmp <- terra::extract(wetland.patches[[yr]], split_ann.low.center.buffer[[i]], exact = T)
    ann.out[ann.out$IDY==i, "pwet.ann.low.center"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(wetland.patches[[yr]], split_ann.hi.center.buffer[[i]], exact = T)
    ann.out[ann.out$IDY==i, "pwet.ann.hi.center"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(wetland.patches[[yr]], ann.akde.polys[[i]], exact = T)
    ann.out[ann.out$IDY==i, "pwet.akde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
  } else {}
  if(i %in% names(ann.kde.polys) == TRUE) {
    tmp <- terra::extract(wetland.patches[[yr]], split_ann.low.center.buffer.kde[[i]], exact = T)
    ann.out[ann.out$IDY==i, "pwet.ann.low.center.kde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(wetland.patches[[yr]], split_ann.hi.center.buffer.kde[[i]], exact = T)
    ann.out[ann.out$IDY==i, "pwet.ann.hi.center.kde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(wetland.patches[[yr]], ann.kde.polys[[i]], exact = T)
    ann.out[ann.out$IDY==i, "pwet.kde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
  } else {}
  if(i %in% names(ann.akde.polys) == TRUE) {
    tmp <- terra::extract(dev.patches[[yr]], split_ann.low.center.buffer[[i]], exact = T)
    ann.out[ann.out$IDY==i, "pdev.ann.low.center"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(dev.patches[[yr]], split_ann.hi.center.buffer[[i]], exact = T)
    ann.out[ann.out$IDY==i, "pdev.ann.hi.center"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(dev.patches[[yr]], ann.akde.polys[[i]], exact = T)
    ann.out[ann.out$IDY==i, "pdev.akde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(ag.patches[[yr]], split_ann.low.center.buffer[[i]], exact = T)
    ann.out[ann.out$IDY==i, "pag.ann.low.center"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(ag.patches[[yr]], split_ann.hi.center.buffer[[i]], exact = T)
    ann.out[ann.out$IDY==i, "pag.ann.hi.center"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(ag.patches[[yr]], ann.akde.polys[[i]], exact = T)
    ann.out[ann.out$IDY==i, "pag.akde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
  } else {}
  if(i %in% names(ann.kde.polys) == TRUE) {
    tmp <- terra::extract(dev.patches[[yr]], split_ann.low.center.buffer.kde[[i]], exact = T)
    ann.out[ann.out$IDY==i, "pdev.ann.low.center.kde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(dev.patches[[yr]], split_ann.hi.center.buffer.kde[[i]], exact = T)
    ann.out[ann.out$IDY==i, "pdev.ann.hi.center.kde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(dev.patches[[yr]], ann.kde.polys[[i]], exact = T)
    ann.out[ann.out$IDY==i, "pdev.kde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(ag.patches[[yr]], split_ann.low.center.buffer.kde[[i]], exact = T)
    ann.out[ann.out$IDY==i, "pag.ann.low.center.kde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(ag.patches[[yr]], split_ann.hi.center.buffer.kde[[i]], exact = T)
    ann.out[ann.out$IDY==i, "pag.ann.hi.center.kde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(ag.patches[[yr]], ann.kde.polys[[i]], exact = T)
    ann.out[ann.out$IDY==i, "pag.kde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
  } else {}
  if(i %in% names(ann.akde.polys) == TRUE) {
    tmp <- terra::crop(lc[[yr]], split_ann.low.center.buffer[[i]], touches = T, mask = T, snap = "ann.out")
    tmp.clump <- lsm_c_clumpy(tmp)
    ann.out[ann.out$IDY==i, "clumpy.ann.low.center"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
    tmp <- terra::crop(lc[[yr]], split_ann.hi.center.buffer[[i]], touches = T, mask = T, snap = "ann.out")
    tmp.clump <- lsm_c_clumpy(tmp)
    ann.out[ann.out$IDY==i, "clumpy.ann.hi.center"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
    tmp <- terra::crop(lc[[yr]], ann.akde.polys[[i]], touches = T, mask = T, snap = "ann.out")
    tmp.clump <- lsm_c_clumpy(tmp)
    ann.out[ann.out$IDY==i, "clumpy.akde"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
  } else {}
  if(i %in% names(ann.kde.polys) == TRUE) {
    tmp <- terra::crop(lc[[yr]], split_ann.low.center.buffer.kde[[i]], touches = T, mask = T, snap = "ann.out")
    tmp.clump <- lsm_c_clumpy(tmp)
    ann.out[ann.out$IDY==i, "clumpy.ann.low.center.kde"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
    tmp <- terra::crop(lc[[yr]], split_ann.hi.center.buffer.kde[[i]], touches = T, mask = T, snap = "ann.out")
    tmp.clump <- lsm_c_clumpy(tmp)
    ann.out[ann.out$IDY==i, "clumpy.ann.hi.center.kde"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
    tmp <- terra::crop(lc[[yr]], ann.kde.polys[[i]], touches = T, mask = T, snap = "ann.out")
    tmp.clump <- lsm_c_clumpy(tmp)
    ann.out[ann.out$IDY==i, "clumpy.kde"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
  } else {}
}

# Repeat for seasonal home ranges

for(i in as.character(unique(seas.out$IDYS))){
  yr <- as.character(seas.out[seas.out$IDYS == i,]$Year)
  dists <- st_distance(split_seas.pts[[i]], split_seas.sf.df[[i]])
  seas.out[seas.out$IDYS==i, "dmax"] <- as.numeric(max(dists))
  seas.out[seas.out$IDYS==i, "dnet"] <- as.numeric(dists[length(dists)])
  if(i %in% names(seas.akde.polys) == TRUE) {
    tmp <- terra::extract(wetland.patches[[yr]], split_seas.low.center.buffer[[i]], exact = T)
    seas.out[seas.out$IDYS==i, "pwet.seas.low.center"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(wetland.patches[[yr]], split_seas.hi.center.buffer[[i]], exact = T)
    seas.out[seas.out$IDYS==i, "pwet.seas.hi.center"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(wetland.patches[[yr]], seas.akde.polys[[i]], exact = T)
    seas.out[seas.out$IDYS==i, "pwet.akde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
  } else {}
  if(i %in% names(seas.kde.polys) == TRUE) {
    tmp <- terra::extract(wetland.patches[[yr]], split_seas.low.center.buffer.kde[[i]], exact = T)
    seas.out[seas.out$IDYS==i, "pwet.seas.low.center.kde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(wetland.patches[[yr]], split_seas.hi.center.buffer.kde[[i]], exact = T)
    seas.out[seas.out$IDYS==i, "pwet.seas.hi.center.kde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(wetland.patches[[yr]], seas.kde.polys[[i]], exact = T)
    seas.out[seas.out$IDYS==i, "pwet.kde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
  } else {}
  if(i %in% names(seas.akde.polys) == TRUE) {
    tmp <- terra::extract(dev.patches[[yr]], split_seas.low.center.buffer[[i]], exact = T)
    seas.out[seas.out$IDYS==i, "pdev.seas.low.center"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(dev.patches[[yr]], split_seas.hi.center.buffer[[i]], exact = T)
    seas.out[seas.out$IDYS==i, "pdev.seas.hi.center"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(dev.patches[[yr]], seas.akde.polys[[i]], exact = T)
    seas.out[seas.out$IDYS==i, "pdev.akde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(ag.patches[[yr]], split_seas.low.center.buffer[[i]], exact = T)
    seas.out[seas.out$IDYS==i, "pag.seas.low.center"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(ag.patches[[yr]], split_seas.hi.center.buffer[[i]], exact = T)
    seas.out[seas.out$IDYS==i, "pag.seas.hi.center"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(ag.patches[[yr]], seas.akde.polys[[i]], exact = T)
    seas.out[seas.out$IDYS==i, "pag.akde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
  } else {}
  if(i %in% names(seas.kde.polys) == TRUE) {
    tmp <- terra::extract(dev.patches[[yr]], split_seas.low.center.buffer.kde[[i]], exact = T)
    seas.out[seas.out$IDYS==i, "pdev.seas.low.center.kde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(dev.patches[[yr]], split_seas.hi.center.buffer.kde[[i]], exact = T)
    seas.out[seas.out$IDYS==i, "pdev.seas.hi.center.kde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(dev.patches[[yr]], seas.kde.polys[[i]], exact = T)
    seas.out[seas.out$IDYS==i, "pdev.kde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(ag.patches[[yr]], split_seas.low.center.buffer.kde[[i]], exact = T)
    seas.out[seas.out$IDYS==i, "pag.seas.low.center.kde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(ag.patches[[yr]], split_seas.hi.center.buffer.kde[[i]], exact = T)
    seas.out[seas.out$IDYS==i, "pag.seas.hi.center.kde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(ag.patches[[yr]], seas.kde.polys[[i]], exact = T)
    seas.out[seas.out$IDYS==i, "pag.kde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
  } else {}
  if(i %in% names(seas.akde.polys) == TRUE) {
    tmp <- terra::crop(lc[[yr]], split_seas.low.center.buffer[[i]], touches = T, mask = T, snap = "seas.out")
    tmp.clump <- lsm_c_clumpy(tmp)
    seas.out[seas.out$IDYS==i, "clumpy.seas.low.center"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
    tmp <- terra::crop(lc[[yr]], split_seas.hi.center.buffer[[i]], touches = T, mask = T, snap = "seas.out")
    tmp.clump <- lsm_c_clumpy(tmp)
    seas.out[seas.out$IDYS==i, "clumpy.seas.hi.center"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
    tmp <- terra::crop(lc[[yr]], seas.akde.polys[[i]], touches = T, mask = T, snap = "seas.out")
    tmp.clump <- lsm_c_clumpy(tmp)
    seas.out[seas.out$IDYS==i, "clumpy.akde"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
  } else {}
  if(i %in% names(seas.kde.polys) == TRUE) {
    tmp <- terra::crop(lc[[yr]], split_seas.low.center.buffer.kde[[i]], touches = T, mask = T, snap = "seas.out")
    tmp.clump <- lsm_c_clumpy(tmp)
    seas.out[seas.out$IDYS==i, "clumpy.seas.low.center.kde"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
    tmp <- terra::crop(lc[[yr]], split_seas.hi.center.buffer.kde[[i]], touches = T, mask = T, snap = "seas.out")
    tmp.clump <- lsm_c_clumpy(tmp)
    seas.out[seas.out$IDYS==i, "clumpy.seas.hi.center.kde"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
    tmp <- terra::crop(lc[[yr]], seas.kde.polys[[i]], touches = T, mask = T, snap = "seas.out")
    tmp.clump <- lsm_c_clumpy(tmp)
    seas.out[seas.out$IDYS==i, "clumpy.kde"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
  } else {}
}

# Calculate landscape variables using NWI rasters

for(i in as.character(unique(full.out$ID))) {
  # Calculate the proportion of wetland within each home range and buffered activity center
  if(i %in% names(akde.polys) == TRUE) {
    tmp <- terra::extract(nwi.all.wetland.patches, split_low.center.buffer[[i]], exact = T)
    full.out[full.out$ID==i, "NWI.pwet.low.center"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(nwi.all.wetland.patches, split_hi.center.buffer[[i]], exact = T)
    full.out[full.out$ID==i, "NWI.pwet.hi.center"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(nwi.all.wetland.patches, akde.polys[[i]], exact = T)
    full.out[full.out$ID==i, "NWI.pwet.akde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
  } else {}
  if(i %in% names(kde.polys) == TRUE) {
    tmp <- terra::extract(nwi.all.wetland.patches, split_low.center.buffer.kde[[i]], exact = T)
    full.out[full.out$ID==i, "NWI.pwet.low.center.kde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(nwi.all.wetland.patches, split_hi.center.buffer.kde[[i]], exact = T)
    full.out[full.out$ID==i, "NWI.pwet.hi.center.kde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(nwi.all.wetland.patches, kde.polys[[i]], exact = T)
    full.out[full.out$ID==i, "NWI.pwet.kde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
  } else {}
  # Calculate clumpiness index for the wetland class
  if(i %in% names(akde.polys) == TRUE) {
    tmp <- terra::crop(binary.nwi.all, split_low.center.buffer[[i]], touches = T, mask = T, snap = "full.out")
    tmp.clump <- lsm_c_clumpy(tmp)
    full.out[full.out$ID==i, "NWI.clumpy.low.center"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
    tmp <- terra::crop(binary.nwi.all, split_hi.center.buffer[[i]], touches = T, mask = T, snap = "full.out")
    tmp.clump <- lsm_c_clumpy(tmp)
    full.out[full.out$ID==i, "NWI.clumpy.hi.center"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
    tmp <- terra::crop(binary.nwi.all, akde.polys[[i]], touches = T, mask = T, snap = "full.out")
    tmp.clump <- lsm_c_clumpy(tmp)
    full.out[full.out$ID==i, "NWI.clumpy.akde"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
  } else {}
  if(i %in% names(kde.polys) == TRUE) {
    tmp <- terra::crop(binary.nwi.all, split_low.center.buffer.kde[[i]], touches = T, mask = T, snap = "full.out")
    tmp.clump <- lsm_c_clumpy(tmp)
    full.out[full.out$ID==i, "NWI.clumpy.low.center.kde"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
    tmp <- terra::crop(binary.nwi.all, split_hi.center.buffer.kde[[i]], touches = T, mask = T, snap = "full.out")
    tmp.clump <- lsm_c_clumpy(tmp)
    full.out[full.out$ID==i, "NWI.clumpy.hi.center.kde"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
    tmp <- terra::crop(binary.nwi.all, kde.polys[[i]], touches = T, mask = T, snap = "full.out")
    tmp.clump <- lsm_c_clumpy(tmp)
    full.out[full.out$ID==i, "NWI.clumpy.kde"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
  } else {}
}

# Repeat for annual home ranges

for(i in as.character(unique(ann.out$IDY))) {
  # Calculate the proportion of wetland within each home range and buffered activity center
  if(i %in% names(ann.akde.polys) == TRUE) {
    tmp <- terra::extract(nwi.all.wetland.patches, split_ann.low.center.buffer[[i]], exact = T)
    ann.out[ann.out$IDY==i, "NWI.pwet.low.center"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(nwi.all.wetland.patches, split_ann.hi.center.buffer[[i]], exact = T)
    ann.out[ann.out$IDY==i, "NWI.pwet.hi.center"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(nwi.all.wetland.patches, ann.akde.polys[[i]], exact = T)
    ann.out[ann.out$IDY==i, "NWI.pwet.akde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
  } else {}
  if(i %in% names(ann.kde.polys) == TRUE) {
    tmp <- terra::extract(nwi.all.wetland.patches, split_ann.low.center.buffer.kde[[i]], exact = T)
    ann.out[ann.out$IDY==i, "NWI.pwet.low.center.kde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(nwi.all.wetland.patches, split_ann.hi.center.buffer.kde[[i]], exact = T)
    ann.out[ann.out$IDY==i, "NWI.pwet.hi.center.kde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(nwi.all.wetland.patches, ann.kde.polys[[i]], exact = T)
    ann.out[ann.out$IDY==i, "NWI.pwet.kde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
  } else {}
  # Calculate clumpiness index for the wetland class
  if(i %in% names(ann.akde.polys) == TRUE) {
    tmp <- terra::crop(binary.nwi.all, split_ann.low.center.buffer[[i]], touches = T, mask = T, snap = "ann.out")
    tmp.clump <- lsm_c_clumpy(tmp)
    ann.out[ann.out$IDY==i, "NWI.clumpy.low.center"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
    tmp <- terra::crop(binary.nwi.all, split_ann.hi.center.buffer[[i]], touches = T, mask = T, snap = "ann.out")
    tmp.clump <- lsm_c_clumpy(tmp)
    ann.out[ann.out$IDY==i, "NWI.clumpy.hi.center"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
    tmp <- terra::crop(binary.nwi.all, ann.akde.polys[[i]], touches = T, mask = T, snap = "ann.out")
    tmp.clump <- lsm_c_clumpy(tmp)
    ann.out[ann.out$IDY==i, "NWI.clumpy.akde"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
  } else {}
  if(i %in% names(ann.kde.polys) == TRUE) {
    tmp <- terra::crop(binary.nwi.all, split_ann.low.center.buffer.kde[[i]], touches = T, mask = T, snap = "ann.out")
    tmp.clump <- lsm_c_clumpy(tmp)
    ann.out[ann.out$IDY==i, "NWI.clumpy.low.center.kde"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
    tmp <- terra::crop(binary.nwi.all, split_ann.hi.center.buffer.kde[[i]], touches = T, mask = T, snap = "ann.out")
    tmp.clump <- lsm_c_clumpy(tmp)
    ann.out[ann.out$IDY==i, "NWI.clumpy.hi.center.kde"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
    tmp <- terra::crop(binary.nwi.all, ann.kde.polys[[i]], touches = T, mask = T, snap = "ann.out")
    tmp.clump <- lsm_c_clumpy(tmp)
    ann.out[ann.out$IDY==i, "NWI.clumpy.kde"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
  } else {}
}

# Repeat for seasonal home ranges

for(i in as.character(unique(seas.out$IDYS))) {
  # Calculate the proportion of wetland within each home range and buffered activity center
  if(i %in% names(seas.akde.polys) == TRUE) {
    tmp <- terra::extract(nwi.all.wetland.patches, split_seas.low.center.buffer[[i]], exact = T)
    seas.out[seas.out$IDYS==i, "NWI.pwet.low.center"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(nwi.all.wetland.patches, split_seas.hi.center.buffer[[i]], exact = T)
    seas.out[seas.out$IDYS==i, "NWI.pwet.hi.center"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(nwi.all.wetland.patches, seas.akde.polys[[i]], exact = T)
    seas.out[seas.out$IDYS==i, "NWI.pwet.akde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
  } else {}
  if(i %in% names(seas.kde.polys) == TRUE) {
    tmp <- terra::extract(nwi.all.wetland.patches, split_seas.low.center.buffer.kde[[i]], exact = T)
    seas.out[seas.out$IDYS==i, "NWI.pwet.low.center.kde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(nwi.all.wetland.patches, split_seas.hi.center.buffer.kde[[i]], exact = T)
    seas.out[seas.out$IDYS==i, "NWI.pwet.hi.center.kde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
    tmp <- terra::extract(nwi.all.wetland.patches, seas.kde.polys[[i]], exact = T)
    seas.out[seas.out$IDYS==i, "NWI.pwet.kde"] <- sum(na.omit(tmp)[,3]) / sum(tmp[,3])
  } else {}
  # Calculate clumpiness index for the wetland class
  if(i %in% names(seas.akde.polys) == TRUE) {
    tmp <- terra::crop(binary.nwi.all, split_seas.low.center.buffer[[i]], touches = T, mask = T, snap = "seas.out")
    tmp.clump <- lsm_c_clumpy(tmp)
    seas.out[seas.out$IDYS==i, "NWI.clumpy.low.center"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
    tmp <- terra::crop(binary.nwi.all, split_seas.hi.center.buffer[[i]], touches = T, mask = T, snap = "seas.out")
    tmp.clump <- lsm_c_clumpy(tmp)
    seas.out[seas.out$IDYS==i, "NWI.clumpy.hi.center"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
    tmp <- terra::crop(binary.nwi.all, seas.akde.polys[[i]], touches = T, mask = T, snap = "seas.out")
    tmp.clump <- lsm_c_clumpy(tmp)
    seas.out[seas.out$IDYS==i, "NWI.clumpy.akde"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
  } else {}
  if(i %in% names(seas.kde.polys) == TRUE) {
    tmp <- terra::crop(binary.nwi.all, split_seas.low.center.buffer.kde[[i]], touches = T, mask = T, snap = "seas.out")
    tmp.clump <- lsm_c_clumpy(tmp)
    seas.out[seas.out$IDYS==i, "NWI.clumpy.low.center.kde"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
    tmp <- terra::crop(binary.nwi.all, split_seas.hi.center.buffer.kde[[i]], touches = T, mask = T, snap = "seas.out")
    tmp.clump <- lsm_c_clumpy(tmp)
    seas.out[seas.out$IDYS==i, "NWI.clumpy.hi.center.kde"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
    tmp <- terra::crop(binary.nwi.all, seas.kde.polys[[i]], touches = T, mask = T, snap = "seas.out")
    tmp.clump <- lsm_c_clumpy(tmp)
    seas.out[seas.out$IDYS==i, "NWI.clumpy.kde"] <- as.numeric(tmp.clump[tmp.clump$class == 1, 6])
  } else {}
}

# For individuals with data across multiple years, see how many out-of-sample points fall within AKDE

for(i in unique(ann.out$ID)){
  if(nrow(ann.out[ann.out$ID == i,]) < 2){
    ann.out[ann.out$ID == i, "AKDEcoverage.in"] <- NA
    ann.out[ann.out$ID == i, "AKDEcoverage.out"] <- NA
  }else{
    for(j in unique(ann.out[ann.out$ID == i,]$IDY)){
      if(j %in% names(ann.akde.polys) == FALSE){
        ann.out[ann.out$IDY == j, "AKDEcoverage.in"] <- NA
        ann.out[ann.out$IDY == j, "AKDEcoverage.out"] <- NA
      }else{
        tmp.df <- sf.df[sf.df$ID == i & sf.df$IDY != j,]
        tmp <- st_intersects(tmp.df, ann.akde.polys[[j]], sparse = F)
        ann.out[ann.out$IDY == j, "AKDEcoverage.in"] <- sum(tmp)
        ann.out[ann.out$IDY == j, "AKDEcoverage.out"] <- nrow(tmp.df) - sum(tmp)
      }
    }
  }
}

# Repeat with KDE estimates

kde.polys <- lapply(seq_along(ann.kdehr), function(i){
  try(temp <- st_transform(as.sf(ann.kdehr[[i]], level = 0.95), crs = crs(lc[[1]]))[2,], silent = T)
})
kde.polys <- Filter(is.list, kde.polys)
kde.names <- lapply(seq_along(kde.polys), function(i){
  temp <- sub(" .*$", "", kde.polys[[i]]$name)
})
names(kde.polys) <- kde.names

for(i in unique(ann.out$ID)){
  if(nrow(ann.out[ann.out$ID == i,]) < 2){
    ann.out[ann.out$ID == i, "KDEcoverage.in"] <- NA
    ann.out[ann.out$ID == i, "KDEcoverage.out"] <- NA
  }else{
    for(j in unique(ann.out[ann.out$ID == i,]$IDY)){
      if(j %in% names(kde.polys) == FALSE){
        ann.out[ann.out$IDY == j, "KDEcoverage.in"] <- NA
        ann.out[ann.out$IDY == j, "KDEcoverage.out"] <- NA
      }else{
        tmp.df <- sf.df[sf.df$ID == i & sf.df$IDY != j,]
        tmp <- st_intersects(tmp.df, kde.polys[[j]], sparse = F)
        ann.out[ann.out$IDY == j, "KDEcoverage.in"] <- sum(tmp)
        ann.out[ann.out$IDY == j, "KDEcoverage.out"] <- nrow(tmp.df) - sum(tmp)
      }
    }
  }
}

# Repeat with 95% MCPs

sf.df[!sf.df$IDY %in% ann.out[ann.out$Points < 5,]$IDY,] -> ade.df
ade.df$IDY <- droplevels(ade.df$IDY)
ade.df <- as(ade.df[,"IDY"], "Spatial")

mcps <- st_transform(st_as_sf(mcp(ade.df, percent = 95)),
                     crs = crs(lc[[1]]))

for(i in unique(ann.out$ID)){
  if(nrow(ann.out[ann.out$ID == i,]) < 2){
    ann.out[ann.out$ID == i, "MCPcoverage.in"] <- NA
    ann.out[ann.out$ID == i, "MCPcoverage.out"] <- NA
  }else{
    for(j in unique(ann.out[ann.out$ID == i,]$IDY)){
      if(j %in% mcps$id == FALSE){
        ann.out[ann.out$IDY == j, "MCPcoverage.in"] <- NA
        ann.out[ann.out$IDY == j, "MCPcoverage.out"] <- NA
      }else{
        tmp.df <- sf.df[sf.df$ID == i & sf.df$IDY != j,]
        tmp <- st_intersects(tmp.df, mcps[mcps$id == j,], sparse = F)
        ann.out[ann.out$IDY == j, "MCPcoverage.in"] <- sum(tmp)
        ann.out[ann.out$IDY == j, "MCPcoverage.out"] <- nrow(tmp.df) - sum(tmp)
      }
    }
  }
}

######################################################################

## STEP THREE: CALCULATE PRECIPITATION AND TEMPERATURE METRICS ##

extra_info <- read.csv("./Extra Info.csv")
extra_info$ID <- as.factor(extra_info$ID)
extra_info$IDY <- interaction(extra_info$ID, 
                              extra_info$Year,
                              sep = "_",
                              drop = T)

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
                                 start = 1980,
                                 end = 2025)

lapply(clim.df,
       function(x) {
         mutate(x$data,
                Site = x$site,
                Date = as.Date(paste(year, yday, sep = "-"), "%Y-%j")) %>%
           mutate(Day = format(Date, "%d"),
                  Month = format(Date, "%m"),
                  Year = as.numeric(format(Date, "%Y")),
                  Precip = prcp..mm.day.,
                  Temp = (tmax..deg.c.+tmin..deg.c.)/2) %>%
           select(Site, Date, Year, Month, Day, Precip, Temp)
       }) %>%
  bind_rows() -> climate.df

# Add long-term precipitation mean and variance, days with precipitation

climate.df %>%
  group_by(Site, Year) %>%
  summarize(Ann.Precip = sum(Precip, na.rm = T),
            PrecipDays = sum(Precip > 0, na.rm = T)) %>%
  group_by(Site) %>%
  summarize(Av.Precip = mean(Ann.Precip, na.rm = T),
            Var.Precip = sd(Ann.Precip, na.rm = T),
            Av.PrecipDays = mean(PrecipDays, na.rm = T)) -> precip.av

full.out %>%
  left_join(y = precip.av, by = "Site") -> full.out

# Add long-term temperature mean and variance, seasonality

climate.df %>%
  group_by(Site, Year) %>%
  summarize(Ann.Temp = mean(Temp, na.rm = T),
            Seasonality = max(Temp) - min(Temp)) %>%
  group_by(Site) %>%
  summarize(Av.Temp = mean(Ann.Temp, na.rm = T),
            Var.Temp = sd(Ann.Temp, na.rm = T),
            Av.Seas = mean(Seasonality, na.rm = T),
            Var.Seas = sd(Seasonality, na.rm = T)) -> temp.av

full.out %>%
  left_join(y = temp.av, by = "Site") -> full.out

# Add new ID column that combines ID and Site name to ensure the IDs are unique

full.out %>%
  mutate(UniqueID = interaction(as.character(Site), 
                                ID,
                                sep = "_",
                                drop = T)) %>%
  select(UniqueID, ID, Site, Year, Sex, SCL, Mass, everything()) -> full.out

# Save variogram plots to pdf and out object to "output.csv"

write.csv(full.out, 
          "full_output.csv",
          row.names = F)

vario.names <- lapply(seq_along(variograms), function(i){
  temp <- sub(" .*$", "", variograms[[i]]@info$identity)
})

names(variograms) <- vario.names

akde.names <- lapply(seq_along(akde.fits), function(i){
  temp <- sub(" .*$", "", akde.fits[[i]]@info$identity)
})

names(akde.fits) <- akde.names

tmp.var <- variograms[names(variograms) %in% names(akde.fits)]

pdf("full_variograms.pdf")

par(mfrow = c(3,2))
lapply(seq_along(tmp.var), function(i){
  plot(tmp.var[[i]],
       CTMM = akde.fits[[i]],
       level = c(0.5, 0.95),
       fraction = 1,
       main = tmp.var[[i]]@info$identity)
})

dev.off()

# Add climate variables to annual home range estimates
# Start with water year (WY) precipitation and number of days with precipitation

climate.df$WY <- ifelse(climate.df$Month %in% c("10","11","12"),
                        as.numeric(climate.df$Year) + 1,
                        as.numeric(climate.df$Year))

wy.df <- climate.df[climate.df$WY <= 2025 & climate.df$WY > 1980,]

wy.df %>%
  group_by(Site, WY) %>%
  summarize(WY.Precip = sum(Precip, na.rm = T)) %>%
  group_by(Site) %>%
  summarize(Av.WY.Precip = mean(WY.Precip, na.rm = T),
            Var.WY.Precip = sd(WY.Precip, na.rm = T)) -> WY.av

wy.df %>%
  group_by(Site, WY) %>%
  summarize(Year = min(WY),
            WY.Precip = sum(Precip, na.rm = T)) %>%
  left_join(y = WY.av, by = "Site") %>%
  select(Site, Year, WY.Precip, Av.WY.Precip, Var.WY.Precip) -> WY.precip

ann.out %>%
  left_join(y = WY.precip, by = c("Site", "Year")) -> ann.out

climate.df %>%
  group_by(Site, Year) %>%
  summarize(PrecipDays = sum(Precip > 0, na.rm = T)) -> precip.days

ann.out %>%
  left_join(y = precip.days, by = c("Site", "Year")) -> ann.out

# Add long-term temperature mean, annual temperature mean

climate.df %>%
  group_by(Site, Year) %>%
  summarize(Ann.Temp = mean(Temp, na.rm = T),
            Seasonality = max(Temp) - min(Temp)) %>%
  group_by(Site) %>%
  summarize(Av.Temp = mean(Ann.Temp, na.rm = T),
            Var.Temp = sd(Ann.Temp, na.rm = T),
            Av.Seas = mean(Seasonality, na.rm = T),
            Var.Seas = sd(Seasonality, na.rm = T)) -> temp.av

climate.df %>%
  group_by(Site, Year) %>%
  summarize(Year = min(Year),
            Ann.Temp = mean(Temp, na.rm = T),
            Ann.Seasonality = max(Temp) - min(Temp)) %>%
  left_join(y = temp.av, by = "Site") -> temp

ann.out %>%
  left_join(y = temp, by = c("Site", "Year")) -> ann.out

# Add new ID column that combines ID and Site name to ensure the IDs are unique

ann.out %>%
  mutate(UniqueID = interaction(as.character(Site), 
                                ID,
                                sep = "_",
                                drop = T)) %>%
  select(UniqueID, ID, Site, Year, IDY, Sex, SCL, Mass, everything()) -> ann.out

# Save variogram plots to pdf and out object to "output.csv"

write.csv(ann.out, 
          "annual_output.csv",
          row.names = F)

vario.names <- lapply(seq_along(ann.variograms), function(i){
  temp <- sub(" .*$", "", ann.variograms[[i]]@info$identity)
})

names(ann.variograms) <- vario.names

ann.akde.names <- lapply(seq_along(ann.akde.fits), function(i){
  temp <- sub(" .*$", "", ann.akde.fits[[i]]@info$identity)
})

names(ann.akde.fits) <- ann.akde.names

tmp.var <- ann.variograms[names(ann.variograms) %in% names(ann.akde.fits)]

pdf("annual_variograms.pdf")

par(mfrow = c(3,2))
lapply(seq_along(tmp.var), function(i){
  plot(tmp.var[[i]],
       CTMM = ann.akde.fits[[i]],
       level = c(0.5, 0.95),
       fraction = 1,
       main = tmp.var[[i]]@info$identity)
})

dev.off()

# Add climate variables to seasonal home range estimates
# Start with precipitation and number of days with precipitation

climate.df$WY <- ifelse(climate.df$Month == "12",
                          as.numeric(climate.df$Year) + 1,
                          as.numeric(climate.df$Year))

climate.df$Season <- ifelse(climate.df$Month %in% c("12", "01", "02"),
                            "Winter",
                            ifelse(climate.df$Month %in% c("03", "04") | climate.df$Month == "05" & climate.df$Day < 16,
                                   "Pre-nesting",
                                   ifelse(climate.df$Month == "05" & climate.df$Day > 15 | climate.df$Month == "06" | climate.df$Month == "07" & climate.df$Day < 16,
                                          "Nesting",
                                          ifelse(climate.df$Month == "07" & climate.df$Day > 15 | climate.df$Month == "08",
                                                 "Summer",
                                                 "Fall"))))

seasonal.df <- climate.df[!(climate.df$WY == 1980 & climate.df$Season == "Winter"),]
seasonal.df <- seasonal.df[!(seasonal.df$WY == 2026 & seasonal.df$Season == "Winter"),]

seasonal.df %>%
  group_by(Site, WY, Season) %>%
  summarize(Seas.Precip = sum(Precip, na.rm = T)) %>%
  group_by(Site, Season) %>%
  summarize(Av.Precip = mean(Seas.Precip, na.rm = T),
            Var.Precip = sd(Seas.Precip, na.rm = T)) -> precip.av

seasonal.df %>%
  group_by(Site, WY, Season) %>%
  summarize(Year = max(WY),
            Precip = sum(Precip, na.rm = T)) %>%
  left_join(y = precip.av, by = c("Site", "Season")) -> seas.precip

seas.out %>%
  left_join(seas.precip, by = c("Site", "Year", "Season")) -> seas.out

seasonal.df %>%
  group_by(Site, WY, Season) %>%
  summarize(Year = max(WY),
            PrecipDays = sum(Precip > 0, na.rm = T)) %>%
  ungroup() %>%
  select(Site, Year, Season, PrecipDays) -> precip.days

seas.out %>%
  left_join(y = precip.days, by = c("Site", "Year", "Season")) -> seas.out

# Add long-term temperature mean and variance, and mean seasonal temp

seasonal.df %>%
  group_by(Site, Year, Season) %>%
  summarize(Seas.Temp = mean(Temp, na.rm = T)) %>%
  group_by(Site, Season) %>%
  summarize(Av.Temp = mean(Seas.Temp, na.rm = T),
            Var.Temp = sd(Seas.Temp, na.rm = T)) -> temp.av

seasonal.df %>%
  group_by(Site, Year, Season) %>%
  summarize(Year = min(Year),
            Temp = mean(Temp, na.rm = T)) %>%
  left_join(y = temp.av, by = c("Site", "Season")) %>%
  select(Site, Year, Season, Temp, Av.Temp, Var.Temp) -> temp

seas.out %>%
  left_join(y = temp, by = c("Site", "Year", "Season")) -> seas.out

# Add new ID column that combines ID and Site name to ensure the IDs are unique

seas.out %>%
  mutate(UniqueID = interaction(as.character(Site), 
                                ID,
                                sep = "_",
                                drop = T)) %>%
  select(UniqueID, ID, Site, Year, Season, Sex, SCL, Mass, everything()) -> seas.out

# Save variogram plots to pdf and out object to "output.csv"

write.csv(seas.out, 
          "seasonal_output.csv",
          row.names = F)

vario.names <- lapply(seq_along(seas.variograms), function(i){
  temp <- sub(" .*$", "", seas.variograms[[i]]@info$identity)
})

names(seas.variograms) <- vario.names

seas.akde.names <- lapply(seq_along(seas.akde.fits), function(i){
  temp <- sub(" .*$", "", seas.akde.fits[[i]]@info$identity)
})

names(seas.akde.fits) <- seas.akde.names

tmp.var <- seas.variograms[names(seas.variograms) %in% names(seas.akde.fits)]

pdf("seasonal_variograms.pdf")

par(mfrow = c(3,2))
lapply(seq_along(tmp.var), function(i){
  plot(tmp.var[[i]],
       CTMM = seas.akde.fits[[i]],
       level = c(0.5, 0.95),
       fraction = 1,
       main = tmp.var[[i]]@info$identity)
})

dev.off()

######################################################################

## STEP FOUR: CALCULATE HOME RANGE OVERLAPS ##

tmp.overlap <- list()

ids <- unique(sub("_[0-9]{4}", "", names(ann.akdehr)))

for(i in unique(ids)){
  tmp.hr <- ann.akdehr[sub("_[0-9]{4}", "", names(ann.akdehr)) == i]
  if(length(tmp.hr) > 1){
    tmp.overlap <- c(tmp.overlap, overlap(tmp.hr))
  } else {}
}

if(length(tmp.overlap) != 0) {
  ann.akde.overlap.df <- as.data.frame.table(tmp.overlap[names(tmp.overlap) == "CI"][[1]])
  
  for(i in 2:(length(tmp.overlap)/2)){
    ann.akde.overlap.df <- rbind(ann.akde.overlap.df,
                                 as.data.frame.table(tmp.overlap[names(tmp.overlap) == "CI"][[i]]))
  }
  
  dofs <- as.data.frame.table(tmp.overlap[names(tmp.overlap) == "DOF"][[1]])
  
  for(i in 2:(length(tmp.overlap)/2)){
    dofs <- rbind(dofs,
                  as.data.frame.table(tmp.overlap[names(tmp.overlap) == "DOF"][[i]]))
  }
  
  ann.akde.overlap.df %>%
    left_join(y = dofs, by = c("Var1", "Var2")) %>%
    filter(Var3 == "est" & Var1 != Var2) %>%
    mutate(ID = sub("_[0-9]{4}", "", Var1),
           IDY = Var1,
           IDY2 = Var2,
           AKDE.Overlap = Freq.x,
           AKDE.DOF = Freq.y) %>%
    select(ID, IDY, IDY2, AKDE.Overlap, AKDE.DOF) -> ann.akde.overlap
  
  tmp.overlap <- list()
  
  ids <- unique(sub("_[0-9]{4}", "", names(ann.kdehr)))
  
  for(i in unique(ids)){
    tmp.hr <- ann.kdehr[sub("_[0-9]{4}", "", names(ann.kdehr)) == i]
    if(length(tmp.hr) > 1){
      tmp.overlap <- c(tmp.overlap, overlap(tmp.hr))
    } else {}
  }
  
  ann.kde.overlap.df <- as.data.frame.table(tmp.overlap[names(tmp.overlap) == "CI"][[1]])
  
  for(i in 2:(length(tmp.overlap)/2)){
    ann.kde.overlap.df <- rbind(ann.kde.overlap.df,
                                as.data.frame.table(tmp.overlap[names(tmp.overlap) == "CI"][[i]]))
  }
  
  dofs <- as.data.frame.table(tmp.overlap[names(tmp.overlap) == "DOF"][[1]])
  
  for(i in 2:(length(tmp.overlap)/2)){
    dofs <- rbind(dofs,
                  as.data.frame.table(tmp.overlap[names(tmp.overlap) == "DOF"][[i]]))
  }
  
  ann.kde.overlap.df %>%
    left_join(y = dofs, by = c("Var1", "Var2")) %>%
    filter(Var3 == "est" & Var1 != Var2) %>%
    mutate(ID = sub("_[0-9]{4}", "", Var1),
           IDY = Var1,
           IDY2 = Var2,
           KDE.Overlap = Freq.x,
           KDE.DOF = Freq.y) %>%
    select(ID, IDY, IDY2, KDE.Overlap, KDE.DOF) -> ann.kde.overlap
  
  ann.akde.overlap %>%
    left_join(y = ann.kde.overlap, by = c("ID", "IDY", "IDY2")) %>%
    left_join(y = ann.out, by = c("ID", "IDY")) -> ann.overlap
  
  write.csv(ann.overlap,
            file = "annual_overlap.csv",
            row.names = FALSE)
} else {}

# Repeat for seasonal home ranges

tmp.overlap <- list()

ids <- unique(sub("_[0-9]{4}_[^_]+$", "", names(seas.akdehr)))

for(i in unique(ids)){
  tmp.hr <- seas.akdehr[sub("_[0-9]{4}_[^_]+$", "", names(seas.akdehr)) == i]
  if(length(tmp.hr) > 1){
    tmp.overlap <- c(tmp.overlap, overlap(tmp.hr))
  } else {}
}

if(length(tmp.overlap) != 0) {
  seas.akde.overlap.df <- as.data.frame.table(tmp.overlap[names(tmp.overlap) == "CI"][[1]])
  
  for(i in 2:(length(tmp.overlap)/2)){
    seas.akde.overlap.df <- rbind(seas.akde.overlap.df,
                                  as.data.frame.table(tmp.overlap[names(tmp.overlap) == "CI"][[i]]))
  }
  
  dofs <- as.data.frame.table(tmp.overlap[names(tmp.overlap) == "DOF"][[1]])
  
  for(i in 2:(length(tmp.overlap)/2)){
    dofs <- rbind(dofs,
                  as.data.frame.table(tmp.overlap[names(tmp.overlap) == "DOF"][[i]]))
  }
  
  seas.akde.overlap.df %>%
    left_join(y = dofs, by = c("Var1", "Var2")) %>%
    filter(Var3 == "est" & Var1 != Var2) %>%
    mutate(ID = sub("_[0-9]{4}_[^_]+$", "", Var1),
           IDYS = Var1,
           IDYS2 = Var2,
           AKDE.Overlap = Freq.x,
           AKDE.DOF = Freq.y) %>%
    select(ID, IDYS, IDYS2, AKDE.Overlap, AKDE.DOF) -> seas.akde.overlap
  
  tmp.overlap <- list()
  
  ids <- unique(sub("_[0-9]{4}_[^_]+$", "", names(seas.kdehr)))
  
  for(i in unique(ids)){
    tmp.hr <- seas.kdehr[sub("_[0-9]{4}_[^_]+$", "", names(seas.kdehr)) == i]
    if(length(tmp.hr) > 1){
      tmp.overlap <- c(tmp.overlap, overlap(tmp.hr))
    } else {}
  }
  
  seas.kde.overlap.df <- as.data.frame.table(tmp.overlap[names(tmp.overlap) == "CI"][[1]])
  
  for(i in 2:(length(tmp.overlap)/2)){
    seas.kde.overlap.df <- rbind(seas.kde.overlap.df,
                                 as.data.frame.table(tmp.overlap[names(tmp.overlap) == "CI"][[i]]))
  }
  
  dofs <- as.data.frame.table(tmp.overlap[names(tmp.overlap) == "DOF"][[1]])
  
  for(i in 2:(length(tmp.overlap)/2)){
    dofs <- rbind(dofs,
                  as.data.frame.table(tmp.overlap[names(tmp.overlap) == "DOF"][[i]]))
  }
  
  seas.kde.overlap.df %>%
    left_join(y = dofs, by = c("Var1", "Var2")) %>%
    filter(Var3 == "est" & Var1 != Var2) %>%
    mutate(ID = sub("_[0-9]{4}_[^_]+$", "", Var1),
           IDYS = Var1,
           IDYS2 = Var2,
           KDE.Overlap = Freq.x,
           KDE.DOF = Freq.y) %>%
    select(ID, IDYS, IDYS2, KDE.Overlap, KDE.DOF) -> seas.kde.overlap
  
  seas.akde.overlap %>%
    left_join(y = seas.kde.overlap, by = c("ID", "IDYS", "IDYS2")) %>%
    left_join(y = seas.out, by = c("ID", "IDYS")) -> seas.overlap
  
  write.csv(seas.overlap,
            file = "seasonal_overlap.csv",
            row.names = FALSE)
} else {}
