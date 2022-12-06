#TODO: 
#Points that lie on the beginning or end of transects (ex/ 2016, ID: 566)
#Force st_intersect with corresponding transect ID (make sure points are assigned to proper transect... maybe st_distance)
#Add Basin Information (let it vary at site level; deal with transect outside of basin area [maybe estuaries, ex/ 1997, ID: 421, 768])
#Bald eagle data (Code: BAEA)
#2015 data seems to be doubled
#Add sum of area by hexagon and year
#Transect 544 exist but not 545. There are points that exist in transect 545... (1994)
#Bathemtry data NA -> 0 ... check 
#Fix how basin is defined. Issues noticed in hexagon.data

#-Libraries-#
library(sf)
library(tidyverse)
library(lubridate)

#-Functions-#

skip.fun <- function(flags = flags){
  flags <- tryCatch(
    {
      skip.site <- flags %>% group_by(L1) %>% summarise(nsite = max(site.group))
      skip.site <- data.frame(L1 = rep(skip.site$L1, skip.site$nsite), site.group = unlist(sapply(skip.site$nsite, function(x) seq(1,x,1))))
      skip.site <- skip.site[which(!do.call(paste, skip.site) %in% do.call(paste, flags %>% select(L1, site.group))),]
      
      skip.site.ID <- mapply(function(x,y) which(flags$site.group == min(flags$site.group[flags$site.group >= x & flags$L1 == y]) & flags$L1 == y), skip.site$site.group, skip.site$L1)  
      
      skip.loci <- flags[skip.site.ID,]
      skip.loci$x <- skip.loci$site.group - skip.site$site.group
      skip.loci <- skip.loci %>% mutate(seg = seg + site.length * x,
                                        site.group = site.group - x)
      skip.loci <- skip.loci %>% select(-x)
      
      rbind(flags, skip.loci) %>%
        arrange(site.group)
      
    },
    error = function(cond){
      return(flags)
    }
    
  )
  return(flags)
}
multi.fun <- function(x){
  x <- as.data.frame(x)
  x <- as.matrix(x)
  x <- st_linestring(x)
  return(x)
}     
site.fun <- function(obj, site.length = 1000){
  geom <- st_geometry(obj)
  att <- st_drop_geometry(obj)
  
  tran.length <- as.numeric(st_length(geom)) #Transect length
  
  coords <- data.frame(st_coordinates(geom))
  coords <- coords %>% group_by(L1) %>% mutate(ID = seq(1:n()))
  coords <- full_join(coords %>% group_by(L1) %>% filter(ID != n()),
                      coords %>% group_by(L1) %>% filter(ID != 1) %>% mutate(ID = seq(1:n())),
                      by = c("L1", "ID"))
  
  coords <- coords %>% mutate(length = sqrt((X.y - X.x)^2 + (Y.y - Y.x)^2),
                              cumsum = cumsum(length))
  
  coords <- coords %>% mutate(site.group = floor(cumsum/site.length))
  
  coords <- coords %>% group_by(L1) %>% mutate(begin.seg = ifelse(duplicated(site.group), 0, 1))
  
  flags <- coords %>% filter(begin.seg == 1) %>%
    filter(site.group != 0) %>%
    group_by(L1) %>%
    ungroup(L1) %>%
    select(!begin.seg)
  
  flags <- flags %>% mutate(diff.X = X.y - X.x,
                            diff.Y = Y.y - Y.x,
                            seg = cumsum - floor(cumsum/site.length)*site.length)
  
  flags <- skip.fun(flags = flags)
  
  flags <- flags %>% mutate(seg.x = cos(atan(abs(diff.Y)/abs(diff.X))) * seg,
                            seg.y = sin(atan(abs(diff.Y)/abs(diff.X))) * seg,
                            X.z = ifelse(is.na(log(diff.X)),
                                         X.y + seg.x,
                                         X.y - seg.x),
                            Y.z = ifelse(is.na(log(diff.Y)),
                                         Y.y + seg.y,
                                         Y.y - seg.y))
  
  
  
  new.loci <- flags %>% select(L1, ID, site.group, X.z, Y.z) %>%
    rename(X = X.z, Y = Y.z)  %>%
    mutate(ID = ID - 0.5)
  
  site <- left_join(data.frame(st_coordinates(geom)) %>%
                      group_by(L1) %>% mutate(ID = seq(1:n()) -  1) %>%
                      ungroup(L1),
                    coords %>% select(L1, ID, site.group),
                    by = c("L1", "ID")) %>%
    mutate(site.group = ifelse(is.na(site.group), 0, site.group))
  
  
  site <- rbind(site, new.loci) %>% 
    arrange(L1, site.group, ID) %>%
    mutate(L1 = factor(L1),
           site.group = factor(site.group)) %>%
    group_by(L1) %>%
    mutate(ID = seq(1:n())) %>%
    ungroup(L1)
  
  site <- left_join(site %>% group_by(L1) %>% filter(ID != n()),
                    site %>% group_by(L1) %>% filter(ID != 1) %>% mutate(ID = seq(1:n())) %>% select(!site.group),
                    by = c("L1", "ID"))
  
  site <- site %>% mutate(length = sqrt((X.y - X.x)^2 + (Y.y - Y.x)^2),
                          cumsum = cumsum(length),
                          site.group = as.numeric(site.group)) %>% 
    group_by(L1, site.group) %>% 
    mutate(site.cumsum = sum(length)) %>% 
    ungroup(L1, site.group)
  
  site <- site %>% group_by(L1) %>% mutate(site.group = ifelse(tran.length[L1] > 1000 & site.cumsum < 500, site.group - 1, site.group))
  
  site <- site %>% 
    group_by(L1, site.group) %>% 
    mutate(site.cumsum = sum(length)) %>% 
    ungroup(L1, site.group)
  
  site <- rbind(site %>% select(X.x, Y.x, L1, site.group, ID) %>% rename(X = X.x, Y = Y.x),
                site %>% select(X.y, Y.y, L1, site.group, ID) %>% rename(X = X.y, Y = Y.y)) %>%
    arrange(L1, site.group, ID) %>%
    select(!ID)
  
  site <- site[!duplicated(site), ]
  
  n.sites <- site %>% group_by(L1) %>% summarise(nsite = n_distinct(site.group)) %>% select(nsite) %>% .$nsite
  site <- st_sfc(lapply(split(site[,1:2], f = ~ site$site.group + site$L1, drop = T), multi.fun), crs = st_crs("EPSG:26910"))
  site <- st_set_geometry(att[rep(1:nrow(att), n.sites), ], site)
  site <- site %>% mutate(Shape_Length = st_length(.))
  return(site)
}
site.id.fun <- function(x, y){
  tmp <- data %>% filter(TransectID == x)
  siteID <- tmp[min(unlist(st_intersects(y, tmp))), "siteID"]
  siteID <- st_drop_geometry(siteID)
  return(siteID)
}

#-Directory-#
dsn <- "/Users/farrm/OneDrive - UW/Projects/Seaduck/Data/MidwinterAerialSeabirdSurveys/MidwinterAerialSeabirdSurveys.gdb"

#-Load data-#
layer.names <- st_layers(dsn)$name

for(i in 1:length(layer.names)){
  object <- with_tz(st_read(dsn, layer.names[i]), "UTC") %>% force_tz("America/Los_Angeles")
  assign(paste(layer.names[i]), object)
  rm(object)
}



#-Species information-#

#Common name of sea ducks
commonnames <- c("Red-breasted merganser", 
                 "Common merganser",
                 "Hooded merganser",
                 "Bufflehead",
                 "Barrow's goldeneye",
                 "Common goldeneye",
                 "White-winged scoter",
                 "American (Black) Scoter",
                 "Surf scoter",
                 "Harlequin duck",
                 "Long-tailed duck"
) 

#Extract species codes
sppcodes <- Species$PSEMP_SpeciesCode[Species$TaxoCommonName %in% commonnames]

#-Observation data-#
# files <- list.files(path = "./PSAMP_MattFarr_Delivery/Data/DataFormattedForStatistics",
#            pattern = "ObsSurvey",
#            full.names = TRUE,
#            recursive = TRUE)
# 
# observation.data <- data.frame()
# 
# for(i in 1:length(files)){
#   observation.data <- rbind(observation.data, read_csv(file = files[i]) %>%
#                               mutate(SurveyYear = as.numeric(substr(files[i], 59, 62))) %>%
#                               select(Species, Count, TransectID, SurveyYear, DateTime, Lat, Lon))
# }
# 
# observation.data <- rbind(st_as_sf(observation.data %>% filter(SurveyYear < 2012), coords = c("Lon", "Lat"), crs = st_crs("EPSG:4267")) %>%
#                              st_transform(., crs = st_crs("EPSG:26910")),
#                            st_as_sf(observation.data %>% filter(SurveyYear >= 2012), coords = c("Lon", "Lat"), crs =st_crs("EPSG:4326")) %>%
#                              st_transform(., crs = st_crs("EPSG:26910")))
# #observation.data <- st_transform(observation.data, crs = st_crs("EPSG:2927"))
# #observation.data <- st_transform(observation.data, crs = st_crs("EPSG:26910"))
# 
# observation.data$Species <- recode(observation.data$Species, OLDS = "LTDU")
# observation.data <- observation.data %>% filter(Species %in% sppcodes)

#-Segment transects-#

year <- PSEMP_SurveyRoutes %>% 
  st_drop_geometry(.) %>%
  summarise(year = unique(SurveyYear)) %>%
  mutate(year = as.numeric(year)) %>%
  arrange(year) %>% .$year

nyears <- length(year)

#Basin information
basin.data <- PSEMP_Analysis_Strata %>%
  st_zm(., drop = TRUE) %>%
  st_transform(., crs = st_crs("EPSG:26910")) %>%
  st_cast(., "MULTIPOLYGON") %>%
  st_make_valid(.) %>%
  group_by(Basin) %>%
  summarise(Shape = st_union(Shape)) %>%
  mutate(Shape_Area = units::set_units(st_area(.), km^2))

#-EPA grid-#
grid <- st_read("/Users/farrm/OneDrive - UW/Projects/Seaduck/Data/EPA40kmHexagons/EPA_Hexagons_40km/EPA_Hexagons_40km.shp")
grid <- st_transform(grid, crs = st_crs("EPSG:26910"))
grid <- st_crop(grid, st_bbox(basin.data %>% summarise(Shape = st_union(Shape)) %>% st_buffer(., dist = 15000)))
#grid <- st_intersection(grid, basin.data %>% summarise(Shape = st_union(Shape)) %>% st_buffer(., dist = 15000))
grid %>% mutate(Shape_Area = units::set_units(st_area(.), km^2))
grid$ID <- seq(1,dim(grid)[1])

#Initiate transect and observation data objects
transect.data <- obs.data <- data.frame()

for(t in 1:nyears){
  obj <- PSEMP_SurveyRoutes %>%
    st_zm(., drop = TRUE) %>%
    filter(SurveyYear== year[t]) %>%
    filter(duplicated(TransectID) == F) %>%
    st_transform(., crs = st_crs("EPSG:26910")) %>%
    st_cast(., "LINESTRING") %>%
    mutate(Shape_Length = st_length(.))
  
  obj <- obj %>% st_buffer(., dist = 88, endCapStyle = "FLAT") %>%
    st_intersection(., basin.data) %>%
    mutate(Shape_Area = st_area(.)) %>%
    arrange(TransectID, -Shape_Area) %>%
    filter(duplicated(TransectID) == F) %>%
    st_drop_geometry(.) %>%
    right_join(., obj, by = c("TransectID", "LAST_DateTime", "FirstDateTime", "TransectCode",
                              "SurveyYear", "Shape_Length")) %>%
    select(-Shape_Area) %>%
    st_set_geometry(., value = "Shape")
  
  #data <- transect.fun(st_geometry(obj), st_drop_geometry(obj))
  #coords <- data.frame(st_coordinates(st_geometry(obj))) %>% group_by(L1) %>% mutate(ID = seq(1:(n())))
  #data <- st_set_geometry(data, st_sfc(mapply(coords.fun, data$tranID, data$begin, data$end)))
  #st_crs(data) <- st_crs("EPSG:26910")
  #data <- data %>% mutate(Shape_Length = st_length(.))
  
  data <- site.fun(obj = obj)
  
  data <- data %>% st_buffer(., dist = 88, 
                             endCapStyle = "FLAT") %>%
    mutate(Shape_Area = units::set_units(st_area(.), km^2))
  
  #data$siteID <- nsite:(nsite + dim(data)[1] - 1)
  #nsite <- nsite + dim(data)[1]
  
  data$siteID <- 1:dim(data)[1]
  
  # data$basin <- st_intersection(data, basin.data) %>% 
  #   mutate(Shape_Area = st_area(.)) %>% 
  #   arrange(siteID, -Shape_Area) %>% 
  #   #filter(duplicated(siteID) == F) %>% 
  #   select(Basin) %>% .$Basin
  
  data$hexID <- st_intersection(data, grid) %>% 
      mutate(Shape_Area = st_area(.)) %>%
      arrange(siteID, -Shape_Area) %>%
      filter(duplicated(siteID) == F) %>%
      select(ID) %>% .$ID
  
  transect.data <- rbind(transect.data, data)
  
  #data$tran.segID <- 1:dim(data)[1]
  
  # obj <- observation.data %>%
  #   filter(SurveyYear == year[t])
  # 
  # obj$tran.segID <- max.col(st_intersects(obj, data), "first")
  # 
  # obj <- obj %>%
  #   left_join(., data %>%
  #              st_drop_geometry(.) %>%
  #              select("SurveyYear", "TransectID", "tran.segID", "siteID", "Basin"),
  #              by = c("SurveyYear", "TransectID", "tran.segID")) %>%
  #   select(-"tran.segID")
  

  

  obj <- PSEMP_SurveyObservations %>%
    st_zm(., drop = TRUE) %>%
    filter(SurveyYear== year[t]) %>%
    st_transform(., crs = st_crs("EPSG:26910"))

  #obj$tran.segID <- max.col(st_intersects(obj, data), "first")
  #obj$siteID <- max.col(st_intersects(obj, data), "first")
  obj <- obj %>% rowwise() %>% mutate(site.id.fun(TransectID, Shape))
  

  obj <- obj %>%
    filter(PSEMP_SpeciesCode %in% sppcodes) %>%
    left_join(., data %>%
                st_drop_geometry(.) %>%
                #select("SurveyYear", "TransectID", "tran.segID", "siteID"),
                select("SurveyYear", "TransectID", "siteID"),
              #by = c("SurveyYear", "TransectID", "tran.segID")) %>%
              by = c("SurveyYear", "TransectID", "siteID")) 
    #select(-"tran.segID")

  #observation.data <- rbind(observation.data, obj)
  
  obs.data <- rbind(obs.data, obj)
  
}

transect.data <- transect.data %>% mutate(Basin = as.numeric(as.factor(Basin)),
                                          SurveyYear = as.numeric(as.factor(SurveyYear)),
                                          hex.factor = as.numeric(as.factor(hexID)))

transect.data <- transect.data %>% drop_na(Basin)

#Hexagon data sampled by transects
hexagon.data <- grid[unique(transect.data$hexID),]
hexagon.data <- hexagon.data %>% arrange(ID)
hexagon.data$basin <- apply(table(transect.data$hexID, transect.data$Basin), 1, which.max)
hexagon.data <- hexagon.data %>% arrange(ID)

st_write(hexagon.data, "./DataFormatting/hexagon.data.shp", delete_layer = T)

nhex <- dim(hexagon.data)[1]

prediction.data <-  st_intersection(hexagon.data, basin.data) %>%
  mutate(km2 = units::set_units(st_area(.), km^2))

prediction.data <- prediction.data %>% arrange(ID)

st_write(prediction.data, "./DataFormatting/prediction.data.shp", delete_layer = T)

npred <- dim(prediction.data)[1]

obs.data <- obs.data %>% drop_na(siteID) #%>% drop_na(Basin)
obs.num <- obs.data %>% #st_drop_geometry(.) %>%
  #select(SurveyYear, siteID, Basin, Species, Count) %>%
  select(SurveyYear, siteID, PSEMP_SpeciesCode, ObservationCount) %>%
  #rename(year = SurveyYear, site = siteID, basin = Basin, species = Species, count = Count) %>%
  rename(year = SurveyYear, site = siteID, species = PSEMP_SpeciesCode, count = ObservationCount) %>%
  mutate(year = as.numeric(factor(year)),
         species = as.numeric(factor(species)))
         #basin = as.numeric(as.factor(basin)))

nsites <- transect.data %>% st_drop_geometry(.) %>% group_by(SurveyYear) %>% summarise(nsite = n()) %>% .$nsite

#-Observation data-#
obs.array <- array(data = NA, dim = c(length(sppcodes), nyears, max(nsites)))

#Set abundances to zero
for(t in 1:nyears){
  obs.array[,t,1:nsites[t]] <- 0
}

#Add observed counts
for(i in 1:dim(obs.num)[1]){
  obs.array[obs.num$species[i], obs.num$year[i], obs.num$site[i]] <- obs.num$count[i]
}

#-Covariate data-#

#Pacific Decadal Oscillation
PDO <- read.table(file = "/Users/farrm/OneDrive - UW/Projects/Seaduck/Data/PDO_data.dat", header = TRUE, skip = 1, fill = T)
PDO <- reshape2::melt(PDO, id.vars = "Year")
colnames(PDO) <- c("Year", "Month", "PDO")
PDO[PDO$Year == "2022" & PDO$Month == "Jul", "PDO"] <- -2.48
PDO$PDO <- as.numeric(PDO$PDO)
PDO <- as.data.frame(PDO)
PDO <- PDO %>% filter(Year > 1992 & Year < 2022) %>% group_by(Year) %>% summarise(Value = mean(PDO, na.rm = T))
PDO <- (PDO$Value - mean(PDO$Value))/sd(PDO$Value)

#Bathymetry
bathymetry <- raster::raster("/Users/farrm/OneDrive - UW/Projects/Seaduck/Data/Marine Environmental Layers/Final_EnvLayers/Bathymetry.tif")
bathymetry <- raster::projectRaster(bathymetry, crs = raster::crs("EPSG:26910"))
bathymetry <- raster::crop(bathymetry, st_bbox(basin.data %>% summarise(Shape = st_union(Shape)) %>% st_buffer(., dist = 15000)))

depthhex <- as.numeric(hexagon.data %>% mutate(depth = raster::extract(bathymetry, hexagon.data, fun = mean, na.rm = T)) %>% .$depth)
depthhex[depthhex>0] <- 0
depthhex <- abs(depthhex)
depth.hex <- (depthhex - mean(depthhex, na.rm = T))/sd(depthhex, na.rm = T)

depthsite <- as.numeric(transect.data %>% mutate(depth = raster::extract(bathymetry, transect.data, fun = mean, na.rm = T)) %>% .$depth)
depthsite[is.na(depthsite)] <- 0
depthsite[depthsite>0] <- 0
depthsite <- abs(depthsite)
#depth.site <- (depth.site - mean(depth.site, na.rm = T))/sd(depth.site, na.rm = T)

depthpred <- as.numeric(prediction.data %>% mutate(depth = raster::extract(bathymetry, prediction.data, fun = mean, na.rm = T)) %>% .$depth)
depthpred[is.na(depthpred)] <- 0
depthpred[depthpred>0] <- 0
depthpred <- abs(depthpred)
#depth.pred <- (depthpred - mean(depthpred, na.rm = T))/sd(depthpred, na.rm = T)

#Basin ID
basin <- array(data = NA, dim = c(nyears, max(nsites)))

#Area sampled per site
depth.site <- hexagon <- area <- array(data = NA, dim = c(nyears, max(nsites)))

for(i in 1:dim(transect.data)[1]){
  t <- transect.data$SurveyYear[i]
  j <- transect.data$siteID[i]
  basin[t,j] <- transect.data$Basin[i]
  area[t,j] <- as.numeric(transect.data$Shape_Area[i])
  hexagon[t,j] <- transect.data$hex.factor[i]
  depth.site[t,j] <- depthsite[i]
}

depth.pred <- (t(matrix(rep(depthpred, nyears), ncol = nyears)) - matrix(rep(apply(depth.site, MARGIN = 1, mean, na.rm = T), npred), ncol = npred, nrow = nyears))/matrix(rep(apply(depth.site, MARGIN = 1, sd, na.rm = T), npred), ncol = npred, nrow = nyears)

depth.site <- t(scale(t(depth.site)))

nbasins <- max(basin, na.rm = T)

#Change of support
COS <- prediction.data$km2
COS <- as.numeric(COS/mean(area, na.rm = T))

#COS2 <- 100 #density per 10 km^2

COS.hex <- hexagon.data$km2
COS.hex <- as.numeric(COS.hex/mean(area, na.rm = T))

#Site area offset
area <- area/mean(area, na.rm = T)

#Hexagon area offset
hex.area.df <- transect.data %>% st_drop_geometry(.) %>% 
  group_by(SurveyYear, hex.factor) %>% 
  summarise(hex.area = sum(Shape_Area))

hex.area <- matrix(0, nrow = nyears, ncol = nhex)
for(i in 1:dim(hex.area.df)[1]){
  t <- hex.area.df$SurveyYear[i]
  h <- hex.area.df$hex.factor[i]
  hex.area[t,h] <- hex.area.df$hex.area[i]
}

hex.area <- (hex.area - mean(hex.area))/sd(hex.area)

#Basin hexagon ID
basinhex <- hexagon.data$basin

#Basin prediction ID
basinpred <- prediction.data$basin

#Hexagon prediction ID
hexpred <- as.numeric(as.factor(prediction.data$ID))

#Effort per year
effort <- nsites/mean(nsites) #MTF: should this be the median?

#-Compile Data-#
data.list <- list(n = obs.array, N = apply(obs.array, MARGIN = c(1,2), sum, na.rm = T))
                  
con.list <- list(nspecies = length(sppcodes), nyears = nyears, nsites = nsites, nbasins = nbasins, nhex = nhex, nped = nped,
                 basin = basin, area = area, effort = effort, hex.area = hex.area, mean.area = mean.area,
                 hexagon = hexagon, basinhex = basinhex, basinpred = basinpred, hexpred = hexpred, PDO = PDO, 
                 depth.site = depth.site, depth.hex = depth.hex, depth.pred = depth.pred, COS = COS, COS.hex = COS.hex)

save(data.list, file = "./DataFormatting/FormattedData.Rds")
save(con.list, file = "./DataFormatting/FormattedCon.Rds")
