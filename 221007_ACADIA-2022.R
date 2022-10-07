### ACADIA 2022
### Date: 24-26 Sept
### Workshop: A data-driven approach for urban design and master planning development

# 1.- Loading/Installing necessary packages:

if ("sp" %in% rownames(installed.packages()) == FALSE)
{
  install.packages("sp")
  library(sp)
} else {
  library(sp)
}

if ("sf" %in% rownames(installed.packages()) == FALSE)
  {
  install.packages("sf")
  library(sf)
  sf::sf_use_s2(FALSE) # This specific command avoids the package "sf" (Simple Features) to use parent package "sp" (Spatial data) spherical geometries 
  } else {
  library(sf)
  sf::sf_use_s2(FALSE)
  }

if ("raster" %in% rownames(installed.packages()) == FALSE)
{
  install.packages("raster")
  library(raster)
} else {
  library(raster)
}

if ("osmar" %in% rownames(installed.packages()) == FALSE)
{
  install.packages("osmar")
  library(osmar)
} else {
  library(osmar)
}

if ("osmdata" %in% rownames(installed.packages()) == FALSE)
{
  install.packages("osmdata")
  library(osmdata)
} else {
  library(osmdata)
}

if ("tmaptools" %in% rownames(installed.packages()) == FALSE)
{
  install.packages("tmaptools")
  library(tmaptools)
} else {
  library(tmaptools)
}

if ("leaflet" %in% rownames(installed.packages()) == FALSE)
{
  install.packages("leaflet")
  library(leaflet)
} else {
  library(leaflet)
}

if ("tidygraph" %in% rownames(installed.packages()) == FALSE)
{
  install.packages("tidygraph")
  library(tidygraph)
} else {
  library(tidygraph)
}

if ("plyr" %in% rownames(installed.packages()) == FALSE)
{
  install.packages("plyr")
  library(plyr)
} else {
  library(plyr)
}

if ("dplyr" %in% rownames(installed.packages()) == FALSE)
{
  install.packages("dplyr")
  library(dplyr)
} else {
  library(dplyr)
}

if ("reshape2" %in% rownames(installed.packages()) == FALSE)
{
  install.packages("reshape2")
  library(reshape2)
} else {
  library(reshape2)
}

if ("jsonify" %in% rownames(installed.packages()) == FALSE)
{
  install.packages("jsonify")
  library(jsonify)
} else {
  library(jsonify)
}

if ("dbscan" %in% rownames(installed.packages()) == FALSE)
{
  install.packages("dbscan")
  library(dbscan)
} else {
  library(dbscan)
}

if ("purrr" %in% rownames(installed.packages()) == FALSE)
{
  install.packages("purrr")
  library(purrr)
} else {
  library(purrr)
}

if ("RColorBrewer" %in% rownames(installed.packages()) == FALSE)
{
  install.packages("RColorBrewer")
  library(RColorBrewer)
} else {
  library(RColorBrewer)
}

if ("nomisr" %in% rownames(installed.packages()) == FALSE)
{
  install.packages("nomisr")
  library(nomisr)
} else {
  library(nomisr)
}

# 2.- Setting up your working folder:
# R language does not recognize back slash as they are reserved for special characters, do please change back slashes on your folder path by a single or double forward slash;
    # The "forward slash" / is actually more common as it used by Unix, Linux, and macOS
    # The "backward slash" \ is actually somewhat painful as it is also an escape character. So whenever you want one, you need to type two in string: "C:\\TEMP".

setwd("C:/Users/Jorge/Documents/Grimshaw/DT/ACADIA 2022")
getwd()

# 3.- Setting up our area of analysis:

longitude <- -0.08897
latitude <- 51.51335 
Analysis_rad <- 1000

BBox_Analysis <- osmar::center_bbox(longitude, latitude, Analysis_rad, Analysis_rad) # defining A bounding Box
BBox_poly <- sf::st_as_sf(tmaptools::bb_poly(BBox_Analysis)) # converting the bounding box to sf polygon 

# Plotting BBox_poly on leaflet
leaflet::leaflet() %>%
  leaflet::addTiles() %>%
  leaflet::addPolygons(data = BBox_poly)

# 4.- Loading all necessary data
  
# OSM_DF <- BBox_Analysis %>% osmdata::opq() # un-comment this line in order to generate a query to the overpass API (Open Street Maps)
# osmdata_xml(OSM_DF, filename = "osm_df.osm") # un-comment this line in order to send the query and receive the desired information on .osm format

# Loading OSM data
OSM_df <- osmdata::osmdata_sf(doc = "osm_df.osm", stringsAsFactors = FALSE)

# Understanding OSM data structure: https://wiki.openstreetmap.org/wiki/Map_features

osm_pol <- OSM_df$osm_polygons # "$" symbol allows you to access different levels on your data structure, in the case of a data-frame (table) it allows you to access different columns
osm_line <- OSM_df$osm_lines
osm_points <- OSM_df$osm_points
osm_mpoly <- OSM_df$osm_multipolygons
osm_mline <- OSM_df$osm_multilines

# Loading LIDAR data

# Digital Surface Model
DSM <- raster::raster("Lidar_National_Program/National-LIDAR-Programme-DSM-2020-TQ38sw/P_12151/DSM_TQ3080_P_12151_20201212_20201212.tif") # loading raster file DIGITAL SURFACE MODEL (DSM)
BBox_poly <- sf::`st_crs<-`(BBox_poly, 4326) # assigning crs to BBox (lat, lng coordinates can be understood as crs mercator, epsg: 4326)
DSM <- raster::crop(DSM, raster::extent(sf::as_Spatial(sf::st_transform(sf::st_as_sf(BBox_poly), 27700)))) # cropping the raster with our bounding box, bear in mind that we are re-projecting the bbox to epsg:27700 crs and migrating the data from sf to sp system as the sp package is the base of the raster package

# Digital Terrain Model
DTM <- raster::raster("Lidar_National_Program/National-LIDAR-Programme-DTM-2020-TQ38sw/P_12151/DTM_TQ3080_P_12151_20201212_20201212.tif") # loading raster file DIGITAL TERRAIN MODEL (DTM)
DTM <- raster::crop(DTM, raster::extent(sf::as_Spatial(sf::st_transform(sf::st_as_sf(BBox_poly), 27700)))) # Cropping raster

# Plotting the raster - please uncomment and run to see plots

#slope_DSM <- raster::terrain(DSM, opt = 'slope') # slope map
#aspect_DSM <- raster::terrain(DSM, opt = 'aspect') # height aspect
#hill_DSM <- raster::hillShade(slope_DSM, aspect_DSM, 40, 270) # hill shade
#raster::plot(hill_DSM, col = grey(0:100/100), legend = FALSE, main = 'Analysis Area') # plot 1
#raster::plot(DSM, col = rainbow(25, alpha = 0.35), add = TRUE) # plot 2
#slope_DTM <- raster::terrain(DTM, opt = 'slope') # slope map
#aspect_DTM <- raster::terrain(DTM, opt = 'aspect') # height aspect
#hill_DTM <- raster::hillShade(slope_DTM, aspect_DTM, 40, 270) # hill shade
#raster::plot(hill_DTM, col = grey(0:100/100), legend = FALSE, main = 'Analysis Area') # plot 3
#raster::plot(DTM, col = rainbow(25, alpha = 0.35), add = TRUE) # plot 4

# Loading Ordnance Survey data: https://osdatahub.os.uk/downloads/open
OS_Buildings <- sf::read_sf("Ordenance_survey/vmdvec_tq/OS VectorMap District (ESRI Shape File) TQ/data/TQ_Building.shp")

# 5.- Processing and visualizing loaded data

# Buildings_sf <- osm_pol[!is.na(osm_pol$building), ] # Accessing generic data frame values in R: DF[Row, Column]
# leaflet::leaflet() %>% leaflet::addTiles() %>% leaflet::addPolygons(data = Buildings_sf, popup = ~ building)

# working with different projection systems
sf::st_crs(OS_Buildings) # checking CRS projection system on ordinance survey spatial data
OS_Buildings <- sf::st_transform(OS_Buildings, 4326) # leaflet and multiple mapping packages only accepts mercator (epsg: 4326) crs system, that is the one used by Google maps and other global map services

sf::st_crs(BBox_poly)
OS_Buildings <- sf::st_intersection(OS_Buildings, BBox_poly) # st_intersection correspond to "clip" command on QGIS, if we do not have both data sets in the same projection system it will not work

OS_Buildings <- sf::st_zm(OS_Buildings, drop = TRUE, what = "ZM") # leaflet only accepts planar geometry (no Z values), rather than mapbox that is capable of plotting 3D spatial data, we use "st_zm" to drop all z axis values 
#leaflet::leaflet() %>% leaflet::addTiles() %>% leaflet::addPolygons(data = OS_Buildings)

# Extracting desired OSM data - Buildings

# plot(osm_pol[unlist(sf::st_contains(BBox_poly, sf::st_centroid(osm_pol))), ][, "geometry"]) # rather than "clipping" lets plot all polygons which it is center point fells within our analysis area (bbox)

OSM_Buildings <- osm_pol[unlist(sf::st_contains(BBox_poly, sf::st_centroid(osm_pol))), ] # we use st_contains to extract only the buildings within the defined building box
OSM_Buildings <- OSM_Buildings[colSums(!is.na(OSM_Buildings)) > 0] # cleaning the data frame so we obtain only the columns with non NA information

# a list of all key and values for the tag building can be found here: https://wiki.openstreetmap.org/wiki/Key:building

keep_columns <- c("osm_id", "addr.city", "addr.country", "addr.housename", "addr.housenumber", "addr.postcode", "addr.street", "addr.unit", "amenity.disused", "disused.amenity", "building", "building.colour", "building.fireproof", "building.flats", "building.levels", "building.architecture", "building.material", "building.min_level", "building.parts", "building.part", "building.soft_storey", "entrance", "height", "max_level", "min_level", "non_existent_levels", "start_date", "amenity", "building.levels.underground", "building.name", "building.use", "addr.buildingnumber", "note.building", "note.building.part", "listed_building", "demolished.building") # list of columns with relevant building information as per OSM wiki page
OSM_Buildings <- OSM_Buildings[, which(colnames(OSM_Buildings) %in% keep_columns)] # filtering columns based on the list above

OSM_Buildings$check <- NA # we create a new column with NA values in order to check if the amenity column have value in it and therefore it could contain a building

# i in seq_along(data.frame$column) : i is going to be the number of the row from row 1 to the final row
# from all the columns we want to operate in all but: the first (osm_id), the last two (geometry and check), so we are going to perform the filtering based on columns: 2 to column length - 2

for (i in seq_along(OSM_Buildings$amenity)) {
  if (sum(!is.na(OSM_Buildings[i, 2:(length(OSM_Buildings) - 2)])) == 1) {OSM_Buildings$check[i] <- NA} else {OSM_Buildings$check[i] <- "YES"}
}

# plot(OSM_Buildings[!is.na(OSM_Buildings$check), "geometry"]) # we can see here that we are still missing some buildings, it they are not on the polygons they should be on the multi-polygons data-frame!

OSM_Buildings <- OSM_Buildings[!is.na(OSM_Buildings$check), ] # final filtering based on check column

# we are going to perform the same operation on the multi-polygon data-frame

OSM_Buildings_m <- osm_mpoly[unlist(sf::st_contains(BBox_poly, sf::st_centroid(osm_mpoly))), ]
OSM_Buildings_m <- OSM_Buildings_m[colSums(!is.na(OSM_Buildings_m)) > 0]
# https://wiki.openstreetmap.org/wiki/Key:building
keep_columns <- c("building", "building.colour", "building.fireproof", "building.flats", "building.levels", "building.architecture", "building.material", "building.min_level", "building.part", "building.soft_storey", "entrance", "height", "max_level", "min_level", "non_existent_levels", "start_date", "amenity")
OSM_Buildings_m <- OSM_Buildings_m[, which(colnames(OSM_Buildings_m) %in% keep_columns)]

OSM_Buildings_m$check <- NA

for (i in seq_along(OSM_Buildings_m$amenity)) {
  if (sum(!is.na(OSM_Buildings_m[i, 1:7])) == 1) {OSM_Buildings_m$check[i] <- NA} else {OSM_Buildings_m$check[i] <- "YES"}
}

OSM_Buildings_m <- OSM_Buildings_m[!is.na(OSM_Buildings_m$check), ]

OSM_Buildings_m <- OSM_Buildings_m[OSM_Buildings_m$building != "", ] # filtering based on building column, those which do not have a value on the building column are going out 
OSM_Buildings_m <- sf::st_cast(OSM_Buildings_m, "POLYGON") # if we want to join both data frames we need to unify the geometry format, with st_cast we can transform multi-polygons into polygons or to other geometry format as sp

OSM_Build_27700 <- plyr::rbind.fill(OSM_Buildings, OSM_Buildings_m) # rbind stands for row bind, so we are going just to join both data frames a fill non matching columns with NA ".fill"
OSM_Build_27700 <- sf::st_transform(sf::st_as_sf(OSM_Build_27700), 27700) # in order to export it to rhino or other cad/3d/GIS software we need to re-project it in its CRS 

#plot(OSM_Build_27700$geometry) now it seems we have all the buildings on the area.

## Extracting Heights from Raster 

osm_hfpMax <- raster::extract(DSM, OSM_Build_27700, fun = max, na.rm = TRUE) # argument fun, can be updated to other mathematical functions as min, mean or other custom calculations
osm_hfpMin <- raster::extract(DTM, OSM_Build_27700, fun = max, na.rm = TRUE)

norm_height <- osm_hfpMax - osm_hfpMin # normalized height (height from DSM and DTM are based on sea level, we should subtract one from the other in order to obtain relative building height or building height over terrain)

OSM_Build_27700$n_height <- round(norm_height, 1) # we can now include the building height results into our building data frame
OSM_Build_27700$t_height <- NA  # we create a new column to check the height data from Open street maps

# unique(OSM_Build_27700$height)
# OSM height data, should be numeric, unfortunately we can see that there are text (strings) on our column, and we need to extract the numeric data out of it

# gsub("([0-9]+)*$", "\\1", OSM_Build_27700$height) # "([0-9]+)*$" <- looking for numeric data, "\\1" <- replace the data by the number if found

for (i in seq_along(OSM_Build_27700$t_height)) {
  if (!is.na(OSM_Build_27700$height[i])) { OSM_Build_27700$t_height[i] <- as.numeric(gsub("([0-9]+)*$", "\\1", OSM_Build_27700$height[i])) } else { OSM_Build_27700$t_height[i] <- OSM_Build_27700$n_height[i] } # if there is data from OSM we keep this data rather the one from the raster
}

OSM_Build_27700$n_height <- unlist(as.list(OSM_Build_27700$n_height))

for (i in seq_along(OSM_Build_27700$t_height)) {
  if (is.na(OSM_Build_27700$t_height[i])) { OSM_Build_27700$t_height[i] <- OSM_Build_27700$n_height[i] } 
}

# Extracting desired OSM data - Street Network

df_lines <- osm_line

Cropped_lines <- sf::st_crop(df_lines, sf::st_bbox(BBox_poly))
Cropped_lines <- Cropped_lines[!is.na(Cropped_lines$highway), ]
St_network <- data.frame(highway = Cropped_lines$highway, geometry = Cropped_lines$geometry)
St_network_27700 <- sf::st_transform(sf::st_as_sf(St_network), 27700)

Cropped_pt <- sf::st_crop(osm_points, sf::st_bbox(BBox_poly))
Cropped_pt <- Cropped_pt[!is.na(Cropped_pt$highway), ]
St_network_Feat <- data.frame(highway = Cropped_pt$highway, geometry = Cropped_pt$geometry)
St_network_Feat_27700 <- sf::st_transform(sf::st_as_sf(St_network_Feat), 27700)

# Extracting desired OSM data - POS

natural_leisure <- osm_pol[, "leisure"] %>% tidygraph::filter(!is.na(leisure))

natural_leisure$natural <- NA
natural_leisure$waterway <- NA
natural_leisure <- natural_leisure[,order(colnames(natural_leisure))]
natural_leisure <- as.data.frame(natural_leisure)
rownames(natural_leisure) <- NULL

natural_natural <- osm_pol[, "natural"] %>% tidygraph::filter(!is.na(natural))
    
natural_natural$leisure <- NA
natural_natural <- natural_natural[,order(colnames(natural_natural))]
natural_natural <- as.data.frame(natural_natural)
rownames(natural_natural) <- NULL
  
natural <- plyr::rbind.fill(natural_leisure, natural_natural)
natural <- natural[!is.na(natural$geometry),]
natural <- as.data.frame(natural)

natural_27700 <- sf::st_transform(st_as_sf(natural), 27700)

# Extracting desired OSM data - POI

df_points <- osm_points

COLS <- colnames(osm_pol)

amenity_list <- c("yes","bar", "biergarten", "cafe", "fast_food", "food_court", "pub", "restaurant", "college","waste_disposal", "driving_school", "kindergarten", "language_school", "library", "toy_library", "music_school", "school", "university", "bank", "clinic", "dentist", "doctors", "hospital", "nursing_home", "pharmacy", "social_facility", "veterinary", "arts_centre", "brothel", "casino", "cinema", "comunity_centre", "events_venue", "gambling", "love_hotel", "nightclub", "planetarium", "social_centre", "studio", "swingerclub", "theatre", "courthouse", "embassy", "fire_station", "police", "post_office", "post_depot", "prison", "ranger_station", "townhall", "shelter", "animal_boarding", "crematorium", "dive_centre", "funeral_hall", "internet_cafe", "kitchen", "monastery", "place_of_worship", "public_bath")
historic_list <- c("yes", "aqueduct", "building", "castle", "castle_wall", "church", "city_gate", "citywalls", "farm", "fort", "manor", "monastery", "monument", "optical_telegraph", "ruins", "tower")

amenity_pol <- osm_pol[, "amenity"] %>% tidygraph::filter(!is.na(amenity))
amenity_filt <- amenity_pol[amenity_pol$amenity %in% amenity_list, ]

amenity_filt$historic <- NA
amenity_filt$building <- NA
amenity_filt <- as.data.frame(amenity_filt)
rownames(amenity_filt) <- NULL
amenity_filt <- amenity_filt[,order(colnames(amenity_filt))]

historic_pol <- osm_pol[, "historic"] %>% tidygraph::filter(!is.na(historic))
historic_filt <- historic_pol[historic_pol$historic %in% historic_list, ]

historic_filt$amenity <- NA
historic_filt$building <- NA
historic_filt <- as.data.frame(historic_filt)
rownames(historic_filt) <- NULL
historic_filt <- historic_filt[,order(colnames(historic_filt))]

#ACTIVITY POINTS: POINTS data set
pts <- sf::st_crop(df_points, sf::st_bbox(BBox_poly))

#ACTIVITY POINTS: POLYGON data set

leisure_pol <- osm_pol[, c("leisure", "geometry")] %>% tidygraph::filter(!is.na(leisure))

office_pol <- osm_pol[, c("office", "geometry")] %>% tidygraph::filter(!is.na(office))

shop_pol <- osm_pol[, c("shop", "geometry")] %>% tidygraph::filter(!is.na(shop))

sport_pol <- osm_pol[, c("sport", "geometry")] %>% tidygraph::filter(!is.na(sport))

tourism_pol <- osm_pol[, c("tourism", "geometry")] %>% tidygraph::filter(!is.na(tourism))

#ACTIVITY POINTS:   PROCESSING DATA
amenity_pt_list <- c("yes", "bar", "biergarten", "cafe", "fast_food", "ice_cream", "pub", "restaurant", "college", "driving_school", "kindergarten", "language_school", "library", "toy_library", "music_school", "school", "university", "atm", "bank", "bureau_de_change", "clinic", "dentist", "doctors", "hospital", "nursing_home", "pharmacy", "social_facility", "veterinary", "arts_centre", "brothel", "casino", "cinema", "comunity_centre", "conference_centre", "events_venue", "gambling", "love_hotel", "nightclub", "social_centre", "stripclubs", "theatre", "couthouse", "embassy", "townhall", "fire_station", "police", "post_office", "prison", "crematorium", "funeral_hall", "grave_yard", "internet_cafe", "marketplace", "place_of_worship", "public_bath")

amenity_pt <- osm_points[, c("amenity", "geometry")] %>% tidygraph::filter(!is.na(amenity))
amenity_pt_filt <- amenity_pt[amenity_pt$amenity %in% amenity_pt_list, ]
amenity_pol_filt <- amenity_pol[amenity_pol$amenity %in% amenity_pt_list, ]
amenity_pol_filt <- sf::st_centroid(amenity_pol_filt)

amenity_pt_total <- as.data.frame(rbind.fill(amenity_pt_filt, amenity_pol_filt))
amenity_pt_total$category <- "amenity"

if (!is.null(amenity_pt_total)) {
  
  amenity_sustenance <- c("bar", "biergarten", "cafe", "fast_food", "ice_cream", "pub", "restaurant")
  amenity_education <- c("college", "driving_school", "kindergarten", "language_school", "library", "toy_library", "music_school", "school", "university")
  amenity_finance <- c("atm", "bank", "bureau_de_change")
  amenity_healthcare <- c("clinic", "dentist", "doctors", "hospital", "nursing_home", "pharmacy", "social_facility", "veterinary")
  amenity_culture <- c("arts_centre", "brothel", "casino", "cinema", "comunity_centre", "conference_centre", "events_venue", "gambling", "love_hotel", "nightclub", "social_centre", "stripclubs", "theatre")
  amenity_public <- c("couthouse", "embassy", "townhall", "fire_station", "police", "post_office", "prison")
  amenity_others <- c("yes", "crematorium", "funeral_hall", "grave_yard", "internet_cafe", "marketplace", "place_of_worship", "public_bath")
  
  for (i in seq_along(amenity_pt_total$amenity)) {
    if (amenity_pt_total$amenity[i] %in% amenity_sustenance) {
      amenity_pt_total$category[i] <- "sustenance"
    } else if (amenity_pt_total$amenity[i] %in% amenity_education) {
      amenity_pt_total$category[i] <- "education"} 
    else if (amenity_pt_total$amenity[i] %in% amenity_finance) {
      amenity_pt_total$category[i] <- "finance"}
    else if (amenity_pt_total$amenity[i] %in% amenity_healthcare) {
      amenity_pt_total$category[i] <- "healthcare"}
    else if (amenity_pt_total$amenity[i] %in% amenity_public) {
      amenity_pt_total$category[i] <- "public"}
    else if (amenity_pt_total$amenity[i] %in% amenity_others) {
      amenity_pt_total$category[i] <- "others"}
  }
  
}

historic_pt_list <- c("church", "monument", "monastery")

historic_pt <- pts[, c("historic", "geometry")] %>% tidygraph::filter(!is.na(historic))

historic_pt_filt <- historic_pt[historic_pt$historic %in% historic_pt_list, ]
historic_pol_filt <- historic_pol[historic_pol$historic %in% historic_pt_list, ]

historic_pol_filt <- st_centroid(historic_pol_filt)

historic_pt_total <- as.data.frame(rbind.fill(historic_pt_filt, historic_pol_filt))
historic_pt_total$category <- "historic"

leisure_pt_list <- c("yes", "adult_gaming_centre", "amusement_arcade", "beach_resort", "disc_golf_course", "escape_game", "fitness_centre", "fitness_station", "horse_riding", "ice_rink", "marina", "miniature_golf", "pitch", "playground", "sports_centre", "stadium", "swimming_area", "swimming_pool", "track", "water_park")
leisure_pt <- pts[, c("leisure", "geometry")] %>% tidygraph::filter(!is.na(leisure))
leisure_pt_filt <- leisure_pt[leisure_pt$leisure %in% leisure_pt_list, ]
leisure_pol <- osm_pol[, c("leisure", "geometry")] %>% tidygraph::filter(!is.na(leisure))
leisure_pol_filt <- leisure_pol[leisure_pol$leisure %in% leisure_pt_list, ]

leisure_pol_filt <- st_centroid(leisure_pol_filt)

leisure_pt_total <- as.data.frame(rbind.fill(leisure_pt_filt, leisure_pol_filt))
leisure_pt_total$category <- "leisure"

craft_pt_list <- c("yes", "atelier", "bakery", "basket_maker", "blackmith", "bookbinder", "brewery", "cabinet_maker", "builder", "car_painter", "carpenter", "caterer", "clockmaker", "confectionery", "dental_technician", "distillery", "dressmaker", "electronic_repair", "electrician", "engraver", "goldsmith", "handcraft", "jeweller", "key_cutter", "locksmith", "musical_instrument", "optician", "pottery", "printer", "printmaker", "sculptor", "shoemaker", "stonemason", "watchmaker", "winery")

craft_pt <- pts[, c("craft", "geometry")] %>% tidygraph::filter(!is.na(craft))
craft_pt_total <- craft_pt[craft_pt$craft %in% craft_pt_list, ]
craft_pt_total$category <- "craft"

office_pt <- pts[, c("office", "geometry")] %>% tidygraph::filter(!is.na(office))
office_pol_filt <- st_centroid(office_pol)

office_pt_total <- as.data.frame(rbind.fill(office_pt, office_pol_filt))
office_pt_total$category <- "office"

shop_pt <- pts[, c("shop", "geometry")] %>% tidygraph::filter(!is.na(shop))
shop_pol_filt <- st_centroid(shop_pol)

shop_pt_total <- as.data.frame(rbind.fill(shop_pt, shop_pol_filt))
shop_pt_total$category <- "shop"

sport_pt <- pts[, c("sport", "geometry")] %>% tidygraph::filter(!is.na(sport))
sport_pol_filt <- st_centroid(sport_pol)
  
sport_pt_total <- as.data.frame(rbind.fill(sport_pt, sport_pol_filt))
sport_pt_total$category <- "sport"

tourism_pt <- pts[, c("tourism", "geometry")] %>% tidygraph::filter(!is.na(tourism))
tourism_pol_filt <- st_centroid(tourism_pol)

tourism_pt_total <- as.data.frame(rbind.fill(tourism_pt, tourism_pol_filt))
tourism_pt_total$category <- "tourism"

#ACTIVITY POINTS:   Rename columns and bind data
names(amenity_pt_total)[names(amenity_pt_total) == "amenity"] <- "type"

names(historic_pt_total)[names(historic_pt_total) == "historic"] <- "type"

names(leisure_pt_total)[names(leisure_pt_total) == "leisure"] <- "type"

names(craft_pt_total)[names(craft_pt_total) == "craft"] <- "type"

names(office_pt_total)[names(office_pt_total) == "office"] <- "type"

names(shop_pt_total)[names(shop_pt_total) == "shop"] <- "type"

names(sport_pt_total)[names(sport_pt_total) == "sport"] <- "type"

names(tourism_pt_total)[names(tourism_pt_total) == "tourism"] <- "type"

activity_df <- rbind.fill(amenity_pt_total, historic_pt_total, leisure_pt_total, craft_pt_total, office_pt_total, shop_pt_total, sport_pt_total, tourism_pt_total)

activity_coord <- as.data.frame(st_coordinates(st_as_sf(activity_df)))
activity_df$lng <- activity_coord$X
activity_df$lat <- activity_coord$Y
  
#ACTIVITY POINTS:   Clean SF features and reset index (distm function does not work if it is not a standard dataframe)
activity_df <- as.data.frame(activity_df)
activity_df <- activity_df[, -2]
row.names(activity_df) <- NULL
  
#ACTIVITY POINTS:  define an EPS distance 
eps_b <- 60 # distance in meters to search for other points to form a cluster, the recommended distance is the standard block size in your analysis area
min_pt <- 5 # the minimum number of points to form a cluster, this parameter will be assigned to certain categories of the POI data set as public, or education
max_pt <- 10 # the minimum number for the rest of the categories as shops, sustenance etc
  
#ACTIVITY POINTS:   Assign cluster function

# we are going to create our own function based on DBSCAN base function, in order to be able to assign different minimum points based on different categories

assign_clusters <- function(poi_df, minPts = NA) {
    
  # for certain categories as (eg) Education and Health we consider that more than 5 points is already a cluster 
  # for the rest of the categories we consider a higher number of points - 10
    
  if (is.na(minPts)) {
    if (poi_df[1, "category"] %in% c("public", "leisure", "education", "healthcare", "finance", "others", "historic", "tourism")) {
      minPts <- min_pt
    } else minPts <- max_pt
  }
    
  eps <- eps_b #  in meters
    
  poi_df[c("lng", "lat")] %>% 
    distm(fun = distHaversine) %>%
    as.dist() %>% 
    dbscan(eps = eps, minPts = minPts) %>% 
    .[["cluster"]] %>% 
    cbind(poi_df, cluster = .)
    
}
  
#ACTIVITY POINTS:   apply functions to data and generate clusters
Clean_data <- activity_df %>%
  split(activity_df$category) %>%  # splits, separate the data frame by categories based on the values of the column category
  purrr::map_df(assign_clusters) # here we are using the map_df() function from purrr package (Apply a function to each element of a list or atomic vector) in order to apply our custom function to Activity_df

get_hull <- function(df) { # another custom function, this one generate a polygon that contains all points in a given cluster
  
  cbind(df$lng, df$lat) %>% 
    as.matrix() %>%
    st_multipoint() %>% 
    st_convex_hull() %>% 
    st_sfc(crs = 4326) %>% 
    {st_sf(category = df$category[1], cluster = df$cluster[1], geom = .)}
}
  
hulls <- function(df) { # custom function which have nested the previous function get_hull() and split the data frame into different cluster and then applies the get_hull() function to each seprate one
  
  df %>%
    split(.$cluster) %>% 
    map(get_hull)
}

hulls_cat <- Clean_data %>%
  dplyr::group_by(category) %>%
  dplyr::summarise() # summarise() creates a new data frame. It will have one (or more) rows for each combination of grouping variables; It will contain one column for each grouping variable and one column for each of the summary statistics that you have specified.

map_cluster_hulls <- Clean_data %>%
  dplyr::filter(cluster != 0) %>%
  dplyr::select(lng, lat, category, cluster) %>% 
  split(.$category) %>% 
  map(hulls)

delete.NULLs  <-  function(x.list){   # delele null/empty entries in a list
  x.list[unlist(lapply(x.list, length) != 0)]
}
  
#ACTIVITY POINTS:   Clean cluster data
mch <- delete.NULLs(map_cluster_hulls)
mch <- reshape2::melt(mch, id = c("category", "cluster", "geom"))
mch <- data.frame(mch$category, mch$cluster, mch$geom)

if (is.null(mch)) {
  
  mch <- NULL
  Cluster <- mch
  
} else if (nrow(mch) == 0) {
  
  mch <- NULL
  Cluster <- mch
  
} else {
  
  Cluster <- mch
  
}

if (is.null(activity_df)) {
  
  activity_df <- NULL
  Activity_Points_df <- activity_df
  
} else if (nrow(activity_df) == 0) {
  
  activity_df <- NULL
  Activity_Points_df <- activity_df
  
} else {
  
  Activity_Points_df <- activity_df
  
}
  
Activity_Points_sf <-  sf::st_as_sf(Activity_Points_df, coords = c("lng", "lat"), crs = 4326, agr = "constant")  
Activity_Points__27700 <- sf::st_transform(Activity_Points_sf, 27700)

Cluster_27700 <- sf::st_transform(sf::st_as_sf(Cluster), 27700)

## Extracting desired OSM Data - Public Transport Routes

df_multilines <- osm_mline

routes <- st_crop(df_multilines, st_bbox(BBox_poly))
routes <- as.data.frame(routes)
routes <- filter(routes, routes$route != "")

routes <- st_transform(st_as_sf(routes), crs = 4326)
route_list <- unique(routes$route)
route_color <- RColorBrewer::brewer.pal(n = length(route_list), name = "Paired")
    
routes_simplified <- routes
    
for (i in seq_along(routes_simplified$route)) {
  a <- which(route_list == routes_simplified$route[i])
  routes_simplified$color[i] <- route_color[a] 
}

Routes <- routes_simplified

Routes_27700 <- data.frame(name = Routes$name, from = Routes$from, to = Routes$to, route = Routes$route, geometry = Routes$geometry)
Routes_27700 <- sf::st_transform(sf::st_as_sf(Routes_27700), 27700)

STOPS <- jsonlite::fromJSON(paste0("https://api.tfl.gov.uk/Stoppoint?lat=", latitude, "&lon=", longitude, "&stoptypes=NaptanMetroStation,NaptanRailStation,NaptanBusCoachStation,NaptanFerryPort,NaptanPublicBusCoachTram&radius=", Analysis_rad))
STOPS_df <- STOPS$stopPoints
STOPS_sf <- sf::st_as_sf(STOPS_df, coords = c("lon", "lat"), crs = 4326)

for (i in seq_along(unique(unlist(STOPS_sf$modes)))) { # create new columns based on unique transport modes within the data frame
  
  STOPS_sf[unique(unlist(STOPS_sf$modes))[i]] <- NA
  
}

for (i in seq_along(STOPS_sf$naptanId)) { # unlisting nested dataframes within STOPS dataframe
  
  for (j in seq_along(STOPS_sf$lineModeGroups[[i]]$modeName)) {

    STOPS_sf[i, STOPS_sf$lineModeGroups[[i]]$modeName[j]] <- length(STOPS_sf$lineModeGroups[[i]]$lineIdentifier[[j]])
    
  }
  
}

STOPS_sf$children <- NULL
STOPS_sf$additionalProperties <- NULL
STOPS_sf$modes <- NULL
STOPS_sf$lines <- NULL
STOPS_sf$lineGroup <- NULL
STOPS_sf$lineModeGroups <- NULL

STOPS_sf <- sf::st_transform(STOPS_sf, 27700)
  
## Census Data

mult <- 1.5
LSOA_QUERY <- data.frame()
  
  while (length(LSOA_QUERY) == 0) {
    mult <- mult + 0.5
    BBox_2 <- center_bbox(longitude, latitude, Analysis_rad*mult, Analysis_rad*mult)
    
    bbox_coord <- as.data.frame(as.list(BBox_2))
    ESRI_query <- paste("https://ons-inspire.esriuk.com/arcgis/rest/services/Census_Boundaries/Lower_Super_Output_Areas_December_2011_Centroids/MapServer/0/query?where=1%3D1&outFields=*&geometry=", as.character(bbox_coord$left), "%2C", as.character(bbox_coord$bottom), "%2C", as.character(bbox_coord$right), '%2C', as.character(bbox_coord$top), "&geometryType=esriGeometryEnvelope&inSR=4326&spatialRel=esriSpatialRelIntersects&outSR=4326&f=json", sep = "")
    
    query_result <- jsonify::from_json(ESRI_query)
    
    LSOA_CODES <- query_result$features$attributes$lsoa11cd
    LSOA_LAT <- query_result$features$geometry$y
    LSOA_LNG <- query_result$features$geometry$x
    
    LSOA_QUERY <- data.frame(lsoacd = LSOA_CODES, lat = LSOA_LAT, lng = LSOA_LNG)
  }
  
  LSOA_df <- LSOA_QUERY
  LSOA_boundaries <- data.frame(geometry = NA, lsoacd = NA)
  
  for (i in seq_along(LSOA_CODES)) {
    q <- paste("https://services1.arcgis.com/ESMARspQHYMw9BZ9/arcgis/rest/services/Lower_Layer_Super_Output_Areas_December_2011_Boundaries_EW_BFE_V2/FeatureServer/0/query?where=LSOA11CD%20%3D%20'", as.character(LSOA_CODES[i]), "'&outFields=*&outSR=4326&f=json", sep = "")
    qr <- jsonify::from_json(q)
    xym <- cbind(qr[["features"]][["geometry"]][["rings"]][[1]][[1]][,1], qr[["features"]][["geometry"]][["rings"]][[1]][[1]][, 2])
    p = sp::Polygon(xym)
    ps = sp::Polygons(list(p),1)
    sps = sp::SpatialPolygons(list(ps))
    sps <- sf::st_as_sf(sps)
    
    LSOA_boundaries[i, "geometry"] <- sps
    LSOA_boundaries[i, "lsoacd"] <- LSOA_CODES[i]
  }
  
  LSOA_boundaries <- st_as_sf(LSOA_boundaries)
  LSOA_sf <- LSOA_boundaries
  
  age_list <- list()
  pop_list <- list()
  resi_list <- list()
  health_list <- list()
  tenure_list <- list()
  household_list <- list()
  travel_list <- list()
  distance_list <- list()
  
  CENSUS_df <- list()
  
  for (i in seq_along(LSOA_CODES)) {
    age <- nomisr::nomis_get_data(id = "NM_159_1", time = "latest", geography = LSOA_CODES[i])
    age_list[[i]] <- age
    resi <- nomisr::nomis_get_data(id = "NM_501_1", time = "latest", geography = LSOA_CODES[i])
    resi_list[[i]] <- resi
    health <- nomisr::nomis_get_data(id = "NM_531_1", time = "latest", geography = LSOA_CODES[i])
    health_list[[i]] <- health
    tenure <- nomisr::nomis_get_data(id = "NM_537_1", time = "latest", geography = LSOA_CODES[i])
    tenure_list[[i]] <- tenure
    household <- nomisr::nomis_get_data(id = "NM_538_1", time = "latest", geography = LSOA_CODES[i])
    household_list[[i]] <- household
    travel <- nomisr::nomis_get_data(id = "NM_568_1", time = "latest", geography = LSOA_CODES[i])
    travel_list[[i]] <- travel
    distance <- nomisr::nomis_get_data(id = "NM_153_1", time = "latest", geography = LSOA_CODES[i])
    distance_list[[i]] <- distance
  }
  
  age_structure <- do.call(rbind, age_list)
  age_structure <- age_structure[, c("GEOGRAPHY_CODE", "CELL_NAME", "CELL_CODE", "CELL_TYPE", "MEASURES_NAME", "OBS_VALUE")]
  CENSUS_df$age_structure <- unique(age_structure)
  residence_type <- do.call(rbind, resi_list)
  residence_type <- residence_type[, c("GEOGRAPHY_CODE", "CELL_NAME", "CELL_CODE", "CELL_TYPE", "MEASURES_NAME", "OBS_VALUE")]
  CENSUS_df$residence_type <- unique(residence_type)
  general_health <- do.call(rbind, health_list)
  general_health <- general_health[, c("GEOGRAPHY_CODE", "C_HEALTH_NAME", "C_HEALTH_CODE", "C_HEALTH_TYPE", "MEASURES_NAME", "OBS_VALUE")]
  CENSUS_df$general_health <- unique(general_health)
  tenure_households <- do.call(rbind, tenure_list)
  tenure_households <- tenure_households[, c("GEOGRAPHY_CODE", "C_TENHUK11_NAME", "C_TENHUK11_CODE", "C_TENHUK11_TYPE", "MEASURES_NAME", "OBS_VALUE")]
  CENSUS_df$tenure_households <- unique(tenure_households)
  household_size <- do.call(rbind, household_list)
  household_size <- household_size[, c("GEOGRAPHY_CODE", "CELL_NAME", "CELL_CODE", "CELL_TYPE", "MEASURES_NAME", "OBS_VALUE")]
  CENSUS_df$household_size <- unique(household_size)
  method_of_travel_to_work <- do.call(rbind, travel_list)
  method_of_travel_to_work <- method_of_travel_to_work[, c("GEOGRAPHY_CODE", "CELL_NAME", "CELL_CODE", "CELL_TYPE", "MEASURES_NAME", "OBS_VALUE")]
  CENSUS_df$method_of_travel_to_work <- unique(method_of_travel_to_work)
  distance_traveled_to_work <- do.call(rbind, distance_list)
  distance_traveled_to_work <- distance_traveled_to_work[, c("GEOGRAPHY_CODE", "CELL_NAME", "CELL_CODE", "CELL_TYPE", "MEASURES_NAME", "OBS_VALUE")]
  CENSUS_df$distance_traveled_to_work <- unique(distance_traveled_to_work)
  
  LSOA_boundaries_27700 <- sf::st_transform(sf::`st_crs<-`(LSOA_boundaries, 4326), 27700)
  
#6.- Saving the Data
  
Save_folder <- "C://Users//Jorge//Documents//Grimshaw//DT//ACADIA 2022//Results//"

DSM_df <- raster::as.data.frame(DSM, xy = TRUE)
DTM_df <- raster::as.data.frame(DTM, xy = TRUE)

DTM_sf <- sf::st_as_sf(DTM_df ,coords = c("x","y"))
DSM_sf <- sf::st_as_sf(DSM_df ,coords = c("x","y"))

sf::st_write(OSM_Build_27700, paste0(Save_folder, "Buildings.shp"), delete_dsn = TRUE)
sf::st_write(St_network_27700, paste0(Save_folder, "St_Network.shp"), delete_dsn = TRUE)
sf::st_write(St_network_Feat_27700, paste0(Save_folder, "St_Network_Features.shp"), delete_dsn = TRUE)
sf::st_write(natural_27700, paste0(Save_folder, "Natural.shp"), delete_dsn = TRUE)
sf::st_write(Activity_Points__27700, paste0(Save_folder, "POI.shp"), delete_dsn = TRUE)
sf::st_write(Cluster_27700, paste0(Save_folder, "Cluster_POI.shp"), delete_dsn = TRUE)
sf::st_write(LSOA_boundaries_27700, paste0(Save_folder, "LSOA.shp"), delete_dsn = TRUE)
sf::st_write(Routes_27700, paste0(Save_folder, "PT_Routes.shp"), delete_dsn = TRUE)
sf::st_write(STOPS_sf, paste0(Save_folder, "PT_STOPS.shp"), delete_dsn = TRUE)
sf::write_sf(DTM_sf, paste0(Save_folder,"DTM_sf.shp"), delete_dsn = TRUE)
sf::write_sf(DSM_sf, paste0(Save_folder,"DSM_sf.shp"), delete_dsn = TRUE)

write.csv(CENSUS_df$age_structure, file = paste0(Save_folder, "CENSUS_age_structure.csv"), sep = ",")
write.csv(CENSUS_df$general_health, file = paste0(Save_folder, "CENSUS_general_health.csv"), sep = ",")
write.csv(CENSUS_df$residence_type, file = paste0(Save_folder, "CENSUS_residence_type.csv"), sep = ",")
write.csv(CENSUS_df$household_size, file = paste0(Save_folder, "CENSUS_household_size.csv"), sep = ",")
write.csv(CENSUS_df$method_of_travel_to_work, file = paste0(Save_folder, "CENSUS_method_of_travel_to_work.csv"), sep = ",")
write.csv(CENSUS_df$distance_traveled_to_work, file = paste0(Save_folder, "CENSUS_distance_traveled_to_work.csv"), sep = ",")