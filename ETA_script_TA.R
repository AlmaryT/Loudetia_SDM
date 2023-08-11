# Load in libraries
library(terra)
library(rgbif)
library(DT)
library(tidyverse)
library(CoordinateCleaner)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(mapview)
library(flexsdm)
library(here)
library(corrplot)
library(SDMtune)
library(plotROC)
library(rasterVis)
library(virtualspecies)
library(blockCV) # spatial block for cross-validation
library(patchwork)
library(ggsci)
library(zeallot)
library(rnaturalearthhires)

#### 1. Download locality data ----
# Name your species
spp <- c("Loudetia simplex")

#### 1.1 Download locality data from GBIF ----

# download GBIF occurrence data for this species; this takes time if there are many data points!
gbif_download <- occ_data(scientificName = spp, hasCoordinate = TRUE, limit = 10000)
# select only the data (ignore the other metadata for now)
gbif_data <- gbif_download$data
# provide a unique ID for this dataset
gbif_data %>% 
  mutate(download_ID = paste0('gbif_', row_number())) %>%
  
  dplyr::select(download_ID, species, decimallongitude = decimalLongitude, decimallatitude = decimalLatitude, year) -> gbif_sel

#### 1.2 Download localilty data from BIEN ----
bien_sel <- read_csv("bien_selected_vars.csv")
# bien_sel <- read_csv("C:/MASTER TCHANA/Loudetia_SDM/bien_selected_vars.csv")

#### Tchana field data ----
field_sel <- read_csv("Field_data.csv")
# field_sel <- read_csv("C:/MASTER TCHANA/Loudetia_SDM/Field_data.csv")
# field_sel


#### 1.3 Combine datasets ----
dat <- rbind(gbif_sel, bien_sel, field_sel)

### 2. Clean GBIF data ----
# using custom filter
dat %>%
  filter(!decimallatitude == 0 | !decimallongitude == 0) %>%
  distinct(decimallongitude,decimallatitude,download_ID,species, year, .keep_all = TRUE) -> gbif_filt

# using cc_*() functions
gbif_filt %>%
  cc_cen(lat = 'decimallatitude', lon = 'decimallongitude', buffer = 2000) %>% # remove country centroids within 2km 
  cc_cap(lat = 'decimallatitude', lon = 'decimallongitude', buffer = 2000) %>% # remove capitals centroids within 2km
  cc_inst(lat = 'decimallatitude', lon = 'decimallongitude', buffer = 2000) -> spp_clean  # remove zoo and herbaria within 2km 

# Downloading country borders
#### Convert the points to a spatial vector
dat_vect <- vect(dat, geom = c('decimallongitude', 'decimallatitude'), crs = '+proj=longlat +datum=WGS84')

# Downloading country borders
# using the rnaturalearth package, download the country border and return it as an sf object
east_afr <- ne_countries(country = c('United Republic of Tanzania', 'Burundi', 'Malawi', 'Mozambique', 'Zambia', 'Zimbabwe', 'Rwanda', 'Uganda', 'Kenya'), scale = 'medium', returnclass = 'sf') %>% dplyr::select(name)

# Download the south east provinces of DRC
drc_states <- ne_states(country = c('Democratic Republic of the Congo'), returnclass = 'sf')
drc_states$name_en

drc_sel_states <- drc_states %>% filter(name_en %in% c('Katanga Province', 'South Kivu', 'Maniema')) %>% dplyr::select(name)
plot(vect(drc_sel_states))
# Note, in 2015 DRC changed the names of its states and these names are no longer used... https://en.wikipedia.org/wiki/Provinces_of_the_Democratic_Republic_of_the_Congo

# Combine the two datasets
east_afr_combine <- rbind(east_afr, drc_sel_states)
plot(vect(east_afr_combine))

#### Mask the points to the shape of the selected countries/provinces
dat_east_afr <- mask(dat_vect, vect(east_afr_combine))

#### Convert back to a dataframe 
coords <- terra::geom(dat_east_afr) %>% as.data.frame() %>% dplyr::select(x, y)
dat_east_afr_df <- cbind(as.data.frame(dat_east_afr), coords) %>% rename(decimallongitude = x, decimallatitude = y)
dat_east_afr_df


# Clean the dataset by selecting only the columns we need, rename the lon/lat variables
dat_east_afr_df %>% dplyr::select(download_ID, species, lon = decimallongitude, lat = decimallatitude, year) -> spp_sel

# Now we want to convert & project our data to the correct Coordinate Reference System
spp_sf <- st_as_sf(spp_sel, coords = c('lon', 'lat'))

#### NEW: filter out points in central S.Afr
spp_sf %>% filter(!download_ID %in% c('bien_2148', 'bien_2131', 'bien_1971', 'bien_2116')) -> spp_sf

spp_sf %>% st_drop_geometry() -> spp_sel2
spp_sel2 <- cbind(spp_sel2, st_coordinates(spp_sf)) %>% rename(lon = X, lat = Y)

### 3. Spatial thinning ----
# load in worldclim locally
worldclim <- rast("all_env_vars.tif")
# worldclim <- rast("C:/MASTER TCHANA/Loudetia_SDM/data/rast/subsah_afr_cov.tiff")

# Crop and mask
ETA <- east_afr_combine
wc_ETA <- worldclim %>% crop(., vect(ETA)) %>% mask(., vect(ETA))

# Thin records by distance
set.seed(42)
spp_filt <- occfilt_geo(
  data = spp_sel2 %>% rename(x = lon, y = lat),
  x = 'x',
  y = 'y',
  method = c('defined', d = 1), # set a 1 km threshold
  env_layer = wc_ETA,
  prj = crs(wc_ETA)
)

# join the result of the filter back onto the metadata and convert to sf
spp_filt_sf <- spp_filt %>% dplyr::select(-c(x, y)) %>% left_join(spp_sel, by = 'download_ID') %>% st_as_sf(coords = c('lon','lat'), crs = st_crs(4326))

### 4. Create pseudo-absences ----
set.seed(42)
pa <-  sample_pseudoabs(
  data = spp_filt,
  x = 'x',
  y = 'y',
  n = nrow(spp_filt),
  method = c('geo_const', width = 10000),
  rlayer = wc_ETA,
  calibarea = vect(ETA)
)
glimpse(pa)

# process our presence points to match pseudoabsences
pres <- spp_filt %>% dplyr::select(x, y) %>% mutate(pr_ab = 1)
glimpse(pres)

# combine presences and pseudoabsences
all_pts <- rbind(pres, pa)
# convert to spatial format
all_pts_sf <- all_pts %>% mutate(pr_ab = as.factor(pr_ab)) %>% st_as_sf(coords = c('x','y'), crs = st_crs(4326))

### 5.  Process environmental data ----
# select variables by hand 
sel_vars <- c('bio01', 'bio02', 'bio04', 'bio12', 'bio17', 'bio15', 'bio19', 'tri', 'ph', 'sand', 'silt')

vars_sel_fin <- wc_ETA[[sel_vars]]
names(vars_sel_fin) <- c('mean_ann_t', 'mean_diurnal_t_range', 't_seas', 'ann_p', 'p_dry_q', 'p_seas', 'p_cold_q', 'tri', 'ph', 'sand', 'silt')

# Re-scale temperature values
vars_sel_fin[[c(1:2)]] <- vars_sel_fin[[c(1:2)]]/10
vars_sel_fin[[c(3)]] <- vars_sel_fin[[c(3)]]/100

#### Check projections 
crs(vars_sel_fin) == crs(vect(all_pts_sf))
all_pts_vect <- project(vect(all_pts_sf), vars_sel_fin)

### 6.  Prepare extracted data ----

#### Extract data
# Prepare a SWD (Sample with Data), which is a class of data specifically used in the SDMtune package
all_pts_df <- as.data.frame(all_pts, geom = 'XY')
all_pts_df

SWDdata <- prepareSWD(
  species = 'LS',
  p = all_pts_df %>% filter(pr_ab == 1) %>% dplyr::select(.,x, y),
  a = all_pts_df %>% filter(pr_ab == 0) %>% dplyr::select(.,x, y),
  env = vars_sel_fin
)

# For some reason ID is included, so we need to remove this...
SWDdata@data <- SWDdata@data[-1]

#### 7.1 cross-validation and hyper paremter tuning ----
#### Spatial CV

# Environmental clustering of the training points into 5 distinct groups
# these groups show clusters with similar environmental conditions based on the environmental layers provided
k = 5
set.seed(42)
env_blocks <- envBlock(rasterLayer = raster::raster(vars_sel_fin),
                       speciesData = all_pts_sf,
                       species = "pr_ab",
                       k = k, 
                       standardization = 'standard', 
                       rasterBlock = FALSE)

env_blocks$foldID
all_pts_sf$fold <- env_blocks$foldID

# visualise the different training dataset clusters
ggplot() + 
  geom_sf(data = ETA, fill = 'gray90') +
  geom_sf(data = all_pts_sf, aes(col = as.factor(fold), shape = as.factor(pr_ab))) +
  scale_color_npg(name = 'Fold') +
  scale_shape(name = 'Presence/Absence') +
  theme_void()

#### Run a random forest model
# when we run the model, we now specify the different environmental blocks
# this means that the model is run 5 times. 
set.seed(42)
cv_model <- train(method = 'RF', data = SWDdata, folds = env_blocks)
cv_model

#### Hyper parameter test
# this sets up a combination of different parameters for us to test. We then look to see which model produces the best AUC values
h <- list(mtry = seq(floor(sqrt(7)), 10, 2),
          ntree = seq(250, 1500, 250),
          nodesize = seq(1, 10, 2))

# we cannot check them all, so we take a look at 50 possible combinations
exp_1 <- randomSearch(cv_model, 
                      hypers = h, 
                      metric = 'auc',
                      pop = 50,
                      seed = 42)

# these are the results of the 50 randomly selected models
exp_1@results

# there isn't a huge difference between them, but we will still select the best possible combination of parameters and train our final model using these parameters

### 7.2  Run & evaluate model ----
#### Final model 
set.seed(42)
final_model <- train('RF',
                     data = exp_1@models[[1]]@data,
                     mtry = exp_1@results[1,1],
                     ntree = exp_1@results[1,2],
                     nodesize = exp_1@results[1,3],
                     folds = env_blocks)

# True Skill Statistic
tss(final_model, test = TRUE)

# Area Under Curve
auc(final_model, test = TRUE)

# Save model
saveRDS(final_model, 'output/models/ETA_RF_model.RDS')

### 8.  Variable insights & predicting habitat suitability ----

#### Predict habitat suitability
vars_sel_fin_stack <- raster::stack(vars_sel_fin)
# vars_sel_fin_stack <- raster::stack('data/rast/subsah_afr_cov_11.tiff')
# names(vars_sel_fin_stack)

pred <- predict(final_model, data = vars_sel_fin_stack)
raster::writeRaster(pred, 'data/rast/ETA_Loudetia_simplex_pred_prob.tiff', overwrite = TRUE)
# pred <- raster('data/rast/AFR_Loudetia_simplex_pred_prob.tiff')

#### END ####

