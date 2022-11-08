
# source("~/00-mount.R")
source("scripts/00_setup.R")
# source("scripts/functions.R")

seq(as_date("1979-01-01"),
    as_date("2021-12-31"),
    "1 month") -> dates

# remove a year for lag + train window (12)
# dates[month(dates) %in% c(3,6,9,12) & 
#         year(dates) > 1980] -> dates_sub

dates %>% 
  days_in_month() %>% 
  unname() -> days_mth



# LAND MASK -----------------------------------------------------------------------------------


c(st_point(c(-180.125, -90.125)),
  st_point(c(179.875, 90.125))) %>% 
  st_bbox() %>% 
  st_set_crs(4326) %>% 
  st_as_stars(dx = 0.05, dy = 0.05, values = -9999) -> rast_reference_0.05

c(st_point(c(-180.125, -90.125)),
  st_point(c(179.875, 90.125))) %>% 
  st_bbox() %>%
  st_set_crs(4326) %>% 
  st_as_stars(dx = 0.25, dy = 0.25, values = -9999) -> rast_reference_0.25

"~/bucket_mine/misc_data/ne_50m_land/ne_50m_land.shp" %>% 
  st_read() %>% 
  mutate(a = 1) %>% 
  select(a) %>% 
  st_rasterize(rast_reference_0.05) -> land

land %>% 
  st_warp(rast_reference_0.25, use_gdal = T, method = "mode") %>% 
  setNames("a") %>% 
  mutate(a = ifelse(a == -9999, "ocean", "land") %>% factor()) -> land

land %>% 
  st_set_dimensions(c(1,2), names = c("longitude", "latitude")) -> land
  
st_get_dimension_values(rast_reference_0.25, "x") -> vect_lon
st_get_dimension_values(rast_reference_0.25, "y") -> vect_lat



# RESPONSE VAR --------------------------------------------------------------------------------

"~/bucket_mine/era/monthly/era5_monthly_mean_daily_precip_1979-2021_shifted.nc" %>% 
  read_ncdf() %>% 
  
  mutate(tp = set_units(tp, mm)) %>%
  
  st_apply(c(1,2), function(x){
    x * days_mth
  },
  .fname = "time") %>%
  aperm(c(2,3,1)) %>%
  st_set_dimensions("time", values = dates) -> precip


# *****************

# each date in "by" is the left bound (included)
# the next date is the right bound (excluded)
precip %>% 
  aggregate(by = dates[month(dates) %in% c(3,6,9,12)], 
            sum) -> precip_seas

# remove last time step (Dec; incomplete quarter)
precip_seas %>% 
  slice(time, -(dim(precip_seas)[1])) %>% 
  aperm(c(2,3,1)) -> precip_seas


# response dataset

#points:

# expand_grid(pos_lon = seq(40, 1420, 30),
#             pos_lat = seq(40, 701, 30)) %>%
#   
#   mutate(longitude = vect_lon[pos_lon],
#          latitude = vect_lat[pos_lat]) %>%
#   
#   left_join(as_tibble(land), by = c("longitude", "latitude")) %>% 
#   filter(a == "land") %>%
#   select(-a) %>% 
#   
#   filter(abs(latitude) < 67) -> tb_resp_locs
#
# ggplot() +
#   geom_stars(data = land) +
#   geom_point(data = tb_resp_locs, aes(x = longitude, y = latitude)) +
#   coord_equal()
# 
# left_join(tb_resp_locs, as_tibble(precip_seas), by = c("longitude", "latitude")) -> resp_precip


# western us:

c(st_point(c(-126.625, 32.625)),
  st_point(c(-115.125, 48.875))) %>% 
  
  st_bbox() %>% 
  st_set_crs(4326) -> w_us

precip_seas %>% 
  st_crop(w_us) -> s_resp

s_resp[st_crop(land, s_resp) == "ocean"] <- NA

s_resp %>% 
  as_tibble() -> tb_resp

expand_grid(longitude = seq(min(tb_resp$longitude), max(tb_resp$longitude), 0.25*4),
            latitude = seq(min(tb_resp$latitude), max(tb_resp$latitude), 0.25*4)) %>% 
  left_join(tb_resp, by = c("longitude", "latitude")) -> tb_resp


# resp_precip %>% 
#   slice(time, 1) %>% 
#   as_tibble() %>%
#   filter(!is.na(tp)) %>% 
#   select(-tp) -> tb_resp_locs


  

# PREDICTORS: LOCAL -------------------------------------------------------------------------------

lg <- 6 # lag
train_window <- 6 # training window width


## PRECIP (MEMORY) ----

# tb_resp_locs %>% 
#   
#   mutate(r = row_number(),
#          r = str_pad(r, 3, "left", "0")) %>%
#   
#   pmap_dfr(function(longitude, latitude, r, ...){
#     
#     print(str_glue("{r} / {nrow(tb_resp_locs)}"))
#     
#     pos_lon <- which.min(abs(vect_lon - longitude))
#     pos_lat <- which.min(abs(vect_lat - latitude))
#     
#     precip[,
#            (pos_lon-1):(pos_lon+1), 
#            (pos_lat-1):(pos_lat+1),
#     ] %>%
#       st_apply(3, mean, na.rm = T) %>% 
#       pull(1) -> x_ts
#     
#     if(!anyNA(x_ts)){
#       
#       dates_sub %>% 
#         map_dfr(function(d){
#           
#           w <- x_ts[(which(dates == d) - lg - train_window):(which(dates == d) - lg - 1)]
#           
#           tibble(
#             
#             pos_lon = pos_lon,
#             pos_lat = pos_lat,
#             tp_avg6mths = mean(w),
#             time = d
#             
#           )
#         })
#     }
#   }) -> pred_precip

# apply SPATIAL moving window (3x3 mean)
foc <- function(x, w) {
  raster::as.matrix(raster::focal(raster::raster(x), w))
}

precip %>% 
  
  # extra rows+cols for 3x3 windows to be complete
  filter(longitude >= (w_us[1]-0.2),
         longitude <= (w_us[3]+0.2),
         latitude >= (w_us[2]-0.2),
         latitude <= (w_us[4]+0.2)) %>%
  
  st_apply(3, foc, w = matrix(1/9, 3, 3)) -> precip_mw

precip_mw[st_crop(land, precip_mw) == "ocean"] <- NA

precip_mw %>% 
  st_crop(w_us) -> precip_mw
  
# apply TEMPORAL moving window
precip_mw %>% 
  st_apply(c(1,2), function(x_ts){
    
    if(!anyNA(x_ts)){
      
      zoo::rollmean(x_ts,
                    k = train_window,
                    na.pad = T,
                    align = "right")
      
    } else {
      rep(NA, dim(precip_mw)[3])
    }
  },
  FUTURE = T,
  .fname = "time") %>% 
  aperm(c(2,3,1)) -> precip_mw_mw

# time = add lag
precip_mw_mw %>% 
  st_set_dimensions("time", values = dates + months(lg+1)) -> precip_mw_mw

# as table
precip_mw_mw %>% 
  filter(time <= st_get_dimension_values(s_resp, "time") %>% last()) %>% 
  as_tibble() %>% 
  filter(!is.na(tp)) -> pred_precip

# join with resp_precip by pos_lon, pos_lat, and time


## SURF TEMP ----

"~/bucket_mine/era/monthly/era5_monthly_mean_meantemp_shifted.nc" %>% 
  read_ncdf() -> pred_temp

# pred_temp %>% 
#   st_apply(3, foc, w = matrix(1,3,3),
#            FUTURE = T) -> pred_temp
# 
# tb_resp_locs %>% 
#   
#   mutate(r = row_number(),
#          r = str_pad(r, 3, "left", "0")) %>%
#   
#   # future_pmap_dfr(function(pos_lon, pos_lat, r, ...){
#   pmap_dfr(function(pos_lon, pos_lat, r, ...){
#     
#     print(str_glue("{r} / {nrow(tb_resp_locs)}"))
#     
#     pred_temp[ , pos_lon, pos_lat, ] %>%
#       pull(1) -> x_ts
#     
#     if(!anyNA(x_ts)){
#       
#       dates_sub %>% 
#         map_dfr(function(d){
#           
#           w <- x_ts[(which(dates == d) - lg - train_window):(which(dates == d) - lg - 1)]
#           
#           tibble(
#             
#             pos_lon = pos_lon,
#             pos_lat = pos_lat,
#             tas_avg6mths = mean(w),
#             time = d
#             
#           )
#         })
#     }
#   }) -> pred_temp

pred_temp %>% 
  
  # extra rows+cols for 3x3 windows to be complete
  filter(longitude >= (w_us[1]-0.2),
         longitude <= (w_us[3]+0.2),
         latitude >= (w_us[2]-0.2),
         latitude <= (w_us[4]+0.2)) %>%
  
  st_apply(3, foc, w = matrix(1/9, 3, 3)) -> temp_mw

precip_mw[st_crop(land, precip_mw) == "ocean"] <- NA

precip_mw %>% 
  st_crop(w_us) -> precip_mw

# apply TEMPORAL moving window
precip_mw %>% 
  st_apply(c(1,2), function(x_ts){
    
    if(!anyNA(x_ts)){
      
      zoo::rollmean(x_ts,
                    k = train_window,
                    na.pad = T,
                    align = "right")
      
    } else {
      rep(NA, dim(precip_mw)[3])
    }
  },
  FUTURE = T,
  .fname = "time") %>% 
  aperm(c(2,3,1)) -> precip_mw_mw

# time = add lag
precip_mw_mw %>% 
  st_set_dimensions("time", values = dates + months(lg+1)) -> precip_mw_mw

# as table
precip_mw_mw %>% 
  filter(time <= st_get_dimension_values(s_resp, "time") %>% last()) %>% 
  as_tibble() %>% 
  filter(!is.na(tp)) -> pred_precip


  
  st_bbox() %>% 
  st_set_crs(4326) %>% 
  
  {st_crop(pred_temp, .)} %>% 
  
  st_apply(3, foc, w = matrix(1/9, 3, 3)) -> temp_agg

temp_agg[st_crop(land, temp_agg) == "ocean"] <- NA

temp_agg %>% 
  st_crop(w_us) -> temp_agg

#aggregate temporally (6 month running mean)
temp_agg %>% 
  st_apply(c(1,2), function(x_ts){
    
    if(!anyNA(x_ts)){
      
      zoo::rollmean(x_ts,
                    k = train_window,
                    na.pad = T,
                    align = "right")
      
    } else {
      rep(NA, dim(temp_agg)[3])
    }
  },
  FUTURE = T,
  .fname = "time") -> temp_agg_agg

# time = add lag
temp_agg_agg %>% 
  aperm(c(2,3,1)) %>% 
  st_set_dimensions("time", values = dates + months(lg+1)) -> temp_agg_agg

# as table
temp_agg_agg %>% 
  filter(time %in% dates_sub) %>% 
  as_tibble() %>% 
  rowwise() %>% 
  mutate(pos_lon = which.min(abs(vect_lon - longitude)),
         pos_lat = which.min(abs(vect_lat - latitude))) %>% 
  select(pos_lon, pos_lat, t2m, time) %>% 
  ungroup() %>% 
  filter(!is.na(t2m)) -> pred_temp

# join with resp_precip by pos_lon, pos_lat, and time


## ELEVATION ---- 

"~/bucket_mine/era/monthly/era5_geopotential_shifted.nc" %>% 
  read_ncdf() %>% 
  adrop() %>% 
  mutate(z = set_units(z, NULL),
         z = z/9.80665,
         z = ifelse(z < 0, 0 , z)) -> pred_elev

st_dimensions(land) <- st_dimensions(pred_elev)

c(pred_elev, land) %>% 
  mutate(z = ifelse(a == "ocean", 0, z)) %>% 
  select(z) -> pred_elev

# tb_resp_locs %>% 
#   left_join(as_tibble(pred_elev), by = c("longitude", "latitude")) -> pred_elev

pred_elev %>% 
  st_crop(w_us) %>% 
  as_tibble() -> pred_elev

# join with resp_precip by longitude and latitude



# PREDICTORS: NON LOCAL ------------------------------------------------------------------------


## SST ----

"~/bucket_mine/era/monthly/era5_monthly_mean_meansst_shifted.nc" %>% 
  read_ncdf() -> pred_sst

pred_sst %>% 
  filter(year(time) <= 2020) -> pred_sst


# *******

# old approach:

# expand_grid(pos_lon = seq(80/2, 1440, 80),
#             # pos_lat = seq(80/2, 721, 80)
#             pos_lat = c(
#               round(cumsum(seq(20, 110, length.out = 5))),
#               -round(cumsum(seq(20, 110, length.out = 5)))
#               ) %>%
#               sort(decreasing = T) %>%
#               {360 - .}
# 
#             ) %>%
# 
#   mutate(longitude = vect_lon[pos_lon],
#          latitude = vect_lat[pos_lat]) %>%
# 
#   mutate(r = row_number(),
#          r = str_pad(r, 3, "left", "0")) -> tb_pos_ref
# 
# ggplot() +
#   geom_stars(data = land) +
#   geom_point(data = tb_pos_ref, aes(x = longitude, y = latitude)) +
#   coord_equal()


# PCA (new approach)

land %>% 
  st_as_sf(as_points = F, merge = T) -> land_pol

set.seed(1)
st_sample(land_pol %>% filter(a == "ocean"), 5000) -> rnd_pts

# aggregate + detrend
plan(multicore, workers = 4)

pred_sst %>% aggregate(by = "years", mean, na.rm = T) -> pred_sst_annual

pred_sst_annual %>%
  st_apply(c("longitude","latitude"), function(x) x - mean(x), 
           .fname = "time",
           FUTURE = T) %>% 
  st_set_dimensions("time", 
                    values = st_get_dimension_values(pred_sst_annual, "time") %>% 
                      as_date()) -> pred_sst_detrended


# calculate pca
pca_maps(pred_sst_detrended, rnd_pts) -> pca

pca$pca_model$sdev^2 -> eigs
(eigs/sum(eigs))[1:10] # proportion var explained

# polygonize
pca$pca_rast %>% 
  select(1:2) -> pca_rast

# pca_rast %>% 
#   select(1) %>% 
#   view_era()

pca_rast %>%
  as_tibble() %>%
  mutate(across(PC1:PC2, ~ifelse(is.na(.x), -999, .x))) -> tb
  
set.seed(1)
kmeans(tb[,-c(1:2)], centers = 9, iter.max = 1000, nstart = 10) -> km
  
km$cluster %>%
  matrix(dim(pca_rast[[1]])) %>% 
  st_as_stars() -> km_rast
  
st_dimensions(km_rast) <- st_dimensions(land)

c(km_rast, land) %>% 
  mutate(A1 = ifelse(a == "land", NA, A1)) -> km_rast

# km_rast %>%
#   view_era()

km_rast %>% 
  st_warp(st_as_stars(st_bbox(), dx = 2)) %>% 
  st_as_sf(as_points = F, merge = T) -> km_pol

# mapview(km_pol)

# create points
km_pol %>% 
  filter(A1 != 7) %>% 
  mutate(area = st_area(.) %>% set_units(km^2)) %>% #mapview()
  filter(area > set_units(14e5, km^2)) %>% 
  group_by(A1) %>% 
  filter(area %in% sort(area, decreasing = T)[c(1:3)]) %>%
  st_point_on_surface() -> km_pts

mapview::mapview(km_pol) + km_pts

# *******

# ASSEMBLE TABLE

assemble_tb_nonlocal(pred_sst, km_pts) -> pred_sst_tb
# join by time



##  GEOPOTENTIAL 500 ----

"~/bucket_mine/era/monthly/era5_monthly_mean_geopot500hpa_shifted.nc" %>% 
  read_ncdf() -> pred_geopot500

pred_geopot500 %>% 
  setNames("geopot500") -> pred_geopot500


# *******

# PCA

set.seed(1)
st_sample(land_pol %>% st_bbox() %>% st_as_sfc(), 5000) -> rnd_pts

# aggregate + detrend
plan(multicore, workers = 4)

pred_geopot500 %>% aggregate(by = "years", mean, na.rm = T) -> pred_geopot500_annual

pred_geopot500_annual %>%
  st_apply(c("longitude","latitude"), function(x) x - mean(x), 
           .fname = "time",
           FUTURE = T) %>% 
  st_set_dimensions("time", 
                    values = st_get_dimension_values(pred_geopot500_annual, "time") %>% 
                      as_date()) -> pred_geopot500_detrended

# calculate pca
pca_maps(pred_geopot500_detrended, rnd_pts) -> pca

pca$pca_model$sdev^2 -> eigs
(eigs/sum(eigs))[1:10] # proportion var explained

# polygonize
pca$pca_rast %>% 
  select(1:3) -> pca_rast

pca_rast %>% 
  select(2) %>% 
  view_era()

pca_rast %>%
  as_tibble() -> tb

set.seed(1)
kmeans(tb[,-c(1:2)], centers = 9, iter.max = 1000, nstart = 10, algorithm = "Lloyd") -> km

km$cluster %>%
  matrix(dim(pca_rast[[1]])) %>% 
  st_as_stars() -> km_rast

st_dimensions(km_rast) <- st_dimensions(land)

km_rast %>%
  view_era()

km_rast %>% 
  st_warp(st_as_stars(st_bbox(), dx = 2)) %>% 
  st_as_sf(as_points = F, merge = T) -> km_pol

mapview::mapview(km_pol)

# create points
km_pol %>% 
  
  mutate(area = st_area(.) %>% set_units(km^2)) %>% #mapview
  filter(area > set_units(5e5, km^2)) -> km_pol
  
km_pol %>% 
  filter(area > set_units(24e6, km^2)) -> km_pol_large

tibble(xmin = seq(-180, 180, 360/3)[-4],
       xmax = seq(-180, 180, 360/3)[-1]) %>% 
  
  pmap(function(xmin, xmax){
    
    km_pol_large %>% 
      st_crop(xmin = xmin,
              xmax = xmax,
              ymin = -90,
              ymax = 90) %>% 
      st_point_on_surface()
    
  }) %>% 
  bind_rows() -> km_pts_large

km_pol %>% 
  filter(area < set_units(24e6, km^2)) %>% 
  group_by(A1) %>% 
  filter(area %in% sort(area, decreasing = T)[c(1:3)]) %>%
  st_point_on_surface() %>% 
  bind_rows(km_pts_large) -> km_pts

mapview::mapview(km_pol) + km_pts


# *******

# ASSEMBLE TABLE

assemble_tb_nonlocal(pred_geopot500, km_pts) -> pred_geopot500_tb

# join by time



# JOIN TABLES ---------------------------------------------------------------------------------

resp_precip %>%
  
  left_join(pred_precip, by = c("pos_lon", "pos_lat", "time")) %>% 
  left_join(pred_temp, by = c("pos_lon", "pos_lat", "time")) %>%
  left_join(pred_elev %>% select(-longitude, -latitude), by = c("pos_lon", "pos_lat")) %>% 
  
  left_join(pred_sst_tb, by = "time") %>% 
  left_join(pred_geopot500_tb, by = "time") %>%
   
  select(-pos_lon, -pos_lat) %>% 
  select(time, tp, everything()) -> tb_f

tb_f[!complete.cases(tb_f),]

saveRDS(tb_f, "output/tb_feat_sel.rds")




# TRAIN MODEL ---------------------------------------------------------------------------------

# source("~/00-mount.R")
source("scripts/00_setup.R")
library(tidymodels)

tb_f <- readRDS("output/tb_feat_sel.rds")

# # add month
# tb_f %>% 
#   mutate(mth = month(time) %>% factor()) -> tb_f

# split data
tb_f %>% 
  pull(time) %>% 
  unique() -> dates_sub

set.seed(1)
sample(dates_sub, length(dates_sub)/5) -> dates_test

tb_f %>% 
  filter(!time %in% dates_test) -> tb_train

tb_f %>% 
  filter(time %in% dates_test) -> tb_test

# recipe
recipe(tp ~ ., data = tb_train) %>% 
  update_role(time, new_role = "time") %>% 
  step_date(time, features = "month") -> rec

# model spec
set.seed(111)
rand_forest(trees = 10) %>% 
  set_engine("randomForest", importance = T) %>%
  set_mode("regression") -> mod

# workflow
workflow() %>% 
  add_model(mod) %>% 
  add_recipe(rec) -> wflow


# fit
wflow %>% 
  fit(data = tb_train) -> rf_fit



augment(rf_fit, tb_test) %>% select(time, longitude, latitude, tp, starts_with(".pre"))





mod %>% 
  fit(tp ~ ., data = tb_f[,-1]) -> rf_fit

predict(rf_fit, 
        tb_f[,-1])












library(randomForest)

readRDS("output/tb_feat_sel.rds") -> tb_f

tb_f %>% 
  mutate(mth = month(time)) -> tb_f

plan(multicore, workers = 18)

tic()
future_map(seq_len(18), function(x){
  randomForest(formula = tp ~ .,
               data = tb_f %>% select(-time),
               importance = T,
               ntree = 100)
}) -> rf_model
toc()

rf_model %>% 
  do.call(combine, .) -> rf_model_comb

saveRDS(rf_model_comb, "output/rf_model_feat_sel.rds")



