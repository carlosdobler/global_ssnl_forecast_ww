
# source("~/00-mount.R")
source("scripts/00_setup.R")
source("scripts/functions.R")

seq(as_date("1979-01-01"),
    as_date("2020-12-31"),
    "1 month") -> dates

dates[month(dates) %in% c(3,6,9,12) & 
        year(dates) > 1980 &
        dates < as_date("2020-12-01")] -> dates_sub

dates %>% 
  days_in_month() %>% 
  unname() -> days_mth



# LAND MASK -----------------------------------------------------------------------------------

rast_reference_0.05 <- 
  c(st_point(c(-180.125, -90.125)),
    st_point(c(179.875, 90.125))) %>% 
  st_bbox() %>% 
  st_set_crs(4326) %>% 
  st_as_stars(dx = 0.05, dy = 0.05, values = -9999)

rast_reference_era <- 
  c(st_point(c(-180.125, -90.125)),
    st_point(c(179.875, 90.125))) %>% 
  st_bbox() %>%
  st_set_crs(4326) %>% 
  st_as_stars(dx = 0.25, dy = 0.25, values = -9999)

land <- 
  "~/bucket_mine/misc_data/ne_50m_land/ne_50m_land.shp" %>% 
  st_read() %>% 
  mutate(a = 1) %>% 
  select(a) %>% 
  st_rasterize(rast_reference_0.05)

land <- 
  land %>% 
  st_warp(rast_reference_era, use_gdal = T, method = "mode") %>% 
  setNames("a") %>% 
  mutate(a = ifelse(a == -9999, "ocean", "land") %>% factor())

land %>% 
  st_set_dimensions(c(1,2), names = c("longitude", "latitude")) -> land
  
st_get_dimension_values(rast_reference_era, "x") -> vect_lon
st_get_dimension_values(rast_reference_era, "y") -> vect_lat



# RESPONSE VAR --------------------------------------------------------------------------------

"~/bucket_mine/era/monthly/era5_monthly_mean_daily_precip_shifted.nc" %>% 
  read_ncdf() %>% 
  
  mutate(tp = set_units(tp, mm)) %>% 
  
  st_apply(c(1,2), function(x){
    x * days_mth
  },
  .fname = "time") %>%
  aperm(c(2,3,1)) %>%
  st_set_dimensions("time", values = dates) -> precip


# *****************

precip %>% 
  aggregate(by = c(dates_sub, as_date("2020-12-01")), sum) -> precip_seas

precip_seas %>% 
  slice(time, -160) -> precip_seas

expand_grid(pos_lon = seq(40, 1420, 40),
            pos_lat = seq(40, 701, 40)) %>%
  
  mutate(longitude = vect_lon[pos_lon],
         latitude = vect_lat[pos_lat]) %>%
  
  left_join(as_tibble(land), by = c("longitude", "latitude")) %>% 
  filter(a == "land") %>%
  select(-a) %>% 
  
  filter(abs(latitude) < 67) -> tb_resp_locs 
  
ggplot() +
  geom_stars(data = land) +
  geom_point(data = tb_resp_locs, aes(x = longitude, y = latitude)) +
  coord_equal()
  
left_join(tb_resp_locs, as_tibble(precip_seas), by = c("longitude", "latitude")) -> resp_precip



# PREDICTORS: LOCAL -------------------------------------------------------------------------------

lg <- 6
train_window <- 6


## PRECIP (MEMORY) ----

# plan(multicore, workers = 2)

tb_resp_locs %>% 
  
  mutate(r = row_number(),
         r = str_pad(r, 3, "left", "0")) %>%
  
  # future_pmap_dfr(function(pos_lon, pos_lat, r, ...){
  pmap_dfr(function(pos_lon, pos_lat, r, ...){
    
    print(str_glue("{r} / {nrow(tb_resp_locs)}"))
    
    precip[,
           (pos_lon-1):(pos_lon+1), 
           (pos_lat-1):(pos_lat+1),
    ] %>%
      st_apply(3, mean, na.rm = T) %>% 
      pull(1) -> x_ts
    
    if(!anyNA(x_ts)){
      
      dates_sub %>% 
        map_dfr(function(d){
          
          w <- x_ts[(which(dates == d) - lg - train_window):(which(dates == d) - lg - 1)]
          
          tibble(
            
            pos_lon = pos_lon,
            pos_lat = pos_lat,
            tp_avg6mths = mean(w),
            time = d
            
          )
        })
    }
  }) -> pred_precip
# join with resp_precip by pos_lon, pos_lat, and time


## SURF TEMP ----

"~/bucket_mine/era/monthly/era5_monthly_mean_meantemp_shifted.nc" %>% 
  read_ncdf() -> pred_temp

plan(multicore, workers = 4)

pred_temp %>% 
  st_apply(3, foc, w = matrix(1,3,3),
           FUTURE = T) -> pred_temp

tb_resp_locs %>% 
  
  mutate(r = row_number(),
         r = str_pad(r, 3, "left", "0")) %>%
  
  # future_pmap_dfr(function(pos_lon, pos_lat, r, ...){
  pmap_dfr(function(pos_lon, pos_lat, r, ...){
    
    print(str_glue("{r} / {nrow(tb_resp_locs)}"))
    
    pred_temp[ , pos_lon, pos_lat, ] %>%
      pull(1) -> x_ts
    
    if(!anyNA(x_ts)){
      
      dates_sub %>% 
        map_dfr(function(d){
          
          w <- x_ts[(which(dates == d) - lg - train_window):(which(dates == d) - lg - 1)]
          
          tibble(
            
            pos_lon = pos_lon,
            pos_lat = pos_lat,
            tp_avg6mths = mean(w),
            time = d
            
          )
        })
    }
  }) -> pred_temp
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

tb_resp_locs %>% 
  left_join(as_tibble(pred_elev), by = c("longitude", "latitude")) -> pred_elev
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

# detrend
plan(multicore, workers = 4)

pred_sst %>%
  st_apply(c("longitude","latitude"), function(x) x - mean(x), 
           .fname = "time",
           FUTURE = T) %>% 
  st_set_dimensions("time", values = dates) -> pred_sst_detrended


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

pca$pca_rast %>% 
  select(1:2) -> pca_rast

pca_rast %>%
  as_tibble() %>%
  mutate(across(PC1:PC2, ~ifelse(is.na(.x), -999, .x))) -> tb
  
kmeans(tb[,-c(1:2)], centers = 8, iter.max = 1000, nstart = 5, algorithm = "Lloyd") -> km
  
km$cluster %>%
  matrix(dim(pca_rast[[1]])) %>% 
  st_as_stars() -> km_rast
  
st_dimensions(km_rast) <- st_dimensions(land)

c(km_rast, land) %>% 
  mutate(A1 = ifelse(a == "land", NA, A1)) -> km_rast

km_rast %>% 
  view_era()

  


pca_rast %>%
  view_era()

pca_rast %>% 
  mutate(PC2 = case_when(PC2 > quantile(PC2, 0.98, na.rm = T) ~ quantile(PC2, 0.95, na.rm = T),
                         TRUE ~ PC2)) %>% 
  view_era()

pca$pca_rast %>% 
  select(1) -> pca_rast
  
pca_rast %>%
  view_era()

pca_rast %>%
  st_warp(st_as_stars(st_bbox(), dx = 0.25)) %>% 
  st_contour(breaks = seq(quantile(.[[1]], 0, na.rm = T), 
                          quantile(.[[1]], 1, na.rm = T), 
                          length.out = 5)) %>% 
  st_cast("POLYGON") -> pca_pol

pca_pol %>%
  mutate(area = st_area(pca_pol) %>% set_units(km^2)) %>% #mapview::mapview()
  filter(area > set_units(35e5, km^2)) -> pca_pol
  
st_point_on_surface(pca_pol) -> pca_pts


# *******

# ASSEMBLE TABLE

# plan(multicore, workers = 8)

pca_pts %>% 
  st_coordinates() %>% 
  as_tibble() %>% 
  rowwise() %>%
  mutate(pos_lon = which.min(abs(vect_lon - X)),
         pos_lat = which.min(abs(vect_lat - Y))) %>% 
  ungroup() %>%
  
  mutate(r = row_number() %>% 
           str_pad(2, "left", "0")) %>% 
  
  # future_pmap_dfc(function(pos_lon, pos_lat, ...){
  pmap_dfc(function(pos_lon, pos_lat, r, ...){
    
    print(str_glue("{r}"))
    
    pred_sst[,
             (pos_lon-2):(pos_lon+2), # 5 x 5 matrix
             (pos_lat-2):(pos_lat+2),
    ] %>%
      st_apply(3, mean, na.rm = T) %>% 
      pull(1) -> x_ts
    
    if(!anyNA(x_ts)){
      
      dates_sub %>% 
        map_dfr(function(d){
          
          w <- x_ts[(which(dates == d) - lg - train_window):(which(dates == d) - lg - 1)]
          
          tibble(
            
            avg6mths = mean(w)
            
          ) -> tb
          
          tb %>% 
            rename_with(~str_glue("{names(pred_sst)}_{names(tb)}_{r}"))
        })
    }
  }) -> pred_sst_tb

pred_sst_tb %>% 
  mutate(time = dates_sub,
  #        mth = month.abb[month(dates_sub)] %>% factor(),
         .before = 1) -> pred_sst_tb

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

# detrend
plan(multicore, workers = 4)
pred_geopot500 %>% 
  st_apply(c(1,2), function(x) x - mean(x), 
           .fname = "time",
           FUTURE = T) %>% 
  st_set_dimensions("time", values = dates) -> pred_geopot500_detrended

# calculate pca
pca_maps(pred_geopot500_detrended, rnd_pts) -> pca

pca$pca_model$sdev^2 -> eigs
(eigs/sum(eigs))[1:2] # proportion var explained

pca$pca_rast %>% 
  select(1) -> pca_rast

pca_rast %>%
  view_era()

pca_rast %>%
  st_warp(st_as_stars(st_bbox(), dx = 0.25)) %>%
  st_contour(breaks = seq(quantile(.[[1]], 0, na.rm = T), 
                          quantile(.[[1]], 1, na.rm = T), 
                          length.out = 7)) %>% 
  st_cast("POLYGON") -> pca_pol

mapview(pca_pol)

pca_pol %>%
  mutate(area = st_area(pca_pol) %>% set_units(km^2)) -> pca_pol
  
pca_pol %>% 
  filter(area < set_units(2e7, km^2)) %>% 
  st_point_on_surface() -> pca_pts_small
  
pca_pol %>% 
  filter(area > set_units(2e7, km^2)) -> pca_pts_large

tibble(xmin = seq(-180, 180, 360/3)[-4],
       xmax = seq(-180, 180, 360/3)[-1]) %>% 
  
  pmap(function(xmin, xmax){
    
    pca_pts_large %>% 
      st_crop(xmin = xmin,
              xmax = xmax,
              ymin = -90,
              ymax = 90) %>% 
      st_point_on_surface()
    
  }) %>% 
  bind_rows() -> pca_pts_large

bind_rows(pca_pts_small,
          pca_pts_large) -> pca_pts

mapview(pca_pol) + pca_pts

# *******

# ASSEMBLE TABLE

# plan(multicore, workers = 8)

pca_pts %>% 
  st_coordinates() %>% 
  as_tibble() %>% 
  rowwise() %>%
  mutate(pos_lon = which.min(abs(vect_lon - X)),
         pos_lat = which.min(abs(vect_lat - Y))) %>% 
  ungroup() %>%
  
  mutate(r = row_number() %>% 
           str_pad(2, "left", "0")) %>% 
  
  # future_pmap_dfc(function(pos_lon, pos_lat, ...){
  pmap_dfc(function(pos_lon, pos_lat, r, ...){
    
    print(str_glue("{r}"))
    
    pred_geopot500[,
                   (pos_lon-2):(pos_lon+2), # 5 x 5 matrix
                   (pos_lat-2):(pos_lat+2),
    ] %>%
      st_apply(3, mean, na.rm = T) %>% 
      pull(1) -> x_ts
    
    if(!anyNA(x_ts)){
      
      dates_sub %>% 
        map_dfr(function(d){
          
          w <- x_ts[(which(dates == d) - lg - train_window):(which(dates == d) - lg - 1)]
          
          tibble(
            
            avg6mths = mean(w)
            
          ) -> tb
          
          tb %>% 
            rename_with(~str_glue("{names(pred_geopot500)}_{names(tb)}_{r}"))
        })
    }
  }) -> pred_geopot500_tb

pred_geopot500_tb %>% 
  mutate(time = dates_sub,
         #        mth = month.abb[month(dates_sub)] %>% factor(),
         .before = 1) -> pred_geopot500_tb

# join by time



# JOIN TABLES ---------------------------------------------------------------------------------

resp_precip %>%
  
  left_join(pred_precip, by = c("pos_lon", "pos_lat", "time")) %>% 
  left_join(pred_temp, by = c("pos_lon", "pos_lat", "time")) %>%
  left_join(pred_elev %>% select(-r, -longitude, -latitude), by = c("pos_lon", "pos_lat")) %>% 
  
  left_join(pred_sst_tb %>% select(-mth), by = "time") %>% 
  left_join(pred_geopot500_tb %>% select(-mth), by = "time") %>% 
   
  select(-pos_lon, -pos_lat) %>% 
  select(time, tp, everything()) -> tb_f

tb_f[!complete.cases(tb_f),]

saveRDS(tb_f, "output/tb_feat_sel.rds")




# TRAIN MODEL ---------------------------------------------------------------------------------

source("scripts/00_setup.R")
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



