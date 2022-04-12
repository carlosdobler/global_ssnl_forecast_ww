
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

expand_grid(pos_lon = seq(20, 1420, 20),
            pos_lat = seq(20, 701, 20)) %>%
  
  mutate(longitude = vect_lon[pos_lon],
         latitude = vect_lat[pos_lat]) %>%
  
  left_join(as_tibble(land), by = c("longitude", "latitude")) %>% 
  filter(a == "land") %>%
  select(-a) %>% 
  
  filter(abs(latitude) < 67) %>% 
  
  # {
  #   ggplot() +
  #     geom_stars(data = land) +
  #     geom_point(data = ., aes(x = longitude, y = latitude)) +
  #     coord_equal()
  # }
  
  left_join(as_tibble(precip_seas), by = c("longitude", "latitude")) -> resp_precip # 547 points / 86,973 obs



# PREDICTOR VAR: PRECIP (MEMORY) --------------------------------------------------------------
lg <- 6
train_window <- 6

# plan(multicore, workers = 2)

expand_grid(pos_lon = seq(20, 1420, 20),
            pos_lat = seq(20, 701, 20)) %>%
  
  mutate(longitude = vect_lon[pos_lon],
         latitude = vect_lat[pos_lat]) %>%
  
  left_join(as_tibble(land), by = c("longitude", "latitude")) %>% 
  filter(a == "land") %>%
  select(-a) %>% 
  
  filter(abs(latitude) < 67) %>% 
  
  mutate(r = row_number(),
         r = str_pad(r, 3, "left", "0")) -> tb_pos_ref_2 
  
tb_pos_ref_2 %>%
  # future_pmap_dfr(function(pos_lon, pos_lat, r, ...){
  pmap_dfr(function(pos_lon, pos_lat, r, ...){
    
    print(str_glue("{r} / {nrow(tb_pos_ref_2)}"))
    
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



# PREDICTOR VARS: TIME INVARIANT ----------------------------------------------------------------

# ELEVATION 

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

tb_pos_ref_2 %>% 
  left_join(as_tibble(pred_elev), by = c("longitude", "latitude")) -> pred_elev

# join with resp_precip by longitude and latitude



# PREDICTOR VARS: SPACE INVARIANT -------------------------------------------------------------

expand_grid(pos_lon = seq(80/2, 1440, 80),
            # pos_lat = seq(80/2, 721, 80)
            pos_lat = c(
              round(cumsum(seq(20, 110, length.out = 5))),
              -round(cumsum(seq(20, 110, length.out = 5)))
              ) %>% 
              sort(decreasing = T) %>% 
              {360 - .}
            
            ) %>%
  
  mutate(longitude = vect_lon[pos_lon],
         latitude = vect_lat[pos_lat]) %>% 
  
  mutate(r = row_number(),
         r = str_pad(r, 3, "left", "0")) -> tb_pos_ref

ggplot() +
  geom_stars(data = land) +
  geom_point(data = tb_pos_ref, aes(x = longitude, y = latitude)) +
  coord_equal()


# *****************


plan(multicore, workers = 10)

lg <- 6
train_window <- 6


# *****************

# SST

"~/bucket_mine/era/monthly/era5_monthly_mean_meansst_shifted.nc" %>% 
  read_ncdf() -> pred_sst

pred_sst %>% 
  filter(year(time) <= 2020) -> sst

sp_inv_pred_tb(pred_sst, tb_pos_ref) -> pred_sst

# join by time


# *******


# GEOPOTENTIAL 500

"~/bucket_mine/era/monthly/era5_monthly_mean_geopot500hpa_shifted.nc" %>% 
  read_ncdf() -> pred_geopot500

pred_geopot500 %>% 
  setNames("geopot500") -> pred_geopot500

sp_inv_pred_tb(pred_geopot500, tb_pos_ref) -> pred_geopot500



# *******


# DIVERGENCE 500

"~/bucket_mine/era/monthly/era5_monthly_mean_divergence500hpa_shifted.nc" %>% 
  read_ncdf() -> pred_div500

pred_div500 %>% 
  setNames("div500") -> pred_div500

sp_inv_pred_tb(pred_div500, tb_pos_ref) -> pred_div500



# JOIN TABLES ---------------------------------------------------------------------------------

resp_precip %>%
  left_join(pred_precip, by = c("pos_lon", "pos_lat", "time")) %>% 
  left_join(pred_elev %>% select(-r, -longitude, -latitude), by = c("pos_lon", "pos_lat")) %>% 
  left_join(pred_sst %>% select(-mth), by = "time") %>% 
  left_join(pred_geopot500 %>% select(-mth), by = "time") %>% 
  left_join(pred_div500 %>% select(-mth), by = "time") %>% 
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



