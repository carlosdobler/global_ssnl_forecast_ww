
# SETUP -------------------------------------------------------------------------------------------

# source("~/00-mount.R")

source("scripts/00_setup.R")

seq(as_date("1979-01-01"),
    as_date("2021-12-31"),
    "1 month") -> dates

dates %>% 
  days_in_month() %>% 
  unname() -> days_mth

# region
c(st_point(c(-126.625, 32.625)),
  st_point(c(-115.125, 48.875))) %>% 
  st_bbox() %>% 
  st_set_crs(4326) -> region

c(st_point(c(region[1]-0.25, region[2]-0.25)),
  st_point(c(region[3]+0.25, region[4]+0.25))) %>% 
  st_bbox() %>% 
  st_set_crs(4326) -> region_plus



# LAND MASK ---------------------------------------------------------------------------------------

region_plus %>% 
  st_as_stars(dx = 0.05, dy = 0.05, values = -9999) -> rast_reference_0.05

region_plus %>% 
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



# RESPONSE VAR --------------------------------------------------------------------------------

"~/bucket_mine/era/monthly/era5_monthly_mean_daily_precip_1979-2021_shifted.nc" %>% 
  read_ncdf(ncsub = cbind(start = c(1,1,1),
                          count = c(NA, NA, 1))) -> s_proxy

which.min(abs(st_get_dimension_values(s_proxy, "longitude", center = F) - region_plus[1])) -> lon_start
which.min(abs(st_get_dimension_values(s_proxy, "longitude", center = F) - region_plus[3])) - lon_start -> lon_count

which.min(abs(st_get_dimension_values(s_proxy, "latitude", center = F) - region_plus[4])) -> lat_start
which.min(abs(st_get_dimension_values(s_proxy, "latitude", center = F) - region_plus[2])) - lat_start -> lat_count


"~/bucket_mine/era/monthly/era5_monthly_mean_daily_precip_1979-2021_shifted.nc" %>% 
  read_ncdf(ncsub = cbind(start = c(lon_start, lat_start, 1),
                          count = c(lon_count, lat_count, NA))) -> precip

st_dimensions(precip)$longitude$offset <- st_get_dimension_values(s_proxy, "longitude", center = F)[lon_start]
st_dimensions(precip)$latitude$offset <- st_get_dimension_values(s_proxy, "latitude", center = F)[lat_start]

precip %>% 
  mutate(tp = set_units(tp, mm)) %>%
  
  st_apply(c(1,2), 
           function(x) x * days_mth,
           .fname = "time") %>%
  aperm(c(2,3,1)) %>%
  st_set_dimensions("time", values = dates) -> precip


# Aggregate precip into seasons
# - each date in "by" is the left bound (included)
# - the next date is the right bound (excluded)
precip %>% 
  aggregate(by = dates[month(dates) %in% c(3,6,9,12)], 
            sum) %>%
  aperm(c(2,3,1)) -> precip_seas

# remove last time step (Dec; incomplete quarter)
precip_seas %>% 
  slice(time, -(dim(precip_seas)[3])) -> precip_seas

# only land
precip_seas[land == "ocean"] <- NA

# response variable table
precip_seas %>%
  st_crop(region) %>% 
  as_tibble() -> tb_resp

# sample (regular)
expand_grid(longitude = seq(min(tb_resp$longitude), max(tb_resp$longitude), 0.25*3),
            latitude = seq(min(tb_resp$latitude), max(tb_resp$latitude), 0.25*3)) %>% 
  left_join(tb_resp, by = c("longitude", "latitude")) -> tb_resp

tb_resp %>% filter(!is.na(tp)) -> tb_resp



# PREDICTOR VARS I: LOCAL --------------------------------------------------------------------------

lg <- 6 # lag
train_window <- 6 # training window width


# 1. PRECIPITATION (MEMORY)
# 6-month accumulated precip of 3x3 window (mean) around sampled cells
# 6-month lag between date to predict and accum precip

# apply SPATIAL moving window (all region; sample later)
foc <- function(x, w) {
  raster::as.matrix(raster::focal(raster::raster(x), w))
}

precip %>%
  st_apply(3, foc, w = matrix(1/9, 3, 3)) -> precip_mw

# apply TEMPORAL moving window
precip_mw[land == "ocean"] <- NA

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
  aperm(c(2,3,1)) -> precip_mw

# add lag (match dates with response dates)
precip_mw %>% 
  st_set_dimensions("time", values = dates + months(lg+1)) -> precip_mw

# table
precip_mw %>%
  setNames("pred_tp") %>% 
  as_tibble() %>%
  {left_join(tb_resp, ., by = c("longitude", "latitude", "time"))} -> tb_f


# 2. TEMPERATURE
# Same approach as above





