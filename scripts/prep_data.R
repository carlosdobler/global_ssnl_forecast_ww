source("scripts/00_setup.R")
library(tsfeatures)


# PRELIM

seq(as_date("1979-01-01"),
    as_date("2020-12-31"),
    "1 month") -> dates

dates %>% 
  days_in_month() %>% 
  unname() -> days_mth

"~/bucket_mine/misc_data/ne_50m_land/" %>% 
  st_read() %>% 
  mutate(a = 1) %>% 
  select(a) -> land

# c(st_point(c(-180.125, -90.125)),
#   st_point(c(179.875, 90.125))) %>% 

c(st_point(c(-126.625, 32.625)),
  st_point(c(-115.125, 47.875))) %>% 
  
  st_bbox() %>% 
  st_set_crs(4326) %>%
  {st_rasterize(land, 
                st_as_stars(.,  dx = 0.25, dy = 0.25, values = NA))} -> land_r


# *****************


# RESPONSE VAR: SEASONAL PRECIPITATION

"~/bucket_mine/era/monthly/era5_monthly_mean_daily_precip_shifted.nc" %>% 
  read_ncdf() -> precip

precip %>% 
  st_crop(land_r) -> precip

precip %>% 
  mutate(tp = set_units(tp, mm)) -> precip

precip %>%
  st_apply(c(1,2), function(x){
    x * days_mth
  },
  .fname = "time") %>%
  aperm(c(2,3,1)) %>%
  st_set_dimensions("time", values = dates) -> precip


dates[month(dates) %in% c(3,6,9,12) & year(dates) > 1980] -> dates_sub

precip %>% 
  aggregate(by = dates_sub, sum) -> precip

# equivalent:
# s_sub_seasonal %>%
#   slice(time, 1)
# 
# s_sub %>%
#   filter(time %in% seq(as_date("1980-03-01"), as_date("1980-05-01"), by = "months")) %>%
#   st_apply(c(1,2), sum)

precip %>% 
  slice(time, -160) -> precip

dates_sub[-length(dates_sub)] -> dates_sub

precip %>% 
  as_tibble() -> tb_response

# set.seed(1)
# tb_response %>% 
#   group_by(time) %>% 
#   nest() %>% 
#   mutate(data_sub = map(data, function(df){
#     df %>% 
#       slice_sample(prop = 0.25)
#   })) %>%
#   select(-data) %>% 
#   unnest(data_sub) %>% 
#   ungroup() -> tb_response


# *****************


# PREDICTOR VARS

# ELEVATION + ASPECT

"~/bucket_mine/era/monthly/era5_geopotential_shifted.nc" %>% 
  read_ncdf() %>% 
  adrop() %>% 
  mutate(z = set_units(z, NULL),
         z = z/9.80665,
         z = ifelse(z < 0, 0 , z)) -> z
  
as(z, "Raster") %>% 
  raster::terrain("aspect", unit = "degrees") %>% 
  st_as_stars() %>% 
  st_set_dimensions(which = c(1,2), names = c("longitude", "latitude")) -> aspect

aspect %>%
  st_crop(land_r) %>%
  setNames("asp") %>% 
  mutate(asp = cut(asp,
                   c(0, seq(22.5, 360, 45), 360),
                   include.lowest = T,
                   labels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW", "N"))) -> aspect

aspect %>% 
  mutate(asp = fct_expand(asp, "no_slope")) -> aspect

aspect[is.na(land_r)] <- "no_slope"

z %>% 
  st_crop(land_r) -> z

z[is.na(land_r)] <- 0

inner_join(as_tibble(z),
           as_tibble(aspect),
           by = c("longitude", "latitude")) -> tb_pred_1


# *****

# SST

lg <- 6 # lag

"~/bucket_mine/era/monthly/era5_monthly_mean_meansst_shifted.nc" %>% 
  read_ncdf() -> sst

sst %>% 
  filter(year(time) <= 2020) -> sst

plan(multicore, workers = 8)

expand_grid(pos_lon = seq(1, 1440, 80),
            pos_lat = seq(1, 721, 80)) %>%
  mutate(r = row_number(),
         r = str_pad(r, 3, "left", "0")) %>%

  # .[1:10,] %>%
  
  future_pmap_dfc(function(pos_lon, pos_lat, r){
    
    # pos_lon <- 81
    # pos_lat <- 321
    
    sst[, pos_lon, pos_lat, ] %>%
      pull(1) %>%
      as.vector() -> x_ts
    
    if(!anyNA(x_ts)){
      
      # print(str_glue("{r} / 180"))
      
      dates_sub %>% 
        map_dfr(function(d){
          
          w <- x_ts[(which(dates == d) - lg - 12): # training period length 
                    (which(dates == d) - lg - 1)]
          
          tibble(
            
            mean1term = mean(w[1:6]),
            mean2term = mean(w[7:12]),
            slope = lm(w ~ seq_len(12)) %>% coefficients() %>% .[2] %>% unname(),
            stl_features(w)[3:6] %>% enframe() %>% pivot_wider()
            
          ) -> tb
          
          tb %>% 
            rename_with(~str_glue("{names(tb)}_{r}"))
        })
    }
  }) -> tb_pred_2

tb_pred_2 %>% 
  mutate(time = dates_sub,
         mth = month.abb[month(dates_sub)] %>% factor(),
         .before = 1) -> tb_pred_2


# *****************


tb_response %>%
  left_join(tb_pred_1, by = c("longitude", "latitude")) %>%
  left_join(tb_pred_2, by = "time") -> tb_f

saveRDS(tb_f, "output/tb_f.rds")
