

# SETUP -----------------------------------------------------------------------

source("scripts/00_setup.R")

plan(multicore)
sf_use_s2(F)


era_grid <- 
  "/mnt/bucket_mine/era/monthly/mean-daily-precip/era5_mon_mean-daily-precip_1959.nc" %>% 
  read_ncdf(ncsub = cbind(start = c(1,1,1),
                          count = c(NA,NA,1))) %>% 
  adrop() %>%
  drop_units() %>% 
  setNames("a") %>%
  mutate(a = NA)



# FUNCTIONS -------------------------------------------------------------------


# Function to copy files to local disk and concatenate into a single file
copy_cat <- function(var_name){
  
  # Copy annual files to local disk
  "gsutil -m cp -r gs://clim_data_reg_useast1/era/monthly/{var_name} {dir_raw_data}" %>%
    str_glue() %>%
    system()
  
  # Concatenate
  ff <- 
    str_glue("{dir_raw_data}/{var_name}") %>%
    list.files(full.names = T) %>% 
    str_flatten(" ")
  
  "cdo cat {ff} {dir_raw_data}/{var_name}.nc" %>% 
    str_glue() %>% 
    system()
  
  # Delete annual files
  "{dir_raw_data}/{var_name}" %>% 
    str_glue() %>% 
    unlink(recursive = T)
  
}



# Function to obtain EOF based on a sample of locations
eof <- function(s_detrended, pts){
  
  # date vector
  d <- 
    st_get_dimension_values(s_detrended, "time") %>% 
    {str_glue("x_{.}")}
  
  # extract values 
  print(str_glue("Extracting values ..."))
  
  pts_val <- 
    aggregate(s_detrended, 
              pts, 
              function(x) x[1], as_points = F) %>% 
    st_as_sf()
  
  pts_val <- 
    bind_cols(
      pts_val %>% 
        st_coordinates() %>% 
        as_tibble(),
      
      pts_val %>% 
        st_drop_geometry() %>% 
        as_tibble() %>% 
        rename_with(~d)
    )
  
  # fit pca
  eof_model <- 
    pts_val %>%
    select(-c(X,Y)) %>% 
    na.omit() %>% 
    prcomp(scale. = T)
  
  # predict eof model to the rest of data
  print(str_glue("Predicting ..."))
  
  s_eof <- 
    s_detrended %>%
    st_set_dimensions("time", 
                      values = d) %>% 
    split("time") %>% 
    predict(eof_model)
  
  list(eof_model = eof_model,
       s_eof = s_eof)
  
}



# LAND MASK -------------------------------------------------------------------

{
  land <- 
    "/mnt/bucket_mine/misc_data/ne_50m_land/ne_50m_land.shp" %>% 
    st_read() %>% 
    mutate(a = 1) %>% 
    select(a) %>% 
    st_rasterize(dx = 0.25)
  
  # Convert to 0-360
  
  land_1 <- 
    land %>% 
    filter(x < 0) %>% 
    st_set_dimensions(which = "x",
                      values = st_get_dimension_values(.,
                                                       "x",
                                                       center = F)+360) %>% 
    st_set_crs(4326)
  
  land_2 <- 
    land %>% 
    filter(x >= 0)
  
  land <- 
    list(land_1, land_2) %>% 
    map(as, "SpatRaster") %>% 
    do.call(terra::merge, .) %>%
    st_as_stars(proxy = F)
  
  rm(land_1, land_2)

  # Warp to ERA grid
  
  land <- 
    land %>% 
    st_warp(era_grid) %>% 
    setNames("a") %>% 
    mutate(a = ifelse(is.na(a), "ocean", "land") %>% factor())
  
  land_pol <- 
    land %>% 
    st_as_sf(as_points = F, merge = T)
  
}



# PROCESS PRECIPITATION -------------------------------------------------------

# read files
precip_raw <-
  str_glue("{dir_raw_data}/mean-daily-precip.nc") %>%
  read_ncdf(proxy = T)

# get dates
dates <-
  st_get_dimension_values(precip_raw, "time") %>%
  as_date()

days_mth <-
  dates %>%
  days_in_month() %>%
  unname()




## Regionalization ********************

# Segment precipitation into homogeneous areas; 
# obtain centroid of those areas; centroids will 
# be used as sample points to train model


# copy_cat(mean-daily-precip)


{
   
  # calculate features to regionalize
  precip_features <-
    precip_raw[, 930:1100, 130:300, ] %>%
    st_apply(c(1,2), function(x){

      x_ts <- ts(x, frequency = 12) # format as ts

      # extract seasonality via decomposition
      x_season <-
        stl(x_ts, s.window = "periodic")$time.series[1:12,1]

      c(avg = mean(x),
        seas = diff(range(x_season)),
        peak_pos = which.max(x_season),
        trough_pos = which.min(x_season))

    },
    FUTURE = T,
    .fname = "stats") %>%
    st_as_stars(proxy = F) %>%
    split("stats")


  # remove ocean
  precip_features[is.na(st_warp(land, precip_features) %>%
                          mutate(a = ifelse(a == "ocean", NA, a)))] <- NA

  # check correlations
  precip_features %>%
    as_tibble() %>%
    select(-(1:2)) %>%
    corrr::correlate() %>%
    corrr::shave() %>%
    corrr::stretch() %>%
    arrange(desc(abs(r)))

  set.seed(1)
  reg <-
    precip_features %>%

    # rescale and weight
    mutate(avg = scales::rescale(avg)*50, # weight
           seas = scales::rescale(seas),
           peak_pos = scales::rescale(peak_pos),
           trough_pos = scales::rescale(trough_pos)
    ) %>%

    supercells::supercells(step = 3, # 4
                           compactness = 3,
                           iter = 20)

  # obtain points
  sample_pts <-
    st_drop_geometry(reg) %>%
    select(2,3) %>% # centroid coords
    st_as_sf(coords = c("x", "y")) %>%
    st_set_crs(4326)


  # st_write(sample_pts, "output/sample_pts.gpkg")
  
}



## Process response variable ***********

# sample_pts <- 
#   st_read("output/sample_pts.gpkg")

# aggregate into seasons/quarters
precip_quarterly <- 
  precip_raw[, 930:1100, 130:300, ] %>% 
  st_apply(c(1,2), function(x){
    
    # monthly totals
    x_mon <- x * days_mth
    
    # 3-month moving window (quarter)
    slider::slide_dbl(x_mon,
                      sum,
                      .complete = T,
                      .after = 2)
    
  },
  FUTURE = T,
  .fname = "time") %>% 
  
  st_as_stars(proxy = F) %>% 
  
  st_set_dimensions("time", values = dates) %>% 
  aperm(c(2,3,1)) %>% 
  
  mutate(tp = tp %>% 
           set_units(m) %>% 
           set_units(mm) %>% 
           set_units(NULL))


# calculate anomalies
# precip_quarterly_anom <- 
#   precip_quarterly %>% 
#   st_apply(c(1,2), function(x){
#     
#     monthly_mean <- 
#       matrix(x, ncol = 12, byrow = T) %>% 
#       apply(2, mean, na.rm = T)
#     
#     anomaly <- x - monthly_mean 
#     
#     return(anomaly)
#     
#   },
#   FUTURE = T,
#   .fname = "time") %>% 
#   aperm(c(2,3,1))




# precip_quarterly_anom <- 
#   precip_quarterly_anom %>% 
#   st_set_dimensions("time", values = dates - months(lead_time + 1)) %>% 
#   setNames("response") %>% 
#   slice(time, -(1:7))

precip_quarterly <- 
  precip_quarterly %>% 
  setNames("response")


# extract from sample points
tb_reponse <-
  aggregate(precip_quarterly, #_anom,
            sample_pts,
            function(x) x[1],
            as_points = FALSE) %>%
  as_tibble() %>% 
  mutate(as_tibble(st_coordinates(geometry))) %>% 
  select(-geometry)




## Process predictor variable *********

# (memory of atmosphere)

precip_predictor <- 
  precip_raw[, 930:1100, 130:300, ] %>% 
  
  # aggregate into 6-month precip
  st_apply(c(1,2), function(x){
    
    x_mon <- x * days_mth
    
    slider::slide_dbl(x_mon,
                      sum,
                      .complete = T,
                      .before = 5) # window size
    
  },
  FUTURE = T,
  .fname = "time") %>% 
  
  st_as_stars(proxy = F) %>% 
  aperm(c(2,3,1)) %>% 
  st_set_dimensions("time", values = dates) %>% 
  mutate(tp = tp %>% 
           set_units(m) %>% 
           set_units(mm) %>% 
           set_units(NULL)) %>% 
  setNames("precip")


# adjust time dimension
# align response to predictors by lagging
lead_time <- 6

precip_predictor <- 
  precip_predictor %>% 
  st_set_dimensions("time", values = dates + months(lead_time + 1))



# extract from sample points
tb_pred_precip <- 
  aggregate(precip_predictor,
            sample_pts,
            function(x) x[1],
            as_points = FALSE) %>%
  as_tibble() %>% 
  mutate(as_tibble(st_coordinates(geometry))) %>% 
  select(-geometry)




# PROCESS TEMPERATURE ---------------------------------------------------------

# cdo -b F32 cat {ff} {dir_raw_data}/mean-tasmean.nc
# copy_cat("mean-tasmean")


# read files 
tasmean_raw <- 
  str_glue("{dir_raw_data}/mean-tasmean.nc") %>% 
  read_ncdf(proxy = T)


tasmean_predictor <- 
  tasmean_raw[, 930:1100, 130:300, ] %>% 
  st_apply(c(1,2), function(x){
    
    # aggregate into 6-month precip
    slider::slide_dbl(x,
                      mean, # could add other stats like range to account for instability
                      .complete = T,
                      .before = 5) # window size
    
  },
  FUTURE = T,
  .fname = "time") %>% 
  
  st_as_stars(proxy = F) %>% 
  aperm(c(2,3,1)) %>% 
  st_set_dimensions("time", values = dates) %>% 
  setNames("tasmean")


tasmean_predictor <- 
  tasmean_predictor %>% 
  st_set_dimensions("time", values = dates + months(lead_time + 1))


# extract from points
tb_pred_tasmean <- 
  aggregate(tasmean_predictor,
            sample_pts,
            function(x) x[1],
            as_points = FALSE) %>%
  as_tibble() %>%
  mutate(as_tibble(st_coordinates(geometry))) %>% 
  select(-geometry)




# PROCESS SEA SURFACE TEMP ----------------------------------------------------

# copy_cat("mean-sst")

# read files
sst_raw <- 
  str_glue("{dir_raw_data}/mean-sst.nc") %>% 
  read_ncdf(proxy = T)

downsampled <-
  sst_raw %>%
  st_as_stars(downsample = c(3,3,0), proxy = F)

## EOF ********************************

{
  # # 5000 sample pts for EOF
  # set.seed(1)
  # rnd_pts <- 
  #   st_sample(land_pol %>% filter(a == "ocean"), 
  #             5000)
  # 
  # 
  # # aggregate and detrend
  # 
  # d <- year(dates)
  # 
  # s_detrended <- 
  #   downsampled %>% 
  #   st_apply(c(1,2), function(x){
  #     
  #     if(any(is.na(x))){
  #       rep(NA, length(unique(d)))
  #     } else {
  #       x_annual_means <- aggregate(x ~ d, FUN = mean)$x
  #       x_annual_means - mean(x_annual_means)
  #     }
  #     
  #   },
  #   FUTURE = T,
  #   .fname = "time") %>% 
  #   st_set_dimensions("time",
  #                     values = d %>% unique() %>% sort())
  # 
  # 
  # # calculate EOFs
  # l_eof <- eof(s_detrended, rnd_pts)
  # 
  # l_eof$eof_model$sdev^2 -> eigs
  # (eigs/sum(eigs))[1:10] # proportion var explained
  # 
  # 
  # # regionalize with K-means
  # eof_rast <- 
  #   l_eof$s_eof %>% 
  #   select(1:2) # first 2 components
  # 
  # tb <- 
  #   eof_rast %>%
  #   as_tibble() %>%
  #   mutate(across(starts_with("PC"), 
  #                 ~ifelse(is.na(.x), -999, .x)))
  # 
  # set.seed(1)
  # km <- 
  #   kmeans(tb[,-c(1:2)], 
  #          centers = 9, 
  #          iter.max = 1000, 
  #          nstart = 10, 
  #          algorithm = "Lloyd")
  # 
  # km_rast <- 
  #   km$cluster %>%
  #   matrix(dim(eof_rast[[1]])) %>% 
  #   st_as_stars()
  # 
  # st_dimensions(km_rast) <- st_dimensions(eof_rast)
  # 
  # # remove land
  # km_rast <- 
  #   c(km_rast, eof_rast) %>% 
  #   mutate(A1 = ifelse(is.na(PC1), NA, A1)) %>% 
  #   select(A1)
  # 
  # # generate points representative of clusters
  # km_pts <-
  #   km_rast %>%
  #   slice(longitude, 2:360) %>% 
  #   slice(latitude, 2:180) %>% 
  #   st_as_sf(as_points = F, merge = T) %>% # polygonize
  #   mutate(area = st_area(.) %>% set_units(km^2)) %>%
  #   group_by(A1) %>% 
  #   filter(area %in% sort(area, decreasing = T)[c(1:3)]) %>% # 3 largest
  #   st_point_on_surface() %>% 
  #   select(A1)
  # 
  # st_write(km_pts, "output/station_pts_sst.gpkg")
}


## Prepare predictor ******************

km_pts <- 
  "output/station_pts_sst.gpkg" %>% st_read()

# extract SST from points
tb_pts <- 
  aggregate(downsampled,
            km_pts,
            function(x) x[1],
            as_points = FALSE) %>%
  as_tibble() %>%
  mutate(as_tibble(st_coordinates(geometry))) %>% 
  select(-geometry) %>% 
  
  mutate(time = as_date(time)) %>% 
  
  group_by(X, Y) %>%
  mutate(id = cur_group_id()) %>% 
  drop_units()
  

# rolling window
tb_pred_sst <- 
  tb_pts %>%
  mutate(sst = slider::slide_dbl(sst,
                                 mean,
                                 .complete = T,
                                 .before = 5)) %>% arrange(id) %>%
  
  # anomalies
  # group_by(mth = month(time), id) %>%
  # mutate(sst = sst - mean(sst, na.rm = T)) %>%
  # select(-mth) %>% 
  
  ungroup() %>% 
  select(-X, -Y) %>%
  pivot_wider(names_from = id, values_from = sst, names_prefix = "sst_")

# adjust date
tb_pred_sst <- 
  tb_pred_sst %>% 
  mutate(time = time + months(lead_time + 1))


# PROCESS GEOPOTENTIAL 500 ----------------------------------------------------

# copy_cat("mean-geopotential-500")

# read files
gpot500_raw <- 
  str_glue("{dir_raw_data}/mean-geopotential-500.nc") %>% 
  read_ncdf(proxy = T)

downsampled <-
  gpot500_raw %>%
  st_as_stars(downsample = c(3,3,0), proxy = F)



## EOF ********************************

{
  
  # # 5000 sample pts for EOF
  # set.seed(1)
  # rnd_pts <- 
  #   st_sample(land_pol, 
  #             5000)
  # 
  # 
  # # aggregate and detrend
  # 
  # 
  # d <- year(dates)
  # 
  # s_detrended <- 
  #   downsampled %>% 
  #   st_apply(c(1,2), function(x){
  #     
  #     if(any(is.na(x))){
  #       rep(NA, length(unique(d)))
  #     } else {
  #       x_annual_means <- aggregate(x ~ d, FUN = mean)$x
  #       x_annual_means - mean(x_annual_means)
  #     }
  #     
  #   },
  #   FUTURE = T,
  #   .fname = "time") %>% 
  #   st_set_dimensions("time",
  #                     values = d %>% unique() %>% sort())
  # 
  # 
  # # calculate EOFs
  # l_eof <- eof(s_detrended, rnd_pts)
  # 
  # l_eof$eof_model$sdev^2 -> eigs
  # (eigs/sum(eigs))[1:10] # proportion var explained
  # 
  # 
  # # regionalize with K-means
  # eof_rast <- 
  #   l_eof$s_eof %>% 
  #   select(1:3) # first 3 components
  # 
  # tb <- 
  #   eof_rast %>%
  #   as_tibble()
  # 
  # set.seed(1)
  # km <- 
  #   kmeans(tb[,-c(1:2)], 
  #          centers = 9, 
  #          iter.max = 1000, 
  #          nstart = 10, 
  #          algorithm = "Lloyd")
  # 
  # km_rast <- 
  #   km$cluster %>%
  #   matrix(dim(eof_rast[[1]])) %>% 
  #   st_as_stars()
  # 
  # st_dimensions(km_rast) <- st_dimensions(eof_rast)
  # 
  # 
  # # generate points representative of clusters
  # km_pol <- 
  #   km_rast %>%
  #   slice(longitude, 2:360) %>% 
  #   slice(latitude, 2:180) %>% 
  #   st_as_sf(as_points = F, merge = T) %>% 
  #   mutate(area = st_area(.) %>% set_units(km^2))
  # 
  # km_pol_large <- 
  #   km_pol %>% 
  #   filter(area > set_units(25e6, km^2))
  # 
  # km_pts_large <- 
  #   tibble(xmin = seq(0, 360, length.out = 4)[-4],
  #          xmax = seq(0, 360, length.out = 4)[-1]) %>% 
  #   
  #   pmap(function(xmin, xmax){
  #     
  #     km_pol_large %>% 
  #       st_crop(xmin = xmin,
  #               xmax = xmax,
  #               ymin = -90,
  #               ymax = 90) %>% 
  #       st_point_on_surface()
  #     
  #   }) %>% 
  #   bind_rows()
  # 
  # km_pts <-
  #   km_pol %>% 
  #   filter(area < set_units(25e6, km^2),
  #          area > set_units(45e4, km^2)) %>%
  #   group_by(A1) %>% 
  #   filter(area %in% sort(area, decreasing = T)[c(1:4)]) %>%
  #   st_point_on_surface() %>% 
  #   bind_rows(km_pts_large) %>% 
  #   select(A1)
  # 
  # st_write(km_pts, "output/station_pts_geopot500.gpkg")
  
}


## Prepare predictor ******************

km_pts <- 
  "output/station_pts_geopot500.gpkg" %>% st_read()

# extract geopot from points
tb_pts <- 
  aggregate(downsampled,
            km_pts,
            function(x) x[1],
            as_points = FALSE) %>%
  as_tibble() %>%
  mutate(as_tibble(st_coordinates(geometry))) %>% 
  select(-geometry) %>% 
  
  mutate(time = as_date(time)) %>% 
  
  group_by(X, Y) %>%
  mutate(id = cur_group_id()) %>% 
  drop_units()


# rolling window
tb_pred_gpot500 <- 
  tb_pts %>%
  rename(gpot500 = z) %>% 
  mutate(gpot500 = slider::slide_dbl(gpot500,
                                     mean,
                                     .complete = T,
                                     .before = 5)) %>% arrange(id) %>%
  
  # anomalies
  # group_by(mth = month(time), id) %>%
  # mutate(gpot500 = gpot500 - mean(gpot500, na.rm = T)) %>%
  # select(-mth) %>% 
  
  ungroup() %>%
  select(-X, -Y) %>%
  pivot_wider(names_from = id, values_from = gpot500, names_prefix = "gpot500_")

# adjust date
tb_pred_gpot500 <- 
  tb_pred_gpot500 %>% 
  mutate(time = time + months(lead_time + 1))



# PROCESS ELEVATION -----------------------------------------------------------

# read file
elev <- 
  "/mnt/bucket_mine/era/era5_geopotential.nc" %>% 
  read_ncdf(proxy = F) %>% 
  # .[, 930:1100, 130:300, ] %>%
  adrop() %>% 
  drop_units()

# set ocean to 0
# elev <- 
#   st_warp(land, elev) %>% 
#   c(elev) %>% 
#   mutate(z = ifelse(a == "ocean" | z < 0, 0, z)) %>% 
#   select(z)

# extract from points
tb_pred_elev <- 
  aggregate(elev,
            sample_pts,
            function(x) x[1],
            as_points = FALSE) %>%
  as_tibble() %>%
  mutate(as_tibble(st_coordinates(geometry))) %>% 
  select(-geometry)



# aspect of slope
terrain <- 
  elev %>% 
  as("SpatRaster") %>% 
  terra::terrain(c("aspect", "slope"), neighbors = 8)

tb_pred_aspect <- 
  terrain %>% 
  st_as_stars() %>%
  # .[, 930:1100, 130:300] %>%
  split("band") %>% 
  mutate(aspect = case_when(is.na(aspect) ~ NA_character_,
                            slope < 3 ~ "flat",
                            aspect > 0 & aspect <= 22 ~ "E",
                            aspect > 22 & aspect <= 67 ~ "SE",
                            aspect > 67 & aspect <= 112 ~ "S",
                            aspect > 112 & aspect <= 157 ~ "SW",
                            aspect > 157 & aspect <= 202 ~ "W",
                            aspect > 202 & aspect <= 247 ~ "NW",
                            aspect > 247 & aspect <= 292 ~ "N",
                            aspect > 292 & aspect <= 337 ~ "NE",
                            aspect > 337 ~ "E"
  )) %>% 
  select(aspect) %>% #plot()
  
  aggregate(sample_pts,
            function(x) x[1],
            as_points = FALSE) %>%
  as_tibble() %>%
  mutate(as_tibble(st_coordinates(geometry))) %>% 
  select(-geometry) %>% 
  mutate(aspect = factor(aspect))




# ASSEMBLE TABLE --------------------------------------------------------------

tb_f <- 
  tb_reponse %>%
  left_join(tb_pred_elev, by = c("X", "Y")) %>% 
  left_join(tb_pred_aspect, by = c("X", "Y")) %>% 
  left_join(tb_pred_precip, by = c("X", "Y", "time")) %>% 
  left_join(tb_pred_tasmean, by = c("X", "Y", "time")) %>% 
  left_join(tb_pred_sst, by = "time") %>% 
  left_join(tb_pred_gpot500, by = "time") %>% 
  na.omit()
  

# tb_f %>% 
#   filter(time >= "1959-06-01") %>% 
#   filter(time <= "2021-03-01")

tb_f %>% 
  write_csv("output/tb_f.csv")



