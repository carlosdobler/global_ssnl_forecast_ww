
# Function to obtain table of space-invariant predictors
sp_inv_pred_tb <- function(s, tble){
  
  tble %>%
    future_pmap_dfc(function(pos_lon, pos_lat, r, ...){
      
      s[,
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
              
              avg6mths = mean(w)
              
            ) -> tb
            
            tb %>% 
              rename_with(~str_glue("{names(s)}_{names(tb)}_{r}"))
          })
      }
    }) %>% 
    
    mutate(time = dates_sub,
           mth = month.abb[month(dates_sub)] %>% factor(),
           .before = 1)
}


# *****************


# Function to visualize ERA slices (mapview)
view_era <- function(s){
  s %>% 
    slice(longitude, 2:1439) %>% 
    slice(latitude, 2:720) %>% 
    mapview::mapview()
}


# *****************


# Function moving window: mean
# foc <- function(x, w) {
#   raster::as.matrix(raster::focal(raster::raster(x), w, mean, na.rm = T, pad = T))
# }


# *****************


# Function to obtain pca
pca_maps <- function(detrended_s, pts){
  
  # extract values 
  print(str_glue("Extracting values ..."))
  
  aggregate(detrended_s, pts, function(x) x[1], as_points = F) %>% 
    st_as_sf() -> pts_val
  
  bind_cols(
    pts_val %>% 
      st_coordinates() %>% 
      as_tibble(),
    
    pts_val %>% 
      st_drop_geometry() %>% 
      as_tibble() %>% 
      rename_with(~ st_get_dimension_values(detrended_s, "time") %>% 
                    {str_c("x_", year(.), "_", month(.))})
  ) -> pts_val
  
  # fit pca
  pts_val %>%
    select(-c(X,Y)) %>% 
    na.omit() %>% 
    prcomp(scale. = T) -> pca_model
  
  # predict pca model to the rest of data
  print(str_glue("Predicting ..."))
  
  detrended_s %>%
    st_set_dimensions("time", 
                      values = st_get_dimension_values(detrended_s, "time") %>% 
                        {str_c("x_", year(.), "_", month(.))}) %>% 
    split("time") %>% 
    predict(pca_model) -> pca_rast
  
  list(pca_model = pca_model,
       pca_rast = pca_rast)
  
}


# *****************


# Function to assemble tb of non-local predictors
assemble_tb_nonlocal <- function(s, pts){
  
  pts %>% 
    st_coordinates() %>% 
    as_tibble() %>% 
    rowwise() %>%
    mutate(pos_lon = which.min(abs(vect_lon - X)),
           pos_lat = which.min(abs(vect_lat - Y))) %>% 
    ungroup() %>%
    
    mutate(r = row_number() %>% 
             str_pad(2, "left", "0")) %>% 
    
    pmap_dfc(function(pos_lon, pos_lat, r, ...){
      
      print(str_glue("{r} / {nrow(pts)}"))
      
      s[,
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
              rename_with(~str_glue("{names(s)}_{names(tb)}_{r}"))
          })
      }
    }) -> tb_pred
  
  tb_pred %>% 
    mutate(time = dates_sub,
           .before = 1) -> tb_pred
  
}






