
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


# Function to visualize ERA slices (mapview)
view_era <- function(s){
  s %>% 
    slice(longitude, 2:1439) %>% 
    slice(latitude, 2:720) %>% 
    mapview::mapview()
}

