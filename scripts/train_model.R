source("scripts/00_setup.R")
library(randomForest)
plan(multicore, workers = availableCores() - 1)

readRDS("output/tb_f.rds") -> tb_f

tb_f %>%
  select(time:mth, 
         starts_with("mean1"),
         starts_with("mean2")) -> tb_f


# tb_f %>%
#   filter(time == "1981-06-01") %>%
#   ggplot(aes(x = longitude, y = latitude, fill = tp)) +
#   geom_raster() +
#   colorspace::scale_fill_continuous_sequential("Viridis", rev = F) +
#   coord_fixed()


# sample 25% of pixels per date
set.seed(1)
tb_f %>% 
  group_by(time) %>% 
  nest() %>% 
  mutate(data_sub = map(data, function(df){
    df %>%
      slice_sample(prop = 0.25)
  })) %>%
  select(-data) %>%
  unnest(data_sub) %>%
  ungroup() -> tb_f_samplepixels
  

# sample 20% of dates for validation
tb_f %>% 
  pull(time) %>% 
  unique() -> dates_sub

set.seed(1)
sample(dates_sub, length(dates_sub)/5) -> dates_test

tb_f_samplepixels %>% 
  filter(!time %in% dates_test) %>% 
  select(-time) -> tb_train

tb_f_samplepixels %>% 
  filter(time %in% dates_test) %>% 
  select(-time) -> tb_test

# tb_train %>%
#   slice_sample(n = 1000) %>%
#   mutate(asp = as.numeric(asp),
#          mth = mth %>% as.factor() %>% as.numeric()) %>% 
#   rowwise(time) %>% 
#   summarise(nas = sum(is.na(c_across(everything())))) %>% 
#   filter(nas > 0)

tic()
future_map(seq_len(11), function(x){
  randomForest(formula = tp ~ .,
               data = tb_train,
               importance = T,
               ntree = 50)
}) -> rf_model
toc()

# 2.96 hrs 
# mean1term, mean2term
# 50 ntree
# 11 forests

rf_model %>% 
  do.call(combine, .) -> rf_model_comb

saveRDS(rf_model_comb, "output/rf_model_sst_2vars_11x50.rds")

