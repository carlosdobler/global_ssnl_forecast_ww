
# SETUP -----------------------------------------------------------------------

source("scripts/00_setup.R")
library(tidymodels)

plan(multicore)
# sf_use_s2(F)


# Read table
tb <- 
  "output/tb_f.csv" %>% 
  read_csv() %>% 
  mutate(time = as_date(time))

tb_f <- 
  tb %>% 
  mutate(yr = year(time)) %>% 
  filter(month(time) == 8)



# split train/test by group
set.seed(1)
data_split <- group_initial_split(tb_f, group = yr, prop = 0.8)

tb_train <- training(data_split)
tb_test <- testing(data_split)


# specify model
model_spec <- 
  rand_forest(trees = 1500) %>% 
  set_engine("ranger", num.threads = 8) %>% 
  set_mode("regression")


# recipe
model_recip <- 
  recipe(response ~ ., data = tb_f) %>% 
  update_role(time, yr, new_role = "id_variable") %>% 
  step_date(time, features = "month")


# workflow
wflow <- 
  workflow() %>% 
  add_model(model_spec) %>% 
  add_recipe(model_recip)


# fit
model_fit <- 
  wflow %>% 
  fit(data = tb_train)


# predict
model_pred <- 
  tb_train %>% 
  select(X, Y, time, response) %>% 
  bind_cols(predict(model_fit, tb_train))

{
  y <- tb_train$yr %>% unique() %>% sample(1)
  # m <- sample(12, 1)
  
  model_pred %>% 
    filter(year(time) == y) %>% 
    # filter(month(time) == m) %>%
    mutate(deviation = response - .pred,
           deviation = case_when(deviation > 100 ~ 100,
                                 deviation < -100 ~ -100,
                                 TRUE ~ deviation)) %>% 
    ggplot(aes(X, Y, color = deviation)) +
    geom_point() +
    colorspace::scale_color_continuous_diverging() +
    coord_equal() #+
    #labs(subtitle = str_glue("Year: {y} / Month: {month.abb[m]}"))
  
}
