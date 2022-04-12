library(colorspace)
library(randomForest)

rf_model <- readRDS("output/rf_model_sst_2vars_11x50.rds")

tb_f %>% 
  filter(time %in% dates_test) -> tb_test_wh

predict(rf_model, newdata = tb_test_wh[,-1]) -> tp_predict

tb_f %>% 
  filter(time %in% dates_test) %>% 
  mutate(tp_pred = unname(tp_predict)) -> tb_f_p

tb_f_p %>% 
  select(time, longitude, latitude, tp, tp_pred) -> tb_f_p

d <- sample(dates_test, 1)

tb_f_p %>% 
  pivot_longer(starts_with("tp"), names_to = "var", values_to = "mm") %>% 
  
  filter(time == d) %>% 
  
  ggplot(aes(x = longitude, y = latitude, fill = mm)) +
  geom_raster() +
  scale_fill_continuous_sequential("Viridis", rev = F) +
  facet_wrap(~var, ncol = 2) +
  coord_fixed() +
  labs(subtitle = d) +
  theme(axis.title = element_blank()) +
  scale_x_continuous(n.breaks = 3)






predict(rf_model, newdata = tb_valid) -> tp_pred

tb_valid %>% 
  mutate(tp_pred = tp_pred) -> tb_valid


ggplot(tb_valid, aes(x = tp, y = tp_pred)) +
  geom_hex(bins = 35, show.legend = F) +
  geom_abline(linetype = "2222") +
  coord_fixed() +
  scale_fill_continuous_sequential("Plasma", rev = F, trans = "log") +
  labs(x = "observed (mm)",
       y = "predicted (mm)",
       title = "All seasons")


tb_valid %>% 
  filter(mth == "Mar") %>% 
  ggplot(aes(x = tp, y = tp_pred)) +
  geom_hex(bins = 35, show.legend = F) +
  geom_abline(linetype = "2222") +
  coord_fixed() +
  scale_fill_continuous_sequential("Plasma", rev = F, trans = "log") +
  labs(x = "observed (mm)",
       y = "predicted (mm)",
       title = "MAM")
  
