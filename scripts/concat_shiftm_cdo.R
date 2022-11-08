# Function to shift central merdidian (ERA5 specific)

func_shift <- function(infile, outfile){
  system(str_glue("cdo sellonlatbox,-180,180,-90.25,90 {infile} {outfile}"))
}


func_cat <- function(infile, outfile){
  
  infile <- str_flatten(infile, " ")
  
  str_glue(
    "cdo cat {infile} {outfile}"
  ) %>%
    system()
  
}





# ***********************


func_cat(list("~/bucket_mine/era/monthly/era5_monthly_mean_daily_precip_1979-2020.nc",
              "~/bucket_mine/era/monthly/era5_monthly_mean_daily_precip_2021.nc"),
         "test.nc")

func_shift(
  "~/bucket_mine/era/monthly/era5_monthly_mean_daily_precip_1979-2021.nc",
  "~/bucket_mine/era/monthly/era5_monthly_mean_daily_precip_1979-2021_shifted.nc"
)

