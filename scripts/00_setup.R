

# PROJECT-WIDE SET-UP

# Load libraries
library(tidyverse)
library(lubridate)
library(tictoc)
library(furrr)
library(stars)
library(units)


# Allow multicore
options(future.fork.enable = T)


# Key directories
dir_raw_data <- "/mnt/pers_disk/raw_data"
