# Use this script when having a .py file
# downloaded from CDS website


library(reticulate)
Sys.unsetenv("RETICULATE_PYTHON")
use_condaenv("r-reticulate", required = T)
source_python("scripts/cds_down_script.py")
