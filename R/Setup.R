# pipeline variables
start.time <- Sys.time()

# Set working directory to source file location
setwd(dirname(getActiveDocumentContext()$path))

# since moving script from local to github - I want to adjust work dir to be main github dir - therefore 
setwd("..")

working.dir <- getwd()

# For reproducibility
set.seed(42)

# Load custom functions
source("R/Functions/scRNAseq_function.R", local = knitr::knit_global())


###############################
# create output directories
###############################

# Common directories 
if(!dir.exists("Exported_RDS_files")){dir.create("Exported_RDS_files", recursive = T)}

# scRNAseq directories
if(!dir.exists("output")){dir.create("output", recursive = T)}
if(!dir.exists("output/figures")){dir.create("output/figures", recursive = T)}
if(!dir.exists("output/tables")){dir.create("output/tables", recursive = T)}
if(!dir.exists("output/QC")){dir.create("output/QC", recursive = T)}
