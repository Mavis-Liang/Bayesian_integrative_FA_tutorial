install.packages("devtools") 
library(devtools)
devtools::install_github("AleAviP/BFR.BE")
library(BFR.BE)

devtools::install_github("Mavis-Liang/MSFA")# Mavis's version for simulation
library(MSFA)
# DeVito's version
# install.packages("remotes")
# remotes::install_github("rdevito/MSFA")


# To use PFA: download the file FBPFA-PFA.R and PFA.cpp from the repository 
# and store them in our working directory
# install.packages("expm") package required in their script
#source("./FBPFA-PFA.R")


# On linux, you might need to install extra dependencies like PROJ, sqlite3 and GDAL onto PATH.
# On Windows, I need to do several updates, particularly for terra. Some trouble shoots are needed.
# We found that the "genedata.rda" is not actually rda format, insteat it's a file that point to a larger file. 
# I remove the file for compilation since it's outputing different errors.
devtools::install_github("noirritchandra/SUFA", build_vignettes = TRUE)
library(SUFA)
