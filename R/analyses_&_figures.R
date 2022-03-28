## --------------------------------------------------------------------------------------------------------- ##
##  Authors:       Juan P. Quimbayo, with contributions in the code from:                                    ##
##                 Fernanda C. Silva, Jonathan S. Lefcheck, Augusto Flores                                   ##
##  Date:          2021-12-12                                                                                ##
##                                                                                                           ##
##  Notes:         1. This file is intended to provide a guide to the basic                                  ##
##                    project workflow, attempting to 'integrate_analysis' the steps                         ##
##                    necessary to conduct the analyses & figures.                                           ##
## --------------------------------------------------------------------------------------------------------- ##

# clean up
rm(list=ls())

# Calling data
data   <- read.csv("data/data.csv", header = TRUE, sep = ";", dec=",")
traits <- 
