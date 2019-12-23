###################################################
### Martone Lab Rocky Intertidal Community Data ###
###     All Data collected on Calvert Island    ###
###################################################
# 
# by Matt Whalen
# updated 20 December 2018

# This script runs the combine script, cleaning script, and generates a species list 

## COMBINE
source( "Data/R code for Data Prep/Martone_Hakai_Intertidal_Data_1_CombineAll.R" )

## CLEAN
source( "Data/R code for Data Prep/Martone_Hakai_Intertidal_Data_2_Clean.R" )

## SPECIES LIST
source( "Data/R code for Data Prep/Martone_Hakai_Intertidal_Data_3_SpeciesList.R" )