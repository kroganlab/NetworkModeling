
library(data.table); library(reshape2)

load("Linked_subdirectories/data/mouse_data_large3")
#load("data/DF_reshape")
Mouse_map <- read.csv(file = "Linked_subdirectories/data/basemap_menche_mouse.csv")
Human_map <- read.csv(file = "Linked_subdirectories/data/basemap_menche.csv")
#mygene_mouse_1to1 <- read.csv(file = "Linked_subdirectories/data/mygene_mouse_1to1.csv")