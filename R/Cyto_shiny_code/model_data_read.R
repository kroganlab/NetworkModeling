
library(data.table); library(reshape2)

load("data/mouse_data_large2")
Mouse_map <- read.csv(file = "data/basemap_menche_mouse.csv")
mygene_mouse_1to1 <- read.csv(file = "data/mygene_mouse_1to1.csv")