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

# Calling data ----------------------------------------------------
data   <- read.csv("data/data.csv", header = TRUE, sep = ";", dec=",")
traits <- read.csv("data/traits.csv", header = TRUE, sep = ";", dec=",")

# Classification the positive/negative records --------------------  
data$captured <- NA
data$captured <- ifelse (data$tot_biomass_kg >0, "yes", data$captured)
data$captured <- ifelse (data$tot_biomass_kg ==0, "no", data$captured)

# Comparing sampling efforts --------------------------------------
library (plyr)
sampling_efforts <- data[,c("year","Id_infraction","captured")]
sampling_efforts <- droplevels(sampling_efforts[!duplicated(sampling_efforts), ])
sampling_efforts <- ddply (sampling_efforts,. (year, captured), summarise, 
                           total=length(Id_infraction))

sampling_efforts$time_covid <- NA
sampling_efforts$time_covid <- ifelse (sampling_efforts$year<=2019, "Before", sampling_efforts$time_covid)
sampling_efforts$time_covid <- ifelse (sampling_efforts$year>2019, "During", sampling_efforts$time_covid)
sampling_efforts

t.test(total ~ time_covid, data = sampling_efforts)

library (ggplot2)
ggplot(sampling_efforts, aes(x = time_covid, y = total, fill = captured, group = captured)) +
  geom_bar(stat = "identity", position = "stack")

# Estimating records distances from the main island ---------------
library (raster)
pts_coordenates <- ddply (data,. (year, Id_infraction, mpa_area),
                          summarise, lat=unique(lat), lon=unique(lon))

Dist_Alcatrazes <- pointDistance (c(-45.693050, -24.102211), 
                                  cbind (lon=pts_coordenates$lon,
                                         lat=pts_coordenates$lat), 
                                  lonlat=TRUE)/1000

Coord_records <- data.frame (Year=pts_coordenates$year, 
                             Id_infraction = pts_coordenates$Id_infraction, 
                             lon=pts_coordenates$lon,
                             lat=pts_coordenates$lat, 
                             Distance=unlist(Dist_Alcatrazes))

# Records Map ------------------------------------------
library (rgdal)
alcatrazes_map <- readOGR(dsn="shapefiles/islands", layer = "ilhas_e_parceis"); crs (alcatrazes_map)
refugia_lim <- readOGR(dsn="shapefiles/limits", layer = "refugio_limite_conforme_decreto")
esec_lim <- readOGR(dsn="shapefiles/limits", layer = "esec__limite_conforme_decreto")

# CRS Projection
crs_projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
coord <-  pts_coordenates
coordenates_loc <- cbind(coord$lon, coord$lat)

prj.coord <- project(coordenates_loc, proj = crs_projection)
coord <- cbind (prj.coord, coord)
names (coord)[1:2] <- c("X.prj", "Y.prj")
coord$COVID19 <- NA
coord$COVID19 <- ifelse(coord$year <=2019, "Before", coord$COVID19)
coord$COVID19 <- ifelse(coord$year >2019, "Current", coord$COVID19)

alcatrazes_map_proj <- spTransform(alcatrazes_map, CRSobj = crs_projection) 
refugia_lim_proj <- spTransform(refugia_lim, CRSobj = crs_projection) 
esec_lim_proj <- spTransform(esec_lim, CRSobj = crs_projection) 

library("RColorBrewer")
vessel_type <- ddply (data,. (Id_infraction, vessel), summarise, presence=sum(presence))
vessel_type$vessel <- gsub("recreationalboat", "Amateur", vessel_type$vessel)
vessel_type$vessel <- gsub("fishingvessel", "Professional", vessel_type$vessel)
vessel_type <- vessel_type [,1:2]

coord <- left_join(coord, vessel_type)
rm (vessel_type)

library (ggspatial)
Alcatrazes_map <- ggplot () +
  geom_polygon(data=refugia_lim_proj, aes(long,lat, group=group), color="black", fill="#8C96C6", size=0.25, linetype=2) +
  geom_polygon(data=esec_lim_proj, aes(long,lat, group=group), color="black", fill="#FEB24C", size=0.25) +
  lims(x=c(-45.9, -45.2), y=c(-24.9,-24)) +
  labs(x="Longitude",y="Latitude")+
  geom_point(inherit.aes = FALSE, 
             data = as.data.frame(coord),
             mapping = aes(x = X.prj, 
                           y = Y.prj, fill = vessel,
                           shape = factor(COVID19)),
             size = 2.5)+
  scale_fill_manual(name = "vessel", values=c("#01579B", "#FEB24C"),
                    breaks = c("Amateur","Professional")) +
  scale_shape_manual(values = c("Before" = 21, "Current" = 23)) +
  geom_polygon(data=alcatrazes_map_proj, aes(long,lat, group=group), color="black", fill="darkgreen", size=0.25) +
  coord_equal() +
  annotation_scale()+
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.03, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering)+
  theme_light()+
  theme (axis.text.y   = element_text(size=11, angle=0, family = "sans"),
         axis.title.y  = element_text(size=15, angle=90, family = "sans"),
         axis.text.x   = element_text(size=11, angle=0, family = "sans"),
         axis.title.x  = element_text(size=15, angle=0, family = "sans"),
         legend.position = "none", 
         legend.title = element_text(colour="black", size=10, face="bold"),
         legend.background = element_blank(),
         panel.background = element_rect(fill = "aliceblue"),
         panel.grid.major = element_line(color = gray(0.5), 
                                         linetype = "dashed", size = 0.05))
Alcatrazes_map

library(geobr)
library (sf)

states <- read_state(code_state = "all", year=2019)
spaulo <- read_state(code_state = "SP", year = 2019)

# South America map
library (maps)
south_amer_map <- borders("world", regions = c("Brazil", "Uruguay", "Argentina", "French Guiana", "Suriname", "Colombia", "Venezuela",
                                               "Bolivia", "Ecuador", "Chile", "Paraguay", "Peru", "Guyana", "Panama"), 
                          fill = "grey70", colour = "grey70")

SAmerica <-  ggplot() + 
  xlim(-90,-30)+
  south_amer_map + 
  geom_sf(data=states, fill="#2D3E50", color="#2D3E50", size=.15, 
          show.legend = FALSE)+
  geom_sf(data = spaulo, fill=NA, color="#FEBF57", size=.2) +
  theme_classic2() + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y   = element_text(size=11, angle=0, family = "sans"),
        axis.title.y  = element_blank(),
        axis.text.x   = element_text(size=11, angle=0, family = "sans"),
        axis.title.x  = element_blank())


# Models ---------------------------------------------------------
library (dplyr)
values_metrics <- ddply (data,. (year, Id_infraction, vessel, mpa_area), summarise,
                         richness=length(unique (species)), biomass=sum(tot_biomass_kg))
values_metrics <- left_join(Coord_records, values_metrics)

values_metrics$richness <- ifelse (values_metrics$biomass==0, 0, values_metrics$richness)

library (brms)
values_metrics$COVID19 <- NA
values_metrics$COVID19 <- ifelse (values_metrics$year <=2019, "before", values_metrics$COVID19)
values_metrics$COVID19 <- ifelse (values_metrics$year >2019, "current", values_metrics$COVID19)
values_metrics$COVID19 <- as.factor(values_metrics$COVID19)

# Species richness model
fit_zinb <- brm (bf(richness ~ COVID19 * vessel + Distance, 
                    zi~ COVID19*vessel), 
                 data = values_metrics, 
                 family = zero_inflated_negbinomial("log"))
summary(fit_zinb)
plot(fit_zinb, ask = FALSE)

data_summary_richness <- values_metrics %>% group_by(COVID19, mpa_area, vessel) %>%
  summarize(mean.value = mean(richness), se.value = plotrix::std.error(richness))


# Biomass model
hist (values_metrics$biomass)
fit_hugamma <- brm (bf(biomass ~ COVID19 * vessel + Distance, 
                       hu ~ COVID19 * vessel), data = values_metrics,
                    family = hurdle_gamma("log"))
summary (fit_hugamma)
plot(fit_hugamma, ask = FALSE)

data_summary_biomass <- values_metrics %>% group_by(COVID19, mpa_area, vessel) %>%
  summarize(mean.value = mean(biomass), se.value = plotrix::std.error(biomass))

# Figure S1
Fig_S1a <- ggplot(data_summary_richness, aes(x = COVID19, y = mean.value, group = mpa_area, col = mpa_area)) +
  geom_line() +
  annotate("rect", xmin=1.5, xmax=Inf, ymin = -Inf, ymax=Inf, alpha = .2)+
  geom_errorbar(aes(ymax = mean.value + se.value, ymin = mean.value - se.value), width = 0.1) +
  geom_point() +
  labs (x="Period", y="Species richness")+
  facet_grid(vessel~ .) +
  scale_x_discrete(labels=c("Before COVID-19", "During COVID-19"))+
  theme_light()+
  theme (axis.text.y   = element_text(size=11, angle=0, family = "sans"),
         axis.title.y  = element_text(size=15, angle=90, family = "sans"),
         axis.text.x   = element_text(size=11, angle=0, family = "sans"),
         axis.title.x  = element_text(size=15, angle=0, family = "sans"),
         legend.title = element_text(colour="black", size=10, face="bold"),
         legend.background = element_blank())


Fig_S1b <- ggplot(data_summary_biomass, aes(x = COVID19, y = mean.value, group = mpa_area, col = mpa_area)) +
  annotate("rect", xmin=1.5, xmax=Inf, ymin = -Inf, ymax=Inf, alpha = .2)+
  geom_line() +
  geom_errorbar(aes(ymax = mean.value + se.value, ymin = mean.value - se.value), width = 0.1) +
  geom_point() +
  labs (x="Period", y="Biomass (Kg)")+
  facet_grid(vessel~ .) +
  scale_x_discrete(labels=c("Before COVID-19", "During COVID-19"))+
  theme_light()+
  theme (axis.text.y   = element_text(size=11, angle=0, family = "sans"),
         axis.title.y  = element_text(size=15, angle=90, family = "sans"),
         axis.text.x   = element_text(size=11, angle=0, family = "sans"),
         axis.title.x  = element_text(size=15, angle=0, family = "sans"),
         legend.position = "none", 
         legend.title = element_text(colour="black", size=10, face="bold"),
         legend.background = element_blank())

library (ggpubr)
ggarrange(Fig_S1a, Fig_S1b, ncol = 2, nrow = 1)


# Functional space -----------------------------------------------
## This functions was extracted from the Maire et al. (2015) - GEB
source ("R/quality_funct_space_fromdist2.R")

mobility  <- sapply(traits$mobility, function(x) {if (x=="sed"){1} else if (x=="mob"){2} else if (x=="vmob"){3}})
schooling <- sapply(traits$size_group,  function(x) {if (x=="sol"){1} else if (x=="pair"){2} else if (x=="smallg"){3} else if (x=="medg"){4} else if (x=="largeg"){5}})
position  <- sapply(traits$position, function(x) {if (x=="benthic"){1} else if (x=="pelagic"){2}})

mobility  <-ordered (mobility); length(mobility)
position  <-ordered (position); length(position)
schooling <-ordered (schooling);length(schooling)

traits_ordered <-data.frame (name=unique(traits$name), 
                             size=traits$body_mass, 
                             mobility=mobility,
                             size_group=schooling, 
                             position=position, 
                             diet=traits$diet, 
                             thermoregulation=traits$thermoregulation)

rownames(traits_ordered) <- traits_ordered$name
traits_ordered           <- traits_ordered[,-1]

rm (mobility, position, schooling)

# Gower Distance
library (cluster)
library (ade4)
gower_matrix <- daisy (data.matrix(traits_ordered), metric=c("gower")) # function daisy compute all dissimilarities distances
pcoa <- dudi.pco (quasieuclid(gower_matrix), scannf=F, nf=4) # here 10 axes are retained
barplot(pcoa$eig)

quality <-quality_funct_space_fromdist (gower_matrix,  nbdim=6,   plot=NA) 
quality$meanSD # the minimal value corresponds to the best space to use, here 6 axes 
Inertia_4axis <- (sum(pcoa$eig[1:4])) /(sum(pcoa$eig)) # percentage of inertia explained by the 4 first axes = 89.4%

# Building the functional space
pcoa1 <- data.frame (pcoa$li)
pcoa1$name   <- rownames(pcoa1)
rownames(pcoa1) <- NULL
pcoa1$name <- gsub (" ", "_", pcoa1$name)

data$presence <- 1
db <- ddply (data,. (year, vessel, species), summarise, 
             presence=sum(presence))
colnames(db)[3] <- "name"

db <- droplevels(db[!db$name=="nofish",])
db <- left_join(db, pcoa1)

db$COVID19 <- NA
db$COVID19 <- ifelse (db$year <=2019, "Before", db$COVID19)
db$COVID19 <- ifelse (db$year >2019, "During", db$COVID19)

# Estimating the total functional volume captured
hull_db <- chull(db[,5:6]) # Using the Axis 1 and 2 from the PCoA
hull_db <- c (hull_db, hull_db[1])

library(geometry)
hull <- convhulln(db[,5:6], options = "FA")
hull$vol

# Split records between professional/amateur before COVID-19 pandemic
# Before
Before_COVID19 <- droplevels (db[!db$COVID19=="During", ])
convhulln(Before_COVID19[,5:6], options = "FA")$vol
hull_Before_COVID19 <- chull(Before_COVID19[,5:6])
hull_Before_COVID19 <- c(hull_Before_COVID19, hull_Before_COVID19[1])
# Before - Professional
Before_COVID19_Professional <- droplevels (subset (Before_COVID19, vessel=="fishingvessel"))
convhulln(Before_COVID19_Professional[,5:6], options = "FA")$vol
hull_Before_COVID19_Professional <- chull(Before_COVID19_Professional[,5:6])
hull_Before_COVID19_Professional <- c(hull_Before_COVID19_Professional, hull_Before_COVID19_Professional[1])
# Before - Amateur
Before_COVID19_Amateur <- droplevels (subset (Before_COVID19, vessel=="recreationalboat"))
convhulln(Before_COVID19_Amateur[,5:6], options = "FA")$vol
hull_Before_COVID19_Amateur <- chull(Before_COVID19_Amateur[,5:6])
hull_Before_COVID19_Amateur <- c(hull_Before_COVID19_Amateur, hull_Before_COVID19_Amateur[1])

# Split records between professional/amateur during COVID-19 pandemic
# During
During_COVID19 <- droplevels (db[!db$COVID19=="Before", ])
convhulln(During_COVID19[,5:6], options = "FA")$vol
hull_During_COVID19 <- chull(During_COVID19[,5:6])
hull_During_COVID19 <- c(hull_During_COVID19, hull_During_COVID19[1])
# During - Professional
During_COVID19_Professional <- droplevels (subset (During_COVID19, vessel=="fishingvessel"))
convhulln(During_COVID19_Professional[,5:6], options = "FA")$vol
hull_During_COVID19_Professional <-chull(During_COVID19_Professional[,5:6])
hull_During_COVID19_Professional <- c(hull_During_COVID19_Professional, hull_During_COVID19_Professional[1])
# During - Amateur
During_COVID19_Amateur <- droplevels (subset (During_COVID19, vessel=="recreationalboat"))
convhulln(During_COVID19_Amateur[,5:6], options = "FA")$vol
hull_During_COVID19_Amateur <- chull(During_COVID19_Amateur[,5:6])
hull_During_COVID19_Amateur <- c(hull_During_COVID19_Amateur, hull_During_COVID19_Amateur[1])

# Figure Functional Space using the Axis 1 and 2 from the PCoA
FSpace_Before_Fig1b <- ggplot(db, aes(x = A1, y = A2)) +
  labs(x="PCoA1", y="PCoA2") +
  theme_classic2(base_family = "sans")+ 
  theme(legend.position = "top", 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=15, angle=90, family = "sans"),
        plot.title   = element_text(size = 19, face = "bold"))+ 
  scale_size_area()+
  geom_polygon(data=db[hull_db,], fill="#CDC7AE", alpha = 2/10, colour="black", linetype = 1) +
  geom_polygon(data=Before_COVID19[hull_Before_COVID19,], fill="grey", alpha = 6/10, colour="black", linetype = 3) +
  geom_polygon(data=Before_COVID19_Professional[hull_Before_COVID19_Professional,], fill="#E69F00", alpha = .7, colour="black", 
               linetype = 1) +
  geom_point(data = Before_COVID19_Professional, aes(x = A1, y = A2), size=6, pch=21, col="black", bg = "#E69F00") +
  #geom_text (data=Before_COVID19_Professional[hull_Before_COVID19_Professional,], aes(label=name), hjust=0, vjust=0) +
  geom_polygon(data=Before_COVID19_Amateur[hull_Before_COVID19_Amateur,], fill="#56B4E9", alpha = .7, colour="black", 
               linetype = 1)+
  geom_point(data = Before_COVID19_Amateur, aes(x = A1, y = A2), size=6, pch=21, col="black", bg = "#56B4E9")
  #geom_text (data=Before_COVID19_Amateur[hull_Before_COVID19_Amateur,], aes(label=code), hjust=0, vjust=0) 

FSpace_Before_Fig1c <-ggplot(db, aes(x = A1, y = A2)) +
  theme_classic2(base_family = "sans")+ 
  labs(x="PCoA1", y="PCoA2") +
  theme(legend.position = "top", 
        axis.title.x = element_text(size=15, angle=0, family = "sans"),
        axis.title.y = element_text(size=15, angle=90, family = "sans"),
        plot.title   = element_text(size = 19, face = "bold"))+ 
  scale_size_area()+
  geom_polygon(data=db[hull_db,], fill="#CDC7AE", alpha = 2/10, colour="black", linetype = 1) +
  geom_polygon(data=During_COVID19[hull_During_COVID19,], fill="grey", alpha = 6/10, colour="black", linetype = 3) +
  geom_polygon(data=During_COVID19_Professional[hull_During_COVID19_Professional,], fill="#E69F00", alpha = .7, colour="black", 
               linetype = 1) +
  geom_point(data = During_COVID19_Professional, aes(x = A1, y = A2), size=6, pch=21, col="black", bg = "#E69F00") +
  #geom_text (data=During_COVID19_Professional[hull_During_COVID19_Professional,], aes(label=name), hjust=0, vjust=0) +
  geom_polygon(data=During_COVID19_Amateur[hull_During_COVID19_Amateur,], fill="#56B4E9", alpha = .7, colour="black", 
               linetype = 1)+
  geom_point(data = During_COVID19_Amateur, aes(x = A1, y = A2), size=6, pch=21, col="black", bg = "#56B4E9")
  #geom_text (data=During_COVID19_Amateur[hull_During_COVID19_Amateur,], aes(label=code), hjust=0, vjust=0) 

ggarrange(FSpace_Before_Fig1b, FSpace_Before_Fig1c, 
          ncol = 1, nrow = 2,
          labels = c("(b)","(c)"),align = "hv")

