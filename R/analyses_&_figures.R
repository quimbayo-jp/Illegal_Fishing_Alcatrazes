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
