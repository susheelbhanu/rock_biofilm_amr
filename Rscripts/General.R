# Packages
library(tidyverse)
library(dplyr)
library(ggplot2)
library(extrafont)
library(artyfarty)
library(RColorBrewer)
library(splitstackshape)
library(ggalluvial)
library(ggsci)
library(cowplot)

# Figure Theme
theme_set(theme_bw() +
            theme(strip.background = element_rect(colour = "black", fill = "white"), strip.text = element_text(size=17, family="Helvetica")) +
            theme(axis.text.x = element_text(angle = 90, size = 18, family = "Helvetica"),
                  axis.text.y = element_text(size = 18, family = "Helvetica"),
                  axis.title.y = element_text(size = 18, face = "bold", family = "Helvetica"),
                  axis.title.x=element_text(size = 18, face = "bold", family = "Helvetica"),
                  axis.ticks.x=element_blank()) +
            theme(plot.title = element_text(hjust = 0.5)) +
            theme(legend.text = element_text(size = 18, family = "Helvetica"),
                  legend.title = element_text(size = 18, family = "Helvetica", face = "bold"), legend.direction = "vertical", legend.box = "vertical"))



theme_set(theme_scientific() +
            theme(strip.background = element_rect(colour = "black", fill = "white"), strip.text = element_text(size=17, family="Helvetica")) +
            theme(axis.text.x = element_text(angle = 90, size = 18, family = "Helvetica"),
                  axis.text.y = element_text(size = 18, family = "Helvetica"),
                  axis.title.y = element_text(size = 18, face = "bold", family = "Helvetica"),
                  axis.title.x=element_text(size = 18, face = "bold", family = "Helvetica"),
                  axis.ticks.x=element_blank()) +
            theme(plot.title = element_text(hjust = 0.5)) +
            theme(legend.text = element_text(size = 18, family = "Helvetica"),
                  legend.title = element_text(size = 18, family = "Helvetica", face = "bold"), legend.direction = "vertical", legend.box = "vertical"))

