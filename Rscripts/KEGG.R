Rock_kegg_antb <- read.delim("~/Documents/Rocks/Rock_kegg_antb.tsv")
Rock_kegg_antb <- Rock_kegg_antb %>% filter(Chr != "Chr")
Rock_kegg_antb$RNum_Gi <- as.numeric(as.character(Rock_kegg_antb$RNum_Gi))
Rock_kegg_antb$Sample <- gsub("_contig.*","",Rock_kegg_antb$Chr)
Rock_kegg_antb <- Rock_kegg_antb %>% filter(Sample != "Gl_R18_GL16_UP_2")
Rock_kegg_antb <- Rock_kegg_antb %>% group_by(PATHWAY, Sample) %>% summarise(RNum_Gi = sum(RNum_Gi))
Rock_kegg_antb$Sample <- gsub("GL_","",Rock_kegg_antb$Sample)
Rock_kegg_antb$Sample <- gsub("_[0-9]","",Rock_kegg_antb$Sample)
Rock_kegg_antb <- separate(Rock_kegg_antb, Sample,  into = c("Rock","Glacier","Location"), sep = "_")

# CIRCLE PLOT
# AMR CAT 
Rock_kegg_circ <- Rock_kegg_antb %>% select(1,2,3,5)
colnames(Rock_kegg_circ) <- c("observation","individual","group","value")
Rock_kegg_circ <- Rock_kegg_circ %>% group_by(observation, individual, group) %>% summarise(value=sum(value))
Rock_kegg_circ <- Rock_kegg_circ %>% spread(observation, value)
Rock_kegg_circ[is.na(Rock_kegg_circ)] <- 0
Rock_kegg_circ<- Rock_kegg_circ %>% gather(key = "observation", value="value", -c(1,2)) 
Rock_kegg_circ <- Rock_kegg_circ %>% arrange(group, individual)
Rock_kegg_circ <- data.frame(Rock_kegg_circ)
Rock_kegg_circ$group <- as.factor(Rock_kegg_circ$group) 
Rock_kegg_circ$individual <- as.factor(Rock_kegg_circ$individual) 

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 1
nObsType <- nlevels(as.factor(Rock_kegg_circ$observation))
to_add <- data.frame( matrix(NA, empty_bar*nlevels(Rock_kegg_circ$group)*nObsType, ncol(Rock_kegg_circ)) )
colnames(to_add) <- colnames(Rock_kegg_circ)
to_add$group <- rep(levels(Rock_kegg_circ$group), each=empty_bar*nObsType )
Rock_kegg_circ <- rbind(Rock_kegg_circ, to_add)
Rock_kegg_circ <- Rock_kegg_circ %>% arrange(group, individual)
Rock_kegg_circ$id <- rep( seq(1, nrow(Rock_kegg_circ)/nObsType) , each=nObsType)

# Get the name and the y position of each label
label_data <- Rock_kegg_circ %>% group_by(id, individual) %>% summarize(tot=sum(value))
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# Make it fancy
Rock_kegg_circ = Rock_kegg_circ %>% arrange(group, value)

## prepare a data frame for base lines
base_data <- Rock_kegg_circ %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

library(RColorBrewer)
nb.cols <- length(unique(Rock_kegg_circ$observation))
myColors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)


# Make the plot
p <- ggplot(Rock_kegg_circ) +      
  geom_bar(aes(x=as.factor(id), y=value, fill=observation), stat="identity", alpha=0.8) +
  geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "black", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0.001, xend = start, yend = 0.001), colour = "black", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  scale_fill_manual(values = myColors) +
  ylim(-0.0008,max(label_data$tot, na.rm=T)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,5), "cm") 
  ) +
  coord_polar() +
  #  ggplot2::annotate("text", x = rep(max(Rock_kegg_circ$id),4), y = c(0,0.00000001, 0.000001, 0.0001), label = c("0", "1e-8", "1e-6", "1e-4") , color="grey", size=6 , angle=0, fontface="bold", hjust=1) +
  geom_text(data=label_data, aes(x=id, y=tot, label=individual, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=3, angle= label_data$angle, inherit.aes = FALSE ) +
  geom_segment(data=base_data, aes(x = start, y = -0.00005, xend = end, yend = -0.00005), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = -0.0002, label=group), colour = "black", alpha=0.8, size=2.5, fontface="bold", inherit.aes = FALSE)

p

ggsave2("~/Documents/Rocks/Figures/Figure3a.tiff", plot = p, device = "tiff", dpi = 350 )
#####################
#    KEGG TAXA      #
#####################
Rock_kegg_taxa <- read.delim("~/Documents/Rocks/Rock_kegg_taxa.tsv")
Rock_kegg_taxa<- Rock_kegg_taxa %>% filter(Chr != "Chr")
Rock_kegg_taxa$RNum_Gi <- as.numeric(as.character(Rock_kegg_taxa$RNum_Gi))
Rock_kegg_taxa$Sample <- gsub("_contig.*","",Rock_kegg_taxa$Chr)
Rock_kegg_taxa$Sample <- gsub("GL_","",Rock_kegg_taxa$Sample)
Rock_kegg_taxa$Sample <- gsub("_[0-9]","",Rock_kegg_taxa$Sample)
Rock_kegg_taxa <- Rock_kegg_taxa %>% filter(Sample != "GL_R18_GL16_UP_2")
Rock_kegg_taxa$classification <- fct_explicit_na(Rock_kegg_taxa$classification, na_level = "unclassified")
Rock_kegg_taxa <- Rock_kegg_taxa %>% group_by(PATHWAY,Sample, classification, full_classification)%>% summarise(RNum_Gi = sum(RNum_Gi)) 

total <- Rock_kegg_taxa %>% group_by(Sample) %>% summarise(sum=sum(RNum_Gi))
Rock_kegg_taxa <- left_join(Rock_kegg_taxa, total)
Rock_kegg_taxa$relab <-(Rock_kegg_taxa$RNum_Gi / Rock_kegg_taxa$sum) * 100
Rock_kegg_taxa$Sample <- gsub("GL_","", Rock_kegg_taxa$Sample)
Rock_kegg_taxa$Sample <- gsub("_[0-9]","", Rock_kegg_taxa$Sample)
Rock_kegg_taxa <- separate(Rock_kegg_taxa, Sample, into = c("Rock","Glacier","Location"), sep = "_")

Rock_kegg_taxa_overall <- Rock_kegg_taxa %>% group_by(Rock,Glacier, classification) %>% summarise(relab = sum(relab)) %>% filter(classification != "unclassified")

# Diverging Barplot
div_kegg_taxa <- Rock_kegg_taxa_overall %>% filter(classification != "Archaea")
div_kegg_taxa <- div_kegg_taxa %>% mutate(relabInv = ifelse(classification == "Bacteria", relab, relab*-1))

p <- div_kegg_taxa %>% 
  ggplot(aes(x=Rock, y=relabInv, fill=classification))+
  geom_bar(stat="identity",position="identity")+
  xlab("samples")+ylab("relative abundance")+
  scale_fill_manual(name="Microbe",values = c("#FFA373","#50486D"))+
  coord_flip()+
  facet_grid(rows = vars(Glacier), scales = "free_y", space = "free_y", switch = "y") +
  theme(strip.background = element_rect(fill = "white")) +
  geom_hline(yintercept=0)+
  scale_y_continuous(breaks = pretty(div_kegg_taxa$relabInv),labels = abs(pretty(div_kegg_taxa$relabInv)))+
  theme_scientific()

ggsave2("~/Documents/Rocks/Figures/Figure3b.tiff", plot = p, device = "tiff", dpi = 350 )
# flow
div_kegg_taxa_2 <- cSplit(Rock_kegg_taxa, "full_classification", ";")
div_kegg_taxa_2 <- div_kegg_taxa_2 %>% filter(classification != "Archaea")
div_kegg_taxa_2 <- div_kegg_taxa_2 %>% group_by(Rock, Glacier, classification, full_classification_3) %>% summarise(relab = sum(relab)) %>% na.omit()

Eukaryota <- div_kegg_taxa_2 %>% filter(classification == "Eukaryota") %>% spread(full_classification_3, relab)
Eukaryota[is.na(Eukaryota)] <- 0
Eukaryota <- Eukaryota %>% gather(full_classification_3, relab, 4:23)

nb.cols <- length(unique(Eukaryota$full_classification_3))
myColors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)

ggplot(Eukaryota,
       aes(x = Rock, stratum = full_classification_3, alluvium = full_classification_3,
           y = relab,
           fill = full_classification_3, label = full_classification_3)) +
  geom_flow() +
  geom_stratum(alpha = .5) +
  scale_fill_manual(values = myColors) +
  facet_grid(~ Glacier, scales = "free_x", space = "free")+
  theme(strip.background = element_rect(fill = "white")) +
  theme_scientific()




div_kegg_taxa_3 <- cSplit(Rock_kegg_taxa, "full_classification", ";")
div_kegg_taxa_3 <- div_kegg_taxa_3 %>% filter(classification != "Archaea")
div_kegg_taxa_3 <- div_kegg_taxa_3 %>% group_by(Rock, Glacier, classification, full_classification_5) %>% summarise(relab = sum(relab)) %>% na.omit()

Bacteria <- div_kegg_taxa_2 %>% filter(classification == "Bacteria") %>% spread(full_classification_3, relab)
Bacteria[is.na(Bacteria)] <- 0
Bacteria <- Bacteria %>% gather(full_classification_3, relab, 4:29)

nb.cols <- length(unique(Bacteria$full_classification_3))
myColors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)

ggplot(Bacteria,
       aes(x = Rock, stratum = full_classification_3, alluvium = full_classification_3,
           y = relab,
           fill = full_classification_3, label = full_classification_3)) +
  geom_flow() +
  geom_stratum(alpha = .5) +
  scale_fill_manual(values = myColors) +
  facet_grid(~ Glacier, scales = "free_x", space = "free")+
  theme(strip.background = element_rect(fill = "white")) +
  theme_scientific()

# BUBLE PLOT
div_kegg_taxa_2 <- cSplit(Rock_kegg_taxa, "full_classification", ";")
div_kegg_taxa_2$full_classification_3 <- as.character(div_kegg_taxa_2$full_classification_3)
div_kegg_taxa_2 <- div_kegg_taxa_2 %>% filter(classification == "Bacteria" | classification == "Eukaryota")
div_kegg_taxa_2 <- div_kegg_taxa_2 %>% group_by(classification, PATHWAY, full_classification_3) %>% summarise_at(vars(RNum_Gi), list(mean=mean))%>% na.omit()

low_occur <- div_kegg_taxa_2 %>% group_by(full_classification_3) %>% summarise(sum_mean=sum(mean))
mean(low_occur$sum_mean) # 13.895803e-06
low_occ <- low_occur %>% filter(sum_mean <= 3.895803e-06)
low_occ <- low_occ$full_classification_3
div_kegg_taxa_2$full_classification_3[div_kegg_taxa_2$full_classification_3 %in% low_occ] = "Others"
div_kegg_taxa_2 <- div_kegg_taxa_2 %>% group_by(classification, PATHWAY, full_classification_3) %>% summarise(mean = sum(mean))

order <-div_kegg_taxa_2 %>% group_by(full_classification_3) %>% tally()
order <- arrange(order, desc(n))

nb.cols <- length(unique(order$full_classification_3))
myColors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)

div_kegg_taxa_2$PATHWAY <- gsub("Biosynthesis of 12-, 14- and 16-membered macrolides","macrolide",div_kegg_taxa_2$PATHWAY)
div_kegg_taxa_2$PATHWAY <- gsub("Biosynthesis of ansamycins","ansamycins",div_kegg_taxa_2$PATHWAY)
div_kegg_taxa_2$PATHWAY <- gsub("Biosynthesis of vancomycin group antibiotics","vancomycin",div_kegg_taxa_2$PATHWAY)
div_kegg_taxa_2$PATHWAY <- gsub("Penicillin and cephalosporin biosynthesis","penicillin/cephalosporin",div_kegg_taxa_2$PATHWAY)
div_kegg_taxa_2$PATHWAY <- gsub("Monobactam biosynthesis","monobactam",div_kegg_taxa_2$PATHWAY)
div_kegg_taxa_2$PATHWAY <- gsub("Streptomycin biosynthesis","streptomcyin",div_kegg_taxa_2$PATHWAY)
div_kegg_taxa_2$PATHWAY <- gsub("Tetracycline biosynthesis","tetracycline",div_kegg_taxa_2$PATHWAY)


div_kegg_taxa_2 %>% mutate(classification = factor(classification, levels = c("Eukaryota","Bacteria"))) %>%
  mutate(full_classification_3 = factor(full_classification_3, levels = rev(order$full_classification_3))) %>%
  ggplot(aes(x=PATHWAY, y=full_classification_3, color = full_classification_3, size=mean)) +
  geom_point(alpha=0.7, stat = "identity") +
  facet_grid(rows = vars(factor(classification, levels = c("Eukaryota","Bacteria"))), scales = "free_y", switch = "y") +
  guides(color = FALSE) +
  scale_size(range = c(.2, 8), name="Abundance") +
  theme(strip.background = element_rect(fill = "white")) +
  theme_scientific() +
  scale_color_manual(values = myColors) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 



