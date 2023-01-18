#####################
#   General AMR     #
#####################

# Load and format files
rock_AMR_FC <- read.delim("~/Documents/Rocks/rock_AMR_FC.tsv")
rock_AMR_FC <- rock_AMR_FC %>% filter(ORF != "ORF")
rock_AMR_FC$RNum_Gi <- as.numeric(as.character(rock_AMR_FC$RNum_Gi))
rock_AMR_FC$AMR_category <- gsub("unclassified","Unclassified", rock_AMR_FC$AMR_category)
rock_AMR_FC$AMR_category <- gsub("fleuroquinolone","fluoroquinolone", rock_AMR_FC$AMR_category)

rock_AMR_FC <- rock_AMR_FC %>% filter(Sample != "GL_R18_GL16_UP_2") %>% na.omit()
rock_AMR_cat <- rock_AMR_FC %>% group_by(AMR_category, Sample) %>% summarise(RNum_Gi=sum(RNum_Gi)) %>% filter(AMR_category != "-")%>% na.omit()
rock_AMR_cat$Sample <- gsub("GL_","", rock_AMR_cat$Sample)
rock_AMR_cat$Sample <- gsub("_[0-9]","", rock_AMR_cat$Sample)
rock_AMR_cat <- separate(rock_AMR_cat, Sample, into = c("Rock","Glacier","Location"), sep = "_")

# Circulized barplot
rock_AMR_circ <- rock_AMR_cat %>% select(1,2,3,5)
colnames(rock_AMR_circ) <- c("observation","individual","group","value")
rock_AMR_circ <- rock_AMR_circ %>% group_by(observation, individual, group) %>% summarise(value=sum(value))
rock_AMR_circ <- rock_AMR_circ %>% spread(observation, value)
rock_AMR_circ[is.na(rock_AMR_circ)] <- 0
rock_AMR_circ<- rock_AMR_circ %>% gather(key = "observation", value="value", -c(1,2)) 
rock_AMR_circ <- rock_AMR_circ %>% arrange(group, individual)
rock_AMR_circ <- data.frame(rock_AMR_circ)
rock_AMR_circ$group <- as.factor(rock_AMR_circ$group) 
rock_AMR_circ$individual <- as.factor(rock_AMR_circ$individual) 

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 1
nObsType <- nlevels(as.factor(rock_AMR_circ$observation))
to_add <- data.frame( matrix(NA, empty_bar*nlevels(rock_AMR_circ$group)*nObsType, ncol(rock_AMR_circ)) )
colnames(to_add) <- colnames(rock_AMR_circ)
to_add$group <- rep(levels(rock_AMR_circ$group), each=empty_bar*nObsType )
rock_AMR_circ <- rbind(rock_AMR_circ, to_add)
rock_AMR_circ <- rock_AMR_circ %>% arrange(group, individual)
rock_AMR_circ$id <- rep( seq(1, nrow(rock_AMR_circ)/nObsType) , each=nObsType)

# Get the name and the y position of each label
label_data <- rock_AMR_circ %>% group_by(id, individual) %>% summarize(tot=sum(value))
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# Make it fancy
rock_AMR_circ = rock_AMR_circ %>% arrange(group, value)

## prepare a data frame for base lines
base_data <- rock_AMR_circ %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

nb.cols <- length(unique(rock_AMR_circ$observation))
myColors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)


# Make the plot
p <- ggplot(rock_AMR_circ) +      
  geom_bar(aes(x=as.factor(id), y=value, fill=observation), stat="identity", alpha=0.8) +
  geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "black", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0.001, xend = start, yend = 0.001), colour = "black", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  scale_fill_manual(
    values = c(
      "acridine_dye"="skyblue3",
      "aminocoumarin"="#E64B35B2",
      "aminoglycoside"="#F7DC6F",
      "aminoglycoside:aminocoumarin"="firebrick",
      "antibacterial_free_fatty_acids"="#00A087B2",
      "bacitracin"="#3C5488B2",
      "beta-lactam"="#F39B7FB2",
      "bicyclomycin"="#8491B4B2",
      "diaminopyrimidine"="#91D1C2B2",
      "elfamycin"="palegreen4",
      "fluoroquinolone"="#DC0000B3",
      "fosfomycin"="#4DBBD5B2",
      "fusidic-acid"="tan3",
      "glycopeptide"="#58D68D",
      "MLS"="#3498DB",
      "multidrug"="#2ECC71",
      "mupirocin"="#F39C12",
      "nitroimidazole"="#EE9A60",
      "nucleoside"="#7A94E9",
      "peptide"="#AED6F1",
      "phenicol"="#C0C0C0",
      "pleuromutilin"="#808000",
      "polyamine:peptide"="#922B21",
      "polymyxin"="#8F283A",
      "rifamycin"="#48ACC2",
      "sulfonamide"="#8F283A",
      "tetracycline"="#236A86",
      "triclosan"="#922B21",
      "Unclassified"="#008080"
    )
  ) +
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
  #  ggplot2::annotate("text", x = rep(max(rock_AMR_circ$id),4), y = c(0,0.00000001, 0.000001, 0.0001), label = c("0", "1e-8", "1e-6", "1e-4") , color="grey", size=6 , angle=0, fontface="bold", hjust=1) +
  geom_text(data=label_data, aes(x=id, y=tot, label=individual, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=3, angle= label_data$angle, inherit.aes = FALSE ) +
  geom_segment(data=base_data, aes(x = start, y = -0.00005, xend = end, yend = -0.00005), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = -0.0002, label=group), colour = "black", alpha=0.8, size=2.5, fontface="bold", inherit.aes = FALSE)

p

ggsave2("~/Documents/Rocks/Figures/Figure1a.tiff", plot = p, device = "tiff", dpi = 350 )

#####################
#     AMR TAXA      #
#####################

# Load and format data
Rock_AMR_taxa <- read.delim("~/Documents/Rocks/Rock_AMR_taxa.tsv")
Rock_AMR_taxa<- Rock_AMR_taxa %>% filter(ORF != "ORF")
Rock_AMR_taxa$RNum_Gi <- as.numeric(as.character(Rock_AMR_taxa$RNum_Gi))
Rock_AMR_taxa <- Rock_AMR_taxa %>% filter(Sample != "GL_R18_GL16_UP_2")
Rock_AMR_taxa$classification <- fct_explicit_na(Rock_AMR_taxa$classification, na_level = "unclassified")
Rock_AMR_taxa$AMR_category <- gsub("unclassified","Unclassified", Rock_AMR_taxa$AMR_category)
Rock_AMR_taxa$AMR_category <- gsub("fleuroquinolone","fluoroquinolone", Rock_AMR_taxa$AMR_category)
Rock_AMR_taxa <- Rock_AMR_taxa %>% group_by(AMR_category, Sample, classification, full_classification) %>% filter(AMR_category != "-") %>% summarise(RNum_Gi = sum(RNum_Gi)) 

total <- Rock_AMR_taxa %>% group_by(Sample) %>% summarise(sum=sum(RNum_Gi))
Rock_AMR_taxa <- left_join(Rock_AMR_taxa, total)
Rock_AMR_taxa$relab <-(Rock_AMR_taxa$RNum_Gi / Rock_AMR_taxa$sum) * 100
Rock_AMR_taxa$Sample <- gsub("GL_","", Rock_AMR_taxa$Sample)
Rock_AMR_taxa$Sample <- gsub("_[0-9]","", Rock_AMR_taxa$Sample)
Rock_AMR_taxa <- separate(Rock_AMR_taxa, Sample, into = c("Rock","Glacier","Location"), sep = "_")
Rock_AMR_taxa_overall <- Rock_AMR_taxa %>% group_by(Rock,Glacier, classification) %>% summarise(relab = sum(relab)) %>% filter(classification != "unclassified")

# Diverging Barplot
div_AMR_taxa <- Rock_AMR_taxa_overall %>% filter(classification != "Archaea")
div_AMR_taxa <- div_AMR_taxa %>% mutate(relabInv = ifelse(classification == "Bacteria", relab, relab*-1))

p <- div_AMR_taxa %>% 
  ggplot(aes(x=Rock, y=relabInv, fill=classification))+
  geom_bar(stat="identity",position="identity")+
  xlab("samples")+ylab("relative abundance")+
  scale_fill_manual(name="Microbe",values = c("#FFA373","#50486D"))+
  coord_flip()+
  facet_grid(rows = vars(Glacier), scales = "free_y", space = "free_y", switch = "y") +
  theme(strip.background = element_rect(fill = "white")) +
  geom_hline(yintercept=0)+
  scale_y_continuous(breaks = pretty(div_AMR_taxa$relabInv),labels = abs(pretty(div_AMR_taxa$relabInv)))+
  theme_scientific()

p

ggsave2("~/Documents/Rocks/Figures/Figure1b.tiff", plot = p, device = "tiff", dpi = 350 )

# flow
div_AMR_taxa_2 <- cSplit(Rock_AMR_taxa, "full_classification", ";")
div_AMR_taxa_2 <- div_AMR_taxa_2 %>% filter(classification != "Archaea")
div_AMR_taxa_2 <- div_AMR_taxa_2 %>% group_by(Rock, Glacier, classification, full_classification_3) %>% summarise(relab = sum(relab)) %>% na.omit()

Eukaryota <- div_AMR_taxa_2 %>% filter(classification == "Eukaryota") %>% spread(full_classification_3, relab)
Eukaryota[is.na(Eukaryota)] <- 0
Eukaryota <- Eukaryota %>% gather(full_classification_3, relab, 4:24)

nb.cols <- length(unique(Eukaryota$full_classification_3))
myColors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)

p <- ggplot(Eukaryota,
      aes(x = Rock, stratum = full_classification_3, alluvium = full_classification_3,
          y = relab,
          fill = full_classification_3, label = full_classification_3)) +
  geom_flow() +
  geom_stratum(alpha = .5) +
  scale_fill_manual(values = myColors) +
  facet_grid(~ Glacier, scales = "free_x", space = "free")+
  theme(strip.background = element_rect(fill = "white")) +
  theme_scientific()

p
ggsave2("~/Documents/Rocks/Figures/Figure2a1.tiff", plot = p, width = 40 , height =15.8, units = "cm", device = "tiff", dpi = 350 )

div_AMR_taxa_3 <- cSplit(Rock_AMR_taxa, "full_classification", ";")
div_AMR_taxa_3 <- div_AMR_taxa_3 %>% filter(classification != "Archaea")
div_AMR_taxa_3 <- div_AMR_taxa_3 %>% group_by(Rock, Glacier, classification, full_classification_5) %>% summarise(relab = sum(relab)) %>% na.omit()

Bacteria <- div_AMR_taxa_2 %>% filter(classification == "Bacteria") %>% spread(full_classification_3, relab)
Bacteria[is.na(Bacteria)] <- 0
Bacteria <- Bacteria %>% gather(full_classification_3, relab, 4:29)

nb.cols <- length(unique(Bacteria$full_classification_3))
myColors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)

p <- ggplot(Bacteria,
       aes(x = Rock, stratum = full_classification_3, alluvium = full_classification_3,
           y = relab,
           fill = full_classification_3, label = full_classification_3)) +
  geom_flow() +
  geom_stratum(alpha = .5) +
  scale_fill_manual(values = myColors) +
  facet_grid(~ Glacier, scales = "free_x", space = "free")+
  theme(strip.background = element_rect(fill = "white")) +
  theme_scientific()

ggsave2("~/Documents/Rocks/Figures/Figure2a2.tiff", plot = p, device = "tiff", dpi = 350 )

# BUBLE PLOT
div_AMR_taxa_2 <- cSplit(Rock_AMR_taxa, "full_classification", ";")
div_AMR_taxa_2$full_classification_3 <- as.character(div_AMR_taxa_2$full_classification_3)
div_AMR_taxa_2 <- div_AMR_taxa_2 %>% filter(classification == "Bacteria" | classification == "Eukaryota")
div_AMR_taxa_2 <- div_AMR_taxa_2 %>% group_by(classification, AMR_category, full_classification_3) %>% summarise_at(vars(RNum_Gi), list(mean=mean))%>% na.omit()

low_occur <- div_AMR_taxa_2 %>% group_by(full_classification_3) %>% summarise(sum_mean=sum(mean))
mean(low_occur$sum_mean) # 1.699499e-06
low_occ <- low_occur %>% filter(sum_mean <= 1.699499e-06)
low_occ <- low_occ$full_classification_3
div_AMR_taxa_2$full_classification_3[div_AMR_taxa_2$full_classification_3 %in% low_occ] = "Others"
div_AMR_taxa_2 <- div_AMR_taxa_2 %>% group_by(classification, AMR_category, full_classification_3) %>% summarise(mean = sum(mean))

order_bac <- div_AMR_taxa_2 %>% filter(classification == "Bacteria")
order <-div_AMR_taxa_2 %>% group_by(full_classification_3) %>% tally()
order <- arrange(order, desc(n))

nb.cols <- length(unique(order$full_classification_3))
myColors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)

p <- div_AMR_taxa_2 %>% mutate(classification = factor(classification, levels = c("Eukaryota","Bacteria"))) %>%
  mutate(full_classification_3 = factor(full_classification_3, levels = rev(order$full_classification_3))) %>%
  ggplot(aes(x=AMR_category, y=full_classification_3, color = full_classification_3, size=mean)) +
  geom_point(alpha=0.7, stat = "identity") +
  facet_grid(rows = vars(factor(classification, levels = c("Eukaryota","Bacteria"))), scales = "free_y", switch = "y") +
  guides(color = FALSE) +
  scale_size(range = c(.2, 14), name="Abundance") +
  theme(strip.background = element_rect(fill = "white")) +
  theme_scientific() +
  scale_color_manual(values = myColors) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 

ggsave2("~/Documents/Rocks/Figures/Figure2b.tiff", plot = p, device = "tiff", dpi = 350 )



