library(circlize)
library(ComplexHeatmap)

ALL_GFF_BGC_Bin <- read.delim("~/Documents/Rocks/ALL_GFF_BGC_Bin.tsv")
ALL_GFF_BGC_Bin <- ALL_GFF_BGC_Bin %>% filter(sequence_id != "sequence_id") %>% filter(detector == "antismash") %>% unique()
ALL_GFF_BGC_Bin[ALL_GFF_BGC_Bin==""]<- NA
ALL_GFF_BGC_Bin$product_activity <- fct_explicit_na(ALL_GFF_BGC_Bin$product_activity, na_level = "unknown")
ALL_GFF_BGC_Bin$product_class <- fct_explicit_na(ALL_GFF_BGC_Bin$product_class, na_level = "unknown")


rock_bacteria <- read.delim("~/Documents/Rocks/rock_gtdbtk.bac120.summary.tsv") %>% select(1,2)
rock_bacteria$user_genome <- gsub(".contig.*","",rock_bacteria$user_genome)
colnames(rock_bacteria)[1] <- "bin"
rock_bacteria$group <- "bacteria"

rock_euk <- read.csv2("~/Documents/Rocks/SummaryMAGs.csv") %>% select(1,4)
colnames(rock_euk) <- c("bin","classification")
rock_euk$group <- "eukaryota"

rock_taxa <- rbind(rock_bacteria, rock_euk)

ALL_GFF_BGC_Bin <- inner_join(ALL_GFF_BGC_Bin, rock_taxa)

# overall 
BGC_overall <- ALL_GFF_BGC_Bin %>% group_by(product_activity,group) %>% tally()

nb.cols <- length(unique(BGC_overall$product_activity))
myColors <- colorRampPalette(brewer.pal(6, "Set1"))(nb.cols)

BGC_overall %>% ggplot(aes(x=group, y=n, fill = product_activity)) +
  geom_bar(stat = "identity", alpha=0.8) +
  xlab("") + ylab("count") +
  ggtitle("BGC") +
  scale_fill_manual(values = myColors) +
  theme_scientific()


BGC_overall_antb <- ALL_GFF_BGC_Bin %>% filter(product_activity == "antibacterial") %>% group_by(product_class,group) %>% tally()
nb.cols <- length(unique(BGC_overall_antb$product_class))
myColors <- colorRampPalette(brewer.pal(6, "Set1"))(nb.cols)

BGC_overall_antb %>% ggplot(aes(x=group, y=n, fill = product_class)) +
  geom_bar(stat = "identity", alpha=0.8) +
  xlab("") + ylab("count") +
  ggtitle("BGC") +
  scale_fill_manual(values = myColors) +
  theme_scientific()

# Taxa_specific
BGC_eukaryota <- ALL_GFF_BGC_Bin %>% filter(group == "eukaryota") %>% filter(product_activity == "antibacterial") %>% group_by(product_class, classification) %>% tally()

BGC_eukaryota%>% ggplot(aes(x=classification, y=n, fill = product_class)) +
  geom_bar(stat = "identity", alpha=0.8) +
  xlab("") + ylab("count") +
  ggtitle("BGC") +
  scale_fill_manual(values = myColors) +
  theme_scientific()


BGC_eukaryota <- ALL_GFF_BGC_Bin %>% filter(group == "eukaryota") %>% group_by(product_activity, classification) %>% tally()

BGC_eukaryota%>% ggplot(aes(x=classification, y=n, fill = product_activity)) +
  geom_bar(stat = "identity", alpha=0.8) +
  xlab("") + ylab("count") +
  ggtitle("BGC") +
  scale_fill_manual(values = myColors) +
  theme_scientific()



# Heatmaps
BGC_bacteria <- ALL_GFF_BGC_Bin %>% filter(group == "bacteria") %>% group_by(product_activity, classification) %>% tally()

BGC_bacteria <- BGC_bacteria %>% separate(classification, into = c("Domain","Phylum","Class","Order","Family","Genus","Species"), sep = ";") %>%
  mutate(Domain = str_replace_all(Domain, "d__", "")) %>%
  mutate(Phylum = str_replace_all(Phylum, "p__", ""))  %>%
  mutate(Class = str_replace_all(Class, "c__", ""))  %>%
  mutate(Order = str_replace_all(Order, "o__", "")) %>%
  mutate(Family = str_replace_all(Family, "f__", ""))  %>%
  mutate(Genus = str_replace_all(Genus, "g__", ""))  %>%
  mutate(Species = str_replace_all(Species, "s__", ""))

BGC_bacteria <- BGC_bacteria %>% spread(product_activity, n)
BGC_bacteria[is.na(BGC_bacteria)] <- 0
single_occ = names(table(BGC_bacteria$Order))[table(BGC_bacteria$Order) == 1]
BGC_bacteria$Order[BGC_bacteria$Order %in% single_occ] = 'Others'
single_occ = names(table(BGC_bacteria$Phylum))[table(BGC_bacteria$Phylum) == 1]
BGC_bacteria$Phylum[BGC_bacteria$Phylum %in% single_occ] = 'Others'

BGC_bacteria[BGC_bacteria == ""] <- 'Others'

col_fun = colorRamp2(c(0,5), c("white", "darkblue"))
col_fun(seq(0, 120))
cols_orders = c('Others' = 'grey',
                "Actinomycetales" = 'steelblue',
                "AKYH767" = '#FC7B6F',
                "Chitinophagales" = '#DB352F',
                "Cytophagales" = '#FF673F',
                "Flavobacteriales" = '#E35F07',
                "Sphingobacteriales" = '#FF1101',
                "Bdellovibrionales" = '#FAD706',
                "Cyanobacteriales" = 'aquamarine',
                "Deinococcales" = '#B37F4F',
                "Polyangiales" = '#68B019',
                "Absconditabacterales" = 'red',
                'Rhizobiales' = 'forestgreen',
                "Rhodobacterales" = 'tomato',
                "Rickettsiales" = '#68D4EB',
                "Sphingomonadales" = '#5A7FC7',
                "Burkholderiales" = '#3AE9CF',
                "Xanthomonadales" = '#2B59E3',
                "Verrucomicrobiales" = '#FF4374')


H1 <- Heatmap(as.matrix(BGC_bacteria[8:13]),
        column_title_rot = 90,
        col = col_fun,
        row_split = as.vector(BGC_bacteria$Phylum),
        row_title_rot = 0,
        row_order = order(BGC_bacteria$Order),
        left_annotation = rowAnnotation(Order = BGC_bacteria$Order, col= list(Order = cols_orders)),
        #        column_split = data_cat$category,
        show_heatmap_legend = F,
        show_column_dend = FALSE,
        show_row_names = FALSE)

# specific
BGC_bacteria <- ALL_GFF_BGC_Bin %>% filter(group == "bacteria") %>% filter(product_activity == "antibacterial") %>% group_by(product_class, classification) %>% tally()

BGC_bacteria <- BGC_bacteria %>% separate(classification, into = c("Domain","Phylum","Class","Order","Family","Genus","Species"), sep = ";") %>%
  mutate(Domain = str_replace_all(Domain, "d__", "")) %>%
  mutate(Phylum = str_replace_all(Phylum, "p__", ""))  %>%
  mutate(Class = str_replace_all(Class, "c__", ""))  %>%
  mutate(Order = str_replace_all(Order, "o__", "")) %>%
  mutate(Family = str_replace_all(Family, "f__", ""))  %>%
  mutate(Genus = str_replace_all(Genus, "g__", ""))  %>%
  mutate(Species = str_replace_all(Species, "s__", ""))

BGC_bacteria <- BGC_bacteria %>% spread(product_class, n)
BGC_bacteria[is.na(BGC_bacteria)] <- 0
single_occ = names(table(BGC_bacteria$Order))[table(BGC_bacteria$Order) == 1]
BGC_bacteria$Order[BGC_bacteria$Order %in% single_occ] = 'Others'
single_occ = names(table(BGC_bacteria$Phylum))[table(BGC_bacteria$Phylum) == 1]
BGC_bacteria$Phylum[BGC_bacteria$Phylum %in% single_occ] = 'Others'

BGC_bacteria[BGC_bacteria == ""] <- 'Others'

col_fun = colorRamp2(c(0,5), c("white", "darkblue"))
col_fun(seq(0, 120))
cols_orders = c('Others' = 'grey',
                "Actinomycetales" = 'steelblue',
                "AKYH767" = '#FC7B6F',
                "Chitinophagales" = '#DB352F',
                "Cytophagales" = '#FF673F',
                "Flavobacteriales" = '#E35F07',
                "Sphingobacteriales" = '#FF1101',
                "Bdellovibrionales" = '#FAD706',
                "Cyanobacteriales" = 'aquamarine',
                "Deinococcales" = '#B37F4F',
                "Polyangiales" = '#68B019',
                "Absconditabacterales" = 'red',
                'Rhizobiales' = 'forestgreen',
                "Rhodobacterales" = 'tomato',
                "Rickettsiales" = '#68D4EB',
                "Sphingomonadales" = '#5A7FC7',
                "Burkholderiales" = '#3AE9CF',
                "Xanthomonadales" = '#2B59E3',
                "Verrucomicrobiales" = '#FF4374')


Heatmap(as.matrix(BGC_bacteria[8:17]),
        column_title_rot = 90,
        col = col_fun,
        row_split = as.vector(BGC_bacteria$Phylum),
        row_title_rot = 0,
        row_order = order(BGC_bacteria$Order),
        left_annotation = rowAnnotation(Order = BGC_bacteria$Order, col= list(Order = cols_orders)),
        #        column_split = data_cat$category,
        show_heatmap_legend = F,
        show_column_dend = FALSE,
        show_row_names = FALSE)




# Heatmap Eukaryota
BGC_eukaryota <- ALL_GFF_BGC_Bin %>% filter(group == "eukaryota") %>% group_by(product_activity, classification) %>% tally()

#BGC_eukaryota <- BGC_eukaryota %>% separate(classification, into = c("Domain","Phylum","Class","Order","Family","Genus","Species"), sep = ";") %>%
#  mutate(Domain = str_replace_all(Domain, "d__", "")) %>%
#  mutate(Phylum = str_replace_all(Phylum, "p__", ""))  %>%
#  mutate(Class = str_replace_all(Class, "c__", ""))  %>%
#  mutate(Order = str_replace_all(Order, "o__", "")) %>%
#  mutate(Family = str_replace_all(Family, "f__", ""))  %>%
#  mutate(Genus = str_replace_all(Genus, "g__", ""))  %>%
#  mutate(Species = str_replace_all(Species, "s__", ""))

BGC_eukaryota <- BGC_eukaryota %>% spread(product_activity, n)
BGC_eukaryota[is.na(BGC_eukaryota)] <- 0
#single_occ = names(table(BGC_eukaryota$Order))[table(BGC_eukaryota$Order) == 1]
#BGC_eukaryota$Order[BGC_eukaryota$Order %in% single_occ] = 'Others'
#single_occ = names(table(BGC_eukaryota$Phylum))[table(BGC_eukaryota$Phylum) == 1]
#BGC_eukaryota$Phylum[BGC_eukaryota$Phylum %in% single_occ] = 'Others'

#BGC_eukaryota[BGC_eukaryota == ""] <- 'Others'

col_fun = colorRamp2(c(0,10), c("white", "darkblue"))
col_fun(seq(0, 150))
cols_orders = c("Bacillariophyta" = 'steelblue',
                "Chytridiomycetes" = '#FC7B6F',
                "Clunio marinus" = '#FAD706',
                "Eukaryota" = '#B37F4F',
                "Fragilariophycidae" = '#68B019',
                "Fungi" = '#68D4EB',
                "Ochromonas" = '#3AE9CF')


H2 <- Heatmap(as.matrix(BGC_eukaryota[2:4]),
        column_title_rot = 90,
        col = col_fun,
        row_split = as.vector(BGC_eukaryota$classification),
        row_title_rot = 0,
        row_order = order(BGC_eukaryota$classification),
        left_annotation = rowAnnotation(Order = BGC_eukaryota$classification, col= list(Order = cols_orders)),
        #        column_split = data_cat$category,
        show_heatmap_legend = F,
        show_column_dend = FALSE,
        show_row_names = FALSE)


# specific
BGC_eukaryota <- ALL_GFF_BGC_Bin %>% filter(group == "eukaryota") %>% filter(product_activity == "antibacterial") %>% group_by(product_class, classification) %>% tally()

BGC_eukaryota <- BGC_eukaryota %>% spread(product_class, n)
BGC_eukaryota[is.na(BGC_eukaryota)] <- 0
#single_occ = names(table(BGC_eukaryota$Order))[table(BGC_eukaryota$Order) == 1]
#BGC_eukaryota$Order[BGC_eukaryota$Order %in% single_occ] = 'Others'
#single_occ = names(table(BGC_eukaryota$Phylum))[table(BGC_eukaryota$Phylum) == 1]
#BGC_eukaryota$Phylum[BGC_eukaryota$Phylum %in% single_occ] = 'Others'

#BGC_eukaryota[BGC_eukaryota == ""] <- 'Others'

col_fun = colorRamp2(c(0,10), c("white", "darkblue"))
col_fun(seq(0, 150))
cols_orders = c("Bacillariophyta" = 'steelblue',
                "Chytridiomycetes" = '#FC7B6F',
                "Clunio marinus" = '#FAD706',
                "Eukaryota" = '#B37F4F',
                "Fragilariophycidae" = '#68B019',
                "Fungi" = '#68D4EB',
                "Ochromonas" = '#3AE9CF')


Heatmap(as.matrix(BGC_eukaryota[2:8]),
        column_title_rot = 90,
        col = col_fun,
        row_split = as.vector(BGC_eukaryota$classification),
        row_title_rot = 0,
        row_order = order(BGC_eukaryota$classification),
        left_annotation = rowAnnotation(Order = BGC_eukaryota$classification, col= list(Order = cols_orders)),
        #        column_split = data_cat$category,
        show_heatmap_legend = F,
        show_column_dend = FALSE,
        show_row_names = FALSE)
