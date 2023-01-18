ALL_GFF_BGC_contig <- read.delim("~/Documents/Rocks/ALL_GFF_BGC_contig.tsv")
ALL_GFF_BGC_contig <- ALL_GFF_BGC_contig %>% filter(sequence_id != "sequence_id") %>% filter(detector == "antismash") %>%  unique()


Rock_AMR_taxa <- read.delim("~/Documents/Rocks/Rock_AMR_taxa.tsv")
Rock_AMR_taxa<- Rock_AMR_taxa %>% filter(ORF != "ORF")
Rock_AMR_taxa$RNum_Gi <- as.numeric(as.character(Rock_AMR_taxa$RNum_Gi))
Rock_AMR_taxa <- Rock_AMR_taxa %>% filter(Sample != "GL_R18_GL16_UP_2")
Rock_AMR_taxa$classification <- fct_explicit_na(Rock_AMR_taxa$classification, na_level = "unclassified")
Rock_AMR_taxa$AMR_category <- gsub("unclassified","Unclassified", Rock_AMR_taxa$AMR_category)
Rock_AMR_taxa$AMR_category <- gsub("fleuroquinolone","fluoroquinolone", Rock_AMR_taxa$AMR_category)
Rock_AMR_select <- Rock_AMR_taxa %>% select(AMR_category,Chr,Sample)

BGC_select <- ALL_GFF_BGC_contig[,c(3,4,5,6)]
BGC_select[BGC_select == ""] <- NA
BGC_select$product_activity <- fct_explicit_na(BGC_select$product_activity, na_level = "unknown")
BGC_select$product_class <- fct_explicit_na(BGC_select$class, na_level = "unknown")

AMR_BGC_select <- merge(Rock_AMR_select, BGC_select, by.x = "Chr", by.y = "contig")
AMR_BGC_select$freq <- 1

rock_bacteria <- read.delim("~/Documents/Rocks/rock_gtdbtk.bac120.summary.tsv") %>% select(1,2)
rock_bacteria$user_genome <- gsub(".contig.*","",rock_bacteria$user_genome)
colnames(rock_bacteria)[1] <- "bin"
rock_bacteria$group <- "bacteria"

AMR_BGC_select <- left_join(AMR_BGC_select,rock_bacteria)

AMR_BGC_select <- AMR_BGC_select %>% separate(classification, into = c("Domain","Phylum","Class","Order","Family","Genus","Species"), sep = ";") %>%
  mutate(Domain = str_replace_all(Domain, "d__", "")) %>%
  mutate(Phylum = str_replace_all(Phylum, "p__", ""))  %>%
  mutate(Class = str_replace_all(Class, "c__", ""))  %>%
  mutate(Order = str_replace_all(Order, "o__", "")) %>%
  mutate(Family = str_replace_all(Family, "f__", ""))  %>%
  mutate(Genus = str_replace_all(Genus, "g__", ""))  %>%
  mutate(Species = str_replace_all(Species, "s__", ""))
AMR_BGC_select$product_class <- fct_explicit_na(as.factor(AMR_BGC_select$product_class), na_level = "unknown")
AMR_BGC_select[AMR_BGC_select == ""] <- "unclassified"

GRAPH <- AMR_BGC_select %>% group_by(AMR_category,Genus, product_activity, product_class)%>% summarise(freq=sum(freq))

nb.cols <- length(unique(GRAPH$Genus))
aggaltColors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)

GRAPH %>% ggplot(aes(y= freq, axis1 = AMR_category, axis2 = Genus, axis3 = product_activity)) +
  geom_alluvium(aes(fill = Genus), reverse = FALSE, knot.pos = 0,) +
  guides(fill = FALSE) +
  geom_stratum(width = 2/8, reverse = FALSE, aes(fill = Genus)) +
  scale_fill_manual(values = aggaltColors, breaks = unique(GRAPH$Genus)) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            reverse = FALSE) +
  scale_x_discrete(limits = c("AMR", "Genus", "BGC activity"), expand = c(.2, .05)) +
  scale_linetype_manual(values = c("blank", "solid"))

ggplot(GRAPH,
       aes(y = freq,
           axis1 = AMR_category, axis2 = Genus, axis3 = product_activity)) +
  geom_alluvium(aes(fill = Genus),
                width = 1/8, reverse = FALSE) +
  scale_fill_manual(values = aggaltColors, breaks = unique(GRAPH$Genus)) +
  guides(fill = FALSE) +
  geom_stratum(alpha = 0, width = 1/8, reverse = FALSE) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            reverse = FALSE) +
  scale_x_continuous(breaks = 1:3, labels = c("AMR", "Genus", "BGC activity"), expand = c(.1, .05)) 
