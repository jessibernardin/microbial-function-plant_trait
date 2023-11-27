#### Microbial function mediates leaf traits in a pitcher plant model system ####
#### Authors: Jessica R. Bernardin, Erica B. Young, Leonora S. Bittleston ####
#### last update : October 16, 2023 ####
#### Metagenome Assembled Genome (MAG) Analysis

# Load and install required packages
packages_to_load <- c(
  "tidyverse", "dplyr", "pheatmap", "qiime2R", "phyloseq", "this.path", "magrittr"
)

for (i in packages_to_load) { #Installs packages if not yet installed
  if (!require(i, character.only = TRUE)) install.packages(i)
}

setwd(this.path::here())

#### MAPPING TRANSCRIPTS AND FUNCTIONS TO MAGS ####
#read in csv blast results of all transcripts blasted against drep_mag blast database
mag_blast <- read.csv("drep_blast_output.csv", header=TRUE)
mag_qual <- read.csv("dRep_checkm2_quality_report.csv", header=TRUE)
#272602

mag_tax <- read.csv("dRep_mag_tax.csv", header=TRUE)
new_bin <- mag_tax[,1]
mag_tax <- mag_tax[,-1]
names(mag_tax)[1] <- "Taxon"
names(mag_tax)[2] <- "Feature.ID"

mag_tax_p <- parse_taxonomy(mag_tax)
mag_tax_p$bin <- rownames(mag_tax_p)
mag_tax_p <- cbind(new_bin, mag_tax_p)
mag_tax_filt<- mag_tax_p[,c(1,8,9)]
split_column <- strsplit(mag_blast$contig_bin, "_", fixed = TRUE)

# Create new columns for contig and bin
mag_blast$contig <- sapply(split_column, function(x) paste(x[1:2], collapse = "_"))
mag_blast$bin <- sapply(split_column, function(x) paste(x[-(1:2)], collapse = "_"))

#functional annotation
KEGG_KO <- read.csv("KEGG_KO.csv", header=TRUE)
KEGG_KO <- KEGG_KO[,-1]

#merge with KEGG_KO
mag_fun <- left_join(mag_blast, KEGG_KO, by = c("DA_contig_dream" = "gene_name"))
#359867

#remove those transcripts that aren't assigned KOs
mag_fun_filt <- na.omit(mag_fun[mag_fun$KO != "NA", ])
#231717

mag_fun_ko_def <- mag_fun_filt[,c(7,8,12)]

filtered_df <- mag_fun_ko_def %>% filter(grepl("aminopeptidase|chit", KO_definition, ignore.case = TRUE))

filtered_df <- merge(filtered_df, mag_tax_p, by.x="bin", by.y="bin")
filtered_df$mag_species <- paste(filtered_df$new_bin, filtered_df$Species, sep = "-")
filtered_df <- filtered_df %>%
  mutate(broad_enzyme = case_when(
    str_detect(KO_definition, "amino") ~ "Protease",
    str_detect(KO_definition, "chit") ~ "Chitinase"
  ))

#Fig7
ggplot(filtered_df, aes(x = mag_species, y = KO_definition, color=broad_enzyme)) +
  geom_jitter(alpha=.3) + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1),
                                                  legend.position = "none")+ scale_color_manual(values = c("black", "black"))

filtered_df_new <- filtered_df %>%
  group_by(mag_species, KO_definition) %>%
  mutate(value = n()) %>%
  ungroup()


ggplot(filtered_df_new, aes(x = mag_species, y = KO_definition, group=KO_definition)) + 
  geom_point(aes(size = value, fill = value), alpha = 0.75, shape = 21)

colours = c( "#A54657",  "#582630", "#F7EE7F", "#4DAA57","#F1A66A","#F26157", "#F9ECCC", "#679289", "#33658A",
             "#F6AE2D","#86BBD8")
colours2 = c( "#F7EE7F","#679289")
ggplot(filtered_df_new, aes(x = mag_species, y = broad_enzyme, group=broad_enzyme)) + 
  geom_point(aes(size = value, fill = value), alpha = 0.75, shape = 21)+ theme_classic()+
   scale_size_continuous(limits = c(1, 300), range = c(1,17), breaks = c(1,10,50,100,200,300)) + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 11), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.title = element_text(size = 12, face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "right")

filtered_df_new$value <- as.integer(filtered_df_new$value)

ggplot(filtered_df_new, aes(x = mag_species, y = broad_enzyme)) + 
  geom_point(aes(size = value, fill = broad_enzyme), alpha = 0.75, shape = 21) + 
  scale_size_continuous(limits = c(1, 300), range = c(1,17), breaks = c(10,50,100,300)) + 
  labs( x= "", y = "", size = "Richness of KOs", fill = "")  + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 11), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.title = element_text(size = 12, face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "right") +  
  scale_fill_manual(values = colours, guide = FALSE) 


filtered_df_new <- filtered_df_new[order(filtered_df_new$value, decreasing = TRUE),]  

ggplot(filtered_df_new, aes(x = mag_species, y = KO_definition)) + 
  geom_point(aes(size = value, fill = KO_definition), alpha = 0.75, shape = 21) + 
  scale_size_continuous(limits = c(1, 300), range = c(1,17), breaks = c(10,50,100,300)) + 
  labs( x= "", y = "", size = "Richness of KOs", fill = "")  + 
          theme(legend.key=element_blank(), 
                axis.text.x = element_text(colour = "black", size = 5, face = "bold", angle = 45, vjust = 0.3, hjust = 1), 
                axis.text.y = element_text(colour = "black", face = "bold", size = 5), 
                legend.text = element_text(size = 5, face ="bold", colour ="black"), 
                legend.title = element_text(size = 5, face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "right") +  
  scale_fill_manual(values = colours, guide = FALSE) 

filtered_df_new$broad_enzyme <- factor(filtered_df_new$broad_enzyme, levels = c("Chitinase", "Protease"))

filtered_df_new <- filtered_df_new %>%
  arrange(broad_enzyme)

filtered_df_new$KO_definition <- as.factor(filtered_df_new$KO_definition)

ggplot(filtered_df_new, aes(x = reorder(mag_species, -value), y = factor(KO_definition, levels=unique(KO_definition)), size = value, fill = KO_definition)) + 
  geom_point(alpha = 0.75, shape = 21) + 
  scale_size_continuous(limits = c(1, 300), range = c(1, 10), breaks = c(10, 50, 100, 300)) + 
  labs(x = "", y = "", size = "# Unique Contigs", fill = "") +theme_classic()+
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 5, vjust = 1.01, hjust = 1, angle=90), 
        axis.text.y = element_text(colour = "black", size = 5), 
        legend.text = element_text(size = 5, colour ="black"), 
        legend.position = "right") +  
  scale_fill_manual(values = colours, guide = FALSE)+ 
  scale_y_discrete(labels = function(x) str_wrap(x, width = 20))



ggplot(filtered_df_new, aes(y = reorder(mag_species, value), x = factor(KO_definition, levels=unique(KO_definition)), size = value, fill = KO_definition)) + 
  geom_point(alpha = 0.75, shape = 21) + 
  scale_size_continuous(limits = c(1, 300), range = c(1, 10), breaks = c(10, 50, 100, 300)) + 
  labs(x = "", y = "", size = "# Unique Contigs", fill = "") +theme_classic()+
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 5, vjust = 1.01, hjust = 1, angle=90), 
        axis.text.y = element_text(colour = "black", size = 5), 
        legend.text = element_text(size = 5, colour ="black"), 
        legend.position = "right") +  
  scale_fill_manual(values = colours, guide = FALSE)+ 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 20))










filtered_df_2 <- filtered_df_new

filtered_df_2 <- filtered_df_2 %>%
  group_by(mag_species, broad_enzyme) %>%
  mutate(value_be = n()) %>%
  ungroup()

ggplot(filtered_df_2, aes(x = reorder(mag_species, -value_be), y = broad_enzyme, size = value_be, fill = broad_enzyme)) + 
  geom_point(alpha = 0.75, shape = 21) + 
  scale_size_continuous(limits = c(1, 300), range = c(1, 17), breaks = c(10, 50, 100, 300)) + 
  labs(x = "", y = "", size = "# Unique Contigs", fill = "") +theme_classic()+
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 5, angle = 45, vjust = 1.01, hjust = 1), 
        axis.text.y = element_text(colour = "black", size = 5, angle=45), 
        legend.text = element_text(size = 5, colour ="black"), 
        legend.position = "right") +  
  scale_fill_manual(values = colours2, guide = FALSE)

