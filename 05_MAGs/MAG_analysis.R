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
  geom_jitter(alpha=.3) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                                                  legend.position = "none")+ scale_color_manual(values = c("black", "black"))
