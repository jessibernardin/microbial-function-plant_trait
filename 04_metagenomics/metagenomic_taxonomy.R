#### Microbial function mediates leaf traits in a pitcher plant model system ####
#### Authors: Jessica R. Bernardin, Erica B. Young, Leonora S. Bittleston ####
#### last update : October 17, 2023 ####
#### Metagenomic Taxonomic Analysis

# Load and install required packages
packages_to_load <- c(
  "tidyverse","janitor","readr", "reshape2", "dplyr", "pheatmap", "qiime2R", "phyloseq", "this.path", "magrittr"
)

for (i in packages_to_load) { #Installs packages if not yet installed
  if (!require(i, character.only = TRUE)) install.packages(i)
}

setwd(this.path::here())

#read in metadata
metadata <- read.csv("metagenomic_metadata.csv", header=TRUE)
metadata <- metadata%>% column_to_rownames("sample_ids")

# read in the sourmash taxonomy results from all samples into a single data frame
sourmash_taxonomy_results <- Sys.glob("sourmash_lineages/*.x.gtdb.with-lineages.csv") %>%
  map_dfr(read_csv, col_types = "ddddddddcccddddcccdc") %>%
  mutate(name = gsub(" .*", "", name))

# We need two tables: a tax table and an "otu" table. 
# The tax table will hold the taxonomic lineages of each of our gather matches.
# To make this, we'll make a table with two columns: the genome match and the lineage of the genome.
# The "otu" table will have the counts of each genome in each sample.
# We'll call this our gather_table.

tax_table <- sourmash_taxonomy_results %>%
  dplyr::select(name, lineage) %>%
  distinct() %>%
  separate(lineage, into = c("domain", "phylum", "class", "order", "family", "genus", "species"), sep = ";") %>%
  column_to_rownames("name")

gather_table <- sourmash_taxonomy_results %>% 
  mutate(n_unique_kmers = (unique_intersect_bp / scaled) * average_abund) %>% # calculate the number of uniquely matched k-mers and multiply it by average k-mer abundance
  dplyr::select(query_name, name, n_unique_kmers) %>% # select only the columns that have information we need
  pivot_wider(id_cols = name, names_from = query_name, values_from = n_unique_kmers) %>% # transform to wide format
  replace(is.na(.), 0) %>% # replace all NAs with 0 
  column_to_rownames("name") # move the metagenome sample name to a rowname

#rearrange
gather_table_t <- as.data.frame(t(gather_table))
a <- gather_table_t[8,]
gather_table_t <- gather_table_t[-8,]

gather_table_t <- gather_table_t %>% add_row(a, .before=13)
rownames(gather_table_t)[rownames(gather_table_t) == "...13"] <- "DNA1_P43_M01_WK1"

b <- gather_table_t[c(13:18),]
gather_table_t <- gather_table_t[-c(13:18),]
gather_table_t <- rbind(b,gather_table_t)
gather_table <- as.matrix(t(gather_table_t))

# Genus level
rownames(tax_table)==rownames(gather_table)
asv16s3tg <- data.frame(genus=tax_table$genus,gather_table)
asv16s3tg$genus[is.na(asv16s3tg$genus)] <- "Unknown"
asv16s3tga <- aggregate(. ~ asv16s3tg$genus, asv16s3tg[,2:ncol(asv16s3tg)], sum) 
row.names(asv16s3tga) <- asv16s3tga[,1]
asv16s3tga <- asv16s3tga[,2:ncol(asv16s3tga)]
asv16s3tgao <- as.matrix(asv16s3tga[order(rowSums(asv16s3tga),decreasing = T),])
asv16s3tgaop <- asv16s3tgao[c(1:20),]
other <- colSums(asv16s3tgao[c(21:nrow(asv16s3tgao)),])
asv16s3tgaopo <- rbind(asv16s3tgaop, other)
data_percentage <- apply(asv16s3tgaopo, 2, function(x){x*100/sum(x,na.rm=T)})
tol21rainbow= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")
#FigS8
barplot(data_percentage,col=tol21rainbow,legend.text=T,axes=F,cex.names = .3,las=2, args.legend = list(x = "topleft", bty = "n", inset=c(-0.05, -0.05), cex=0.4))

tol21rainbow2= c("black",  "#777711",  "#AA7744",  "#44AA77",
                 "#771122",  "#114477",  "#4477AA",  "#117744",  "#AA4455",
                 "#DDDD77",  "#AA4488", "#DDAA77",  "#44AAAA",  "#774411", "#AAAA44",
                 "#77AADD",  "#DD7788",  "#77CCCC",  "#CC99BB", "#771155" ,  "#117777")

#plot using ggplot
ave.gen_gg <- as.data.frame(asv16s3tgaopo)
ave.gen_gg <- adorn_percentages(ave.gen_gg,, 1:18, denominator="col")
other <- ave.gen_gg[21,]
ave.gen_gg <- ave.gen_gg[-21,]
ave.gen_gg <- ave.gen_gg[order(rowSums(ave.gen_gg),decreasing = F),]
ave.gen_gg <- rbind(other, ave.gen_gg)
ave.gen_gg[ "Genus" ] <- rownames(ave.gen_gg)
ave.gen_gg$Genus <- factor(ave.gen_gg$Genus, levels = ave.gen_gg$Genus)
ave.gen_gg.rmelt<- melt( ave.gen_gg, id.vars="Genus", value.name="Relative_Abundance", variable.name="Sample")
ggplot(ave.gen_gg.rmelt, aes( x = Sample, y = Relative_Abundance, fill = Genus)) + geom_bar( position = "fill", stat = "identity", width=1) + theme_classic()+ theme(axis.text.x = element_text(angle = 90, size=5))+
  scale_fill_manual(values = tol21rainbow2)+ guides(fill = guide_legend(ncol = 1))