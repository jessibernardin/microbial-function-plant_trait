#### Microbial function mediates leaf traits in a pitcher plant model system ####
#### Authors: Jessica R. Bernardin, Erica B. Young, Leonora S. Bittleston ####
#### last update : November 27, 2023 ####
#### Metatranscriptomic Functional Analysis

# Load and install required packages
packages_to_load <- c(
  "tidyr", "dplyr", "compositions", "ggpubr", "variancePartition", "edgeR", "BiocParallel",
  "stringr", "pheatmap", "qiime2R","this.path", "limma", "ggrepel", "vegan"
)

for (i in packages_to_load) { #Installs packages if not yet installed
  if (!require(i, character.only = TRUE)) install.packages(i)
}

setwd(this.path::here())

#kallisto TPM table (gene by sample) from metatranscriptomic analysis
tpm <- read.table("kallisto_contigs_transcript_tpms_all_samples.tsv", header =T, row.names = 1)# 145,861 contigs

#filter all contigs less than 500, remove length column, round, change column names, remove neg control
tpm_500 <- filter(tpm, length >= 500) #83,977 contigs
tpm_500 <- tpm_500[,-1]
tpm_500 <- round(tpm_500, digits = 0)
# Remove the ".aligned" from all column names
new_column_names <- sub("\\.aligned$", "", colnames(tpm_500))
colnames(tpm_500) <- new_column_names
#remove gene name and neg control from tpm 
tpm_500_filt <- subset(tpm_500, select = !grepl("neg", names(tpm_500)))
tpm_500_filt <- tpm_500_filt_adjusted[,-9]
write.csv(tpm_500_filt, "tpm_500_filt.csv")
tpm_500_filt_adjusted <- tpm_500_filt + 0.1 #this is the input for the DREAM model

#read in metadata
meta_all <- read.csv("Exp_1_function_metadata_2021.csv", header=TRUE)
meta_rna <- meta_all %>% filter(nucleic_acid == "RNA")
meta_rna <- meta_rna[,-1]
rownames(meta_rna) <- meta_rna[,1]
meta_rna$treatment <- as.factor(meta_rna$treatment)
#read in results from kofamscan (gene names are different bc there can be several ko per gene)
#this file is very large and could not be added to the github repo
#ko_tsv <- read.table("../../../2021_Comm_Exp1/meta_RNASeq/metaT.kofamscan_detail.tsv", sep = "\t", header = FALSE) #122,695 genes/kos

# simplifying_KofamScan_output KEGG_table
kegg.f <- ko_tsv
kegg.f$V2 <- sapply(strsplit(as.character(kegg.f$V2), "_"), function(x) paste(x[1:2], collapse="_"))
kegg.f<-subset(kegg.f, V6 < 0.001)[,c(2:7)] #e value less than .001
kegg.f<-subset(kegg.f, V5 > 99)[,c(1:6)] #bit score value more than 99

#find the KO with the best match for each contig (using bit score)
kegg.f <- kegg.f %>%
  group_by(V2) %>%
  filter(V5 == max(V5)) %>%
  ungroup()
colnames(kegg.f) <- c("gene_name", "KO", "thrshld", "score", "E_value", "KO_definition")
#82,732 unique contigs

#read in KEGG orthology from website, has pathway info for the KEGGs
KEGG_hier <- read.table("KO_Orthology_ko00001.txt", sep = "\t", header = FALSE) #https://merenlab.org/2018/01/17/importing-ghostkoala-annotations/
names(KEGG_hier) <- c("L1", "L2", "L3", "L4")
KEGG_hier$KO <- word(KEGG_hier$L4, 1) #17,045 KOs
KEGG_KO <- left_join(kegg.f, KEGG_hier, by="KO")

#remove any KOs that don't have functional pathways
KEGG_KO <- KEGG_KO[complete.cases(KEGG_KO$L1), ]

#normalize
dge2 <- DGEList(tpm_500_filt_adjusted) ##the matrix with gene names as row names and counts in columns, 83977 genes
norm.mat2 <- calcNormFactors(dge2)

# Specify parallel processing parameters
param <- SnowParam(4, "SOCK", progressbar=TRUE)

#find differentially abundant genes between treatments using week as a random effect
form <- ~ 0 + treatment + (1|plant_number)

# estimate weights using linear mixed model of dream (no intercept baseline)
vobjDream2 <- voomWithDreamWeights(norm.mat2, form, meta_rna, BPPARAM=param )

#make contrasts
L1 <- makeContrastsDream(form, meta_rna, 
                         contrasts = c(TestCommA = "treatmentCommA - (treatmentCommB + treatmentCommC)/2",
                                       TestCommB = "treatmentCommB - (treatmentCommA + treatmentCommC)/2",
                                       TestCommC = "treatmentCommC - (treatmentCommB + treatmentCommA)/2"))

# evaluate model with contrasts
fitmm3 <- dream(vobjDream2, form, meta_rna, L=L1, BPPARAM=param)
fitmm3 <- eBayes(fitmm3)#83977 genes
saveRDS(fitmm3, "dream_model_ave_contrasts.RDS")
head(fitmm3$coef, 3)

##### Comparing each treatment to the meat of the other treatments
topA <- topTable(fitmm3, coef = "TestCommA", n=Inf)
topB <- topTable(fitmm3, coef = "TestCommB", n=Inf)
topC <- topTable(fitmm3, coef = "TestCommC", n=Inf)

### CommA
ta <- topA %>% 
  mutate(status=factor(case_when(logFC>0 & adj.P.Val<0.05 ~ "up",
                                 logFC<0 & adj.P.Val<0.05 ~ "down",
                                 TRUE ~ "not.signif"),
                       levels=c("not.signif","up","down")))
#make rowname a column called gene_name
ta$gene_name <- rownames(ta)
ta$delabel <- NA
ta$delabel[ta$status != "not.signif"] <- rownames(ta)[ta$status != "not.signif"]
ta$outliers <- NA

# Create a logical vector based on the specified conditions
outlier_conditions <- abs(ta$logFC) >= 6 | ta$AveExpr > 6

# Update the outliers column using the conditions
ta$outliers[outlier_conditions & ta$status != "not.signif"] <- rownames(ta)[outlier_conditions & ta$status != "not.signif"]

#combine the kegg.f so we can see what the function of the outliers are
ta_modified_KO <- left_join(ta, kegg.f, by = "gene_name")
na_rows <- is.na(ta_modified_KO$outliers)

# Set columns 12 to 20 to NA for the identified rows
ta_modified_KO[na_rows, 12:16] <- NA
ta_modified_KO <- ta_modified_KO %>% arrange(status)

#FIG5A
ggplot(ta_modified_KO, aes(x=AveExpr, y=logFC, color=status)) +
  geom_point(size=.5, alpha = .7) +
  scale_color_manual(values=c("grey","firebrick3",  "navy")) +
  ggtitle("MA plot CommA") + theme_classic()+ geom_text_repel(aes(label = KO_definition), size = 3, box.padding = 0.5, point.padding = 0.2, max.overlaps = Inf)

### CommB
tb <- topB %>% 
  mutate(status=factor(case_when(logFC>0 & adj.P.Val<0.05 ~ "up",
                                 logFC<0 & adj.P.Val<0.05 ~ "down",
                                 TRUE ~ "not.signif"),
                       levels=c("not.signif","up","down")))
#make rowname a column called gene_name
tb$gene_name <- rownames(tb)
tb$delabel <- NA
tb$delabel[tb$status != "not.signif"] <- rownames(tb)[tb$status != "not.signif"]
tb$outliers <- NA

# Create a logical vector based on the specified conditions
outlier_conditions <- abs(tb$logFC) >= 6| tb$AveExpr > 6

# Update the outliers column using the conditions
tb$outliers[outlier_conditions & tb$status != "not.signif"] <- rownames(tb)[outlier_conditions & tb$status != "not.signif"]

#combine the kegg.f so we can see what the function of the outliers are
tb_modified_KO <- left_join(tb, kegg.f, by = "gene_name")
na_rows <- is.na(tb_modified_KO$outliers)

# Set columns 12 to 20 to NA for the identified rows
tb_modified_KO[na_rows, 12:16] <- NA
tb_modified_KO <- tb_modified_KO %>% arrange(status)

#FIG5A
ggplot(tb_modified_KO, aes(x=AveExpr, y=logFC, color=status)) +
  geom_point(size=.5, alpha = .7) +
  scale_color_manual(values=c("grey","firebrick3",  "navy")) +
  ggtitle("MA plot CommB - ave others") + theme_classic()+ geom_text_repel(aes(label = KO_definition), size = 3, box.padding = 0.5, point.padding = 0.2, max.overlaps = Inf)

### CommC
tc <- topC %>% 
  mutate(status=factor(case_when(logFC>0 & adj.P.Val<0.05 ~ "up",
                                 logFC<0 & adj.P.Val<0.05 ~ "down",
                                 TRUE ~ "not.signif"),
                       levels=c("not.signif","up","down")))
#make rowname a column called gene_name
tc$gene_name <- rownames(tc)
tc$delabel <- NA
tc$delabel[tc$status != "not.signif"] <- rownames(tc)[tc$status != "not.signif"]
tc$outliers <- NA

# Create a logical vector based on the specified conditions
outlier_conditions <- abs(tc$logFC) >= 6 | tc$AveExpr > 6

# Update the outliers column using the conditions
tc$outliers[outlier_conditions & tc$status != "not.signif"] <- rownames(tc)[outlier_conditions & tc$status != "not.signif"]

#combine the kegg.f so we can see what the function of the outliers are
tc_modified_KO <- left_join(tc, kegg.f, by = "gene_name")
na_rows <- is.na(tc_modified_KO$outliers)

# Set columns 12 to 20 to NA for the identified rows
tc_modified_KO[na_rows, 12:16] <- NA
tc_modified_KO <- tc_modified_KO %>% arrange(status)

#FIG5A
ggplot(tc_modified_KO, aes(x=AveExpr, y=logFC, color=status)) +
  geom_point(size=.5, alpha = .7) +
  scale_color_manual(values=c("grey","firebrick3",  "navy")) +
  ggtitle("MA plot CommC") + theme_classic()+ geom_text_repel(aes(label = KO_definition), size = 3, box.padding = 0.5, point.padding = 0.2, max.overlaps = Inf)

#get DEG
tab1_filtered <- topA[topA$adj.P.Val < 0.05, ] 
tab2_filtered <- topB[topB$adj.P.Val < 0.05, ] 
tab3_filtered <- topC[topC$adj.P.Val < 0.05, ] 

tab1_filtered$treatment_comp <- "CommA"
tab2_filtered$treatment_comp <- "CommB"
tab3_filtered$treatment_comp <- "CommC"

tab1_filtered$gene_name <- rownames(tab1_filtered)
tab2_filtered$gene_name <- rownames(tab2_filtered)
tab3_filtered$gene_name <- rownames(tab3_filtered)

tab_filtered <- rbind(tab1_filtered, tab2_filtered, tab3_filtered)#6567

#make different df of DEG based on logFC
DEG_lfc4 <- tab_filtered %>%
  filter(abs(logFC) > 4)
DEG_lfc4_unique <- unique(DEG_lfc4$gene_name)#675

#match sig genes with TPM and ko
tpm_500_filt$gene_name <- rownames(tpm_500_filt)

tpm_sig_lfc4 <- tpm_500_filt[tpm_500_filt$gene_name %in% DEG_lfc4_unique, ]
tpm_sig_lfc4_KO <- left_join(tpm_sig_lfc4, kegg.f, by="gene_name")

tpm_sig_lfc5 <- tpm_500_filt[tpm_500_filt$gene_name %in% DEG_lfc5_unique, ]
tpm_sig_lfc5_KO <- left_join(tpm_sig_lfc5, kegg.f, by="gene_name")

#prep dataset for looking at specific functions
unique_genes <- unique(tab_filtered$gene_name)#5732

# Filter transcript abundance table based on unique genes
tpm_gene_sig <- tpm_500_filt[tpm_500_filt$gene_name %in% unique_genes, ]
tpm_gene_sig_KO <- left_join(tpm_gene_sig, kegg.f, by="gene_name")

melted_sig_genes <- tpm_gene_sig_KO %>%
  gather(key = Sample, value = TPM, starts_with("RNA"))

separated_df <- melted_sig_genes %>%
  separate(Sample, into = c("sample", "plant", "treatment", "week"), sep = "_")

separated_df$treatment <- ifelse(grepl("M01", separated_df$treatment) & grepl("M01", separated_df$treatment), "CommA",
                              ifelse(grepl("M06", separated_df$treatment) & grepl("M06", separated_df$treatment), "CommB",
                                     ifelse(grepl("M09", separated_df$treatment) & grepl("M09", separated_df$treatment), "CommC",NA)))


separated_df$week <- ifelse(grepl("1", separated_df$week) & grepl("1", separated_df$week), "0",
                                 ifelse(grepl("4", separated_df$week) & grepl("4", separated_df$week), "3",
                                        ifelse(grepl("9", separated_df$week) & grepl("9", separated_df$week), "8",NA)))
separated_df <- separated_df %>% mutate(Day =
                     case_when(week <= 0 ~ 1, 
                               week <= 3 ~ 22,
                               week >= 8 ~ 55))
                     
separated_df <- as.data.frame(separated_df)

#filter whole dataset by chitinase and protease
filtered_whole_chit <- separated_df %>%
  filter(str_detect(KO_definition, regex("chit|aminopeptidase", ignore_case = TRUE)))

filtered_whole_chit$week <- as.factor(filtered_whole_chit$week)
filtered_whole_chit$KO_definition <- as.factor(filtered_whole_chit$KO_definition)
col_treat <- c("#014b7a", "#ffaf49", "#44b7c2")

#Fig3C
filtered_whole_chit %>% filter(week != 3) %>% 
  mutate(log2_TPM = log2(TPM+0.5)) %>%
  ggplot(aes(x = KO_definition, y = log2_TPM, group_by=treatment,fill = treatment)) +
  geom_boxplot(outlier.shape=NA, size=1, alpha=.7) + geom_jitter(aes(color=treatment),position = position_jitterdodge(0.2),cex=1) + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_x_discrete(labels = function(KO_definition) str_wrap(KO_definition, width = 10))+
  facet_wrap(~ Day, nrow = 2)+
  scale_color_manual(values=col_treat)+
  scale_fill_manual(values=col_treat)+ylab("log2(TPM+.5)")


#IAA = "indol|tryptophan|phenylpyruvate|tyramine oxidase"
#pectinesterase (cytokinins)="pectinesterase | zeatin | isopentenyl|dihydrozeatin"
#ACC=non-proteinogenic amino acid ACC is the precursor and means of long-distance transport of ethylene, a plant hormone associated with growth arrest
#One of the most important genes associated with the increase in plant biomass and stress resistance is acdS, which encodes a 1-aminocyclopropane-1-carboxylate- or ACC-deaminase

#filter whole dataset by all important hormones
filtered_combi_horm <- separated_df %>%
  filter(str_detect(KO_definition, regex("aminocyclopropane|pectinesterase | zeatin | isopentenyl|dihydrozeatin|indol|tryptophan|phenylpyruvate|tyramine oxidase", ignore_case = TRUE)))

#FigS8
filtered_combi_horm %>% filter(week != 3) %>% 
  mutate(log2_TPM = log2(TPM+0.5)) %>%
  ggplot(aes(x = KO_definition, y = log2_TPM, group_by=treatment,fill = treatment)) +
  geom_boxplot(outlier.shape=NA, size=1, alpha=.7) + geom_jitter(aes(color=treatment),position = position_jitterdodge(0.2),cex=1) + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_x_discrete(labels = function(KO_definition) str_wrap(KO_definition, width = 10))+
  facet_wrap(~ Day, nrow = 2)+
  scale_color_manual(values=col_treat)+
  scale_fill_manual(values=col_treat)+ylab("log2(TPM+.5)")


#Dominant ecoplate substrates
filtered_eco <- separated_df %>%
  filter(str_detect(KO, regex("K05349|K00677|K05350|K02535|K16363|K02851|K24300|K24301|K00790", ignore_case = TRUE)))

#Fig3F
filtered_eco %>% filter(week != 3) %>% 
  mutate(log2_TPM = log2(TPM+0.5)) %>%
  ggplot(aes(x = KO_definition, y = log2_TPM, group_by=treatment,fill = treatment)) +
  geom_boxplot(outlier.shape=NA, size=1, alpha=.7) + geom_jitter(aes(color=treatment),position = position_jitterdodge(0.2),cex=1) + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_x_discrete(labels = function(KO_definition) str_wrap(KO_definition, width = 10))+
  facet_wrap(~ Day, nrow = 2)+
  scale_color_manual(values=col_treat)+
  scale_fill_manual(values=col_treat)+ylab("log2(TPM+.5)")

#########################

#LOGFC4
tran_L4sig <- tpm_sig_lfc4_KO[,c(1:26,32)]
NAs <- tran_L4sig[!complete.cases(tran_L4sig),] 
tran_L4sig <- tran_L4sig[complete.cases(tran_L4sig),] 
tran_L4sig_agg <- aggregate(. ~ KO_definition, data = tran_L4sig, FUN = sum)
rownames(tran_L4sig_agg) <- tran_L4sig_agg[,1]
tran_L4sig_agg <- tran_L4sig_agg[,-1]
Bac.Log2.agg_df <- log2(tran_L4sig_agg + 0.5)
Bac.Log2.agg_df <- as.matrix(Bac.Log2.agg_df[order(rowSums(Bac.Log2.agg_df),decreasing = T),])
# Remove columns with names containing "_WK4"
filtered_Bac.Log2.agg_df <- as.data.frame(Bac.Log2.agg_df) %>%
  dplyr::select(-matches("_WK4"))
m01 <- filtered_Bac.Log2.agg_df[,c(7,13:17)]
m06_9 <- filtered_Bac.Log2.agg_df[,-c(7,13:17)]
all <- cbind(m01, m06_9)

#FIG5B
pheatmap(all,
         cluster_cols = FALSE,
         clustering_method = "average", 
         col = colorRampPalette(c("navy", "white", "firebrick3"))(50), border_color=NA, cex=1, fontsize=4)

all_norownames<- all
rownames(all_norownames) <- NULL

#FIG5B
pheatmap(all_norownames,
         cluster_cols = FALSE,
         clustering_method = "average", 
         col = colorRampPalette(c("navy", "white", "firebrick3"))(50), border_color=NA, cex=1, fontsize=4)
