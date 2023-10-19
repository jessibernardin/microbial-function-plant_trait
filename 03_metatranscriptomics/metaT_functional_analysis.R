#### Microbial function mediates leaf traits in a pitcher plant model system ####
#### Authors: Jessica R. Bernardin, Erica B. Young, Leonora S. Bittleston ####
#### last update : October 18, 2023 ####
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
tpm_500_filt_adjusted <- tpm_500_filt + 0.1 #this is the input for the DREAM model

#read in metadata
meta_all <- read.csv("Exp_1_function_metadata_2021.csv", header=TRUE)
meta_rna <- meta_all %>% filter(nucleic_acid == "RNA")
meta_rna <- meta_rna[,-1]
rownames(meta_rna) <- meta_rna[,1]

#read in results from kofamscan (gene names are different bc there can be several ko per gene)
#this file is very large and could not be added to the github repo
#ko_tsv <- read.table("metaT.kofamscan_detail.tsv", sep = "\t", header = FALSE) #122,695 genes/kos

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

#find differentially abundant genes between treatments using week as a random effect
# Specify parallel processing parameters
param <- SnowParam(4, "SOCK", progressbar=TRUE)

# The variable to be tested must be a fixed effect
form <- ~ treatment + (1|plant_number) 

# estimate weights using linear mixed model of dream (CommA is baseline)
vobjDream2 <- voomWithDreamWeights(norm.mat2, form, meta_rna, BPPARAM=param )
fitmm2 <- dream( vobjDream2, form, meta_rna, BPPARAM=param )
fitmm2 <- eBayes(fitmm2)#83977 genes

# Explicitly specify contrasts using desired baseline (B to C)
# This evaluates CommB - CommC
L1 <- makeContrastsDream(form, meta_rna, contrasts = c(compareB_C = "treatmentCommB - treatmentCommC"))

# evaluate model with contrasts
fitmm3 <- dream(vobjDream2, form, meta_rna, L=L1, BPPARAM=param)
fitmm3 <- eBayes(fitmm3)

#fitmm2<- readRDS("~/Dropbox/2021_Comm_Exp1/meta_RNASeq/dream_fitmm2_allcontigs.RDS")
#fitmm3<- readRDS("~/Dropbox/2021_Comm_Exp1/meta_RNASeq/dream_fitmm3_allcontigs.RDS")

#filter the data by significant differentially abundant genes
#A sample mean with a z-score greater than or equal to the critical value of 1.645 is significant at the 0.05 level.
#can't compare t-statistics here
head(fitmm3$coef, 3) #shows 3 samples
colnames(fitmm2)
colnames(fitmm3)

topB_A <- topTable(fitmm2, coef=c("treatmentCommB"), n=Inf)
topC_A <- topTable(fitmm2, coef=c("treatmentCommC"), n=Inf)
topB_C <- topTable(fitmm3, coef=c("compareB_C"), n=Inf)

###volcano
ggplot(topB_A, aes(x=logFC, y=-log10(P.Value))) + ggtitle("CommB - CommA")+geom_point()
ggplot(topC_A, aes(x=logFC, y=-log10(P.Value))) + ggtitle("CommC - CommA")+geom_point()
ggplot(topB_C, aes(x=logFC, y=-log10(P.Value))) + ggtitle("CommB - CommC")+geom_point()

###mean adjusted
ggplot(topB_A, aes(x=AveExpr, y=logFC)) + ggtitle("CommB - CommA")+geom_point(size=.5, alpha=.5)
ggplot(topC_A, aes(x=AveExpr, y=logFC)) + ggtitle("CommC - CommA")+geom_point()
ggplot(topB_C, aes(x=AveExpr, y=logFC)) + ggtitle("CommB - CommC")+geom_point()

#### B to A ####
# Classify genes into significantly up and down
tt_modified <- topB_A %>% 
  mutate(status=factor(case_when(logFC>0 & adj.P.Val<0.05 ~ "up",
                                 logFC<0 & adj.P.Val<0.05 ~ "down",
                                 TRUE ~ "not.signif"),
                       levels=c("not.signif","up","down")))
#make rowname a column called gene_name
tt_modified$gene_name <- rownames(tt_modified)
tt_modified$delabel <- NA
tt_modified$delabel[tt_modified$status != "not.signif"] <- rownames(tt_modified)[tt_modified$status != "not.signif"]
tt_modified$outliers <- NA

# Create a logical vector based on the specified conditions
outlier_conditions <- abs(tt_modified$logFC) >= 7 | tt_modified$AveExpr > 6

# Update the outliers column using the conditions
tt_modified$outliers[outlier_conditions & tt_modified$status != "not.signif"] <- rownames(tt_modified)[outlier_conditions & tt_modified$status != "not.signif"]

#combine the kegg.f so we can see what the function of the outliers are
#if want L3, or L2, L1, need to merge to KEGG_KO
tt_modified_KO <- left_join(tt_modified, kegg.f, by = "gene_name")
na_rows <- is.na(tt_modified_KO$outliers)

# Set columns 12 to 20 to NA for the identified rows
tt_modified_KO[na_rows, 12:16] <- NA
tt_modified_KO <- tt_modified_KO %>% arrange(status)

# MA-plot
#Fig6B
ggplot(tt_modified_KO, aes(x=AveExpr, y=logFC, color=status)) +
  geom_point(size=.5, alpha = .7) +
  scale_color_manual(values=c("grey","firebrick3",  "navy")) +
  ggtitle("MA plot CommB - CommA") + theme_classic()+ geom_text_repel(aes(label = KO_definition), size = 3, box.padding = 0.5, point.padding = 0.2, max.overlaps = Inf)


########
tt_modified3 <- topC_A %>% 
  mutate(status=factor(case_when(logFC>0 & adj.P.Val<0.05 ~ "up",
                                 logFC<0 & adj.P.Val<0.05 ~ "down",
                                 TRUE ~ "not.signif"),
                       levels=c("not.signif","up","down")))
tt_modified3$gene_name <- rownames(tt_modified3)
tt_modified3$delabel <- NA
tt_modified3$delabel[tt_modified3$status != "not.signif"] <- rownames(tt_modified3)[tt_modified3$status != "not.signif"]
tt_modified3$outliers <- NA
outlier_conditions <- abs(tt_modified3$logFC) >= 7 | tt_modified3$AveExpr > 6
tt_modified3$outliers[outlier_conditions & tt_modified3$status != "not.signif"] <- rownames(tt_modified3)[outlier_conditions & tt_modified3$status != "not.signif"]
tt_modified3_KO <- left_join(tt_modified3, kegg.f, by = "gene_name")
na_rows <- is.na(tt_modified3_KO$outliers)
tt_modified3_KO[na_rows, 12:16] <- NA
tt_modified3_KO <- tt_modified3_KO %>% arrange(status)
ggplot(tt_modified3_KO, aes(x=AveExpr, y=logFC, color=status)) +
  geom_point(size=.5, alpha=.7) +
  scale_color_manual(values=c( "grey", "firebrick","dodgerblue")) +
  ggtitle("MA plot CommC - CommA") + theme_classic()+ geom_text_repel(aes(label = KO_definition), size = 3, box.padding = 0.5, point.padding = 0.2, max.overlaps = Inf)


#### B to C ####
tt_modified4 <- topB_C %>% 
  mutate(status=factor(case_when(logFC>0 & adj.P.Val<0.05 ~ "up",
                                 logFC<0 & adj.P.Val<0.05 ~ "down",
                                 TRUE ~ "not.signif"),
                       levels=c("not.signif","up",  "down")))
tt_modified4$gene_name <- rownames(tt_modified4)
tt_modified4$delabel <- NA
tt_modified4$delabel[tt_modified4$status != "not.signif"] <- rownames(tt_modified4)[tt_modified4$status != "not.signif"]
tt_modified4$outliers <- NA
outlier_conditions <- abs(tt_modified4$logFC) >= 7 | tt_modified4$AveExpr > 6
tt_modified4$outliers[outlier_conditions & tt_modified4$status != "not.signif"] <- rownames(tt_modified4)[outlier_conditions & tt_modified4$status != "not.signif"]
tt_modified4_KO <- left_join(tt_modified4, kegg.f, by = "gene_name")
na_rows <- is.na(tt_modified4_KO$outliers)
tt_modified4_KO[na_rows, 12:16] <- NA
tt_modified4_KO <- tt_modified4_KO %>% arrange(status)

#Fig6A
ggplot(tt_modified4_KO, aes(x=AveExpr, y=logFC, color=status)) +
  geom_point(size=.5, alpha=.7) +
  scale_color_manual(values=c("grey", "firebrick3", "navy")) +
  ggtitle("MA plot CommB - CommC") + theme_classic()+ geom_text_repel(aes(label = KO_definition), size = 3, box.padding = 0.5, point.padding = 0.2, max.overlaps = Inf)

tab1_filtered <- tt_modified4[tt_modified4$adj.P.Val < 0.05, ] 
tab2_filtered <- tt_modified3[tt_modified3$adj.P.Val < 0.05, ] 
tab3_filtered <- tt_modified[tt_modified$adj.P.Val < 0.05, ] 

tab1_filtered$treatment_comp <- "B_C"
tab2_filtered$treatment_comp <- "C_A"
tab3_filtered$treatment_comp <- "B_A"

tab_filtered <- rbind(tab1_filtered, tab2_filtered, tab3_filtered)

#match sig genes with TPM and ko
#make gene name a new column from the rownames
tpm_500_filt$gene_name <- row.names(tpm_500_filt) #community matrix for genes and samples83977
tab_filtered$gene_name <- row.names(tab_filtered) #significant genes according to dream7455

sig_genes <- left_join(tab_filtered, tpm_500_filt,by="gene_name") #7455
sig_genes_description_all <- left_join(sig_genes, KEGG_KO, by="gene_name")
sig_genes_description_all_naomit <- sig_genes_description_all[complete.cases(sig_genes_description_all$L1), ]
sig_genes_dedup <- sig_genes_description_all_naomit %>% distinct(treatment_comp, gene_name, .keep_all = TRUE)

melted_sig_genes <- sig_genes_dedup %>%
  gather(key = Sample, value = TPM, starts_with("RNA"))

melted_sig_genes <- melted_sig_genes %>%
  arrange(desc(abs(logFC)))

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
                     
colcb <- c("#DDAA33", "#BB5566", "#004488")
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
  scale_fill_manual(values=col_treat)


#IAA = "indol|tryptophan|phenylpyruvate|tyramine oxidase"
#pectinesterasem (cytokinins)="pectinesterase | zeatin | isopentenyl|dihydrozeatin"
#ACC=non-proteinogenic amino acid ACC is the precursor and means of long-distance transport of ethylene, a plant hormone associated with growth arrest
#One of the most important genes associated with the increase in plant biomass and stress resistance is acdS, which encodes a 1-aminocyclopropane-1-carboxylate- or ACC-deaminase


#filter whole dataset by all important hormones
filtered_combi_horm <- separated_df %>%
  filter(str_detect(KO_definition, regex("aminocyclopropane|pectinesterase | zeatin | isopentenyl|dihydrozeatin|indol|tryptophan|phenylpyruvate|tyramine oxidase", ignore_case = TRUE)))

#FigS7
filtered_combi_horm %>% filter(week != 3) %>% 
  mutate(log2_TPM = log2(TPM+0.5)) %>%
  ggplot(aes(x = KO_definition, y = log2_TPM, group_by=treatment,fill = treatment)) +
  geom_boxplot(outlier.shape=NA, size=1, alpha=.7) + geom_jitter(aes(color=treatment),position = position_jitterdodge(0.2),cex=1) + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_x_discrete(labels = function(KO_definition) str_wrap(KO_definition, width = 10))+
  facet_wrap(~ Day, nrow = 2)+
  scale_color_manual(values=col_treat)+
  scale_fill_manual(values=col_treat)

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
  scale_fill_manual(values=col_treat)

#### metatranscriptomic function heatmap ####
sig_genes_dedup.no18 <- sig_genes_dedup[,-21]#no reads in this sample

##filter the sig genes by the largest |LFC|
abs_log_fold_change <- abs(sig_genes_dedup.no18$logFC)

#histogram
hist(abs_log_fold_change, breaks = 10, main = "Absolute Log-Fold Change Histogram", xlab = "Absolute Log-Fold Change")
filtered_df <- sig_genes_dedup.no18 %>%
  filter(abs(logFC) > 4) #about half
tran_L3sig <- filtered_df[,c(13:38,46)]
agg_df <- aggregate(. ~ L3, data = tran_L3sig, FUN = sum)
rownames(agg_df) <- agg_df[,1]
agg_df <- agg_df[,-1]

#rarify matrix
colSums(agg_df)
agg_df.t <- t(agg_df)
agg_df.r <- rrarefy(agg_df.t, 4217)
agg_df.rt <- t(agg_df.r)
Bac.Log2.agg_df <- log2(agg_df.rt + 0.5)
Bac.Log2.agg_df <- as.matrix(Bac.Log2.agg_df[order(rowSums(Bac.Log2.agg_df),decreasing = T),])

# Remove columns with names containing "_WK4"
filtered_Bac.Log2.agg_df <- as.data.frame(Bac.Log2.agg_df) %>%
  dplyr::select(-matches("_WK4"))

m01 <- filtered_Bac.Log2.agg_df[,c(7,13:17)]

m06_9 <- filtered_Bac.Log2.agg_df[,-c(7,13:17)]
all <- cbind(m01, m06_9)

#Fig6C
pheatmap(all,
         cluster_cols = FALSE,
         clustering_method = "average", 
         col = colorRampPalette(c("navy", "white", "firebrick3"))(50), border_color=NA, cex=1, fontsize=4)

