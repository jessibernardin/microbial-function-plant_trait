#### Microbial function mediates leaf traits in a pitcher plant model system ####
#### Authors: Jessica R. Bernardin, Erica B. Young, Leonora S. Bittleston ####
#### last update : November 27, 2023 ####
#### 16S Metabarcoding Analysis

# List of package names to load
packages_to_load <- c(
  "ggplot2", "vegan", "picante", "dplyr", "lme4", "effects", "coin", "usedist",
  "metacoder", "tibble", "qiime2R", "phyloseq", "RColorBrewer", "decontam", "ANCOMBC",
  "microbiome", "caret", "DT", "pairwiseAdonis", "brms", "reshape2", "gganimate", "gapminder",
  "ggthemes", "readr", "tidyr", "ggpubr", "janitor", "bayesplot", "ggeffects",
  "marginaleffects", "multcompView", "pals", "ggplot2", "dplyr", "gapminder", "ggthemes",
  "readr", "tidyr", "viridis", "reshape", "forcats", "Rlda", "pheatmap", "this.path"
)

# Load and install required packages
for (i in packages_to_load) { #Installs packages if not yet installed
  if (!require(i, character.only = TRUE)) install.packages(i)
}

setwd(this.path::here())

#### Read in asv data, metadata, taxonomy ####
#samples by asv
asv16s <- read_tsv("asv_table_dada2.txt")
asv16s$ASV <- paste("ASV", 1:2633, sep="")
asv_name_match <- asv16s
asv_name_match <- asv_name_match[,c(1,318)]
#write.csv(asv_name_match, "DADA_ASV_Match_newname.csv")

asv16s <- asv16s %>% remove_rownames %>% column_to_rownames(var="ASV")
asv16s <- asv16s[,-1]
summary(rowSums(asv16s))#1
summary(colSums(asv16s))#3406 ASV are row names and samples are columns

#### metadata ####
#this has all the physiology added to the metadata
#samples are rownames
meta <- read_csv("../01_physiology_in_R/micro_function_plant_trait_metadata.csv")
meta <- meta %>% drop_na(Description)
meta$ID <- meta$"#SampleID"
meta <- meta %>% remove_rownames %>% column_to_rownames(var="#SampleID") 

#### read in taxonomy ####
#taxonomy for each ASV
tax.16s <- read_tsv("taxonomy.tsv")
asv_name_match_tax <- merge(tax.16s, asv_name_match, by.x="Feature ID", by.y = "OTU ID")
#write.csv(asv_name_match_tax, "DADA_ASV_Match_newname_tax.csv")
tax.16s <- asv_name_match_tax %>% remove_rownames %>% column_to_rownames(var="ASV")
tax.16s <- tax.16s[,-1]
tax.16s <- tax.16s[order(rownames(tax.16s)), ]
asv16s <- asv16s[order(rownames(asv16s)), ]

row.names(asv16s) == row.names(tax.16s) # sanity check
nrow(asv16s) #=2633 ASVs
nrow(tax.16s) #=2633 ASVs

#### Tree #### 
asv16s.tree <- read.tree("tree.nwk")
asv16s.tree <-root(asv16s.tree, "5eb10ceecb6da365335a61e120997a28")## root with an archaeon

#rename the tree asv to the new names
# Convert the dataframe into a named vector for easier mapping
rename_vec <- setNames(asv_name_match_tax$ASV, asv_name_match_tax$`Feature ID`)

# Rename the ASVs in the phyloseq tree
asv16s.tree$tip.label <- rename_vec[asv16s.tree$tip.label]

#### asv16s.physeq3 = USING all data (raw) ####
row.names(asv16s) == row.names(tax.16s)
asv16s.tax1raw <- data.frame(Feature.ID=row.names(tax.16s),Taxon=tax.16s[,1])
asv16s.tax2raw <- parse_taxonomy(asv16s.tax1raw)
asv16s.tax3raw <- cbind(ASVs=row.names(asv16s.tax2raw),asv16s.tax2raw)

tree.16s.root.raw <-root(asv16s.tree, "ASV1324")## root with an archaeon
new_treeraw <- ape::multi2di(tree.16s.root.raw)

TAX.16sraw <- tax_table(as.matrix(asv16s.tax3raw))
OTU.16sraw <- otu_table(as.matrix(asv16s), taxa_are_rows = TRUE)
SAM.16sraw <- sample_data(meta)
asv16s.physeq3 <-merge_phyloseq(phyloseq(OTU.16sraw),SAM.16sraw,TAX.16sraw,new_treeraw)

#### DECONTAM PACKAGE for identifying contaminants####
#https://benjjneb.github.io/decontam/vignettes/decontam_intro.html
### Prevalance Method (this is raw data, no filtering)
sample_data(asv16s.physeq3)$is.neg <- sample_data(asv16s.physeq3)$plant_number == "NEG"
contamdf.prev <- isContaminant(asv16s.physeq3, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant) ### 128 ASVs are identified as contaminants, 2505 as not
head(which(contamdf.prev$contaminant), n=20)

ggplot(data=contamdf.prev, aes(x=p)) + 
  labs(x = 'decontam-prevalence Score', y='Number of ASVs') + 
  geom_histogram(binwidth=0.02)

###filter the df to samples id'd as contaminants and find out home many samples they are present in and the read abundance for these
cont.prev.true <-subset(contamdf.prev, contaminant=="TRUE")
abund <- asv16s
abund$reads <- rowSums(asv16s)
cont.prev.true.reads <- subset(abund, row.names(abund) %in% rownames(cont.prev.true)) 
cont.true.reads.only <- cont.prev.true.reads %>% dplyr::select(reads)  
cont.prev.true.reads <- cont.prev.true.reads[,-317]

#count up how many samples each ASV is present in
count_non_zero <- function(row) {
  sum(row != 0)
}
cont.prev.true.reads$NonZeroCount <- apply(cont.prev.true.reads[, -1], MARGIN = 1, count_non_zero)

#filter tax to cont asv to check abundance
cont.tax <- subset(tax.16s, row.names(tax.16s) %in% rownames(cont.prev.true.reads)) 
#because these 128 ASV have low frequency and low abundance and were closely associated with the negative controls, they will be removed from the dataset

#### Clean up data ####
# remove chloroplasts
tax.16s2 <- tax.16s[grep("Chloroplast",tax.16s$Taxon, invert = T),]
nrow(tax.16s2) #2624, removed 9 ASVs

# remove mitochondria
tax.16s2 <- tax.16s2[grep("Mitochondria",tax.16s2$Taxon, invert = T),]
nrow(tax.16s2) #2552, removed 72 ASVs

# remove unassigned
tax.16s2 <- tax.16s2[grep("Unassigned",tax.16s2$Taxon, invert = T),]
nrow(tax.16s2) #2544, removed 8 ASVs

#remove contaminants identified above
cont.tax <- cont.tax %>% rownames_to_column(var = "ASV")
tax.16s2 <- tax.16s2 %>% rownames_to_column(var = "ASV")
tax.16s2 <- anti_join(tax.16s2, cont.tax, by="ASV") #2421    3
rownames(tax.16s2) <- tax.16s2$ASV
tax.16s2$ASV<- NULL

#remove the filtered taxa from the ASV table
asv16s_filt <- asv16s[,-c(311:316)]#2633  310
asv16s_filt <- subset(asv16s_filt, row.names(asv16s_filt) %in% row.names(tax.16s2)) #2421  310

#### asv16s.physeq2 = decontaminated data, negative controls removed###
#### also removed experimental controls
row.names(asv16s_filt) == row.names(tax.16s2)
asv16s.tax1raw <- data.frame(Feature.ID=row.names(tax.16s2),Taxon=tax.16s2[,1])
asv16s.tax2raw <- parse_taxonomy(asv16s.tax1raw)
asv16s.tax3raw <- cbind(ASVs=row.names(asv16s.tax2raw),asv16s.tax2raw)
tree.16s <- picante::prune.sample(asv16s_filt, asv16s.tree)# Subset tree to relevant samples
tree.16s.root.raw <-root(asv16s.tree, "ASV1324")## root with an archaeon
new_treeraw <- ape::multi2di(tree.16s.root.raw)
TAX.16sraw <- tax_table(as.matrix(asv16s.tax3raw))
OTU.16sraw <- otu_table(as.matrix(asv16s_filt), taxa_are_rows = TRUE)
SAM.16sraw <- sample_data(meta)
asv16s.physeq2 <-merge_phyloseq(phyloseq(OTU.16sraw),SAM.16sraw,TAX.16sraw,new_treeraw)
asv16s.physeq2 <- subset_samples(asv16s.physeq2, treatment != "ACM")
asv16s.physeq2 <- subset_samples(asv16s.physeq2, treatment != "WATER")

#### Rarify ####
#take out asv less than 10
asv16s_filt[asv16s_filt < 10] <- 0# think about this 
asv16s_filt <- asv16s_filt[as.logical(rowSums(asv16s_filt != 0)), ] # removed 1385 ASVs, 1036 ASV left
tax.16s_filt <- subset(tax.16s2, row.names(tax.16s2) %in% row.names(asv16s_filt))

asv16s_filt.t <- t(asv16s_filt)
min(rowSums(asv16s_filt.t)) #min number of reads per sample = 1669
max(rowSums(asv16s_filt.t)) #max number of reads per sample = 57927
set.seed(1115)
asv16s.rt <- rrarefy(asv16s_filt.t, 1669) ## rarefy at 1669, samples are rows
asv16s.rt <- asv16s.rt[,colSums(asv16s.rt) > 0] ## remove ASVs no longer present
asv16s.rt <- asv16s.rt[order(row.names(asv16s.rt)),] # order samples alphabetically
asv16s.rt <- asv16s.rt[,order(colnames(asv16s.rt))] # order asvs alphabetically
dim(asv16s.rt) #310 samples, 917 ASVs
summary(rowSums(asv16s.rt)) #=1669
summary(colSums(asv16s.rt)) #1

asv16s.r <- t(asv16s.rt) #rows are asvs and columns are samples

#filter taxonomy to match the 917 ASVs
tax.16s2.r <- subset(tax.16s2, row.names(tax.16s2) %in% row.names(asv16s.r)) #filter tax table to match asvs
nrow(asv16s.r)#917
nrow(tax.16s2.r)#917

meta <- meta[order(row.names(meta)),] # order samples alphabetically
meta.r <- subset(meta, row.names(meta) %in% colnames(asv16s.r)) #filter meta to relevant samples (no week 0, just week 1 - 8)

#now asv, tax, meta all the same size
row.names(asv16s.rt) == row.names(meta.r)
tax.16s2.r <- tax.16s2.r[order(row.names(tax.16s2.r)),] # order ASVs alphabetically
asv16s.r <- asv16s.r[order(row.names(asv16s.r)),] # order ASVs alphabetically
row.names(asv16s.r) == row.names(tax.16s2.r)

#### 16S alpha diversity ####
shannon.16s <- diversity(asv16s.r, index="shannon")
ef.16s <- exp(shannon.16s)
colnames(ef.16s)[1] ="effective_species"
summary(ef.16s)
ef.16s <- round(ef.16s)
richness <- colSums(asv16s.r !=0)
md16s <- cbind(meta.r, richness)
md16s <- cbind(md16s, shannon.16s)
md16s <- cbind(md16s, ef.16s)

#relevel the predictor variables to ACM
md16s$treatment <- as.factor(md16s$treatment)
md16s <- md16s %>% arrange(desc(treatment))
md16s$treatment <- relevel(md16s$treatment, ref = "ACM")

#look at effective number of species (ASVs) day 8, 55 for treatments only
md_filt <- md16s %>% subset(.,treatment != "ACM")
md_filt <- md_filt %>% subset(.,treatment != "WATER")
md_filt <- md_filt %>% subset(.,treatment != "WATER")
md_filt <- md_filt %>% subset(.,week %in% c("1", "8"))
md_filt$week <- as.factor(md_filt$week)

#Fig4B
ggplot(data=md_filt,aes(x = treatment, y = effective_species, group_by=treatment,fill = treatment)) +
geom_boxplot(outlier.shape=NA, size=1, alpha=.7) + geom_jitter(aes(color=treatment),cex=3) + theme_classic() +
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
facet_wrap(~ day)+
scale_color_manual(values=c("#014b7a", "#ffaf49", "#44b7c2"))+
scale_fill_manual(values=c("#014b7a", "#ffaf49", "#44b7c2"))

anova.ens_wk18<- aov(effective_species ~ treatment*week, data= md_filt)
summary(anova.ens_wk18)
ths_ens_wk18 <- TukeyHSD(anova.ens_wk18)
ths_ens_wk18

glm1 <- brm(effective_species ~ treatment*week + (1|plant_number), family = poisson, data=md_filt, iter=10000)
pp_check(glm1)
summary(glm1)
posterior_interval(glm1, prob=.95)
mcmc_areas(glm1,
           prob = 0.8, # 80% intervals
           prob_outer = 0.95, # 95%
           point_est = "mean") +
  theme_classic() + xaxis_text(color="black")+ yaxis_text(color="black") +
  geom_vline(xintercept=0, linetype="dotted", color="darkgray")
posterior6 <- mcmc_intervals_data(glm1, 
                                  prob_outer=0.95,
                                  prob=0.5)

posterior6$nonzero <- NA
posterior6$nonzero[posterior6$ll>0 & posterior6$hh>0] <- "nonzero"
posterior6$nonzero[posterior6$ll<0 & posterior6$hh<0] <- "nonzero"
posterior6$nonzero[is.na(posterior6$nonzero)] <- "zero"
posterior6<- posterior6[1:6,]

#FIGS9
ggplot(posterior6, aes(x = parameter,
                       shape=nonzero)) +
  geom_hline(yintercept = 0, linetype = 3, 
             size=1, color = "#b0b5b3") +
  geom_pointrange(aes(ymin = ll, ymax = hh, y = m),
                  position= position_dodge(width=0.75),
                  size = 3/4) +
  scale_color_manual(name="",
                     values = c("grey60", "#484c8d")) +
  scale_shape_manual(values=c(16, 17), 
                     labels=c("95% CI does\nnot contain zero", 
                              "95% CI\ncontains zero"))+
  coord_flip() +theme_bw() + 
  theme(axis.text.y = element_text( size=7), 
        axis.text.x=element_text(size=7),
        axis.title = element_text(size=7), 
        legend.text = element_text(size=7)) +
  xlab(NULL) +
  ylab("Probability of effect of\nTreatment on Pitcher Length (cm)")+
  guides(linetype=FALSE) 

#### MAKE PHYLOSEQ OBJECTS ####
#asv16s.physeq1 = RARIFIED, experimental controls removed
row.names(asv16s.r) == row.names(tax.16s2.r)
r.asv16s.tax1 <- data.frame(Feature.ID=row.names(tax.16s2.r),Taxon=tax.16s2.r[,1])
r.asv16s.tax2 <- parse_taxonomy(r.asv16s.tax1)
r.asv16s.tax3 <- cbind(ASVs=row.names(r.asv16s.tax2),r.asv16s.tax2)

r.tree.16s <- picante::prune.sample(asv16s.rt, asv16s.tree)# Subset tree to relevant samples
r.tree.16s.root <-root(r.tree.16s, "ASV1324")## root with an archaeon
r.new_tree <- ape::multi2di(r.tree.16s.root)

r.TAX.16s <- tax_table(as.matrix(r.asv16s.tax3))
r.OTU.16s <- otu_table(as.matrix(asv16s.r), taxa_are_rows = TRUE)
r.SAM.16s <- sample_data(meta.r)
asv16s.physeq1 <-merge_phyloseq(phyloseq(r.OTU.16s),r.SAM.16s,r.TAX.16s,r.new_tree)
asv16s.physeq1 
asv16s.physeq1 <- subset_samples(asv16s.physeq1, treatment != "ACM")
asv16s.physeq1 <- subset_samples(asv16s.physeq1, treatment != "WATER")
#917 ASVs, 217 samples

#### Unweighted UniFrac ####
ps_treat_only_enzymena_uu <- subset_samples(asv16s.physeq1, chitinase_rate_uM_per_min != "NA")
ps_treat_only_enzymena_uu <- subset_samples(ps_treat_only_enzymena_uu, protease_rate_nM_per_min != "NA")
biomass_uu <- subset_samples(asv16s.physeq1, drymass_leaf1_grams != "NA")

## Calculate unweighted Unifrac distance and run NMDS
uu.dist.16s.treat <- distance(asv16s.physeq1,"uUniFrac")
ps_treat_only_enzymena_uu.dist <- distance(ps_treat_only_enzymena_uu,"uUniFrac")
mass.dist.uu <- distance(biomass_uu,"uUniFrac")

set.seed(123)
uu.nmds.16s.treat <- metaMDS(uu.dist.16s.treat,
                             k = 2, 
                             trymax = 1000,
                             wascores = TRUE)
## Plot NMDS
data.scores2 <- as.data.frame(scores(uu.nmds.16s.treat$points[,1:2]))
data.scores2$ID <- rownames(data.scores2)
md16s.filt <- subset(md16s, treatment != "WATER")
md16s.filt <- subset(md16s.filt, treatment != "ACM")
md16s.drop <- droplevels(md16s.filt)
data_merge2 <- merge(data.scores2, md16s.drop, by = c("ID"))
data_merge2$week <- as.factor(data_merge2$week)
data_merge2$treament <- as.factor(data_merge2$treatment)

#Fig4C
plot(uu.nmds.16s.treat$points[,1:2], xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     col= c("#004488", "#ffaf49", "#44b7c2")[data_merge2$treament],
     pch=19)
legend("topleft", 
       legend=c("CommA","CommB","CommC"),
       col= c("#004488", "#ffaf49", "#44b7c2"),
       pch=19,
       cex=1,
       bty = "n")
ordispider(uu.nmds.16s.treat,groups = data_merge2$treatment, show.groups = "CommA", col = "#004488")
ordispider(uu.nmds.16s.treat,groups = data_merge2$treatment, show.groups = "CommB", col = "#ffaf49")
ordispider(uu.nmds.16s.treat,groups = data_merge2$treatment, show.groups = "CommC", col = "#44b7c2")

#### ggplot ####
ggplot(data_merge2, aes(x = MDS1, y = MDS2)) + 
  geom_point(size = 4, aes( color = treatment))+ scale_color_manual(values=c(c("#004488", "#ffaf49", "#44b7c2"))) +
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", y = "NMDS2", color = "treatment")

#### 16s adonis and mantel tests ####set margin to 2
ad.16s.treat <- adonis2(uu.dist.16s.treat ~ data_merge2$treatment, by="margin")
ad.16s.treat

#### Pairwise adonis ####
pw.ad.16s.treat <- pairwise.adonis2(uu.dist.16s.treat ~ treatment/week, data=data_merge2)
pw.ad.16s.treat 

#### Mantel Test (filtered to three treatments) ####
data_merge_uu <- data_merge2 %>% drop_na(chitinase_rate_uM_per_min)
data_merge_uu <- data_merge_uu %>% drop_na(protease_rate_nM_per_min)

chit.dist.treat <- vegdist(data_merge_uu$chitinase_rate_uM_per_min,method="euclidean", na.rm=TRUE)
prot.dist.treat <- vegdist(data_merge_uu$protease_rate_nM_per_min,method="euclidean", na.rm=TRUE)

treat.asv16s.chit.uu.man <- mantel(ps_treat_only_enzymena_uu.dist,chit.dist.treat, method = "spearman", permutations=999)
treat.asv16s.chit.uu.man 

treat.asv16s.prot.uu.man <- mantel(ps_treat_only_enzymena_uu.dist,prot.dist.treat, method = "spearman", permutations=999)
treat.asv16s.prot.uu.man

data_merge_biomass <- data_merge2 %>% drop_na(drymass_leaf1_grams)
biomass.dist.treat <- vegdist(data_merge_biomass$drymass_leaf1_grams, method="euclidean", na.rm=TRUE)

treat.asv16s.biomass.uu.man <- mantel(mass.dist.uu,biomass.dist.treat, method = "spearman", permutations=999)
treat.asv16s.biomass.uu.man 

#### 16s Barchart ####
asv <- as.data.frame(otu_table(asv16s.physeq1)) #917 ASV and 217 samples
tax <- as.data.frame(tax_table(asv16s.physeq1)) #917 ASV
meta <- as.data.frame(sample_data(asv16s.physeq1)) #917 ASV

tol21rainbow2= c("black",  "#777711",  "#AA7744",  "#44AA77",
                 "#771122",  "#114477",  "#4477AA",  "#117744",  "#AA4455",
                 "#DDDD77",  "#AA4488", "#DDAA77",  "#44AAAA",  "#774411", "#AAAA44",
                 "#77AADD",  "#DD7788",  "#77CCCC",  "#CC99BB", "#771155" ,  "#117777")

##Average biological reps
asv.avg <- asv
colnames(asv.avg) <- sub("\\.R$", "", colnames(asv.avg))

# Convert row names to a new column and melt the data frame
asv.avg.melt <- reshape2::melt(cbind(RowNames = rownames(asv.avg), asv.avg), id.vars = "RowNames")

# Rename the columns
colnames(asv.avg.melt) <- c("ASV", "Sample", "Reads")

asv.melt.sep <- asv.avg.melt %>%
  separate(Sample, into = c("plant", "treatment", "week"), sep = "\\.")

#average the plants within a treatment
asv.avg.melt.sep.ave <- asv.melt.sep %>%
  group_by(ASV,treatment, week) %>%
  mutate(average_reads = mean(Reads, na.rm = TRUE))

asv.avg.melt.sep.ave <- asv.avg.melt.sep.ave[,-c(2,5)]
asv.avg.melt.sep.ave <- unique(asv.avg.melt.sep.ave)
asv.avg.melt.sep.ave$treatment <- factor(asv.avg.melt.sep.ave$treatment, levels = c("M01", "M06", "M09"))
asv.avg.melt.sep.ave <- asv.avg.melt.sep.ave %>%
  arrange(treatment, week)

asv.avg.melt.sep.ave$treatment_week <- paste(asv.avg.melt.sep.ave$treatment, asv.avg.melt.sep.ave$week, sep = ".")
asv.avg.melt.sep.ave <- asv.avg.melt.sep.ave[,-c(2,3)]
pivot_asv.avg <- pivot_wider(asv.avg.melt.sep.ave, names_from = treatment_week, values_from = average_reads)
pivot_asv.avg <- as.data.frame(pivot_asv.avg)
rownames(pivot_asv.avg) <- pivot_asv.avg$ASV
pivot_asv.avg <- pivot_asv.avg[-1]  
rownames(tax)==rownames(pivot_asv.avg)
asv16s3tg <- data.frame(Genus=tax$Genus,pivot_asv.avg)
asv16s3tg$Genus[is.na(asv16s3tg$Genus)] <- "Unknown"
asv16s3tga <- aggregate(. ~ asv16s3tg$Genus, asv16s3tg[,2:ncol(asv16s3tg)], sum)
row.names(asv16s3tga) <- asv16s3tga[,1]
asv16s3tga <- asv16s3tga[,2:ncol(asv16s3tga)]
asv16s3tgao <- as.matrix(asv16s3tga[order(rowSums(asv16s3tga),decreasing = T),])
asv16s3tgaop <- asv16s3tgao[c(1,2,4:21),]
other <- colSums(asv16s3tgao[c(3,22:nrow(asv16s3tgao)),])
asv16s3tgaopo <- rbind(asv16s3tgaop, other)
ave.gen_gg <- as.data.frame(asv16s3tgaopo)
ave.gen_gg <- adorn_percentages(ave.gen_gg,, 1:24, denominator="col")
other <- ave.gen_gg[21,]
ave.gen_gg <- ave.gen_gg[-21,]
ave.gen_gg <- ave.gen_gg[order(rowSums(ave.gen_gg),decreasing = F),]
ave.gen_gg <- rbind(other, ave.gen_gg)
ave.gen_gg[ "Genus" ] <- rownames(ave.gen_gg)
ave.gen_gg$Genus <- factor(ave.gen_gg$Genus, levels = ave.gen_gg$Genus)
ave.gen_gg.rmelt<- reshape2::melt( ave.gen_gg, id.vars="Genus", value.name="Relative_Abundance", variable.name="Sample")

#Fig4A
ggplot(ave.gen_gg.rmelt, aes( x = Sample, y = Relative_Abundance, fill = Genus)) + geom_bar( position = "fill", stat = "identity", width=1) + theme_classic()+ theme(axis.text.x = element_text(angle = 90, size=5))+
  scale_fill_manual(values = tol21rainbow2)+ guides(fill = guide_legend(ncol = 1))+ theme(axis.text.x = element_blank(),axis.ticks = element_blank())

#week 1 only
asv16s3tga.wk1 <- asv16s3tga[,c(1,2,10,18)]
names(asv16s3tga.wk1)[1] <- "Genus"
asv16s3tga <- aggregate(. ~ asv16s3tga.wk1$Genus, asv16s3tga.wk1[,2:ncol(asv16s3tga.wk1)], sum)
row.names(asv16s3tga) <- asv16s3tga[,1]
asv16s3tga <- asv16s3tga[,2:ncol(asv16s3tga)]
asv16s3tgao <- as.matrix(asv16s3tga[order(rowSums(asv16s3tga),decreasing = T),])
asv16s3tgaop <- asv16s3tgao[c(1:3,5:16),]
other <- colSums(asv16s3tgao[c(4,17:nrow(asv16s3tgao)),])
asv16s3tgaopo <- rbind(asv16s3tgaop, other)
tol16rainbow2= c("black", "#44AA77", "#AA4455", "#DDDD77",
                 "#4477AA", "#AAAA44",  "#DD7788",
                 "#44AAAA","#117744", "#CC99BB",
                 "#114477","#77CCCC", "#AA4488",  
                 "#77AADD", "#117777", "#771155")
ave.gen_gg <- as.data.frame(asv16s3tgaopo)
ave.gen_gg <- adorn_percentages(ave.gen_gg,, 1:3, denominator="col")
other <- ave.gen_gg[16,]
ave.gen_gg <- ave.gen_gg[-16,]
ave.gen_gg <- ave.gen_gg[order(rowSums(ave.gen_gg),decreasing = F),]
ave.gen_gg <- rbind(other, ave.gen_gg)
ave.gen_gg[ "Genus" ] <- rownames(ave.gen_gg)
ave.gen_gg$Genus <- factor(ave.gen_gg$Genus, levels = ave.gen_gg$Genus)
ave.gen_gg.rmelt<- reshape2::melt( ave.gen_gg, id.vars="Genus", value.name="Relative_Abundance", variable.name="Sample")

#Fig1F
ggplot(ave.gen_gg.rmelt, aes( x = Sample, y = Relative_Abundance, fill = Genus)) + geom_bar( position = "fill", stat = "identity", width=1) + theme_classic()+ theme(axis.text.x = element_text(angle = 90, size=5))+
  scale_fill_manual(values = tol16rainbow2)+ guides(fill = guide_legend(ncol = 1))

#### Top 16S ASVs, with taxonomy #####
asv16s.top <- data.frame(sort(colSums(asv16s.rt),decreasing = T))
names(asv16s.top) <- "Sequences"
asv16s.top.tax <- subset(tax.16s2, row.names(tax.16s2) %in% row.names(asv16s.top))
asv16s.top.tax <- asv16s.top.tax[match(rownames(asv16s.top), rownames(asv16s.top.tax)), ]
row.names(asv16s.top.tax) == row.names(asv16s.top.tax)
asv16s.top <- data.frame(asv16s.top,asv16s.top.tax)

####look at top 20 asv abundance taxonomy in 16S
sum(asv16s.top$Sequences) #517390 16S reads post clean up
asv16s.top.20 <- asv16s.top[1:20,]
sum(asv16s.top.20$Sequences) #325149 16S reads in top 20

#top 20 asv reads divided by total reads
(325149/517390)*100 #62.84408

Feature.ID <- rownames(asv16s.top.20)
asv16s.top.20 <- cbind(asv16s.top.20, Feature.ID)
asv16s.top.20.tax <- parse_taxonomy(asv16s.top.20)
rownames(asv16s.top.20.tax) == rownames(asv16s.top.20)
asv16s.top.20.tax.abund <- cbind(asv16s.top.20.tax, asv16s.top.20)

#look at top 16S families
asv.16S.mv <- subset(asv16s.top.20.tax.abund, Genus == "Microvirgula")
((sum(asv.16S.mv$Sequences)/325149)*100) #14.21594%
((sum(asv.16S.mv$Sequences)/517390)*100) #8.93388%

asv.16S.st <- subset(asv16s.top.20.tax.abund, Genus == "Stenotrophomonas")
((sum(asv.16S.st$Sequences)/325149)*100) #12.3977%
((sum(asv.16S.st$Sequences)/517390)*100) #7.791221%

asv.16S.pd <- subset(asv16s.top.20.tax.abund, Genus == "Pedobacter")
((sumasv.16S.pd/325149)*100) #11.32527%
((sum(asv.16S.pd$Sequences)/517390)*100) #7.117262%

#### TAX SUMMARY ####
tax <- as.data.frame(tax_table(asv16s.physeq1))
print(length(unique(tax$Kingdom)))#2
print(length(unique(tax$Phylum)))#26
print(length(unique(tax$Class)))#55
print(length(unique(tax$Order)))#121
print(length(unique(tax$Family)))#192
print(length(unique(tax$Genus)))#286

#### ANCOMBC DIFFERENTIAL ABUNDANCE ####
#subset to treatment and only the first and second week
asv16s.physeq1_filt <- asv16s.physeq1 %>% subset_samples(week %in% c("1", "8"))
asv16s.physeq1_WK1 <- asv16s.physeq1_filt %>% subset_samples(week == 1)
asv16s.physeq1_WK8 <- asv16s.physeq1_filt %>% subset_samples(week == 8)

outputwk1_8 <- ancombc2(data = asv16s.physeq1_filt, assay_name = "counts",
                      fix_formula = "treatment", rand_formula = NULL, pairwise = TRUE,
                      p_adj_method = "fdr",
                      prv_cut = 0.3, 
                      group = "treatment",
                      alpha = 0.05,
                      global = TRUE)


outputwk1 <- ancombc2(data = asv16s.physeq1_WK1, assay_name = "counts",
                  fix_formula = "treatment", rand_formula = NULL, pairwise = TRUE,
                  p_adj_method = "fdr",
                  prv_cut = 0.3, 
                  group = "treatment",
                  alpha = 0.05,
                  global = TRUE)

outputwk8 <- ancombc2(data = asv16s.physeq1_WK8, assay_name = "counts",
                      fix_formula = "treatment", rand_formula = NULL, pairwise = TRUE,
                      p_adj_method = "fdr",
                      prv_cut = 0.3, 
                      group = "treatment",
                      alpha = 0.05,
                      global = TRUE)

outputwk1_res <- outputwk1$res
outputwk8_res <- outputwk8$res
outputwk18_res <- outputwk1_8$res

outputwk1_res_sig <- subset(outputwk1_res, `diff_(Intercept)`==T | diff_treatmentCommB==T |
                                              diff_treatmentCommC==T)

dim(outputwk1_res)
dim(outputwk1_res_sig)
#if a taxon has nonzero counts presented in less than 30% (prv_cut=.3)samples, it will not be further analyzed
# 52 ASVs were present in at least 30% of pitcher samples at week 1
# Of these 52, 38 were DA when comparing the three treatments

outputwk8_res_sig <- subset(outputwk8_res, `diff_(Intercept)`==T | diff_treatmentCommB==T |
                              diff_treatmentCommC==T)
dim(outputwk8_res)
dim(outputwk8_res_sig)
# 70 ASVs were present in at least 30% of pitcher samples at week 8
# Of these 70, 19 were DA when comparing the three treatments

outputwk18_res_sig <- subset(outputwk18_res, `diff_(Intercept)`==T | diff_treatmentCommB==T |
                              diff_treatmentCommC==T)
dim(outputwk18_res)
dim(outputwk18_res_sig)
# 61 ASVs were present in at least 30% of pitcher samples at week 1/8
# Of these 61, 61 were DA when comparing the three treatments

#taxonomy
tax <- as.data.frame(tax.16s2.r)
tax <- cbind(Feature.ID=rownames(tax),tax)
tax.p <- parse_taxonomy(tax)
tax.p$ASV <- row.names(tax.p)

#asvs
asv16s.anc <- asv16s.r.t
asv16s.anc.t <- t(asv16s.anc)
asv16s.anc.t <- as.data.frame(asv16s.anc.t)
asv16s.anc.t$ASV <- row.names(asv16s.anc.t)

#ancombc
names(ancombc_wk1_sig)[2] <- "ASV"
names(ancombc_wk8_sig)[2] <- "ASV"
names(ancombc_wk18_sig)[2] <- "ASV"

wk1_ancom_tax <- merge(ancombc_wk1_sig, tax.p, by="ASV")
wk8_ancom_tax <- merge(ancombc_wk8_sig, tax.p, by="ASV")
wk18_ancom_tax <- merge(ancombc_wk18_sig, tax.p, by="ASV")

asv16s.rt.1 <- asv16s.anc.t[, (grep(pattern="*\\.1", colnames(asv16s.anc.t)))]
asv16s.rt.8 <- asv16s.anc.t[, (grep(pattern="*\\.8", colnames(asv16s.anc.t)))]
asv16s.rt.18 <- asv16s.anc.t[, (grep(pattern="*\\.8|*\\.1", colnames(asv16s.anc.t)))]

asv16s.rt.1$ASV <- row.names(asv16s.rt.1)
asv16s.rt.8$ASV <- row.names(asv16s.rt.8)
asv16s.rt.18$ASV <- row.names(asv16s.rt.18)

wk1_ancom_tax <- merge(wk1_ancom_tax, asv16s.rt.1, by="ASV")
wk8_ancom_tax <- merge(wk8_ancom_tax, asv16s.rt.8, by="ASV")
wk18_ancom_tax <- merge(wk18_ancom_tax, asv16s.rt.18, by="ASV")

#combined
wk18_ancom_tax <- wk18_ancom_tax[, -(grep(pattern="ACM", colnames(wk18_ancom_tax)))]
wk18_ancom_tax <- wk18_ancom_tax[, -(grep(pattern="WATER", colnames(wk18_ancom_tax)))]
wk18_ancom_tax.filt <- wk18_ancom_tax[,-c(2:24,27)]

wk1_ancom_tax <- wk1_ancom_tax[, -(grep(pattern="ACM", colnames(wk1_ancom_tax)))]
wk1_ancom_tax <- wk1_ancom_tax[, -(grep(pattern="WATER", colnames(wk1_ancom_tax)))]
wk1_ancom_tax.filt <- wk1_ancom_tax[,-c(2:24,27)]

wk8_ancom_tax <- wk8_ancom_tax[, -(grep(pattern="ACM", colnames(wk8_ancom_tax)))]
wk8_ancom_tax <- wk8_ancom_tax[, -(grep(pattern="WATER", colnames(wk8_ancom_tax)))]
wk8_ancom_tax.filt <- wk8_ancom_tax[,-c(2:24, 27)]

#### 16S Sig ASV HEATMAP ####
#combined wk18_ancom_tax.filt
heat_asv_wk18 <- wk18_ancom_tax.filt
rownames(heat_asv_wk18) <- heat_asv_wk18[,1]
heat_asv_wk18 <- heat_asv_wk18[,-c(1:3)]
Bac.Log2.counts500 <- log2(heat_asv_wk18 + 0.5)


meta.r$week <- as.factor(meta.r$week)
Bac.factorsDS <- dplyr::select(meta.r, treatment, week)
Bac.factorsDS <- Bac.factorsDS %>%
  filter(!(treatment %in% c("WATER", "ACM"))) %>%
  filter(!(week %in% 2:7))

# Reorder treatment
Bac.factorsDS$treatment <- factor(Bac.factorsDS$treatment, levels = c("CommA", "CommB", "CommC"))
DensityCol <- c("#014b7a", "#ffaf49", "#44b7c2")
names(DensityCol) <- levels(Bac.factorsDS$treatment)

# Reorder week
Bac.factorsDS$week <- factor(Bac.factorsDS$week, levels = c("1", "8"))
SpeciesCol <- c("white", "black")
names(SpeciesCol) <- levels(Bac.factorsDS$week)

# Add to a list, where names match those in factors dataframe
AnnColour <- list(
  treatment = DensityCol,
  week = SpeciesCol
)
SampleOrder <- order(Bac.factorsDS$treatment, Bac.factorsDS$week)

#Fig4D
pheatmap(Bac.Log2.counts500,
         cluster_cols = FALSE,
         clustering_method = "average", annotation_colors = AnnColour, annotation_col = Bac.factorsDS, 
         col = colorRampPalette(c("navy", "white", "firebrick3"))(50), border_color=NA)

#### LATENT DIRICHLET ALLOCATION ####
#Hyperparameter associated with the Dirichlet Phi matrix, hyperparameter for membership of species (ie. features) in latent communities)
#degree of species mixing within communities, lower beta minimal mixing and few high abundance species
#increasing beta high species mixing, higher probability of more uniform distribution of species abundances
asv_table <- otu_table(asv16s.physeq1)
asv_table <- as.data.frame(asv_table)
asv16s.treatments.t<- t(asv_table)
asv16s.treatments.t <- as.data.frame(asv16s.treatments.t)
beta <- rep(1,ncol(asv16s.treatments.t)) #samples are rows
meta<- sample_data(asv16s.physeq1)

#site mixing hyperparameter (sometimes set to 0) which would result in no mixing across communities, high species turnover
gamma <- 1 #Hyperparameter associated with the Stick-Breaking prior.
ngibbs <- 20000 #Total number of Gibbs Samples.
n_comm <- 4 #Total number of communities to return
set.seed(125)

# model convergence takes a while
#LDA_Exp1_16S <- rlda.multinomial(data=asv16s.treatments.t, n_community = n_comm, beta=beta, gamma=gamma, n_gibbs=ngibbs, ll_prior=TRUE, display_progress = TRUE)

#The vector of Log-Likelihoods compute for each Gibbs Sample
ll <- LDA_Exp1_16S$logLikelihood[18000:20000]
plot(ll, type="l", xlab="Iterations",ylab="Log(likel.)+log(prior)") # convergence

#community membership in sample units
#theta is the individual probability for each observation (ex: location) belong in each SG (ex: community). It is a matrix with dimension equal n_gibbs by nrow(data) * n_community
theta <- summary(LDA_Exp1_16S, burnin=.95)$Theta.mean
colnames(theta)=paste("SG",1:ncol(theta))

theta %>%melt()%>%
  ggplot(aes(x=Var2,y=value))+geom_boxplot()+
  labs(y="Probability",x="Sub-Group")+theme_bw()+
  theme(axis.text.x=element_text(angle = 45,hjust=1),
        axis.text.y=element_text(size=11),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

theta <- as.data.frame(theta)
meta_filt <- meta[,c(4,7)]
theta_meta <- cbind(meta_filt,theta)
theta_meta_melt <- melt(theta_meta, id.vars=c("treatment","week"))
names(theta_meta_melt) <- c("Treatment","Week","SG","Probability")

# show community mixture in each site
cols=viridis(4)

#Fig5A
theta_meta_melt%>%filter(SG%in%paste("SG",c(1:4)),Treatment%in%c("CommA","CommB", "CommC"),
                         Week%in%c("1","5","8"))%>%
  group_by(SG)%>%
  ggplot(aes(factor(SG), Probability,fill=SG))+
  geom_bar(stat="identity",show.legend=FALSE)+ylim(0,1) +
  scale_fill_manual(values=cols)+ylim(0,1) +
  facet_grid(Treatment~Week)+
  theme_bw()+
  theme(axis.text.x=element_text(angle = 45,hjust=1),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#feature membership in communities
#The individual probability for each variable (ex: Species) belong in each cluster (ex: community). It is a matrix with dimension equal n_gibbs by ncol(data) * n_community
phi <- summary(LDA_Exp1_16S, burnin=.9)$Phi.mean
phi <- as.data.frame(phi)
rownames(phi) <- paste("SG",1:nrow(phi))
phi <- cbind(comm=rownames(phi),phi)
comm <- melt(phi, id.vars="comm")
names(comm) <- c("SG","ASV","Probability")

tax <- tax_table(asv16s.physeq1)
tax <- as.data.frame(tax)
tax$ASV <- row.names(tax)

phi_tax <- left_join(comm, tax, by = c("ASV" = "ASVs"))
grouped_data <- phi_tax %>%
  group_by(SG, Genus) %>%
  summarize(Total_Probability = sum(Probability))
c <- grouped_data %>% group_by(SG)%>%top_n(5,Total_Probability)%>%arrange(SG,Total_Probability)%>%ungroup()%>%mutate(order=-row_number())

#Fig5B
ggplot(c,aes(x=order,y=Total_Probability,fill=factor(SG)))+geom_bar(stat="identity",show.legend=FALSE)+
  facet_wrap(~SG, scales="free")+
  scale_x_continuous(breaks=c$order,labels=c$Genus,expand=c(0,0))+
  scale_fill_manual(values=cols)+
  theme_bw()+
  theme(axis.text.x=element_text(angle = 45,hjust=1, size=6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  xlab("Genus")
