#### Microbial function mediates leaf traits in a pitcher plant model system ####
#### Authors: Jessica R. Bernardin, Erica B. Young, Leonora S. Bittleston ####
#### last update : October 16, 2023 ####
#### Physiological Functions Analysis

#### Load Required Packages ####
packages_to_load <- c(
  "ggplot2", "vegan", "lme4", "tidyverse", "effects", "growthcurver",
  "plyr", "dplyr", "reshape", "reshape2", "ape", "DiagrammeR", "tidyverse",
  "tidybayes", "coefplot", "standardize", "bayesplot", "MCMCvis", "car",
  "patchwork", "ggpubr", "corrr", "ggcorrplot", "factoextra", "MASS",
  "pairwiseAdonis", "plotrix", "gridExtra", "multcompView", "ggeffects", "this.path"
)

# Load and install required packages
for (i in packages_to_load) { #Installs packages if not yet installed
  if (!require(i, character.only = TRUE)) install.packages(i)
}

setwd(this.path::here())

#### Read in Data
#bacterial and plant metadata
data <- read.csv("micro_function_plant_trait_metadata.csv", header = TRUE, check.names = FALSE)

#filter out plants that had damage
data_filt <- data[!(data$plant_number %in% c("25", "44", "51", "9", "29", "NEG")), ]
data_filt$treatment <- as.factor(data_filt$treatment)
data_filt$day <- as.factor(data_filt$day)
data_filt$bacterial_respiration_rate_at_24hr <- as.numeric(data_filt$bacterial_respiration_rate_at_24hr)

#ecoplate data
ecodat <- read.csv("eco_all_2021.csv", header = TRUE)

#original cultures (day0) data
orig_resp <- read.csv("original_cultures_respiration.csv", header = TRUE)
orig_enz_rkn <- read.csv("original_cultures_enzyme_growth_data.csv", header = TRUE)

#### Microbial Functions of cultures prior to inoculation in planta ####
#growth rate original cultures
ggplot(orig_enz_rkn, aes(x = treatment, y = growth_r, color=treatment)) +  geom_boxplot() +
  scale_fill_manual(values=c("#014b7a", "#ffaf49", "#44b7c2"))+ 
  scale_color_manual(values=c("#014b7a", "#ffaf49", "#44b7c2")) +
  theme_classic() + ylab("Growth Rate") + xlab("Original cultures used as inoculant")+
  theme(text = element_text(size = 10)) +geom_jitter()

anova.growth_cult<- aov(growth_r ~ treatment, data= orig_enz_rkn)
summary(anova.growth_cult)

#chitinase original cultures
anova.chit.org<- aov(chitinase_7_29 ~ treatment, data= orig_enz_rkn)
summary(anova.chit.org)

tukey.test_chit <- TukeyHSD(anova.chit.org)
tukey.test_chit
cld <- multcompLetters4(anova.chit.org, tukey.test_chit)
detach(package:plyr)
dt <- orig_enz_rkn %>% na.omit(orig_enz_rkn$chitinase_7_29) %>% dplyr::group_by(treatment) %>%
  summarise(c_mean=mean(chitinase_7_29), sd=sd(chitinase_7_29)) %>%
  arrange(desc(c_mean))
cld <- as.data.frame.list(cld$treatment)
dt$tukey.test_chit <- cld$Letters
color <- c("#ea7000", "#bdbf02", "#39888f")

#Fig1C
ggplot(dt, aes(x = treatment, y = c_mean, fill = treatment)) +
  geom_bar(stat = "identity", position = "dodge", show.legend = FALSE) +
  geom_errorbar(aes(ymax = c_mean + sd, ymin = c_mean - sd),
                position = position_dodge(0.9), width = 0.25, color = "Gray25") +
  scale_fill_manual(values=c("#014b7a", "#ffaf49", "#44b7c2")) +theme_classic()+
  geom_text(aes(label=tukey.test_chit, y = c_mean + sd*1.2), vjust = -.1, size = 10, color = "Gray25", show.legend = FALSE, position = position_dodge(0.9)) + ylab("Chitinase Rate")

#protease original cultures
anova.prot.org<- aov(protease_7_29 ~ treatment, data= orig_enz_rkn)
summary(anova.prot.org)

tukey.test_prot <- TukeyHSD(anova.prot.org)
tukey.test_prot
cld <- multcompLetters4(anova.prot.org, tukey.test_prot)
dt <- orig_enz_rkn %>% na.omit(orig_enz_rkn$protease_7_29) %>% dplyr::group_by(treatment) %>%
  summarise(p_mean=mean(protease_7_29), sd=sd(protease_7_29)) %>%
  arrange(desc(p_mean))
cld <- as.data.frame.list(cld$treatment)
dt$tukey.test_prot <- cld$Letters
color <- c("#ea7000", "#bdbf02", "#39888f")

#Fig1D
ggplot(dt, aes(x = treatment, y = p_mean, fill = treatment)) +
  geom_bar(stat = "identity", position = "dodge", show.legend = FALSE) +
  geom_errorbar(aes(ymax = p_mean + sd, ymin = p_mean - sd),
                position = position_dodge(0.9), width = 0.25, color = "Gray25") +
  scale_fill_manual(values=c("#014b7a", "#ffaf49", "#44b7c2")) +theme_classic()+
  geom_text(aes(label=tukey.test_prot, y = p_mean + sd*1.2), vjust = -.1, size = 10, color = "Gray25", show.legend = FALSE, position = position_dodge(0.9)) + ylab("Protease Rate")

#respiration
# Microresp for each of the 5 weeks of inoculating (CommA, CommB, CommC)
#set as factor
orig_resp$week <- as.factor(orig_resp$week)
orig_resp$rep <- as.factor(orig_resp$rep)
orig_resp$treatment <- as.factor(orig_resp$treatment)

#calculate rate fo all the data
resp_ppm_rate <- orig_resp %>% mutate(rate24 = (orig_resp$T24 - orig_resp$T0) / 24)
ggplot(resp_ppm_rate, aes(x=treatment, y=rate24, color=week)) + geom_jitter()

#remove week 1 data, something not quite right with microresp plate
resp_ppm_rate_filt <- resp_ppm_rate %>% filter(!(week == 1))
resp_ppm_rate_filt <- resp_ppm_rate_filt[,-7]
ggplot(resp_ppm_rate_filt, aes(x=treatment, y=rate24, color=week)) + geom_jitter()

#average reps for each week
resp_ppm_rate_filt_mean <- plyr::ddply(resp_ppm_rate_filt, c("treatment", "week"), summarise,
                                       N    = sum(!is.na(rate24)),
                                       ppm_mean = mean(rate24, na.rm = TRUE),
                                       sd   = sd(rate24, na.rm = TRUE),
                                       se   = sd / sqrt(N))

ggplot(resp_ppm_rate_filt_mean, aes(x=treatment, y=ppm_mean, color=week)) + geom_jitter()

#average each week
resp_ave <- plyr::ddply(resp_ppm_rate_filt_mean, c("treatment"), summarise,
                        N    = sum(!is.na(ppm_mean)),
                        mean = mean(ppm_mean, na.rm = TRUE),
                        sd   = sd(ppm_mean, na.rm = TRUE),
                        se   = sd / sqrt(N))

ggplot(data = resp_ave, aes(x=treatment, y=mean, colour=treatment, group = treatment)) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1) +
  geom_line() + geom_point() + theme_minimal()

##stats###
anova_treatment_resp <- aov(rate24 ~ treatment, data = resp_ppm_rate_filt)
summary(anova_treatment_resp)
tukey.test_resp <- TukeyHSD(anova_treatment_resp)
tukey.test_resp
cld <- multcompLetters4(anova_treatment_resp, tukey.test_resp)
# Table with the mean, the standard deviation and the letters indications significant differences for each treatment
resp_ppm_rate_filt$treatment <- as.factor(resp_ppm_rate_filt$treatment)
levels(resp_ppm_rate_filt$treatment)
dt <- resp_ppm_rate_filt %>% dplyr::group_by(treatment) %>%
  summarise(ppm_mean=mean(rate24), sd=sd(rate24)) %>%
  arrange(desc(ppm_mean))
cld <- as.data.frame.list(cld$treatment)
dt$tukey.test_resp <- cld$Letters
color <- c("#014b7a", "#ffaf49", "#44b7c2")

#FigS4A
ggplot(dt, aes(x = treatment, y = ppm_mean, fill = treatment)) +
  geom_bar(stat = "identity", position = "dodge", show.legend = FALSE) +
  geom_errorbar(aes(ymax = ppm_mean + sd, ymin = ppm_mean - sd),
                position = position_dodge(0.9), width = 0.25, color = "Gray25") +
  scale_fill_manual(values=c("#014b7a", "#ffaf49", "#44b7c2")) +
  geom_text(aes(label=tukey.test_resp, y = ppm_mean + sd*1.2), vjust = -.1, size = 10, color = "Gray25",
            show.legend = FALSE, position = position_dodge(0.9)) + theme_classic() +
  labs(y=expression(CO[2]~Production~(ppm/hr))) + xlab("Treatment Community")


#### Chitinase Activity ####
#Day 1 in plant for 24 hrs
data_day1 <- data_filt %>%
  filter(day == 1, !(treatment %in% c("ACM", "WATER")))

#FigS3A
ggplot(data=data_day1,aes(x = treatment, y = chitinase_rate_uM_per_min, group_by=treatment,fill = treatment)) +
  geom_boxplot(outlier.shape=NA, size=1, alpha=.7) + geom_jitter(aes(color=treatment),position = position_jitterdodge(0.2),cex=1) + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_color_manual(values=c("#014b7a", "#ffaf49", "#44b7c2"))+
  scale_fill_manual(values=c("#014b7a", "#ffaf49", "#44b7c2"))

anova_treatment_chit <- aov(chitinase_rate_uM_per_min ~ treatment, data = data_day1)
summary(anova_treatment_chit)

#chitinase activity day1 through day55
#FigS3C
ggplot(data=data_filt, aes(x=treatment, y=chitinase_rate_uM_per_min, color = treatment))+ 
  facet_wrap(~week, nrow=3) + geom_boxplot() + geom_jitter()+
  ylab("Chitinase Rate") + theme_classic() +
  scale_color_manual(values=c("#5C4033","#014b7a", "#ffaf49", "#44b7c2", "gray"))+
  scale_fill_manual(values=c("#5C4033", "#014b7a", "#ffaf49", "#44b7c2", "gray"))

#Subset to only those samples in metatranscriptomic analysis
col_treat <- c("#014b7a", "#ffaf49", "#44b7c2")

#Fig3A
data_filt %>%  subset(., plant_number %in% c("5", "18", "45", "19", "43", "1", "13", "3", "33")) %>% 
ggplot(aes(x = treatment, y = cumulative_chitinase_uM_min, group_by=treatment,fill = treatment)) +
  geom_boxplot(outlier.shape=NA, size=1, alpha=.7) + geom_jitter(aes(color=treatment),position = position_jitterdodge(0.2),cex=1) + theme_classic() +
  scale_color_manual(values=col_treat)+
  scale_fill_manual(values=col_treat)

#### Protease Activity ####
#Day 1 in plant for 24 hrs
#FigS3B
ggplot(data=data_day1,aes(x = treatment, y = protease_rate_nM_per_min, group_by=treatment,fill = treatment)) +
  geom_boxplot(outlier.shape=NA, size=1, alpha=.7) + geom_jitter(aes(color=treatment),position = position_jitterdodge(0.2),cex=1) + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_color_manual(values=c("#014b7a", "#ffaf49", "#44b7c2"))+
  scale_fill_manual(values=c("#014b7a", "#ffaf49", "#44b7c2"))

anova_treatment_prot <- aov(protease_rate_nM_per_min ~ treatment, data = data_day1)
summary(anova_treatment_prot)

#protease activity day1 through day55
#FigS3D
ggplot(data=data_filt, aes(x=treatment, y=protease_rate_nM_per_min, color = treatment))+ 
  facet_wrap(~week, nrow=3) + geom_boxplot()+ geom_jitter()+
  ylab("Protease Rate") + theme_classic()  +
  scale_color_manual(values=c("#5C4033","#014b7a", "#ffaf49", "#44b7c2", "gray"))+
  scale_fill_manual(values=c("#5C4033", "#014b7a", "#ffaf49", "#44b7c2", "gray"))

#Subset to only those samples in metatranscriptomic analysis
#Fig3B
data_filt %>%  subset(., plant_number %in% c("5", "18", "45", "19", "43", "1", "13", "3", "33")) %>% 
  ggplot(aes(x = treatment, y = cumulative_protease_nM_min, group_by=treatment,fill = treatment)) +
  geom_boxplot(outlier.shape=NA, size=1, alpha=.7) + geom_jitter(aes(color=treatment),position = position_jitterdodge(0.2),cex=1) + theme_classic() +
  scale_color_manual(values=col_treat)+
  scale_fill_manual(values=col_treat)

#### RESPIRATION ####
#FigS4B
ggplot(data=data_filt,aes(x = treatment, y = bacterial_respiration_rate_at_24hr, group_by=treatment,fill = treatment)) +
  geom_boxplot(outlier.shape=NA, size=1, alpha=.7) + geom_jitter(aes(color=treatment),position = position_jitterdodge(0.2),cex=1, size=3) + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_color_manual(values=c("#5C4033", "#014b7a", "#ffaf49", "#44b7c2", "gray"))+
  scale_fill_manual(values=c("#5C4033", "#014b7a", "#ffaf49", "#44b7c2", "gray"))

resp.aov <- aov(bacterial_respiration_rate_at_24hr ~ treatment, data=data_filt)
summary(resp.aov)
tukey.anova_respall <- TukeyHSD(resp.aov)
tukey.anova_respall

#### GROWTH RATE ####
anova_growth <- aov(wk8_r ~ treatment, data= data_filt)
summary(anova_growth) 

#### Plant Traits ####
#FigS6A
ggplot(data=data_filt,aes(x = day, y = pitcher_length_cm, group_by=treatment,fill = treatment)) +
 geom_jitter(aes(color=treatment),position = position_jitterdodge(0.2),cex=1, size=3) + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_color_manual(values=c("#5C4033", "#014b7a", "#ffaf49", "#44b7c2", "gray"))+
  scale_fill_manual(values=c("#5C4033", "#014b7a", "#ffaf49", "#44b7c2", "gray"))+geom_smooth()+
  facet_wrap(~treatment, nrow=1)

#FigS6B
ggplot(data=data_filt,aes(x = treatment, y = photosynthetic_rate_leaf1_μmol_per_m2sec, group_by=treatment,fill = treatment)) +geom_boxplot(outlier.shape=NA, size=1, alpha=.7)+
  geom_jitter(aes(color=treatment),size=3) + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_color_manual(values=c("gray","#5C4033", "#014b7a", "#ffaf49", "#44b7c2"))+
  scale_fill_manual(values=c("gray","#5C4033", "#014b7a", "#ffaf49", "#44b7c2"))

leaf1_photo_anova <- aov(photosynthetic_rate_leaf1_μmol_per_m2sec ~treatment, data = data_filt)
summary(leaf1_photo_anova)

#FigS6C
ggplot(data=data_filt,aes(x = treatment, y = new_leaves, group_by=treatment,fill = treatment)) +geom_boxplot(outlier.shape=NA, size=1, alpha=.7)+
  geom_jitter(aes(color=treatment),size=3) + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_color_manual(values=c("gray","#5C4033", "#014b7a", "#ffaf49", "#44b7c2"))+
  scale_fill_manual(values=c("gray","#5C4033", "#014b7a", "#ffaf49", "#44b7c2"))

#### Correlation between plant traits ####
#subset data to just plant traits
plant_trait <- data_filt[,c(13,14,17,18,19, 23,24,25,26,29,30,31,32,39,40,41)]

plant_trait_naomit <- na.omit(plant_trait)
data_normalized <- scale(plant_trait_naomit)
corr_matrix <- cor(data_normalized)
ggcorrplot(corr_matrix)
data.pca <- princomp(corr_matrix)
summary(data.pca)
data.pca$loadings[, 1:2]
fviz_eig(data.pca, addlabels = TRUE)
fviz_pca_var(data.pca, col.var = "black")
fviz_cos2(data.pca, choice = "var", axes = 1:2)

#Fig2A
fviz_pca_var(data.pca, col.var = "x",
             gradient.cols = c("navy","black", "maroon"),
             repel = TRUE)

#### ECOPLATES ####
# no blanks, experimental controls, or damaged plants
eco_filt <- subset(ecodat, treatment != "WATER" & treatment != "ACM" & 
                 !((day == "55") & (plant %in% c("15", "21", "41", "51", "25", "44", "27", "34", "28", "42", "11", "9", "29"))))

com5 <- eco_filt

com5_meta <- subset(com5, select= c(plant, treatment, day))
com5 <- subset(com5, select= -c(plant, treatment, day))
com5[com5 < 0] <- 0

#run NMDS on ecoplate data, day 1 and day 55
set.seed(123)
nmds5 <- metaMDS(com5, k=2, trymax=500)

#get the data scores from the nmds5
all.scores <- as.data.frame(scores(nmds5, "sites"))
#add metadata to the scores info
all.scores <- cbind(all.scores, com5_meta)
all.scores$day <- as.factor(all.scores$day)

# make a new column that has unique names for week and treatment
all.scores <- all.scores %>% 
  unite(hull.id, treatment, day, remove = FALSE)

hullsall <- all.scores %>%
  group_by(hull.id) %>%
  slice(chull(NMDS1, NMDS2))

#FigS2
ggplot(all.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes(colour = treatment, shape = day)) +
  theme(axis.text.y = element_text(colour = "black", size = 16, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 16), 
        legend.text = element_text(size = 16, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 16), 
        axis.title.x = element_text(face = "bold", size = 16, colour = "black"), 
        legend.title = element_text(size = 16, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "treatment", y = "NMDS2", shape = "day")  + 
  scale_colour_manual(values = c("#014b7a", "#ffaf49", "#44b7c2"))  + aes(fill = factor(hull.id)) + geom_polygon(data = hullsall, alpha = 0.2) + guides(fill=FALSE)+ 
  scale_fill_manual(values = c("#014b7a","#014b7a", "#ffaf49", "#ffaf49", "#44b7c2","#44b7c2"))

#day 1 only
com1 <- subset(ecodat, treatment != "WATER" & treatment != "ACM" & day != "55")
com1_meta <- subset(com1, select= c(plant, treatment, day))
com1 <- subset(com1, select= -c(plant, treatment, day))
com1[com1 < 0] <- 0

#run NMDS on day 1 data only
nmds1 <- metaMDS(com1, k=2, trymax=500)
wk1.scores <- as.data.frame(scores(nmds1, "sites"))
wk1.scores <- cbind(wk1.scores, com1_meta)
wk1.scores$day <- as.factor(wk1.scores$day)

# make a new column that has unique names for week and treatment
wk1.scores <- wk1.scores %>% 
  unite(hull.id, treatment, day, remove = FALSE)

hullswk1 <- wk1.scores %>%
  group_by(hull.id) %>%
  slice(chull(NMDS1, NMDS2))

#Fig1E
ggplot(wk1.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes(colour = treatment, shape = day)) +
  theme(axis.text.y = element_text(colour = "black", size = 16, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 16), 
        legend.text = element_text(size = 16, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 16), 
        axis.title.x = element_text(face = "bold", size = 16, colour = "black"), 
        legend.title = element_text(size = 16, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "treatment", y = "NMDS2", shape = "day")  + 
  scale_colour_manual(values = c("#014b7a", "#ffaf49", "#44b7c2"))  + aes(fill = factor(hull.id)) + geom_polygon(data = hullswk1, alpha = 0.2) + guides(fill=FALSE)+ 
  scale_fill_manual(values = c("#014b7a", "#ffaf49", "#44b7c2"))

## Ecoplate Statistical Analysis
# beta disper to check dispersion day 1 and day 55
ecoall.bd <- betadisper(vegdist(com5), com5_meta$treatment)
ecoall.bd
boxplot(ecoall.bd)
anova(ecoall.bd)# p=0.5602

ecoall.bdday <- betadisper(vegdist(com5), com5_meta$day)
ecoall.bdday
boxplot(ecoall.bdday)
anova(ecoall.bdday)# p=0.06473

# beta disper to check dispersion day 8 and day 55
ecoall.day1 <- betadisper(vegdist(com1), com1_meta$treatment)
ecoall.day1
boxplot(ecoall.day1)
anova(ecoall.day1)# p=0.02602

# PERMANOVA for categorical variables (factors)
#week1
adonis.eco.w1 <- adonis2(com1 ~ treatment, data=wk1.scores, method = "bray")
adonis.eco.w1

#day 8 and day 55
# set number of permutations
perm <- how(nperm = 999)
#specify a random variable (plant).
setBlocks(perm) <- with(all.scores, plant)
adonis <- adonis2(com5 ~ treatment * day, data=all.scores, permutations = perm, method = "bray")
adonis

# pairwise adonis
ecoall.ad.pw <- pairwise.adonis2(com5 ~ treatment * day, data = com5_meta, strata = 'plant')
ecoall.ad.pw

#### ECOPLATES THROUGH TIME ####
#week one doesn't have plant 35 because that pitcher had low volume
dat.1 <- eco_filt[!(eco_filt$plant=="35"),]
dat.1 <- pmax(dat.1, 0)

dat.carbons <- melt(dat.1[, c("plant","treatment", "day", 
                              "Pyruvic.acid.methyl.ester","Tween.40","Tween.80","Alpha.cyclodextrin",
                              "Glycogen","D.Cellobiose","alpha.D.Lactose","Beta.Methyl.D.Glucoside",
                              "D.Xylose","I.Erythritol","D.Mannitol","N.Acetyl.D.Glucosamine",
                              "D.Glucosaminic.acid","Glucose.1.phosphate","D.L.alpha.glycerol.phosphate",
                              "D.Galactonic.acid.gamma.Lactone","D.Galacturonic.acid","X2.Hydroxy.benzoic.acid",
                              "X4.Hydroxy.benzoic.acid","gamma.Hydroxybutyric.acid","Itaconic.acid",
                              "alpha.Ketobutyric.acid","D.Malic.acid","L.Arginine","L.Asparagine",
                              "L.Phenylalanine","L.Serine","L.Threonine","Glycyl.L.glutamic.acid",
                              "Phenylethylamine","Putrescine")], 
                    measure.vars = c("Pyruvic.acid.methyl.ester","Tween.40","Tween.80","Alpha.cyclodextrin",
                                     "Glycogen","D.Cellobiose","alpha.D.Lactose","Beta.Methyl.D.Glucoside",
                                     "D.Xylose","I.Erythritol","D.Mannitol","N.Acetyl.D.Glucosamine",
                                     "D.Glucosaminic.acid","Glucose.1.phosphate","D.L.alpha.glycerol.phosphate",
                                     "D.Galactonic.acid.gamma.Lactone","D.Galacturonic.acid","X2.Hydroxy.benzoic.acid",
                                     "X4.Hydroxy.benzoic.acid","gamma.Hydroxybutyric.acid","Itaconic.acid",
                                     "alpha.Ketobutyric.acid","D.Malic.acid","L.Arginine","L.Asparagine",
                                     "L.Phenylalanine","L.Serine","L.Threonine","Glycyl.L.glutamic.acid",
                                     "Phenylethylamine","Putrescine"))
dat.carbons <- within(dat.carbons, {
  c2 <- as.integer(variable == "Pyruvic.acid.methyl.ester")
  c3<- as.integer(variable == "Tween.40")
  c4<- as.integer(variable == "Tween.80")
  c5<- as.integer(variable == "Alpha.cyclodextrin")
  c6<- as.integer(variable == "Glycogen")
  c7<- as.integer(variable == "D.Cellobiose")
  c8<- as.integer(variable == "alpha.D.Lactose")
  c9<- as.integer(variable == "Beta.Methyl.D.Glucoside")
  c10<- as.integer(variable == "D.Xylose")
  c11<- as.integer(variable == "I.Erythritol")
  c12<- as.integer(variable == "D.Mannitol")
  c13<- as.integer(variable == "N.Acetyl.D.Glucosamine")
  c14<- as.integer(variable == "D.Glucosaminic.acid")
  c15<- as.integer(variable == "Glucose.1.Phosphate")
  c16<- as.integer(variable == "D.L.alpha.glycerol.phosphate")
  c17<- as.integer(variable == "D.Galactonic.acid.gamma.Lactone")
  c18<- as.integer(variable == "D.Galacturonic.acid")
  c19<- as.integer(variable == "X2.Hydroxy.benzoic.acid")
  c20<- as.integer(variable == "X4.Hydroxy.benzoic.acid")
  c21<- as.integer(variable == "gamma.Hydroxybutyric.acid")
  c22<- as.integer(variable == "Itaconic.acid")
  c23<- as.integer(variable == "alpha.Ketobutyric.acid")
  c24<- as.integer(variable == "D.Malic.acid")
  c25<- as.integer(variable == "L.Arginine")
  c26<- as.integer(variable == "L.Asparagine")
  c27<- as.integer(variable == "L.Phenylalanine")
  c28<- as.integer(variable == "L.Serine")
  c29<- as.integer(variable == "L.Threonine")
  c30<- as.integer(variable == "Glycyl.L.glutamic.acid")
  c31<- as.integer(variable == "Phenylethylamine")
  c32<- as.integer(variable == "Putrescine")
})

## VISUALIZATIONS
# Plot average value and se for each substrate within each treatment
# x-axis is the day1 timepoint and y-axis is the day55 timepoint.

#day 1 data
dat.carbons1<-dat.carbons[!(dat.carbons$day=="55"),]
avg.subst1 <- with(dat.carbons1 , aggregate(value, list("Treatment"=as.factor(treatment),"Substrate"=as.factor(variable)), 
                                            mean))
avg.subst1$se <- with(dat.carbons1, aggregate(value, list("Treatment"=as.factor(treatment),"Substrate"=as.factor(variable)),
                                              function(x) std.error(x)))[,3]

#day 55 data
dat.carbons2<-dat.carbons[!(dat.carbons$day=="1"),]
avg.subst2 <- with(dat.carbons2 , aggregate(value, list("Treatment"=as.factor(treatment),"Substrate"=as.factor(variable)), 
                                            mean))
avg.subst2$se <- with(dat.carbons2, aggregate(value, list("Treatment"=as.factor(treatment),"Substrate"=as.factor(variable)),
                                              function(x) std.error(x)))[,3]
avg.subst1$day8 <-avg.subst2$x
avg.subst1$day8.se <- avg.subst2$se



agg <- plyr::ddply(dat.carbons, c("treatment", "day", "variable"), summarise,
             N    = sum(!is.na(value)),
             mean = mean(value, na.rm = TRUE),
             sd   = sd(value, na.rm = TRUE),
             se   = sd / sqrt(N))

pd <- position_dodge(0.4)

pairs <- c(0:30, 0:30, 31:61, 31:61, 62:92, 62:92)

agg <- cbind(agg, pairs)

#lets look at the treatments individually
treatCommA <- subset(agg, treatment == "CommA")
treatCommB <- subset(agg, treatment == "CommB")
treatCommC <- subset(agg, treatment == "CommC")

### correlations
wide <- cbind(treatCommA, treatCommB)
wide <- cbind(wide, treatCommC)

# Define new column names
new_colnames <- c(
  "CommA", "CommA_day", "CommA_variable", "CommA_N", "CommA_mean", "CommA_sd", "CommA_se", "CommA_pairs",
  "CommB", "CommB_day", "CommB_variable", "CommB_N", "CommB_mean", "CommB_sd", "CommB_se", "CommB_pairs",
  "CommC", "CommC_day", "CommC_variable", "CommC_N", "CommC_mean", "CommC_sd", "CommC_se", "CommC_pairs"
)

# Assign the new column names
colnames(wide) <- new_colnames

wide$CommA_day <- as.factor(wide$CommA_day)
wide$CommB_day <- as.factor(wide$CommB_day)
wide$CommC_day <- as.factor(wide$CommC_day)

#Fig3E
##CommBvsCommA
ggplot(data=wide, aes(x=CommA_mean, y=CommB_mean, color = CommA_day, label=CommA_variable))+geom_point()+
  geom_smooth(method="lm",aes(x=CommA_mean,y=CommB_mean))+geom_errorbar(aes(ymax=CommB_mean+CommB_se, ymin=CommB_mean-CommB_se,x=CommA_mean), stat="identity", width=0.01) +
  geom_errorbar(aes(xmax=CommA_mean+CommA_se, xmin=CommA_mean-CommA_se,y=CommB_mean), stat="identity", width=0.01) +
  labs(x = "CommA", y = "CommB") + geom_text(size=3) + theme_classic() + scale_color_manual(values=c("#393D47", "#990099"))

#Fig3D
##CommBvsCommC
ggplot(data=wide, aes(x=CommC_mean, 
                      y=CommB_mean, color = CommC_day, label=CommC_variable))+ geom_point()+
  geom_smooth(method="lm",aes(x=CommC_mean,y=CommB_mean))+
  geom_errorbar(aes(ymax=CommB_mean+CommB_se, ymin=CommB_mean-CommB_se,x=CommC_mean), stat="identity", width=0.01) +
  geom_errorbar(aes(xmax=CommC_mean+CommC_se, xmin=CommC_mean-CommC_se,y=CommB_mean), stat="identity", width=0.01) +
  labs(x = "CommC", y = "CommB") + geom_text(size=3) + theme_classic()+ scale_color_manual(values=c("#393D47", "#990099"))

#### GLMMs ####
## standardize the predictor variables using standardize package
data_filt$cumulative_chitinase_uM_min_scaled <- scale(data_filt$cumulative_chitinase_uM_min)[, 1]
data_filt$cumulative_protease_nM_min_scaled <- scale(data_filt$cumulative_protease_nM_min)[, 1]

#relevel the predictor variables to ACM
data_filt <- data_filt %>% arrange(desc(treatment))
data_filt$treatment <- relevel(data_filt$treatment, ref = "ACM")

#treatment ~ leaf 1 biomass
mbiomass <- brm(drymass_leaf1_grams ~ treatment,data=data_filt, family=Gamma(link="log"), iter = 10000, chains = 4, cores = 4)
saveRDS(mbiomass, file = "brms_mbiomass_ACMlevel.RDS")
mbiomass <- readRDS("brms_mbiomass_ACMlevel.RDS")
summary(mbiomass)
pp_check(mbiomass)
#transformed by exp
fixef(mbiomass)
#CommA
exp(0.15)*exp(-1.29)
#0.319819
#CommB
exp(0.53)*exp(-1.29)
#0.4676664
#CommC
exp(0.17)*exp(-1.29)
#0.3262798
#Water
exp(-0.16)*exp(-1.29)
#0.2345703

mcmc_areas(mbiomass, pars=c("b_Intercept", "b_treatmentCommA","b_treatmentCommB", "b_treatmentCommC", "b_treatmentWATER"),
           prob = 0.80, # 80% intervals
           prob_outer = 1, # 99%
           point_est = "mean") +
  theme_classic() + xaxis_text(color="black")+ yaxis_text(color="black") +
  geom_vline(xintercept=0, linetype="dotted", color="darkgray")


pmbiomass <- ggeffects::ggpredict(mbiomass, "treatment")
colnames(pmbiomass)[1] <- "treatment"
colnames(pmbiomass)[2] <- "drymass_leaf1_grams"

#Fig2B
ggplot(data_filt, aes(x = treatment, y = drymass_leaf1_grams, color=treatment)) + geom_jitter(size=6)+
  scale_color_manual(values=c("#5C4033","#004488", "#ffaf49", "#44b7c2", "gray"))+
  theme(legend.position="none") + geom_point(data=pmbiomass, aes(x=treatment, y=drymass_leaf1_grams), color="black", size=5) +
  geom_linerange(data=pmbiomass,aes(ymin=conf.low, ymax=conf.high), size=2, color="black",
                 position=position_dodge(width = 0.5)) + labs(
                   x = "Treatment", 
                   y = "Pitcher Biomass (g)") + theme_classic()+
  theme(legend.position="none")

posterior_mbiomass <- as.data.frame(mbiomass)

#Directional Effects
CommA <- posterior_mbiomass %>% filter(b_treatmentCommA >0)
nrow(CommA)/nrow(posterior_mbiomass) #the probability of direction 0.8104
CommB <- posterior_mbiomass %>% filter(b_treatmentCommB >0)
nrow(CommB)/nrow(posterior_mbiomass) #the probability of direction .99725
CommC <- posterior_mbiomass %>% filter(b_treatmentCommC >0)
nrow(CommC)/nrow(posterior_mbiomass) #the probability of direction 0.8432
WATER <- posterior_mbiomass %>% filter(b_treatmentWATER <0)
nrow(WATER)/nrow(posterior_mbiomass) #the probability of direction 0.8125

posterior_mbiomass_melt <- posterior_mbiomass[,1:5]
posterior_mbiomass_melt <- reshape2::melt(posterior_mbiomass_melt)

#FigS5A
ggplot(posterior_mbiomass_melt, aes(x = value, y=variable,
                                    fill = stat(x < 0))) +
  stat_halfeye() +
  scale_fill_manual(values=c("#F7D95C", "gray"))+
  geom_vline(aes(xintercept=0), 
             color="black", size=1, linetype="dashed")+
  ylab("Probability density") +
  theme_classic() + guides(fill="none") + xlab("")



## leaf nutrients ~ enzyme activity or treeatment ####

m1 <- brm(leaf1_percent_N ~ cumulative_chitinase_uM_min_scaled + cumulative_protease_nM_min_scaled, data=data_filt, family=Gamma(link="log"), iter = 10000, chains = 4, cores = 4) 
m2 <- brm(leaf1_percent_N ~ treatment, data=data_filt, family=Gamma(link="log"), iter = 10000, chains = 4, cores = 4) 
m3 <- brm(leaf1_percent_C ~ cumulative_chitinase_uM_min_scaled + cumulative_protease_nM_min_scaled, data=data_filt, family=Gamma(link="log"), iter = 10000, chains = 4, cores = 4) 
m4 <- brm(leaf1_percent_C ~ treatment, data=data_filt, family=Gamma(link="log"), iter = 10000, chains = 4, cores = 4) 
m5 <- brm(leaf1_percent_N ~  number_of_pitchers + new_leaves, data=data_filt, family=Gamma(link="log"), iter = 10000, chains = 4, cores = 4)

summary(m1)
summary(m2)
summary(m3)
summary(m4)
summary(m5)

# model 2 estimate effects of treatment on leaf nitrogen
#acm only 1.682028
exp(0.52)
#water 0.9704455
exp(-0.55)*exp(0.52)
#commA 1.568312
exp(-0.07)*exp(0.52)
#commB 1.138828
exp(-0.37)*exp(0.5)
#commC 1.491825
exp(-0.10)*exp(0.5)


# Visualize the credible intervals for our model's parameters:
mcmc_areas(
  m1, 
  pars = c("b_Intercept", "b_cumulative_chitinase_uM_min_scaled", "b_cumulative_protease_nM_min_scaled"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean") + theme_classic() + xaxis_text(color="black")+ yaxis_text(color="black") +
  geom_vline(xintercept=0, linetype="dotted", color="darkgray")+ theme(text = element_text(size = 15))

mcmc_areas(
  m2, 
  pars = c("b_Intercept", "b_treatmentWATER", "b_treatmentCommA", "b_treatmentCommB", "b_treatmentCommC"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean")+ theme_classic() + xaxis_text(color="black")+ yaxis_text(color="black") +
  geom_vline(xintercept=0, linetype="dotted", color="darkgray")

mcmc_areas(
  m3, 
  pars = c("b_Intercept", "b_cumulative_chitinase_uM_min_scaled", "b_cumulative_protease_nM_min_scaled"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean")+ theme_classic() + xaxis_text(color="black")+ yaxis_text(color="black") +
  geom_vline(xintercept=0, linetype="dotted", color="darkgray")

mcmc_areas(
  m4, 
  pars = c("b_Intercept", "b_treatmentWATER", "b_treatmentCommA", "b_treatmentCommB", "b_treatmentCommC"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean")+ theme_classic() + xaxis_text(color="black")+ yaxis_text(color="black") +
  geom_vline(xintercept=0, linetype="dotted", color="darkgray")+ theme(text = element_text(size = 30))

mcmc_areas(
  m5, 
  pars = c("b_Intercept", "b_number_of_pitchers", "b_new_leaves"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean")+ theme_classic() + xaxis_text(color="black")+ yaxis_text(color="black") +
  geom_vline(xintercept=0, linetype="dotted", color="darkgray")

#marginal effect plots for C and N ~ treatment
pm2 <- ggpredict(m2, "treatment")
colnames(pm2)[1] <- "treatment"
colnames(pm2)[2] <- "leaf1_percent_N"

#Fig2C
ggplot(data_filt, aes(x = treatment, y = leaf1_percent_N, color=treatment)) + geom_jitter(size=6, alpha=5)+
  scale_color_manual(values=c("#5C4033","#004488", "#ffaf49", "#44b7c2", "gray"))+
  theme(legend.position="none") + geom_point(data=pm2, aes(x=treatment, y=leaf1_percent_N), color="black", size=5) +
  geom_linerange(data=pm2,aes(ymin=conf.low, ymax=conf.high), size=2, color="black",
                 position=position_dodge(width = 0.5)) + labs(
                   x = "Treatment", 
                   y = "Predicted Leaf Nitrogen (%)") + theme(text = element_text(size =30)) + theme_classic()+
  theme(legend.position="none")

posterior_nit_treat <- as.data.frame(m2)
posterior_nit_treat_melt <- posterior_nit_treat[,1:5]
posterior_nit_treat_melt <- reshape2::melt(posterior_nit_treat_melt)

#FigS5B
ggplot(posterior_nit_treat_melt, aes(x = value, y=variable,
                                    fill = stat(x > 0))) +
  stat_halfeye() +
  scale_fill_manual(values=c("#F7D95C", "gray"))+
  geom_vline(aes(xintercept=0), 
             color="black", size=1, linetype="dashed")+
  ylab("Probability density") +
  theme_classic() + guides(fill="none") + xlab("")

#Check model fit
pp_check(m1)
pp_check(m2)
pp_check(m3)
pp_check(m4)
pp_check(m5)


#positive or negative directional effect
posterior.m1 <- as.data.frame(m1)
protsum_pos <- posterior.m1 %>% filter(b_cumulative_protease_nM_min_scaled >0)
nrow(protsum_pos)/nrow(posterior.m1) 

chitsum_pos <- posterior.m1 %>% filter(b_cumulative_chitinase_uM_min_scaled <0)
nrow(chitsum_pos)/nrow(posterior.m1) 

posterior.m1_melt <- posterior.m1[,1:3]
posterior.m1_melt <- reshape2::melt(posterior.m1_melt)

#FigS5C
ggplot(posterior.m1_melt, aes(x = value, y=variable,
                                     fill = stat(x < 0))) +
  stat_halfeye() +
  scale_fill_manual(values=c("#F7D95C", "gray"))+
  geom_vline(aes(xintercept=0), 
             color="black", size=1, linetype="dashed")+
  ylab("Probability density") +
  theme_classic() + guides(fill="none") + xlab("")

#PROTEASE

#Fig2D
plot(data_filt$leaf1_percent_N ~ data_filt$cumulative_protease_nM_min_scaled, 
     pch=21, bg=c("black"), cex.lab = 2, cex.axis = 1.5,
     xlab="Cumulative protease rate (scaled)", 
     ylab="Leaf Nitrogen (%)", par(mar=c(6,6,4,4)))

# predictive distribution
preds <- posterior_predict(m1, 
                           ndraws = 1000
)
#This "for loop" will plot all of the predictions from our posterior (for this set of values of our IVs), generated by each row (sample) from our posterior.
# The first part (for(i in...)) will repeat the subsequent action for each row [i] in our posterior (so each  MCMC sample):

x <- for(i in 1:nrow(preds)){ 
  curve(exp(      #Inverse Link function for gamma
    posterior.m1$b_Intercept[i]+ #Posterior for intercept
      posterior.m1$b_cumulative_chitinase_uM_min_scaled[i]*0 + 
      posterior.m1$b_cumulative_protease_nM_min_scaled[i]*x),
    add=T, # Add to our existing plot
    col=rgb(0,0,0,0.01)) }

#CHITINASE
# predictive distribution
plot(data_filt$leaf1_percent_N~data_filt$cumulative_chitinase_uM_min_scaled, 
     pch=21, bg=c("black"),cex.lab = 2, cex.axis = 1.5,
     xlab="Cumulative chitinase rate (scaled)",
     ylab="Leaf Nitrogen (%)", par(mar=c(6,6,4,4))) 

x <- for(i in 1:nrow(preds)){ 
  curve(exp(      #Inverse Link function for gamma
    posterior.m1$b_Intercept[i]+ #Posterior for intercept
      posterior.m1$b_cumulative_protease_nM_min_scaled[i]*0 + # holding protease density at mean
      posterior.m1$b_cumulative_chitinase_uM_min_scaled[i]*x),
    add=T, # Add to our existing plot
    col=rgb(0,0,0,0.01)) }
