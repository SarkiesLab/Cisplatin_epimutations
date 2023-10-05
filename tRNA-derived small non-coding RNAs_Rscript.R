# Cisplatin analysis

#if needed install packages
#BiocManager::install("DESeq2")

########Call the libraries########
library(tibble)
library(DESeq2)
library("writexl")
library(ggplot2)
library(openxlsx)
library(forcats)
library(datasets)
library(RColorBrewer)
library(FSA)
library(dplyr)
library(ggpubr)
library(nortest)
library(lme4)
library(arm)
library(stringr)
library("survival")
library("survminer")
library(factoextra)
library("gplots")
library(tidyr)
library(VennDiagram)
library(tidyverse)
library(stats)
library(openxlsx)
library(vioplot)
#--------------------------
#####Sup.Fig.2.A.#####

control <- c(199, 177, 184, 217)
with_cisplatin_150 <- c(182, 165, 148, 109)

t_mean <- mean(with_cisplatin_150)
c_mean <- mean(control)

means <- c(c_mean, t_mean)

names <- c("control", "control", "control", "control", "with_cisplatin_150", "with_cisplatin_150", "with_cisplatin_150", "with_cisplatin_150")

table_cis <- c(control, with_cisplatin_150)

table_cis <- data.frame(table_cis, names)
colnames(table_cis) <- c("values", "names")

boxplot(values ~ names, data = table_cis, col = "white", main = "Worm counts: Control vs Cisplatin 150 uM", xlab = "Condition", ylab = "Number of worms")

# Points
stripchart(values ~ names, data = table_cis,
           method = "jitter",
           pch = 19,
           col = 4:1,
           vertical = TRUE,
           add = TRUE)

points(means,col="red",pch=18)

var.test(values ~ names, data = table_cis, 
         alternative = "two.sided")


t.test(control, with_cisplatin_150, alternative = "two.sided", var.equal = T)
wilcox.test(control, with_cisplatin_150, alternative = c("two.sided"))
library(tidyverse)

table_cis$names <- str_replace(table_cis$names, "with_cisplatin_150", "High dose")
table_cis$names <- str_replace(table_cis$names, "control", "Control")

write.xlsx(table_cis,"Data_Sup.Fig.2.A.xlsx")

#--------------------------
#####Work on DNA mutations, Fig.2.A. & Sup.Fig.3#####
#Background F0 CONTROL total#
#Create F0 control total files = mutation already present in the original population#
#Upload data from varscan
setwd("~/Papers/Cisplatin paper/Varscan outputs 0.5")
C1_F0_indels<-read.xlsx("Cisp_A_0.sorted.bam_pileup.vcfF.indel.xlsx")
C2_F0_indels<-read.xlsx("Cisp_B_0.sorted.bam_pileup.vcfF.indel.xlsx")
All_C_indels<-rbind(C2_F0_indels,C1_F0_indels)
C1_F0_snps<-read.xlsx("Cisp_A_0.sorted.bam_pileup.vcfF.snp.xlsx")
C2_F0_snps<-read.xlsx("Cisp_B_0.sorted.bam_pileup.vcfF.snp.xlsx")
All_C_snps<-rbind(C2_F0_snps,C1_F0_snps)
#DNA mutations C1
setwd("~/Papers/Cisplatin paper/Varscan outputs 0.5/C1")
#Indels
C1_intels_list <- list.files(path = "~/Papers/Cisplatin paper/Varscan outputs 0.5/C1", pattern = ".indel.xlsx")
lst <- lapply(C1_intels_list, read.xlsx)
C1_gene0_mut<-lst[[1]]
C1_gene2_mut<-lst[[7]]
C1_gene4_mut<-lst[[9]]
C1_gene6_mut<-lst[[10]]
C1_gene10_mut<-lst[[2]]
C1_gene12_mut<-lst[[3]]
C1_gene14_mut<-lst[[4]]
C1_gene16_mut<-lst[[5]]
C1_gene18_mut<-lst[[6]]
C1_gene20_mut<-lst[[8]]
#Discard background mutations
C1_gene0_mut_wt_bg<-C1_gene0_mut[ !(C1_gene0_mut$Position %in% All_C_indels$Position), ]
C1_gene2_mut_wt_bg<-C1_gene2_mut[ !(C1_gene2_mut$Position %in% All_C_indels$Position), ]
C1_gene4_mut_wt_bg<-C1_gene4_mut[ !(C1_gene4_mut$Position %in% All_C_indels$Position), ]
C1_gene6_mut_wt_bg<-C1_gene6_mut[ !(C1_gene6_mut$Position %in% All_C_indels$Position), ]
C1_gene10_mut_wt_bg<-C1_gene10_mut[ !(C1_gene10_mut$Position %in% All_C_indels$Position), ]
C1_gene12_mut_wt_bg<-C1_gene12_mut[ !(C1_gene12_mut$Position %in% All_C_indels$Position), ]
C1_gene14_mut_wt_bg<-C1_gene14_mut[ !(C1_gene14_mut$Position %in% All_C_indels$Position), ]
C1_gene16_mut_wt_bg<-C1_gene16_mut[ !(C1_gene16_mut$Position %in% All_C_indels$Position), ]
C1_gene18_mut_wt_bg<-C1_gene18_mut[ !(C1_gene18_mut$Position %in% All_C_indels$Position), ]
C1_gene20_mut_wt_bg<-C1_gene20_mut[ !(C1_gene20_mut$Position %in% All_C_indels$Position), ]
#C1_mutations_from_0_to_20
df_mergebg <- merge(C1_gene0_mut_wt_bg,C1_gene2_mut_wt_bg,by="Position")
df_mergebg <- merge(df_mergebg,C1_gene4_mut_wt_bg,by="Position")
df_mergebg <- merge(df_mergebg,C1_gene6_mut_wt_bg,by="Position")
df_mergebg <- merge(df_mergebg,C1_gene10_mut_wt_bg,by="Position")
df_mergebg <- merge(df_mergebg,C1_gene12_mut_wt_bg,by="Position")
df_mergebg <- merge(df_mergebg,C1_gene14_mut_wt_bg,by="Position")
df_mergebg <- merge(df_mergebg,C1_gene16_mut_wt_bg,by="Position")
df_mergebg <- merge(df_mergebg,C1_gene18_mut_wt_bg,by="Position")
C1_indels_from_0_to_20 <- merge(df_mergebg,C1_gene20_mut_wt_bg,by="Position")
write.xlsx(C1_indels_from_0_to_20,"C1_indels_from_0_to_20.xlsx")
#Discard F0 mutations
C1_gene2_mut_wt_gene0<-C1_gene2_mut_wt_bg[ !(C1_gene2_mut_wt_bg$Position %in% C1_gene0_mut_wt_bg$Position), ]
C1_gene4_mut_wt_gene0<-C1_gene4_mut_wt_bg[ !(C1_gene4_mut_wt_bg$Position %in% C1_gene0_mut_wt_bg$Position), ]
C1_gene6_mut_wt_gene0<-C1_gene6_mut_wt_bg[ !(C1_gene6_mut_wt_bg$Position %in% C1_gene0_mut_wt_bg$Position), ]
C1_gene10_mut_wt_gene0<-C1_gene10_mut_wt_bg[ !(C1_gene10_mut_wt_bg$Position %in% C1_gene0_mut_wt_bg$Position), ]
C1_gene12_mut_wt_gene0<-C1_gene12_mut_wt_bg[ !(C1_gene12_mut_wt_bg$Position %in% C1_gene0_mut_wt_bg$Position), ]
C1_gene14_mut_wt_gene0<-C1_gene14_mut_wt_bg[ !(C1_gene14_mut_wt_bg$Position %in% C1_gene0_mut_wt_bg$Position), ]
C1_gene16_mut_wt_gene0<-C1_gene16_mut_wt_bg[ !(C1_gene16_mut_wt_bg$Position %in% C1_gene0_mut_wt_bg$Position), ]
C1_gene18_mut_wt_gene0<-C1_gene18_mut_wt_bg[ !(C1_gene18_mut_wt_bg$Position %in% C1_gene0_mut_wt_bg$Position), ]
C1_gene20_mut_wt_gene0<-C1_gene20_mut_wt_bg[ !(C1_gene20_mut_wt_bg$Position %in% C1_gene0_mut_wt_bg$Position), ]
#C1_indels_from_2_to_20
df_merge <- merge(C1_gene2_mut_wt_gene0,C1_gene4_mut_wt_gene0,by="Position")
df_merge <- merge(df_merge,C1_gene6_mut_wt_gene0,by="Position")
df_merge <- merge(df_merge,C1_gene10_mut_wt_gene0,by="Position")
df_merge <- merge(df_merge,C1_gene12_mut_wt_gene0,by="Position")
df_merge <- merge(df_merge,C1_gene14_mut_wt_gene0,by="Position")
df_merge <- merge(df_merge,C1_gene16_mut_wt_gene0,by="Position")
df_merge <- merge(df_merge,C1_gene18_mut_wt_gene0,by="Position")
C1_indels_from_2_to_20 <- merge(df_merge,C1_gene20_mut_wt_gene0,by="Position")
write.xlsx(C1_indels_from_2_to_20,"C1_indels_from_2_to_20.xlsx")
#Discard F2 mutations
C1_gene4_mut_wt_gene2<-C1_gene4_mut_wt_gene0[ !(C1_gene4_mut_wt_gene0$Position %in% C1_gene2_mut_wt_gene0$Position), ]
C1_gene6_mut_wt_gene2<-C1_gene6_mut_wt_gene0[ !(C1_gene6_mut_wt_gene0$Position %in% C1_gene2_mut_wt_gene0$Position), ]
C1_gene10_mut_wt_gene2<-C1_gene10_mut_wt_gene0[ !(C1_gene10_mut_wt_gene0$Position %in% C1_gene2_mut_wt_gene0$Position), ]
C1_gene12_mut_wt_gene2<-C1_gene12_mut_wt_gene0[ !(C1_gene12_mut_wt_gene0$Position %in% C1_gene2_mut_wt_gene0$Position), ]
C1_gene14_mut_wt_gene2<-C1_gene14_mut_wt_gene0[ !(C1_gene14_mut_wt_gene0$Position %in% C1_gene2_mut_wt_gene0$Position), ]
C1_gene16_mut_wt_gene2<-C1_gene16_mut_wt_gene0[ !(C1_gene16_mut_wt_gene0$Position %in% C1_gene2_mut_wt_gene0$Position), ]
C1_gene18_mut_wt_gene2<-C1_gene18_mut_wt_gene0[ !(C1_gene18_mut_wt_gene0$Position %in% C1_gene2_mut_wt_gene0$Position), ]
C1_gene20_mut_wt_gene2<-C1_gene20_mut_wt_gene0[ !(C1_gene20_mut_wt_gene0$Position %in% C1_gene2_mut_wt_gene0$Position), ]
#C1_mutations_from_4_to_20
df_merge4 <- merge(C1_gene4_mut_wt_gene2,C1_gene6_mut_wt_gene2,by="Position")
df_merge4 <- merge(df_merge4,C1_gene10_mut_wt_gene2,by="Position")
df_merge4 <- merge(df_merge4,C1_gene12_mut_wt_gene2,by="Position")
df_merge4 <- merge(df_merge4,C1_gene14_mut_wt_gene2,by="Position")
df_merge4 <- merge(df_merge4,C1_gene16_mut_wt_gene2,by="Position")
df_merge4 <- merge(df_merge4,C1_gene18_mut_wt_gene2,by="Position")
C1_indels_from_4_to_20 <- merge(df_merge4,C1_gene20_mut_wt_gene2,by="Position")
write.xlsx(C1_indels_from_4_to_20,"C1_indels_from_4_to_20.xlsx")
#Discard F4 mutations
C1_gene6_mut_wt_gene4<-C1_gene6_mut_wt_gene2[ !(C1_gene6_mut_wt_gene2$Position %in% C1_gene4_mut_wt_gene2$Position), ]
C1_gene10_mut_wt_gene4<-C1_gene10_mut_wt_gene2[ !(C1_gene10_mut_wt_gene2$Position %in% C1_gene4_mut_wt_gene2$Position), ]
C1_gene12_mut_wt_gene4<-C1_gene12_mut_wt_gene2[ !(C1_gene12_mut_wt_gene2$Position %in% C1_gene4_mut_wt_gene2$Position), ]
C1_gene14_mut_wt_gene4<-C1_gene14_mut_wt_gene2[ !(C1_gene14_mut_wt_gene2$Position %in% C1_gene4_mut_wt_gene2$Position), ]
C1_gene16_mut_wt_gene4<-C1_gene16_mut_wt_gene2[ !(C1_gene16_mut_wt_gene2$Position %in% C1_gene4_mut_wt_gene2$Position), ]
C1_gene18_mut_wt_gene4<-C1_gene18_mut_wt_gene2[ !(C1_gene18_mut_wt_gene2$Position %in% C1_gene4_mut_wt_gene2$Position), ]
C1_gene20_mut_wt_gene4<-C1_gene20_mut_wt_gene2[ !(C1_gene20_mut_wt_gene2$Position %in% C1_gene4_mut_wt_gene2$Position), ]
#C1_mutations_from_6_to_20
df_merge6 <- merge(C1_gene6_mut_wt_gene4,C1_gene10_mut_wt_gene4,by="Position")
df_merge6 <- merge(df_merge6,C1_gene12_mut_wt_gene4,by="Position")
df_merge6 <- merge(df_merge6,C1_gene14_mut_wt_gene4,by="Position")
df_merge6 <- merge(df_merge6,C1_gene16_mut_wt_gene4,by="Position")
df_merge6 <- merge(df_merge6,C1_gene18_mut_wt_gene4,by="Position")
C1_indels_from_6_to_20 <- merge(df_merge6,C1_gene20_mut_wt_gene4,by="Position")
write.xlsx(C1_indels_from_6_to_20,"C1_indels_from_6_to_20.xlsx")
#Discard F6 mutations
C1_gene10_mut_wt_gene6<-C1_gene10_mut_wt_gene4[ !(C1_gene10_mut_wt_gene4$Position %in% C1_gene6_mut_wt_gene4$Position), ]
C1_gene12_mut_wt_gene6<-C1_gene12_mut_wt_gene4[ !(C1_gene12_mut_wt_gene4$Position %in% C1_gene6_mut_wt_gene4$Position), ]
C1_gene14_mut_wt_gene6<-C1_gene14_mut_wt_gene4[ !(C1_gene14_mut_wt_gene4$Position %in% C1_gene6_mut_wt_gene4$Position), ]
C1_gene16_mut_wt_gene6<-C1_gene16_mut_wt_gene4[ !(C1_gene16_mut_wt_gene4$Position %in% C1_gene6_mut_wt_gene4$Position), ]
C1_gene18_mut_wt_gene6<-C1_gene18_mut_wt_gene4[ !(C1_gene18_mut_wt_gene4$Position %in% C1_gene6_mut_wt_gene4$Position), ]
C1_gene20_mut_wt_gene6<-C1_gene20_mut_wt_gene4[ !(C1_gene20_mut_wt_gene4$Position %in% C1_gene6_mut_wt_gene4$Position), ]
#C1_mutations_from_10_to_20
df_merge10 <- merge(C1_gene10_mut_wt_gene6,C1_gene12_mut_wt_gene6,by="Position")
df_merge10 <- merge(df_merge10,C1_gene14_mut_wt_gene6,by="Position")
df_merge10 <- merge(df_merge10,C1_gene16_mut_wt_gene6,by="Position")
df_merge10 <- merge(df_merge10,C1_gene18_mut_wt_gene6,by="Position")
C1_indels_from_10_to_20 <- merge(df_merge10,C1_gene20_mut_wt_gene6,by="Position")
write.xlsx(C1_indels_from_10_to_20,"C1_indels_from_10_to_20.xlsx")
#Discard F10 mutations
C1_gene12_mut_wt_gene10<-C1_gene12_mut_wt_gene6[ !(C1_gene12_mut_wt_gene6$Position %in% C1_gene10_mut_wt_gene6$Position), ]
C1_gene14_mut_wt_gene10<-C1_gene14_mut_wt_gene6[ !(C1_gene14_mut_wt_gene6$Position %in% C1_gene10_mut_wt_gene6$Position), ]
C1_gene16_mut_wt_gene10<-C1_gene16_mut_wt_gene6[ !(C1_gene16_mut_wt_gene6$Position %in% C1_gene10_mut_wt_gene6$Position), ]
C1_gene18_mut_wt_gene10<-C1_gene18_mut_wt_gene6[ !(C1_gene18_mut_wt_gene6$Position %in% C1_gene10_mut_wt_gene6$Position), ]
C1_gene20_mut_wt_gene10<-C1_gene20_mut_wt_gene6[ !(C1_gene20_mut_wt_gene6$Position %in% C1_gene10_mut_wt_gene6$Position), ]
#C1_mutations_from_12_to_20
df_merge12 <- merge(C1_gene12_mut_wt_gene10,C1_gene14_mut_wt_gene10,by="Position")
df_merge12 <- merge(df_merge12,C1_gene16_mut_wt_gene10,by="Position")
df_merge12 <- merge(df_merge12,C1_gene18_mut_wt_gene10,by="Position")
C1_indels_from_12_to_20 <- merge(df_merge12,C1_gene20_mut_wt_gene10,by="Position")
write.xlsx(C1_indels_from_12_to_20,"C1_indels_from_12_to_20.xlsx")
#Discard F12 mutations
C1_gene14_mut_wt_gene12<-C1_gene14_mut_wt_gene10[ !(C1_gene14_mut_wt_gene10$Position %in% C1_gene12_mut_wt_gene10$Position), ]
C1_gene16_mut_wt_gene12<-C1_gene16_mut_wt_gene10[ !(C1_gene16_mut_wt_gene10$Position %in% C1_gene12_mut_wt_gene10$Position), ]
C1_gene18_mut_wt_gene12<-C1_gene18_mut_wt_gene10[ !(C1_gene18_mut_wt_gene10$Position %in% C1_gene12_mut_wt_gene10$Position), ]
C1_gene20_mut_wt_gene12<-C1_gene20_mut_wt_gene10[ !(C1_gene20_mut_wt_gene10$Position %in% C1_gene12_mut_wt_gene10$Position), ]
#C1_mutations_from_14_to_20
df_merge14 <- merge(C1_gene14_mut_wt_gene12,C1_gene16_mut_wt_gene12,by="Position")
df_merge14 <- merge(df_merge14,C1_gene18_mut_wt_gene12,by="Position")
C1_indels_from_14_to_20 <- merge(df_merge14,C1_gene20_mut_wt_gene12,by="Position")
write.xlsx(C1_indels_from_14_to_20,"C1_indels_from_14_to_20.xlsx")
#Discard F14 mutations
C1_gene16_mut_wt_gene14<-C1_gene16_mut_wt_gene12[ !(C1_gene16_mut_wt_gene12$Position %in% C1_gene14_mut_wt_gene12$Position), ]
C1_gene18_mut_wt_gene14<-C1_gene18_mut_wt_gene12[ !(C1_gene18_mut_wt_gene12$Position %in% C1_gene14_mut_wt_gene12$Position), ]
C1_gene20_mut_wt_gene14<-C1_gene20_mut_wt_gene12[ !(C1_gene20_mut_wt_gene12$Position %in% C1_gene14_mut_wt_gene12$Position), ]
#C1_mutations_from_16_to_20
df_merge16 <- merge(C1_gene16_mut_wt_gene14,C1_gene18_mut_wt_gene14,by="Position")
C1_indels_from_16_to_20 <- merge(df_merge16,C1_gene20_mut_wt_gene14,by="Position")
write.xlsx(C1_indels_from_16_to_20,"C1_indels_from_16_to_20.xlsx")
#Discard F16 mutations
C1_gene18_mut_wt_gene16<-C1_gene18_mut_wt_gene14[ !(C1_gene18_mut_wt_gene14$Position %in% C1_gene16_mut_wt_gene14$Position), ]
C1_gene20_mut_wt_gene16<-C1_gene20_mut_wt_gene14[ !(C1_gene20_mut_wt_gene14$Position %in% C1_gene16_mut_wt_gene14$Position), ]
#C1_mutations_from_18_to_20
C1_indels_from_18_to_20 <- merge(C1_gene18_mut_wt_gene16,C1_gene20_mut_wt_gene16,by="Position")
write.xlsx(C1_indels_from_18_to_20,"C1_indels_from_18_to_20.xlsx")
#Discard F18 mutations
C1_gene20_mut_wt_gene18<-C1_gene20_mut_wt_gene16[ !(C1_gene20_mut_wt_gene16$Position %in% C1_gene18_mut_wt_gene16$Position), ]
#C1_mutations_from_20_to_20
C1_indels_from_20_to_20 <- C1_gene20_mut_wt_gene18
write.xlsx(C1_indels_from_20_to_20,"C1_indels_from_20_to_20.xlsx")
######Snps
C1_snps_list <- list.files(path = "~/Papers/Cisplatin paper/Varscan outputs 0.5/C1", pattern = ".snp.xlsx")
lst <- lapply(C1_snps_list, read.xlsx)
C1_gene0_snps<-lst[[1]]
C1_gene2_snps<-lst[[7]]
C1_gene4_snps<-lst[[9]]
C1_gene6_snps<-lst[[10]]
C1_gene10_snps<-lst[[2]]
C1_gene12_snps<-lst[[3]]
C1_gene14_snps<-lst[[4]]
C1_gene16_snps<-lst[[5]]
C1_gene18_snps<-lst[[6]]
C1_gene20_snps<-lst[[8]]
#Discard background mutations
C1_gene0_snps_wt_bg<-C1_gene0_snps[ !(C1_gene0_snps$Position %in% All_C_snps$Position), ]
C1_gene2_snps_wt_bg<-C1_gene2_snps[ !(C1_gene2_snps$Position %in% All_C_snps$Position), ]
C1_gene4_snps_wt_bg<-C1_gene4_snps[ !(C1_gene4_snps$Position %in% All_C_snps$Position), ]
C1_gene6_snps_wt_bg<-C1_gene6_snps[ !(C1_gene6_snps$Position %in% All_C_snps$Position), ]
C1_gene10_snps_wt_bg<-C1_gene10_snps[ !(C1_gene10_snps$Position %in% All_C_snps$Position), ]
C1_gene12_snps_wt_bg<-C1_gene12_snps[ !(C1_gene12_snps$Position %in% All_C_snps$Position), ]
C1_gene14_snps_wt_bg<-C1_gene14_snps[ !(C1_gene14_snps$Position %in% All_C_snps$Position), ]
C1_gene16_snps_wt_bg<-C1_gene16_snps[ !(C1_gene16_snps$Position %in% All_C_snps$Position), ]
C1_gene18_snps_wt_bg<-C1_gene18_snps[ !(C1_gene18_snps$Position %in% All_C_snps$Position), ]
C1_gene20_snps_wt_bg<-C1_gene20_snps[ !(C1_gene20_snps$Position %in% All_C_snps$Position), ]
#C1_mutations_from_0_to_20
df_merge0 <- merge(C1_gene0_snps_wt_bg,C1_gene2_snps_wt_bg,by="Position")
df_merge0 <- merge(df_merge0,C1_gene4_snps_wt_bg,by="Position")
df_merge0 <- merge(df_merge0,C1_gene6_snps_wt_bg,by="Position")
df_merge0 <- merge(df_merge0,C1_gene10_snps_wt_bg,by="Position")
df_merge0 <- merge(df_merge0,C1_gene12_snps_wt_bg,by="Position")
df_merge0 <- merge(df_merge0,C1_gene14_snps_wt_bg,by="Position")
df_merge0 <- merge(df_merge0,C1_gene16_snps_wt_bg,by="Position")
df_merge0 <- merge(df_merge0,C1_gene18_snps_wt_bg,by="Position")
C1_snps_from_0_to_20 <- merge(df_merge0,C1_gene20_snps_wt_bg,by="Position")
write.xlsx(C1_snps_from_0_to_20,"C1_snps_from_0_to_20.xlsx")
#Discard F0 mutations
C1_gene2_snps_wt_gene0<-C1_gene2_snps_wt_bg[ !(C1_gene2_snps_wt_bg$Position %in% C1_gene0_snps_wt_bg$Position), ]
C1_gene4_snps_wt_gene0<-C1_gene4_snps_wt_bg[ !(C1_gene4_snps_wt_bg$Position %in% C1_gene0_snps_wt_bg$Position), ]
C1_gene6_snps_wt_gene0<-C1_gene6_snps_wt_bg[ !(C1_gene6_snps_wt_bg$Position %in% C1_gene0_snps_wt_bg$Position), ]
C1_gene10_snps_wt_gene0<-C1_gene10_snps_wt_bg[ !(C1_gene10_snps_wt_bg$Position %in% C1_gene0_snps_wt_bg$Position), ]
C1_gene12_snps_wt_gene0<-C1_gene12_snps_wt_bg[ !(C1_gene12_snps_wt_bg$Position %in% C1_gene0_snps_wt_bg$Position), ]
C1_gene14_snps_wt_gene0<-C1_gene14_snps_wt_bg[ !(C1_gene14_snps_wt_bg$Position %in% C1_gene0_snps_wt_bg$Position), ]
C1_gene16_snps_wt_gene0<-C1_gene16_snps_wt_bg[ !(C1_gene16_snps_wt_bg$Position %in% C1_gene0_snps_wt_bg$Position), ]
C1_gene18_snps_wt_gene0<-C1_gene18_snps_wt_bg[ !(C1_gene18_snps_wt_bg$Position %in% C1_gene0_snps_wt_bg$Position), ]
C1_gene20_snps_wt_gene0<-C1_gene20_snps_wt_bg[ !(C1_gene20_snps_wt_bg$Position %in% C1_gene0_snps_wt_bg$Position), ]
#C1_mutations_from_2_to_20
df_merge <- merge(C1_gene2_snps_wt_gene0,C1_gene4_snps_wt_gene0,by="Position")
df_merge <- merge(df_merge,C1_gene6_snps_wt_gene0,by="Position")
df_merge <- merge(df_merge,C1_gene10_snps_wt_gene0,by="Position")
df_merge <- merge(df_merge,C1_gene12_snps_wt_gene0,by="Position")
df_merge <- merge(df_merge,C1_gene14_snps_wt_gene0,by="Position")
df_merge <- merge(df_merge,C1_gene16_snps_wt_gene0,by="Position")
df_merge <- merge(df_merge,C1_gene18_snps_wt_gene0,by="Position")
C1_snps_from_2_to_20 <- merge(df_merge,C1_gene20_snps_wt_gene0,by="Position")
write.xlsx(C1_snps_from_2_to_20,"C1_snps_from_2_to_20.xlsx")
#Discard F2 snps
C1_gene4_snps_wt_gene2<-C1_gene4_snps_wt_gene0[ !(C1_gene4_snps_wt_gene0$Position %in% C1_gene2_snps_wt_gene0$Position), ]
C1_gene6_snps_wt_gene2<-C1_gene6_snps_wt_gene0[ !(C1_gene6_snps_wt_gene0$Position %in% C1_gene2_snps_wt_gene0$Position), ]
C1_gene10_snps_wt_gene2<-C1_gene10_snps_wt_gene0[ !(C1_gene10_snps_wt_gene0$Position %in% C1_gene2_snps_wt_gene0$Position), ]
C1_gene12_snps_wt_gene2<-C1_gene12_snps_wt_gene0[ !(C1_gene12_snps_wt_gene0$Position %in% C1_gene2_snps_wt_gene0$Position), ]
C1_gene14_snps_wt_gene2<-C1_gene14_snps_wt_gene0[ !(C1_gene14_snps_wt_gene0$Position %in% C1_gene2_snps_wt_gene0$Position), ]
C1_gene16_snps_wt_gene2<-C1_gene16_snps_wt_gene0[ !(C1_gene16_snps_wt_gene0$Position %in% C1_gene2_snps_wt_gene0$Position), ]
C1_gene18_snps_wt_gene2<-C1_gene18_snps_wt_gene0[ !(C1_gene18_snps_wt_gene0$Position %in% C1_gene2_snps_wt_gene0$Position), ]
C1_gene20_snps_wt_gene2<-C1_gene20_snps_wt_gene0[ !(C1_gene20_snps_wt_gene0$Position %in% C1_gene2_snps_wt_gene0$Position), ]
#C1_snps_from_4_to_20
df_merge4 <- merge(C1_gene4_snps_wt_gene2,C1_gene6_snps_wt_gene2,by="Position")
df_merge4 <- merge(df_merge4,C1_gene10_snps_wt_gene2,by="Position")
df_merge4 <- merge(df_merge4,C1_gene12_snps_wt_gene2,by="Position")
df_merge4 <- merge(df_merge4,C1_gene14_snps_wt_gene2,by="Position")
df_merge4 <- merge(df_merge4,C1_gene16_snps_wt_gene2,by="Position")
df_merge4 <- merge(df_merge4,C1_gene18_snps_wt_gene2,by="Position")
C1_snps_from_4_to_20 <- merge(df_merge4,C1_gene20_snps_wt_gene2,by="Position")
write.xlsx(C1_snps_from_4_to_20,"C1_snps_from_4_to_20.xlsx")
#Discard F4 snps
C1_gene6_snps_wt_gene4<-C1_gene6_snps_wt_gene2[ !(C1_gene6_snps_wt_gene2$Position %in% C1_gene4_snps_wt_gene2$Position), ]
C1_gene10_snps_wt_gene4<-C1_gene10_snps_wt_gene2[ !(C1_gene10_snps_wt_gene2$Position %in% C1_gene4_snps_wt_gene2$Position), ]
C1_gene12_snps_wt_gene4<-C1_gene12_snps_wt_gene2[ !(C1_gene12_snps_wt_gene2$Position %in% C1_gene4_snps_wt_gene2$Position), ]
C1_gene14_snps_wt_gene4<-C1_gene14_snps_wt_gene2[ !(C1_gene14_snps_wt_gene2$Position %in% C1_gene4_snps_wt_gene2$Position), ]
C1_gene16_snps_wt_gene4<-C1_gene16_snps_wt_gene2[ !(C1_gene16_snps_wt_gene2$Position %in% C1_gene4_snps_wt_gene2$Position), ]
C1_gene18_snps_wt_gene4<-C1_gene18_snps_wt_gene2[ !(C1_gene18_snps_wt_gene2$Position %in% C1_gene4_snps_wt_gene2$Position), ]
C1_gene20_snps_wt_gene4<-C1_gene20_snps_wt_gene2[ !(C1_gene20_snps_wt_gene2$Position %in% C1_gene4_snps_wt_gene2$Position), ]
#C1_snps_from_6_to_20
df_merge6 <- merge(C1_gene6_snps_wt_gene4,C1_gene10_snps_wt_gene4,by="Position")
df_merge6 <- merge(df_merge6,C1_gene12_snps_wt_gene4,by="Position")
df_merge6 <- merge(df_merge6,C1_gene14_snps_wt_gene4,by="Position")
df_merge6 <- merge(df_merge6,C1_gene16_snps_wt_gene4,by="Position")
df_merge6 <- merge(df_merge6,C1_gene18_snps_wt_gene4,by="Position")
C1_snps_from_6_to_20 <- merge(df_merge6,C1_gene20_snps_wt_gene4,by="Position")
write.xlsx(C1_snps_from_6_to_20,"C1_snps_from_6_to_20.xlsx")
#Discard F6 snps
C1_gene10_snps_wt_gene6<-C1_gene10_snps_wt_gene4[ !(C1_gene10_snps_wt_gene4$Position %in% C1_gene6_snps_wt_gene4$Position), ]
C1_gene12_snps_wt_gene6<-C1_gene12_snps_wt_gene4[ !(C1_gene12_snps_wt_gene4$Position %in% C1_gene6_snps_wt_gene4$Position), ]
C1_gene14_snps_wt_gene6<-C1_gene14_snps_wt_gene4[ !(C1_gene14_snps_wt_gene4$Position %in% C1_gene6_snps_wt_gene4$Position), ]
C1_gene16_snps_wt_gene6<-C1_gene16_snps_wt_gene4[ !(C1_gene16_snps_wt_gene4$Position %in% C1_gene6_snps_wt_gene4$Position), ]
C1_gene18_snps_wt_gene6<-C1_gene18_snps_wt_gene4[ !(C1_gene18_snps_wt_gene4$Position %in% C1_gene6_snps_wt_gene4$Position), ]
C1_gene20_snps_wt_gene6<-C1_gene20_snps_wt_gene4[ !(C1_gene20_snps_wt_gene4$Position %in% C1_gene6_snps_wt_gene4$Position), ]
#C1_snps_from_10_to_20
df_merge10 <- merge(C1_gene10_snps_wt_gene6,C1_gene12_snps_wt_gene6,by="Position")
df_merge10 <- merge(df_merge10,C1_gene14_snps_wt_gene6,by="Position")
df_merge10 <- merge(df_merge10,C1_gene16_snps_wt_gene6,by="Position")
df_merge10 <- merge(df_merge10,C1_gene18_snps_wt_gene6,by="Position")
C1_snps_from_10_to_20 <- merge(df_merge10,C1_gene20_snps_wt_gene6,by="Position")
write.xlsx(C1_snps_from_10_to_20,"C1_snps_from_10_to_20.xlsx")
#Discard F10 snps
C1_gene12_snps_wt_gene10<-C1_gene12_snps_wt_gene6[ !(C1_gene12_snps_wt_gene6$Position %in% C1_gene10_snps_wt_gene6$Position), ]
C1_gene14_snps_wt_gene10<-C1_gene14_snps_wt_gene6[ !(C1_gene14_snps_wt_gene6$Position %in% C1_gene10_snps_wt_gene6$Position), ]
C1_gene16_snps_wt_gene10<-C1_gene16_snps_wt_gene6[ !(C1_gene16_snps_wt_gene6$Position %in% C1_gene10_snps_wt_gene6$Position), ]
C1_gene18_snps_wt_gene10<-C1_gene18_snps_wt_gene6[ !(C1_gene18_snps_wt_gene6$Position %in% C1_gene10_snps_wt_gene6$Position), ]
C1_gene20_snps_wt_gene10<-C1_gene20_snps_wt_gene6[ !(C1_gene20_snps_wt_gene6$Position %in% C1_gene10_snps_wt_gene6$Position), ]
#C1_snps_from_12_to_20
df_merge12 <- merge(C1_gene12_snps_wt_gene10,C1_gene14_snps_wt_gene10,by="Position")
df_merge12 <- merge(df_merge12,C1_gene16_snps_wt_gene10,by="Position")
df_merge12 <- merge(df_merge12,C1_gene18_snps_wt_gene10,by="Position")
C1_snps_from_12_to_20 <- merge(df_merge12,C1_gene20_snps_wt_gene10,by="Position")
write.xlsx(C1_snps_from_12_to_20,"C1_snps_from_12_to_20.xlsx")
#Discard F12 snps
C1_gene14_snps_wt_gene12<-C1_gene14_snps_wt_gene10[ !(C1_gene14_snps_wt_gene10$Position %in% C1_gene12_snps_wt_gene10$Position), ]
C1_gene16_snps_wt_gene12<-C1_gene16_snps_wt_gene10[ !(C1_gene16_snps_wt_gene10$Position %in% C1_gene12_snps_wt_gene10$Position), ]
C1_gene18_snps_wt_gene12<-C1_gene18_snps_wt_gene10[ !(C1_gene18_snps_wt_gene10$Position %in% C1_gene12_snps_wt_gene10$Position), ]
C1_gene20_snps_wt_gene12<-C1_gene20_snps_wt_gene10[ !(C1_gene20_snps_wt_gene10$Position %in% C1_gene12_snps_wt_gene10$Position), ]
#C1_snps_from_14_to_20
df_merge14 <- merge(C1_gene14_snps_wt_gene12,C1_gene16_snps_wt_gene12,by="Position")
df_merge14 <- merge(df_merge14,C1_gene18_snps_wt_gene12,by="Position")
C1_snps_from_14_to_20 <- merge(df_merge14,C1_gene20_snps_wt_gene12,by="Position")
write.xlsx(C1_snps_from_14_to_20,"C1_snps_from_14_to_20.xlsx")
#Discard F14 snps
C1_gene16_snps_wt_gene14<-C1_gene16_snps_wt_gene12[ !(C1_gene16_snps_wt_gene12$Position %in% C1_gene14_snps_wt_gene12$Position), ]
C1_gene18_snps_wt_gene14<-C1_gene18_snps_wt_gene12[ !(C1_gene18_snps_wt_gene12$Position %in% C1_gene14_snps_wt_gene12$Position), ]
C1_gene20_snps_wt_gene14<-C1_gene20_snps_wt_gene12[ !(C1_gene20_snps_wt_gene12$Position %in% C1_gene14_snps_wt_gene12$Position), ]
#C1_snps_from_16_to_20
df_merge16 <- merge(C1_gene16_snps_wt_gene14,C1_gene18_snps_wt_gene14,by="Position")
C1_snps_from_16_to_20 <- merge(df_merge16,C1_gene20_snps_wt_gene14,by="Position")
write.xlsx(C1_snps_from_16_to_20,"C1_snps_from_16_to_20.xlsx")
#Discard F16 snps
C1_gene18_snps_wt_gene16<-C1_gene18_snps_wt_gene14[ !(C1_gene18_snps_wt_gene14$Position %in% C1_gene16_snps_wt_gene14$Position), ]
C1_gene20_snps_wt_gene16<-C1_gene20_snps_wt_gene14[ !(C1_gene20_snps_wt_gene14$Position %in% C1_gene16_snps_wt_gene14$Position), ]
#C1_snps_from_18_to_20
C1_snps_from_18_to_20 <- merge(C1_gene18_snps_wt_gene16,C1_gene20_snps_wt_gene16,by="Position")
write.xlsx(C1_snps_from_18_to_20,"C1_snps_from_18_to_20.xlsx")
#Discard F18 snps
C1_gene20_snps_wt_gene18<-C1_gene20_snps_wt_gene16[ !(C1_gene20_snps_wt_gene16$Position %in% C1_gene18_snps_wt_gene16$Position), ]
#C1_snps_from_18_to_20
C1_snps_from_20_to_20 <- C1_gene20_snps_wt_gene18
write.xlsx(C1_snps_from_20_to_20,"C1_snps_from_20_to_20.xlsx")
#####DNA mutations C2
setwd("~/Papers/Cisplatin paper/Varscan outputs 0.5/C2")
#Indels#
C2_intels_list <- list.files(path = "~/Papers/Cisplatin paper/Varscan outputs 0.5/C2", pattern = ".indel.xlsx")
lst <- lapply(C2_intels_list, read.xlsx)
C2_gene0_mut<-lst[[1]]
C2_gene2_mut<-lst[[5]]
C2_gene4_mut<-lst[[7]]
C2_gene6_mut<-lst[[8]]
C2_gene12_mut<-lst[[2]]
C2_gene16_mut<-lst[[3]]
C2_gene18_mut<-lst[[4]]
C2_gene20_mut<-lst[[6]]
#Discard background mutations
C2_gene0_mut_wt_bg<-C2_gene0_mut[ !(C2_gene0_mut$Position %in% All_C_indels$Position), ]
C2_gene2_mut_wt_bg<-C2_gene2_mut[ !(C2_gene2_mut$Position %in% All_C_indels$Position), ]
C2_gene4_mut_wt_bg<-C2_gene4_mut[ !(C2_gene4_mut$Position %in% All_C_indels$Position), ]
C2_gene6_mut_wt_bg<-C2_gene6_mut[ !(C2_gene6_mut$Position %in% All_C_indels$Position), ]
C2_gene12_mut_wt_bg<-C2_gene12_mut[ !(C2_gene12_mut$Position %in% All_C_indels$Position), ]
C2_gene16_mut_wt_bg<-C2_gene16_mut[ !(C2_gene16_mut$Position %in% All_C_indels$Position), ]
C2_gene18_mut_wt_bg<-C2_gene18_mut[ !(C2_gene18_mut$Position %in% All_C_indels$Position), ]
C2_gene20_mut_wt_bg<-C2_gene20_mut[ !(C2_gene20_mut$Position %in% All_C_indels$Position), ]
#C2_mutations_from_0_to_20
df_merge0 <- merge(C2_gene0_mut_wt_bg,C2_gene2_mut_wt_bg,by="Position")
df_merge0 <- merge(df_merge0,C2_gene4_mut_wt_bg,by="Position")
df_merge0 <- merge(df_merge0,C2_gene6_mut_wt_bg,by="Position")
df_merge0 <- merge(df_merge0,C2_gene12_mut_wt_bg,by="Position")
df_merge0 <- merge(df_merge0,C2_gene16_mut_wt_bg,by="Position")
df_merge0 <- merge(df_merge0,C2_gene18_mut_wt_bg,by="Position")
C2_indels_from_0_to_20 <- merge(df_merge0,C2_gene20_mut_wt_bg,by="Position")
write.xlsx(C2_indels_from_0_to_20,"C2_indels_from_0_to_20.xlsx")
#Discard F0 mutations
C2_gene2_mut_wt_gene0<-C2_gene2_mut_wt_bg[ !(C2_gene2_mut_wt_bg$Position %in% C2_gene0_mut_wt_bg$Position), ]
C2_gene4_mut_wt_gene0<-C2_gene4_mut_wt_bg[ !(C2_gene4_mut_wt_bg$Position %in% C2_gene0_mut_wt_bg$Position), ]
C2_gene6_mut_wt_gene0<-C2_gene6_mut_wt_bg[ !(C2_gene6_mut_wt_bg$Position %in% C2_gene0_mut_wt_bg$Position), ]
C2_gene12_mut_wt_gene0<-C2_gene12_mut_wt_bg[ !(C2_gene12_mut_wt_bg$Position %in% C2_gene0_mut_wt_bg$Position), ]
C2_gene16_mut_wt_gene0<-C2_gene16_mut_wt_bg[ !(C2_gene16_mut_wt_bg$Position %in% C2_gene0_mut_wt_bg$Position), ]
C2_gene18_mut_wt_gene0<-C2_gene18_mut_wt_bg[ !(C2_gene18_mut_wt_bg$Position %in% C2_gene0_mut_wt_bg$Position), ]
C2_gene20_mut_wt_gene0<-C2_gene20_mut_wt_bg[ !(C2_gene20_mut_wt_bg$Position %in% C2_gene0_mut_wt_bg$Position), ]
#C2_mutations_from_2_to_20
df_merge <- merge(C2_gene2_mut_wt_gene0,C2_gene4_mut_wt_gene0,by="Position")
df_merge <- merge(df_merge,C2_gene6_mut_wt_gene0,by="Position")
df_merge <- merge(df_merge,C2_gene12_mut_wt_gene0,by="Position")
df_merge <- merge(df_merge,C2_gene16_mut_wt_gene0,by="Position")
df_merge <- merge(df_merge,C2_gene18_mut_wt_gene0,by="Position")
C2_indels_from_2_to_20 <- merge(df_merge,C2_gene20_mut_wt_gene0,by="Position")
write.xlsx(C2_indels_from_2_to_20,"C2_indels_from_2_to_20.xlsx")
#Discard F2 mutations
C2_gene4_mut_wt_gene2<-C2_gene4_mut_wt_gene0[ !(C2_gene4_mut_wt_gene0$Position %in% C2_gene2_mut_wt_gene0$Position), ]
C2_gene6_mut_wt_gene2<-C2_gene6_mut_wt_gene0[ !(C2_gene6_mut_wt_gene0$Position %in% C2_gene2_mut_wt_gene0$Position), ]
C2_gene12_mut_wt_gene2<-C2_gene12_mut_wt_gene0[ !(C2_gene12_mut_wt_gene0$Position %in% C2_gene2_mut_wt_gene0$Position), ]
C2_gene16_mut_wt_gene2<-C2_gene16_mut_wt_gene0[ !(C2_gene16_mut_wt_gene0$Position %in% C2_gene2_mut_wt_gene0$Position), ]
C2_gene18_mut_wt_gene2<-C2_gene18_mut_wt_gene0[ !(C2_gene18_mut_wt_gene0$Position %in% C2_gene2_mut_wt_gene0$Position), ]
C2_gene20_mut_wt_gene2<-C2_gene20_mut_wt_gene0[ !(C2_gene20_mut_wt_gene0$Position %in% C2_gene2_mut_wt_gene0$Position), ]
#C2_mutations_from_4_to_20
df_merge4 <- merge(C2_gene4_mut_wt_gene2,C2_gene6_mut_wt_gene2,by="Position")
df_merge4 <- merge(df_merge4,C2_gene12_mut_wt_gene2,by="Position")
df_merge4 <- merge(df_merge4,C2_gene16_mut_wt_gene2,by="Position")
df_merge4 <- merge(df_merge4,C2_gene18_mut_wt_gene2,by="Position")
C2_indels_from_4_to_20 <- merge(df_merge4,C2_gene20_mut_wt_gene2,by="Position")
write.xlsx(C2_indels_from_4_to_20,"C2_indels_from_4_to_20.xlsx")
#Discard F4 mutations
C2_gene6_mut_wt_gene4<-C2_gene6_mut_wt_gene2[ !(C2_gene6_mut_wt_gene2$Position %in% C2_gene4_mut_wt_gene2$Position), ]
C2_gene12_mut_wt_gene4<-C2_gene12_mut_wt_gene2[ !(C2_gene12_mut_wt_gene2$Position %in% C2_gene4_mut_wt_gene2$Position), ]
C2_gene16_mut_wt_gene4<-C2_gene16_mut_wt_gene2[ !(C2_gene16_mut_wt_gene2$Position %in% C2_gene4_mut_wt_gene2$Position), ]
C2_gene18_mut_wt_gene4<-C2_gene18_mut_wt_gene2[ !(C2_gene18_mut_wt_gene2$Position %in% C2_gene4_mut_wt_gene2$Position), ]
C2_gene20_mut_wt_gene4<-C2_gene20_mut_wt_gene2[ !(C2_gene20_mut_wt_gene2$Position %in% C2_gene4_mut_wt_gene2$Position), ]
#C2_mutations_from_6_to_20
df_merge6 <- merge(C2_gene6_mut_wt_gene4,C2_gene12_mut_wt_gene4,by="Position")
df_merge6 <- merge(df_merge6,C2_gene16_mut_wt_gene4,by="Position")
df_merge6 <- merge(df_merge6,C2_gene18_mut_wt_gene4,by="Position")
C2_indels_from_6_to_20 <- merge(df_merge6,C2_gene20_mut_wt_gene4,by="Position")
write.xlsx(C2_indels_from_6_to_20,"C2_indels_from_6_to_20.xlsx")
#Discard F6 mutations
C2_gene12_mut_wt_gene6<-C2_gene12_mut_wt_gene4[ !(C2_gene12_mut_wt_gene4$Position %in% C2_gene6_mut_wt_gene4$Position), ]
C2_gene16_mut_wt_gene6<-C2_gene16_mut_wt_gene4[ !(C2_gene16_mut_wt_gene4$Position %in% C2_gene6_mut_wt_gene4$Position), ]
C2_gene18_mut_wt_gene6<-C2_gene18_mut_wt_gene4[ !(C2_gene18_mut_wt_gene4$Position %in% C2_gene6_mut_wt_gene4$Position), ]
C2_gene20_mut_wt_gene6<-C2_gene20_mut_wt_gene4[ !(C2_gene20_mut_wt_gene4$Position %in% C2_gene6_mut_wt_gene4$Position), ]
#C2_mutations_from_12_to_20
df_merge12 <- merge(C2_gene12_mut_wt_gene6,C2_gene16_mut_wt_gene6,by="Position")
df_merge12 <- merge(df_merge12,C2_gene18_mut_wt_gene6,by="Position")
C2_indels_from_12_to_20 <- merge(df_merge12,C2_gene20_mut_wt_gene6,by="Position")
write.xlsx(C2_indels_from_12_to_20,"C2_indels_from_12_to_20.xlsx")
#Discard F12 mutations
C2_gene16_mut_wt_gene12<-C2_gene16_mut_wt_gene6[ !(C2_gene16_mut_wt_gene6$Position %in% C2_gene12_mut_wt_gene6$Position), ]
C2_gene18_mut_wt_gene12<-C2_gene18_mut_wt_gene6[ !(C2_gene18_mut_wt_gene6$Position %in% C2_gene12_mut_wt_gene6$Position), ]
C2_gene20_mut_wt_gene12<-C2_gene20_mut_wt_gene6[ !(C2_gene20_mut_wt_gene6$Position %in% C2_gene12_mut_wt_gene6$Position), ]
#C2_mutations_from_16_to_20
df_merge16 <- merge(C2_gene16_mut_wt_gene12,C2_gene18_mut_wt_gene12,by="Position")
C2_indels_from_16_to_20 <- merge(df_merge16,C2_gene20_mut_wt_gene12,by="Position")
write.xlsx(C2_indels_from_16_to_20,"C2_indels_from_16_to_20.xlsx")
#Discard F16 mutations
C2_gene18_mut_wt_gene16<-C2_gene18_mut_wt_gene12[ !(C2_gene18_mut_wt_gene12$Position %in% C2_gene16_mut_wt_gene12$Position), ]
C2_gene20_mut_wt_gene16<-C2_gene20_mut_wt_gene12[ !(C2_gene20_mut_wt_gene12$Position %in% C2_gene16_mut_wt_gene12$Position), ]
#C2_mutations_from_18_to_20
C2_indels_from_18_to_20 <- merge(C2_gene18_mut_wt_gene16,C2_gene20_mut_wt_gene16,by="Position")
write.xlsx(C2_indels_from_18_to_20,"C2_indels_from_18_to_20.xlsx")
#Discard F18 mutations
C2_gene20_mut_wt_gene18<-C2_gene20_mut_wt_gene16[ !(C2_gene20_mut_wt_gene16$Position %in% C2_gene18_mut_wt_gene16$Position), ]
#C2_mutations_from_20_to_20
C2_indels_from_20_to_20 <- C2_gene20_mut_wt_gene18
write.xlsx(C2_indels_from_20_to_20,"C2_indels_from_20_to_20.xlsx")
#####Snps
C2_snps_list <- list.files(path = "~/Papers/Cisplatin paper/Varscan outputs 0.5/C2", pattern = ".snp.xlsx")
lst <- lapply(C2_snps_list, read.xlsx)
C2_gene0_snps<-lst[[1]]
C2_gene2_snps<-lst[[5]]
C2_gene4_snps<-lst[[7]]
C2_gene6_snps<-lst[[8]]
C2_gene12_snps<-lst[[2]]
C2_gene16_snps<-lst[[3]]
C2_gene18_snps<-lst[[4]]
C2_gene20_snps<-lst[[6]]
#Discard background mutations
C2_gene0_snps_wt_bg<-C2_gene0_snps[ !(C2_gene0_snps$Position %in% All_C_snps$Position), ]
C2_gene2_snps_wt_bg<-C2_gene2_snps[ !(C2_gene2_snps$Position %in% All_C_snps$Position), ]
C2_gene4_snps_wt_bg<-C2_gene4_snps[ !(C2_gene4_snps$Position %in% All_C_snps$Position), ]
C2_gene6_snps_wt_bg<-C2_gene6_snps[ !(C2_gene6_snps$Position %in% All_C_snps$Position), ]
C2_gene12_snps_wt_bg<-C2_gene12_snps[ !(C2_gene12_snps$Position %in% All_C_snps$Position), ]
C2_gene16_snps_wt_bg<-C2_gene16_snps[ !(C2_gene16_snps$Position %in% All_C_snps$Position), ]
C2_gene18_snps_wt_bg<-C2_gene18_snps[ !(C2_gene18_snps$Position %in% All_C_snps$Position), ]
C2_gene20_snps_wt_bg<-C2_gene20_snps[ !(C2_gene20_snps$Position %in% All_C_snps$Position), ]
#C2_mutations_from_0_to_20
df_merge0 <- merge(C2_gene0_snps_wt_bg,C2_gene2_snps_wt_bg,by="Position")
df_merge0 <- merge(df_merge0,C2_gene4_snps_wt_bg,by="Position")
df_merge0 <- merge(df_merge0,C2_gene6_snps_wt_bg,by="Position")
df_merge0 <- merge(df_merge0,C2_gene12_snps_wt_bg,by="Position")
df_merge0 <- merge(df_merge0,C2_gene16_snps_wt_bg,by="Position")
df_merge0 <- merge(df_merge0,C2_gene18_snps_wt_bg,by="Position")
C2_snps_from_0_to_20 <- merge(df_merge0,C2_gene20_snps_wt_bg,by="Position")
write.xlsx(C2_snps_from_0_to_20,"C2_snps_from_0_to_20.xlsx")
#Discard F0 mutations
C2_gene2_snps_wt_gene0<-C2_gene2_snps_wt_bg[ !(C2_gene2_snps_wt_bg$Position %in% C2_gene0_snps_wt_bg$Position), ]
C2_gene4_snps_wt_gene0<-C2_gene4_snps_wt_bg[ !(C2_gene4_snps_wt_bg$Position %in% C2_gene0_snps_wt_bg$Position), ]
C2_gene6_snps_wt_gene0<-C2_gene6_snps_wt_bg[ !(C2_gene6_snps_wt_bg$Position %in% C2_gene0_snps_wt_bg$Position), ]
C2_gene12_snps_wt_gene0<-C2_gene12_snps_wt_bg[ !(C2_gene12_snps_wt_bg$Position %in% C2_gene0_snps_wt_bg$Position), ]
C2_gene16_snps_wt_gene0<-C2_gene16_snps_wt_bg[ !(C2_gene16_snps_wt_bg$Position %in% C2_gene0_snps_wt_bg$Position), ]
C2_gene18_snps_wt_gene0<-C2_gene18_snps_wt_bg[ !(C2_gene18_snps_wt_bg$Position %in% C2_gene0_snps_wt_bg$Position), ]
C2_gene20_snps_wt_gene0<-C2_gene20_snps_wt_bg[ !(C2_gene20_snps_wt_bg$Position %in% C2_gene0_snps_wt_bg$Position), ]
#C2_mutations_from_2_to_20
df_merge <- merge(C2_gene2_snps_wt_gene0,C2_gene4_snps_wt_gene0,by="Position")
df_merge <- merge(df_merge,C2_gene6_snps_wt_gene0,by="Position")
df_merge <- merge(df_merge,C2_gene12_snps_wt_gene0,by="Position")
df_merge <- merge(df_merge,C2_gene16_snps_wt_gene0,by="Position")
df_merge <- merge(df_merge,C2_gene18_snps_wt_gene0,by="Position")
C2_snps_from_2_to_20 <- merge(df_merge,C2_gene20_snps_wt_gene0,by="Position")
write.xlsx(C2_snps_from_2_to_20,"C2_snps_from_2_to_20.xlsx")
#Discard F2 snps
C2_gene4_snps_wt_gene2<-C2_gene4_snps_wt_gene0[ !(C2_gene4_snps_wt_gene0$Position %in% C2_gene2_snps_wt_gene0$Position), ]
C2_gene6_snps_wt_gene2<-C2_gene6_snps_wt_gene0[ !(C2_gene6_snps_wt_gene0$Position %in% C2_gene2_snps_wt_gene0$Position), ]
C2_gene12_snps_wt_gene2<-C2_gene12_snps_wt_gene0[ !(C2_gene12_snps_wt_gene0$Position %in% C2_gene2_snps_wt_gene0$Position), ]
C2_gene16_snps_wt_gene2<-C2_gene16_snps_wt_gene0[ !(C2_gene16_snps_wt_gene0$Position %in% C2_gene2_snps_wt_gene0$Position), ]
C2_gene18_snps_wt_gene2<-C2_gene18_snps_wt_gene0[ !(C2_gene18_snps_wt_gene0$Position %in% C2_gene2_snps_wt_gene0$Position), ]
C2_gene20_snps_wt_gene2<-C2_gene20_snps_wt_gene0[ !(C2_gene20_snps_wt_gene0$Position %in% C2_gene2_snps_wt_gene0$Position), ]
#C2_snps_from_4_to_20
df_merge4 <- merge(C2_gene4_snps_wt_gene2,C2_gene6_snps_wt_gene2,by="Position")
df_merge4 <- merge(df_merge4,C2_gene12_snps_wt_gene2,by="Position")
df_merge4 <- merge(df_merge4,C2_gene16_snps_wt_gene2,by="Position")
df_merge4 <- merge(df_merge4,C2_gene18_snps_wt_gene2,by="Position")
C2_snps_from_4_to_20 <- merge(df_merge4,C2_gene20_snps_wt_gene2,by="Position")
write.xlsx(C2_snps_from_4_to_20,"C2_snps_from_4_to_20.xlsx")
#Discard F4 snps
C2_gene6_snps_wt_gene4<-C2_gene6_snps_wt_gene2[ !(C2_gene6_snps_wt_gene2$Position %in% C2_gene4_snps_wt_gene2$Position), ]
C2_gene12_snps_wt_gene4<-C2_gene12_snps_wt_gene2[ !(C2_gene12_snps_wt_gene2$Position %in% C2_gene4_snps_wt_gene2$Position), ]
C2_gene16_snps_wt_gene4<-C2_gene16_snps_wt_gene2[ !(C2_gene16_snps_wt_gene2$Position %in% C2_gene4_snps_wt_gene2$Position), ]
C2_gene18_snps_wt_gene4<-C2_gene18_snps_wt_gene2[ !(C2_gene18_snps_wt_gene2$Position %in% C2_gene4_snps_wt_gene2$Position), ]
C2_gene20_snps_wt_gene4<-C2_gene20_snps_wt_gene2[ !(C2_gene20_snps_wt_gene2$Position %in% C2_gene4_snps_wt_gene2$Position), ]
#C2_snps_from_6_to_20
df_merge6 <- merge(C2_gene6_snps_wt_gene4,C2_gene12_snps_wt_gene4,by="Position")
df_merge6 <- merge(df_merge6,C2_gene16_snps_wt_gene4,by="Position")
df_merge6 <- merge(df_merge6,C2_gene18_snps_wt_gene4,by="Position")
C2_snps_from_6_to_20 <- merge(df_merge6,C2_gene20_snps_wt_gene4,by="Position")
write.xlsx(C2_snps_from_6_to_20,"C2_snps_from_6_to_20.xlsx")
#Discard F6 snps
C2_gene12_snps_wt_gene6<-C2_gene12_snps_wt_gene4[ !(C2_gene12_snps_wt_gene4$Position %in% C2_gene6_snps_wt_gene4$Position), ]
C2_gene16_snps_wt_gene6<-C2_gene16_snps_wt_gene4[ !(C2_gene16_snps_wt_gene4$Position %in% C2_gene6_snps_wt_gene4$Position), ]
C2_gene18_snps_wt_gene6<-C2_gene18_snps_wt_gene4[ !(C2_gene18_snps_wt_gene4$Position %in% C2_gene6_snps_wt_gene4$Position), ]
C2_gene20_snps_wt_gene6<-C2_gene20_snps_wt_gene4[ !(C2_gene20_snps_wt_gene4$Position %in% C2_gene6_snps_wt_gene4$Position), ]
#C2_snps_from_12_to_20
df_merge12 <- merge(C2_gene12_snps_wt_gene6,C2_gene16_snps_wt_gene6,by="Position")
df_merge12 <- merge(df_merge12,C2_gene18_snps_wt_gene6,by="Position")
C2_snps_from_12_to_20 <- merge(df_merge12,C2_gene20_snps_wt_gene6,by="Position")
write.xlsx(C2_snps_from_12_to_20,"C2_snps_from_12_to_20.xlsx")
#Discard F12 snps
C2_gene16_snps_wt_gene12<-C2_gene16_snps_wt_gene6[ !(C2_gene16_snps_wt_gene6$Position %in% C2_gene12_snps_wt_gene6$Position), ]
C2_gene18_snps_wt_gene12<-C2_gene18_snps_wt_gene6[ !(C2_gene18_snps_wt_gene6$Position %in% C2_gene12_snps_wt_gene6$Position), ]
C2_gene20_snps_wt_gene12<-C2_gene20_snps_wt_gene6[ !(C2_gene20_snps_wt_gene6$Position %in% C2_gene12_snps_wt_gene6$Position), ]
#C2_snps_from_16_to_20
df_merge16 <- merge(C2_gene16_snps_wt_gene12,C2_gene18_snps_wt_gene12,by="Position")
C2_snps_from_16_to_20 <- merge(df_merge16,C2_gene20_snps_wt_gene12,by="Position")
write.xlsx(C2_snps_from_16_to_20,"C2_snps_from_16_to_20.xlsx")
#Discard F16 snps
C2_gene18_snps_wt_gene16<-C2_gene18_snps_wt_gene12[ !(C2_gene18_snps_wt_gene12$Position %in% C2_gene16_snps_wt_gene12$Position), ]
C2_gene20_snps_wt_gene16<-C2_gene20_snps_wt_gene12[ !(C2_gene20_snps_wt_gene12$Position %in% C2_gene16_snps_wt_gene12$Position), ]
#C2_snps_from_18_to_20
C2_snps_from_18_to_20 <- merge(C2_gene18_snps_wt_gene16,C2_gene20_snps_wt_gene16,by="Position")
write.xlsx(C2_snps_from_18_to_20,"C2_snps_from_18_to_20.xlsx")
#Discard F18 snps
C2_gene20_snps_wt_gene18<-C2_gene20_snps_wt_gene16[ !(C2_gene20_snps_wt_gene16$Position %in% C2_gene18_snps_wt_gene16$Position), ]
#C2_snps_from_20_to_20
C2_snps_from_20_to_20 <- C2_gene20_snps_wt_gene18
write.xlsx(C2_snps_from_20_to_20,"C2_snps_from_20_to_20.xlsx")
#####DNA mutations LD1
setwd("~/Papers/Cisplatin paper/Varscan outputs 0.5/LD1")
######Indels
LD1_intels_list <- list.files(path = "~/Papers/Cisplatin paper/Varscan outputs 0.5/LD1", pattern = ".indel.xlsx")
lst <- lapply(LD1_intels_list, read.xlsx)
LD1_gene0_mut<-lst[[1]]
LD1_gene2_mut<-lst[[5]]
LD1_gene4_mut<-lst[[6]]
LD1_gene8_mut<-lst[[7]]
LD1_gene10_mut<-lst[[2]]
LD1_gene12_mut<-lst[[3]]
LD1_gene16_mut<-lst[[4]]
#Discard background mutations
LD1_gene0_mut_wt_bg<-LD1_gene0_mut[ !(LD1_gene0_mut$Position %in% All_C_indels$Position), ]
LD1_gene2_mut_wt_bg<-LD1_gene2_mut[ !(LD1_gene2_mut$Position %in% All_C_indels$Position), ]
LD1_gene4_mut_wt_bg<-LD1_gene4_mut[ !(LD1_gene4_mut$Position %in% All_C_indels$Position), ]
LD1_gene8_mut_wt_bg<-LD1_gene8_mut[ !(LD1_gene8_mut$Position %in% All_C_indels$Position), ]
LD1_gene10_mut_wt_bg<-LD1_gene10_mut[ !(LD1_gene10_mut$Position %in% All_C_indels$Position), ]
LD1_gene12_mut_wt_bg<-LD1_gene12_mut[ !(LD1_gene12_mut$Position %in% All_C_indels$Position), ]
LD1_gene16_mut_wt_bg<-LD1_gene16_mut[ !(LD1_gene16_mut$Position %in% All_C_indels$Position), ]
#LD1_mutations_from_0_to_20
df_merge0 <- merge(LD1_gene0_mut_wt_bg,LD1_gene2_mut_wt_bg,by="Position")
df_merge0 <- merge(df_merge0,LD1_gene4_mut_wt_bg,by="Position")
df_merge0 <- merge(df_merge0,LD1_gene8_mut_wt_bg,by="Position")
df_merge0 <- merge(df_merge0,LD1_gene10_mut_wt_bg,by="Position")
df_merge0 <- merge(df_merge0,LD1_gene12_mut_wt_bg,by="Position")
LD1_indels_from_0_to_16 <- merge(df_merge0,LD1_gene16_mut_wt_bg,by="Position")
write.xlsx(LD1_indels_from_0_to_16,"LD1_indels_from_0_to_16.xlsx")
#Discard F0 mutations
LD1_gene2_mut_wt_gene0<-LD1_gene2_mut_wt_bg[ !(LD1_gene2_mut_wt_bg$Position %in% LD1_gene0_mut_wt_bg$Position), ]
LD1_gene4_mut_wt_gene0<-LD1_gene4_mut_wt_bg[ !(LD1_gene4_mut_wt_bg$Position %in% LD1_gene0_mut_wt_bg$Position), ]
LD1_gene8_mut_wt_gene0<-LD1_gene8_mut_wt_bg[ !(LD1_gene8_mut_wt_bg$Position %in% LD1_gene0_mut_wt_bg$Position), ]
LD1_gene10_mut_wt_gene0<-LD1_gene10_mut_wt_bg[ !(LD1_gene10_mut_wt_bg$Position %in% LD1_gene0_mut_wt_bg$Position), ]
LD1_gene12_mut_wt_gene0<-LD1_gene12_mut_wt_bg[ !(LD1_gene12_mut_wt_bg$Position %in% LD1_gene0_mut_wt_bg$Position), ]
LD1_gene16_mut_wt_gene0<-LD1_gene16_mut_wt_bg[ !(LD1_gene16_mut_wt_bg$Position %in% LD1_gene0_mut_wt_bg$Position), ]
#LD1_mutations_from_2_to_16
df_merge <- merge(LD1_gene2_mut_wt_gene0,LD1_gene4_mut_wt_gene0,by="Position")
df_merge <- merge(df_merge,LD1_gene8_mut_wt_gene0,by="Position")
df_merge <- merge(df_merge,LD1_gene10_mut_wt_gene0,by="Position")
df_merge <- merge(df_merge,LD1_gene12_mut_wt_gene0,by="Position")
LD1_indels_from_2_to_16 <- merge(df_merge,LD1_gene16_mut_wt_gene0,by="Position")
write.xlsx(LD1_indels_from_2_to_16,"LD1_indels_from_2_to_16.xlsx")
#Discard F2 mutations
LD1_gene4_mut_wt_gene2<-LD1_gene4_mut_wt_gene0[ !(LD1_gene4_mut_wt_gene0$Position %in% LD1_gene2_mut_wt_gene0$Position), ]
LD1_gene8_mut_wt_gene2<-LD1_gene8_mut_wt_gene0[ !(LD1_gene8_mut_wt_gene0$Position %in% LD1_gene2_mut_wt_gene0$Position), ]
LD1_gene10_mut_wt_gene2<-LD1_gene10_mut_wt_gene0[ !(LD1_gene10_mut_wt_gene0$Position %in% LD1_gene2_mut_wt_gene0$Position), ]
LD1_gene12_mut_wt_gene2<-LD1_gene12_mut_wt_gene0[ !(LD1_gene12_mut_wt_gene0$Position %in% LD1_gene2_mut_wt_gene0$Position), ]
LD1_gene16_mut_wt_gene2<-LD1_gene16_mut_wt_gene0[ !(LD1_gene16_mut_wt_gene0$Position %in% LD1_gene2_mut_wt_gene0$Position), ]
#LD1_mutations_from_4_to_16
df_merge4 <- merge(LD1_gene4_mut_wt_gene2,LD1_gene8_mut_wt_gene2,by="Position")
df_merge4 <- merge(df_merge4,LD1_gene10_mut_wt_gene2,by="Position")
df_merge4 <- merge(df_merge4,LD1_gene12_mut_wt_gene2,by="Position")
LD1_indels_from_4_to_16 <- merge(df_merge4,LD1_gene16_mut_wt_gene2,by="Position")
write.xlsx(LD1_indels_from_4_to_16,"LD1_indels_from_4_to_16.xlsx")
#Discard F4 mutations
LD1_gene8_mut_wt_gene4<-LD1_gene8_mut_wt_gene2[ !(LD1_gene8_mut_wt_gene2$Position %in% LD1_gene4_mut_wt_gene2$Position), ]
LD1_gene10_mut_wt_gene4<-LD1_gene10_mut_wt_gene2[ !(LD1_gene10_mut_wt_gene2$Position %in% LD1_gene4_mut_wt_gene2$Position), ]
LD1_gene12_mut_wt_gene4<-LD1_gene12_mut_wt_gene2[ !(LD1_gene12_mut_wt_gene2$Position %in% LD1_gene4_mut_wt_gene2$Position), ]
LD1_gene16_mut_wt_gene4<-LD1_gene16_mut_wt_gene2[ !(LD1_gene16_mut_wt_gene2$Position %in% LD1_gene4_mut_wt_gene2$Position), ]
#LD1_mutations_from_8_to_16
df_merge8 <- merge(LD1_gene8_mut_wt_gene4,LD1_gene10_mut_wt_gene4,by="Position")
df_merge8 <- merge(df_merge8,LD1_gene12_mut_wt_gene4,by="Position")
LD1_indels_from_8_to_16 <- merge(df_merge8,LD1_gene16_mut_wt_gene4,by="Position")
write.xlsx(LD1_indels_from_8_to_16,"LD1_indels_from_8_to_16.xlsx")
#Discard F8 mutations
LD1_gene10_mut_wt_gene8<-LD1_gene10_mut_wt_gene4[ !(LD1_gene10_mut_wt_gene4$Position %in% LD1_gene8_mut_wt_gene4$Position), ]
LD1_gene12_mut_wt_gene6<-LD1_gene12_mut_wt_gene4[ !(LD1_gene12_mut_wt_gene4$Position %in% LD1_gene8_mut_wt_gene4$Position), ]
LD1_gene16_mut_wt_gene6<-LD1_gene16_mut_wt_gene4[ !(LD1_gene16_mut_wt_gene4$Position %in% LD1_gene8_mut_wt_gene4$Position), ]
#LD1_mutations_from_10_to_16
df_merge10 <- merge(LD1_gene10_mut_wt_gene8,LD1_gene12_mut_wt_gene6,by="Position")
LD1_indels_from_10_to_16 <- merge(df_merge10,LD1_gene16_mut_wt_gene6,by="Position")
write.xlsx(LD1_indels_from_10_to_16,"LD1_indels_from_10_to_16.xlsx")
#Discard F10 mutations
LD1_gene12_mut_wt_gene10<-LD1_gene12_mut_wt_gene6[ !(LD1_gene12_mut_wt_gene6$Position %in% LD1_gene10_mut_wt_gene8$Position), ]
LD1_gene16_mut_wt_gene10<-LD1_gene16_mut_wt_gene6[ !(LD1_gene16_mut_wt_gene6$Position %in% LD1_gene10_mut_wt_gene8$Position), ]
#LD1_mutations_from_12_to_16
LD1_indels_from_12_to_16 <- merge(LD1_gene12_mut_wt_gene10,LD1_gene16_mut_wt_gene10,by="Position")
write.xlsx(LD1_indels_from_12_to_16,"LD1_indels_from_12_to_16.xlsx")
#Discard F12 mutations
LD1_gene16_mut_wt_gene12<-LD1_gene16_mut_wt_gene10[ !(LD1_gene16_mut_wt_gene10$Position %in% LD1_gene12_mut_wt_gene10$Position), ]
#LD1_mutations_from_16_to_16
LD1_indels_from_16_to_16 <- LD1_gene16_mut_wt_gene12
write.xlsx(LD1_indels_from_16_to_16,"LD1_indels_from_16_to_16.xlsx")
####Snps
LD1_snps_list <- list.files(path = "~/Papers/Cisplatin paper/Varscan outputs 0.5/LD1", pattern = ".snp.xlsx")
lst <- lapply(LD1_snps_list, read.xlsx)
LD1_gene0_mut<-lst[[1]]
LD1_gene2_mut<-lst[[5]]
LD1_gene4_mut<-lst[[6]]
LD1_gene8_mut<-lst[[7]]
LD1_gene10_mut<-lst[[2]]
LD1_gene12_mut<-lst[[3]]
LD1_gene16_mut<-lst[[4]]
#Discard background mutations
LD1_gene0_mut_wt_bg<-LD1_gene0_mut[ !(LD1_gene0_mut$Position %in% All_C_snps$Position), ]
LD1_gene2_mut_wt_bg<-LD1_gene2_mut[ !(LD1_gene2_mut$Position %in% All_C_snps$Position), ]
LD1_gene4_mut_wt_bg<-LD1_gene4_mut[ !(LD1_gene4_mut$Position %in% All_C_snps$Position), ]
LD1_gene8_mut_wt_bg<-LD1_gene8_mut[ !(LD1_gene8_mut$Position %in% All_C_snps$Position), ]
LD1_gene10_mut_wt_bg<-LD1_gene10_mut[ !(LD1_gene10_mut$Position %in% All_C_snps$Position), ]
LD1_gene12_mut_wt_bg<-LD1_gene12_mut[ !(LD1_gene12_mut$Position %in% All_C_snps$Position), ]
LD1_gene16_mut_wt_bg<-LD1_gene16_mut[ !(LD1_gene16_mut$Position %in% All_C_snps$Position), ]
#LD1_mutations_from_0_to_20
df_merge0 <- merge(LD1_gene0_mut_wt_bg,LD1_gene2_mut_wt_bg,by="Position")
df_merge0 <- merge(df_merge0,LD1_gene4_mut_wt_bg,by="Position")
df_merge0 <- merge(df_merge0,LD1_gene8_mut_wt_bg,by="Position")
df_merge0 <- merge(df_merge0,LD1_gene10_mut_wt_bg,by="Position")
df_merge0 <- merge(df_merge0,LD1_gene12_mut_wt_bg,by="Position")
LD1_snps_from_0_to_16 <- merge(df_merge0,LD1_gene16_mut_wt_bg,by="Position")
write.xlsx(LD1_snps_from_0_to_16,"LD1_snps_from_0_to_16.xlsx")
#Discard F0 mutations
LD1_gene2_mut_wt_gene0<-LD1_gene2_mut_wt_bg[ !(LD1_gene2_mut_wt_bg$Position %in% LD1_gene0_mut_wt_bg$Position), ]
LD1_gene4_mut_wt_gene0<-LD1_gene4_mut_wt_bg[ !(LD1_gene4_mut_wt_bg$Position %in% LD1_gene0_mut_wt_bg$Position), ]
LD1_gene8_mut_wt_gene0<-LD1_gene8_mut_wt_bg[ !(LD1_gene8_mut_wt_bg$Position %in% LD1_gene0_mut_wt_bg$Position), ]
LD1_gene10_mut_wt_gene0<-LD1_gene10_mut_wt_bg[ !(LD1_gene10_mut_wt_bg$Position %in% LD1_gene0_mut_wt_bg$Position), ]
LD1_gene12_mut_wt_gene0<-LD1_gene12_mut_wt_bg[ !(LD1_gene12_mut_wt_bg$Position %in% LD1_gene0_mut_wt_bg$Position), ]
LD1_gene16_mut_wt_gene0<-LD1_gene16_mut_wt_bg[ !(LD1_gene16_mut_wt_bg$Position %in% LD1_gene0_mut_wt_bg$Position), ]
#LD1_mutations_from_2_to_16
df_merge <- merge(LD1_gene2_mut_wt_gene0,LD1_gene4_mut_wt_gene0,by="Position")
df_merge <- merge(df_merge,LD1_gene8_mut_wt_gene0,by="Position")
df_merge <- merge(df_merge,LD1_gene10_mut_wt_gene0,by="Position")
df_merge <- merge(df_merge,LD1_gene12_mut_wt_gene0,by="Position")
LD1_snps_from_2_to_16 <- merge(df_merge,LD1_gene16_mut_wt_gene0,by="Position")
write.xlsx(LD1_snps_from_2_to_16,"LD1_snps_from_2_to_16.xlsx")
#Discard F2 mutations
LD1_gene4_mut_wt_gene2<-LD1_gene4_mut_wt_gene0[ !(LD1_gene4_mut_wt_gene0$Position %in% LD1_gene2_mut_wt_gene0$Position), ]
LD1_gene8_mut_wt_gene2<-LD1_gene8_mut_wt_gene0[ !(LD1_gene8_mut_wt_gene0$Position %in% LD1_gene2_mut_wt_gene0$Position), ]
LD1_gene10_mut_wt_gene2<-LD1_gene10_mut_wt_gene0[ !(LD1_gene10_mut_wt_gene0$Position %in% LD1_gene2_mut_wt_gene0$Position), ]
LD1_gene12_mut_wt_gene2<-LD1_gene12_mut_wt_gene0[ !(LD1_gene12_mut_wt_gene0$Position %in% LD1_gene2_mut_wt_gene0$Position), ]
LD1_gene16_mut_wt_gene2<-LD1_gene16_mut_wt_gene0[ !(LD1_gene16_mut_wt_gene0$Position %in% LD1_gene2_mut_wt_gene0$Position), ]
#LD1_mutations_from_4_to_16
df_merge4 <- merge(LD1_gene4_mut_wt_gene2,LD1_gene8_mut_wt_gene2,by="Position")
df_merge4 <- merge(df_merge4,LD1_gene10_mut_wt_gene2,by="Position")
df_merge4 <- merge(df_merge4,LD1_gene12_mut_wt_gene2,by="Position")
LD1_snps_from_4_to_16 <- merge(df_merge4,LD1_gene16_mut_wt_gene2,by="Position")
write.xlsx(LD1_snps_from_4_to_16,"LD1_snps_from_4_to_16.xlsx")
#Discard F4 mutations
LD1_gene8_mut_wt_gene4<-LD1_gene8_mut_wt_gene2[ !(LD1_gene8_mut_wt_gene2$Position %in% LD1_gene4_mut_wt_gene2$Position), ]
LD1_gene10_mut_wt_gene4<-LD1_gene10_mut_wt_gene2[ !(LD1_gene10_mut_wt_gene2$Position %in% LD1_gene4_mut_wt_gene2$Position), ]
LD1_gene12_mut_wt_gene4<-LD1_gene12_mut_wt_gene2[ !(LD1_gene12_mut_wt_gene2$Position %in% LD1_gene4_mut_wt_gene2$Position), ]
LD1_gene16_mut_wt_gene4<-LD1_gene16_mut_wt_gene2[ !(LD1_gene16_mut_wt_gene2$Position %in% LD1_gene4_mut_wt_gene2$Position), ]
#LD1_mutations_from_8_to_16
df_merge8 <- merge(LD1_gene8_mut_wt_gene4,LD1_gene10_mut_wt_gene4,by="Position")
df_merge8 <- merge(df_merge8,LD1_gene12_mut_wt_gene4,by="Position")
LD1_snps_from_8_to_16 <- merge(df_merge8,LD1_gene16_mut_wt_gene4,by="Position")
write.xlsx(LD1_snps_from_8_to_16,"LD1_snps_from_8_to_16.xlsx")
#Discard F8 mutations
LD1_gene10_mut_wt_gene8<-LD1_gene10_mut_wt_gene4[ !(LD1_gene10_mut_wt_gene4$Position %in% LD1_gene8_mut_wt_gene4$Position), ]
LD1_gene12_mut_wt_gene6<-LD1_gene12_mut_wt_gene4[ !(LD1_gene12_mut_wt_gene4$Position %in% LD1_gene8_mut_wt_gene4$Position), ]
LD1_gene16_mut_wt_gene6<-LD1_gene16_mut_wt_gene4[ !(LD1_gene16_mut_wt_gene4$Position %in% LD1_gene8_mut_wt_gene4$Position), ]
#LD1_mutations_from_10_to_16
df_merge10 <- merge(LD1_gene10_mut_wt_gene8,LD1_gene12_mut_wt_gene6,by="Position")
LD1_snps_from_10_to_16 <- merge(df_merge10,LD1_gene16_mut_wt_gene6,by="Position")
write.xlsx(LD1_snps_from_10_to_16,"LD1_snps_from_10_to_16.xlsx")
#Discard F10 mutations
LD1_gene12_mut_wt_gene10<-LD1_gene12_mut_wt_gene6[ !(LD1_gene12_mut_wt_gene6$Position %in% LD1_gene10_mut_wt_gene8$Position), ]
LD1_gene16_mut_wt_gene10<-LD1_gene16_mut_wt_gene6[ !(LD1_gene16_mut_wt_gene6$Position %in% LD1_gene10_mut_wt_gene8$Position), ]
#LD1_mutations_from_12_to_16
LD1_snps_from_12_to_16 <- merge(LD1_gene12_mut_wt_gene10,LD1_gene16_mut_wt_gene10,by="Position")
write.xlsx(LD1_snps_from_12_to_16,"LD1_snps_from_12_to_16.xlsx")
#Discard F12 mutations
LD1_gene16_mut_wt_gene12<-LD1_gene16_mut_wt_gene10[ !(LD1_gene16_mut_wt_gene10$Position %in% LD1_gene12_mut_wt_gene10$Position), ]
#LD1_mutations_from_16_to_16
LD1_snps_from_16_to_16 <- LD1_gene16_mut_wt_gene12
write.xlsx(LD1_snps_from_16_to_16,"LD1_snps_from_16_to_16.xlsx")
####DNA mutations LD2
setwd("~/Papers/Cisplatin paper/Varscan outputs 0.5/LD2")
####Indels
LD2_intels_list <- list.files(path = "~/Papers/Cisplatin paper/Varscan outputs 0.5/LD2", pattern = ".indel.xlsx")
lst <- lapply(LD2_intels_list, read.xlsx)
LD2_gene0_mut<-lst[[1]]
LD2_gene4_mut<-lst[[5]]
LD2_gene8_mut<-lst[[6]]
LD2_gene10_mut<-lst[[2]]
LD2_gene12_mut<-lst[[3]]
LD2_gene16_mut<-lst[[4]]
#Discard background mutations
LD2_gene0_mut_wt_bg<-LD2_gene0_mut[ !(LD2_gene0_mut$Position %in% All_C_indels$Position), ]
LD2_gene4_mut_wt_bg<-LD2_gene4_mut[ !(LD2_gene4_mut$Position %in% All_C_indels$Position), ]
LD2_gene8_mut_wt_bg<-LD2_gene8_mut[ !(LD2_gene8_mut$Position %in% All_C_indels$Position), ]
LD2_gene10_mut_wt_bg<-LD2_gene10_mut[ !(LD2_gene10_mut$Position %in% All_C_indels$Position), ]
LD2_gene12_mut_wt_bg<-LD2_gene12_mut[ !(LD2_gene12_mut$Position %in% All_C_indels$Position), ]
LD2_gene16_mut_wt_bg<-LD2_gene16_mut[ !(LD2_gene16_mut$Position %in% All_C_indels$Position), ]
#LD2_mutations_from_0_to_26
df_merge0 <- merge(LD2_gene0_mut_wt_bg,LD2_gene4_mut_wt_bg,by="Position")
df_merge0 <- merge(df_merge0,LD2_gene8_mut_wt_bg,by="Position")
df_merge0 <- merge(df_merge0,LD2_gene10_mut_wt_bg,by="Position")
df_merge0 <- merge(df_merge0,LD2_gene12_mut_wt_bg,by="Position")
LD2_indels_from_0_to_16 <- merge(df_merge0,LD2_gene16_mut_wt_bg,by="Position")
write.xlsx(LD2_indels_from_0_to_16,"LD2_indels_from_0_to_16.xlsx")
#Discard F0 mutations
LD2_gene4_mut_wt_gene0<-LD2_gene4_mut_wt_bg[ !(LD2_gene4_mut_wt_bg$Position %in% LD2_gene0_mut_wt_bg$Position), ]
LD2_gene8_mut_wt_gene0<-LD2_gene8_mut_wt_bg[ !(LD2_gene8_mut_wt_bg$Position %in% LD2_gene0_mut_wt_bg$Position), ]
LD2_gene10_mut_wt_gene0<-LD2_gene10_mut_wt_bg[ !(LD2_gene10_mut_wt_bg$Position %in% LD2_gene0_mut_wt_bg$Position), ]
LD2_gene12_mut_wt_gene0<-LD2_gene12_mut_wt_bg[ !(LD2_gene12_mut_wt_bg$Position %in% LD2_gene0_mut_wt_bg$Position), ]
LD2_gene16_mut_wt_gene0<-LD2_gene16_mut_wt_bg[ !(LD2_gene16_mut_wt_bg$Position %in% LD2_gene0_mut_wt_bg$Position), ]
#LD2_mutations_from_4_to_16
df_merge <- merge(LD2_gene4_mut_wt_gene0,LD2_gene8_mut_wt_gene0,by="Position")
df_merge <- merge(df_merge,LD2_gene10_mut_wt_gene0,by="Position")
df_merge <- merge(df_merge,LD2_gene12_mut_wt_gene0,by="Position")
LD2_indels_from_4_to_16 <- merge(df_merge,LD2_gene16_mut_wt_gene0,by="Position")
write.xlsx(LD2_indels_from_4_to_16,"LD2_indels_from_4_to_16.xlsx")
#Discard F4 mutations
LD2_gene8_mut_wt_gene4<-LD2_gene8_mut_wt_gene0[ !(LD2_gene8_mut_wt_gene0$Position %in% LD2_gene4_mut_wt_gene0$Position), ]
LD2_gene10_mut_wt_gene4<-LD2_gene10_mut_wt_gene0[ !(LD2_gene10_mut_wt_gene0$Position %in% LD2_gene4_mut_wt_gene0$Position), ]
LD2_gene12_mut_wt_gene4<-LD2_gene12_mut_wt_gene0[ !(LD2_gene12_mut_wt_gene0$Position %in% LD2_gene4_mut_wt_gene0$Position), ]
LD2_gene16_mut_wt_gene4<-LD2_gene16_mut_wt_gene0[ !(LD2_gene16_mut_wt_gene0$Position %in% LD2_gene4_mut_wt_gene0$Position), ]
#LD2_mutations_from_8_to_16
df_merge8 <- merge(LD2_gene8_mut_wt_gene4,LD2_gene10_mut_wt_gene4,by="Position")
df_merge8 <- merge(df_merge8,LD2_gene12_mut_wt_gene4,by="Position")
LD2_indels_from_8_to_16 <- merge(df_merge8,LD2_gene16_mut_wt_gene4,by="Position")
write.xlsx(LD2_indels_from_8_to_16,"LD2_indels_from_8_to_16.xlsx")
#Discard F8 mutations
LD2_gene10_mut_wt_gene8<-LD2_gene10_mut_wt_gene4[ !(LD2_gene10_mut_wt_gene4$Position %in% LD2_gene8_mut_wt_gene4$Position), ]
LD2_gene12_mut_wt_gene8<-LD2_gene12_mut_wt_gene4[ !(LD2_gene12_mut_wt_gene4$Position %in% LD2_gene8_mut_wt_gene4$Position), ]
LD2_gene16_mut_wt_gene8<-LD2_gene16_mut_wt_gene4[ !(LD2_gene16_mut_wt_gene4$Position %in% LD2_gene8_mut_wt_gene4$Position), ]
#LD2_mutations_from_10_to_16
df_merge10 <- merge(LD2_gene10_mut_wt_gene8,LD2_gene12_mut_wt_gene8,by="Position")
LD2_indels_from_10_to_16 <- merge(df_merge10,LD2_gene16_mut_wt_gene8,by="Position")
write.xlsx(LD2_indels_from_10_to_16,"LD2_indels_from_10_to_16.xlsx")
#Discard F10 mutations
LD2_gene12_mut_wt_gene10<-LD2_gene12_mut_wt_gene8[ !(LD2_gene12_mut_wt_gene8$Position %in% LD2_gene10_mut_wt_gene8$Position), ]
LD2_gene16_mut_wt_gene10<-LD2_gene16_mut_wt_gene8[ !(LD2_gene16_mut_wt_gene8$Position %in% LD2_gene10_mut_wt_gene8$Position), ]
#LD1_mutations_from_12_to_16
LD2_indels_from_12_to_16 <- merge(LD2_gene12_mut_wt_gene10,LD2_gene16_mut_wt_gene10,by="Position")
write.xlsx(LD2_indels_from_12_to_16,"LD2_indels_from_12_to_16.xlsx")
#Discard F12 mutations
LD2_gene16_mut_wt_gene12<-LD2_gene16_mut_wt_gene10[ !(LD2_gene16_mut_wt_gene10$Position %in% LD2_gene12_mut_wt_gene10$Position), ]
#LD1_mutations_from_16_to_16
LD2_indels_from_16_to_16 <- LD2_gene16_mut_wt_gene12
write.xlsx(LD2_indels_from_16_to_16,"LD2_indels_from_16_to_16.xlsx")
####Snps
LD2_snps_list <- list.files(path = "~/Papers/Cisplatin paper/Varscan outputs 0.5/LD2", pattern = ".snp.xlsx")
lst <- lapply(LD2_snps_list, read.xlsx)
LD2_gene0_mut<-lst[[1]]
LD2_gene4_mut<-lst[[5]]
LD2_gene8_mut<-lst[[6]]
LD2_gene10_mut<-lst[[2]]
LD2_gene12_mut<-lst[[3]]
LD2_gene16_mut<-lst[[4]]
#Discard background mutations
LD2_gene0_mut_wt_bg<-LD2_gene0_mut[ !(LD2_gene0_mut$Position %in% All_C_snps$Position), ]
LD2_gene4_mut_wt_bg<-LD2_gene4_mut[ !(LD2_gene4_mut$Position %in% All_C_snps$Position), ]
LD2_gene8_mut_wt_bg<-LD2_gene8_mut[ !(LD2_gene8_mut$Position %in% All_C_snps$Position), ]
LD2_gene10_mut_wt_bg<-LD2_gene10_mut[ !(LD2_gene10_mut$Position %in% All_C_snps$Position), ]
LD2_gene12_mut_wt_bg<-LD2_gene12_mut[ !(LD2_gene12_mut$Position %in% All_C_snps$Position), ]
LD2_gene16_mut_wt_bg<-LD2_gene16_mut[ !(LD2_gene16_mut$Position %in% All_C_snps$Position), ]
#LD2_mutations_from_0_to_20
df_merge0 <- merge(LD2_gene0_mut_wt_bg,LD2_gene4_mut_wt_bg,by="Position")
df_merge0 <- merge(df_merge0,LD2_gene8_mut_wt_bg,by="Position")
df_merge0 <- merge(df_merge0,LD2_gene10_mut_wt_bg,by="Position")
df_merge0 <- merge(df_merge0,LD2_gene12_mut_wt_bg,by="Position")
LD2_snps_from_0_to_16 <- merge(df_merge0,LD2_gene16_mut_wt_bg,by="Position")
write.xlsx(LD2_snps_from_0_to_16,"LD2_snps_from_0_to_16.xlsx")
#Discard F0 mutations
LD2_gene4_mut_wt_gene0<-LD2_gene4_mut_wt_bg[ !(LD2_gene4_mut_wt_bg$Position %in% LD2_gene0_mut_wt_bg$Position), ]
LD2_gene8_mut_wt_gene0<-LD2_gene8_mut_wt_bg[ !(LD2_gene8_mut_wt_bg$Position %in% LD2_gene0_mut_wt_bg$Position), ]
LD2_gene10_mut_wt_gene0<-LD2_gene10_mut_wt_bg[ !(LD2_gene10_mut_wt_bg$Position %in% LD2_gene0_mut_wt_bg$Position), ]
LD2_gene12_mut_wt_gene0<-LD2_gene12_mut_wt_bg[ !(LD2_gene12_mut_wt_bg$Position %in% LD2_gene0_mut_wt_bg$Position), ]
LD2_gene16_mut_wt_gene0<-LD2_gene16_mut_wt_bg[ !(LD2_gene16_mut_wt_bg$Position %in% LD2_gene0_mut_wt_bg$Position), ]
#LD2_mutations_from_4_to_16
df_merge <- merge(LD2_gene4_mut_wt_gene0,LD2_gene8_mut_wt_gene0,by="Position")
df_merge <- merge(df_merge,LD2_gene10_mut_wt_gene0,by="Position")
df_merge <- merge(df_merge,LD2_gene12_mut_wt_gene0,by="Position")
LD2_snps_from_4_to_16 <- merge(df_merge,LD2_gene16_mut_wt_gene0,by="Position")
write.xlsx(LD2_snps_from_4_to_16,"LD2_snps_from_4_to_16.xlsx")
#Discard F4 mutations
LD2_gene8_mut_wt_gene4<-LD2_gene8_mut_wt_gene0[ !(LD2_gene8_mut_wt_gene0$Position %in% LD2_gene4_mut_wt_gene0$Position), ]
LD2_gene10_mut_wt_gene4<-LD2_gene10_mut_wt_gene0[ !(LD2_gene10_mut_wt_gene0$Position %in% LD2_gene4_mut_wt_gene0$Position), ]
LD2_gene12_mut_wt_gene4<-LD2_gene12_mut_wt_gene0[ !(LD2_gene12_mut_wt_gene0$Position %in% LD2_gene4_mut_wt_gene0$Position), ]
LD2_gene16_mut_wt_gene4<-LD2_gene16_mut_wt_gene0[ !(LD2_gene16_mut_wt_gene0$Position %in% LD2_gene4_mut_wt_gene0$Position), ]
#LD2_mutations_from_8_to_16
df_merge8 <- merge(LD2_gene8_mut_wt_gene4,LD2_gene10_mut_wt_gene4,by="Position")
df_merge8 <- merge(df_merge8,LD2_gene12_mut_wt_gene4,by="Position")
LD2_snps_from_8_to_16 <- merge(df_merge8,LD2_gene16_mut_wt_gene4,by="Position")
write.xlsx(LD2_snps_from_8_to_16,"LD2_snps_from_8_to_16.xlsx")
#Discard F8 mutations
LD2_gene10_mut_wt_gene8<-LD2_gene10_mut_wt_gene4[ !(LD2_gene10_mut_wt_gene4$Position %in% LD2_gene8_mut_wt_gene4$Position), ]
LD2_gene12_mut_wt_gene6<-LD2_gene12_mut_wt_gene4[ !(LD2_gene12_mut_wt_gene4$Position %in% LD2_gene8_mut_wt_gene4$Position), ]
LD2_gene16_mut_wt_gene6<-LD2_gene16_mut_wt_gene4[ !(LD2_gene16_mut_wt_gene4$Position %in% LD2_gene8_mut_wt_gene4$Position), ]
#LD1_mutations_from_10_to_16
df_merge10 <- merge(LD2_gene10_mut_wt_gene8,LD2_gene12_mut_wt_gene6,by="Position")
LD2_snps_from_10_to_16 <- merge(df_merge10,LD2_gene16_mut_wt_gene6,by="Position")
write.xlsx(LD2_snps_from_10_to_16,"LD2_snps_from_10_to_16.xlsx")
#Discard F10 mutations
LD2_gene12_mut_wt_gene10<-LD2_gene12_mut_wt_gene6[ !(LD2_gene12_mut_wt_gene6$Position %in% LD2_gene10_mut_wt_gene8$Position), ]
LD2_gene16_mut_wt_gene10<-LD2_gene16_mut_wt_gene6[ !(LD2_gene16_mut_wt_gene6$Position %in% LD2_gene10_mut_wt_gene8$Position), ]
#LD1_mutations_from_12_to_16
LD2_snps_from_12_to_16 <- merge(LD2_gene12_mut_wt_gene10,LD2_gene16_mut_wt_gene10,by="Position")
write.xlsx(LD2_snps_from_12_to_16,"LD2_snps_from_12_to_16.xlsx")
#Discard F12 mutations
LD2_gene16_mut_wt_gene12<-LD2_gene16_mut_wt_gene10[ !(LD2_gene16_mut_wt_gene10$Position %in% LD2_gene12_mut_wt_gene10$Position), ]
#LD1_mutations_from_16_to_16
LD2_snps_from_16_to_16 <- LD2_gene16_mut_wt_gene12
write.xlsx(LD2_snps_from_16_to_16,"LD2_snps_from_16_to_16.xlsx")
####DNA mutations HD1
setwd("~/Papers/Cisplatin paper/Varscan outputs 0.5/HD1")
#####Indels
HD1_intels_list <- list.files(path = "~/Papers/Cisplatin paper/Varscan outputs 0.5/HD1", pattern = ".indel.xlsx")
lst <- lapply(HD1_intels_list, read.xlsx)
HD1_gene0_mut<-lst[[1]]
HD1_gene2_mut<-lst[[4]]
HD1_gene4_mut<-lst[[5]]
HD1_gene6_mut<-lst[[6]]
HD1_gene8_mut<-lst[[7]]
HD1_gene10_mut<-lst[[2]]
HD1_gene18_mut<-lst[[3]]
#Discard background mutations
HD1_gene0_mut_wt_bg<-HD1_gene0_mut[ !(HD1_gene0_mut$Position %in% All_C_indels$Position), ]
HD1_gene2_mut_wt_bg<-HD1_gene2_mut[ !(HD1_gene2_mut$Position %in% All_C_indels$Position), ]
HD1_gene4_mut_wt_bg<-HD1_gene4_mut[ !(HD1_gene4_mut$Position %in% All_C_indels$Position), ]
HD1_gene6_mut_wt_bg<-HD1_gene6_mut[ !(HD1_gene6_mut$Position %in% All_C_indels$Position), ]
HD1_gene8_mut_wt_bg<-HD1_gene8_mut[ !(HD1_gene8_mut$Position %in% All_C_indels$Position), ]
HD1_gene10_mut_wt_bg<-HD1_gene10_mut[ !(HD1_gene10_mut$Position %in% All_C_indels$Position), ]
HD1_gene18_mut_wt_bg<-HD1_gene18_mut[ !(HD1_gene18_mut$Position %in% All_C_indels$Position), ]
#HD1_mutations_from_0_to_18
df_merge0 <- merge(HD1_gene0_mut_wt_bg,HD1_gene2_mut_wt_bg,by="Position")
df_merge0 <- merge(df_merge0,HD1_gene4_mut_wt_bg,by="Position")
df_merge0 <- merge(df_merge0,HD1_gene6_mut_wt_bg,by="Position")
df_merge0 <- merge(df_merge0,HD1_gene8_mut_wt_bg,by="Position")
df_merge0 <- merge(df_merge0,HD1_gene10_mut_wt_bg,by="Position")
HD1_indels_from_0_to_18 <- merge(df_merge0,HD1_gene18_mut_wt_bg,by="Position")
write.xlsx(HD1_indels_from_0_to_18,"HD1_indels_from_0_to_18.xlsx")
#Discard F0 mutations
HD1_gene2_mut_wt_gene0<-HD1_gene2_mut_wt_bg[ !(HD1_gene2_mut_wt_bg$Position %in% HD1_gene0_mut_wt_bg$Position), ]
HD1_gene4_mut_wt_gene0<-HD1_gene4_mut_wt_bg[ !(HD1_gene4_mut_wt_bg$Position %in% HD1_gene0_mut_wt_bg$Position), ]
HD1_gene6_mut_wt_gene0<-HD1_gene6_mut_wt_bg[ !(HD1_gene6_mut_wt_bg$Position %in% HD1_gene0_mut_wt_bg$Position), ]
HD1_gene8_mut_wt_gene0<-HD1_gene8_mut_wt_bg[ !(HD1_gene8_mut_wt_bg$Position %in% HD1_gene0_mut_wt_bg$Position), ]
HD1_gene10_mut_wt_gene0<-HD1_gene10_mut_wt_bg[ !(HD1_gene10_mut_wt_bg$Position %in% HD1_gene0_mut_wt_bg$Position), ]
HD1_gene18_mut_wt_gene0<-HD1_gene18_mut_wt_bg[ !(HD1_gene18_mut_wt_bg$Position %in% HD1_gene0_mut_wt_bg$Position), ]
#HD1_mutations_from_2_to_18
df_merge <- merge(HD1_gene2_mut_wt_gene0,HD1_gene4_mut_wt_gene0,by="Position")
df_merge <- merge(df_merge,HD1_gene6_mut_wt_gene0,by="Position")
df_merge <- merge(df_merge,HD1_gene8_mut_wt_gene0,by="Position")
df_merge <- merge(df_merge,HD1_gene10_mut_wt_gene0,by="Position")
HD1_indels_from_2_to_18 <- merge(df_merge,HD1_gene18_mut_wt_gene0,by="Position")
write.xlsx(HD1_indels_from_2_to_18,"HD1_indels_from_2_to_18.xlsx")
#Discard F2 mutations
HD1_gene4_mut_wt_gene2<-HD1_gene4_mut_wt_gene0[ !(HD1_gene4_mut_wt_gene0$Position %in% HD1_gene2_mut_wt_gene0$Position), ]
HD1_gene6_mut_wt_gene2<-HD1_gene6_mut_wt_gene0[ !(HD1_gene6_mut_wt_gene0$Position %in% HD1_gene2_mut_wt_gene0$Position), ]
HD1_gene8_mut_wt_gene2<-HD1_gene8_mut_wt_gene0[ !(HD1_gene8_mut_wt_gene0$Position %in% HD1_gene2_mut_wt_gene0$Position), ]
HD1_gene10_mut_wt_gene2<-HD1_gene10_mut_wt_gene0[ !(HD1_gene10_mut_wt_gene0$Position %in% HD1_gene2_mut_wt_gene0$Position), ]
HD1_gene18_mut_wt_gene2<-HD1_gene18_mut_wt_gene0[ !(HD1_gene18_mut_wt_gene0$Position %in% HD1_gene2_mut_wt_gene0$Position), ]
#C1_mutations_from_4_to_18
df_merge4 <- merge(HD1_gene4_mut_wt_gene2,HD1_gene6_mut_wt_gene2,by="Position")
df_merge4 <- merge(df_merge4,HD1_gene8_mut_wt_gene2,by="Position")
df_merge4 <- merge(df_merge4,HD1_gene10_mut_wt_gene2,by="Position")
HD1_indels_from_4_to_18 <- merge(df_merge4,HD1_gene18_mut_wt_gene2,by="Position")
write.xlsx(HD1_indels_from_4_to_18,"HD1_indels_from_4_to_18.xlsx")
#Discard F4 mutations
HD1_gene6_mut_wt_gene4<-HD1_gene6_mut_wt_gene2[ !(HD1_gene6_mut_wt_gene2$Position %in% HD1_gene4_mut_wt_gene2$Position), ]
HD1_gene8_mut_wt_gene4<-HD1_gene8_mut_wt_gene2[ !(HD1_gene8_mut_wt_gene2$Position %in% HD1_gene4_mut_wt_gene2$Position), ]
HD1_gene10_mut_wt_gene4<-HD1_gene10_mut_wt_gene2[ !(HD1_gene10_mut_wt_gene2$Position %in% HD1_gene4_mut_wt_gene2$Position), ]
HD1_gene18_mut_wt_gene4<-HD1_gene18_mut_wt_gene2[ !(HD1_gene18_mut_wt_gene2$Position %in% HD1_gene4_mut_wt_gene2$Position), ]
#C1_mutations_from_6_to_18
df_merge6 <- merge(HD1_gene6_mut_wt_gene4,HD1_gene8_mut_wt_gene4,by="Position")
df_merge6 <- merge(df_merge6,HD1_gene10_mut_wt_gene4,by="Position")
HD1_indels_from_6_to_18 <- merge(df_merge6,HD1_gene18_mut_wt_gene4,by="Position")
write.xlsx(HD1_indels_from_6_to_18,"HD1_indels_from_6_to_18.xlsx")
#Discard F6 mutations
HD1_gene8_mut_wt_gene6<-HD1_gene8_mut_wt_gene4[ !(HD1_gene8_mut_wt_gene4$Position %in% HD1_gene6_mut_wt_gene4$Position), ]
HD1_gene10_mut_wt_gene6<-HD1_gene10_mut_wt_gene4[ !(HD1_gene10_mut_wt_gene4$Position %in% HD1_gene6_mut_wt_gene4$Position), ]
HD1_gene18_mut_wt_gene6<-HD1_gene18_mut_wt_gene4[ !(HD1_gene18_mut_wt_gene4$Position %in% HD1_gene6_mut_wt_gene4$Position), ]
#C1_mutations_from_8_to_18
df_merge8 <- merge(HD1_gene8_mut_wt_gene6,HD1_gene10_mut_wt_gene6,by="Position")
HD1_indels_from_8_to_18 <- merge(df_merge8,HD1_gene18_mut_wt_gene6,by="Position")
write.xlsx(HD1_indels_from_8_to_18,"HD1_indels_from_8_to_18.xlsx")
#Discard F8 mutations
HD1_gene10_mut_wt_gene8<-HD1_gene10_mut_wt_gene6[ !(HD1_gene10_mut_wt_gene6$Position %in% HD1_gene8_mut_wt_gene6$Position), ]
HD1_gene18_mut_wt_gene8<-HD1_gene18_mut_wt_gene6[ !(HD1_gene18_mut_wt_gene6$Position %in% HD1_gene8_mut_wt_gene6$Position), ]
#HD1_mutations_from_10_to_18
HD1_indels_from_10_to_18 <- merge(HD1_gene10_mut_wt_gene8,HD1_gene18_mut_wt_gene8,by="Position")
write.xlsx(HD1_indels_from_10_to_18,"HD1_indels_from_10_to_18.xlsx")
#Discard F10 mutations
HD1_gene18_mut_wt_gene10<-HD1_gene18_mut_wt_gene8[ !(HD1_gene18_mut_wt_gene8$Position %in% HD1_gene10_mut_wt_gene8$Position), ]
#HD1_mutations_from_18_to_18
HD1_indels_from_18_to_18 <- HD1_gene18_mut_wt_gene10
write.xlsx(HD1_indels_from_18_to_18,"HD1_indels_from_18_to_18.xlsx")
####Snps
HD1_snps_list <- list.files(path = "~/Papers/Cisplatin paper/Varscan outputs 0.5/HD1", pattern = ".snp.xlsx")
lst <- lapply(HD1_snps_list, read.xlsx)
HD1_gene0_snps<-lst[[1]]
HD1_gene2_snps<-lst[[4]]
HD1_gene4_snps<-lst[[5]]
HD1_gene6_snps<-lst[[6]]
HD1_gene8_snps<-lst[[7]]
HD1_gene10_snps<-lst[[2]]
HD1_gene18_snps<-lst[[3]]
#Discard background mutations
HD1_gene0_snps_wt_bg<-HD1_gene0_snps[ !(HD1_gene0_snps$Position %in% All_C_snps$Position), ]
HD1_gene2_snps_wt_bg<-HD1_gene2_snps[ !(HD1_gene2_snps$Position %in% All_C_snps$Position), ]
HD1_gene4_snps_wt_bg<-HD1_gene4_snps[ !(HD1_gene4_snps$Position %in% All_C_snps$Position), ]
HD1_gene6_snps_wt_bg<-HD1_gene6_snps[ !(HD1_gene6_snps$Position %in% All_C_snps$Position), ]
HD1_gene8_snps_wt_bg<-HD1_gene8_snps[ !(HD1_gene8_snps$Position %in% All_C_snps$Position), ]
HD1_gene10_snps_wt_bg<-HD1_gene10_snps[ !(HD1_gene10_snps$Position %in% All_C_snps$Position), ]
HD1_gene18_snps_wt_bg<-HD1_gene18_snps[ !(HD1_gene18_snps$Position %in% All_C_snps$Position), ]
#HD1_snps_from_0_to_18
df_merge <- merge(HD1_gene0_snps_wt_bg,HD1_gene2_snps_wt_bg,by="Position")
df_merge <- merge(df_merge,HD1_gene4_snps_wt_bg,by="Position")
df_merge <- merge(df_merge,HD1_gene6_snps_wt_bg,by="Position")
df_merge <- merge(df_merge,HD1_gene8_snps_wt_bg,by="Position")
df_merge <- merge(df_merge,HD1_gene10_snps_wt_bg,by="Position")
HD1_snps_from_0_to_18 <- merge(df_merge,HD1_gene18_snps_wt_bg,by="Position")
write.xlsx(HD1_snps_from_0_to_18,"HD1_snps_from_0_to_18.xlsx")
#Discard F0 mutations
HD1_gene2_snps_wt_gene0<-HD1_gene2_snps_wt_bg[ !(HD1_gene2_snps_wt_bg$Position %in% HD1_gene0_snps_wt_bg$Position), ]
HD1_gene4_snps_wt_gene0<-HD1_gene4_snps_wt_bg[ !(HD1_gene4_snps_wt_bg$Position %in% HD1_gene0_snps_wt_bg$Position), ]
HD1_gene6_snps_wt_gene0<-HD1_gene6_snps_wt_bg[ !(HD1_gene6_snps_wt_bg$Position %in% HD1_gene0_snps_wt_bg$Position), ]
HD1_gene8_snps_wt_gene0<-HD1_gene8_snps_wt_bg[ !(HD1_gene8_snps_wt_bg$Position %in% HD1_gene0_snps_wt_bg$Position), ]
HD1_gene10_snps_wt_gene0<-HD1_gene10_snps_wt_bg[ !(HD1_gene10_snps_wt_bg$Position %in% HD1_gene0_snps_wt_bg$Position), ]
HD1_gene18_snps_wt_gene0<-HD1_gene18_snps_wt_bg[ !(HD1_gene18_snps_wt_bg$Position %in% HD1_gene0_snps_wt_bg$Position), ]
#HD1_snps_from_2_to_18
df_merge <- merge(HD1_gene2_snps_wt_gene0,HD1_gene4_snps_wt_gene0,by="Position")
df_merge <- merge(df_merge,HD1_gene6_snps_wt_gene0,by="Position")
df_merge <- merge(df_merge,HD1_gene8_snps_wt_gene0,by="Position")
df_merge <- merge(df_merge,HD1_gene10_snps_wt_gene0,by="Position")
HD1_snps_from_2_to_18 <- merge(df_merge,HD1_gene18_snps_wt_gene0,by="Position")
write.xlsx(HD1_snps_from_2_to_18,"HD1_snps_from_2_to_18.xlsx")
#Discard F2 snps
HD1_gene4_snps_wt_gene2<-HD1_gene4_snps_wt_gene0[ !(HD1_gene4_snps_wt_gene0$Position %in% HD1_gene2_snps_wt_gene0$Position), ]
HD1_gene6_snps_wt_gene2<-HD1_gene6_snps_wt_gene0[ !(HD1_gene6_snps_wt_gene0$Position %in% HD1_gene2_snps_wt_gene0$Position), ]
HD1_gene8_snps_wt_gene2<-HD1_gene8_snps_wt_gene0[ !(HD1_gene8_snps_wt_gene0$Position %in% HD1_gene2_snps_wt_gene0$Position), ]
HD1_gene10_snps_wt_gene2<-HD1_gene10_snps_wt_gene0[ !(HD1_gene10_snps_wt_gene0$Position %in% HD1_gene2_snps_wt_gene0$Position), ]
HD1_gene18_snps_wt_gene2<-HD1_gene18_snps_wt_gene0[ !(HD1_gene18_snps_wt_gene0$Position %in% HD1_gene2_snps_wt_gene0$Position), ]
#C1_snps_from_4_to_18
df_merge4 <- merge(HD1_gene4_snps_wt_gene2,HD1_gene6_snps_wt_gene2,by="Position")
df_merge4 <- merge(df_merge4,HD1_gene8_snps_wt_gene2,by="Position")
df_merge4 <- merge(df_merge4,HD1_gene10_snps_wt_gene2,by="Position")
HD1_snps_from_4_to_18 <- merge(df_merge4,HD1_gene18_snps_wt_gene2,by="Position")
write.xlsx(HD1_snps_from_4_to_18,"HD1_snps_from_4_to_18.xlsx")
#Discard F4 snps
HD1_gene6_snps_wt_gene4<-HD1_gene6_snps_wt_gene2[ !(HD1_gene6_snps_wt_gene2$Position %in% HD1_gene4_snps_wt_gene2$Position), ]
HD1_gene8_snps_wt_gene4<-HD1_gene8_snps_wt_gene2[ !(HD1_gene8_snps_wt_gene2$Position %in% HD1_gene4_snps_wt_gene2$Position), ]
HD1_gene10_snps_wt_gene4<-HD1_gene10_snps_wt_gene2[ !(HD1_gene10_snps_wt_gene2$Position %in% HD1_gene4_snps_wt_gene2$Position), ]
HD1_gene18_snps_wt_gene4<-HD1_gene18_snps_wt_gene2[ !(HD1_gene18_snps_wt_gene2$Position %in% HD1_gene4_snps_wt_gene2$Position), ]
#C1_snps_from_6_to_18
df_merge6 <- merge(HD1_gene6_snps_wt_gene4,HD1_gene8_snps_wt_gene4,by="Position")
df_merge6 <- merge(df_merge6,HD1_gene10_snps_wt_gene4,by="Position")
HD1_snps_from_6_to_18 <- merge(df_merge6,HD1_gene18_snps_wt_gene4,by="Position")
write.xlsx(HD1_snps_from_6_to_18,"HD1_snps_from_6_to_18.xlsx")
#Discard F6 snps
HD1_gene8_snps_wt_gene6<-HD1_gene8_snps_wt_gene4[ !(HD1_gene8_snps_wt_gene4$Position %in% HD1_gene6_snps_wt_gene4$Position), ]
HD1_gene10_snps_wt_gene6<-HD1_gene10_snps_wt_gene4[ !(HD1_gene10_snps_wt_gene4$Position %in% HD1_gene6_snps_wt_gene4$Position), ]
HD1_gene18_snps_wt_gene6<-HD1_gene18_snps_wt_gene4[ !(HD1_gene18_snps_wt_gene4$Position %in% HD1_gene6_snps_wt_gene4$Position), ]
#C1_snps_from_8_to_18
df_merge8 <- merge(HD1_gene8_snps_wt_gene6,HD1_gene10_snps_wt_gene6,by="Position")
HD1_snps_from_8_to_18 <- merge(df_merge8,HD1_gene18_snps_wt_gene6,by="Position")
write.xlsx(HD1_snps_from_8_to_18,"HD1_snps_from_8_to_18.xlsx")
#Discard F8 snps
HD1_gene10_snps_wt_gene8<-HD1_gene10_snps_wt_gene6[ !(HD1_gene10_snps_wt_gene6$Position %in% HD1_gene8_snps_wt_gene6$Position), ]
HD1_gene18_snps_wt_gene8<-HD1_gene18_snps_wt_gene6[ !(HD1_gene18_snps_wt_gene6$Position %in% HD1_gene8_snps_wt_gene6$Position), ]
#HD1_snps_from_10_to_18
HD1_snps_from_10_to_18 <- merge(HD1_gene10_snps_wt_gene8,HD1_gene18_snps_wt_gene8,by="Position")
write.xlsx(HD1_snps_from_10_to_18,"HD1_snps_from_10_to_18.xlsx")
#Discard F10 snps
HD1_gene18_snps_wt_gene10<-HD1_gene18_snps_wt_gene8[ !(HD1_gene18_snps_wt_gene8$Position %in% HD1_gene10_snps_wt_gene8$Position), ]
#HD1_snps_from_10_to_18
HD1_snps_from_18_to_18 <- HD1_gene18_snps_wt_gene10
write.xlsx(HD1_snps_from_18_to_18,"HD1_snps_from_18_to_18.xlsx")
#####DNA mutations HD2
setwd("~/Papers/Cisplatin paper/Varscan outputs 0.5/HD2")
####Indels
HD2_intels_list <- list.files(path = "~/Papers/Cisplatin paper/Varscan outputs 0.5/HD2", pattern = ".indel.xlsx")
lst <- lapply(HD2_intels_list, read.xlsx)
HD2_gene0_mut<-lst[[1]]
HD2_gene2_mut<-lst[[6]]
HD2_gene4_mut<-lst[[7]]
HD2_gene8_mut<-lst[[8]]
HD2_gene10_mut<-lst[[2]]
HD2_gene12_mut<-lst[[3]]
HD2_gene16_mut<-lst[[4]]
HD2_gene18_mut<-lst[[5]]
#Discard background mutations
HD2_gene0_mut_wt_bg<-HD2_gene0_mut[ !(HD2_gene0_mut$Position %in% All_C_indels$Position), ]
HD2_gene2_mut_wt_bg<-HD2_gene2_mut[ !(HD2_gene2_mut$Position %in% All_C_indels$Position), ]
HD2_gene4_mut_wt_bg<-HD2_gene4_mut[ !(HD2_gene4_mut$Position %in% All_C_indels$Position), ]
HD2_gene8_mut_wt_bg<-HD2_gene8_mut[ !(HD2_gene8_mut$Position %in% All_C_indels$Position), ]
HD2_gene10_mut_wt_bg<-HD2_gene10_mut[ !(HD2_gene10_mut$Position %in% All_C_indels$Position), ]
HD2_gene12_mut_wt_bg<-HD2_gene12_mut[ !(HD2_gene12_mut$Position %in% All_C_indels$Position), ]
HD2_gene16_mut_wt_bg<-HD2_gene16_mut[ !(HD2_gene16_mut$Position %in% All_C_indels$Position), ]
HD2_gene18_mut_wt_bg<-HD2_gene18_mut[ !(HD2_gene18_mut$Position %in% All_C_indels$Position), ]
#HD2_mutations_from_0_to_18
df_merge0 <- merge(HD2_gene0_mut_wt_bg,HD2_gene2_mut_wt_bg,by="Position")
df_merge0 <- merge(df_merge0,HD2_gene4_mut_wt_bg,by="Position")
df_merge0 <- merge(df_merge0,HD2_gene8_mut_wt_bg,by="Position")
df_merge0 <- merge(df_merge0,HD2_gene10_mut_wt_bg,by="Position")
df_merge0 <- merge(df_merge0,HD2_gene12_mut_wt_bg,by="Position")
df_merge0 <- merge(df_merge0,HD2_gene16_mut_wt_bg,by="Position")
HD2_indels_from_0_to_18 <- merge(df_merge0,HD2_gene18_mut_wt_bg,by="Position")
write.xlsx(HD2_indels_from_0_to_18,"HD2_indels_from_0_to_18.xlsx")
#Discard F0 mutations
HD2_gene2_mut_wt_gene0<-HD2_gene2_mut_wt_bg[ !(HD2_gene2_mut_wt_bg$Position %in% HD2_gene0_mut_wt_bg$Position), ]
HD2_gene4_mut_wt_gene0<-HD2_gene4_mut_wt_bg[ !(HD2_gene4_mut_wt_bg$Position %in% HD2_gene0_mut_wt_bg$Position), ]
HD2_gene8_mut_wt_gene0<-HD2_gene8_mut_wt_bg[ !(HD2_gene8_mut_wt_bg$Position %in% HD2_gene0_mut_wt_bg$Position), ]
HD2_gene10_mut_wt_gene0<-HD2_gene10_mut_wt_bg[ !(HD2_gene10_mut_wt_bg$Position %in% HD2_gene0_mut_wt_bg$Position), ]
HD2_gene12_mut_wt_gene0<-HD2_gene12_mut_wt_bg[ !(HD2_gene12_mut_wt_bg$Position %in% HD2_gene0_mut_wt_bg$Position), ]
HD2_gene16_mut_wt_gene0<-HD2_gene16_mut_wt_bg[ !(HD2_gene16_mut_wt_bg$Position %in% HD2_gene0_mut_wt_bg$Position), ]
HD2_gene18_mut_wt_gene0<-HD2_gene18_mut_wt_bg[ !(HD2_gene18_mut_wt_bg$Position %in% HD2_gene0_mut_wt_bg$Position), ]
#HD2_mutations_from_2_to_18
df_merge <- merge(HD2_gene2_mut_wt_gene0,HD2_gene4_mut_wt_gene0,by="Position")
df_merge <- merge(df_merge,HD2_gene8_mut_wt_gene0,by="Position")
df_merge <- merge(df_merge,HD2_gene10_mut_wt_gene0,by="Position")
df_merge <- merge(df_merge,HD2_gene12_mut_wt_gene0,by="Position")
df_merge <- merge(df_merge,HD2_gene16_mut_wt_gene0,by="Position")
HD2_indels_from_2_to_18 <- merge(df_merge,HD2_gene18_mut_wt_gene0,by="Position")
write.xlsx(HD2_indels_from_2_to_18,"HD2_indels_from_2_to_18.xlsx")
#Discard F2 mutations
HD2_gene4_mut_wt_gene2<-HD2_gene4_mut_wt_gene0[ !(HD2_gene4_mut_wt_gene0$Position %in% HD2_gene2_mut_wt_gene0$Position), ]
HD2_gene8_mut_wt_gene2<-HD2_gene8_mut_wt_gene0[ !(HD2_gene8_mut_wt_gene0$Position %in% HD2_gene2_mut_wt_gene0$Position), ]
HD2_gene10_mut_wt_gene2<-HD2_gene10_mut_wt_gene0[ !(HD2_gene10_mut_wt_gene0$Position %in% HD2_gene2_mut_wt_gene0$Position), ]
HD2_gene12_mut_wt_gene2<-HD2_gene12_mut_wt_gene0[ !(HD2_gene12_mut_wt_gene0$Position %in% HD2_gene2_mut_wt_gene0$Position), ]
HD2_gene16_mut_wt_gene2<-HD2_gene16_mut_wt_gene0[ !(HD2_gene16_mut_wt_gene0$Position %in% HD2_gene2_mut_wt_gene0$Position), ]
HD2_gene18_mut_wt_gene2<-HD2_gene18_mut_wt_gene0[ !(HD2_gene18_mut_wt_gene0$Position %in% HD2_gene2_mut_wt_gene0$Position), ]
#HD2_mutations_from_4_to_18
df_merge4 <- merge(HD2_gene4_mut_wt_gene2,HD2_gene8_mut_wt_gene2,by="Position")
df_merge4 <- merge(df_merge4,HD2_gene10_mut_wt_gene2,by="Position")
df_merge4 <- merge(df_merge4,HD2_gene12_mut_wt_gene2,by="Position")
df_merge4 <- merge(df_merge4,HD2_gene16_mut_wt_gene2,by="Position")
HD2_indels_from_4_to_18 <- merge(df_merge4,HD2_gene18_mut_wt_gene2,by="Position")
write.xlsx(HD2_indels_from_4_to_18,"HD2_indels_from_4_to_18.xlsx")
#Discard F4 mutations
HD2_gene8_mut_wt_gene4<-HD2_gene8_mut_wt_gene2[ !(HD2_gene8_mut_wt_gene2$Position %in% HD2_gene4_mut_wt_gene2$Position), ]
HD2_gene10_mut_wt_gene4<-HD2_gene10_mut_wt_gene2[ !(HD2_gene10_mut_wt_gene2$Position %in% HD2_gene4_mut_wt_gene2$Position), ]
HD2_gene12_mut_wt_gene4<-HD2_gene12_mut_wt_gene2[ !(HD2_gene12_mut_wt_gene2$Position %in% HD2_gene4_mut_wt_gene2$Position), ]
HD2_gene16_mut_wt_gene4<-HD2_gene16_mut_wt_gene2[ !(HD2_gene16_mut_wt_gene2$Position %in% HD2_gene4_mut_wt_gene2$Position), ]
HD2_gene18_mut_wt_gene4<-HD2_gene18_mut_wt_gene2[ !(HD2_gene18_mut_wt_gene2$Position %in% HD2_gene4_mut_wt_gene2$Position), ]
#HD2_mutations_from_8_to_18
df_merge8 <- merge(HD2_gene8_mut_wt_gene4,HD2_gene10_mut_wt_gene4,by="Position")
df_merge8 <- merge(df_merge8,HD2_gene12_mut_wt_gene4,by="Position")
df_merge8 <- merge(df_merge8,HD2_gene16_mut_wt_gene4,by="Position")
HD2_indels_from_8_to_18 <- merge(df_merge8,HD2_gene18_mut_wt_gene4,by="Position")
write.xlsx(HD2_indels_from_8_to_18,"HD2_indels_from_8_to_18.xlsx")
#Discard F8 mutations
HD2_gene10_mut_wt_gene6<-HD2_gene10_mut_wt_gene4[ !(HD2_gene10_mut_wt_gene4$Position %in% HD2_gene8_mut_wt_gene4$Position), ]
HD2_gene12_mut_wt_gene6<-HD2_gene12_mut_wt_gene4[ !(HD2_gene12_mut_wt_gene4$Position %in% HD2_gene8_mut_wt_gene4$Position), ]
HD2_gene16_mut_wt_gene6<-HD2_gene16_mut_wt_gene4[ !(HD2_gene16_mut_wt_gene4$Position %in% HD2_gene8_mut_wt_gene4$Position), ]
HD2_gene18_mut_wt_gene6<-HD2_gene18_mut_wt_gene4[ !(HD2_gene18_mut_wt_gene4$Position %in% HD2_gene8_mut_wt_gene4$Position), ]
#HD2_mutations_from_10_to_18
df_merge10 <- merge(HD2_gene10_mut_wt_gene6,HD2_gene12_mut_wt_gene6,by="Position")
df_merge10 <- merge(df_merge10,HD2_gene16_mut_wt_gene6,by="Position")
HD2_indels_from_10_to_18 <- merge(df_merge10,HD2_gene18_mut_wt_gene6,by="Position")
write.xlsx(HD2_indels_from_10_to_18,"HD2_indels_from_10_to_18.xlsx")
#Discard F10 mutations
HD2_gene12_mut_wt_gene10<-HD2_gene12_mut_wt_gene6[ !(HD2_gene12_mut_wt_gene6$Position %in% HD2_gene10_mut_wt_gene6$Position), ]
HD2_gene16_mut_wt_gene10<-HD2_gene16_mut_wt_gene6[ !(HD2_gene16_mut_wt_gene6$Position %in% HD2_gene10_mut_wt_gene6$Position), ]
HD2_gene18_mut_wt_gene10<-HD2_gene18_mut_wt_gene6[ !(HD2_gene18_mut_wt_gene6$Position %in% HD2_gene10_mut_wt_gene6$Position), ]
#HD2_mutations_from_12_to_18
df_merge12 <- merge(HD2_gene12_mut_wt_gene10,HD2_gene16_mut_wt_gene10,by="Position")
HD2_indels_from_12_to_18 <- merge(df_merge12,HD2_gene18_mut_wt_gene10,by="Position")
write.xlsx(HD2_indels_from_12_to_18,"HD2_indels_from_12_to_18.xlsx")
#Discard F12 mutations
HD2_gene16_mut_wt_gene12<-HD2_gene16_mut_wt_gene10[ !(HD2_gene16_mut_wt_gene10$Position %in% HD2_gene12_mut_wt_gene10$Position), ]
HD2_gene18_mut_wt_gene12<-HD2_gene18_mut_wt_gene10[ !(HD2_gene18_mut_wt_gene10$Position %in% HD2_gene12_mut_wt_gene10$Position), ]
#HD2_mutations_from_16_to_18
HD2_indels_from_16_to_18 <- merge(HD2_gene16_mut_wt_gene12,HD2_gene18_mut_wt_gene12,by="Position")
write.xlsx(HD2_indels_from_16_to_18,"HD2_indels_from_16_to_18.xlsx")
#Discard F16 mutations
HD2_gene18_mut_wt_gene16<-HD2_gene18_mut_wt_gene12[ !(HD2_gene18_mut_wt_gene12$Position %in% HD2_gene16_mut_wt_gene12$Position), ]
#HD2_mutations_from_18_to_18
HD2_indels_from_18_to_18 <- HD2_gene18_mut_wt_gene16
write.xlsx(HD2_indels_from_18_to_18,"HD2_indels_from_18_to_18.xlsx")
#####Snps
HD2_snps_list <- list.files(path = "~/Papers/Cisplatin paper/Varscan outputs 0.5/HD2", pattern = ".snp.xlsx")
lst <- lapply(HD2_snps_list, read.xlsx)
HD2_gene0_snps<-lst[[1]]
HD2_gene2_snps<-lst[[6]]
HD2_gene4_snps<-lst[[7]]
HD2_gene8_snps<-lst[[8]]
HD2_gene10_snps<-lst[[2]]
HD2_gene12_snps<-lst[[3]]
HD2_gene16_snps<-lst[[4]]
HD2_gene18_snps<-lst[[5]]
#Discard background snps
HD2_gene0_snps_wt_bg<-HD2_gene0_snps[ !(HD2_gene0_snps$Position %in% All_C_snps$Position), ]
HD2_gene2_snps_wt_bg<-HD2_gene2_snps[ !(HD2_gene2_snps$Position %in% All_C_snps$Position), ]
HD2_gene4_snps_wt_bg<-HD2_gene4_snps[ !(HD2_gene4_snps$Position %in% All_C_snps$Position), ]
HD2_gene8_snps_wt_bg<-HD2_gene8_snps[ !(HD2_gene8_snps$Position %in% All_C_snps$Position), ]
HD2_gene10_snps_wt_bg<-HD2_gene10_snps[ !(HD2_gene10_snps$Position %in% All_C_snps$Position), ]
HD2_gene12_snps_wt_bg<-HD2_gene12_snps[ !(HD2_gene12_snps$Position %in% All_C_snps$Position), ]
HD2_gene16_snps_wt_bg<-HD2_gene16_snps[ !(HD2_gene16_snps$Position %in% All_C_snps$Position), ]
HD2_gene18_snps_wt_bg<-HD2_gene18_snps[ !(HD2_gene18_snps$Position %in% All_C_snps$Position), ]
#HD2_snps_from_0_to_18
df_merge0 <- merge(HD2_gene0_snps_wt_bg,HD2_gene2_snps_wt_bg,by="Position")
df_merge0 <- merge(df_merge0,HD2_gene4_snps_wt_bg,by="Position")
df_merge0 <- merge(df_merge0,HD2_gene8_snps_wt_bg,by="Position")
df_merge0 <- merge(df_merge0,HD2_gene10_snps_wt_bg,by="Position")
df_merge0 <- merge(df_merge0,HD2_gene12_snps_wt_bg,by="Position")
df_merge0 <- merge(df_merge0,HD2_gene16_snps_wt_bg,by="Position")
HD2_snps_from_0_to_18 <- merge(df_merge0,HD2_gene18_snps_wt_bg,by="Position")
write.xlsx(HD2_snps_from_0_to_18,"HD2_snps_from_0_to_18.xlsx")
#Discard F0 snps
HD2_gene2_snps_wt_gene0<-HD2_gene2_snps_wt_bg[ !(HD2_gene2_snps_wt_bg$Position %in% HD2_gene0_snps_wt_bg$Position), ]
HD2_gene4_snps_wt_gene0<-HD2_gene4_snps_wt_bg[ !(HD2_gene4_snps_wt_bg$Position %in% HD2_gene0_snps_wt_bg$Position), ]
HD2_gene8_snps_wt_gene0<-HD2_gene8_snps_wt_bg[ !(HD2_gene8_snps_wt_bg$Position %in% HD2_gene0_snps_wt_bg$Position), ]
HD2_gene10_snps_wt_gene0<-HD2_gene10_snps_wt_bg[ !(HD2_gene10_snps_wt_bg$Position %in% HD2_gene0_snps_wt_bg$Position), ]
HD2_gene12_snps_wt_gene0<-HD2_gene12_snps_wt_bg[ !(HD2_gene12_snps_wt_bg$Position %in% HD2_gene0_snps_wt_bg$Position), ]
HD2_gene16_snps_wt_gene0<-HD2_gene16_snps_wt_bg[ !(HD2_gene16_snps_wt_bg$Position %in% HD2_gene0_snps_wt_bg$Position), ]
HD2_gene18_snps_wt_gene0<-HD2_gene18_snps_wt_bg[ !(HD2_gene18_snps_wt_bg$Position %in% HD2_gene0_snps_wt_bg$Position), ]
#HD2_snps_from_2_to_18
df_merge <- merge(HD2_gene2_snps_wt_gene0,HD2_gene4_snps_wt_gene0,by="Position")
df_merge <- merge(df_merge,HD2_gene8_snps_wt_gene0,by="Position")
df_merge <- merge(df_merge,HD2_gene10_snps_wt_gene0,by="Position")
df_merge <- merge(df_merge,HD2_gene12_snps_wt_gene0,by="Position")
df_merge <- merge(df_merge,HD2_gene16_snps_wt_gene0,by="Position")
HD2_snps_from_2_to_18 <- merge(df_merge,HD2_gene18_snps_wt_gene0,by="Position")
write.xlsx(HD2_snps_from_2_to_18,"HD2_snps_from_2_to_18.xlsx")
#Discard F2 snps
HD2_gene4_snps_wt_gene2<-HD2_gene4_snps_wt_gene0[ !(HD2_gene4_snps_wt_gene0$Position %in% HD2_gene2_snps_wt_gene0$Position), ]
HD2_gene8_snps_wt_gene2<-HD2_gene8_snps_wt_gene0[ !(HD2_gene8_snps_wt_gene0$Position %in% HD2_gene2_snps_wt_gene0$Position), ]
HD2_gene10_snps_wt_gene2<-HD2_gene10_snps_wt_gene0[ !(HD2_gene10_snps_wt_gene0$Position %in% HD2_gene2_snps_wt_gene0$Position), ]
HD2_gene12_snps_wt_gene2<-HD2_gene12_snps_wt_gene0[ !(HD2_gene12_snps_wt_gene0$Position %in% HD2_gene2_snps_wt_gene0$Position), ]
HD2_gene16_snps_wt_gene2<-HD2_gene16_snps_wt_gene0[ !(HD2_gene16_snps_wt_gene0$Position %in% HD2_gene2_snps_wt_gene0$Position), ]
HD2_gene18_snps_wt_gene2<-HD2_gene18_snps_wt_gene0[ !(HD2_gene18_snps_wt_gene0$Position %in% HD2_gene2_snps_wt_gene0$Position), ]
#HD2_snps_from_4_to_18
df_merge4 <- merge(HD2_gene4_snps_wt_gene2,HD2_gene8_snps_wt_gene2,by="Position")
df_merge4 <- merge(df_merge4,HD2_gene10_snps_wt_gene2,by="Position")
df_merge4 <- merge(df_merge4,HD2_gene12_snps_wt_gene2,by="Position")
df_merge4 <- merge(df_merge4,HD2_gene16_snps_wt_gene2,by="Position")
HD2_snps_from_4_to_18 <- merge(df_merge4,HD2_gene18_snps_wt_gene2,by="Position")
write.xlsx(HD2_snps_from_4_to_18,"HD2_snps_from_4_to_18.xlsx")
#Discard F4 snps
HD2_gene8_snps_wt_gene4<-HD2_gene8_snps_wt_gene2[ !(HD2_gene8_snps_wt_gene2$Position %in% HD2_gene4_snps_wt_gene2$Position), ]
HD2_gene10_snps_wt_gene4<-HD2_gene10_snps_wt_gene2[ !(HD2_gene10_snps_wt_gene2$Position %in% HD2_gene4_snps_wt_gene2$Position), ]
HD2_gene12_snps_wt_gene4<-HD2_gene12_snps_wt_gene2[ !(HD2_gene12_snps_wt_gene2$Position %in% HD2_gene4_snps_wt_gene2$Position), ]
HD2_gene16_snps_wt_gene4<-HD2_gene16_snps_wt_gene2[ !(HD2_gene16_snps_wt_gene2$Position %in% HD2_gene4_snps_wt_gene2$Position), ]
HD2_gene18_snps_wt_gene4<-HD2_gene18_snps_wt_gene2[ !(HD2_gene18_snps_wt_gene2$Position %in% HD2_gene4_snps_wt_gene2$Position), ]
#HD2_snps_from_8_to_18
df_merge8 <- merge(HD2_gene8_snps_wt_gene4,HD2_gene10_snps_wt_gene4,by="Position")
df_merge8 <- merge(df_merge8,HD2_gene12_snps_wt_gene4,by="Position")
df_merge8 <- merge(df_merge8,HD2_gene16_snps_wt_gene4,by="Position")
HD2_snps_from_8_to_18 <- merge(df_merge8,HD2_gene18_snps_wt_gene4,by="Position")
write.xlsx(HD2_snps_from_8_to_18,"HD2_snps_from_8_to_18.xlsx")
#Discard F8 snps
HD2_gene10_snps_wt_gene6<-HD2_gene10_snps_wt_gene4[ !(HD2_gene10_snps_wt_gene4$Position %in% HD2_gene8_snps_wt_gene4$Position), ]
HD2_gene12_snps_wt_gene6<-HD2_gene12_snps_wt_gene4[ !(HD2_gene12_snps_wt_gene4$Position %in% HD2_gene8_snps_wt_gene4$Position), ]
HD2_gene16_snps_wt_gene6<-HD2_gene16_snps_wt_gene4[ !(HD2_gene16_snps_wt_gene4$Position %in% HD2_gene8_snps_wt_gene4$Position), ]
HD2_gene18_snps_wt_gene6<-HD2_gene18_snps_wt_gene4[ !(HD2_gene18_snps_wt_gene4$Position %in% HD2_gene8_snps_wt_gene4$Position), ]
#HD2_snps_from_10_to_18
df_merge10 <- merge(HD2_gene10_snps_wt_gene6,HD2_gene12_snps_wt_gene6,by="Position")
df_merge10 <- merge(df_merge10,HD2_gene16_snps_wt_gene6,by="Position")
HD2_snps_from_10_to_18 <- merge(df_merge10,HD2_gene18_snps_wt_gene6,by="Position")
write.xlsx(HD2_snps_from_10_to_18,"HD2_snps_from_10_to_18.xlsx")
#Discard F10 snps
HD2_gene12_snps_wt_gene10<-HD2_gene12_snps_wt_gene6[ !(HD2_gene12_snps_wt_gene6$Position %in% HD2_gene10_snps_wt_gene6$Position), ]
HD2_gene16_snps_wt_gene10<-HD2_gene16_snps_wt_gene6[ !(HD2_gene16_snps_wt_gene6$Position %in% HD2_gene10_snps_wt_gene6$Position), ]
HD2_gene18_snps_wt_gene10<-HD2_gene18_snps_wt_gene6[ !(HD2_gene18_snps_wt_gene6$Position %in% HD2_gene10_snps_wt_gene6$Position), ]
#HD2_snps_from_12_to_18
df_merge12 <- merge(HD2_gene12_snps_wt_gene10,HD2_gene16_snps_wt_gene10,by="Position")
HD2_snps_from_12_to_18 <- merge(df_merge12,HD2_gene18_snps_wt_gene10,by="Position")
write.xlsx(HD2_snps_from_12_to_18,"HD2_snps_from_12_to_18.xlsx")
#Discard F12 snps
HD2_gene16_snps_wt_gene12<-HD2_gene16_snps_wt_gene10[ !(HD2_gene16_snps_wt_gene10$Position %in% HD2_gene12_snps_wt_gene10$Position), ]
HD2_gene18_snps_wt_gene12<-HD2_gene18_snps_wt_gene10[ !(HD2_gene18_snps_wt_gene10$Position %in% HD2_gene12_snps_wt_gene10$Position), ]
#HD2_snps_from_16_to_18
HD2_snps_from_16_to_18 <- merge(HD2_gene16_snps_wt_gene12,HD2_gene18_snps_wt_gene12,by="Position")
write.xlsx(HD2_snps_from_16_to_18,"HD2_snps_from_16_to_18.xlsx")
#Discard F16 snps
HD2_gene18_snps_wt_gene16<-HD2_gene18_snps_wt_gene12[ !(HD2_gene18_snps_wt_gene12$Position %in% HD2_gene16_snps_wt_gene12$Position), ]
#HD2_snps_from_18_to_18
HD2_snps_from_18_to_18 <- HD2_gene18_snps_wt_gene16
write.xlsx(HD2_snps_from_18_to_18,"HD2_snps_from_18_to_18.xlsx")

setwd("~/Documents/Cisplatin project/DNA mutations")
####All mutations
#Import table manually created from previous generated tables (exp:HD2_snps_from_18_to_18)
summary<-read.xlsx("Summary_all_mut_wt_last_gen.xlsx")
#Identification of significant differences
ggdensity(summary$New_lasting_snps, 
          main = "Density plot of SNPs epimutations_number",
          xlab = "Total number of SNPs epimutations_number")
shapiro.test(summary$New_lasting_snps)
kruskal.test(New_lasting_snps ~ Condition, data = summary)
#--> no difference in SNPS

ggdensity(summary$New_lasting_indels, 
          main = "Density plot of indels epimutations_number",
          xlab = "Total number of indels epimutations_number")
shapiro.test(summary$New_lasting_indels)
kruskal.test(New_lasting_indels ~ Condition, data = summary)
library(conover.test)
conover.test  (summary$New_lasting_indels, g=summary$Condition, method=p.adjustment.methods, kw=TRUE, label=TRUE, 
               wrap=FALSE, table=TRUE, list=FALSE, rmc=FALSE, alpha=0.05, altp=FALSE)

summary$Condition <- fct_relevel(summary$Condition, c("Control", "Low dose", "High dose"))
New_indels_mut_gen_condition <- ggplot(summary, aes(x=Condition, y=New_lasting_indels, group=Condition, color = Condition)) + 
  geom_boxplot(fatten = 1, lwd = 1)+
  scale_color_manual(values=c("cornflowerblue", "darkgreen", "red"))+
  labs(y = "Nb of indels", x = "\nCondition")+
  geom_dotplot(binaxis='y', binwidth = 0.25, stackdir='center', dotsize=1.2)+
  theme_bw()+
  theme_linedraw()+
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(color = "black", size = 24, angle = 0, face = "plain"),
        axis.text.y = element_text(color = "black", size = 24, angle = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 20, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, face = "bold"))+
  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))+  
  ggtitle(paste(""))
#Fig.2.A
New_indels_mut_gen_condition

New_indels_mut_gen_condition_detailed <- ggplot(summary, aes(x=Generation, y=New_lasting_indels, color = Condition, group= Generation)) + 
  geom_boxplot(fatten = 1, lwd = 1)+
  scale_color_manual(values=c("cornflowerblue", "darkgreen", "red"))+
  labs(y = "Nb of indels", x = "\nGeneration")+
  geom_jitter(aes(shape=Lineage, size=8),
              position=position_jitter(width = 0.2,
                                       height = 0.2))+
  theme_bw()+
  theme(axis.text.x = element_text(color = "black", size = 24, angle = 0, face = "plain"),
        axis.text.y = element_text(color = "black", size = 24, angle = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 20, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, face = "bold"))+
  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))+  
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=16), #change legend title font size
        legend.text = element_text(size=14))+
  ggtitle(paste(""))
#Sup.Fig.3
New_indels_mut_gen_condition_detailed

#Get table with raw data used for figures
write.xlsx(summary,"Data_Fig_2_A_&_Sup_Fig_3.xlsx")
#--------------------------
#####Fig.2.B.#####
#Set working directory
setwd("~/Documents/Cisplatin project/Data_analysis/Count data")
#Data prep
#Import raw data RNAseq
RNA_list <- list.files(path = "~/Documents/Cisplatin project/Data_analysis/Count data", pattern = "_Alt_Gen_RNA_")
RNA_counts_total <- c()
for(i in 1:length(RNA_list)){
  Tmp <- read.table(RNA_list[i], sep="\t", stringsAsFactors=F)
  counts <- Tmp[,6]
  RNA_counts_total <- cbind(RNA_counts_total, counts)
}
#Import one sample table to extract genes names
C_elegans_gene_names <- read.xlsx("C_elegans_gene_names2.xlsx")
row.names(RNA_counts_total)<-paste(C_elegans_gene_names[,1], C_elegans_gene_names[,2], C_elegans_gene_names[,3], C_elegans_gene_names[,4],C_elegans_gene_names[,5], sep = ":")
#Give names to columns 
colnames(RNA_counts_total)<-c("C1_0","C2_0","L1_0","L2_0","H1_0","H2_0",
                              "C1_10","C2_10","L1_10","L2_10","H1_10","H2_10",
                              "C1_12","L1_12","L2_12","H1_12","H2_12",
                              "C1_14","C2_14","L1_14","L2_14","H1_14","H2_14",
                              "C1_16","C2_16","L1_16","L2_16","H1_16","H2_16",
                              "C1_18","C2_18","L1_18","L2_18","H1_18","H2_18",
                              "C1_2","C2_2","L1_2","L2_2","H1_2","H2_2",
                              "C1_20","C2_20","L1_20","L2_20","H2_20",
                              "C1_4","C2_4","L1_4","L2_4","H1_4","H2_4",
                              "C1_6","C2_6","L1_6","L2_6","H1_6","H2_6",
                              "C1_8","C2_8","L1_8","L2_8","H1_8","H2_8")
#We remove the 21-URS
Names <- table(C_elegans_gene_names[,4]) 
Names <- Names[-grep("21u", names(Names))]
#Keep only the longest isoform per gene
longest<-c()
for(i in 1:length(Names)){q<-which(C_elegans_gene_names[,4]==names(Names)[i]); L <- C_elegans_gene_names[q,3] - C_elegans_gene_names[q,2]; longest<-c(longest,q[order(L, decreasing=T)][1])}
RNA_counts_total_longest_isoform <- RNA_counts_total[longest,]
# The final step is to select out entries where there were >1 count in at least one column
# This removes lines that do not have substantial reads in at least one sample.
RNA_counts_total2 <- RNA_counts_total_longest_isoform[which(apply(RNA_counts_total_longest_isoform,1,max)>1),]
RNA_countData <- RNA_counts_total2
#Import raw data tinyRNAs
alltiny<-read.xlsx("tinyrna_2022-10-19_10-36-23_feature_counts.xlsx")
colnames(alltiny)<-c("Feature.ID","Tag","Feature.Name","Feature.Class","C1_0","C1_18","C1_20","C1_2","C1_4","C1_6","C1_8","C1_10","C1_12","C1_14","C1_16",
                     "C2_0","C2_20","C2_2","C2_4","C2_8","C2_10","C2_12","C2_14","C2_16","C2_18",
                     "L1_0","L1_20","L1_2","L1_4","L1_6","L1_8","L1_12","L1_14","L1_16","L1_18",
                     "L2_0","L2_20","L2_2","L2_4","L2_8","L2_10","L2_12","L2_14","L2_16","L2_18",
                     "H1_0","H1_20","H1_2","H1_4","H1_8","H1_10","H1_12","H1_14","H1_16","H1_18",
                     "H2_0","H2_18","H2_20","H2_2","H2_4","H2_6","H2_8","H2_10","H2_12","H2_14","H2_16")
alltiny_genes<-read.xlsx("alltiny_genes.xlsx")
alltiny<-cbind(alltiny,alltiny_genes)
#Import Ahringer table
load("/Users/manonfallet/Documents/Cisplatin project/Data_analysis/Count data/Ahringer_single_gene_ref_table.Rdata")
Ahringer<-Ahringer_single_gene_ref_table
write.table(Ahringer, file = "Ahringer.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE)
#Change working directory
setwd("~/Documents/Cisplatin project/Data_analysis/Analysis")
#Normalisation
#Deseq on RNA to compare conditions
###Create data frame with control and low dose 
RNA_counts_totalCvsL<-cbind(RNA_countData[,1], RNA_countData[,36], RNA_countData[,47], RNA_countData[,53], RNA_countData[,59], RNA_countData[,7], RNA_countData[,13], RNA_countData[,18], RNA_countData[,24], RNA_countData[,30], RNA_countData[,42],
                            RNA_countData[,2], RNA_countData[,37],RNA_countData[,48],RNA_countData[,54],RNA_countData[,60],RNA_countData[,8], RNA_countData[,19], RNA_countData[,25],  RNA_countData[,31],  RNA_countData[,43],
                            RNA_countData[,3], RNA_countData[,38],RNA_countData[,49],RNA_countData[,55],RNA_countData[,61],RNA_countData[,9], RNA_countData[,14], RNA_countData[,20],  RNA_countData[,26], RNA_countData[,32],RNA_countData[,44],
                            RNA_countData[,4], RNA_countData[,39],RNA_countData[,50],RNA_countData[,56],RNA_countData[,62],RNA_countData[,10], RNA_countData[,15], RNA_countData[,21],  RNA_countData[,27], RNA_countData[,33],RNA_countData[,45])
colnames(RNA_counts_totalCvsL) <- c("C1_0", "C1_2", "C1_4", "C1_6", "C1_8", "C1_10", "C1_12", "C1_14", "C1_16", "C1_18","C1_20",
                                    "C2_1","C2|_2","C2_4", "C2_6", "C2_8", "C2_10", "C2_14", "C2_16", "C2_18","C2_20",
                                    "L1_0", "L1_2", "L1_4", "L1_6", "L1_8", "L1_10", "L1_12", "L1_14", "L1_16", "L1_18","L1_20",
                                    "L2_0", "L2_2", "L2_4", "L2_6", "L2_8", "L2_10", "L2_12", "L2_14", "L2_16", "L2_18","L2_20")
#Create the list of conditions to compare
conditionCL<-c("Control","Control","Control","Control","Control","Control",
               "Control","Control","Control","Control","Control","Control",
               "Control","Control","Control","Control","Control",
               "Control","Control","Control","Control","Low_dose","Low_dose",
               "Low_dose","Low_dose","Low_dose","Low_dose","Low_dose","Low_dose",
               "Low_dose","Low_dose","Low_dose","Low_dose","Low_dose","Low_dose",
               "Low_dose","Low_dose","Low_dose","Low_dose","Low_dose","Low_dose",
               "Low_dose","Low_dose")
# generate the DESeqDataSet
RNA_DESeq.dsCvsL<-DESeqDataSetFromMatrix(countData = RNA_counts_totalCvsL,
                                         colData=DataFrame(conditionCL),
                                         design=~conditionCL)
ddsCvsL <- DESeq(RNA_DESeq.dsCvsL)
resCvsL <- results(ddsCvsL)
#order results by pval
resOrderedCvsL <- resCvsL[order(resCvsL$pvalue),]
#How many genes are dif exp between low dose and control?
sum(resCvsL$padj < 0.05, na.rm=TRUE)
#get the significant results
resSigCvsL <- subset(resOrderedCvsL, padj < 0.05)
resSigCvsL<-as.data.frame(resSigCvsL)
library(tibble)
resSigCvsL_table <- tibble::rownames_to_column(resSigCvsL, "WB")
write.xlsx(resSigCvsL_table,"resSigCvsL.xlsx")
#Volcano plot
resOrderedCvsL<-as.data.frame(resOrderedCvsL)
resOrderedCvsL <- tibble::rownames_to_column(resOrderedCvsL, "WB")
resOrderedCvsL$Dif.exp <- "NO"
resOrderedCvsL$Dif.exp[resOrderedCvsL$log2FoldChange > 0.5 & resOrderedCvsL$padj < 0.05] <- "UP"
resOrderedCvsL$Dif.exp[resOrderedCvsL$log2FoldChange < 0.5 & resOrderedCvsL$padj < 0.05] <- "DOWN"
resOrderedCvsL$WB <-sapply(strsplit(resOrderedCvsL$WB,":"), `[`, 4)
resOrderedCvsL$delabel <- NA
resOrderedCvsL$delabel[resOrderedCvsL$Dif.exp != "NO"] <- resOrderedCvsL$WB[resOrderedCvsL$Dif.exp != "NO"]
#Final plot-Fig.1.B
library(ggrepel)
ggplot(data=resOrderedCvsL, aes(x=log2FoldChange, y=-log10(padj), col=Dif.exp, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  theme(axis.text.x=element_text(size=30),axis.text.y=element_text(size=30), axis.title=element_text(size=32,face="bold")) +
  geom_text_repel(size=5) +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.5, 0.5), col="red") +
  geom_hline(yintercept=-log10(0.05), col="purple")
#Get raw data table
write.xlsx(resOrderedCvsL,"Data_Fig_2_B.xlsx")
#--------------------------
#####Fig.2.C.#####
#Create data frame with control and high dose 
RNA_counts_totalCvsH<-cbind(RNA_countData[,1], RNA_countData[,36], RNA_countData[,47], RNA_countData[,53], RNA_countData[,59], RNA_countData[,7], RNA_countData[,13], RNA_countData[,18], RNA_countData[,24], RNA_countData[,30], RNA_countData[,42],
                            RNA_countData[,2], RNA_countData[,37],RNA_countData[,48],RNA_countData[,54],RNA_countData[,60],RNA_countData[,8], RNA_countData[,19], RNA_countData[,25],  RNA_countData[,31],  RNA_countData[,43],
                            RNA_countData[,5], RNA_countData[,40],RNA_countData[,51],RNA_countData[,57],RNA_countData[,63],RNA_countData[,11], RNA_countData[,16], RNA_countData[,22],  RNA_countData[,28], RNA_countData[,34],
                            RNA_countData[,6], RNA_countData[,41],RNA_countData[,52],RNA_countData[,58],RNA_countData[,64],RNA_countData[,12], RNA_countData[,17], RNA_countData[,23],  RNA_countData[,29], RNA_countData[,35],RNA_countData[,46])
colnames(RNA_counts_totalCvsH) <- c("C1_0", "C1_2", "C1_4", "C1_6", "C1_8", "C1_10", "C1_12", "C1_14", "C1_16", "C1_18","C1_20",
                                    "C2_1","C2|_2","C2_4", "C2_6", "C2_8", "C2_10", "C2_14", "C2_16", "C2_18","C2_20",
                                    "H1_0", "H1_2", "H1_4", "H1_6", "H1_8", "H1_10", "H1_12", "H1_14", "H1_16", "H1_18",
                                    "H2_0", "H2_2", "H2_4", "H2_6", "H2_8", "H2_10", "H2_12", "H2_14", "H2_16", "H2_18","H2_20")
#Create the list of conditions to compare
conditionCH<-c("Control","Control","Control","Control","Control","Control",
               "Control","Control","Control","Control","Control","Control",
               "Control","Control","Control","Control","Control",
               "Control","Control","Control","Control","High_dose","High_dose",
               "High_dose","High_dose","High_dose","High_dose","High_dose","High_dose",
               "High_dose","High_dose","High_dose","High_dose","High_dose","High_dose",
               "High_dose","High_dose","High_dose","High_dose","High_dose","High_dose",
               "High_dose")
# generate the DESeqDataSet
RNA_DESeq.dsCvsH<-DESeqDataSetFromMatrix(countData = RNA_counts_totalCvsH,
                                         colData=DataFrame(conditionCH),
                                         design=~conditionCH)

ddsCvsH <- DESeq(RNA_DESeq.dsCvsH)
resCvsH <- results(ddsCvsH)
#order results by pval
resOrderedCvsH <- resCvsH[order(resCvsH$pvalue),]
#How many genes are dif exp between low dose and control?
sum(resCvsH$padj < 0.05, na.rm=TRUE)
#get the significant results
resSigCvsH <- subset(resOrderedCvsH, padj < 0.05)
resSigCvsH<-as.data.frame(resSigCvsH)
library(tibble)
resSigCvsH_table <- tibble::rownames_to_column(resSigCvsH, "WB")
write.xlsx(resSigCvsH_table,"resSigCvsH.xlsx")
#Volcano plot
resOrderedCvsH<-as.data.frame(resOrderedCvsH)
resOrderedCvsH <- tibble::rownames_to_column(resOrderedCvsH, "WB")
resOrderedCvsH$Dif.exp <- "NO"
resOrderedCvsH$Dif.exp[resOrderedCvsH$log2FoldChange > 0.5 & resOrderedCvsH$padj < 0.05] <- "UP"
resOrderedCvsH$Dif.exp[resOrderedCvsH$log2FoldChange < -0.5 & resOrderedCvsH$padj < 0.05] <- "DOWN"
resOrderedCvsH$WB <-sapply(strsplit(resOrderedCvsH$WB,":"), `[`, 4)
resOrderedCvsH$delabel <- NA
resOrderedCvsH$delabel[resOrderedCvsH$Dif.exp != "NO"] <- resOrderedCvsH$WB[resOrderedCvsH$Dif.exp != "NO"]

#Final plot-Fig.2.C
library(ggrepel)
ggplot(data=resOrderedCvsH, aes(x=log2FoldChange, y=-log10(padj), col=Dif.exp, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  theme(axis.text.x=element_text(size=30),axis.text.y=element_text(size=30), axis.title=element_text(size=32,face="bold")) +
  geom_text_repel(size=5) +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.5, 0.5), col="red") +
  geom_hline(yintercept=-log10(0.05), col="purple")

#Get the raw data
write.xlsx(resOrderedCvsH,"Data_Fig_2_C.xlsx")
#--------------------------
#####Fig.2.D.#####
#Shared genes with DEGs
DegL<-row.names(resSigCvsL)
DegH<-row.names(resSigCvsH)
venn.diagram(
  x = list(DegL, DegH),
  category.names = c("DegL", "DegH"),
  col = "transparent",
  fill = c("darkgreen","red"),
  alpha = 0.30,
  print.mode=c("raw","percent"),
  filename = "venn_diagramm_DEGL_vs_DEGH.png",
  imagetype="png",
  output=TRUE,
)
#Fisher test
dfvennDEGLvsDEGH <- data.frame("DEGL" = c(5, 1), "DEGL'" = c(1587, 24433), row.names = c("DEGH", "DEGH'"))
fisher.test(dfvennDEGLvsDEGH)
#Get the list of shared genes 
Shared_genes_names_list<-intersect(DegH, DegL)
Shared_genes_L<-resSigCvsL[rownames(resSigCvsL) %in% Shared_genes_names_list,]
Log2FC_shared_genes_L<-as.data.frame(Shared_genes_L[,2])
row.names(Log2FC_shared_genes_L)<-row.names(Shared_genes_L)
colnames(Log2FC_shared_genes_L)<-c("LLog2FC")
Shared_genes_H<-resSigCvsH[rownames(resSigCvsH) %in% Shared_genes_names_list,]
Log2FC_shared_genes_H<-as.data.frame(Shared_genes_H[,2])
row.names(Log2FC_shared_genes_H)<-row.names(Shared_genes_H)
colnames(Log2FC_shared_genes_H)<-c("HLog2FC")
Log2FC<-cbind(Log2FC_shared_genes_L,Log2FC_shared_genes_H)

#--------------------------
#####Sup.Table.4 & Sup.Fig.2.B#####
#EnrichR preparation
library(enrichR)
# tell the enrichR package to use worms
setEnrichrSite("WormEnrichr")
# find gene set databases available
dbs <- listEnrichrDbs()
# We need to tell enrichR which databases (from the selection in dbs) we would like to query.
# We can start with KEGG and the three 2018 GO databases
chosendbs <- c("KEGG_2019",
               "GO_Cellular_Component_2018",
               "GO_Molecular_Function_2018",
               "GO_Biological_Process_2018", 
               "InterPro_Domains_2019")
#High dose gene expression vs all genes
# The test gene list
DegH<-row.names(resSigCvsH)
Sim_test_list <- unique(DegH)
Sim_test_list <-sapply(strsplit(Sim_test_list,":"), `[`, 4)
# The background list
Sim_bg_list <- unique(C_elegans_gene_names[,4])
intersect <- intersect(Sim_bg_list, Sim_test_list)
Sim_bg_list <- Sim_bg_list[Sim_bg_list %in% intersect == F]
# Send both lists to enrichr
Sim_epimutated_enriched_list  <- enrichr(Sim_test_list, chosendbs)
Sim_epimutated_bg_list  <- enrichr(Sim_bg_list, chosendbs)
Sim_test_bg_list <- list(Sim_epimutated_enriched_list, Sim_epimutated_bg_list)
# Annotate with name of libraries
Sim_test_bg_results <- lapply(Sim_test_bg_list, function(x){
  lapply(1:5, function(j){
    if(nrow(x[[j]]) > 0){
      x[[j]][1:nrow(x[[j]]), "library"] <- names(Sim_test_bg_list[[1]])[j]
      return(x[[j]])
    } 
  })
})
# then we combine them and put them in order of significance using the Old.Adjusted.P.value
Sim_test_bg_RESULTS <- lapply(Sim_test_bg_results, function(x){
  tempdf <- do.call(rbind, x)
  tempdf[order(tempdf$Old.Adjusted.P.value), ]
})
split_value_list <- list()
for(i in 1:length(Sim_test_bg_RESULTS)){
  split <- strsplit(Sim_test_bg_RESULTS[[i]][, 2], "/")
  split_value <- c()
  for(j in 1:length(split)){
    save <- split[[j]][1]
    split_value <- as.numeric(c(split_value, save))
  }
  Sim_test_bg_RESULTS[[i]]$Overlap <- split_value
}
write.xlsx(Sim_test_bg_results[[1]],"Sim_epimutated_enriched_list.xlsx")
write.xlsx(Sim_test_bg_results[[2]],"Sim_epimutated_bg_list.xlsx")
# get the test_term, the test_not, the mean(sample_term) and the sample_not
# test_term is overlap
# test_not is length of gene set minus overlap
# then add 1 to all if either test_term OR sample_term = 0
# then only pursue if test_term &/OR sample_term > 4
# let the terms of interest be the terms for the enrichr result for the test list, i.e. with simultaneously epimutated genes
terms_of_interest <- Sim_test_bg_RESULTS[[1]]$Term
# and if the term is not present in the other gene list then assume there is no enrichment for that term in those lists
Test_List <- Sim_test_bg_RESULTS[[1]]
BG_List <- Sim_test_bg_RESULTS[[2]]
table_comp <- c()
for(i in 1:length(terms_of_interest)){
  Term_select <- terms_of_interest[[i]]
  sample_term <- 0
  if(Term_select %in% BG_List$Term == T){
    sample_term <- BG_List[which(BG_List$Term %in% Term_select), 2]}
    sample_not <- length(Sim_bg_list) - sample_term
  test_term <- Test_List[which(Test_List$Term %in% Term_select), 2]
  test_not <- length(Sim_test_list) - test_term
  # save the components in a row in a new table
  component_row <- data.frame(test_term, test_not, sample_term, sample_not)
  rownames(component_row) <- Term_select
  table_comp <- rbind(table_comp, component_row)
  colnames(table_comp) <- c("test_term", "test_not", "sample_term", "sample_not")
}
# Now we will adjust  all values by adding 1 if there are any zeros
zeros <- colSums(table_comp == 0)
if(sum(zeros) > 0){
  table_comp <-  table_comp + 1
}
# Now modify the table to remove any rows where test term or mean sample term are < 5
# Now apply the function only if at least one category has 5 or more genes
# It is only reasonable to assess for relative enrichment/depletion if either the test list or sample list has > 5 genes in that ontology category
table_comp <- table_comp[(table_comp$test_term < 5 & table_comp$sample_term < 5)==F, ]
# Now we have found terms where there are at least 5 genes in 1 test or sample category
# we have adjusted the data for each gene set for any 0 values relative to the specific test list
# Now we can do the analysis but only on terms present in each table_comp list
Pval_col <- c()
OR_col <- c()
for(q in 1:nrow(table_comp)){
  test_term <- table_comp[q, 1] 
  test_not <- table_comp[q, 2]
  mean_Sample_term <- table_comp[q, 3]
  sample_not <- table_comp[q, 4]
  contingency_table <-
    rbind(c(test_term, test_not),
          c(mean_Sample_term, sample_not))
  FT_out <- fisher.test(contingency_table)
  pval <- FT_out$p.value
  OR <- FT_out$estimate
  Pval_col <- c(Pval_col, pval)
  OR_col<- c(OR_col, OR)
}    
OR_pval_frame <- data.frame(Pval_col, OR_col)
colnames(OR_pval_frame) <- c("Pval", "OR")
rownames(OR_pval_frame) <- rownames(table_comp)
# Bonferroni correction
OR_pval_frame[, 1] <- p.adjust(OR_pval_frame[, 1], method="bonferroni", n = nrow(OR_pval_frame)) 
# Make plots for this
# i) In each plot order them by the largest significant enrichment out of all of the comparisons 
OR_pval_frame_ordered <- OR_pval_frame[order(OR_pval_frame[, 2], decreasing = T), ]
# Make a bubble plot 
P_1 <- as.numeric(OR_pval_frame_ordered[, 1])
OR_1 <- log10(as.numeric(OR_pval_frame_ordered[, 2]))
neg_log_p_1 <- -log(P_1)
bubble_table_1 <- data.frame(rownames(OR_pval_frame_ordered), OR_1,  P_1, neg_log_p_1)
colnames(bubble_table_1) <- c("Term", "log_OR", "p_val", "neg_log_P")
bubble_table_1$Term <- factor(bubble_table_1$Term, levels = rev(bubble_table_1$Term))
find_size <- c()
find_alpha <- c()
for(i in 1:nrow(bubble_table_1)){
  alpha <- 0.3
  size_point <- 14
  if(as.numeric(bubble_table_1[i, 3]) < 0.1){
    alpha <- 1
    size_point <- 1*bubble_table_1[i,4]
    }
  find_alpha <- c(find_alpha, alpha)  
  find_size <- c(find_size, size_point)
}
bubble_table <- cbind(bubble_table_1, find_alpha, find_size)
#Sup.Table.4
write.xlsx(bubble_table,"bubble_table_HD_GO.xlsx")
# just plot first 10
allsignifHD_DEGs<-subset(bubble_table,bubble_table$p_val<0.05)
write.xlsx(allsignifHD_DEGs,"allsignifHD_DEGs.xlsx")
bubble_table <- bubble_table[1:10, ]
trunc_terms <- as.character(bubble_table$Term)
# Trunc terms manually entered as Y axis labels in Adobe Illustrator
trunc_terms <- c(
  "nucleus",
  "integral component of plasma membrane",                                               
  "mitochondrion",                
  "cytosol",                                  
  "embryo development ending in birth or egg hatching",                            
  "proteasome-mediated ubiquitin-dependent protein catabolic process",              
  "ubiquitin-dependent protein catabolic process",                                       
  "chemical synaptic transmission",     
  "proteasomal protein catabolic process",                 
  "Spliceosome"
)
guides_merge <- function(gdefs) {
  gdefs <- lapply(gdefs, function(g) { g$hash <- paste(g$order, g$hash, sep = "z"); g})
  tapply(gdefs, sapply(gdefs, function(g)g$hash), function(gs)Reduce(guide_merge, gs))
}
environment(guides_merge) <- environment(ggplot)
assignInNamespace("guides_merge", guides_merge, pos = "package:ggplot2")
bubble_Sim_RNA_chrom <-
  ggplot(bubble_table, aes(y=Term, x=as.numeric(log_OR)))+
  geom_point(aes(color=Term, size=as.numeric(find_size), alpha=find_alpha))+
  scale_x_continuous(breaks = c(1.40,1.45), limits = c(1.40,1.46)) +
  geom_vline(xintercept = 0, colour = "grey")+
  theme_bw()+
  scale_size_continuous(range = c(2, 15), breaks = c(200,100,50), 
                        limits = c(12, 300))+  
  scale_alpha(name = paste("Transparency indicates \nsignificance of enrichment"), range = c(0.3, 1),
              breaks = seq(1, 0.3, length = 2),
              limits = c(0.29, 1.05), 
              labels = c("p value significant < 0.05", "p value not significant > 0.05"))+
  ggtitle(paste("Gene ontology High dose DEGs"))
#Sup.Fig.2.B
bubble_Sim_RNA_chrom <- bubble_Sim_RNA_chrom + labs(y="Gene Ontology Terms\n", x = "log10(Odds Ratio for enrichment)")+
  theme(axis.title.x=element_text(face = "bold", size=14))+
  theme(axis.title.y=element_text(face = "bold", size=14))+
  theme(plot.title = element_text(face = "bold", size=16))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_discrete(labels= rev(trunc_terms))+
  # scale_y_discrete(labels= "")+
  theme(axis.text.x = element_text(color="#000000", size=10))+
  theme(axis.text.y = element_text(color="#000000", size = 17))+
  guides(size = guide_legend(order = 3), colour = "none", alpha = guide_legend(order = 2))+
  theme(legend.text=element_text(color="#000000", size=12))
bubble_Sim_RNA_chrom

#Get the raw data
write.xlsx(bubble_table,"Data_Sup_Fig_2_B.xlsx")
#--------------------------
#####Table 1 - Association between DNA mutations & gene exp. changes#####
#RNA normalization
sample_info<-data.frame(condition=colnames(RNA_countData),row.names=colnames(RNA_countData))
RNA_DESeq.dsall<-DESeqDataSetFromMatrix(countData = RNA_countData,
                                        colData=sample_info,
                                        design=~condition)

# obtain regularized log - transformed values
vstall<-vst(RNA_DESeq.dsall,blind = TRUE)
RNA_vst.norm.count.all<-assay(vstall)

#calculate the size factor and add it to the data set
RNA_DESeq.dsall<-estimateSizeFactors(RNA_DESeq.dsall)
sizeFactors(RNA_DESeq.dsall)
colData(RNA_DESeq.dsall)

# retrieve the _ normalized _ read counts
RNA_countData_all<-counts(RNA_DESeq.dsall,normalized=TRUE)
save(RNA_countData_all,file="RNA_countData_all.Rdata")

# transform size - factor normalized read counts to log2 scale using a pseudocount of 1
RNA_log.norm.counts<-log2(RNA_countData_all+1)
boxplot(RNA_log.norm.counts,notch=TRUE,
        main="log2-transformed read counts",
        ylab="log2(read counts)")

#Detection of epimutations in each lineage
#Create data frame per lineage
RNA_Lineage_C1 <- cbind(RNA_log.norm.counts[,1], RNA_log.norm.counts[,36], RNA_log.norm.counts[,47], RNA_log.norm.counts[,53], RNA_log.norm.counts[,59], RNA_log.norm.counts[,7], RNA_log.norm.counts[,13],RNA_log.norm.counts[,18],RNA_log.norm.counts[,24],RNA_log.norm.counts[,30],RNA_log.norm.counts[,42])
colnames(RNA_Lineage_C1) <- c("C1_0", "C1_2", "C1_4", "C1_6", "C1_8", "C1_10", "C1_12", "C1_14", "C1_16", "C1_18","C1_20")
RNA_Lineage_C2 <- cbind(RNA_log.norm.counts[,2], RNA_log.norm.counts[,37],RNA_log.norm.counts[,48],RNA_log.norm.counts[,54],RNA_log.norm.counts[,60],RNA_log.norm.counts[,8], RNA_log.norm.counts[,19],RNA_log.norm.counts[,25],RNA_log.norm.counts[,31], RNA_log.norm.counts[,43])
colnames(RNA_Lineage_C2) <- c("C2_0", "C2_2", "C2_4", "C2_6", "C2_8", "C2_10","C2_14", "C2_16", "C2_18","C2_20")
RNA_Lineage_L1 <- cbind(RNA_log.norm.counts[,3], RNA_log.norm.counts[,38],RNA_log.norm.counts[,49],RNA_log.norm.counts[,55],RNA_log.norm.counts[,61],RNA_log.norm.counts[,9], RNA_log.norm.counts[,14],RNA_log.norm.counts[,20],RNA_log.norm.counts[,26],RNA_log.norm.counts[,32],RNA_log.norm.counts[,44])
colnames(RNA_Lineage_L1) <- c("L1_0", "L1_2", "L1_4", "L1_6", "L1_8", "L1_10", "L1_12", "L1_14", "L1_16", "L1_18","L1_20")
RNA_Lineage_L2 <- cbind(RNA_log.norm.counts[,4], RNA_log.norm.counts[,39],RNA_log.norm.counts[,50],RNA_log.norm.counts[,56],RNA_log.norm.counts[,62],RNA_log.norm.counts[,10], RNA_log.norm.counts[,15],RNA_log.norm.counts[,21], RNA_log.norm.counts[,27],RNA_log.norm.counts[,33],RNA_log.norm.counts[,45])
colnames(RNA_Lineage_L2) <- c("L2_0", "L2_2", "L2_4", "L2_6", "L2_8", "L2_10", "L2_12", "L2_14", "L2_16", "L2_18","L2_20")
RNA_Lineage_H1 <- cbind(RNA_log.norm.counts[,5], RNA_log.norm.counts[,40],RNA_log.norm.counts[,51],RNA_log.norm.counts[,57],RNA_log.norm.counts[,63],RNA_log.norm.counts[,11], RNA_log.norm.counts[,16], RNA_log.norm.counts[,22],RNA_log.norm.counts[,28],RNA_log.norm.counts[,34])
colnames(RNA_Lineage_H1) <- c("H1_0", "H1_2", "H1_4", "H1_6", "H1_8", "H1_10", "H1_12", "H1_14", "H1_16", "H1_18")
RNA_Lineage_H2 <- cbind(RNA_log.norm.counts[,6], RNA_log.norm.counts[,41],RNA_log.norm.counts[,52],RNA_log.norm.counts[,58],RNA_log.norm.counts[,64],RNA_log.norm.counts[,12], RNA_log.norm.counts[,17], RNA_log.norm.counts[,23],  RNA_log.norm.counts[,29], RNA_log.norm.counts[,35],RNA_log.norm.counts[,46])
colnames(RNA_Lineage_H2) <- c("H2_0", "H2_2", "H2_4", "H2_6", "H2_8", "H2_10", "H2_12", "H2_14", "H2_16", "H2_18","H2_20")
#C1
#Identify epimutations using linear model in each lineage 
RNA_Lineage_C1_2_to_20<-RNA_Lineage_C1[,2:11]
RNA_z_scores_table_C1<-matrix(0,ncol=10, nrow = nrow(RNA_Lineage_C1))

for (i in 1:ncol(RNA_Lineage_C1_2_to_20)) {
  EpimutC1<-
    loess(RNA_Lineage_C1[,1]~RNA_Lineage_C1_2_to_20[,i])
  ResidC1<-EpimutC1$residuals
  RNA_z_scores_table_C1[,i]<-(ResidC1-mean(ResidC1))/sd(ResidC1)
}

#keep only relevant epimutations and binarised data
RNA_binarised_z_scores_table_C1<-matrix(0,ncol=ncol(RNA_z_scores_table_C1),nrow=nrow(RNA_z_scores_table_C1))

for (i in 1:ncol(RNA_z_scores_table_C1)) {
  for(j in 1:nrow(RNA_z_scores_table_C1)){
    
    if(RNA_z_scores_table_C1[j,i]> (-2.25) & RNA_z_scores_table_C1[j,i]< 2.25) {
      RNA_binarised_z_scores_table_C1[j,i]<-0
    }
    if(RNA_z_scores_table_C1[j,i] > 2.25){
      RNA_binarised_z_scores_table_C1[j,i]<-1
    }
    if(RNA_z_scores_table_C1[j,i] < (-2.25)){
      RNA_binarised_z_scores_table_C1[j,i]<- (-1)
    }
  }  
}
RNA_C1epimutations<-RNA_binarised_z_scores_table_C1 %>% as.data.frame() 
colnames(RNA_C1epimutations) <- c('F2','4','6','8','10','12','14','16','18','20')
RNA_C1epimutations <- RNA_C1epimutations %>%
  # Creating an empty column:
  add_column(F0 = 0, .before='F2')
RNA_C1epimutations <- RNA_C1epimutations %>%
  # Creating a condition column:
  add_column(condition = "C1", .before="F0")
RNA_C1epimutations <- RNA_C1epimutations %>%
  # Creating a genes column:
  add_column(genes = row.names(RNA_Lineage_C1), .before="condition")
colnames(RNA_C1epimutations) <- c('genes','condition','0','2','4','6','8','10','12','14','16','18','20')
save(RNA_C1epimutations,file="RNA_C1epimutations.Rdata")
#C2
RNA_Lineage_C2_2_to_20<-RNA_Lineage_C2[,2:10]
RNA_z_scores_table_C2<-matrix(0,ncol=9, nrow = nrow(RNA_Lineage_C2))

for (i in 1:ncol(RNA_Lineage_C2_2_to_20)) {
  EpimutC2<-
    loess(RNA_Lineage_C2[,1]~RNA_Lineage_C2_2_to_20[,i])
  ResidC2<-EpimutC2$residuals
  RNA_z_scores_table_C2[,i]<-(ResidC2-mean(ResidC2))/sd(ResidC2)
}

#keep only relevant epimutations and binarised data
RNA_binarised_z_scores_table_C2<-matrix(0,ncol=ncol(RNA_z_scores_table_C2),nrow=nrow(RNA_z_scores_table_C2))

for (i in 1:ncol(RNA_z_scores_table_C2)) {
  for(j in 1:nrow(RNA_z_scores_table_C2)){
    
    if(RNA_z_scores_table_C2[j,i]> (-2.25) & RNA_z_scores_table_C2[j,i]< 2.25) {
      RNA_binarised_z_scores_table_C2[j,i]<-0
    }
    if(RNA_z_scores_table_C2[j,i] > 2.25){
      RNA_binarised_z_scores_table_C2[j,i]<-1
    }
    if(RNA_z_scores_table_C2[j,i] < (-2.25)){
      RNA_binarised_z_scores_table_C2[j,i]<- (-1)
    }
  }  
}
RNA_C2epimutations<-RNA_binarised_z_scores_table_C2 %>% as.data.frame() 
colnames(RNA_C2epimutations) <- c('F2','4','6','8','10','14','16','18','20')
RNA_C2epimutations <- RNA_C2epimutations %>%
  # Creating an empty column:
  add_column(F0 = 0, .before='F2')
RNA_C2epimutations <- RNA_C2epimutations %>%
  # Creating a condition column:
  add_column(condition = "C2", .before="F0")
RNA_C2epimutations <- RNA_C2epimutations %>%
  # Creating a genes column:
  add_column(genes = row.names(RNA_Lineage_C2), .before="condition")
RNA_C2epimutations <- RNA_C2epimutations %>%
  # Creating a NA column:
  add_column(F12 = "NA", .before="14")
colnames(RNA_C2epimutations) <- c('genes','condition','0','2','4','6','8','10','12','14','16','18','20')
save(RNA_C2epimutations,file="RNA_C2epimutations.Rdata")
#low1
RNA_Lineage_L1_2_to_20<-RNA_Lineage_L1[,2:11]
RNA_z_scores_table_L1<-matrix(0,ncol=10, nrow = nrow(RNA_Lineage_L1))

for (i in 1:ncol(RNA_Lineage_L1_2_to_20)) {
  EpimutL1<-
    loess(RNA_Lineage_L1[,1]~RNA_Lineage_L1_2_to_20[,i])
  ResidL1<-EpimutL1$residuals
  RNA_z_scores_table_L1[,i]<-(ResidL1-mean(ResidL1))/sd(ResidL1)
}

#keep only relevant epimutations and binarised data
RNA_binarised_z_scores_table_L1<-matrix(0,ncol=ncol(RNA_z_scores_table_L1),nrow=nrow(RNA_z_scores_table_L1))

for (i in 1:ncol(RNA_z_scores_table_L1)) {
  for(j in 1:nrow(RNA_z_scores_table_L1)){
    
    if(RNA_z_scores_table_L1[j,i]> (-2.25) & RNA_z_scores_table_L1[j,i]< 2.25) {
      RNA_binarised_z_scores_table_L1[j,i]<-0
    }
    if(RNA_z_scores_table_L1[j,i] > 2.25){
      RNA_binarised_z_scores_table_L1[j,i]<-1
    }
    if(RNA_z_scores_table_L1[j,i] < (-2.25)){
      RNA_binarised_z_scores_table_L1[j,i]<- (-1)
    }
  }  
}

RNA_L1epimutations <- as.data.frame(RNA_binarised_z_scores_table_L1)
colnames(RNA_L1epimutations) <- c('F2','4','6','8','10','12','14','16','18','20')
RNA_L1epimutations <- RNA_L1epimutations %>%
  # Creating an empty column:
  add_column(F0 = 0, .before="F2")
RNA_L1epimutations <- RNA_L1epimutations %>%
  # Creating an empty column:
  add_column(condition = 'L1', .before="F0")
RNA_L1epimutations <- RNA_L1epimutations %>%
  # Creating a genes column:
  add_column(genes = row.names(RNA_Lineage_L1), .before="condition")
colnames(RNA_L1epimutations) <- c('genes','condition','0','2','4','6','8','10','12','14','16','18','20')
save(RNA_L1epimutations,file="RNA_L1epimutations.Rdata")
#low2
RNA_Lineage_L2_2_to_20<-RNA_Lineage_L2[,2:11]
RNA_z_scores_table_L2<-matrix(0,ncol=10, nrow = nrow(RNA_Lineage_L2))

for (i in 1:ncol(RNA_Lineage_L2_2_to_20)) {
  EpimutL2<-
    loess(RNA_Lineage_L2[,1]~RNA_Lineage_L2_2_to_20[,i])
  ResidL2<-EpimutL2$residuals
  RNA_z_scores_table_L2[,i]<-(ResidL2-mean(ResidL2))/sd(ResidL2)
}

#keep only relevant epimutations and binarised data
RNA_binarised_z_scores_table_L2<-matrix(0,ncol=ncol(RNA_z_scores_table_L2),nrow=nrow(RNA_z_scores_table_L2))

for (i in 1:ncol(RNA_z_scores_table_L2)) {
  for(j in 1:nrow(RNA_z_scores_table_L2)){
    
    if(RNA_z_scores_table_L2[j,i]> (-2.25) & RNA_z_scores_table_L2[j,i]< 2.25) {
      RNA_binarised_z_scores_table_L2[j,i]<-0
    }
    if(RNA_z_scores_table_L2[j,i] > 2.25){
      RNA_binarised_z_scores_table_L2[j,i]<-1
    }
    if(RNA_z_scores_table_L2[j,i] < (-2.25)){
      RNA_binarised_z_scores_table_L2[j,i]<- (-1)
    }
  }  
}

RNA_L2epimutations <- as.data.frame(RNA_binarised_z_scores_table_L2)
colnames(RNA_L2epimutations) <- c('F2','4','6','8','10','12','14','16','18','20')
RNA_L2epimutations <- RNA_L2epimutations %>%
  # Creating an empty column:
  add_column(F0 = 0, .before="F2")
RNA_L2epimutations <- RNA_L2epimutations %>%
  # Creating an empty column:
  add_column(condition = 'L2', .before="F0")
RNA_L2epimutations <- RNA_L2epimutations %>%
  # Creating a genes column:
  add_column(genes = row.names(RNA_Lineage_L2), .before="condition")
colnames(RNA_L2epimutations) <- c('genes','condition','0','2','4','6','8','10','12','14','16','18','20')
save(RNA_L2epimutations,file="RNA_L2epimutations.Rdata")

table(RNA_L2epimutations$`2`)
table(RNA_L2epimutations$`4`)
table(RNA_L2epimutations$`6`)
table(RNA_L2epimutations$`8`)
table(RNA_L2epimutations$`10`)
table(RNA_L2epimutations$`12`)
table(RNA_L2epimutations$`14`)
table(RNA_L2epimutations$`16`)
table(RNA_L2epimutations$`18`)
table(RNA_L2epimutations$`20`)

#high1
#Identify epimutations using linear model in each lineage 
RNA_Lineage_H1_2_to_20<-RNA_Lineage_H1[,2:10]

RNA_z_scores_table_H1<-matrix(0,ncol=9, nrow = nrow(RNA_Lineage_H1))

for (i in 1:ncol(RNA_Lineage_H1_2_to_20)){
  EpimutH1<-
    loess(RNA_Lineage_H1[,1]~RNA_Lineage_H1_2_to_20[,i])
  ResidH1<-EpimutH1$residuals
  RNA_z_scores_table_H1[,i]<-(ResidH1-mean(ResidH1))/sd(ResidH1)
}

RNA_binarised_z_scores_table_H1<-matrix(0,ncol=9, nrow = nrow(RNA_z_scores_table_H1))
for (i in 1:ncol(RNA_z_scores_table_H1)){
  for (j in 1:nrow(RNA_z_scores_table_H1)){
    
    if(RNA_z_scores_table_H1[j,i]> (-2.25) & RNA_z_scores_table_H1[j,i]< 2.25) {
      RNA_binarised_z_scores_table_H1[j,i]<-0
    }
    if(RNA_z_scores_table_H1[j,i] > 2.25){
      RNA_binarised_z_scores_table_H1[j,i]<-1
    }
    if(RNA_z_scores_table_H1[j,i] < (-2.25)){
      RNA_binarised_z_scores_table_H1[j,i]<- (-1)
    }
  }
}
RNA_H1epimutations<-RNA_binarised_z_scores_table_H1 %>% as.data.frame() 
colnames(RNA_H1epimutations) <- c('F2','4','6','8','10','12','14','16','18')
RNA_H1epimutations <- RNA_H1epimutations %>%
  # Creating an empty column:
  add_column(F0 = 0, .before='F2')
RNA_H1epimutations <- RNA_H1epimutations %>%
  # Creating a condition column:
  add_column(condition = "H1", .before="F0")
RNA_H1epimutations <- RNA_H1epimutations %>%
  # Creating a genes column:
  add_column(genes = row.names(RNA_Lineage_H1), .before="condition")
RNA_H1epimutations <- RNA_H1epimutations %>%
  # Creating a NA column for the missing generation:
  add_column(F20 = "NA", .after="18")
colnames(RNA_H1epimutations) <- c('genes','condition','0','2','4','6','8','10','12','14','16','18','20')
save(RNA_H1epimutations,file="RNA_H1epimutations.Rdata")
#high2
RNA_Lineage_H2_2_to_20<-RNA_Lineage_H2[,2:11]
RNA_z_scores_table_H2<-matrix(0,ncol=10, nrow = nrow(RNA_Lineage_H2))

for (i in 1:ncol(RNA_Lineage_H2_2_to_20)) {
  EpimutH2<-loess(RNA_Lineage_H2[,1]~RNA_Lineage_H2_2_to_20[,i])
  ResidH2<-EpimutH2$residuals
  RNA_z_scores_table_H2[,i]<-(ResidH2-mean(ResidH2))/sd(ResidH2)
}

#keep only relevant epimutations and binarised data
RNA_binarised_z_scores_table_H2<-matrix(0,ncol=ncol(RNA_z_scores_table_H2),nrow=nrow(RNA_z_scores_table_H2))

for (i in 1:ncol(RNA_z_scores_table_H2)) {
  for(j in 1:nrow(RNA_z_scores_table_H2)){
    
    if(RNA_z_scores_table_H2[j,i]> (-2.25) & RNA_z_scores_table_H2[j,i]< 2.25) {
      RNA_binarised_z_scores_table_H2[j,i]<-0
    }
    if(RNA_z_scores_table_H2[j,i] > 2.25){
      RNA_binarised_z_scores_table_H2[j,i]<-1
    }
    if(RNA_z_scores_table_H2[j,i] < (-2.25)){
      RNA_binarised_z_scores_table_H2[j,i]<- (-1)
    }
  }  
}

RNA_H2epimutations <- as.data.frame(RNA_binarised_z_scores_table_H2)
colnames(RNA_H2epimutations) <- c('F2','4','6','8','10','12','14','16','18','20')
RNA_H2epimutations <- RNA_H2epimutations %>%
  # Creating an empty column:
  add_column(F0 = 0, .before="F2")
RNA_H2epimutations <- RNA_H2epimutations %>%
  # Creating an empty column:
  add_column(condition = 'H2', .before="F0")
RNA_H2epimutations <- RNA_H2epimutations %>%
  # Creating a genes column:
  add_column(genes = row.names(RNA_Lineage_H2), .before="condition")
colnames(RNA_H2epimutations) <- c('genes','condition','0','2','4','6','8','10','12','14','16','18','20')
save(RNA_H2epimutations,file="RNA_H2epimutations.Rdata")
#Integration RNA/DNA mutations
#An RNA centric table with RNA genes, RNA expression change and DNA mutations
#Data prep
colnames(RNA_C2epimutations) <- c('genes','condition','0','2','4','6','8','10','F12','14','16','18','20')
RNA_C2epimutationsbis <- subset(RNA_C2epimutations,  select = -F12)
RNA_C2epimutationsbis<-cbind(RNA_C2epimutationsbis, rep(RNA_C2epimutationsbis[9],1))
RNA_C2epimutationsbis<-RNA_C2epimutationsbis[,c(1,2,3,4,5,6,7,8,13,9,10,11,12)]
colnames(RNA_C2epimutationsbis) <- c('genes','condition','0','2','4','6','8','10','12','14','16','18','20')
colnames(RNA_H1epimutations) <- c('genes','condition','0','2','4','6','8','10','12','14','16','18','F20')
RNA_H1epimutationsbis <- subset(RNA_H1epimutations,  select = -F20)

rownames(RNA_C1epimutations)<-RNA_C1epimutations$genes
RNA_C1epimutations <- subset(RNA_C1epimutations,  select = -genes)
RNA_C1epimutations <- subset(RNA_C1epimutations,  select = -condition)

rownames(RNA_C2epimutationsbis)<-RNA_C2epimutationsbis$genes
RNA_C2epimutationsbis <- subset(RNA_C2epimutationsbis,  select = -genes)
RNA_C2epimutationsbis <- subset(RNA_C2epimutationsbis,  select = -condition)

rownames(RNA_L1epimutations)<-RNA_L1epimutations$genes
RNA_L1epimutations <- subset(RNA_L1epimutations,  select = -genes)
RNA_L1epimutations <- subset(RNA_L1epimutations,  select = -condition)

rownames(RNA_L2epimutations)<-RNA_L2epimutations$genes
RNA_L2epimutations <- subset(RNA_L2epimutations,  select = -genes)
RNA_L2epimutations <- subset(RNA_L2epimutations,  select = -condition)

rownames(RNA_H1epimutationsbis)<-RNA_H1epimutationsbis$genes
RNA_H1epimutationsbis <- subset(RNA_H1epimutationsbis,  select = -genes)
RNA_H1epimutationsbis <- subset(RNA_H1epimutationsbis,  select = -condition)

rownames(RNA_H2epimutations)<-RNA_H2epimutations$genes
RNA_H2epimutations <- subset(RNA_H2epimutations,  select = -genes)
RNA_H2epimutations <- subset(RNA_H2epimutations,  select = -condition)

# We will need piRNA cluster genes for this

# Coordinates for piRNA cluster genes (5 - 7.5 & 13 - 17 MB from Ruby et al. 2006) have been intersected with RNA coordinates 

RNA_piRNA_domains <- read.table("RNA_Coords_intersected_piRNA_c.bed", header = FALSE)

piRNA_domains <- RNA_piRNA_domains[RNA_piRNA_domains$V5==1, ]

# When we are considering the genes with chromatin domain annotations we can only consider the genes which map to Ahringer as this is where we get the domain annotations from. 

all_Ahr_RNA_genes <- unique(Ahringer_single_gene_ref_table$Gene)

piRNA_cluster_genes <- intersect(unique(piRNA_domains$V4), all_Ahr_RNA_genes)

#Control
Control_RNA_list <- list(RNA_C1epimutations, RNA_C2epimutationsbis)
DNA_mut_C1<-read.xlsx("DNA_mut_C1_in_genes.xlsx")
DNA_mut_C2<-read.xlsx("DNA_mut_C2_in_genes.xlsx")
row.names(DNA_mut_C1)<-DNA_mut_C1$WB
DNA_mut_C1<-DNA_mut_C1[,2:12]
row.names(DNA_mut_C2)<-DNA_mut_C2$WB
DNA_mut_C2<-DNA_mut_C2[,2:12]
Control_DNA_mut_genes <- list(DNA_mut_C1, DNA_mut_C2)

lin <- c("C1", "C2")
Control_RNA_and_DNA_mut_integrated_table <- c()

for(x in 1:length(Control_RNA_list)){
  RNA_bin <- Control_RNA_list[[x]]
  DNA_mut_bin <- Control_DNA_mut_genes[[x]]
  Lineage <- lin[[x]]
  integrated_table <- c()
  ep_table <- c()
  
  for(e in 1:length(rownames(RNA_bin))){
    coord <- rownames(RNA_bin)[e]
    gene <- strsplit(coord, ":")[[1]][4]
    RNA_UP <- 0
    RNA_DOWN <- 0
    RNA_MixedUPDOWN <- 0
    
    RNA_mut <- 0
    
    if(sum(abs(RNA_bin[coord, ])) >0){
      RNA_mut <- 1
      # Direction of RNA expression changes  
      # Determine the original direction of the RNA events    
      events <- sum(abs(RNA_bin[coord, ]))
      sum <- sum(RNA_bin[coord, ])
      
      if(sum == events){
        RNA_UP <- 1
      }     
      
      if(sum == -1*(events)){
        RNA_DOWN <- 1
      }     
      
      if(!abs(sum) == events){
        RNA_MixedUPDOWN <- 1
      }  
    }
    RNA_gens <- "0"
    RNA_inherited <- 0
    
    if(RNA_mut ==1){
      RNA_gens <-colnames(RNA_bin[(which(abs(RNA_bin[coord,])>0))])
      
      if(length(RNA_gens) > 1){
        RNA_gens <- paste(RNA_gens, collapse = "_")
        
        if(2 %in% diff(as.numeric(strsplit(RNA_gens, "_")[[1]]))){
          RNA_inherited <- 1
        } 
      }
    }
    # The RNA coding locus has to a) be present in the small RNA genes so as to get a small RNA signal data point
    # b) be present in the Ahringer data set so as to get associated annotations on chromatin domain distribution
    # First assume the coding locus with RNAseq data is not present in small RNA 
    maps_to_Ahr <- 0
    maps_to_DNA_mut <- 0
    DNA_mut <- 0
    DNA_mut_inherited <- 0
    DNA_mut_gens <- 0
    multiple_DNA_mut_gens  <- 0
    
    if(gene %in% Control_DNA_mut_genes == T){
      maps_to_DNA_mut <- 1
      
      if(sum(abs(DNA_mut_bin[gene, ]))>0){
        DNA_mut <- 1
        target_row <- DNA_mut_bin[gene, ]
        events <- sum(abs(DNA_mut_bin[gene, ]))
      }
      if(DNA_mut ==1){
        # Generations
        DNA_mut_gens <- colnames(DNA_mut_bin[(which(abs(DNA_mut_bin[gene, ]) > 0))])
        
        if(length(DNA_mut_gens)> 1){
          DNA_mut_gens <- paste(DNA_mut_gens, collapse = "_")
          multiple_DNA_mut_gens <- 1
        }
        # Inherited
        inherited <- c()
        inherit <- 0 
        each_row <- colnames(DNA_mut_bin[(which(abs(DNA_mut_bin[gene, ]) > 0))])
        
        if(2 %in% diff(as.numeric(each_row))){
          inherit <- 1}
        inherited <- c(inherited, inherit) 
        
        if(sum(inherited > 0)){
          DNA_mut_inherited <- 1
        }
      }}
    time_matched_to_RNA <- 0
    # Determine time matching of small RNA to RNA
    
    if(RNA_mut==1&DNA_mut==1){
      DNA_mut_gens <- unlist(str_split(DNA_mut_gens, "_"))
      rna_gens <- unlist(str_split(RNA_gens, "_"))
      matched_gens <- DNA_mut_gens[DNA_mut_gens %in% rna_gens]
      matched <-  length(which(DNA_mut_gens %in% rna_gens))
      
      if(matched > 0){
        time_matched_to_RNA <- 1
      }
    }
    save <- data.frame(Lineage, gene, coord,RNA_gens, maps_to_Ahr, maps_to_DNA_mut, RNA_mut, RNA_inherited, 
                       DNA_mut, 
                       DNA_mut_inherited, time_matched_to_RNA,DNA_mut_gens,
                       RNA_UP, RNA_DOWN, RNA_MixedUPDOWN 
    )
    integrated_table <- rbind(integrated_table, save)
  }
  Control_RNA_and_DNA_mut_integrated_table <- rbind(Control_RNA_and_DNA_mut_integrated_table, integrated_table)
}
save(Control_RNA_and_DNA_mut_integrated_table,file="Control_RNA_and_DNA_mut_integrated_table.Rdata")
write.xlsx(Control_RNA_and_DNA_mut_integrated_table,"Control_RNA_and_DNA_mut_integrated_table.xlsx")
#Low dose
Low_dose_RNA_list <- list(RNA_L1epimutations, RNA_L2epimutations)
DNA_mut_L1<-read.xlsx("DNA_mut_L1_in_genes.xlsx")
DNA_mut_L2<-read.xlsx("DNA_mut_L2_in_genes.xlsx")
row.names(DNA_mut_L1)<-DNA_mut_L1$WB
DNA_mut_L1<-DNA_mut_L1[,2:12]
row.names(DNA_mut_L2)<-DNA_mut_L2$WB
DNA_mut_L2<-DNA_mut_L2[,2:12]
Low_dose_DNA_mut_genes <- list(DNA_mut_L1, DNA_mut_L2)

lin <- c("L1", "L2")
Low_dose_RNA_and_DNA_mut_integrated_table <- c()

for(x in 1:length(Low_dose_RNA_list)){
  RNA_bin <- Low_dose_RNA_list[[x]]
  DNA_mut_bin <- Low_dose_DNA_mut_genes[[x]]
  Lineage <- lin[[x]]
  integrated_table <- c()
  ep_table <- c()
  
  for(e in 1:length(rownames(RNA_bin))){
    coord <- rownames(RNA_bin)[e]
    gene <- strsplit(coord, ":")[[1]][4]
    RNA_UP <- 0
    RNA_DOWN <- 0
    RNA_MixedUPDOWN <- 0
    
    RNA_mut <- 0
    
    if(sum(abs(RNA_bin[coord, ])) >0){
      RNA_mut <- 1
      # Direction of RNA expression changes  
      # Determine the original direction of the RNA events    
      events <- sum(abs(RNA_bin[coord, ]))
      sum <- sum(RNA_bin[coord, ])
      
      if(sum == events){
        RNA_UP <- 1
      }     
      
      if(sum == -1*(events)){
        RNA_DOWN <- 1
      }     
      
      if(!abs(sum) == events){
        RNA_MixedUPDOWN <- 1
      }  
    }
    RNA_gens <- "0"
    RNA_inherited <- 0
    
    if(RNA_mut ==1){
      RNA_gens <-colnames(RNA_bin[(which(abs(RNA_bin[coord,])>0))])
      
      if(length(RNA_gens) > 1){
        RNA_gens <- paste(RNA_gens, collapse = "_")
        
        if(2 %in% diff(as.numeric(strsplit(RNA_gens, "_")[[1]]))){
          RNA_inherited <- 1
        } 
      }
    }
    # The RNA coding locus has to a) be present in the small RNA genes so as to get a small RNA signal data point
    # b) be present in the Ahringer data set so as to get associated annotations on chromatin domain distribution
    # First assume the coding locus with RNAseq data is not present in small RNA 
    maps_to_Ahr <- 0
    maps_to_DNA_mut <- 0
    DNA_mut <- 0
    DNA_mut_inherited <- 0
    DNA_mut_gens <- 0
    multiple_DNA_mut_gens  <- 0
    
    if(gene %in% Low_dose_DNA_mut_genes == T){
      maps_to_DNA_mut <- 1
      
      if(sum(abs(DNA_mut_bin[gene, ]))>0){
        DNA_mut <- 1
        target_row <- DNA_mut_bin[gene, ]
        events <- sum(abs(DNA_mut_bin[gene, ]))
      }
      if(DNA_mut ==1){
        # Generations
        DNA_mut_gens <- colnames(DNA_mut_bin[(which(abs(DNA_mut_bin[gene, ]) > 0))])
        
        if(length(DNA_mut_gens)> 1){
          DNA_mut_gens <- paste(DNA_mut_gens, collapse = "_")
          multiple_DNA_mut_gens <- 1
        }
        # Inherited
        inherited <- c()
        inherit <- 0 
        each_row <- colnames(DNA_mut_bin[(which(abs(DNA_mut_bin[gene, ]) > 0))])
        
        if(2 %in% diff(as.numeric(each_row))){
          inherit <- 1}
        inherited <- c(inherited, inherit) 
        
        if(sum(inherited > 0)){
          DNA_mut_inherited <- 1
        }
      }}
    time_matched_to_RNA <- 0
    # Determine time matching of small RNA to RNA
    
    if(RNA_mut==1&DNA_mut==1){
      DNA_mut_gens <- unlist(str_split(DNA_mut_gens, "_"))
      rna_gens <- unlist(str_split(RNA_gens, "_"))
      matched_gens <- DNA_mut_gens[DNA_mut_gens %in% rna_gens]
      matched <-  length(which(DNA_mut_gens %in% rna_gens))
      
      if(matched > 0){
        time_matched_to_RNA <- 1
      }
    }
    save <- data.frame(Lineage, gene, coord,RNA_gens, maps_to_Ahr, maps_to_DNA_mut, RNA_mut, RNA_inherited, 
                       DNA_mut, 
                       DNA_mut_inherited, time_matched_to_RNA,DNA_mut_gens,
                       RNA_UP, RNA_DOWN, RNA_MixedUPDOWN 
    )
    integrated_table <- rbind(integrated_table, save)
  }
  Low_dose_RNA_and_DNA_mut_integrated_table <- rbind(Low_dose_RNA_and_DNA_mut_integrated_table, integrated_table)
}
save(Low_dose_RNA_and_DNA_mut_integrated_table,file="Low_dose_RNA_and_DNA_mut_integrated_table.Rdata")
write.xlsx(Low_dose_RNA_and_DNA_mut_integrated_table,"Low_dose_RNA_and_DNA_mut_integrated_table.xlsx")
#High dose
High_dose_RNA_list <- list(RNA_H1epimutationsbis, RNA_H2epimutations)
DNA_mut_H1<-read.xlsx("DNA_mut_H1_in_genes.xlsx")
DNA_mut_H2<-read.xlsx("DNA_mut_H2_in_genes.xlsx")
row.names(DNA_mut_H1)<-DNA_mut_H1$WB
DNA_mut_H1<-DNA_mut_H1[,2:12]
row.names(DNA_mut_H2)<-DNA_mut_H2$WB
DNA_mut_H2<-DNA_mut_H2[,2:12]
High_dose_DNA_mut_genes <- list(DNA_mut_H1, DNA_mut_H2)

lin <- c("H1", "H2")
High_dose_RNA_and_DNA_mut_integrated_table <- c()

for(x in 1:length(High_dose_RNA_list)){
  RNA_bin <- High_dose_RNA_list[[x]]
  DNA_mut_bin <- High_dose_DNA_mut_genes[[x]]
  Lineage <- lin[[x]]
  integrated_table <- c()
  ep_table <- c()
  
  for(e in 1:length(rownames(RNA_bin))){
    coord <- rownames(RNA_bin)[e]
    gene <- strsplit(coord, ":")[[1]][4]
    RNA_UP <- 0
    RNA_DOWN <- 0
    RNA_MixedUPDOWN <- 0
    
    RNA_mut <- 0
    
    if(sum(abs(RNA_bin[coord, ])) >0){
      RNA_mut <- 1
      # Direction of RNA expression changes  
      # Determine the original direction of the RNA events    
      events <- sum(abs(RNA_bin[coord, ]))
      sum <- sum(RNA_bin[coord, ])
      
      if(sum == events){
        RNA_UP <- 1
      }     
      
      if(sum == -1*(events)){
        RNA_DOWN <- 1
      }     
      
      if(!abs(sum) == events){
        RNA_MixedUPDOWN <- 1
      }  
    }
    RNA_gens <- "0"
    RNA_inherited <- 0
    
    if(RNA_mut ==1){
      RNA_gens <-colnames(RNA_bin[(which(abs(RNA_bin[coord,])>0))])
      
      if(length(RNA_gens) > 1){
        RNA_gens <- paste(RNA_gens, collapse = "_")
        
        if(2 %in% diff(as.numeric(strsplit(RNA_gens, "_")[[1]]))){
          RNA_inherited <- 1
        } 
      }
    }
    # The RNA coding locus has to a) be present in the small RNA genes so as to get a small RNA signal data point
    # b) be present in the Ahringer data set so as to get associated annotations on chromatin domain distribution
    # First assume the coding locus with RNAseq data is not present in small RNA 
    maps_to_Ahr <- 0
    maps_to_DNA_mut <- 0
    DNA_mut <- 0
    DNA_mut_inherited <- 0
    DNA_mut_gens <- 0
    multiple_DNA_mut_gens  <- 0
    
    if(gene %in% High_dose_DNA_mut_genes == T){
      maps_to_DNA_mut <- 1
      
      if(sum(abs(DNA_mut_bin[gene, ]))>0){
        DNA_mut <- 1
        target_row <- DNA_mut_bin[gene, ]
        events <- sum(abs(DNA_mut_bin[gene, ]))
      }
      if(DNA_mut ==1){
        # Generations
        DNA_mut_gens <- colnames(DNA_mut_bin[(which(abs(DNA_mut_bin[gene, ]) > 0))])
        
        if(length(DNA_mut_gens)> 1){
          DNA_mut_gens <- paste(DNA_mut_gens, collapse = "_")
          multiple_DNA_mut_gens <- 1
        }
        # Inherited
        inherited <- c()
        inherit <- 0 
        each_row <- colnames(DNA_mut_bin[(which(abs(DNA_mut_bin[gene, ]) > 0))])
        
        if(2 %in% diff(as.numeric(each_row))){
          inherit <- 1}
        inherited <- c(inherited, inherit) 
        
        if(sum(inherited > 0)){
          DNA_mut_inherited <- 1
        }
      }}
    time_matched_to_RNA <- 0
    # Determine time matching of small RNA to RNA
    
    if(RNA_mut==1&DNA_mut==1){
      DNA_mut_gens <- unlist(str_split(DNA_mut_gens, "_"))
      rna_gens <- unlist(str_split(RNA_gens, "_"))
      matched_gens <- DNA_mut_gens[DNA_mut_gens %in% rna_gens]
      matched <-  length(which(DNA_mut_gens %in% rna_gens))
      
      if(matched > 0){
        time_matched_to_RNA <- 1
      }
    }
    save <- data.frame(Lineage, gene, coord,RNA_gens, maps_to_Ahr, maps_to_DNA_mut, RNA_mut, RNA_inherited, 
                       DNA_mut, 
                       DNA_mut_inherited, time_matched_to_RNA,DNA_mut_gens,
                       RNA_UP, RNA_DOWN, RNA_MixedUPDOWN 
    )
    integrated_table <- rbind(integrated_table, save)
  }
  High_dose_RNA_and_DNA_mut_integrated_table <- rbind(High_dose_RNA_and_DNA_mut_integrated_table, integrated_table)
}
save(High_dose_RNA_and_DNA_mut_integrated_table,file="High_dose_RNA_and_DNA_mut_integrated_table.Rdata")
write.xlsx(High_dose_RNA_and_DNA_mut_integrated_table,"High_dose_RNA_and_DNA_mut_integrated_table.xlsx")

#The final Table.1 has been written using the data from Control_RNA_and_DNA_mut_integrated_table, Low_dose_RNA_and_DNA_mut_integrated_table and High_dose_RNA_and_DNA_mut_integrated_table
#--------------------------
#####Fig.3.A & Sup.Fig.4.A#####
load("RNA_C1epimutations.Rdata")
load("RNA_C2epimutations.Rdata")
load("RNA_L1epimutations.Rdata")
load("RNA_L2epimutations.Rdata")
load("RNA_H1epimutations.Rdata")
load("RNA_H2epimutations.Rdata")
#Calculate the new epimutations each generation
#Function creation
# Defining the number of transitions UP/DOWN per generation for each data type
# This code counts the number of transitions in chromatin/gene expression state in each generation

# UP transition function
#define function
UP_transition_func<-function(vector_in,input_name){
  Gen_ON <- vector()
  for(i in 2:length(vector_in)){
    if(vector_in[i]==1&vector_in[i-1]== 0){
      Gen_ON[[i]] <-1 # there is an UP transition to ON at generation[i]
    }
    
    if(vector_in[i]== 1&vector_in[i-1]== -1){
      Gen_ON[[i]] <-1 # there is an UP transition to ON at generation[i]
    }
    
    if(vector_in[i]== 1&vector_in[i-1]=="NA"){
      Gen_ON[[i]] <-1 # there is an UP transition to ON at generation[i]
    }
    
    if(vector_in[i]== -1&vector_in[i-1]== 0){
      Gen_ON[[i]] <- 0 # this is a DOWN transition
    }
    
    if(vector_in[i]== -1&vector_in[i-1]== 1){
      Gen_ON[[i]] <- 0 # this is a DOWN epimutation
    }
    
    if(vector_in[i]== -1&vector_in[i-1]== "NA"){
      Gen_ON[[i]] <- 0 # this is a DOWN epimutation
    }
    
    if(vector_in[i]==1&vector_in[i-1]==1){
      Gen_ON[[i]] <- 0 # inside a run of 1 1 epimutations at generation[i]
    }
    
    if(vector_in[i]==-1&vector_in[i-1]==-1){
      Gen_ON[[i]] <- 0 # inside a run of -1 -1 epimutations at generation[i]
    }
    
    if(vector_in[i]==0&vector_in[i-1]==0){
      Gen_ON[[i]] <- 0 # between epimutations at generation[i]
    }
    
    if(vector_in[i]==0&vector_in[i-1]==1){
      Gen_ON[[i]] <- 0 # there is a transition to OFF at generation[i]
    }
    
    if(vector_in[i]==0&vector_in[i-1]==-1){
      Gen_ON[[i]] <- 0 # there is a transition to OFF at generation[i]
    }
    
    if(vector_in[i]==0&vector_in[i-1]=="NA"){
      Gen_ON[[i]] <- 0 # there is a transition to OFF at generation[i]
    } 
    output <- Gen_ON
  }
  
  return(output)
}

# DOWN transition function

#define function
DOWN_transition_func<-function(vector_in,input_name){
  Gen_ON <- vector()
  for(i in 2:length(vector_in)){
    if(vector_in[i]== -1&vector_in[i-1]== 0){
      Gen_ON[[i]] <-1 # there is a DOWN transition to ON at generation[i]
    }
    
    if(vector_in[i]== -1&vector_in[i-1]== 1){
      Gen_ON[[i]] <-1 # there is a DOWN transition to ON at generation[i]
    }
    
    if(vector_in[i]== -1&vector_in[i-1]== "NA"){
      Gen_ON[[i]] <-1 # there is a DOWN transition to ON at generation[i]
    }
    
    if(vector_in[i]== 1&vector_in[i-1]== 0){
      Gen_ON[[i]] <- 0 # this is an UP transition
    }
    
    if(vector_in[i]== 1&vector_in[i-1]== -1){
      Gen_ON[[i]] <- 0 # this is an UP epimutation
    }
    
    if(vector_in[i]==1&vector_in[i-1]==1){
      Gen_ON[[i]] <- 0 # inside a run of 1 1 epimutations at generation[i]
    }
    
    if(vector_in[i]==1&vector_in[i-1]=="NA"){
      Gen_ON[[i]] <- 0 # this is an UP epimutation
    }
    
    if(vector_in[i]==-1&vector_in[i-1]==-1){
      Gen_ON[[i]] <- 0 # inside a run of -1 -1 epimutations at generation[i]
    }
    
    if(vector_in[i]==0&vector_in[i-1]==0){
      Gen_ON[[i]] <- 0 # between epimutations at generation[i]
    }
    
    if(vector_in[i]==0&vector_in[i-1]==1){
      Gen_ON[[i]] <- 0 # there is a transition to OFF at generation[i]
    }
    
    if(vector_in[i]==0&vector_in[i-1]==-1){
      Gen_ON[[i]] <- 0 # there is a transition to OFF at generation[i]
    }
    
    if(vector_in[i]==0&vector_in[i-1]=="NA"){
      Gen_ON[[i]] <- 0 # there is a transition to OFF at generation[i]
    }
    
    output <- Gen_ON
  }
  
  return(output)
}
#Control1 new epi
#UP
UP_output_RNA_C1<- c()
row.names(RNA_C1epimutations)<-RNA_C1epimutations$genes
RNA_C1epimutations<-RNA_C1epimutations[,c(3:13)]
for(i in 1:nrow(RNA_C1epimutations)){
  UP_output_RNA_C1 <-rbind(UP_output_RNA_C1, UP_transition_func(RNA_C1epimutations[i,], input_name=row.names(RNA_C1epimutations)[i]))}
colnames(UP_output_RNA_C1) <- c("0", "2", "4", "6", "8", "10", "12", "14", "16", "18", "20")
row.names(UP_output_RNA_C1) <- rownames(RNA_C1epimutations)
# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 
transitions_at_2 <- sum(UP_output_RNA_C1[,2])
transitions_at_4 <- sum(UP_output_RNA_C1[,3])
transitions_at_6 <- sum(UP_output_RNA_C1[,4])
transitions_at_8 <- sum(UP_output_RNA_C1[,5])
transitions_at_10 <- sum(UP_output_RNA_C1[,6])
transitions_at_12 <- sum(UP_output_RNA_C1[,7])
transitions_at_14 <- sum(UP_output_RNA_C1[,8])
transitions_at_16 <- sum(UP_output_RNA_C1[,9])
transitions_at_18 <- sum(UP_output_RNA_C1[,10])
transitions_at_20 <- sum(UP_output_RNA_C1[,11])

UP_RNA_C1_Table_of_new_epimutations <- rbind(transitions_at_2,
                                             transitions_at_4, 
                                             transitions_at_6,
                                             transitions_at_8, 
                                             transitions_at_10, 
                                             transitions_at_12, 
                                             transitions_at_14, 
                                             transitions_at_16, 
                                             transitions_at_18, 
                                             transitions_at_20)

#DOWN
DOWN_output_RNA_C1<- c()
for(i in 1:nrow(RNA_C1epimutations)){
  DOWN_output_RNA_C1 <-rbind(DOWN_output_RNA_C1, DOWN_transition_func(RNA_C1epimutations[i,], input_name=row.names(RNA_C1epimutations)[i]))}
colnames(DOWN_output_RNA_C1) <- c("0", "2", "4", "6", "8", "10", "12", "14", "16", "18", "20")
row.names(DOWN_output_RNA_C1) <- rownames(RNA_C1epimutations)
# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 
transitions_at_2 <- sum(DOWN_output_RNA_C1[,2])
transitions_at_4 <- sum(DOWN_output_RNA_C1[,3])
transitions_at_6 <- sum(DOWN_output_RNA_C1[,4])
transitions_at_8 <- sum(DOWN_output_RNA_C1[,5])
transitions_at_10 <- sum(DOWN_output_RNA_C1[,6])
transitions_at_12 <- sum(DOWN_output_RNA_C1[,7])
transitions_at_14 <- sum(DOWN_output_RNA_C1[,8])
transitions_at_16 <- sum(DOWN_output_RNA_C1[,9])
transitions_at_18 <- sum(DOWN_output_RNA_C1[,10])
transitions_at_20 <- sum(DOWN_output_RNA_C1[,11])
DOWN_RNA_C1_Table_of_new_epimutations <- rbind(transitions_at_2,
                                               transitions_at_4, 
                                               transitions_at_6,
                                               transitions_at_8, 
                                               transitions_at_10, 
                                               transitions_at_12, 
                                               transitions_at_14, 
                                               transitions_at_16, 
                                               transitions_at_18, 
                                               transitions_at_20)

#Control2 new epi
#UP
UP_output_RNA_C2<- c()
row.names(RNA_C2epimutations)<-RNA_C2epimutations$genes
RNA_C2epimutations<-RNA_C2epimutations[,c(3:13)]
for(i in 1:nrow(RNA_C2epimutations)){
  UP_output_RNA_C2 <-rbind(UP_output_RNA_C2, UP_transition_func(RNA_C2epimutations[i,], input_name=row.names(RNA_C2epimutations)[i]))}
colnames(UP_output_RNA_C2) <- c("0", "2", "4", "6", "8", "10", "12","14", "16", "18", "20")
row.names(UP_output_RNA_C2) <- rownames(RNA_C2epimutations)
# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 
transitions_at_2 <- sum(UP_output_RNA_C2[,2])
transitions_at_4 <- sum(UP_output_RNA_C2[,3])
transitions_at_6 <- sum(UP_output_RNA_C2[,4])
transitions_at_8 <- sum(UP_output_RNA_C2[,5])
transitions_at_10 <- sum(UP_output_RNA_C2[,6])
transitions_at_12 <- NA
transitions_at_14 <- sum(UP_output_RNA_C2[,8])
transitions_at_16 <- sum(UP_output_RNA_C2[,9])
transitions_at_18 <- sum(UP_output_RNA_C2[,10])
transitions_at_20 <- sum(UP_output_RNA_C2[,11])
UP_RNA_C2_Table_of_new_epimutations <- rbind(transitions_at_2,
                                             transitions_at_4, 
                                             transitions_at_6,
                                             transitions_at_8, 
                                             transitions_at_10, 
                                             transitions_at_12, 
                                             transitions_at_14, 
                                             transitions_at_16, 
                                             transitions_at_18, 
                                             transitions_at_20)

#DOWN
DOWN_output_RNA_C2<- c()
for(i in 1:nrow(RNA_C2epimutations)){
  DOWN_output_RNA_C2 <-rbind(DOWN_output_RNA_C2, DOWN_transition_func(RNA_C2epimutations[i,], input_name=row.names(RNA_C2epimutations)[i]))}
colnames(DOWN_output_RNA_C2) <- c("0", "2", "4", "6", "8", "10", "12", "14", "16", "18", "20")
row.names(DOWN_output_RNA_C2) <- rownames(RNA_C2epimutations)
# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 
transitions_at_2 <- sum(DOWN_output_RNA_C2[,2])
transitions_at_4 <- sum(DOWN_output_RNA_C2[,3])
transitions_at_6 <- sum(DOWN_output_RNA_C2[,4])
transitions_at_8 <- sum(DOWN_output_RNA_C2[,5])
transitions_at_10 <- sum(DOWN_output_RNA_C2[,6])
transitions_at_12 <- NA
transitions_at_14 <- sum(DOWN_output_RNA_C2[,8])
transitions_at_16 <- sum(DOWN_output_RNA_C2[,9])
transitions_at_18 <- sum(DOWN_output_RNA_C2[,10])
transitions_at_20 <- sum(DOWN_output_RNA_C2[,11])
DOWN_RNA_C2_Table_of_new_epimutations <- rbind(transitions_at_2,
                                               transitions_at_4, 
                                               transitions_at_6,
                                               transitions_at_8, 
                                               transitions_at_10, 
                                               transitions_at_12, 
                                               transitions_at_14, 
                                               transitions_at_16, 
                                               transitions_at_18, 
                                               transitions_at_20)

#Low1 new epi
#UP
UP_output_RNA_L1<- c()
row.names(RNA_L1epimutations)<-RNA_L1epimutations$genes
RNA_L1epimutations<-RNA_L1epimutations[,c(3:13)]
for(i in 1:nrow(RNA_L1epimutations)){
  UP_output_RNA_L1 <-rbind(UP_output_RNA_L1, UP_transition_func(RNA_L1epimutations[i,], input_name=row.names(RNA_L1epimutations)[i]))}
colnames(UP_output_RNA_L1) <- c("0", "2", "4", "6", "8", "10", "12","14", "16", "18", "20")
row.names(UP_output_RNA_L1) <- rownames(RNA_L1epimutations)
# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 
transitions_at_2 <- sum(UP_output_RNA_L1[,2])
transitions_at_4 <- sum(UP_output_RNA_L1[,3])
transitions_at_6 <- sum(UP_output_RNA_L1[,4])
transitions_at_8 <- sum(UP_output_RNA_L1[,5])
transitions_at_10 <- sum(UP_output_RNA_L1[,6])
transitions_at_12 <- sum(UP_output_RNA_L1[,7])
transitions_at_14 <- sum(UP_output_RNA_L1[,8])
transitions_at_16 <- sum(UP_output_RNA_L1[,9])
transitions_at_18 <- sum(UP_output_RNA_L1[,10])
transitions_at_20 <- sum(UP_output_RNA_L1[,11])
UP_RNA_L1_Table_of_new_epimutations <- rbind(transitions_at_2,
                                             transitions_at_4, 
                                             transitions_at_6,
                                             transitions_at_8, 
                                             transitions_at_10, 
                                             transitions_at_12, 
                                             transitions_at_14, 
                                             transitions_at_16, 
                                             transitions_at_18, 
                                             transitions_at_20)

#DOWN
DOWN_output_RNA_L1<- c()
for(i in 1:nrow(RNA_L1epimutations)){
  DOWN_output_RNA_L1 <-rbind(DOWN_output_RNA_L1, DOWN_transition_func(RNA_L1epimutations[i,], input_name=row.names(RNA_L1epimutations)[i]))}
colnames(DOWN_output_RNA_L1) <- c("0", "2", "4", "6", "8", "10", "12", "14", "16", "18", "20")
row.names(DOWN_output_RNA_L1) <- rownames(RNA_L1epimutations)
# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 
transitions_at_2 <- sum(DOWN_output_RNA_L1[,2])
transitions_at_4 <- sum(DOWN_output_RNA_L1[,3])
transitions_at_6 <- sum(DOWN_output_RNA_L1[,4])
transitions_at_8 <- sum(DOWN_output_RNA_L1[,5])
transitions_at_10 <- sum(DOWN_output_RNA_L1[,6])
transitions_at_12 <- sum(DOWN_output_RNA_L1[,7])
transitions_at_14 <- sum(DOWN_output_RNA_L1[,8])
transitions_at_16 <- sum(DOWN_output_RNA_L1[,9])
transitions_at_18 <- sum(DOWN_output_RNA_L1[,10])
transitions_at_20 <- sum(DOWN_output_RNA_L1[,11])
DOWN_RNA_L1_Table_of_new_epimutations <- rbind(transitions_at_2,
                                               transitions_at_4, 
                                               transitions_at_6,
                                               transitions_at_8, 
                                               transitions_at_10, 
                                               transitions_at_12, 
                                               transitions_at_14, 
                                               transitions_at_16, 
                                               transitions_at_18, 
                                               transitions_at_20)

#Low2 new epi
#UP
UP_output_RNA_L2<- c()
row.names(RNA_L2epimutations)<-RNA_L2epimutations$genes
RNA_L2epimutations<-RNA_L2epimutations[,c(3:13)]
for(i in 1:nrow(RNA_L2epimutations)){
  UP_output_RNA_L2 <-rbind(UP_output_RNA_L2, UP_transition_func(RNA_L2epimutations[i,], input_name=row.names(RNA_L2epimutations)[i]))}
colnames(UP_output_RNA_L2) <- c("0", "2", "4", "6", "8", "10", "12","14", "16", "18", "20")
row.names(UP_output_RNA_L2) <- rownames(RNA_L2epimutations)
# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 
transitions_at_2 <- sum(UP_output_RNA_L2[,2])
transitions_at_4 <- sum(UP_output_RNA_L2[,3])
transitions_at_6 <- sum(UP_output_RNA_L2[,4])
transitions_at_8 <- sum(UP_output_RNA_L2[,5])
transitions_at_10 <- sum(UP_output_RNA_L2[,6])
transitions_at_12 <- sum(UP_output_RNA_L2[,7])
transitions_at_14 <- sum(UP_output_RNA_L2[,8])
transitions_at_16 <- sum(UP_output_RNA_L2[,9])
transitions_at_18 <- sum(UP_output_RNA_L2[,10])
transitions_at_20 <- sum(UP_output_RNA_L2[,11])
UP_RNA_L2_Table_of_new_epimutations <- rbind(transitions_at_2,
                                             transitions_at_4, 
                                             transitions_at_6,
                                             transitions_at_8, 
                                             transitions_at_10, 
                                             transitions_at_12, 
                                             transitions_at_14, 
                                             transitions_at_16, 
                                             transitions_at_18, 
                                             transitions_at_20)

#DOWN
DOWN_output_RNA_L2<- c()
for(i in 1:nrow(RNA_L2epimutations)){
  DOWN_output_RNA_L2 <-rbind(DOWN_output_RNA_L2, DOWN_transition_func(RNA_L2epimutations[i,], input_name=row.names(RNA_L2epimutations)[i]))}
colnames(DOWN_output_RNA_L2) <- c("0", "2", "4", "6", "8", "10", "12", "14", "16", "18", "20")
row.names(DOWN_output_RNA_L2) <- rownames(RNA_L2epimutations)
# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 
transitions_at_2 <- sum(DOWN_output_RNA_L2[,2])
transitions_at_4 <- sum(DOWN_output_RNA_L2[,3])
transitions_at_6 <- sum(DOWN_output_RNA_L2[,4])
transitions_at_8 <- sum(DOWN_output_RNA_L2[,5])
transitions_at_10 <- sum(DOWN_output_RNA_L2[,6])
transitions_at_12 <- sum(DOWN_output_RNA_L2[,7])
transitions_at_14 <- sum(DOWN_output_RNA_L2[,8])
transitions_at_16 <- sum(DOWN_output_RNA_L2[,9])
transitions_at_18 <- sum(DOWN_output_RNA_L2[,10])
transitions_at_20 <- sum(DOWN_output_RNA_L2[,11])
DOWN_RNA_L2_Table_of_new_epimutations <- rbind(transitions_at_2,
                                               transitions_at_4, 
                                               transitions_at_6,
                                               transitions_at_8, 
                                               transitions_at_10, 
                                               transitions_at_12, 
                                               transitions_at_14, 
                                               transitions_at_16, 
                                               transitions_at_18, 
                                               transitions_at_20)

#High1 new epi
#UP
UP_output_RNA_H1<- c()
row.names(RNA_H1epimutations)<-RNA_H1epimutations$genes
RNA_H1epimutations<-RNA_H1epimutations[,c(3:13)]
for(i in 1:nrow(RNA_H1epimutations)){
  UP_output_RNA_H1 <-rbind(UP_output_RNA_H1, UP_transition_func(RNA_H1epimutations[i,], input_name=row.names(RNA_H1epimutations)[i]))}
colnames(UP_output_RNA_H1) <- c("0", "2", "4", "6", "8", "10", "12","14", "16", "18")
row.names(UP_output_RNA_H1) <- rownames(RNA_H1epimutations)
# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 
transitions_at_2 <- sum(UP_output_RNA_H1[,2])
transitions_at_4 <- sum(UP_output_RNA_H1[,3])
transitions_at_6 <- sum(UP_output_RNA_H1[,4])
transitions_at_8 <- sum(UP_output_RNA_H1[,5])
transitions_at_10 <- sum(UP_output_RNA_H1[,6])
transitions_at_12 <- sum(UP_output_RNA_H1[,7])
transitions_at_14 <- sum(UP_output_RNA_H1[,8])
transitions_at_16 <- sum(UP_output_RNA_H1[,9])
transitions_at_18 <- sum(UP_output_RNA_H1[,10])
UP_RNA_H1_Table_of_new_epimutations <- rbind(transitions_at_2,
                                             transitions_at_4, 
                                             transitions_at_6,
                                             transitions_at_8, 
                                             transitions_at_10, 
                                             transitions_at_12, 
                                             transitions_at_14, 
                                             transitions_at_16, 
                                             transitions_at_18)

#DOWN
DOWN_output_RNA_H1<- c()
for(i in 1:nrow(RNA_H1epimutations)){
  DOWN_output_RNA_H1 <-rbind(DOWN_output_RNA_H1, DOWN_transition_func(RNA_H1epimutations[i,], input_name=row.names(RNA_H1epimutations)[i]))}
colnames(DOWN_output_RNA_H1) <- c("0", "2", "4", "6", "8", "10", "12", "14", "16", "18")
row.names(DOWN_output_RNA_H1) <- rownames(RNA_H1epimutations)
# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 
transitions_at_2 <- sum(DOWN_output_RNA_H1[,2])
transitions_at_4 <- sum(DOWN_output_RNA_H1[,3])
transitions_at_6 <- sum(DOWN_output_RNA_H1[,4])
transitions_at_8 <- sum(DOWN_output_RNA_H1[,5])
transitions_at_10 <- sum(DOWN_output_RNA_H1[,6])
transitions_at_12 <- sum(DOWN_output_RNA_H1[,7])
transitions_at_14 <- sum(DOWN_output_RNA_H1[,8])
transitions_at_16 <- sum(DOWN_output_RNA_H1[,9])
transitions_at_18 <- sum(DOWN_output_RNA_H1[,10])
DOWN_RNA_H1_Table_of_new_epimutations <- rbind(transitions_at_2,
                                               transitions_at_4, 
                                               transitions_at_6,
                                               transitions_at_8, 
                                               transitions_at_10, 
                                               transitions_at_12, 
                                               transitions_at_14, 
                                               transitions_at_16, 
                                               transitions_at_18)
#High2 new epi
#UP
UP_output_RNA_H2<- c()
row.names(RNA_H2epimutations)<-RNA_H2epimutations$genes
RNA_H2epimutations<-RNA_H2epimutations[,c(3:13)]
for(i in 1:nrow(RNA_H2epimutations)){
  UP_output_RNA_H2 <-rbind(UP_output_RNA_H2, UP_transition_func(RNA_H2epimutations[i,], input_name=row.names(RNA_H2epimutations)[i]))}
colnames(UP_output_RNA_H2) <- c("0", "2", "4", "6", "8", "10", "12","14", "16", "18","20")
row.names(UP_output_RNA_H2) <- rownames(RNA_H2epimutations)
# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 
transitions_at_2 <- sum(UP_output_RNA_H2[,2])
transitions_at_4 <- sum(UP_output_RNA_H2[,3])
transitions_at_6 <- sum(UP_output_RNA_H2[,4])
transitions_at_8 <- sum(UP_output_RNA_H2[,5])
transitions_at_10 <- sum(UP_output_RNA_H2[,6])
transitions_at_12 <- sum(UP_output_RNA_H2[,7])
transitions_at_14 <- sum(UP_output_RNA_H2[,8])
transitions_at_16 <- sum(UP_output_RNA_H2[,9])
transitions_at_18 <- sum(UP_output_RNA_H2[,10])
transitions_at_20 <- sum(UP_output_RNA_H2[,11])
UP_RNA_H2_Table_of_new_epimutations <- rbind(transitions_at_2,
                                             transitions_at_4, 
                                             transitions_at_6,
                                             transitions_at_8, 
                                             transitions_at_10, 
                                             transitions_at_12, 
                                             transitions_at_14, 
                                             transitions_at_16, 
                                             transitions_at_18,
                                             transitions_at_20)

#DOWN
DOWN_output_RNA_H2<- c()
for(i in 1:nrow(RNA_H2epimutations)){
  DOWN_output_RNA_H2 <-rbind(DOWN_output_RNA_H2, DOWN_transition_func(RNA_H2epimutations[i,], input_name=row.names(RNA_H2epimutations)[i]))}
colnames(DOWN_output_RNA_H2) <- c("0", "2", "4", "6", "8", "10", "12", "14", "16", "18","20")
row.names(DOWN_output_RNA_H2) <- rownames(RNA_H2epimutations)
# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 
transitions_at_2 <- sum(DOWN_output_RNA_H2[,2])
transitions_at_4 <- sum(DOWN_output_RNA_H2[,3])
transitions_at_6 <- sum(DOWN_output_RNA_H2[,4])
transitions_at_8 <- sum(DOWN_output_RNA_H2[,5])
transitions_at_10 <- sum(DOWN_output_RNA_H2[,6])
transitions_at_12 <- sum(DOWN_output_RNA_H2[,7])
transitions_at_14 <- sum(DOWN_output_RNA_H2[,8])
transitions_at_16 <- sum(DOWN_output_RNA_H2[,9])
transitions_at_18 <- sum(DOWN_output_RNA_H2[,10])
transitions_at_20 <- sum(DOWN_output_RNA_H2[,11])
DOWN_RNA_H2_Table_of_new_epimutations <- rbind(transitions_at_2,
                                               transitions_at_4, 
                                               transitions_at_6,
                                               transitions_at_8, 
                                               transitions_at_10, 
                                               transitions_at_12, 
                                               transitions_at_14, 
                                               transitions_at_16, 
                                               transitions_at_18,
                                               transitions_at_20)
#Allnewepi
Transitions<-row.names(UP_RNA_C1_Table_of_new_epimutations)
colnames(UP_RNA_C1_Table_of_new_epimutations)<-c("UP")
colnames(DOWN_RNA_C1_Table_of_new_epimutations)<-c("DOWN")
RNA_allepiC1 <- cbind(Transitions,UP_RNA_C1_Table_of_new_epimutations,DOWN_RNA_C1_Table_of_new_epimutations)
RNA_allepiC1<-as.data.frame(RNA_allepiC1)
RNA_allepiC1 <- RNA_allepiC1 %>%
  # Creating an empty column:
  add_column(Lineage = "C1", .before="UP")
colnames(UP_RNA_C2_Table_of_new_epimutations)<-c("UP")
colnames(DOWN_RNA_C2_Table_of_new_epimutations)<-c("DOWN")
RNA_allepiC2 <- cbind(Transitions,UP_RNA_C2_Table_of_new_epimutations,DOWN_RNA_C2_Table_of_new_epimutations)
RNA_allepiC2<-as.data.frame(RNA_allepiC2)
RNA_allepiC2 <- RNA_allepiC2 %>%
  # Creating an empty column:
  add_column(Lineage = "C2", .before="UP")
RNA_allC <- rbind(RNA_allepiC1, RNA_allepiC2)
RNA_allC <- RNA_allC %>%
  # Creating an empty column:
  add_column(Condition = "Control", .before="Lineage")

colnames(UP_RNA_L1_Table_of_new_epimutations)<-c("UP")
colnames(DOWN_RNA_L1_Table_of_new_epimutations)<-c("DOWN")
RNA_allepiL1 <- cbind(Transitions,UP_RNA_L1_Table_of_new_epimutations,DOWN_RNA_L1_Table_of_new_epimutations)
RNA_allepiL1<-as.data.frame(RNA_allepiL1)
RNA_allepiL1 <- RNA_allepiL1 %>%
  # Creating an empty column:
  add_column(Lineage = "L1", .before="UP")
colnames(UP_RNA_L2_Table_of_new_epimutations)<-c("UP")
colnames(DOWN_RNA_L2_Table_of_new_epimutations)<-c("DOWN")
RNA_allepiL2 <- cbind(Transitions,UP_RNA_L2_Table_of_new_epimutations,DOWN_RNA_L2_Table_of_new_epimutations)
RNA_allepiL2<-as.data.frame(RNA_allepiL2)
RNA_allepiL2 <- RNA_allepiL2 %>%
  # Creating an empty column:
  add_column(Lineage = "L2", .before="UP")
RNA_allL <- rbind(RNA_allepiL1, RNA_allepiL2)
RNA_allL <- RNA_allL %>%
  # Creating an empty column:
  add_column(Condition = "Low dose", .before="Lineage")

colnames(UP_RNA_H1_Table_of_new_epimutations)<-c("UP")
colnames(DOWN_RNA_H1_Table_of_new_epimutations)<-c("DOWN")
RNA_allepiH1 <- cbind(Transitions,UP_RNA_H1_Table_of_new_epimutations,DOWN_RNA_H1_Table_of_new_epimutations)
RNA_allepiH1<-as.data.frame(RNA_allepiH1)
RNA_allepiH1 <- RNA_allepiH1 %>%
  # Creating an empty column:
  add_column(Lineage = "H1", .before="UP")
colnames(UP_RNA_H2_Table_of_new_epimutations)<-c("UP")
colnames(DOWN_RNA_H2_Table_of_new_epimutations)<-c("DOWN")
RNA_allepiH2 <- cbind(Transitions,UP_RNA_H2_Table_of_new_epimutations,DOWN_RNA_H2_Table_of_new_epimutations)
RNA_allepiH2<-as.data.frame(RNA_allepiH2)
RNA_allepiH2 <- RNA_allepiH2 %>%
  # Creating an empty column:
  add_column(Lineage = "H2", .before="UP")
RNA_allH <- rbind(RNA_allepiH1, RNA_allepiH2)
RNA_allH <- RNA_allH %>%
  # Creating an empty column:
  add_column(Condition = "High dose", .before="UP")
RNA_allnewepi<-rbind(RNA_allC,RNA_allL,RNA_allH)
rownames(RNA_allnewepi)<-NULL

RNA_allnewepi$UP <- as.numeric(RNA_allnewepi$UP)
RNA_allnewepi$DOWN <- as.numeric(RNA_allnewepi$DOWN)
RNA_allnewepi$Total=rowSums(cbind(RNA_allnewepi$UP,RNA_allnewepi$DOWN),na.rm=FALSE)

RNA_allnewepi$Condition <- fct_relevel(RNA_allnewepi$Condition, c("Control", "Low dose", "High dose"))

#All RNA epimutations
ggdensity(RNA_allnewepi$Total, 
          main = "Density plot of epimutations_number",
          xlab = "Total number of epimutations_number")
shapiro.test(RNA_allnewepi$Total)
kruskal.test(Total ~ Condition, data = RNA_allnewepi)
dunnTest(Total ~ Condition, data = RNA_allnewepi)

#Fig.3.A
RNA_allLineage_rate <- ggplot(RNA_allnewepi, aes(x=Condition, y=Total, color = Condition,group=Condition)) + 
  geom_boxplot(fatten = 1, lwd = 1, width=0.5)+
  scale_color_manual(values=c("cornflowerblue", "darkgreen", "red"))+
  labs(y = "Number of new epimutations", x = "\nCondition")+
  geom_dotplot(binaxis='y', binwidth = 0.25, stackdir='center',  dotsize=150, binpositions = "all", stackgroups = TRUE)+
  theme_bw()+
  theme_linedraw()+
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(color = "black", size = 20, angle = 0, face = "plain"),
        axis.text.y = element_text(color = "black", size = 20, angle = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 16, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 16, angle = 90, face = "bold"))+
  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))+  
  scale_fill_brewer(palette="Pastel2")+
  ggtitle(paste("RNA_All"))
RNA_allLineage_rate 

#Sup.Fig.4.A
Generation<-c("2","4","6","8","10","12","14","16","18","20","2","4","6","8","10","12","14","16","18","20","2","4","6","8","10","12","14","16","18","20","2","4","6","8","10","12","14","16","18","20","2","4","6","8","10","12","14","16","18","2","4","6","8","10","12","14","16","18","20")
RNA_allnewepi_detailed<-cbind(RNA_allnewepi,Generation)
RNA_allnewepi_detailed$Generation <- fct_relevel(RNA_allnewepi_detailed$Generation, c("2", "4", "6","8","10","12","14","16","18","20"))
RNA_allLineage_rate_detailed <- ggplot(RNA_allnewepi_detailed, aes(x=Generation , y=Total, color = Condition)) + 
  geom_boxplot(fatten = 1, lwd = 1, width=0.5)+
  scale_color_manual(values=c("cornflowerblue", "darkgreen", "red"))+
  labs(y = "Number of new epimutations", x = "\nGeneration")+
  geom_jitter(aes(shape=Lineage, size=8),
              position=position_jitter(width = 0.2,
                                       height = 0.2))+
  theme_bw()+
  theme(axis.text.x = element_text(color = "black", size = 20, angle = 0, face = "plain"),
        axis.text.y = element_text(color = "black", size = 20, angle = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 16, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 16, angle = 90, face = "bold"))+
  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))+  
  scale_fill_brewer(palette="Pastel2")+
  ggtitle(paste(""))
RNA_allLineage_rate_detailed 

#Get the raw data
write.xlsx(RNA_allnewepi,"Data_Fig_3_A.xlsx")
write.xlsx(RNA_allnewepi_detailed,"Data_Sup_Fig_4_A.xlsx")
#--------------------------
#####Fig.3.B & Sup.Fig.4.B#####
##RNA-centric table
load("RNA_C1epimutations.Rdata")
load("RNA_C2epimutations.Rdata")
load("RNA_L1epimutations.Rdata")
load("RNA_L2epimutations.Rdata")
load("RNA_H1epimutations.Rdata")
load("RNA_H2epimutations.Rdata")
# Here we show for each gene for which we have RNAseq data, whether they have enhancer/promoter coordinates as per Ahringer data set, if RNA expression changed and additional data 

# Produce a data table of all genes that have enhancer/promoter coordinates and their RNA expression change 
#Data prep
colnames(RNA_C2epimutations) <- c('genes','condition','0','2','4','6','8','10','F12','14','16','18','20')
RNA_C2epimutationsbis <- subset(RNA_C2epimutations,  select = -F12)
RNA_C2epimutationsbis<-cbind(RNA_C2epimutationsbis, rep(RNA_C2epimutationsbis[9],1))
RNA_C2epimutationsbis<-RNA_C2epimutationsbis[,c(1,2,3,4,5,6,7,8,13,9,10,11,12)]
colnames(RNA_C2epimutationsbis) <- c('genes','condition','0','2','4','6','8','10','12','14','16','18','20')
colnames(RNA_H1epimutations) <- c('genes','condition','0','2','4','6','8','10','12','14','16','18','F20')
RNA_H1epimutationsbis <- subset(RNA_H1epimutations,  select = -F20)

rownames(RNA_C1epimutations)<-RNA_C1epimutations$genes
RNA_C1epimutations <- subset(RNA_C1epimutations,  select = -genes)
RNA_C1epimutations <- subset(RNA_C1epimutations,  select = -condition)

rownames(RNA_C2epimutationsbis)<-RNA_C2epimutationsbis$genes
RNA_C2epimutationsbis <- subset(RNA_C2epimutationsbis,  select = -genes)
RNA_C2epimutationsbis <- subset(RNA_C2epimutationsbis,  select = -condition)

rownames(RNA_L1epimutations)<-RNA_L1epimutations$genes
RNA_L1epimutations <- subset(RNA_L1epimutations,  select = -genes)
RNA_L1epimutations <- subset(RNA_L1epimutations,  select = -condition)

rownames(RNA_L2epimutations)<-RNA_L2epimutations$genes
RNA_L2epimutations <- subset(RNA_L2epimutations,  select = -genes)
RNA_L2epimutations <- subset(RNA_L2epimutations,  select = -condition)

rownames(RNA_H1epimutationsbis)<-RNA_H1epimutationsbis$genes
RNA_H1epimutationsbis <- subset(RNA_H1epimutationsbis,  select = -genes)
RNA_H1epimutationsbis <- subset(RNA_H1epimutationsbis,  select = -condition)

rownames(RNA_H2epimutations)<-RNA_H2epimutations$genes
RNA_H2epimutations <- subset(RNA_H2epimutations,  select = -genes)
RNA_H2epimutations <- subset(RNA_H2epimutations,  select = -condition)

# We will need piRNA cluster genes for this

# Coordinates for piRNA cluster genes (5 - 7.5 & 13 - 17 MB from Ruby et al. 2006) have been intersected with RNA coordinates 

RNA_piRNA_domains <- read.table("RNA_Coords_intersected_piRNA_c.bed", header = FALSE)

piRNA_domains <- RNA_piRNA_domains[RNA_piRNA_domains$V5==1, ]

# When we are considering the genes with chromatin domain annotations we can only consider the genes which map to Ahringer as this is where we get the domain annotations from. 

all_Ahr_RNA_genes <- unique(Ahringer_single_gene_ref_table$Gene)

piRNA_cluster_genes <- intersect(unique(piRNA_domains$V4), all_Ahr_RNA_genes)

#Control
RNA_list <- list(RNA_C1epimutations, RNA_C2epimutationsbis)
lin <- c("C1", "C2")

RNA_Control_integrated_table <- c()

for(x in 1:length(RNA_list)){
  RNA_bin <- RNA_list[[x]]
  Lineage <- lin[[x]]
  integrated_table <- c()
  ep_table <- c()
  
  for(e in 1:length(rownames(RNA_bin))){
    coord <- rownames(RNA_bin)[e]
    gene <- strsplit(coord, ":")[[1]][4]
    AllUP <- 0
    AllDOWN <- 0
    MixedUPDOWN <- 0
    RNA_mut <- 0
    
    if(sum(abs(RNA_bin[coord, ])) >0){
      RNA_mut <- 1
      # Direction of RNA expression changes  
      # Determine the original direction of the RNA events    
      events <- sum(abs(RNA_bin[coord, ]))
      sum <- sum(RNA_bin[coord, ])
      
      if(sum == events){
        AllUP <- 1
      }     
      
      if(sum == -1*(events)){
        AllDOWN <- 1
      }     
      
      if(!abs(sum) == events){
        MixedUPDOWN <- 1
      }  
    }
    
    RNA_gens <- "0"
    RNA_inherited <- 0
    
    if(RNA_mut ==1){
      RNA_gens<-colnames(RNA_bin[(which(abs(RNA_bin[coord,])>0))])
      
      if(length(RNA_gens) > 1){
        RNA_gens <- paste(RNA_gens, collapse = "_")
        
        if(2 %in% diff(as.numeric(strsplit(RNA_gens, "_")[[1]]))){
          RNA_inherited <- 1
        } 
      }
    }
    save <- data.frame(Lineage, gene, coord, RNA_mut, RNA_gens, RNA_inherited, 
                       AllUP, AllDOWN, MixedUPDOWN)
    integrated_table <- rbind(integrated_table, save)
  }
  RNA_Control_integrated_table <- rbind(RNA_Control_integrated_table, integrated_table)
}
colnames(RNA_Control_integrated_table) <- c(
  "Lineage", "gene", "RNA_coord", "is_RNA_exp_change",
  "RNA_mut_gens", "is_RNA_inherited","RNA_All_UP", "RNA_All_DOWN", "RNA_Mixed_Muts"
)
save(RNA_Control_integrated_table,file="RNA_Control_integrated_table.Rdata")
#Low dose
RNA_list <- list(RNA_L1epimutations, RNA_L2epimutations)
lin <- c("L1", "L2")

RNA_Low_dose_integrated_table <- c()

for(x in 1:length(RNA_list)){
  RNA_bin <- RNA_list[[x]]
  Lineage <- lin[[x]]
  integrated_table <- c()
  ep_table <- c()
  
  for(e in 1:length(rownames(RNA_bin))){
    coord <- rownames(RNA_bin)[e]
    gene <- strsplit(coord, ":")[[1]][4]
    AllUP <- 0
    AllDOWN <- 0
    MixedUPDOWN <- 0
    RNA_mut <- 0
    
    if(sum(abs(RNA_bin[coord, ])) >0){
      RNA_mut <- 1
      # Direction of RNA expression changes  
      # Determine the original direction of the RNA events    
      events <- sum(abs(RNA_bin[coord, ]))
      sum <- sum(RNA_bin[coord, ])
      
      if(sum == events){
        AllUP <- 1
      }     
      
      if(sum == -1*(events)){
        AllDOWN <- 1
      }     
      
      if(!abs(sum) == events){
        MixedUPDOWN <- 1
      }  
    }
    
    RNA_gens <- "0"
    RNA_inherited <- 0
    
    if(RNA_mut ==1){
      RNA_gens<-colnames(RNA_bin[(which(abs(RNA_bin[coord,])>0))])
      
      if(length(RNA_gens) > 1){
        RNA_gens <- paste(RNA_gens, collapse = "_")
        
        if(2 %in% diff(as.numeric(strsplit(RNA_gens, "_")[[1]]))){
          RNA_inherited <- 1
        } 
      }
    }
    save <- data.frame(Lineage, gene, coord, RNA_mut, RNA_gens, RNA_inherited, 
                       AllUP, AllDOWN, MixedUPDOWN)
    integrated_table <- rbind(integrated_table, save)
  }
  RNA_Low_dose_integrated_table <- rbind(RNA_Low_dose_integrated_table, integrated_table)
}
colnames(RNA_Low_dose_integrated_table) <- c("Lineage", "gene", "RNA_coord", "is_RNA_exp_change","RNA_mut_gens", "is_RNA_inherited",  "RNA_All_UP", "RNA_All_DOWN", "RNA_Mixed_Muts")
save(RNA_Low_dose_integrated_table,file="RNA_Low_dose_integrated_table.Rdata")
#High dose
RNA_list <- list(RNA_H1epimutationsbis, RNA_H2epimutations)
lin <- c("H1", "H2")

RNA_High_dose_integrated_table <- c()

for(x in 1:length(RNA_list)){
  RNA_bin <- RNA_list[[x]]
  Lineage <- lin[[x]]
  integrated_table <- c()
  ep_table <- c()
  
  for(e in 1:length(rownames(RNA_bin))){
    coord <- rownames(RNA_bin)[e]
    gene <- strsplit(coord, ":")[[1]][4]
    AllUP <- 0
    AllDOWN <- 0
    MixedUPDOWN <- 0
    RNA_mut <- 0
    
    if(sum(abs(RNA_bin[coord, ])) >0){
      RNA_mut <- 1
      # Direction of RNA expression changes  
      # Determine the original direction of the RNA events    
      events <- sum(abs(RNA_bin[coord, ]))
      sum <- sum(RNA_bin[coord, ])
      
      if(sum == events){
        AllUP <- 1
      }     
      
      if(sum == -1*(events)){
        AllDOWN <- 1
      }     
      
      if(!abs(sum) == events){
        MixedUPDOWN <- 1
      }  
    }
    
    RNA_gens <- "0"
    RNA_inherited <- 0
    
    if(RNA_mut ==1){
      RNA_gens<-colnames(RNA_bin[(which(abs(RNA_bin[coord,])>0))])
      if(length(RNA_gens) > 1){
        RNA_gens <- paste(RNA_gens, collapse = "_")
        
        if(2 %in% diff(as.numeric(strsplit(RNA_gens, "_")[[1]]))){
          RNA_inherited <- 1
        } 
      }
    }
    save <- data.frame(Lineage, gene, coord, RNA_mut, RNA_gens, RNA_inherited, 
                       AllUP, AllDOWN, MixedUPDOWN)
    integrated_table <- rbind(integrated_table, save)
  }
  RNA_High_dose_integrated_table <- rbind(RNA_High_dose_integrated_table, integrated_table)
}
colnames(RNA_High_dose_integrated_table) <- c("Lineage", "gene", "RNA_coord", "is_RNA_exp_change","RNA_mut_gens", "is_RNA_inherited",  "RNA_All_UP", "RNA_All_DOWN", "RNA_Mixed_Muts")
save(RNA_High_dose_integrated_table,file="RNA_High_dose_integrated_table.Rdata")
#Calc duration epimut

RNA_all_epi<-list(RNA_C1epimutations,RNA_C2epimutationsbis,RNA_L1epimutations,RNA_L2epimutations,RNA_H1epimutationsbis,RNA_H2epimutations)
names(RNA_all_epi)<-c("C1","C2","L1","L2","H1","H2")
#Calculation
RNA_All_output_list_with_missing_as_muts <- list()

for(e in 1:length(RNA_all_epi)){
  
  main_frame <- as.data.frame(RNA_all_epi[[e]])
  
  Overall_output <- c()
  
  for(t in 1:nrow(main_frame)){
    
    vector_in <- main_frame[t, ]
    gen_names <- as.numeric(colnames(vector_in))
    
    number_transitions <- 0
    onset_gen <- 0
    tempL <- 0
    output_frame <- c()
    
    is_up <- 0
    is_down <- 0
    complete <- 0
    
    # first deal with no epimutations in the lineage for a locus/gene
    
    if(sum(abs(vector_in)) == 0){
      
      length = 0
      is_up <- 0
      is_down <- 0
      gene_name <- rownames(vector_in)
      gene <- strsplit(gene_name, ":")[[1]][4]
      
      number_transitions <- 0 
      onset_gen <- 0
      
      Lineage <- names(RNA_all_epi[e])
      complete <- 0
      
      save_mut <- data.frame(gene_name, gene, number_transitions, onset_gen, length,  complete,  Lineage, is_up, is_down)   
      
      output_frame <- rbind(output_frame, save_mut)
      
    }
    
    else
      
      if(sum(abs(vector_in)) > 0){
        
        for(i in 2:length(vector_in)){
          
          
          # if it is an UP transition, turn off any down transitions and start new epimutation
          
          if(vector_in[i]==1&vector_in[i-1]== 0|vector_in[i]== 1&vector_in[i-1]== -1) {
            
            if(is_down == 1){ # a down epimutation that is now turning off
              
              gene_name <- rownames(vector_in)
              gene <- strsplit(gene_name, ":")[[1]][4]
              
              number_transitions <- number_transitions + 1 
              
              length <-tempL
              
              Lineage <- names(RNA_all_epi[e])
              
              complete <- 1
              
              save_mut <- data.frame(gene_name, gene, number_transitions, onset_gen, length,  complete, Lineage, is_up, is_down)   
              
              output_frame <- rbind(output_frame, save_mut)
              
            }
            
            # start the UP transition
            
            is_up <- 1
            is_down <- 0
            tempL <- 1
            
            number_transitions <- 1 # there is a transition to an UP epimutation at generation[i]
            onset_gen <- gen_names[i]
            
            onset_check <- i -1
          }
          
          # if it is a new DOWN transition, turn off any up transitions and start new epimutation
          
          if(vector_in[i]==-1&vector_in[i-1]== 0|vector_in[i]== -1&vector_in[i-1]== 1) {
            
            if(is_up == 1){ # an up epimutation that is now turning off
              
              gene_name <- rownames(vector_in)
              gene <- strsplit(gene_name, ":")[[1]][4]
              
              number_transitions <- number_transitions + 1
              
              length <- tempL
              
              Lineage <- names(RNA_all_epi[e])
              
              complete <- 1
              
              save_mut <- data.frame(gene_name, gene, number_transitions, onset_gen, length,  complete, Lineage, is_up, is_down)   
              
              output_frame <- rbind(output_frame, save_mut)
              
            }
            
            # start the down transition
            
            is_down <- 1
            is_up <- 0
            tempL <- 1
            number_transitions <- 1 # there is a transition to a DOWN epimutation at generation[i]
            onset_gen <- gen_names[i]
            onset_check <- i-1
          }
          
          # if it is an UP epimutation, does it continue, does it turn off, does it extend to last gen?
          
          if(is_up==1){
            
            
            if(vector_in[i]==1&vector_in[i-1]==1)
              
            {tempL<-tempL+(gen_names[i]-gen_names[i-1])}  # the epimutation continues
            
            if(vector_in[i]==0&vector_in[i-1]==1){
              
              gene_name <- rownames(vector_in)
              gene <- strsplit(gene_name, ":")[[1]][4]
              
              number_transitions <- number_transitions +1  # there is a transition to OFF from an UP epimutation at generation[i]
              
              length <- tempL
              
              Lineage <- names(RNA_all_epi[e])
              
              complete <- 1
              
              save_mut <- data.frame(gene_name, gene, number_transitions, onset_gen, length,  complete, Lineage, is_up, is_down)   
              
              output_frame <- rbind(output_frame, save_mut)
              
              is_up <- 0
              is_down <- 0
              tempL <- 0
              length <- 0
              complete <- 0
              
            }
            
            # deal with an up epimutation extending to end of lineage
            
            if(vector_in[i]==1&i==length(vector_in)){
              
              gene_name <- rownames(vector_in)
              
              gene <- strsplit(gene_name, ":")[[1]][4]
              
              number_transitions <- number_transitions + 0
              
              length <- tempL
              
              Lineage <- names(RNA_all_epi[e])
              
              complete <- 0
              
              save_mut <- data.frame(gene_name, gene, number_transitions, onset_gen, length,  complete, Lineage, is_up, is_down)   
              
              output_frame <- rbind(output_frame, save_mut)
              
              is_up <- 0
              is_down <- 0
              length <- 0
              tempL <- 0
              complete <- 0
              
            }
            
          }
          
          # if it is a DOWN epimutation, does it continue, does it turn off, does it extend to last gen?
          
          if(is_down==1){
            
            if(vector_in[i]==-1&vector_in[i-1]==-1)
              
            {tempL<-tempL+(gen_names[i]-gen_names[i-1])}
            
            
            if(vector_in[i]==0&vector_in[i-1]==-1){
              
              gene_name <- rownames(vector_in)
              
              gene <- strsplit(gene_name, ":")[[1]][4]
              
              number_transitions <- number_transitions +1  # there is a transition to OFF from a DOWN epimutation at generation[i]
              
              length <- tempL
              
              Lineage <- names(RNA_all_epi[e])
              
              complete <- 1
              
              save_mut <- data.frame(gene_name, gene, number_transitions, onset_gen, length,  complete, Lineage, is_up, is_down)   
              
              output_frame <- rbind(output_frame, save_mut)
              
              is_up <- 0
              is_down <- 0
              tempL <- 0
              complete <- 0
              
            }
            
            # deal with a down epimutation extending to end of lineage
            
            if(vector_in[i]==-1&i==length(vector_in)){
              
              gene_name <- rownames(vector_in)
              
              gene <- strsplit(gene_name, ":")[[1]][4]
              
              number_transitions <- number_transitions + 0
              
              length <- tempL
              
              Lineage <- names(RNA_all_epi[e])
              
              complete <- 0
              
              save_mut <- data.frame(gene_name, gene, number_transitions, onset_gen, length,  complete, Lineage, is_up, is_down)   
              output_frame <- rbind(output_frame, save_mut)
              
              is_up <- 0
              is_down <- 0
              tempL <- 0
              complete <- 0
              
            }
            
          }    
          
        }} 
    
    Overall_output <- rbind(Overall_output, output_frame)  
    
  }
  
  RNA_All_output_list_with_missing_as_muts[[e]] <- Overall_output
}

RNA_allepiduration<-c()
names(RNA_All_output_list_with_missing_as_muts) <- names(RNA_all_epi)
RNA_All_output_list_with_missing_as_muts$C1 <- RNA_All_output_list_with_missing_as_muts$C1 %>%
  # Creating an empty column:
  add_column(Condition = "Control")
RNA_allepiduration<-RNA_All_output_list_with_missing_as_muts$C1
RNA_All_output_list_with_missing_as_muts$C2 <- RNA_All_output_list_with_missing_as_muts$C2 %>%
  # Creating an empty column:
  add_column(Condition = "Control")
RNA_allepiduration<-rbind(RNA_allepiduration, RNA_All_output_list_with_missing_as_muts$C2)
RNA_All_output_list_with_missing_as_muts$L1 <- RNA_All_output_list_with_missing_as_muts$L1 %>%
  # Creating an empty column:
  add_column(Condition = "Low dose")
RNA_allepiduration<-rbind(RNA_allepiduration, RNA_All_output_list_with_missing_as_muts$L1)
RNA_All_output_list_with_missing_as_muts$L2 <- RNA_All_output_list_with_missing_as_muts$L2 %>%
  # Creating an empty column:
  add_column(Condition = "Low dose")
RNA_allepiduration<-rbind(RNA_allepiduration, RNA_All_output_list_with_missing_as_muts$L2)
RNA_All_output_list_with_missing_as_muts$H1 <- RNA_All_output_list_with_missing_as_muts$H1 %>%
  # Creating an empty column:
  add_column(Condition = "High dose")
RNA_allepiduration<-rbind(RNA_allepiduration, RNA_All_output_list_with_missing_as_muts$H1)
RNA_All_output_list_with_missing_as_muts$H2 <- RNA_All_output_list_with_missing_as_muts$H2 %>%
  # Creating an empty column:
  add_column(Condition = "High dose")
RNA_allepiduration<-rbind(RNA_allepiduration, RNA_All_output_list_with_missing_as_muts$H2)
save(RNA_allepiduration,file = "RNA_allepiduration.Rdata")

RNA_epimut<-subset(RNA_allepiduration,RNA_allepiduration$length !=0)
save(RNA_epimut,file="RNA_epimut.Rdata")

########Survival analysis##
#Look at condition effect
RNA_epimut_C<-subset(RNA_epimut,RNA_epimut$Condition=="Control")
RNA_epimut_L<-subset(RNA_epimut,RNA_epimut$Condition=="Low dose")
RNA_epimut_H<-subset(RNA_epimut,RNA_epimut$Condition=="High dose")

fit <- survfit(Surv(length, complete) ~ Condition, data = RNA_epimut)

CoxMod<-coxph(Surv(length,complete)~Condition,data=RNA_epimut)
ggforest(model=CoxMod,data=RNA_epimut,fontsize = 0.8, noDigits = 2)
#Fig 3. B
ggsurvplot <- ggsurvplot(fit,  RNA_epimut, censor = T, break.time.by= 1, pval=TRUE,pval.coord = c(0, 0.03),
                         font.main = c(16, "bold"),
                         font.x = c(20,"bold"),
                         palette = c("#3399FF", "#006600","2E9FDF"),
                         font.y = c(20,"bold"),
                         font.legend = 20,
                         legend.labs=c("Control","Low dose","High dose"),
                         font.tickslab = 18)+
  ggtitle("Changes survival all conditions")+
  xlab("Time (generations)")
ggsurvplot$plot +theme(plot.title = element_text(hjust = 0.5)) 
#Sup Fig 4. B
CoxMod<-coxph(Surv(length,complete)~Lineage,data=RNA_epimut)
ggforest(model=CoxMod,data=RNA_epimut,fontsize = 0.8, noDigits = 2)
fit_2 <- survfit(Surv(length, complete) ~ Lineage, data = RNA_epimut)
ggsurvplot <- ggsurvplot(fit_2,  RNA_epimut, censor = T, break.time.by= 1, pval=TRUE,pval.coord = c(0, 0.03),                         font.main = c(16, "bold"),
                         font.x = c(20,"bold"),
                         palette = c("#3399FF","blue","2E9FDF","red", "green","#006600"),
                         font.y = c(20,"bold"),
                         font.legend = 20,
                         legend.labs=c("C1","C2","H1","H2","L1","L2"),
                         font.tickslab = 18)+
  ggtitle("Changes survival all conditions")+
  xlab("Time (generations)")
ggsurvplot$plot +theme(plot.title = element_text(hjust = 0.5)) 
#Get raw data
write.xlsx(RNA_epimut,"Data_Fig_3_B_&_Sup_Fig_4_B.xlsx")
#--------------------------
#####Fig.3.C & Sup.Fig.4.C - Focus on TRUE epimutations (lasting >1 generation)#####
#Look at epimutations starting at each generation for control
tot_C<-nrow(RNA_epimut_C)
tot_C_2<-nrow(subset(RNA_epimut_C,RNA_epimut_C$onset_gen==2))
lg_epimut_C_2<-nrow(subset(RNA_epimut_C,RNA_epimut_C$onset_gen==2 & RNA_epimut_C$length>1))
per_lg_epimut_C_2<-(lg_epimut_C_2/tot_C_2)*100
tot_C_4<-nrow(subset(RNA_epimut_C,RNA_epimut_C$onset_gen==4))
lg_epimut_C_4<-nrow(subset(RNA_epimut_C,RNA_epimut_C$onset_gen==4 & RNA_epimut_C$length>1))
per_lg_epimut_C_4<-(lg_epimut_C_4/tot_C_4)*100
tot_C_6<-nrow(subset(RNA_epimut_C,RNA_epimut_C$onset_gen==6))
lg_epimut_C_6<-nrow(subset(RNA_epimut_C,RNA_epimut_C$onset_gen==6 & RNA_epimut_C$length>1))
per_lg_epimut_C_6<-(lg_epimut_C_6/tot_C_6)*100
tot_C_8<-nrow(subset(RNA_epimut_C,RNA_epimut_C$onset_gen==8))
lg_epimut_C_8<-nrow(subset(RNA_epimut_C,RNA_epimut_C$onset_gen==8 & RNA_epimut_C$length>1))
per_lg_epimut_C_8<-(lg_epimut_C_8/tot_C_8)*100
tot_C_10<-nrow(subset(RNA_epimut_C,RNA_epimut_C$onset_gen==10))
lg_epimut_C_10<-nrow(subset(RNA_epimut_C,RNA_epimut_C$onset_gen==10 & RNA_epimut_C$length>1))
per_lg_epimut_C_10<-(lg_epimut_C_10/tot_C_10)*100
tot_C_12<-nrow(subset(RNA_epimut_C,RNA_epimut_C$onset_gen==12))
lg_epimut_C_12<-nrow(subset(RNA_epimut_C,RNA_epimut_C$onset_gen==12 & RNA_epimut_C$length>1))
per_lg_epimut_C_12<-(lg_epimut_C_12/tot_C_12)*100
tot_C_14<-nrow(subset(RNA_epimut_C,RNA_epimut_C$onset_gen==14))
lg_epimut_C_14<-nrow(subset(RNA_epimut_C,RNA_epimut_C$onset_gen==14 & RNA_epimut_C$length>1))
per_lg_epimut_C_14<-(lg_epimut_C_14/tot_C_14)*100
tot_C_16<-nrow(subset(RNA_epimut_C,RNA_epimut_C$onset_gen==16))
lg_epimut_C_16<-nrow(subset(RNA_epimut_C,RNA_epimut_C$onset_gen==16 & RNA_epimut_C$length>1))
per_lg_epimut_C_16<-(lg_epimut_C_16/tot_C_16)*100
tot_C_18<-nrow(subset(RNA_epimut_C,RNA_epimut_C$onset_gen==18))
lg_epimut_C_18<-nrow(subset(RNA_epimut_C,RNA_epimut_C$onset_gen==18 & RNA_epimut_C$length>1))
per_lg_epimut_C_18<-(lg_epimut_C_18/tot_C_18)*100

#Look at epimutations starting at each generation for low dose
tot_L<-nrow(RNA_epimut_L)
tot_L_2<-nrow(subset(RNA_epimut_C,RNA_epimut_L$onset_gen==2))
lg_epimut_L_2<-nrow(subset(RNA_epimut_L,RNA_epimut_L$onset_gen==2 & RNA_epimut_L$length>1))
per_lg_epimut_L_2<-(lg_epimut_L_2/tot_L_2)*100
tot_L_4<-nrow(subset(RNA_epimut_L,RNA_epimut_L$onset_gen==4))
lg_epimut_L_4<-nrow(subset(RNA_epimut_L,RNA_epimut_L$onset_gen==4 & RNA_epimut_L$length>1))
per_lg_epimut_L_4<-(lg_epimut_L_4/tot_L_4)*100
tot_L_6<-nrow(subset(RNA_epimut_L,RNA_epimut_L$onset_gen==6))
lg_epimut_L_6<-nrow(subset(RNA_epimut_L,RNA_epimut_L$onset_gen==6 & RNA_epimut_L$length>1))
per_lg_epimut_L_6<-(lg_epimut_L_6/tot_L_6)*100
tot_L_8<-nrow(subset(RNA_epimut_L,RNA_epimut_L$onset_gen==8))
lg_epimut_L_8<-nrow(subset(RNA_epimut_L,RNA_epimut_L$onset_gen==8 & RNA_epimut_L$length>1))
per_lg_epimut_L_8<-(lg_epimut_L_8/tot_L_8)*100
tot_L_10<-nrow(subset(RNA_epimut_L,RNA_epimut_L$onset_gen==10))
lg_epimut_L_10<-nrow(subset(RNA_epimut_L,RNA_epimut_L$onset_gen==10 & RNA_epimut_L$length>1))
per_lg_epimut_L_10<-(lg_epimut_L_10/tot_L_10)*100
tot_L_12<-nrow(subset(RNA_epimut_L,RNA_epimut_L$onset_gen==12))
lg_epimut_L_12<-nrow(subset(RNA_epimut_L,RNA_epimut_L$onset_gen==12 & RNA_epimut_L$length>1))
per_lg_epimut_L_12<-(lg_epimut_L_12/tot_L_12)*100
tot_L_14<-nrow(subset(RNA_epimut_L,RNA_epimut_L$onset_gen==14))
lg_epimut_L_14<-nrow(subset(RNA_epimut_L,RNA_epimut_L$onset_gen==14 & RNA_epimut_L$length>1))
per_lg_epimut_L_14<-(lg_epimut_L_14/tot_L_14)*100
tot_L_16<-nrow(subset(RNA_epimut_L,RNA_epimut_L$onset_gen==16))
lg_epimut_L_16<-nrow(subset(RNA_epimut_L,RNA_epimut_L$onset_gen==16 & RNA_epimut_L$length>1))
per_lg_epimut_L_16<-(lg_epimut_L_16/tot_L_16)*100
tot_L_18<-nrow(subset(RNA_epimut_L,RNA_epimut_L$onset_gen==18))
lg_epimut_L_18<-nrow(subset(RNA_epimut_L,RNA_epimut_L$onset_gen==18 & RNA_epimut_L$length>1))
per_lg_epimut_L_18<-(lg_epimut_L_18/tot_L_18)*100

#Look at epimutations starting at each generation for high dose
tot_H<-nrow(RNA_epimut_H)
tot_H_2<-nrow(subset(RNA_epimut_H,RNA_epimut_H$onset_gen==2))
lg_epimut_H_2<-nrow(subset(RNA_epimut_H,RNA_epimut_H$onset_gen==2 & RNA_epimut_H$length>1))
per_lg_epimut_H_2<-(lg_epimut_H_2/tot_H_2)*100
tot_H_4<-nrow(subset(RNA_epimut_H,RNA_epimut_H$onset_gen==4))
lg_epimut_H_4<-nrow(subset(RNA_epimut_H,RNA_epimut_H$onset_gen==4 & RNA_epimut_H$length>1))
per_lg_epimut_H_4<-(lg_epimut_H_4/tot_H_4)*100
tot_H_6<-nrow(subset(RNA_epimut_H,RNA_epimut_H$onset_gen==6))
lg_epimut_H_6<-nrow(subset(RNA_epimut_H,RNA_epimut_H$onset_gen==6 & RNA_epimut_H$length>1))
per_lg_epimut_H_6<-(lg_epimut_H_6/tot_H_6)*100
tot_H_8<-nrow(subset(RNA_epimut_H,RNA_epimut_H$onset_gen==8))
lg_epimut_H_8<-nrow(subset(RNA_epimut_H,RNA_epimut_H$onset_gen==8 & RNA_epimut_H$length>1))
per_lg_epimut_H_8<-(lg_epimut_H_8/tot_H_8)*100
tot_H_10<-nrow(subset(RNA_epimut_H,RNA_epimut_H$onset_gen==10))
lg_epimut_H_10<-nrow(subset(RNA_epimut_H,RNA_epimut_H$onset_gen==10 & RNA_epimut_H$length>1))
per_lg_epimut_H_10<-(lg_epimut_H_10/tot_H_10)*100
tot_H_12<-nrow(subset(RNA_epimut_H,RNA_epimut_H$onset_gen==12))
lg_epimut_H_12<-nrow(subset(RNA_epimut_H,RNA_epimut_H$onset_gen==12 & RNA_epimut_H$length>1))
per_lg_epimut_H_12<-(lg_epimut_H_12/tot_H_12)*100
tot_H_14<-nrow(subset(RNA_epimut_H,RNA_epimut_H$onset_gen==14))
lg_epimut_H_14<-nrow(subset(RNA_epimut_H,RNA_epimut_H$onset_gen==14 & RNA_epimut_H$length>1))
per_lg_epimut_H_14<-(lg_epimut_H_14/tot_H_14)*100
tot_H_16<-nrow(subset(RNA_epimut_H,RNA_epimut_H$onset_gen==16))
lg_epimut_H_16<-nrow(subset(RNA_epimut_H,RNA_epimut_H$onset_gen==16 & RNA_epimut_H$length>1))
per_lg_epimut_H_16<-(lg_epimut_H_16/tot_H_16)*100
tot_H_18<-nrow(subset(RNA_epimut_H,RNA_epimut_H$onset_gen==18))
lg_epimut_H_18<-nrow(subset(RNA_epimut_H,RNA_epimut_H$onset_gen==18 & RNA_epimut_H$length>1))
per_lg_epimut_H_18<-(lg_epimut_H_18/tot_H_18)*100

#Create a table with generations as replicates
Lg_RNA_epimut_table <- data.frame (Condition  = c("Control", "Control","Control","Control","Control","Control","Control","Control","Control",
                                                  "Low dose","Low dose","Low dose","Low dose","Low dose","Low dose","Low dose","Low dose","Low dose",
                                                  "High dose","High dose","High dose","High dose","High dose","High dose","High dose","High dose","High dose"),
                                   Generation = c("2", "4","6","8","10","12","14","16","18",
                                                  "2", "4","6","8","10","12","14","16","18",
                                                  "2", "4","6","8","10","12","14","16","18"),
                                   Total_new_epimut = c(tot_C_2,tot_C_4,tot_C_6,tot_C_8,tot_C_10,tot_C_12,tot_C_14,tot_C_16,tot_C_18,
                                                        tot_L_2,tot_L_4,tot_L_6,tot_L_8,tot_L_10,tot_L_12,tot_L_14,tot_L_16,tot_L_18,
                                                        tot_H_2,tot_H_4,tot_H_6,tot_H_8,tot_H_10,tot_H_12,tot_H_14,tot_H_16,tot_H_18),
                                   Lg_epimut = c(lg_epimut_C_2,lg_epimut_C_4,lg_epimut_C_6,lg_epimut_C_8,lg_epimut_C_10,lg_epimut_C_12,lg_epimut_C_14,lg_epimut_C_16,lg_epimut_C_18,
                                                 lg_epimut_L_2,lg_epimut_L_4,lg_epimut_L_6,lg_epimut_L_8,lg_epimut_L_10,lg_epimut_L_12,lg_epimut_L_14,lg_epimut_L_16,lg_epimut_L_18,
                                                 lg_epimut_H_2,lg_epimut_H_4,lg_epimut_H_6,lg_epimut_H_8,lg_epimut_H_10,lg_epimut_H_12,lg_epimut_H_14,lg_epimut_H_16,lg_epimut_H_18),
                                   Percentage=c(per_lg_epimut_C_2,per_lg_epimut_C_4,per_lg_epimut_C_6,per_lg_epimut_C_8,per_lg_epimut_C_10,per_lg_epimut_C_12,per_lg_epimut_C_14,per_lg_epimut_C_16,per_lg_epimut_C_18,
                                                per_lg_epimut_L_2,per_lg_epimut_L_4,per_lg_epimut_L_6,per_lg_epimut_L_8,per_lg_epimut_L_10,per_lg_epimut_L_12,per_lg_epimut_L_14,per_lg_epimut_L_16,per_lg_epimut_L_18,
                                                per_lg_epimut_H_2,per_lg_epimut_H_4,per_lg_epimut_H_6,per_lg_epimut_H_8,per_lg_epimut_H_10,per_lg_epimut_H_12,per_lg_epimut_H_14,per_lg_epimut_H_16,per_lg_epimut_H_18))

my_sum <- Lg_RNA_epimut_table %>%
  group_by(Condition) %>% 
  dplyr::summarize(count = n(),
      mean=mean(Percentage),
      sd=sd(Percentage)
  )

#Fig.3.C
my_sum$Condition <- fct_relevel(my_sum$Condition, c("Control","Low dose","High dose"))
ggplot(my_sum) +
  geom_bar( aes(x=Condition, y=mean), stat="identity", fill=c("cornflowerblue","red","darkgreen"), alpha=0.5) +
  geom_errorbar( aes(x=Condition, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="black", alpha=0.9, size=1.5) +
  labs(x="Condition", y="Mean % of new epimutations arising each \n generation  and lasting >1 generation") +
  theme(axis.text.x = element_text(color = "black", size = 20, angle = 0, face = "plain"),
        axis.text.y = element_text(color = "black", size = 20, angle = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 16, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 16, angle = 90, face = "bold"))
ggtitle("")
#Get the raw data
write.xlsx(my_sum,"Data_Fig_3_C.xlsx")
##Sup.Fig.4.C
#Look at epimutations starting at each generation for C1
RNA_epimut_C1<-subset(RNA_epimut_C,Lineage == "C1")
tot_C1<-nrow(RNA_epimut_C1)
tot_C1_2<-nrow(subset(RNA_epimut_C1,RNA_epimut_C1$onset_gen==2))
lg_epimut_C1_2<-nrow(subset(RNA_epimut_C1,RNA_epimut_C1$onset_gen==2 & RNA_epimut_C1$length>1))
per_lg_epimut_C1_2<-(lg_epimut_C1_2/tot_C1_2)*100
tot_C1_4<-nrow(subset(RNA_epimut_C1,RNA_epimut_C1$onset_gen==4))
lg_epimut_C1_4<-nrow(subset(RNA_epimut_C1,RNA_epimut_C1$onset_gen==4 & RNA_epimut_C1$length>1))
per_lg_epimut_C1_4<-(lg_epimut_C1_4/tot_C1_4)*100
tot_C1_6<-nrow(subset(RNA_epimut_C1,RNA_epimut_C1$onset_gen==6))
lg_epimut_C1_6<-nrow(subset(RNA_epimut_C1,RNA_epimut_C1$onset_gen==6 & RNA_epimut_C1$length>1))
per_lg_epimut_C1_6<-(lg_epimut_C1_6/tot_C1_6)*100
tot_C1_8<-nrow(subset(RNA_epimut_C1,RNA_epimut_C1$onset_gen==8))
lg_epimut_C1_8<-nrow(subset(RNA_epimut_C1,RNA_epimut_C1$onset_gen==8 & RNA_epimut_C1$length>1))
per_lg_epimut_C1_8<-(lg_epimut_C1_8/tot_C1_8)*100
tot_C1_10<-nrow(subset(RNA_epimut_C1,RNA_epimut_C1$onset_gen==10))
lg_epimut_C1_10<-nrow(subset(RNA_epimut_C1,RNA_epimut_C1$onset_gen==10 & RNA_epimut_C1$length>1))
per_lg_epimut_C1_10<-(lg_epimut_C1_10/tot_C1_10)*100
tot_C1_12<-nrow(subset(RNA_epimut_C1,RNA_epimut_C1$onset_gen==12))
lg_epimut_C1_12<-nrow(subset(RNA_epimut_C1,RNA_epimut_C1$onset_gen==12 & RNA_epimut_C1$length>1))
per_lg_epimut_C1_12<-(lg_epimut_C1_12/tot_C1_12)*100
tot_C1_14<-nrow(subset(RNA_epimut_C1,RNA_epimut_C1$onset_gen==14))
lg_epimut_C1_14<-nrow(subset(RNA_epimut_C1,RNA_epimut_C1$onset_gen==14 & RNA_epimut_C1$length>1))
per_lg_epimut_C1_14<-(lg_epimut_C1_14/tot_C1_14)*100
tot_C1_16<-nrow(subset(RNA_epimut_C1,RNA_epimut_C1$onset_gen==16))
lg_epimut_C1_16<-nrow(subset(RNA_epimut_C1,RNA_epimut_C1$onset_gen==16 & RNA_epimut_C1$length>1))
per_lg_epimut_C1_16<-(lg_epimut_C1_16/tot_C1_16)*100
tot_C1_18<-nrow(subset(RNA_epimut_C1,RNA_epimut_C1$onset_gen==18))
lg_epimut_C1_18<-nrow(subset(RNA_epimut_C1,RNA_epimut_C1$onset_gen==18 & RNA_epimut_C1$length>1))
per_lg_epimut_C1_18<-(lg_epimut_C1_18/tot_C1_18)*100
#Look at epimutations starting at each generation for C2
RNA_epimut_C2<-subset(RNA_epimut_C,Lineage == "C2")
tot_C2<-nrow(RNA_epimut_C2)
tot_C2_2<-nrow(subset(RNA_epimut_C2,RNA_epimut_C2$onset_gen==2))
lg_epimut_C2_2<-nrow(subset(RNA_epimut_C2,RNA_epimut_C2$onset_gen==2 & RNA_epimut_C2$length>1))
per_lg_epimut_C2_2<-(lg_epimut_C2_2/tot_C2_2)*100
tot_C2_4<-nrow(subset(RNA_epimut_C2,RNA_epimut_C2$onset_gen==4))
lg_epimut_C2_4<-nrow(subset(RNA_epimut_C2,RNA_epimut_C2$onset_gen==4 & RNA_epimut_C2$length>1))
per_lg_epimut_C2_4<-(lg_epimut_C2_4/tot_C2_4)*100
tot_C2_6<-nrow(subset(RNA_epimut_C2,RNA_epimut_C2$onset_gen==6))
lg_epimut_C2_6<-nrow(subset(RNA_epimut_C2,RNA_epimut_C2$onset_gen==6 & RNA_epimut_C2$length>1))
per_lg_epimut_C2_6<-(lg_epimut_C2_6/tot_C2_6)*100
tot_C2_8<-nrow(subset(RNA_epimut_C2,RNA_epimut_C2$onset_gen==8))
lg_epimut_C2_8<-nrow(subset(RNA_epimut_C2,RNA_epimut_C2$onset_gen==8 & RNA_epimut_C2$length>1))
per_lg_epimut_C2_8<-(lg_epimut_C2_8/tot_C2_8)*100
tot_C2_10<-nrow(subset(RNA_epimut_C2,RNA_epimut_C2$onset_gen==10))
lg_epimut_C2_10<-nrow(subset(RNA_epimut_C2,RNA_epimut_C2$onset_gen==10 & RNA_epimut_C2$length>1))
per_lg_epimut_C2_10<-(lg_epimut_C2_10/tot_C2_10)*100
tot_C2_12<-nrow(subset(RNA_epimut_C2,RNA_epimut_C2$onset_gen==12))
lg_epimut_C2_12<-nrow(subset(RNA_epimut_C2,RNA_epimut_C2$onset_gen==12 & RNA_epimut_C2$length>1))
per_lg_epimut_C2_12<-(lg_epimut_C2_12/tot_C2_12)*100
tot_C2_14<-nrow(subset(RNA_epimut_C2,RNA_epimut_C2$onset_gen==14))
lg_epimut_C2_14<-nrow(subset(RNA_epimut_C2,RNA_epimut_C2$onset_gen==14 & RNA_epimut_C2$length>1))
per_lg_epimut_C2_14<-(lg_epimut_C2_14/tot_C2_14)*100
tot_C2_16<-nrow(subset(RNA_epimut_C2,RNA_epimut_C2$onset_gen==16))
lg_epimut_C2_16<-nrow(subset(RNA_epimut_C2,RNA_epimut_C2$onset_gen==16 & RNA_epimut_C2$length>1))
per_lg_epimut_C2_16<-(lg_epimut_C2_16/tot_C2_16)*100
tot_C2_18<-nrow(subset(RNA_epimut_C2,RNA_epimut_C2$onset_gen==18))
lg_epimut_C2_18<-nrow(subset(RNA_epimut_C2,RNA_epimut_C2$onset_gen==18 & RNA_epimut_C2$length>1))
per_lg_epimut_C2_18<-(lg_epimut_C2_18/tot_C2_18)*100

#Look at epimutations starting at each generation for low dose 1
RNA_epimut_L1<-subset(RNA_epimut_L,Lineage == "L1")
tot_L1<-nrow(RNA_epimut_L1)
tot_L1_2<-nrow(subset(RNA_epimut_L1,RNA_epimut_L1$onset_gen==2))
lg_epimut_L1_2<-nrow(subset(RNA_epimut_L1,RNA_epimut_L1$onset_gen==2 & RNA_epimut_L1$length>1))
per_lg_epimut_L1_2<-(lg_epimut_L1_2/tot_L1_2)*100
tot_L1_4<-nrow(subset(RNA_epimut_L1,RNA_epimut_L1$onset_gen==4))
lg_epimut_L1_4<-nrow(subset(RNA_epimut_L1,RNA_epimut_L1$onset_gen==4 & RNA_epimut_L1$length>1))
per_lg_epimut_L1_4<-(lg_epimut_L1_4/tot_L1_4)*100
tot_L1_6<-nrow(subset(RNA_epimut_L1,RNA_epimut_L1$onset_gen==6))
lg_epimut_L1_6<-nrow(subset(RNA_epimut_L1,RNA_epimut_L1$onset_gen==6 & RNA_epimut_L1$length>1))
per_lg_epimut_L1_6<-(lg_epimut_L1_6/tot_L1_6)*100
tot_L1_8<-nrow(subset(RNA_epimut_L1,RNA_epimut_L1$onset_gen==8))
lg_epimut_L1_8<-nrow(subset(RNA_epimut_L1,RNA_epimut_L1$onset_gen==8 & RNA_epimut_L1$length>1))
per_lg_epimut_L1_8<-(lg_epimut_L1_8/tot_L1_8)*100
tot_L1_10<-nrow(subset(RNA_epimut_L1,RNA_epimut_L1$onset_gen==10))
lg_epimut_L1_10<-nrow(subset(RNA_epimut_L1,RNA_epimut_L1$onset_gen==10 & RNA_epimut_L1$length>1))
per_lg_epimut_L1_10<-(lg_epimut_L1_10/tot_L1_10)*100
tot_L1_12<-nrow(subset(RNA_epimut_L1,RNA_epimut_L1$onset_gen==12))
lg_epimut_L1_12<-nrow(subset(RNA_epimut_L1,RNA_epimut_L1$onset_gen==12 & RNA_epimut_L1$length>1))
per_lg_epimut_L1_12<-(lg_epimut_L1_12/tot_L1_12)*100
tot_L1_14<-nrow(subset(RNA_epimut_L1,RNA_epimut_L1$onset_gen==14))
lg_epimut_L1_14<-nrow(subset(RNA_epimut_L1,RNA_epimut_L1$onset_gen==14 & RNA_epimut_L1$length>1))
per_lg_epimut_L1_14<-(lg_epimut_L1_14/tot_L1_14)*100
tot_L1_16<-nrow(subset(RNA_epimut_L1,RNA_epimut_L1$onset_gen==16))
lg_epimut_L1_16<-nrow(subset(RNA_epimut_L1,RNA_epimut_L1$onset_gen==16 & RNA_epimut_L1$length>1))
per_lg_epimut_L1_16<-(lg_epimut_L1_16/tot_L1_16)*100
tot_L1_18<-nrow(subset(RNA_epimut_L1,RNA_epimut_L1$onset_gen==18))
lg_epimut_L1_18<-nrow(subset(RNA_epimut_L1,RNA_epimut_L1$onset_gen==18 & RNA_epimut_L1$length>1))
per_lg_epimut_L1_18<-(lg_epimut_L1_18/tot_L1_18)*100
#Look at epimutations starting at each generation for L2
RNA_epimut_L2<-subset(RNA_epimut_L,Lineage == "L2")
tot_L2<-nrow(RNA_epimut_L2)
tot_L2_2<-nrow(subset(RNA_epimut_L2,RNA_epimut_L2$onset_gen==2))
lg_epimut_L2_2<-nrow(subset(RNA_epimut_L2,RNA_epimut_L2$onset_gen==2 & RNA_epimut_L2$length>1))
per_lg_epimut_L2_2<-(lg_epimut_L2_2/tot_L2_2)*100
tot_L2_4<-nrow(subset(RNA_epimut_L2,RNA_epimut_L2$onset_gen==4))
lg_epimut_L2_4<-nrow(subset(RNA_epimut_L2,RNA_epimut_L2$onset_gen==4 & RNA_epimut_L2$length>1))
per_lg_epimut_L2_4<-(lg_epimut_L2_4/tot_L2_4)*100
tot_L2_6<-nrow(subset(RNA_epimut_L2,RNA_epimut_L2$onset_gen==6))
lg_epimut_L2_6<-nrow(subset(RNA_epimut_L2,RNA_epimut_L2$onset_gen==6 & RNA_epimut_L2$length>1))
per_lg_epimut_L2_6<-(lg_epimut_L2_6/tot_L2_6)*100
tot_L2_8<-nrow(subset(RNA_epimut_L2,RNA_epimut_L2$onset_gen==8))
lg_epimut_L2_8<-nrow(subset(RNA_epimut_L2,RNA_epimut_L2$onset_gen==8 & RNA_epimut_L2$length>1))
per_lg_epimut_L2_8<-(lg_epimut_L2_8/tot_L2_8)*100
tot_L2_10<-nrow(subset(RNA_epimut_L2,RNA_epimut_L2$onset_gen==10))
lg_epimut_L2_10<-nrow(subset(RNA_epimut_L2,RNA_epimut_L2$onset_gen==10 & RNA_epimut_L2$length>1))
per_lg_epimut_L2_10<-(lg_epimut_L2_10/tot_L2_10)*100
tot_L2_12<-nrow(subset(RNA_epimut_L2,RNA_epimut_L2$onset_gen==12))
lg_epimut_L2_12<-nrow(subset(RNA_epimut_L2,RNA_epimut_L2$onset_gen==12 & RNA_epimut_L2$length>1))
per_lg_epimut_L2_12<-(lg_epimut_L2_12/tot_L2_12)*100
tot_L2_14<-nrow(subset(RNA_epimut_L2,RNA_epimut_L2$onset_gen==14))
lg_epimut_L2_14<-nrow(subset(RNA_epimut_L2,RNA_epimut_L2$onset_gen==14 & RNA_epimut_L2$length>1))
per_lg_epimut_L2_14<-(lg_epimut_L2_14/tot_L2_14)*100
tot_L2_16<-nrow(subset(RNA_epimut_L2,RNA_epimut_L2$onset_gen==16))
lg_epimut_L2_16<-nrow(subset(RNA_epimut_L2,RNA_epimut_L2$onset_gen==16 & RNA_epimut_L2$length>1))
per_lg_epimut_L2_16<-(lg_epimut_L2_16/tot_L2_16)*100
tot_L2_18<-nrow(subset(RNA_epimut_L2,RNA_epimut_L2$onset_gen==18))
lg_epimut_L2_18<-nrow(subset(RNA_epimut_L2,RNA_epimut_L2$onset_gen==18 & RNA_epimut_L2$length>1))
per_lg_epimut_L2_18<-(lg_epimut_L2_18/tot_L2_18)*100

#Look at epimutations starting at each generation for high dose 1
RNA_epimut_H1<-subset(RNA_epimut_H,Lineage == "H1")
tot_H1<-nrow(RNA_epimut_H1)
tot_H1_2<-nrow(subset(RNA_epimut_H1,RNA_epimut_H1$onset_gen==2))
lg_epimut_H1_2<-nrow(subset(RNA_epimut_H1,RNA_epimut_H1$onset_gen==2 & RNA_epimut_H1$length>1))
per_lg_epimut_H1_2<-(lg_epimut_H1_2/tot_H1_2)*100
tot_H1_4<-nrow(subset(RNA_epimut_H1,RNA_epimut_H1$onset_gen==4))
lg_epimut_H1_4<-nrow(subset(RNA_epimut_H1,RNA_epimut_H1$onset_gen==4 & RNA_epimut_H1$length>1))
per_lg_epimut_H1_4<-(lg_epimut_H1_4/tot_H1_4)*100
tot_H1_6<-nrow(subset(RNA_epimut_H1,RNA_epimut_H1$onset_gen==6))
lg_epimut_H1_6<-nrow(subset(RNA_epimut_H1,RNA_epimut_H1$onset_gen==6 & RNA_epimut_H1$length>1))
per_lg_epimut_H1_6<-(lg_epimut_H1_6/tot_H1_6)*100
tot_H1_8<-nrow(subset(RNA_epimut_H1,RNA_epimut_H1$onset_gen==8))
lg_epimut_H1_8<-nrow(subset(RNA_epimut_H1,RNA_epimut_H1$onset_gen==8 & RNA_epimut_H1$length>1))
per_lg_epimut_H1_8<-(lg_epimut_H1_8/tot_H1_8)*100
tot_H1_10<-nrow(subset(RNA_epimut_H1,RNA_epimut_H1$onset_gen==10))
lg_epimut_H1_10<-nrow(subset(RNA_epimut_H1,RNA_epimut_H1$onset_gen==10 & RNA_epimut_H1$length>1))
per_lg_epimut_H1_10<-(lg_epimut_H1_10/tot_H1_10)*100
tot_H1_12<-nrow(subset(RNA_epimut_H1,RNA_epimut_H1$onset_gen==12))
lg_epimut_H1_12<-nrow(subset(RNA_epimut_H1,RNA_epimut_H1$onset_gen==12 & RNA_epimut_H1$length>1))
per_lg_epimut_H1_12<-(lg_epimut_H1_12/tot_H1_12)*100
tot_H1_14<-nrow(subset(RNA_epimut_H1,RNA_epimut_H1$onset_gen==14))
lg_epimut_H1_14<-nrow(subset(RNA_epimut_H1,RNA_epimut_H1$onset_gen==14 & RNA_epimut_H1$length>1))
per_lg_epimut_H1_14<-(lg_epimut_H1_14/tot_H1_14)*100
tot_H1_16<-nrow(subset(RNA_epimut_H1,RNA_epimut_H1$onset_gen==16))
lg_epimut_H1_16<-nrow(subset(RNA_epimut_H1,RNA_epimut_H1$onset_gen==16 & RNA_epimut_H1$length>1))
per_lg_epimut_H1_16<-(lg_epimut_H1_16/tot_H1_16)*100
tot_H1_18<-nrow(subset(RNA_epimut_H1,RNA_epimut_H1$onset_gen==18))
lg_epimut_H1_18<-nrow(subset(RNA_epimut_H1,RNA_epimut_H1$onset_gen==18 & RNA_epimut_H1$length>1))
per_lg_epimut_H1_18<-(lg_epimut_H1_18/tot_H1_18)*100
#Look at epimutations starting at each generation for H2
RNA_epimut_H2<-subset(RNA_epimut_H,Lineage == "H2")
tot_H2<-nrow(RNA_epimut_H2)
tot_H2_2<-nrow(subset(RNA_epimut_H2,RNA_epimut_H2$onset_gen==2))
lg_epimut_H2_2<-nrow(subset(RNA_epimut_H2,RNA_epimut_H2$onset_gen==2 & RNA_epimut_H2$length>1))
per_lg_epimut_H2_2<-(lg_epimut_H2_2/tot_H2_2)*100
tot_H2_4<-nrow(subset(RNA_epimut_H2,RNA_epimut_H2$onset_gen==4))
lg_epimut_H2_4<-nrow(subset(RNA_epimut_H2,RNA_epimut_H2$onset_gen==4 & RNA_epimut_H2$length>1))
per_lg_epimut_H2_4<-(lg_epimut_H2_4/tot_H2_4)*100
tot_H2_6<-nrow(subset(RNA_epimut_H2,RNA_epimut_H2$onset_gen==6))
lg_epimut_H2_6<-nrow(subset(RNA_epimut_H2,RNA_epimut_H2$onset_gen==6 & RNA_epimut_H2$length>1))
per_lg_epimut_H2_6<-(lg_epimut_H2_6/tot_H2_6)*100
tot_H2_8<-nrow(subset(RNA_epimut_H2,RNA_epimut_H2$onset_gen==8))
lg_epimut_H2_8<-nrow(subset(RNA_epimut_H2,RNA_epimut_H2$onset_gen==8 & RNA_epimut_H2$length>1))
per_lg_epimut_H2_8<-(lg_epimut_H2_8/tot_H2_8)*100
tot_H2_10<-nrow(subset(RNA_epimut_H2,RNA_epimut_H2$onset_gen==10))
lg_epimut_H2_10<-nrow(subset(RNA_epimut_H2,RNA_epimut_H2$onset_gen==10 & RNA_epimut_H2$length>1))
per_lg_epimut_H2_10<-(lg_epimut_H2_10/tot_H2_10)*100
tot_H2_12<-nrow(subset(RNA_epimut_H2,RNA_epimut_H2$onset_gen==12))
lg_epimut_H2_12<-nrow(subset(RNA_epimut_H2,RNA_epimut_H2$onset_gen==12 & RNA_epimut_H2$length>1))
per_lg_epimut_H2_12<-(lg_epimut_H2_12/tot_H2_12)*100
tot_H2_14<-nrow(subset(RNA_epimut_H2,RNA_epimut_H2$onset_gen==14))
lg_epimut_H2_14<-nrow(subset(RNA_epimut_H2,RNA_epimut_H2$onset_gen==14 & RNA_epimut_H2$length>1))
per_lg_epimut_H2_14<-(lg_epimut_H2_14/tot_H2_14)*100
tot_H2_16<-nrow(subset(RNA_epimut_H2,RNA_epimut_H2$onset_gen==16))
lg_epimut_H2_16<-nrow(subset(RNA_epimut_H2,RNA_epimut_H2$onset_gen==16 & RNA_epimut_H2$length>1))
per_lg_epimut_H2_16<-(lg_epimut_H2_16/tot_H2_16)*100
tot_H2_18<-nrow(subset(RNA_epimut_H2,RNA_epimut_H2$onset_gen==18))
lg_epimut_H2_18<-nrow(subset(RNA_epimut_H2,RNA_epimut_H2$onset_gen==18 & RNA_epimut_H2$length>1))
per_lg_epimut_H2_18<-(lg_epimut_H2_18/tot_H2_18)*100

#Create a table with generations as replicates
Lg_RNA_epimut_table_lineage <- data.frame (Lineage = c("C1", "C1","C1","C1","C1","C1","C1","C1","C1",
                                                       "C2", "C2","C2","C2","C2","C2","C2","C2",
                                                       "L1", "L1","L1","L1","L1","L1","L1","L1","L1",
                                                       "L2", "L2","L2","L2","L2","L2","L2","L2","L2",
                                                       "H1", "H1","H1","H1","H1","H1","H1","H1","H1",
                                                       "H2", "H2","H2","H2","H2","H2","H2","H2","H2"),
                                   Condition  = c("Control", "Control","Control","Control","Control","Control","Control","Control","Control",
                                                  "Control", "Control","Control","Control","Control","Control","Control","Control",
                                                  "Low dose","Low dose","Low dose","Low dose","Low dose","Low dose","Low dose","Low dose","Low dose",
                                                  "Low dose","Low dose","Low dose","Low dose","Low dose","Low dose","Low dose","Low dose","Low dose",
                                                  "High dose","High dose","High dose","High dose","High dose","High dose","High dose","High dose","High dose",
                                                  "High dose","High dose","High dose","High dose","High dose","High dose","High dose","High dose","High dose"),
                                   Generation = c("2", "4","6","8","10","12","14","16","18",
                                                  "2", "4","6","8","10","12","16","18",
                                                  "2", "4","6","8","10","12","14","16","18",
                                                  "2", "4","6","8","10","12","14","16","18",
                                                  "2", "4","6","8","10","12","14","16","18",
                                                  "2", "4","6","8","10","12","14","16","18"),
                                   Total_new_epimut = c(tot_C1_2,tot_C1_4,tot_C1_6,tot_C1_8,tot_C1_10,tot_C1_12,tot_C1_14,tot_C1_16,tot_C1_18,
                                                        tot_C2_2,tot_C2_4,tot_C2_6,tot_C2_8,tot_C2_10,tot_C2_12,tot_C2_16,tot_C2_18,
                                                        tot_L1_2,tot_L1_4,tot_L1_6,tot_L1_8,tot_L1_10,tot_L1_12,tot_L1_14,tot_L1_16,tot_L1_18,
                                                        tot_L2_2,tot_L2_4,tot_L2_6,tot_L2_8,tot_L2_10,tot_L2_12,tot_L2_14,tot_L2_16,tot_L2_18,
                                                        tot_H1_2,tot_H1_4,tot_H1_6,tot_H1_8,tot_H1_10,tot_H1_12,tot_H1_14,tot_H1_16,tot_H1_18,
                                                        tot_H2_2,tot_H2_4,tot_H2_6,tot_H2_8,tot_H2_10,tot_H2_12,tot_H2_14,tot_H2_16,tot_H2_18),
                                   Lg_epimut = c(lg_epimut_C1_2,lg_epimut_C1_4,lg_epimut_C1_6,lg_epimut_C1_8,lg_epimut_C1_10,lg_epimut_C1_12,lg_epimut_C1_14,lg_epimut_C1_16,lg_epimut_C1_18,
                                                 lg_epimut_C2_2,lg_epimut_C2_4,lg_epimut_C2_6,lg_epimut_C2_8,lg_epimut_C2_10,lg_epimut_C2_12,lg_epimut_C2_16,lg_epimut_C2_18,
                                                 lg_epimut_L1_2,lg_epimut_L1_4,lg_epimut_L1_6,lg_epimut_L1_8,lg_epimut_L1_10,lg_epimut_L1_12,lg_epimut_L1_14,lg_epimut_L1_16,lg_epimut_L1_18,
                                                 lg_epimut_L2_2,lg_epimut_L2_4,lg_epimut_L2_6,lg_epimut_L2_8,lg_epimut_L2_10,lg_epimut_L2_12,lg_epimut_L2_14,lg_epimut_L2_16,lg_epimut_L2_18,
                                                 lg_epimut_H1_2,lg_epimut_H1_4,lg_epimut_H1_6,lg_epimut_H1_8,lg_epimut_H1_10,lg_epimut_H1_12,lg_epimut_H1_14,lg_epimut_H1_16,lg_epimut_H1_18,
                                                 lg_epimut_H2_2,lg_epimut_H2_4,lg_epimut_H2_6,lg_epimut_H2_8,lg_epimut_H2_10,lg_epimut_H2_12,lg_epimut_H2_14,lg_epimut_H2_16,lg_epimut_H2_18),
                                   Percentage=c(per_lg_epimut_C1_2,per_lg_epimut_C1_4,per_lg_epimut_C1_6,per_lg_epimut_C1_8,per_lg_epimut_C1_10,per_lg_epimut_C1_12,per_lg_epimut_C1_14,per_lg_epimut_C1_16,per_lg_epimut_C1_18,
                                                per_lg_epimut_C2_2,per_lg_epimut_C2_4,per_lg_epimut_C2_6,per_lg_epimut_C2_8,per_lg_epimut_C2_10,per_lg_epimut_C2_12,per_lg_epimut_C2_16,per_lg_epimut_C2_18,
                                                per_lg_epimut_L1_2,per_lg_epimut_L1_4,per_lg_epimut_L1_6,per_lg_epimut_L1_8,per_lg_epimut_L1_10,per_lg_epimut_L1_12,per_lg_epimut_L1_14,per_lg_epimut_L1_16,per_lg_epimut_L1_18,
                                                per_lg_epimut_L2_2,per_lg_epimut_L2_4,per_lg_epimut_L2_6,per_lg_epimut_L2_8,per_lg_epimut_L2_10,per_lg_epimut_L2_12,per_lg_epimut_L2_14,per_lg_epimut_L2_16,per_lg_epimut_L2_18,
                                                per_lg_epimut_H1_2,per_lg_epimut_H1_4,per_lg_epimut_H1_6,per_lg_epimut_H1_8,per_lg_epimut_H1_10,per_lg_epimut_H1_12,per_lg_epimut_H1_14,per_lg_epimut_H1_16,per_lg_epimut_H1_18,
                                                per_lg_epimut_H2_2,per_lg_epimut_H2_4,per_lg_epimut_H2_6,per_lg_epimut_H2_8,per_lg_epimut_H2_10,per_lg_epimut_H2_12,per_lg_epimut_H2_14,per_lg_epimut_H2_16,per_lg_epimut_H2_18))

my_sum_2 <- Lg_RNA_epimut_table_lineage %>%
  group_by(Lineage) %>% 
  dplyr::summarize(count = n(),
                   mean=mean(Percentage),
                   sd=sd(Percentage)
  )

#Sup.Fig.4.C
my_sum_2$Lineage <- fct_relevel(my_sum_2$Lineage, c("C1","C2","L1","L2","H1","H2"))
ggplot(my_sum_2) +
  geom_bar( aes(x=Lineage, y=mean), stat="identity", fill=c("#3399FF","blue","red","2E9FDF","green","#006600"), alpha=0.5) +
  geom_errorbar( aes(x=Lineage, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="black", alpha=0.9, size=1.5) +
  labs(x="Lineage", y="Mean % of new epimutations arising each \n generation  and lasting >1 generation") +
  theme(axis.text.x = element_text(color = "black", size = 20, angle = 0, face = "plain"),
        axis.text.y = element_text(color = "black", size = 20, angle = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 16, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 16, angle = 90, face = "bold"))

#Get the raw data
write.xlsx(my_sum_2,"Data_Sup_Fig_4_C.xlsx")
#--------------------------
#####Fig.3.D - GO term enrichment analysis for RNA epimutations####
#Use of EnrichR
library(enrichR)

# tell the enrichR package to use worms
setEnrichrSite("WormEnrichr")

# find gene set databases available
dbs <- listEnrichrDbs()

# We need to tell enrichR which databases (from the selection in dbs) we would like to query.
# We can start with KEGG and the three 2018 GO databases

chosendbs <- c("KEGG_2019",
               "GO_Cellular_Component_2018",
               "GO_Molecular_Function_2018",
               "GO_Biological_Process_2018", 
               "InterPro_Domains_2019")
#selection of TRUE epimutations
tot_C_sup_1<-subset(RNA_epimut_C,RNA_epimut_C$length>1)
tot_L_sup_1<-subset(RNA_epimut_L,RNA_epimut_L$length>1)
tot_H_sup_1<-subset(RNA_epimut_H,RNA_epimut_H$length>1)
#RNA Control epimut>1 vs all genes
# The test gene list
Sim_test_list <- unique(tot_C_sup_1$gene)

# The background list
Sim_bg_list <- unique(C_elegans_gene_names[,4])

# Remove any genes from the background list that feature in the test list

intersect <- intersect(Sim_bg_list, Sim_test_list)
Sim_bg_list <- Sim_bg_list[Sim_bg_list %in% intersect == F]

# Send both lists to enrichr
Sim_epimutated_enriched_list  <- enrichr(Sim_test_list, chosendbs)
Sim_epimutated_bg_list  <- enrichr(Sim_bg_list, chosendbs)

Sim_test_bg_list <- list(Sim_epimutated_enriched_list, Sim_epimutated_bg_list)

# Annotate with name of libraries
#write.xlsx(Sim_test_bg_results[[2]],"Sim_test_bg_results_overlap.xlsx")
Sim_test_bg_results <- lapply(Sim_test_bg_list, function(x){
  
  lapply(1:5, function(j){
    
    if(nrow(x[[j]]) > 0){
      x[[j]][1:nrow(x[[j]]), "library"] <- names(Sim_test_bg_list[[1]])[j]
      return(x[[j]])
    } 
    
  })
  
})

# then we combine them and put them in order of significance using the Old.Adjusted.P.value
Sim_test_bg_RESULTS <- lapply(Sim_test_bg_results, function(x){
  
  tempdf <- do.call(rbind, x)
  tempdf[order(tempdf$Old.Adjusted.P.value), ]
  
})

split_value_list <- list()

for(i in 1:length(Sim_test_bg_RESULTS)){
  
  split <- strsplit(Sim_test_bg_RESULTS[[i]][, 2], "/")
  
  split_value <- c()
  
  for(j in 1:length(split)){
    
    save <- split[[j]][1]
    
    split_value <- as.numeric(c(split_value, save))
  }
  
  Sim_test_bg_RESULTS[[i]]$Overlap <- split_value
  
}

# get the test_term, the test_not, the mean(sample_term) and the sample_not

# test_term is overlap
# test_not is length of gene set minus overlap

# then add 1 to all if either test_term OR sample_term = 0
# then only pursue if test_term &/OR sample_term > 4
# let the terms of interest be the terms for the enrichr result for the test list, i.e. with simultaneously epimutated genes

terms_of_interest <- Sim_test_bg_RESULTS[[1]]$Term

# and if the term is not present in the other gene list then assume there is no enrichment for that term in those lists

Test_List <- Sim_test_bg_RESULTS[[1]]
BG_List <- Sim_test_bg_RESULTS[[2]]

table_comp <- c()

for(i in 1:length(terms_of_interest)){
  
  Term_select <- terms_of_interest[[i]]
  sample_term <- 0
  
  if(Term_select %in% BG_List$Term == T){
    
    sample_term <- BG_List[which(BG_List$Term %in% Term_select), 2]}
  sample_not <- length(Sim_bg_list) - sample_term
  
  test_term <- Test_List[which(Test_List$Term %in% Term_select), 2]
  test_not <- length(Sim_test_list) - test_term
  
  # save the components in a row in a new table
  
  component_row <- data.frame(test_term, test_not, sample_term, sample_not)
  
  rownames(component_row) <- Term_select
  
  table_comp <- rbind(table_comp, component_row)
  colnames(table_comp) <- c("test_term", "test_not", "sample_term", "sample_not")
  
}

# Now we will adjust  all values by adding 1 if there are any zeros

zeros <- colSums(table_comp == 0)

if(sum(zeros) > 0){
  
  table_comp <-  table_comp + 1
}

# Now modify the table to remove any rows where test term or mean sample term are < 5
# Now apply the function only if at least one category has 5 or more genes
# It is only reasonable to assess for relative enrichment/depletion if either the test list or sample list has > 5 genes in that ontology category

table_comp <- table_comp[(table_comp$test_term < 5 & table_comp$sample_term < 5)==F, ]

# Now we have found terms where there are at least 5 genes in 1 test or sample category
# we have adjusted the data for each gene set for any 0 values relative to the specific test list
# Now we can do the analysis but only on terms present in each table_comp list

Pval_col <- c()
OR_col <- c()

for(q in 1:nrow(table_comp)){
  
  test_term <- table_comp[q, 1] 
  test_not <- table_comp[q, 2]
  mean_Sample_term <- table_comp[q, 3]
  sample_not <- table_comp[q, 4]
  
  contingency_table <-
    rbind(c(test_term, test_not),
          c(mean_Sample_term, sample_not))
  FT_out <- fisher.test(contingency_table)
  pval <- FT_out$p.value
  OR <- FT_out$estimate
  
  Pval_col <- c(Pval_col, pval)
  OR_col<- c(OR_col, OR)
  
}    

OR_pval_frame <- data.frame(Pval_col, OR_col)

colnames(OR_pval_frame) <- c("Pval", "OR")
rownames(OR_pval_frame) <- rownames(table_comp)

# Bonferroni correction

OR_pval_frame[, 1] <- p.adjust(OR_pval_frame[, 1], method="bonferroni", n = nrow(OR_pval_frame)) 

# Make plots for this

# i) In each plot order them by the largest significant enrichment out of all of the comparisons 

OR_pval_frame_ordered <- OR_pval_frame[order(OR_pval_frame[, 2], decreasing = T), ]

# Make a bubble plot 

P_1 <- as.numeric(OR_pval_frame_ordered[, 1])

OR_1 <- log10(as.numeric(OR_pval_frame_ordered[, 2]))

neg_log_p_1 <- -log(P_1)

bubble_table_1 <- data.frame(rownames(OR_pval_frame_ordered), OR_1,  P_1, neg_log_p_1)

colnames(bubble_table_1) <- c("Term", "log_OR", "p_val", "neg_log_P")

bubble_table_1$Term <- factor(bubble_table_1$Term, levels = rev(bubble_table_1$Term))

find_size <- c()
find_alpha <- c()

for(i in 1:nrow(bubble_table_1)){
  
  alpha <- 0.3
  size_point <- bubble_table_1[i,4]
  
  if(as.numeric(bubble_table_1[i, 3]) < 0.1){
    alpha <- 1
    size_point <- bubble_table_1[i, 4]
    
  }
  find_alpha <- c(find_alpha, alpha)  
  find_size <- c(find_size, size_point)
}

bubble_table <- cbind(bubble_table_1, find_alpha, find_size)

# just plot first 10

allsignifCepi<-subset(bubble_table,bubble_table$p_val<0.05)
write.xlsx(allsignifCepi,"allsignifCepi.xlsx")
bubble_table_C <- bubble_table[1:10, ]
trunc_terms_C <- as.character(bubble_table$Term)

#RNA Low dose epimut>1 vs all genes
# The test gene list
Sim_test_list <- unique(tot_L_sup_1$gene)

# The background list
Sim_bg_list <- unique(C_elegans_gene_names[,4])

# Remove any genes from the background list that feature in the test list

intersect <- intersect(Sim_bg_list, Sim_test_list)
Sim_bg_list <- Sim_bg_list[Sim_bg_list %in% intersect == F]

# Send both lists to enrichr
Sim_epimutated_enriched_list  <- enrichr(Sim_test_list, chosendbs)
Sim_epimutated_bg_list  <- enrichr(Sim_bg_list, chosendbs)

Sim_test_bg_list <- list(Sim_epimutated_enriched_list, Sim_epimutated_bg_list)

# Annotate with name of libraries
#write.xlsx(Sim_test_bg_results[[2]],"Sim_test_bg_results_overlap.xlsx")
Sim_test_bg_results <- lapply(Sim_test_bg_list, function(x){
  
  lapply(1:5, function(j){
    
    if(nrow(x[[j]]) > 0){
      x[[j]][1:nrow(x[[j]]), "library"] <- names(Sim_test_bg_list[[1]])[j]
      return(x[[j]])
    } 
    
  })
  
})

# then we combine them and put them in order of significance using the Old.Adjusted.P.value
Sim_test_bg_RESULTS <- lapply(Sim_test_bg_results, function(x){
  
  tempdf <- do.call(rbind, x)
  tempdf[order(tempdf$Old.Adjusted.P.value), ]
  
})

split_value_list <- list()

for(i in 1:length(Sim_test_bg_RESULTS)){
  
  split <- strsplit(Sim_test_bg_RESULTS[[i]][, 2], "/")
  
  split_value <- c()
  
  for(j in 1:length(split)){
    
    save <- split[[j]][1]
    
    split_value <- as.numeric(c(split_value, save))
  }
  
  Sim_test_bg_RESULTS[[i]]$Overlap <- split_value
  
}

# get the test_term, the test_not, the mean(sample_term) and the sample_not

# test_term is overlap
# test_not is length of gene set minus overlap

# then add 1 to all if either test_term OR sample_term = 0
# then only pursue if test_term &/OR sample_term > 4
# let the terms of interest be the terms for the enrichr result for the test list, i.e. with simultaneously epimutated genes

terms_of_interest <- Sim_test_bg_RESULTS[[1]]$Term

# and if the term is not present in the other gene list then assume there is no enrichment for that term in those lists

Test_List <- Sim_test_bg_RESULTS[[1]]
BG_List <- Sim_test_bg_RESULTS[[2]]

table_comp <- c()

for(i in 1:length(terms_of_interest)){
  
  Term_select <- terms_of_interest[[i]]
  sample_term <- 0
  
  if(Term_select %in% BG_List$Term == T){
    
    sample_term <- BG_List[which(BG_List$Term %in% Term_select), 2]}
  sample_not <- length(Sim_bg_list) - sample_term
  
  test_term <- Test_List[which(Test_List$Term %in% Term_select), 2]
  test_not <- length(Sim_test_list) - test_term
  
  # save the components in a row in a new table
  
  component_row <- data.frame(test_term, test_not, sample_term, sample_not)
  
  rownames(component_row) <- Term_select
  
  table_comp <- rbind(table_comp, component_row)
  colnames(table_comp) <- c("test_term", "test_not", "sample_term", "sample_not")
  
}

# Now we will adjust  all values by adding 1 if there are any zeros

zeros <- colSums(table_comp == 0)

if(sum(zeros) > 0){
  
  table_comp <-  table_comp + 1
}

# Now modify the table to remove any rows where test term or mean sample term are < 5
# Now apply the function only if at least one category has 5 or more genes
# It is only reasonable to assess for relative enrichment/depletion if either the test list or sample list has > 5 genes in that ontology category

table_comp <- table_comp[(table_comp$test_term < 5 & table_comp$sample_term < 5)==F, ]

# Now we have found terms where there are at least 5 genes in 1 test or sample category
# we have adjusted the data for each gene set for any 0 values relative to the specific test list
# Now we can do the analysis but only on terms present in each table_comp list

Pval_col <- c()
OR_col <- c()

for(q in 1:nrow(table_comp)){
  
  test_term <- table_comp[q, 1] 
  test_not <- table_comp[q, 2]
  mean_Sample_term <- table_comp[q, 3]
  sample_not <- table_comp[q, 4]
  
  contingency_table <-
    rbind(c(test_term, test_not),
          c(mean_Sample_term, sample_not))
  FT_out <- fisher.test(contingency_table)
  pval <- FT_out$p.value
  OR <- FT_out$estimate
  
  Pval_col <- c(Pval_col, pval)
  OR_col<- c(OR_col, OR)
  
}    

OR_pval_frame <- data.frame(Pval_col, OR_col)

colnames(OR_pval_frame) <- c("Pval", "OR")
rownames(OR_pval_frame) <- rownames(table_comp)

# Bonferroni correction

OR_pval_frame[, 1] <- p.adjust(OR_pval_frame[, 1], method="bonferroni", n = nrow(OR_pval_frame)) 

# Make plots for this

# i) In each plot order them by the largest significant enrichment out of all of the comparisons 

OR_pval_frame_ordered <- OR_pval_frame[order(OR_pval_frame[, 2], decreasing = T), ]

# Make a bubble plot 

P_1 <- as.numeric(OR_pval_frame_ordered[, 1])

OR_1 <- log10(as.numeric(OR_pval_frame_ordered[, 2]))

neg_log_p_1 <- -log(P_1)

bubble_table_1 <- data.frame(rownames(OR_pval_frame_ordered), OR_1,  P_1, neg_log_p_1)

colnames(bubble_table_1) <- c("Term", "log_OR", "p_val", "neg_log_P")

bubble_table_1$Term <- factor(bubble_table_1$Term, levels = rev(bubble_table_1$Term))

find_size <- c()
find_alpha <- c()

for(i in 1:nrow(bubble_table_1)){
  
  alpha <- 0.3
  size_point <- bubble_table_1[i,4]
  
  if(as.numeric(bubble_table_1[i, 3]) < 0.1){
    alpha <- 1
    size_point <- bubble_table_1[i, 4]
    
  }
  find_alpha <- c(find_alpha, alpha)  
  find_size <- c(find_size, size_point)
}

bubble_table <- cbind(bubble_table_1, find_alpha, find_size)

# just plot first 10

allsignifLepi<-subset(bubble_table,bubble_table$p_val<0.05)
write.xlsx(allsignifLepi,"allsignifLepi.xlsx")
bubble_table_L <- bubble_table[1:10, ]
trunc_terms_L <- as.character(bubble_table$Term)

#RNA High dose epimut>1 vs all genes
# The test gene list
Sim_test_list <- unique(tot_H_sup_1$gene)

# The background list
Sim_bg_list <- unique(C_elegans_gene_names[,4])

# Remove any genes from the background list that feature in the test list

intersect <- intersect(Sim_bg_list, Sim_test_list)
Sim_bg_list <- Sim_bg_list[Sim_bg_list %in% intersect == F]

# Send both lists to enrichr
Sim_epimutated_enriched_list  <- enrichr(Sim_test_list, chosendbs)
Sim_epimutated_bg_list  <- enrichr(Sim_bg_list, chosendbs)

Sim_test_bg_list <- list(Sim_epimutated_enriched_list, Sim_epimutated_bg_list)

# Annotate with name of libraries
#write.xlsx(Sim_test_bg_results[[2]],"Sim_test_bg_results_overlap.xlsx")
Sim_test_bg_results <- lapply(Sim_test_bg_list, function(x){
  
  lapply(1:5, function(j){
    
    if(nrow(x[[j]]) > 0){
      x[[j]][1:nrow(x[[j]]), "library"] <- names(Sim_test_bg_list[[1]])[j]
      return(x[[j]])
    } 
    
  })
  
})

# then we combine them and put them in order of significance using the Old.Adjusted.P.value
Sim_test_bg_RESULTS <- lapply(Sim_test_bg_results, function(x){
  
  tempdf <- do.call(rbind, x)
  tempdf[order(tempdf$Old.Adjusted.P.value), ]
  
})

split_value_list <- list()

for(i in 1:length(Sim_test_bg_RESULTS)){
  
  split <- strsplit(Sim_test_bg_RESULTS[[i]][, 2], "/")
  
  split_value <- c()
  
  for(j in 1:length(split)){
    
    save <- split[[j]][1]
    
    split_value <- as.numeric(c(split_value, save))
  }
  
  Sim_test_bg_RESULTS[[i]]$Overlap <- split_value
  
}

# get the test_term, the test_not, the mean(sample_term) and the sample_not

# test_term is overlap
# test_not is length of gene set minus overlap

# then add 1 to all if either test_term OR sample_term = 0
# then only pursue if test_term &/OR sample_term > 4
# let the terms of interest be the terms for the enrichr result for the test list, i.e. with simultaneously epimutated genes

terms_of_interest <- Sim_test_bg_RESULTS[[1]]$Term

# and if the term is not present in the other gene list then assume there is no enrichment for that term in those lists

Test_List <- Sim_test_bg_RESULTS[[1]]
BG_List <- Sim_test_bg_RESULTS[[2]]

table_comp <- c()

for(i in 1:length(terms_of_interest)){
  
  Term_select <- terms_of_interest[[i]]
  sample_term <- 0
  
  if(Term_select %in% BG_List$Term == T){
    
    sample_term <- BG_List[which(BG_List$Term %in% Term_select), 2]}
  sample_not <- length(Sim_bg_list) - sample_term
  
  test_term <- Test_List[which(Test_List$Term %in% Term_select), 2]
  test_not <- length(Sim_test_list) - test_term
  
  # save the components in a row in a new table
  
  component_row <- data.frame(test_term, test_not, sample_term, sample_not)
  
  rownames(component_row) <- Term_select
  
  table_comp <- rbind(table_comp, component_row)
  colnames(table_comp) <- c("test_term", "test_not", "sample_term", "sample_not")
  
}

# Now we will adjust  all values by adding 1 if there are any zeros

zeros <- colSums(table_comp == 0)

if(sum(zeros) > 0){
  
  table_comp <-  table_comp + 1
}

# Now modify the table to remove any rows where test term or mean sample term are < 5
# Now apply the function only if at least one category has 5 or more genes
# It is only reasonable to assess for relative enrichment/depletion if either the test list or sample list has > 5 genes in that ontology category

table_comp <- table_comp[(table_comp$test_term < 5 & table_comp$sample_term < 5)==F, ]

# Now we have found terms where there are at least 5 genes in 1 test or sample category
# we have adjusted the data for each gene set for any 0 values relative to the specific test list
# Now we can do the analysis but only on terms present in each table_comp list

Pval_col <- c()
OR_col <- c()

for(q in 1:nrow(table_comp)){
  
  test_term <- table_comp[q, 1] 
  test_not <- table_comp[q, 2]
  mean_Sample_term <- table_comp[q, 3]
  sample_not <- table_comp[q, 4]
  
  contingency_table <-
    rbind(c(test_term, test_not),
          c(mean_Sample_term, sample_not))
  FT_out <- fisher.test(contingency_table)
  pval <- FT_out$p.value
  OR <- FT_out$estimate
  
  Pval_col <- c(Pval_col, pval)
  OR_col<- c(OR_col, OR)
  
}    

OR_pval_frame <- data.frame(Pval_col, OR_col)

colnames(OR_pval_frame) <- c("Pval", "OR")
rownames(OR_pval_frame) <- rownames(table_comp)

# Bonferroni correction

OR_pval_frame[, 1] <- p.adjust(OR_pval_frame[, 1], method="bonferroni", n = nrow(OR_pval_frame)) 

# Make plots for this

# i) In each plot order them by the largest significant enrichment out of all of the comparisons 

OR_pval_frame_ordered <- OR_pval_frame[order(OR_pval_frame[, 2], decreasing = T), ]

# Make a bubble plot 

P_1 <- as.numeric(OR_pval_frame_ordered[, 1])

OR_1 <- log10(as.numeric(OR_pval_frame_ordered[, 2]))

neg_log_p_1 <- -log(P_1)

bubble_table_1 <- data.frame(rownames(OR_pval_frame_ordered), OR_1,  P_1, neg_log_p_1)

colnames(bubble_table_1) <- c("Term", "log_OR", "p_val", "neg_log_P")

bubble_table_1$Term <- factor(bubble_table_1$Term, levels = rev(bubble_table_1$Term))

find_size <- c()
find_alpha <- c()

for(i in 1:nrow(bubble_table_1)){
  
  alpha <- 0.3
  size_point <- bubble_table_1[i,4]
  
  if(as.numeric(bubble_table_1[i, 3]) < 0.1){
    alpha <- 1
    size_point <- bubble_table_1[i, 4]
    
  }
  find_alpha <- c(find_alpha, alpha)  
  find_size <- c(find_size, size_point)
}

bubble_table <- cbind(bubble_table_1, find_alpha, find_size)

# just plot first 10

allsignifHepi<-subset(bubble_table,bubble_table$p_val<0.05)
write.xlsx(allsignifHepi,"allsignifHepi.xlsx")
bubble_table_H <- bubble_table[1:10, ]
trunc_terms_H <- as.character(bubble_table$Term)


#Plot all conditions
bubble_table_C<-bubble_table_C  %>%
  mutate(Condition = "Control")
bubble_table_L<-bubble_table_L %>%
  mutate(Condition = "Low dose")
bubble_table_H<-bubble_table_H %>%
  mutate(Condition = "High dose")

df_list <- list(bubble_table_C, bubble_table_L, bubble_table_H)
Bubble_table_all<-Reduce(function(x, y) merge(x, y, all=TRUE), df_list, accumulate=FALSE)

# Trunc terms manually entered as Y axis labels in Adobe Illustrator

trunc_terms <- c(
  
  "detection of chemical stimulus involved \n in sensory perception (GO:0050907)",
  "sensory perception of \n chemical stimulus (GO:0007606)" ,                                               
  "integral component of \n plasma membrane (GO:0005887)" ,                
  "nucleus (GO:0005634)",                                  
  "olfactory receptor activity (GO:0004984)",                            
  "G-protein coupled receptor \n activity (GO:0004930)" ,              
  "G-protein coupled olfactory \n receptor activity (GO:0038022)",                                       
  "protein kinase activity (GO:0004672)",  
  "protein serine/threonine \n kinase activity (GO:0004674)",
  "protein phosphorylation (GO:0006468)",
  "peptidyl-serine modification (GO:0018209)"
)

guides_merge <- function(gdefs) {
  gdefs <- lapply(gdefs, function(g) { g$hash <- paste(g$order, g$hash, sep = "z"); g})
  tapply(gdefs, sapply(gdefs, function(g)g$hash), function(gs)Reduce(guide_merge, gs))
}
environment(guides_merge) <- environment(ggplot)
assignInNamespace("guides_merge", guides_merge, pos = "package:ggplot2")

Bubble_table_all_plot <-
  
  ggplot(Bubble_table_all, aes(y=Term, x=as.numeric(log_OR)))+
  geom_point(aes(color=Condition, size=as.numeric(find_size), alpha=find_alpha))+
  scale_x_continuous(breaks = c(1.1,1.2,1.3), limits = c(1.1, 1.3)) +
  geom_vline(xintercept = 0, colour = "grey")+
  theme_bw()+
  
  
  scale_size_continuous(range = c(2, 15), breaks = c(25,50,100), 
                        limits = c(25, 135))+  
  
  scale_alpha(name = paste("Condition"), range = c(0.3, 1),
              breaks = seq(1, 0.3, length = 3),
              limits = c(0.29, 1.05), 
              labels = c("Control", "Low dose", "High dose"))

Bubble_table_all_plot <- Bubble_table_all_plot + labs(y="Gene Ontology Terms\n", x = "log10(Odds Ratio for enrichment)")+
  theme(axis.title.x=element_text(face = "bold", size=18))+
  theme(axis.title.y=element_text(face = "bold", size=18))+
  theme(plot.title = element_text(face = "bold", size=16))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_discrete(labels= rev(trunc_terms))+
  scale_colour_manual(values = c("cornflowerblue", "red","darkgreen"), guide=F) +
  theme(axis.text.x = element_text(color="#000000", size=16))+
  theme(axis.text.y = element_text(color="#000000", size = 16))+
  guides(size = guide_legend(order = 3), colour = "none", alpha = guide_legend(order = 2))+
  theme(legend.text=element_text(color="#000000", size=16))

Bubble_table_all_plot

#Get raw data
write.xlsx(Bubble_table_all,"Data_Fig_3_D.xlsx")

#--------------------------
#####Sup.Fig.4.D - GO term enrichment analysis for RNA epimutations####
tot_C1_sup_1<-subset(RNA_epimut_C,RNA_epimut_C$Lineage=="C1" & RNA_epimut_C$length>1)
tot_C2_sup_1<-subset(RNA_epimut_C,RNA_epimut_C$Lineage=="C2" & RNA_epimut_C$length>1)
tot_L1_sup_1<-subset(RNA_epimut_L,RNA_epimut_L$Lineage=="L1" & RNA_epimut_L$length>1)
tot_L2_sup_1<-subset(RNA_epimut_L,RNA_epimut_L$Lineage=="L2" & RNA_epimut_L$length>1)
tot_H1_sup_1<-subset(RNA_epimut_H,RNA_epimut_H$Lineage=="H1" & RNA_epimut_H$length>1)
tot_H2_sup_1<-subset(RNA_epimut_H,RNA_epimut_H$Lineage=="H2" & RNA_epimut_H$length>1)
#RNA Control 1 epimut>1 vs all genes
# The test gene list
Sim_test_list <- unique(tot_C1_sup_1$gene)

# The background list
Sim_bg_list <- unique(C_elegans_gene_names[,4])

# Remove any genes from the background list that feature in the test list

intersect <- intersect(Sim_bg_list, Sim_test_list)
Sim_bg_list <- Sim_bg_list[Sim_bg_list %in% intersect == F]

# Send both lists to enrichr
Sim_epimutated_enriched_list  <- enrichr(Sim_test_list, chosendbs)
Sim_epimutated_bg_list  <- enrichr(Sim_bg_list, chosendbs)

Sim_test_bg_list <- list(Sim_epimutated_enriched_list, Sim_epimutated_bg_list)

# Annotate with name of libraries
#write.xlsx(Sim_test_bg_results[[2]],"Sim_test_bg_results_overlap.xlsx")
Sim_test_bg_results <- lapply(Sim_test_bg_list, function(x){
  
  lapply(1:5, function(j){
    
    if(nrow(x[[j]]) > 0){
      x[[j]][1:nrow(x[[j]]), "library"] <- names(Sim_test_bg_list[[1]])[j]
      return(x[[j]])
    } 
    
  })
  
})

# then we combine them and put them in order of significance using the Old.Adjusted.P.value
Sim_test_bg_RESULTS <- lapply(Sim_test_bg_results, function(x){
  
  tempdf <- do.call(rbind, x)
  tempdf[order(tempdf$Old.Adjusted.P.value), ]
  
})

split_value_list <- list()

for(i in 1:length(Sim_test_bg_RESULTS)){
  
  split <- strsplit(Sim_test_bg_RESULTS[[i]][, 2], "/")
  
  split_value <- c()
  
  for(j in 1:length(split)){
    
    save <- split[[j]][1]
    
    split_value <- as.numeric(c(split_value, save))
  }
  
  Sim_test_bg_RESULTS[[i]]$Overlap <- split_value
  
}

# get the test_term, the test_not, the mean(sample_term) and the sample_not

# test_term is overlap
# test_not is length of gene set minus overlap

# then add 1 to all if either test_term OR sample_term = 0
# then only pursue if test_term &/OR sample_term > 4
# let the terms of interest be the terms for the enrichr result for the test list, i.e. with simultaneously epimutated genes

terms_of_interest <- Sim_test_bg_RESULTS[[1]]$Term

# and if the term is not present in the other gene list then assume there is no enrichment for that term in those lists

Test_List <- Sim_test_bg_RESULTS[[1]]
BG_List <- Sim_test_bg_RESULTS[[2]]

table_comp <- c()

for(i in 1:length(terms_of_interest)){
  
  Term_select <- terms_of_interest[[i]]
  sample_term <- 0
  
  if(Term_select %in% BG_List$Term == T){
    
    sample_term <- BG_List[which(BG_List$Term %in% Term_select), 2]}
  sample_not <- length(Sim_bg_list) - sample_term
  
  test_term <- Test_List[which(Test_List$Term %in% Term_select), 2]
  test_not <- length(Sim_test_list) - test_term
  
  # save the components in a row in a new table
  
  component_row <- data.frame(test_term, test_not, sample_term, sample_not)
  
  rownames(component_row) <- Term_select
  
  table_comp <- rbind(table_comp, component_row)
  colnames(table_comp) <- c("test_term", "test_not", "sample_term", "sample_not")
  
}

# Now we will adjust  all values by adding 1 if there are any zeros

zeros <- colSums(table_comp == 0)

if(sum(zeros) > 0){
  
  table_comp <-  table_comp + 1
}

# Now modify the table to remove any rows where test term or mean sample term are < 5
# Now apply the function only if at least one category has 5 or more genes
# It is only reasonable to assess for relative enrichment/depletion if either the test list or sample list has > 5 genes in that ontology category

table_comp <- table_comp[(table_comp$test_term < 5 & table_comp$sample_term < 5)==F, ]

# Now we have found terms where there are at least 5 genes in 1 test or sample category
# we have adjusted the data for each gene set for any 0 values relative to the specific test list
# Now we can do the analysis but only on terms present in each table_comp list

Pval_col <- c()
OR_col <- c()

for(q in 1:nrow(table_comp)){
  
  test_term <- table_comp[q, 1] 
  test_not <- table_comp[q, 2]
  mean_Sample_term <- table_comp[q, 3]
  sample_not <- table_comp[q, 4]
  
  contingency_table <-
    rbind(c(test_term, test_not),
          c(mean_Sample_term, sample_not))
  FT_out <- fisher.test(contingency_table)
  pval <- FT_out$p.value
  OR <- FT_out$estimate
  
  Pval_col <- c(Pval_col, pval)
  OR_col<- c(OR_col, OR)
  
}    

OR_pval_frame <- data.frame(Pval_col, OR_col)

colnames(OR_pval_frame) <- c("Pval", "OR")
rownames(OR_pval_frame) <- rownames(table_comp)

# Bonferroni correction

OR_pval_frame[, 1] <- p.adjust(OR_pval_frame[, 1], method="bonferroni", n = nrow(OR_pval_frame)) 

# Make plots for this

# i) In each plot order them by the largest significant enrichment out of all of the comparisons 

OR_pval_frame_ordered <- OR_pval_frame[order(OR_pval_frame[, 2], decreasing = T), ]

# Make a bubble plot 

P_1 <- as.numeric(OR_pval_frame_ordered[, 1])

OR_1 <- log10(as.numeric(OR_pval_frame_ordered[, 2]))

neg_log_p_1 <- -log(P_1)

bubble_table_1 <- data.frame(rownames(OR_pval_frame_ordered), OR_1,  P_1, neg_log_p_1)

colnames(bubble_table_1) <- c("Term", "log_OR", "p_val", "neg_log_P")

bubble_table_1$Term <- factor(bubble_table_1$Term, levels = rev(bubble_table_1$Term))

find_size <- c()
find_alpha <- c()

for(i in 1:nrow(bubble_table_1)){
  
  alpha <- 0.3
  size_point <- bubble_table_1[i,4]
  
  if(as.numeric(bubble_table_1[i, 3]) < 0.1){
    alpha <- 1
    size_point <- bubble_table_1[i, 4]
    
  }
  find_alpha <- c(find_alpha, alpha)  
  find_size <- c(find_size, size_point)
}

bubble_table <- cbind(bubble_table_1, find_alpha, find_size)

# just plot first 10

allsignifC1epi<-subset(bubble_table,bubble_table$p_val<0.05)
write.xlsx(allsignifC1epi,"allsignifC1epi.xlsx")
bubble_table_C1 <- bubble_table[1:10, ]
trunc_terms_C1 <- as.character(bubble_table$Term)

#RNA Control 2 epimut>1 vs all genes
# The test gene list
Sim_test_list <- unique(tot_C2_sup_1$gene)

# The background list
Sim_bg_list <- unique(C_elegans_gene_names[,4])

# Remove any genes from the background list that feature in the test list

intersect <- intersect(Sim_bg_list, Sim_test_list)
Sim_bg_list <- Sim_bg_list[Sim_bg_list %in% intersect == F]

# Send both lists to enrichr
Sim_epimutated_enriched_list  <- enrichr(Sim_test_list, chosendbs)
Sim_epimutated_bg_list  <- enrichr(Sim_bg_list, chosendbs)

Sim_test_bg_list <- list(Sim_epimutated_enriched_list, Sim_epimutated_bg_list)

# Annotate with name of libraries
#write.xlsx(Sim_test_bg_results[[2]],"Sim_test_bg_results_overlap.xlsx")
Sim_test_bg_results <- lapply(Sim_test_bg_list, function(x){
  
  lapply(1:5, function(j){
    
    if(nrow(x[[j]]) > 0){
      x[[j]][1:nrow(x[[j]]), "library"] <- names(Sim_test_bg_list[[1]])[j]
      return(x[[j]])
    } 
    
  })
  
})

# then we combine them and put them in order of significance using the Old.Adjusted.P.value
Sim_test_bg_RESULTS <- lapply(Sim_test_bg_results, function(x){
  
  tempdf <- do.call(rbind, x)
  tempdf[order(tempdf$Old.Adjusted.P.value), ]
  
})

split_value_list <- list()

for(i in 1:length(Sim_test_bg_RESULTS)){
  
  split <- strsplit(Sim_test_bg_RESULTS[[i]][, 2], "/")
  
  split_value <- c()
  
  for(j in 1:length(split)){
    
    save <- split[[j]][1]
    
    split_value <- as.numeric(c(split_value, save))
  }
  
  Sim_test_bg_RESULTS[[i]]$Overlap <- split_value
  
}

# get the test_term, the test_not, the mean(sample_term) and the sample_not

# test_term is overlap
# test_not is length of gene set minus overlap

# then add 1 to all if either test_term OR sample_term = 0
# then only pursue if test_term &/OR sample_term > 4
# let the terms of interest be the terms for the enrichr result for the test list, i.e. with simultaneously epimutated genes

terms_of_interest <- Sim_test_bg_RESULTS[[1]]$Term

# and if the term is not present in the other gene list then assume there is no enrichment for that term in those lists

Test_List <- Sim_test_bg_RESULTS[[1]]
BG_List <- Sim_test_bg_RESULTS[[2]]

table_comp <- c()

for(i in 1:length(terms_of_interest)){
  
  Term_select <- terms_of_interest[[i]]
  sample_term <- 0
  
  if(Term_select %in% BG_List$Term == T){
    
    sample_term <- BG_List[which(BG_List$Term %in% Term_select), 2]}
  sample_not <- length(Sim_bg_list) - sample_term
  
  test_term <- Test_List[which(Test_List$Term %in% Term_select), 2]
  test_not <- length(Sim_test_list) - test_term
  
  # save the components in a row in a new table
  
  component_row <- data.frame(test_term, test_not, sample_term, sample_not)
  
  rownames(component_row) <- Term_select
  
  table_comp <- rbind(table_comp, component_row)
  colnames(table_comp) <- c("test_term", "test_not", "sample_term", "sample_not")
  
}

# Now we will adjust  all values by adding 1 if there are any zeros

zeros <- colSums(table_comp == 0)

if(sum(zeros) > 0){
  
  table_comp <-  table_comp + 1
}

# Now modify the table to remove any rows where test term or mean sample term are < 5
# Now apply the function only if at least one category has 5 or more genes
# It is only reasonable to assess for relative enrichment/depletion if either the test list or sample list has > 5 genes in that ontology category

table_comp <- table_comp[(table_comp$test_term < 5 & table_comp$sample_term < 5)==F, ]

# Now we have found terms where there are at least 5 genes in 1 test or sample category
# we have adjusted the data for each gene set for any 0 values relative to the specific test list
# Now we can do the analysis but only on terms present in each table_comp list

Pval_col <- c()
OR_col <- c()

for(q in 1:nrow(table_comp)){
  
  test_term <- table_comp[q, 1] 
  test_not <- table_comp[q, 2]
  mean_Sample_term <- table_comp[q, 3]
  sample_not <- table_comp[q, 4]
  
  contingency_table <-
    rbind(c(test_term, test_not),
          c(mean_Sample_term, sample_not))
  FT_out <- fisher.test(contingency_table)
  pval <- FT_out$p.value
  OR <- FT_out$estimate
  
  Pval_col <- c(Pval_col, pval)
  OR_col<- c(OR_col, OR)
  
}    

OR_pval_frame <- data.frame(Pval_col, OR_col)

colnames(OR_pval_frame) <- c("Pval", "OR")
rownames(OR_pval_frame) <- rownames(table_comp)

# Bonferroni correction

OR_pval_frame[, 1] <- p.adjust(OR_pval_frame[, 1], method="bonferroni", n = nrow(OR_pval_frame)) 

# Make plots for this

# i) In each plot order them by the largest significant enrichment out of all of the comparisons 

OR_pval_frame_ordered <- OR_pval_frame[order(OR_pval_frame[, 2], decreasing = T), ]

# Make a bubble plot 

P_1 <- as.numeric(OR_pval_frame_ordered[, 1])

OR_1 <- log10(as.numeric(OR_pval_frame_ordered[, 2]))

neg_log_p_1 <- -log(P_1)

bubble_table_1 <- data.frame(rownames(OR_pval_frame_ordered), OR_1,  P_1, neg_log_p_1)

colnames(bubble_table_1) <- c("Term", "log_OR", "p_val", "neg_log_P")

bubble_table_1$Term <- factor(bubble_table_1$Term, levels = rev(bubble_table_1$Term))

find_size <- c()
find_alpha <- c()

for(i in 1:nrow(bubble_table_1)){
  
  alpha <- 0.3
  size_point <- bubble_table_1[i,4]
  
  if(as.numeric(bubble_table_1[i, 3]) < 0.1){
    alpha <- 1
    size_point <- bubble_table_1[i, 4]
    
  }
  find_alpha <- c(find_alpha, alpha)  
  find_size <- c(find_size, size_point)
}

bubble_table <- cbind(bubble_table_1, find_alpha, find_size)

# just plot first 10

allsignifC2epi<-subset(bubble_table,bubble_table$p_val<0.05)
write.xlsx(allsignifC2epi,"allsignifC2epi.xlsx")
bubble_table_C2 <- bubble_table[1:10, ]
trunc_terms_C2 <- as.character(bubble_table$Term)

#RNA Low dose 1 epimut>1 vs all genes
# The test gene list
Sim_test_list <- unique(tot_L1_sup_1$gene)

# The background list
Sim_bg_list <- unique(C_elegans_gene_names[,4])

# Remove any genes from the background list that feature in the test list

intersect <- intersect(Sim_bg_list, Sim_test_list)
Sim_bg_list <- Sim_bg_list[Sim_bg_list %in% intersect == F]

# Send both lists to enrichr
Sim_epimutated_enriched_list  <- enrichr(Sim_test_list, chosendbs)
Sim_epimutated_bg_list  <- enrichr(Sim_bg_list, chosendbs)

Sim_test_bg_list <- list(Sim_epimutated_enriched_list, Sim_epimutated_bg_list)

# Annotate with name of libraries
#write.xlsx(Sim_test_bg_results[[2]],"Sim_test_bg_results_overlap.xlsx")
Sim_test_bg_results <- lapply(Sim_test_bg_list, function(x){
  
  lapply(1:5, function(j){
    
    if(nrow(x[[j]]) > 0){
      x[[j]][1:nrow(x[[j]]), "library"] <- names(Sim_test_bg_list[[1]])[j]
      return(x[[j]])
    } 
    
  })
  
})

# then we combine them and put them in order of significance using the Old.Adjusted.P.value
Sim_test_bg_RESULTS <- lapply(Sim_test_bg_results, function(x){
  
  tempdf <- do.call(rbind, x)
  tempdf[order(tempdf$Old.Adjusted.P.value), ]
  
})

split_value_list <- list()

for(i in 1:length(Sim_test_bg_RESULTS)){
  
  split <- strsplit(Sim_test_bg_RESULTS[[i]][, 2], "/")
  
  split_value <- c()
  
  for(j in 1:length(split)){
    
    save <- split[[j]][1]
    
    split_value <- as.numeric(c(split_value, save))
  }
  
  Sim_test_bg_RESULTS[[i]]$Overlap <- split_value
  
}

# get the test_term, the test_not, the mean(sample_term) and the sample_not

# test_term is overlap
# test_not is length of gene set minus overlap

# then add 1 to all if either test_term OR sample_term = 0
# then only pursue if test_term &/OR sample_term > 4
# let the terms of interest be the terms for the enrichr result for the test list, i.e. with simultaneously epimutated genes

terms_of_interest <- Sim_test_bg_RESULTS[[1]]$Term

# and if the term is not present in the other gene list then assume there is no enrichment for that term in those lists

Test_List <- Sim_test_bg_RESULTS[[1]]
BG_List <- Sim_test_bg_RESULTS[[2]]

table_comp <- c()

for(i in 1:length(terms_of_interest)){
  
  Term_select <- terms_of_interest[[i]]
  sample_term <- 0
  
  if(Term_select %in% BG_List$Term == T){
    
    sample_term <- BG_List[which(BG_List$Term %in% Term_select), 2]}
  sample_not <- length(Sim_bg_list) - sample_term
  
  test_term <- Test_List[which(Test_List$Term %in% Term_select), 2]
  test_not <- length(Sim_test_list) - test_term
  
  # save the components in a row in a new table
  
  component_row <- data.frame(test_term, test_not, sample_term, sample_not)
  
  rownames(component_row) <- Term_select
  
  table_comp <- rbind(table_comp, component_row)
  colnames(table_comp) <- c("test_term", "test_not", "sample_term", "sample_not")
  
}

# Now we will adjust  all values by adding 1 if there are any zeros

zeros <- colSums(table_comp == 0)

if(sum(zeros) > 0){
  
  table_comp <-  table_comp + 1
}

# Now modify the table to remove any rows where test term or mean sample term are < 5
# Now apply the function only if at least one category has 5 or more genes
# It is only reasonable to assess for relative enrichment/depletion if either the test list or sample list has > 5 genes in that ontology category

table_comp <- table_comp[(table_comp$test_term < 5 & table_comp$sample_term < 5)==F, ]

# Now we have found terms where there are at least 5 genes in 1 test or sample category
# we have adjusted the data for each gene set for any 0 values relative to the specific test list
# Now we can do the analysis but only on terms present in each table_comp list

Pval_col <- c()
OR_col <- c()

for(q in 1:nrow(table_comp)){
  
  test_term <- table_comp[q, 1] 
  test_not <- table_comp[q, 2]
  mean_Sample_term <- table_comp[q, 3]
  sample_not <- table_comp[q, 4]
  
  contingency_table <-
    rbind(c(test_term, test_not),
          c(mean_Sample_term, sample_not))
  FT_out <- fisher.test(contingency_table)
  pval <- FT_out$p.value
  OR <- FT_out$estimate
  
  Pval_col <- c(Pval_col, pval)
  OR_col<- c(OR_col, OR)
  
}    

OR_pval_frame <- data.frame(Pval_col, OR_col)

colnames(OR_pval_frame) <- c("Pval", "OR")
rownames(OR_pval_frame) <- rownames(table_comp)

# Bonferroni correction

OR_pval_frame[, 1] <- p.adjust(OR_pval_frame[, 1], method="bonferroni", n = nrow(OR_pval_frame)) 

# Make plots for this

# i) In each plot order them by the largest significant enrichment out of all of the comparisons 

OR_pval_frame_ordered <- OR_pval_frame[order(OR_pval_frame[, 2], decreasing = T), ]

# Make a bubble plot 

P_1 <- as.numeric(OR_pval_frame_ordered[, 1])

OR_1 <- log10(as.numeric(OR_pval_frame_ordered[, 2]))

neg_log_p_1 <- -log(P_1)

bubble_table_1 <- data.frame(rownames(OR_pval_frame_ordered), OR_1,  P_1, neg_log_p_1)

colnames(bubble_table_1) <- c("Term", "log_OR", "p_val", "neg_log_P")

bubble_table_1$Term <- factor(bubble_table_1$Term, levels = rev(bubble_table_1$Term))

find_size <- c()
find_alpha <- c()

for(i in 1:nrow(bubble_table_1)){
  
  alpha <- 0.3
  size_point <- bubble_table_1[i,4]
  
  if(as.numeric(bubble_table_1[i, 3]) < 0.1){
    alpha <- 1
    size_point <- bubble_table_1[i, 4]
    
  }
  find_alpha <- c(find_alpha, alpha)  
  find_size <- c(find_size, size_point)
}

bubble_table <- cbind(bubble_table_1, find_alpha, find_size)

# just plot first 10

allsignifL1epi<-subset(bubble_table,bubble_table$p_val<0.05)
write.xlsx(allsignifL1epi,"allsignifL1epi.xlsx")
bubble_table_L1 <- bubble_table[1:10, ]
trunc_terms_L1 <- as.character(bubble_table$Term)

#RNA Low dose 2 epimut>1 vs all genes
# The test gene list
Sim_test_list <- unique(tot_L2_sup_1$gene)

# The background list
Sim_bg_list <- unique(C_elegans_gene_names[,4])

# Remove any genes from the background list that feature in the test list

intersect <- intersect(Sim_bg_list, Sim_test_list)
Sim_bg_list <- Sim_bg_list[Sim_bg_list %in% intersect == F]

# Send both lists to enrichr
Sim_epimutated_enriched_list  <- enrichr(Sim_test_list, chosendbs)
Sim_epimutated_bg_list  <- enrichr(Sim_bg_list, chosendbs)

Sim_test_bg_list <- list(Sim_epimutated_enriched_list, Sim_epimutated_bg_list)

# Annotate with name of libraries
#write.xlsx(Sim_test_bg_results[[2]],"Sim_test_bg_results_overlap.xlsx")
Sim_test_bg_results <- lapply(Sim_test_bg_list, function(x){
  
  lapply(1:5, function(j){
    
    if(nrow(x[[j]]) > 0){
      x[[j]][1:nrow(x[[j]]), "library"] <- names(Sim_test_bg_list[[1]])[j]
      return(x[[j]])
    } 
    
  })
  
})

# then we combine them and put them in order of significance using the Old.Adjusted.P.value
Sim_test_bg_RESULTS <- lapply(Sim_test_bg_results, function(x){
  
  tempdf <- do.call(rbind, x)
  tempdf[order(tempdf$Old.Adjusted.P.value), ]
  
})

split_value_list <- list()

for(i in 1:length(Sim_test_bg_RESULTS)){
  
  split <- strsplit(Sim_test_bg_RESULTS[[i]][, 2], "/")
  
  split_value <- c()
  
  for(j in 1:length(split)){
    
    save <- split[[j]][1]
    
    split_value <- as.numeric(c(split_value, save))
  }
  
  Sim_test_bg_RESULTS[[i]]$Overlap <- split_value
  
}

# get the test_term, the test_not, the mean(sample_term) and the sample_not

# test_term is overlap
# test_not is length of gene set minus overlap

# then add 1 to all if either test_term OR sample_term = 0
# then only pursue if test_term &/OR sample_term > 4
# let the terms of interest be the terms for the enrichr result for the test list, i.e. with simultaneously epimutated genes

terms_of_interest <- Sim_test_bg_RESULTS[[1]]$Term

# and if the term is not present in the other gene list then assume there is no enrichment for that term in those lists

Test_List <- Sim_test_bg_RESULTS[[1]]
BG_List <- Sim_test_bg_RESULTS[[2]]

table_comp <- c()

for(i in 1:length(terms_of_interest)){
  
  Term_select <- terms_of_interest[[i]]
  sample_term <- 0
  
  if(Term_select %in% BG_List$Term == T){
    
    sample_term <- BG_List[which(BG_List$Term %in% Term_select), 2]}
  sample_not <- length(Sim_bg_list) - sample_term
  
  test_term <- Test_List[which(Test_List$Term %in% Term_select), 2]
  test_not <- length(Sim_test_list) - test_term
  
  # save the components in a row in a new table
  
  component_row <- data.frame(test_term, test_not, sample_term, sample_not)
  
  rownames(component_row) <- Term_select
  
  table_comp <- rbind(table_comp, component_row)
  colnames(table_comp) <- c("test_term", "test_not", "sample_term", "sample_not")
  
}

# Now we will adjust  all values by adding 1 if there are any zeros

zeros <- colSums(table_comp == 0)

if(sum(zeros) > 0){
  
  table_comp <-  table_comp + 1
}

# Now modify the table to remove any rows where test term or mean sample term are < 5
# Now apply the function only if at least one category has 5 or more genes
# It is only reasonable to assess for relative enrichment/depletion if either the test list or sample list has > 5 genes in that ontology category

table_comp <- table_comp[(table_comp$test_term < 5 & table_comp$sample_term < 5)==F, ]

# Now we have found terms where there are at least 5 genes in 1 test or sample category
# we have adjusted the data for each gene set for any 0 values relative to the specific test list
# Now we can do the analysis but only on terms present in each table_comp list

Pval_col <- c()
OR_col <- c()

for(q in 1:nrow(table_comp)){
  
  test_term <- table_comp[q, 1] 
  test_not <- table_comp[q, 2]
  mean_Sample_term <- table_comp[q, 3]
  sample_not <- table_comp[q, 4]
  
  contingency_table <-
    rbind(c(test_term, test_not),
          c(mean_Sample_term, sample_not))
  FT_out <- fisher.test(contingency_table)
  pval <- FT_out$p.value
  OR <- FT_out$estimate
  
  Pval_col <- c(Pval_col, pval)
  OR_col<- c(OR_col, OR)
  
}    

OR_pval_frame <- data.frame(Pval_col, OR_col)

colnames(OR_pval_frame) <- c("Pval", "OR")
rownames(OR_pval_frame) <- rownames(table_comp)

# Bonferroni correction

OR_pval_frame[, 1] <- p.adjust(OR_pval_frame[, 1], method="bonferroni", n = nrow(OR_pval_frame)) 

# Make plots for this

# i) In each plot order them by the largest significant enrichment out of all of the comparisons 

OR_pval_frame_ordered <- OR_pval_frame[order(OR_pval_frame[, 2], decreasing = T), ]

# Make a bubble plot 

P_1 <- as.numeric(OR_pval_frame_ordered[, 1])

OR_1 <- log10(as.numeric(OR_pval_frame_ordered[, 2]))

neg_log_p_1 <- -log(P_1)

bubble_table_1 <- data.frame(rownames(OR_pval_frame_ordered), OR_1,  P_1, neg_log_p_1)

colnames(bubble_table_1) <- c("Term", "log_OR", "p_val", "neg_log_P")

bubble_table_1$Term <- factor(bubble_table_1$Term, levels = rev(bubble_table_1$Term))

find_size <- c()
find_alpha <- c()

for(i in 1:nrow(bubble_table_1)){
  
  alpha <- 0.3
  size_point <- bubble_table_1[i,4]
  
  if(as.numeric(bubble_table_1[i, 3]) < 0.1){
    alpha <- 1
    size_point <- bubble_table_1[i, 4]
    
  }
  find_alpha <- c(find_alpha, alpha)  
  find_size <- c(find_size, size_point)
}

bubble_table <- cbind(bubble_table_1, find_alpha, find_size)

# just plot first 10

allsignifL2epi<-subset(bubble_table,bubble_table$p_val<0.05)
write.xlsx(allsignifL2epi,"allsignifL2epi.xlsx")
bubble_table_L2 <- bubble_table[1:10, ]
trunc_terms_L2 <- as.character(bubble_table$Term)

#RNA High dose 1 epimut>1 vs all genes
# The test gene list
Sim_test_list <- unique(tot_H1_sup_1$gene)

# The background list
Sim_bg_list <- unique(C_elegans_gene_names[,4])

# Remove any genes from the background list that feature in the test list

intersect <- intersect(Sim_bg_list, Sim_test_list)
Sim_bg_list <- Sim_bg_list[Sim_bg_list %in% intersect == F]

# Send both lists to enrichr
Sim_epimutated_enriched_list  <- enrichr(Sim_test_list, chosendbs)
Sim_epimutated_bg_list  <- enrichr(Sim_bg_list, chosendbs)

Sim_test_bg_list <- list(Sim_epimutated_enriched_list, Sim_epimutated_bg_list)

# Annotate with name of libraries
#write.xlsx(Sim_test_bg_results[[2]],"Sim_test_bg_results_overlap.xlsx")
Sim_test_bg_results <- lapply(Sim_test_bg_list, function(x){
  
  lapply(1:5, function(j){
    
    if(nrow(x[[j]]) > 0){
      x[[j]][1:nrow(x[[j]]), "library"] <- names(Sim_test_bg_list[[1]])[j]
      return(x[[j]])
    } 
    
  })
  
})

# then we combine them and put them in order of significance using the Old.Adjusted.P.value
Sim_test_bg_RESULTS <- lapply(Sim_test_bg_results, function(x){
  
  tempdf <- do.call(rbind, x)
  tempdf[order(tempdf$Old.Adjusted.P.value), ]
  
})

split_value_list <- list()

for(i in 1:length(Sim_test_bg_RESULTS)){
  
  split <- strsplit(Sim_test_bg_RESULTS[[i]][, 2], "/")
  
  split_value <- c()
  
  for(j in 1:length(split)){
    
    save <- split[[j]][1]
    
    split_value <- as.numeric(c(split_value, save))
  }
  
  Sim_test_bg_RESULTS[[i]]$Overlap <- split_value
  
}

# get the test_term, the test_not, the mean(sample_term) and the sample_not

# test_term is overlap
# test_not is length of gene set minus overlap

# then add 1 to all if either test_term OR sample_term = 0
# then only pursue if test_term &/OR sample_term > 4
# let the terms of interest be the terms for the enrichr result for the test list, i.e. with simultaneously epimutated genes

terms_of_interest <- Sim_test_bg_RESULTS[[1]]$Term

# and if the term is not present in the other gene list then assume there is no enrichment for that term in those lists

Test_List <- Sim_test_bg_RESULTS[[1]]
BG_List <- Sim_test_bg_RESULTS[[2]]

table_comp <- c()

for(i in 1:length(terms_of_interest)){
  
  Term_select <- terms_of_interest[[i]]
  sample_term <- 0
  
  if(Term_select %in% BG_List$Term == T){
    
    sample_term <- BG_List[which(BG_List$Term %in% Term_select), 2]}
  sample_not <- length(Sim_bg_list) - sample_term
  
  test_term <- Test_List[which(Test_List$Term %in% Term_select), 2]
  test_not <- length(Sim_test_list) - test_term
  
  # save the components in a row in a new table
  
  component_row <- data.frame(test_term, test_not, sample_term, sample_not)
  
  rownames(component_row) <- Term_select
  
  table_comp <- rbind(table_comp, component_row)
  colnames(table_comp) <- c("test_term", "test_not", "sample_term", "sample_not")
  
}

# Now we will adjust  all values by adding 1 if there are any zeros

zeros <- colSums(table_comp == 0)

if(sum(zeros) > 0){
  
  table_comp <-  table_comp + 1
}

# Now modify the table to remove any rows where test term or mean sample term are < 5
# Now apply the function only if at least one category has 5 or more genes
# It is only reasonable to assess for relative enrichment/depletion if either the test list or sample list has > 5 genes in that ontology category

table_comp <- table_comp[(table_comp$test_term < 5 & table_comp$sample_term < 5)==F, ]

# Now we have found terms where there are at least 5 genes in 1 test or sample category
# we have adjusted the data for each gene set for any 0 values relative to the specific test list
# Now we can do the analysis but only on terms present in each table_comp list

Pval_col <- c()
OR_col <- c()

for(q in 1:nrow(table_comp)){
  
  test_term <- table_comp[q, 1] 
  test_not <- table_comp[q, 2]
  mean_Sample_term <- table_comp[q, 3]
  sample_not <- table_comp[q, 4]
  
  contingency_table <-
    rbind(c(test_term, test_not),
          c(mean_Sample_term, sample_not))
  FT_out <- fisher.test(contingency_table)
  pval <- FT_out$p.value
  OR <- FT_out$estimate
  
  Pval_col <- c(Pval_col, pval)
  OR_col<- c(OR_col, OR)
  
}    

OR_pval_frame <- data.frame(Pval_col, OR_col)

colnames(OR_pval_frame) <- c("Pval", "OR")
rownames(OR_pval_frame) <- rownames(table_comp)

# Bonferroni correction

OR_pval_frame[, 1] <- p.adjust(OR_pval_frame[, 1], method="bonferroni", n = nrow(OR_pval_frame)) 

# Make plots for this

# i) In each plot order them by the largest significant enrichment out of all of the comparisons 

OR_pval_frame_ordered <- OR_pval_frame[order(OR_pval_frame[, 2], decreasing = T), ]

# Make a bubble plot 

P_1 <- as.numeric(OR_pval_frame_ordered[, 1])

OR_1 <- log10(as.numeric(OR_pval_frame_ordered[, 2]))

neg_log_p_1 <- -log(P_1)

bubble_table_1 <- data.frame(rownames(OR_pval_frame_ordered), OR_1,  P_1, neg_log_p_1)

colnames(bubble_table_1) <- c("Term", "log_OR", "p_val", "neg_log_P")

bubble_table_1$Term <- factor(bubble_table_1$Term, levels = rev(bubble_table_1$Term))

find_size <- c()
find_alpha <- c()

for(i in 1:nrow(bubble_table_1)){
  
  alpha <- 0.3
  size_point <- bubble_table_1[i,4]
  
  if(as.numeric(bubble_table_1[i, 3]) < 0.1){
    alpha <- 1
    size_point <- bubble_table_1[i, 4]
    
  }
  find_alpha <- c(find_alpha, alpha)  
  find_size <- c(find_size, size_point)
}

bubble_table <- cbind(bubble_table_1, find_alpha, find_size)

# just plot first 10

allsignifH1epi<-subset(bubble_table,bubble_table$p_val<0.05)
write.xlsx(allsignifH1epi,"allsignifH1epi.xlsx")
bubble_table_H1 <- bubble_table[1:10, ]
trunc_terms_H1 <- as.character(bubble_table$Term)

#RNA High dose 2 epimut>1 vs all genes
# The test gene list
Sim_test_list <- unique(tot_H2_sup_1$gene)

# The background list
Sim_bg_list <- unique(C_elegans_gene_names[,4])

# Remove any genes from the background list that feature in the test list

intersect <- intersect(Sim_bg_list, Sim_test_list)
Sim_bg_list <- Sim_bg_list[Sim_bg_list %in% intersect == F]

# Send both lists to enrichr
Sim_epimutated_enriched_list  <- enrichr(Sim_test_list, chosendbs)
Sim_epimutated_bg_list  <- enrichr(Sim_bg_list, chosendbs)

Sim_test_bg_list <- list(Sim_epimutated_enriched_list, Sim_epimutated_bg_list)

# Annotate with name of libraries
#write.xlsx(Sim_test_bg_results[[2]],"Sim_test_bg_results_overlap.xlsx")
Sim_test_bg_results <- lapply(Sim_test_bg_list, function(x){
  
  lapply(1:5, function(j){
    
    if(nrow(x[[j]]) > 0){
      x[[j]][1:nrow(x[[j]]), "library"] <- names(Sim_test_bg_list[[1]])[j]
      return(x[[j]])
    } 
    
  })
  
})

# then we combine them and put them in order of significance using the Old.Adjusted.P.value
Sim_test_bg_RESULTS <- lapply(Sim_test_bg_results, function(x){
  
  tempdf <- do.call(rbind, x)
  tempdf[order(tempdf$Old.Adjusted.P.value), ]
  
})

split_value_list <- list()

for(i in 1:length(Sim_test_bg_RESULTS)){
  
  split <- strsplit(Sim_test_bg_RESULTS[[i]][, 2], "/")
  
  split_value <- c()
  
  for(j in 1:length(split)){
    
    save <- split[[j]][1]
    
    split_value <- as.numeric(c(split_value, save))
  }
  
  Sim_test_bg_RESULTS[[i]]$Overlap <- split_value
  
}

# get the test_term, the test_not, the mean(sample_term) and the sample_not

# test_term is overlap
# test_not is length of gene set minus overlap

# then add 1 to all if either test_term OR sample_term = 0
# then only pursue if test_term &/OR sample_term > 4
# let the terms of interest be the terms for the enrichr result for the test list, i.e. with simultaneously epimutated genes

terms_of_interest <- Sim_test_bg_RESULTS[[1]]$Term

# and if the term is not present in the other gene list then assume there is no enrichment for that term in those lists

Test_List <- Sim_test_bg_RESULTS[[1]]
BG_List <- Sim_test_bg_RESULTS[[2]]

table_comp <- c()

for(i in 1:length(terms_of_interest)){
  
  Term_select <- terms_of_interest[[i]]
  sample_term <- 0
  
  if(Term_select %in% BG_List$Term == T){
    
    sample_term <- BG_List[which(BG_List$Term %in% Term_select), 2]}
  sample_not <- length(Sim_bg_list) - sample_term
  
  test_term <- Test_List[which(Test_List$Term %in% Term_select), 2]
  test_not <- length(Sim_test_list) - test_term
  
  # save the components in a row in a new table
  
  component_row <- data.frame(test_term, test_not, sample_term, sample_not)
  
  rownames(component_row) <- Term_select
  
  table_comp <- rbind(table_comp, component_row)
  colnames(table_comp) <- c("test_term", "test_not", "sample_term", "sample_not")
  
}

# Now we will adjust  all values by adding 1 if there are any zeros

zeros <- colSums(table_comp == 0)

if(sum(zeros) > 0){
  
  table_comp <-  table_comp + 1
}

# Now modify the table to remove any rows where test term or mean sample term are < 5
# Now apply the function only if at least one category has 5 or more genes
# It is only reasonable to assess for relative enrichment/depletion if either the test list or sample list has > 5 genes in that ontology category

table_comp <- table_comp[(table_comp$test_term < 5 & table_comp$sample_term < 5)==F, ]

# Now we have found terms where there are at least 5 genes in 1 test or sample category
# we have adjusted the data for each gene set for any 0 values relative to the specific test list
# Now we can do the analysis but only on terms present in each table_comp list

Pval_col <- c()
OR_col <- c()

for(q in 1:nrow(table_comp)){
  
  test_term <- table_comp[q, 1] 
  test_not <- table_comp[q, 2]
  mean_Sample_term <- table_comp[q, 3]
  sample_not <- table_comp[q, 4]
  
  contingency_table <-
    rbind(c(test_term, test_not),
          c(mean_Sample_term, sample_not))
  FT_out <- fisher.test(contingency_table)
  pval <- FT_out$p.value
  OR <- FT_out$estimate
  
  Pval_col <- c(Pval_col, pval)
  OR_col<- c(OR_col, OR)
  
}    

OR_pval_frame <- data.frame(Pval_col, OR_col)

colnames(OR_pval_frame) <- c("Pval", "OR")
rownames(OR_pval_frame) <- rownames(table_comp)

# Bonferroni correction

OR_pval_frame[, 1] <- p.adjust(OR_pval_frame[, 1], method="bonferroni", n = nrow(OR_pval_frame)) 

# Make plots for this

# i) In each plot order them by the largest significant enrichment out of all of the comparisons 

OR_pval_frame_ordered <- OR_pval_frame[order(OR_pval_frame[, 2], decreasing = T), ]

# Make a bubble plot 

P_1 <- as.numeric(OR_pval_frame_ordered[, 1])

OR_1 <- log10(as.numeric(OR_pval_frame_ordered[, 2]))

neg_log_p_1 <- -log(P_1)

bubble_table_1 <- data.frame(rownames(OR_pval_frame_ordered), OR_1,  P_1, neg_log_p_1)

colnames(bubble_table_1) <- c("Term", "log_OR", "p_val", "neg_log_P")

bubble_table_1$Term <- factor(bubble_table_1$Term, levels = rev(bubble_table_1$Term))

find_size <- c()
find_alpha <- c()

for(i in 1:nrow(bubble_table_1)){
  
  alpha <- 0.3
  size_point <- bubble_table_1[i,4]
  
  if(as.numeric(bubble_table_1[i, 3]) < 0.1){
    alpha <- 1
    size_point <- bubble_table_1[i, 4]
    
  }
  find_alpha <- c(find_alpha, alpha)  
  find_size <- c(find_size, size_point)
}

bubble_table <- cbind(bubble_table_1, find_alpha, find_size)

# just plot first 10

allsignifH2epi<-subset(bubble_table,bubble_table$p_val<0.05)
write.xlsx(allsignifH2epi,"allsignifH2epi.xlsx")
bubble_table_H2 <- bubble_table[1:10, ]
trunc_terms_H2 <- as.character(bubble_table$Term)

#Plot all conditions
bubble_table_C1<-bubble_table_C1  %>%
  mutate(Lineage = "C1")
bubble_table_C2<-bubble_table_C2  %>%
  mutate(Lineage = "C2")
bubble_table_L1<-bubble_table_L1 %>%
  mutate(Lineage = "L1")
bubble_table_L2<-bubble_table_L2 %>%
  mutate(Lineage = "L2")
bubble_table_H1<-bubble_table_H1 %>%
  mutate(Lineage = "H1")
bubble_table_H2<-bubble_table_H2 %>%
  mutate(Lineage = "H2")

df_list <- list(bubble_table_C1, bubble_table_C2, bubble_table_L1,bubble_table_L2,bubble_table_H1,bubble_table_H2)
Bubble_table_all_lineage<-Reduce(function(x, y) merge(x, y, all=TRUE), df_list, accumulate=FALSE)

trunc_terms <- c("nuclear lumen (GO:0031981)",
  "RNA binding (GO:0003723)",
  "Ribosome biogenesis in eukaryotes",
  "membrane raft (GO:0045121)",
  "detection of chemical stimulus \n involved in sensory perception (GO:0050907)",
  "sensory perception of \n chemical stimulus (GO:0007606)",
  "integral component of plasma membrane (GO:0005887)",       
  "olfactory receptor activity (GO:0004984)",
  "G-protein coupled receptor \n activity (GO:0004930)",  
  "nucleus (GO:0005634)", 
  "G-protein coupled olfactory \n receptor activity (GO:0038022)" ,  
  "protein kinase activity",
  "defense response to bacterium (GO:0042742)",  
  "protein phosphorylation (GO:0006468)",     
  "protein serine/threonine \n kinase activity (GO:0004674)" ,
  "peptidyl-serine modification (GO:0018209)" ,
  "peptidyl-serine phosphorylation (GO:0018105)",
)

guides_merge <- function(gdefs) {
  gdefs <- lapply(gdefs, function(g) { g$hash <- paste(g$order, g$hash, sep = "z"); g})
  tapply(gdefs, sapply(gdefs, function(g)g$hash), function(gs)Reduce(guide_merge, gs))
}
environment(guides_merge) <- environment(ggplot)
assignInNamespace("guides_merge", guides_merge, pos = "package:ggplot2")

Bubble_table_all_lineage_plot <-
  
  ggplot(Bubble_table_all_lineage, aes(y=Term, x=as.numeric(log_OR)))+
  geom_point(aes(color=Lineage, size=as.numeric(find_size), alpha=find_alpha))+
  scale_x_continuous(breaks = c(1.3,1.4,1.5,1.6), limits = c(1.3, 1.6)) +
  geom_vline(xintercept = 0, colour = "grey")+
  theme_bw()+
  labs(y="Gene Ontology Terms\n", x = "log10(Odds Ratio for enrichment)")+
  scale_size_continuous(range = c(2, 15), breaks = c(25,50,100), 
                        limits = c(20, 110))+  
  scale_colour_manual(values = c("cornflowerblue","blue", "red","2E9FDF","green","darkgreen"), guide=F) +
  scale_alpha(name = paste("Lineage"), range = c(0.3, 1),
              breaks = seq(1, 0.3, length = 6),
              limits = c(0.29, 1.05), 
              labels = c("C1", "C2","L1", "L2", "H1", "H2"))+
  theme(axis.title.x=element_text(face = "bold", size=18))+
  theme(axis.title.y=element_text(face = "bold", size=18))+
  theme(plot.title = element_text(face = "bold", size=16))+
  theme(axis.text.x = element_text(color="#000000", size=16))+
  theme(axis.text.y = element_text(color="#000000", size = 16))
  
Bubble_table_all_lineage_plot

#Get raw data
write.xlsx(Bubble_table_all_lineage,"Data_Sup_Fig_4_D.xlsx")
#--------------------------
#####Fig.3.E - RNA epimutations according to chromatin domain####
Ahringer_Chromatin_domain<-Ahringer[ , c("Gene", "Chromatin_domain","Tissue_specificity")]   
colnames(RNA_epimut)<-c("Position", "Gene","number_transitions","onset_gen","length","complete","Lineage","is_up","is_down","Condition")
RNA_allepiduration_Ahringer <- merge(RNA_epimut,Ahringer_Chromatin_domain,by="Gene")
RNA_all_true_epiduration_Ahringer <- subset(RNA_allepiduration_Ahringer,RNA_allepiduration_Ahringer$length >1)
colnames(C_elegans_gene_names)<-c("V1","V2","V3","Gene","V5","V6")
all_genes_Ahringer <- merge(C_elegans_gene_names,Ahringer_Chromatin_domain,by="Gene")
all_genes_Ahringer_genes<-unique(all_genes_Ahringer$Gene)

# The lists of length categorised expression changes must be restricted to genes which have Ahringer annotations 

Ahr_RNA_Control <- subset(RNA_all_true_epiduration_Ahringer,RNA_all_true_epiduration_Ahringer$Condition=="Control")
Ahr_RNA_Control_genes<-unique(Ahr_RNA_Control$Gene)
Ahr_RNA_Low_dose <- subset(RNA_all_true_epiduration_Ahringer,RNA_all_true_epiduration_Ahringer$Condition=="Low dose")
Ahr_RNA_Low_dose_genes<-unique(Ahr_RNA_Low_dose$Gene)
Ahr_RNA_High_dose <- subset(RNA_all_true_epiduration_Ahringer,RNA_all_true_epiduration_Ahringer$Condition=="High dose")
Ahr_RNA_High_dose_genes<-unique(Ahr_RNA_High_dose$Gene)

# RNA_background are genes which map to Ahringer

Active_genes <- unique(all_genes_Ahringer)[all_genes_Ahringer$Chromatin_domain == "A", 1]
Active_genes<-unique(Active_genes)
Regulated_genes <- unique(all_genes_Ahringer[all_genes_Ahringer$Chromatin_domain == "R", 1])
Regulated_genes<-unique(Regulated_genes)
Chr_X_genes <- unique(all_genes_Ahringer[all_genes_Ahringer$Chromatin_domain == ".", 1])
Chr_X_genes<-unique(Chr_X_genes)

RNA_test_list <- list(Ahr_RNA_Control_genes, Ahr_RNA_Low_dose_genes, Ahr_RNA_High_dose_genes)
RNA_Features <- list(Regulated_genes, Active_genes, piRNA_cluster_genes, Chr_X_genes)

# Define the function

epiIntersect <-
  function(common,
           Epimutation_total,
           Feature_total,
           overall_total) {
    contingency_table <-
      rbind(c(common, Epimutation_total - common),
            c(
              Feature_total - common,
              overall_total - (Epimutation_total + Feature_total - common)
            ))
    FT_out <- fisher.test(contingency_table)
    pval <- FT_out$p.value
    OddsRatio <- FT_out$estimate
    return(c(pval, OddsRatio))
  }

Non_Nested_RNA_OR <- c()
Non_Nested_RNA_pval <- c()

for(d in 1:length(RNA_test_list)){
  
  test <- RNA_test_list[[d]]
  
  temp_pval <- matrix(0, ncol = 1, nrow = 4)
  temp_oR <- matrix(0, ncol = 1, nrow = 4)
  
  for (i in 1:length(RNA_Features)) {
    common_in <- length(intersect(RNA_Features[[i]], test))
    Epimutation_total_in <- length(test)
    Feature_total_in <- length(RNA_Features[[i]])
    
    output <- epiIntersect(
      common = common_in,
      Epimutation_total = Epimutation_total_in,
      Feature_total = Feature_total_in,
      overall_total = 12623  # Total number of RNA genes mapping to Ahringer
    )
    temp_pval[i, ] <- output[1]
    temp_oR[i, ] <- output[2]
  }
  
  temp_pval <- as.matrix(p.adjust(temp_pval, "bonferroni", 4))
  
  colnames(temp_pval) <-
    c("pval")
  rownames(temp_pval) <-
    c("Regulated", "Active", "piRNA", "ChrX")
  
  colnames(temp_oR) <-
    c("oR")
  rownames(temp_oR) <-
    c("Regulated", "Active", "piRNA", "ChrX")
  
  Non_Nested_RNA_OR <- cbind(Non_Nested_RNA_OR, temp_oR)
  
  Non_Nested_RNA_pval <- cbind(Non_Nested_RNA_pval, temp_pval)
}

colnames(Non_Nested_RNA_OR) <- c("Control", "Low dose", "High dose")
colnames(Non_Nested_RNA_pval) <- c("Control", "Low dose", "High dose")

# RNA epimutations distributions in distinct chromatin domains

ordered_sequence <- c("piRNA", "ChrX", "Regulated", "Active")
colnames(Non_Nested_RNA_OR) <- c("Test_1", "Test_2", "Test_3")
colnames(Non_Nested_RNA_pval) <- c("Test_1", "Test_2", "Test_3")

Non_Nested_RNA_OR <- as.data.frame(Non_Nested_RNA_OR[ordered_sequence, ])
Non_Nested_RNA_pval <- as.data.frame(Non_Nested_RNA_pval[ordered_sequence, ])

P_1 <- as.numeric(Non_Nested_RNA_pval$Test_1)
log_Odds_1 <- log2(as.numeric(Non_Nested_RNA_OR$Test_1))
neg_log_p_1 <- -log2(P_1)
bubble_table_1 <- cbind(rownames(Non_Nested_RNA_OR), Non_Nested_RNA_pval$Test_1, neg_log_p_1, log_Odds_1)
bubble_table_1 <- data.frame(rep("Control"), bubble_table_1)
bubble_table_1[, 2] <- factor(bubble_table_1[, 2], levels = bubble_table_1[, 2])
colnames(bubble_table_1) <- c("Analysis",  "Description", "p_val", "neg_log_P", "log_odds")

P_2 <- as.numeric(Non_Nested_RNA_pval$Test_2)
log_Odds_2 <- log2(as.numeric(Non_Nested_RNA_OR$Test_2))
neg_log_p_2 <- -log2(P_2)
bubble_table_2 <- cbind(rownames(Non_Nested_RNA_OR), Non_Nested_RNA_pval$Test_2, neg_log_p_2, log_Odds_2)
bubble_table_2 <- data.frame(rep("Low dose"), bubble_table_2)
bubble_table_2[, 2] <- factor(bubble_table_2[, 2], levels = bubble_table_2[, 2])
colnames(bubble_table_2) <- c("Analysis",  "Description", "p_val", "neg_log_P", "log_odds")

P_3 <- as.numeric(Non_Nested_RNA_pval$Test_3)
log_Odds_3 <- log2(as.numeric(Non_Nested_RNA_OR$Test_3))
neg_log_p_3 <- -log2(P_3)
bubble_table_3 <- cbind(rownames(Non_Nested_RNA_OR), Non_Nested_RNA_pval$Test_3, neg_log_p_3, log_Odds_3)
bubble_table_3 <- data.frame(rep("High dose"), bubble_table_3)
bubble_table_3[, 2] <- factor(bubble_table_3[, 2], levels = bubble_table_3[, 2])
colnames(bubble_table_3) <- c("Analysis",  "Description", "p_val", "neg_log_P", "log_odds")

c_tests <- rbind(bubble_table_1, bubble_table_2, bubble_table_3)

find_size <- c()
find_alpha <- c()

for(i in 1:nrow(c_tests)){
  
  alpha <- 0.3
  size_point <- abs(log2(as.numeric(c_tests[i, 4]))+log2(as.numeric(c_tests[i, 4])))
  
  if(as.numeric(c_tests[i, 3]) < 0.05){
    alpha <- 1
    size_point <- abs(log2(as.numeric(c_tests[i, 4]))+log2(as.numeric(c_tests[i, 4])))
  }
  
  if(as.numeric(c_tests[i, 3]) == 1){
    alpha <- 0.3
    size_point <- 1
  }
  
  find_alpha <- c(find_alpha, alpha)  
  find_size <- c(find_size, size_point)
}

RNA_c_combined_tests <- cbind(c_tests, find_alpha, find_size)

guides_merge <- function(gdefs) {
  gdefs <- lapply(gdefs, function(g) { g$hash <- paste(g$order, g$hash, sep = "z"); g})
  tapply(gdefs, sapply(gdefs, function(g)g$hash), function(gs)Reduce(guide_merge, gs))
}
environment(guides_merge) <- environment(ggplot)
assignInNamespace("guides_merge", guides_merge, pos = "package:ggplot2")

RNA_c_combined_tests$Analysis <- factor(RNA_c_combined_tests$Analysis, levels = unique(RNA_c_combined_tests$Analysis))

# Bubble plot to show distribution of RNAs in chromatin domains according to length

RNA_Lengths_Domain <-
  
  ggplot(RNA_c_combined_tests, aes(y=Description, x=as.numeric(log_odds)))+
  geom_point(aes(color=Analysis, size=find_size, alpha=find_alpha))+
  expand_limits(x=c(-1, 2))+
  geom_vline(xintercept = 0, colour = "grey")+
  theme_bw()+
  scale_color_manual(values = c("cornflowerblue", "darkgreen", "red"))+
  
  scale_size_continuous(range = c(0, 70), breaks = c(1, 5, 10,15), 
                        limits = c(0, 70))+
  
  scale_alpha(name = paste("Transparency of bubble represents \nsignificance of enrichment"), range = c(0.3, 1),
              breaks = seq(1, 0.3, length = 2),
              limits = c(0.29, 1.05), 
              labels = c("p value < 0.05", "p value > 0.05"))

RNA_Lengths_Domain <- 
  
  RNA_Lengths_Domain + labs(y="Chromatin Domains\n", x=paste('log2(Odds Ratio) for enrichment'))+
  theme(axis.title=element_text(face = "bold", size=14))+
  theme(plot.title = element_text(face = "bold", size=16))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(color="#000000", size=10))+
  theme(axis.text.y = element_text(color="#000000", size = 14))+
  guides(size = guide_legend(order = 3), colour = guide_legend(order = 1), alpha = guide_legend(order = 2))+
  theme(legend.text=element_text(color="#000000", size=12))

RNA_Lengths_Domain 

#Get raw data
write.xlsx(RNA_c_combined_tests,"Data_Fig_3_E.xlsx")

#--------------------------
#####Sup.Fig.4.E - RNA epimutations according to chromatin domain####
# The lists of length categorised expression changes must be restricted to genes which have Ahringer annotations 

Ahr_RNA_Control_1 <- subset(RNA_all_true_epiduration_Ahringer,RNA_all_true_epiduration_Ahringer$Condition=="Control" & RNA_all_true_epiduration_Ahringer$Lineage=="C1")
Ahr_RNA_Control_1_genes<-unique(Ahr_RNA_Control_1$Gene)
Ahr_RNA_Control_2 <- subset(RNA_all_true_epiduration_Ahringer,RNA_all_true_epiduration_Ahringer$Condition=="Control" & RNA_all_true_epiduration_Ahringer$Lineage=="C2")
Ahr_RNA_Control_2_genes<-unique(Ahr_RNA_Control_2$Gene)

Ahr_RNA_Low_dose_1 <- subset(RNA_all_true_epiduration_Ahringer,RNA_all_true_epiduration_Ahringer$Condition=="Low dose" & RNA_all_true_epiduration_Ahringer$Lineage=="L1")
Ahr_RNA_Low_dose_1_genes<-unique(Ahr_RNA_Low_dose_1$Gene)
Ahr_RNA_Low_dose_2 <- subset(RNA_all_true_epiduration_Ahringer,RNA_all_true_epiduration_Ahringer$Condition=="Low dose" & RNA_all_true_epiduration_Ahringer$Lineage=="L2")
Ahr_RNA_Low_dose_2_genes<-unique(Ahr_RNA_Low_dose_2$Gene)

Ahr_RNA_High_dose_1 <- subset(RNA_all_true_epiduration_Ahringer,RNA_all_true_epiduration_Ahringer$Condition=="High dose" & RNA_all_true_epiduration_Ahringer$Lineage=="H1")
Ahr_RNA_High_dose_1_genes<-unique(Ahr_RNA_High_dose_1$Gene)
Ahr_RNA_High_dose_2 <- subset(RNA_all_true_epiduration_Ahringer,RNA_all_true_epiduration_Ahringer$Condition=="High dose" & RNA_all_true_epiduration_Ahringer$Lineage=="H2")
Ahr_RNA_High_dose_2_genes<-unique(Ahr_RNA_High_dose_2$Gene)

# RNA_background are genes which map to Ahringer

Active_genes <- unique(all_genes_Ahringer)[all_genes_Ahringer$Chromatin_domain == "A", 1]
Active_genes<-unique(Active_genes)
Regulated_genes <- unique(all_genes_Ahringer[all_genes_Ahringer$Chromatin_domain == "R", 1])
Regulated_genes<-unique(Regulated_genes)
Chr_X_genes <- unique(all_genes_Ahringer[all_genes_Ahringer$Chromatin_domain == ".", 1])
Chr_X_genes<-unique(Chr_X_genes)

RNA_test_list <- list(Ahr_RNA_Control_1_genes, Ahr_RNA_Control_2_genes, Ahr_RNA_Low_dose_1_genes, Ahr_RNA_Low_dose_2_genes, Ahr_RNA_High_dose_1_genes, Ahr_RNA_High_dose_2_genes)
RNA_Features <- list(Regulated_genes, Active_genes, piRNA_cluster_genes, Chr_X_genes)

# Define the function

epiIntersect <-
  function(common,
           Epimutation_total,
           Feature_total,
           overall_total) {
    contingency_table <-
      rbind(c(common, Epimutation_total - common),
            c(
              Feature_total - common,
              overall_total - (Epimutation_total + Feature_total - common)
            ))
    FT_out <- fisher.test(contingency_table)
    pval <- FT_out$p.value
    OddsRatio <- FT_out$estimate
    return(c(pval, OddsRatio))
  }

Non_Nested_RNA_OR <- c()
Non_Nested_RNA_pval <- c()

for(d in 1:length(RNA_test_list)){
  
  test <- RNA_test_list[[d]]
  
  temp_pval <- matrix(0, ncol = 1, nrow = 4)
  temp_oR <- matrix(0, ncol = 1, nrow = 4)
  
  for (i in 1:length(RNA_Features)) {
    common_in <- length(intersect(RNA_Features[[i]], test))
    Epimutation_total_in <- length(test)
    Feature_total_in <- length(RNA_Features[[i]])
    
    output <- epiIntersect(
      common = common_in,
      Epimutation_total = Epimutation_total_in,
      Feature_total = Feature_total_in,
      overall_total = 12623  # Total number of RNA genes mapping to Ahringer
    )
    temp_pval[i, ] <- output[1]
    temp_oR[i, ] <- output[2]
  }
  
  temp_pval <- as.matrix(p.adjust(temp_pval, "bonferroni", 4))
  
  colnames(temp_pval) <-
    c("pval")
  rownames(temp_pval) <-
    c("Regulated", "Active", "piRNA", "ChrX")
  
  colnames(temp_oR) <-
    c("oR")
  rownames(temp_oR) <-
    c("Regulated", "Active", "piRNA", "ChrX")
  
  Non_Nested_RNA_OR <- cbind(Non_Nested_RNA_OR, temp_oR)
  
  Non_Nested_RNA_pval <- cbind(Non_Nested_RNA_pval, temp_pval)
}

colnames(Non_Nested_RNA_OR) <- c("C1", "C2", "L1","L2", "H1", "H2")
colnames(Non_Nested_RNA_pval) <- c("C1", "C2", "L1","L2", "H1", "H2")

# RNA epimutations distributions in distinct chromatin domains

ordered_sequence <- c("piRNA", "ChrX", "Regulated", "Active")
colnames(Non_Nested_RNA_OR) <- c("Test_1", "Test_2", "Test_3", "Test_4", "Test_5", "Test_6")
colnames(Non_Nested_RNA_pval) <- c("Test_1", "Test_2", "Test_3", "Test_4", "Test_5", "Test_6")

Non_Nested_RNA_OR <- as.data.frame(Non_Nested_RNA_OR[ordered_sequence, ])
Non_Nested_RNA_pval <- as.data.frame(Non_Nested_RNA_pval[ordered_sequence, ])

P_1 <- as.numeric(Non_Nested_RNA_pval$Test_1)
log_Odds_1 <- log2(as.numeric(Non_Nested_RNA_OR$Test_1))
neg_log_p_1 <- -log2(P_1)
bubble_table_1 <- cbind(rownames(Non_Nested_RNA_OR), Non_Nested_RNA_pval$Test_1, neg_log_p_1, log_Odds_1)
bubble_table_1 <- data.frame(rep("C1"), bubble_table_1)
bubble_table_1[, 2] <- factor(bubble_table_1[, 2], levels = bubble_table_1[, 2])
colnames(bubble_table_1) <- c("Analysis",  "Description", "p_val", "neg_log_P", "log_odds")

P_2 <- as.numeric(Non_Nested_RNA_pval$Test_2)
log_Odds_2 <- log2(as.numeric(Non_Nested_RNA_OR$Test_2))
neg_log_p_2 <- -log2(P_2)
bubble_table_2 <- cbind(rownames(Non_Nested_RNA_OR), Non_Nested_RNA_pval$Test_2, neg_log_p_2, log_Odds_2)
bubble_table_2 <- data.frame(rep("C2"), bubble_table_2)
bubble_table_2[, 2] <- factor(bubble_table_2[, 2], levels = bubble_table_2[, 2])
colnames(bubble_table_2) <- c("Analysis",  "Description", "p_val", "neg_log_P", "log_odds")

P_3 <- as.numeric(Non_Nested_RNA_pval$Test_3)
log_Odds_3 <- log2(as.numeric(Non_Nested_RNA_OR$Test_3))
neg_log_p_3 <- -log2(P_3)
bubble_table_3 <- cbind(rownames(Non_Nested_RNA_OR), Non_Nested_RNA_pval$Test_3, neg_log_p_3, log_Odds_3)
bubble_table_3 <- data.frame(rep("L1"), bubble_table_3)
bubble_table_3[, 2] <- factor(bubble_table_3[, 2], levels = bubble_table_3[, 2])
colnames(bubble_table_3) <- c("Analysis",  "Description", "p_val", "neg_log_P", "log_odds")

P_4 <- as.numeric(Non_Nested_RNA_pval$Test_4)
log_Odds_4 <- log2(as.numeric(Non_Nested_RNA_OR$Test_4))
neg_log_p_4 <- -log2(P_4)
bubble_table_4 <- cbind(rownames(Non_Nested_RNA_OR), Non_Nested_RNA_pval$Test_4, neg_log_p_4, log_Odds_4)
bubble_table_4 <- data.frame(rep("L2"), bubble_table_4)
bubble_table_4[, 2] <- factor(bubble_table_4[, 2], levels = bubble_table_4[, 2])
colnames(bubble_table_4) <- c("Analysis",  "Description", "p_val", "neg_log_P", "log_odds")

P_5 <- as.numeric(Non_Nested_RNA_pval$Test_5)
log_Odds_5 <- log2(as.numeric(Non_Nested_RNA_OR$Test_5))
neg_log_p_5 <- -log2(P_5)
bubble_table_5 <- cbind(rownames(Non_Nested_RNA_OR), Non_Nested_RNA_pval$Test_5, neg_log_p_5, log_Odds_5)
bubble_table_5 <- data.frame(rep("H1"), bubble_table_5)
bubble_table_5[, 2] <- factor(bubble_table_5[, 2], levels = bubble_table_5[, 2])
colnames(bubble_table_5) <- c("Analysis",  "Description", "p_val", "neg_log_P", "log_odds")

P_6 <- as.numeric(Non_Nested_RNA_pval$Test_6)
log_Odds_6 <- log2(as.numeric(Non_Nested_RNA_OR$Test_6))
neg_log_p_6 <- -log2(P_6)
bubble_table_6 <- cbind(rownames(Non_Nested_RNA_OR), Non_Nested_RNA_pval$Test_6, neg_log_p_6, log_Odds_6)
bubble_table_6 <- data.frame(rep("H2"), bubble_table_6)
bubble_table_6[, 2] <- factor(bubble_table_6[, 2], levels = bubble_table_6[, 2])
colnames(bubble_table_6) <- c("Analysis",  "Description", "p_val", "neg_log_P", "log_odds")

c_tests <- rbind(bubble_table_1, bubble_table_2, bubble_table_3, bubble_table_4, bubble_table_5, bubble_table_6)

find_size <- c()
find_alpha <- c()

for(i in 1:nrow(c_tests)){
  
  alpha <- 0.3
  size_point <- abs(log2(as.numeric(c_tests[i, 4]))+log2(as.numeric(c_tests[i, 4])))
  
  if(as.numeric(c_tests[i, 3]) < 0.05){
    alpha <- 1
    size_point <- abs(log2(as.numeric(c_tests[i, 4]))+log2(as.numeric(c_tests[i, 4])))
  }
  
  if(as.numeric(c_tests[i, 3]) == 1){
    alpha <- 0.3
    size_point <- 1
  }
  
  find_alpha <- c(find_alpha, alpha)  
  find_size <- c(find_size, size_point)
}

RNA_c_combined_tests <- cbind(c_tests, find_alpha, find_size)

guides_merge <- function(gdefs) {
  gdefs <- lapply(gdefs, function(g) { g$hash <- paste(g$order, g$hash, sep = "z"); g})
  tapply(gdefs, sapply(gdefs, function(g)g$hash), function(gs)Reduce(guide_merge, gs))
}
environment(guides_merge) <- environment(ggplot)
assignInNamespace("guides_merge", guides_merge, pos = "package:ggplot2")

RNA_c_combined_tests$Analysis <- factor(RNA_c_combined_tests$Analysis, levels = unique(RNA_c_combined_tests$Analysis))

# Bubble plot to show distribution of RNAs in chromatin domains according to length

RNA_Lengths_Domain <-
  
  ggplot(RNA_c_combined_tests, aes(y=Description, x=as.numeric(log_odds)))+
  geom_point(aes(color=Analysis, size=find_size, alpha=find_alpha))+
  expand_limits(x=c(-1, 2))+
  geom_vline(xintercept = 0, colour = "grey")+
  theme_bw()+
  scale_color_manual(values = c("cornflowerblue","blue","green","darkgreen", "red","2E9FDF"))+
  
  scale_size_continuous(range = c(0, 50), breaks = c(1, 7 ,15), 
                        limits = c(0, 50))+
  
  scale_alpha(name = paste("Transparency of bubble represents \nsignificance of enrichment"), range = c(0.3, 1),
              breaks = seq(1, 0.3, length = 2),
              limits = c(0.29, 1.05), 
              labels = c("p value < 0.05", "p value > 0.05"))

RNA_Lengths_Domain <- 
  
  RNA_Lengths_Domain + labs(y="Chromatin Domains\n", x=paste('log2(Odds Ratio) for enrichment'))+
  theme(axis.title=element_text(face = "bold", size=14))+
  theme(plot.title = element_text(face = "bold", size=16))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(color="#000000", size=10))+
  theme(axis.text.y = element_text(color="#000000", size = 14))+
  guides(size = guide_legend(order = 3), colour = guide_legend(order = 1), alpha = guide_legend(order = 2))+
  theme(legend.text=element_text(color="#000000", size=12))

RNA_Lengths_Domain 

#Get raw data
write.xlsx(RNA_c_combined_tests,"Data_Sup_Fig_4_E.xlsx")
#--------------------------
#####Fig.4.A & Sup.Fig.5.A - 22G-RNAs epimutations####
#Importation of sncRNAs data
setwd("~/Documents/Cisplatin project/Data_analysis/Count data")
alltiny<-read.xlsx("tinyrna_2022-10-19_10-36-23_feature_counts.xlsx")
colnames(alltiny)<-c("Feature.ID","Tag","Feature.Name","Feature.Class","C1_0","C1_18","C1_20","C1_2","C1_4","C1_6","C1_8","C1_10","C1_12","C1_14","C1_16",
                     "C2_0","C2_20","C2_2","C2_4","C2_8","C2_10","C2_12","C2_14","C2_16","C2_18",
                     "L1_0","L1_20","L1_2","L1_4","L1_6","L1_8","L1_12","L1_14","L1_16","L1_18",
                     "L2_0","L2_20","L2_2","L2_4","L2_8","L2_10","L2_12","L2_14","L2_16","L2_18",
                     "H1_0","H1_20","H1_2","H1_4","H1_8","H1_10","H1_12","H1_14","H1_16","H1_18",
                     "H2_0","H2_18","H2_20","H2_2","H2_4","H2_6","H2_8","H2_10","H2_12","H2_14","H2_16")
alltiny_genes<-read.xlsx("alltiny_genes.xlsx")
alltiny<-cbind(alltiny,alltiny_genes)
setwd("~/Documents/Cisplatin project/Data_analysis/Analysis")
count_matrix_tinyRNA<-alltiny[,5:66]
sample_info_tinyRNA<-data.frame(condition=colnames(count_matrix_tinyRNA),row.names=colnames(count_matrix_tinyRNA))

DESeq.ds_tinyRNA<-DESeqDataSetFromMatrix(countData = round(count_matrix_tinyRNA),
                                         colData=sample_info_tinyRNA,
                                         design=~condition)

# obtain regularized log - transformed values
vst_tinyRNA<-vst(DESeq.ds_tinyRNA,blind = TRUE)
vst.norm.count.timyRNA<-assay(vst_tinyRNA)

#calculate the size factor and add it to the data set
DESeq.ds_tinyRNA<-estimateSizeFactors(DESeq.ds_tinyRNA)
DESeq.ds_tinyRNA_sf<-sizeFactors(DESeq.ds_tinyRNA)
colData(DESeq.ds_tinyRNA)

# retrieve the _ normalized _ read counts
countData_normalized_tinyRNA<-counts(DESeq.ds_tinyRNA,normalized=TRUE)
sizeFactOut<-sizeFactors(DESeq.ds_tinyRNA)

# transform size - factor normalized read counts to log2 scale using a pseudocount of 1
log.norm.counts_tinyRNA<-log2(countData_normalized_tinyRNA+1)
boxplot(log.norm.counts_tinyRNA,notch=TRUE,
        main="log2-transformed read counts",
        ylab="log2(read counts)")

#Detection of small-RNAs epimutations
#Create data frame per lineage
log.norm.counts_tinyRNA<-as.data.frame(log.norm.counts_tinyRNA)
log.norm.counts_tinyRNA_Lineage_C1 <- cbind(log.norm.counts_tinyRNA$C1_0,log.norm.counts_tinyRNA$C1_2,log.norm.counts_tinyRNA$C1_4,log.norm.counts_tinyRNA$C1_6,log.norm.counts_tinyRNA$C1_8,log.norm.counts_tinyRNA$C1_10,log.norm.counts_tinyRNA$C1_12,log.norm.counts_tinyRNA$C1_14,log.norm.counts_tinyRNA$C1_16,log.norm.counts_tinyRNA$C1_18,log.norm.counts_tinyRNA$C1_20)
colnames(log.norm.counts_tinyRNA_Lineage_C1) <- c("C1_0", "C1_2", "C1_4", "C1_6", "C1_8", "C1_10","C1_12","C1_14", "C1_16", "C1_18","C1_20")
log.norm.counts_tinyRNA_Lineage_C2 <- cbind(log.norm.counts_tinyRNA$C2_0,log.norm.counts_tinyRNA$C2_2,log.norm.counts_tinyRNA$C2_4,log.norm.counts_tinyRNA$C2_8,log.norm.counts_tinyRNA$C2_10,log.norm.counts_tinyRNA$C2_12,log.norm.counts_tinyRNA$C2_14,log.norm.counts_tinyRNA$C2_16,log.norm.counts_tinyRNA$C2_18,log.norm.counts_tinyRNA$C2_20)
colnames(log.norm.counts_tinyRNA_Lineage_C2) <- c("C2_0", "C2_2", "C2_4", "C2_8", "C2_10","C2_12","C2_14", "C2_16", "C2_18","C2_20")
log.norm.counts_tinyRNA_Lineage_L1 <- cbind(log.norm.counts_tinyRNA$L1_0,log.norm.counts_tinyRNA$L1_2,log.norm.counts_tinyRNA$L1_4,log.norm.counts_tinyRNA$L1_6,log.norm.counts_tinyRNA$L1_8,log.norm.counts_tinyRNA$L1_12,log.norm.counts_tinyRNA$L1_14,log.norm.counts_tinyRNA$L1_16,log.norm.counts_tinyRNA$L1_18,log.norm.counts_tinyRNA$L1_20)
colnames(log.norm.counts_tinyRNA_Lineage_L1) <- c("L1_0", "L1_2", "L1_4", "L1_6", "L1_8", "L1_12", "L1_14", "L1_16", "L1_18","L1_20")
log.norm.counts_tinyRNA_Lineage_L2 <- cbind(log.norm.counts_tinyRNA$L2_0,log.norm.counts_tinyRNA$L2_2,log.norm.counts_tinyRNA$L2_4,log.norm.counts_tinyRNA$L2_8,log.norm.counts_tinyRNA$L2_10,log.norm.counts_tinyRNA$L2_12,log.norm.counts_tinyRNA$L2_14,log.norm.counts_tinyRNA$L2_16,log.norm.counts_tinyRNA$L2_18,log.norm.counts_tinyRNA$L2_20)
colnames(log.norm.counts_tinyRNA_Lineage_L2) <- c("L2_0", "L2_2", "L2_4", "L2_8", "L2_10", "L2_12", "L2_14", "L2_16", "L2_18","L2_20")
log.norm.counts_tinyRNA_Lineage_H1 <- cbind(log.norm.counts_tinyRNA$H1_0,log.norm.counts_tinyRNA$H1_2,log.norm.counts_tinyRNA$H1_4,log.norm.counts_tinyRNA$H1_8,log.norm.counts_tinyRNA$H1_10,log.norm.counts_tinyRNA$H1_12,log.norm.counts_tinyRNA$H1_14,log.norm.counts_tinyRNA$H1_16,log.norm.counts_tinyRNA$H1_18,log.norm.counts_tinyRNA$H1_20)
colnames(log.norm.counts_tinyRNA_Lineage_H1) <- c("H1_0", "H1_2", "H1_4", "H1_8", "H1_10", "H1_12", "H1_14", "H1_16", "H1_18", "H1_20")
log.norm.counts_tinyRNA_Lineage_H2 <- cbind(log.norm.counts_tinyRNA$H2_0,log.norm.counts_tinyRNA$H2_2,log.norm.counts_tinyRNA$H2_4,log.norm.counts_tinyRNA$H2_6,log.norm.counts_tinyRNA$H2_8,log.norm.counts_tinyRNA$H2_10,log.norm.counts_tinyRNA$H2_12,log.norm.counts_tinyRNA$H2_14,log.norm.counts_tinyRNA$H2_16,log.norm.counts_tinyRNA$H2_18,log.norm.counts_tinyRNA$H2_20)
colnames(log.norm.counts_tinyRNA_Lineage_H2) <- c("H2_0", "H2_2", "H2_4", "H2_6", "H2_8", "H2_10", "H2_12", "H2_14", "H2_16", "H2_18","H2_20")

#C1
#Identify epimutations in tinyRNA using linear model in each lineage 
tinyRNA_Lineage_C1_2_to_20<-log.norm.counts_tinyRNA_Lineage_C1[,2:11]
tinyRNA_z_scores_table_C1<-matrix(0,ncol=10, nrow = nrow(log.norm.counts_tinyRNA_Lineage_C1))

for (i in 1:ncol(tinyRNA_Lineage_C1_2_to_20)) {
  tinyRNA_EpimutC1<-
    lm(log.norm.counts_tinyRNA_Lineage_C1[,1]~tinyRNA_Lineage_C1_2_to_20[,i])
  tinyRNA_ResidC1<-tinyRNA_EpimutC1$residuals
  tinyRNA_z_scores_table_C1[,i]<-(tinyRNA_ResidC1-mean(tinyRNA_ResidC1))/sd(tinyRNA_ResidC1)
}

#keep only relevant epimutations and binarised data
tiniyRNA_binarised_z_scores_table_C1<-matrix(0,ncol=ncol(tinyRNA_z_scores_table_C1),nrow=nrow(tinyRNA_z_scores_table_C1))

for (i in 1:ncol(tinyRNA_z_scores_table_C1)) {
  for(j in 1:nrow(tinyRNA_z_scores_table_C1)){
    
    if(tinyRNA_z_scores_table_C1[j,i]> (-2.25) & tinyRNA_z_scores_table_C1[j,i]< 2.25) {
      tiniyRNA_binarised_z_scores_table_C1[j,i]<-0
    }
    if(tinyRNA_z_scores_table_C1[j,i] > 2.25){
      tiniyRNA_binarised_z_scores_table_C1[j,i]<-1
    }
    if(tinyRNA_z_scores_table_C1[j,i] < (-2.25)){
      tiniyRNA_binarised_z_scores_table_C1[j,i]<- (-1)
    }
  }  
}
tinyRNA_C1epimutations<-tiniyRNA_binarised_z_scores_table_C1 %>% as.data.frame() 
colnames(tinyRNA_C1epimutations) <- c('F2','4','6','8','10','12','14','16','18','20')
tinyRNA_C1epimutations <- tinyRNA_C1epimutations %>%
  # Creating an empty column:
  add_column(F0 = 0, .before='F2')
tinyRNA_C1epimutations <- tinyRNA_C1epimutations %>%
  # Creating a condition column:
  add_column(condition = "C1", .before="F0")
tinyRNA_C1epimutations <- cbind(alltiny[,1:4],alltiny[,67],tinyRNA_C1epimutations)
colnames(tinyRNA_C1epimutations) <- c('Feature.ID','Tag','Feature.Name','Feature.Class','genes','condition','0','2','4','6','8','10','12','14','16','18','20')
save(tinyRNA_C1epimutations,file="tinyRNA_C1epimutations.Rdata")

#see number of epimutations per generations
tinyRNA_C1epimutations_22G<-subset(tinyRNA_C1epimutations,tinyRNA_C1epimutations$Tag=="22G")
tinyRNA_C1epimutations_26G<-subset(tinyRNA_C1epimutations,tinyRNA_C1epimutations$Tag=="26G")
tinyRNA_C1epimutations_piRNA<-subset(tinyRNA_C1epimutations,tinyRNA_C1epimutations$Feature.Class=="piRNA")
tinyRNA_C1epimutations_miRNA<-subset(tinyRNA_C1epimutations,tinyRNA_C1epimutations$Feature.Class=="miRNA")
save(tinyRNA_C1epimutations_22G,file="tinyRNA_C1epimutations_22G.Rdata")

#C2
#Identify epimutations in tinyRNA using linear model in each lineage 
tinyRNA_Lineage_C2_2_to_20<-log.norm.counts_tinyRNA_Lineage_C2[,2:10]
tinyRNA_z_scores_table_C2<-matrix(0,ncol=9, nrow = nrow(log.norm.counts_tinyRNA_Lineage_C2))

for (i in 1:ncol(tinyRNA_Lineage_C2_2_to_20)) {
  tinyRNA_EpimutC2<-
    lm(log.norm.counts_tinyRNA_Lineage_C2[,1]~tinyRNA_Lineage_C2_2_to_20[,i])
  tinyRNA_ResidC2<-tinyRNA_EpimutC2$residuals
  tinyRNA_z_scores_table_C2[,i]<-(tinyRNA_ResidC2-mean(tinyRNA_ResidC2))/sd(tinyRNA_ResidC2)
}

#keep only relevant epimutations and binarised data
tiniyRNA_binarised_z_scores_table_C2<-matrix(0,ncol=ncol(tinyRNA_z_scores_table_C2),nrow=nrow(tinyRNA_z_scores_table_C2))

for (i in 1:ncol(tinyRNA_z_scores_table_C2)) {
  for(j in 1:nrow(tinyRNA_z_scores_table_C2)){
    
    if(tinyRNA_z_scores_table_C2[j,i]> (-2.25) & tinyRNA_z_scores_table_C2[j,i]< 2.25) {
      tiniyRNA_binarised_z_scores_table_C2[j,i]<-0
    }
    if(tinyRNA_z_scores_table_C2[j,i] > 2.25){
      tiniyRNA_binarised_z_scores_table_C2[j,i]<-1
    }
    if(tinyRNA_z_scores_table_C2[j,i] < (-2.25)){
      tiniyRNA_binarised_z_scores_table_C2[j,i]<- (-1)
    }
  }  
}
tinyRNA_C2epimutations<-tiniyRNA_binarised_z_scores_table_C2 %>% as.data.frame() 
colnames(tinyRNA_C2epimutations) <- c('F2','4','8','10','12','14','16','18','20')
tinyRNA_C2epimutations <- tinyRNA_C2epimutations %>%
  # Creating an empty column:
  add_column(F0 = 0, .before='F2')
tinyRNA_C2epimutations <- tinyRNA_C2epimutations %>%
  # Creating a condition column:
  add_column(condition = "C2", .before="F0")
tinyRNA_C2epimutations <- cbind(alltiny[,1:4],alltiny[,67],tinyRNA_C2epimutations)
colnames(tinyRNA_C2epimutations) <- c('Feature.ID','Tag','Feature.Name','Feature.Class','genes','condition','0','2','4','8','10','12','14','16','18','20')
save(tinyRNA_C2epimutations,file="tinyRNA_C2epimutations.Rdata")

#see number of epimutations per generations
tinyRNA_C2epimutations_22G<-subset(tinyRNA_C2epimutations,tinyRNA_C2epimutations$Tag=="22G")
tinyRNA_C2epimutations_26G<-subset(tinyRNA_C2epimutations,tinyRNA_C2epimutations$Tag=="26G")
tinyRNA_C2epimutations_piRNA<-subset(tinyRNA_C2epimutations,tinyRNA_C2epimutations$Feature.Class=="piRNA")
tinyRNA_C2epimutations_miRNA<-subset(tinyRNA_C2epimutations,tinyRNA_C2epimutations$Feature.Class=="miRNA")
save(tinyRNA_C2epimutations_22G,file="tinyRNA_C2epimutations_22G.Rdata")

#L1
#Identify epimutations in tinyRNA using linear model in each lineage 
# Replace all values of 0 with 
tinyRNA_Lineage_L1_2_to_20<-log.norm.counts_tinyRNA_Lineage_L1[,2:10]
tinyRNA_z_scores_table_L1<-matrix(0,ncol=9, nrow = nrow(log.norm.counts_tinyRNA_Lineage_L1))

for (i in 1:ncol(tinyRNA_Lineage_L1_2_to_20)) {
  tinyRNA_EpimutL1<-
    lm(log.norm.counts_tinyRNA_Lineage_L1[,1]~tinyRNA_Lineage_L1_2_to_20[,i])
  tinyRNA_ResidL1<-tinyRNA_EpimutL1$residuals
  tinyRNA_z_scores_table_L1[,i]<-(tinyRNA_ResidL1-mean(tinyRNA_ResidL1))/sd(tinyRNA_ResidL1)
}

#keep only relevant epimutations and binarised data
tiniyRNA_binarised_z_scores_table_L1<-matrix(0,ncol=ncol(tinyRNA_z_scores_table_L1),nrow=nrow(tinyRNA_z_scores_table_L1))

for (i in 1:ncol(tinyRNA_z_scores_table_L1)) {
  for(j in 1:nrow(tinyRNA_z_scores_table_L1)){
    
    if(tinyRNA_z_scores_table_L1[j,i]> (-2.25) & tinyRNA_z_scores_table_L1[j,i]< 2.25) {
      tiniyRNA_binarised_z_scores_table_L1[j,i]<-0
    }
    if(tinyRNA_z_scores_table_L1[j,i] > 2.25){
      tiniyRNA_binarised_z_scores_table_L1[j,i]<-1
    }
    if(tinyRNA_z_scores_table_L1[j,i] < (-2.25)){
      tiniyRNA_binarised_z_scores_table_L1[j,i]<- (-1)
    }
  }  
}
tinyRNA_L1epimutations<-tiniyRNA_binarised_z_scores_table_L1 %>% as.data.frame() 
colnames(tinyRNA_L1epimutations) <- c('F2','4','6','8','12','14','16','18','20')
tinyRNA_L1epimutations <- tinyRNA_L1epimutations %>%
  # Creating an empty column:
  add_column(F0 = 0, .before='F2')
tinyRNA_L1epimutations <- tinyRNA_L1epimutations %>%
  # Creating a condition column:
  add_column(condition = "L1", .before="F0")
tinyRNA_L1epimutations <- cbind(alltiny[,1:4],alltiny[,67],tinyRNA_L1epimutations)
colnames(tinyRNA_L1epimutations) <- c('Feature.ID','Tag','Feature.Name','Feature.Class','genes', 'condition','0','2','4','6','8','12','14','16','18','20')
save(tinyRNA_L1epimutations,file="tinyRNA_L1epimutations.Rdata")

#see number of epimutations per generations
tinyRNA_L1epimutations_22G<-subset(tinyRNA_L1epimutations,tinyRNA_L1epimutations$Tag=="22G")
tinyRNA_L1epimutations_26G<-subset(tinyRNA_L1epimutations,tinyRNA_L1epimutations$Tag=="26G")
tinyRNA_L1epimutations_piRNA<-subset(tinyRNA_L1epimutations,tinyRNA_L1epimutations$Feature.Class=="piRNA")
tinyRNA_L1epimutations_miRNA<-subset(tinyRNA_L1epimutations,tinyRNA_L1epimutations$Feature.Class=="miRNA")
save(tinyRNA_L1epimutations_22G,file="tinyRNA_L1epimutations_22G.Rdata")

#L2
#Identify epimutations in tinyRNA using linear model in each lineage 
tinyRNA_Lineage_L2_2_to_20<-log.norm.counts_tinyRNA_Lineage_L2[,2:10]
tinyRNA_z_scores_table_L2<-matrix(0,ncol=9, nrow = nrow(log.norm.counts_tinyRNA_Lineage_L2))

for (i in 1:ncol(tinyRNA_Lineage_L2_2_to_20)) {
  tinyRNA_EpimutL2<-
    lm(log.norm.counts_tinyRNA_Lineage_L2[,1]~tinyRNA_Lineage_L2_2_to_20[,i])
  tinyRNA_ResidL2<-tinyRNA_EpimutL2$residuals
  tinyRNA_z_scores_table_L2[,i]<-(tinyRNA_ResidL2-mean(tinyRNA_ResidL2))/sd(tinyRNA_ResidL2)
}

#keep only relevant epimutations and binarised data
tiniyRNA_binarised_z_scores_table_L2<-matrix(0,ncol=ncol(tinyRNA_z_scores_table_L2),nrow=nrow(tinyRNA_z_scores_table_L2))

for (i in 1:ncol(tinyRNA_z_scores_table_L2)) {
  for(j in 1:nrow(tinyRNA_z_scores_table_L2)){
    
    if(tinyRNA_z_scores_table_L2[j,i]> (-2.25) & tinyRNA_z_scores_table_L2[j,i]< 2.25) {
      tiniyRNA_binarised_z_scores_table_L2[j,i]<-0
    }
    if(tinyRNA_z_scores_table_L2[j,i] > 2.25){
      tiniyRNA_binarised_z_scores_table_L2[j,i]<-1
    }
    if(tinyRNA_z_scores_table_L2[j,i] < (-2.25)){
      tiniyRNA_binarised_z_scores_table_L2[j,i]<- (-1)
    }
  }  
}
tinyRNA_L2epimutations<-tiniyRNA_binarised_z_scores_table_L2 %>% as.data.frame() 
colnames(tinyRNA_L2epimutations) <- c('F2','4','8','10','12','14','16','18','20')
tinyRNA_L2epimutations <- tinyRNA_L2epimutations %>%
  # Creating an empty column:
  add_column(F0 = 0, .before='F2')
tinyRNA_L2epimutations <- tinyRNA_L2epimutations %>%
  # Creating a condition column:
  add_column(condition = "L2", .before="F0")
tinyRNA_L2epimutations <- cbind(alltiny[,1:4],alltiny[,67],tinyRNA_L2epimutations)
colnames(tinyRNA_L2epimutations) <- c('Feature.ID','Tag','Feature.Name','Feature.Class','genes','condition','0','2','4','8','10','12','14','16','18','20')
save(tinyRNA_L2epimutations,file="tinyRNA_L2epimutations.Rdata")

#see number of epimutations per generations
tinyRNA_L2epimutations_22G<-subset(tinyRNA_L2epimutations,tinyRNA_L2epimutations$Tag=="22G")
tinyRNA_L2epimutations_26G<-subset(tinyRNA_L2epimutations,tinyRNA_L2epimutations$Tag=="26G")
tinyRNA_L2epimutations_piRNA<-subset(tinyRNA_L2epimutations,tinyRNA_L2epimutations$Feature.Class=="piRNA")
tinyRNA_L2epimutations_miRNA<-subset(tinyRNA_L2epimutations,tinyRNA_L2epimutations$Feature.Class=="miRNA")
save(tinyRNA_L2epimutations_22G,file="tinyRNA_L2epimutations_22G.Rdata")

#H1
#Identify epimutations in tinyRNA using linear model in each lineage 
tinyRNA_Lineage_H1_2_to_20<-log.norm.counts_tinyRNA_Lineage_H1[,2:10]
tinyRNA_z_scores_table_H1<-matrix(0,ncol=9, nrow = nrow(log.norm.counts_tinyRNA_Lineage_H1))

for (i in 1:ncol(tinyRNA_Lineage_H1_2_to_20)) {
  tinyRNA_EpimutH1<-
    lm(log.norm.counts_tinyRNA_Lineage_H1[,1]~tinyRNA_Lineage_H1_2_to_20[,i])
  tinyRNA_ResidH1<-tinyRNA_EpimutH1$residuals
  tinyRNA_z_scores_table_H1[,i]<-(tinyRNA_ResidH1-mean(tinyRNA_ResidH1))/sd(tinyRNA_ResidH1)
}

#keep only relevant epimutations and binarised data
tiniyRNA_binarised_z_scores_table_H1<-matrix(0,ncol=ncol(tinyRNA_z_scores_table_H1),nrow=nrow(tinyRNA_z_scores_table_H1))

for (i in 1:ncol(tinyRNA_z_scores_table_H1)) {
  for(j in 1:nrow(tinyRNA_z_scores_table_H1)){
    
    if(tinyRNA_z_scores_table_H1[j,i]> (-2.25) & tinyRNA_z_scores_table_H1[j,i]< 2.25) {
      tiniyRNA_binarised_z_scores_table_H1[j,i]<-0
    }
    if(tinyRNA_z_scores_table_H1[j,i] > 2.25){
      tiniyRNA_binarised_z_scores_table_H1[j,i]<-1
    }
    if(tinyRNA_z_scores_table_H1[j,i] < (-2.25)){
      tiniyRNA_binarised_z_scores_table_H1[j,i]<- (-1)
    }
  }  
}
tinyRNA_H1epimutations<-tiniyRNA_binarised_z_scores_table_H1 %>% as.data.frame() 
colnames(tinyRNA_H1epimutations) <- c('F2','4','8','10','12','14','16','18','20')
tinyRNA_H1epimutations <- tinyRNA_H1epimutations %>%
  # Creating an empty column:
  add_column(F0 = 0, .before='F2')
tinyRNA_H1epimutations <- tinyRNA_H1epimutations %>%
  # Creating a condition column:
  add_column(condition = "H1", .before="F0")
tinyRNA_H1epimutations <- cbind(alltiny[,1:4],alltiny[,67],tinyRNA_H1epimutations)
colnames(tinyRNA_H1epimutations) <- c('Feature.ID','Tag','Feature.Name','Feature.Class','genes','condition','0','2','4','8','10','12','14','16','18','20')
save(tinyRNA_H1epimutations,file="tinyRNA_H1epimutations.Rdata")

#see number of epimutations per generations
tinyRNA_H1epimutations_22G<-subset(tinyRNA_H1epimutations,tinyRNA_H1epimutations$Tag=="22G")
tinyRNA_H1epimutations_26G<-subset(tinyRNA_H1epimutations,tinyRNA_H1epimutations$Tag=="26G")
tinyRNA_H1epimutations_piRNA<-subset(tinyRNA_H1epimutations,tinyRNA_H1epimutations$Feature.Class=="piRNA")
tinyRNA_H1epimutations_miRNA<-subset(tinyRNA_H1epimutations,tinyRNA_H1epimutations$Feature.Class=="miRNA")
save(tinyRNA_H1epimutations_22G,file="tinyRNA_H1epimutations_22G.Rdata")

#H2
#Identify epimutations in tinyRNA using linear model in each lineage 
tinyRNA_Lineage_H2_2_to_20<-log.norm.counts_tinyRNA_Lineage_H2[,2:11]
tinyRNA_z_scores_table_H2<-matrix(0,ncol=10, nrow = nrow(log.norm.counts_tinyRNA_Lineage_H2))

for (i in 1:ncol(tinyRNA_Lineage_H2_2_to_20)) {
  tinyRNA_EpimutH2<-
    lm(log.norm.counts_tinyRNA_Lineage_H2[,1]~tinyRNA_Lineage_H2_2_to_20[,i])
  tinyRNA_ResidH2<-tinyRNA_EpimutH2$residuals
  tinyRNA_z_scores_table_H2[,i]<-(tinyRNA_ResidH2-mean(tinyRNA_ResidH2))/sd(tinyRNA_ResidH2)
}

#keep only relevant epimutations and binarised data
tiniyRNA_binarised_z_scores_table_H2<-matrix(0,ncol=ncol(tinyRNA_z_scores_table_H2),nrow=nrow(tinyRNA_z_scores_table_H2))

for (i in 1:ncol(tinyRNA_z_scores_table_H2)) {
  for(j in 1:nrow(tinyRNA_z_scores_table_H2)){
    
    if(tinyRNA_z_scores_table_H2[j,i]> (-2.25) & tinyRNA_z_scores_table_H2[j,i]< 2.25) {
      tiniyRNA_binarised_z_scores_table_H2[j,i]<-0
    }
    if(tinyRNA_z_scores_table_H2[j,i] > 2.25){
      tiniyRNA_binarised_z_scores_table_H2[j,i]<-1
    }
    if(tinyRNA_z_scores_table_H2[j,i] < (-2.25)){
      tiniyRNA_binarised_z_scores_table_H2[j,i]<- (-1)
    }
  }  
}
tinyRNA_H2epimutations<-tiniyRNA_binarised_z_scores_table_H2 %>% as.data.frame() 
colnames(tinyRNA_H2epimutations) <- c('F2','4','6','8','10','12','14','16','18','20')
tinyRNA_H2epimutations <- tinyRNA_H2epimutations %>%
  # Creating an empty column:
  add_column(F0 = 0, .before='F2')
tinyRNA_H2epimutations <- tinyRNA_H2epimutations %>%
  # Creating a condition column:
  add_column(condition = "H2", .before="F0")
tinyRNA_H2epimutations <- cbind(alltiny[,1:4],alltiny[,67],tinyRNA_H2epimutations)
colnames(tinyRNA_H2epimutations) <- c('Feature.ID','Tag','Feature.Name','Feature.Class','genes','condition','0','2','4','6','8','10','12','14','16','18','20')
save(tinyRNA_H2epimutations,file="tinyRNA_H2epimutations.Rdata")

#see number of epimutations per generations
tinyRNA_H2epimutations_22G<-subset(tinyRNA_H2epimutations,tinyRNA_H2epimutations$Tag=="22G")
tinyRNA_H2epimutations_26G<-subset(tinyRNA_H2epimutations,tinyRNA_H2epimutations$Tag=="26G")
tinyRNA_H2epimutations_piRNA<-subset(tinyRNA_H2epimutations,tinyRNA_H2epimutations$Feature.Class=="piRNA")
tinyRNA_H2epimutations_miRNA<-subset(tinyRNA_H2epimutations,tinyRNA_H2epimutations$Feature.Class=="miRNA")
save(tinyRNA_H2epimutations_22G,file="tinyRNA_H2epimutations_22G.Rdata")

#Calculate new epimutations
#Control1
tinyRNA_C1epimutations_22G<-subset(tinyRNA_C1epimutations,tinyRNA_C1epimutations$Tag=="22G")
#UP
tinyRNA_22G_UP_output_C1<- c()
tinyRNA_ID <- paste(tinyRNA_C1epimutations_22G$Feature.ID, tinyRNA_C1epimutations_22G$Tag,tinyRNA_C1epimutations_22G$Feature.Name,tinyRNA_C1epimutations_22G$Feature.Class, sep="_")
row.names(tinyRNA_C1epimutations_22G)<-tinyRNA_ID
tinyRNA_C1epimutations_22G<-tinyRNA_C1epimutations_22G[,c(7:17)]
for(i in 1:nrow(tinyRNA_C1epimutations_22G)){
  tinyRNA_22G_UP_output_C1 <-rbind(tinyRNA_22G_UP_output_C1, UP_transition_func(tinyRNA_C1epimutations_22G[i,], input_name=row.names(tinyRNA_C1epimutations_22G)[i]))}
colnames(tinyRNA_22G_UP_output_C1) <- c("0", "2", "4", "6", "8", "10", "12", "14", "16", "18", "20")
row.names(tinyRNA_22G_UP_output_C1) <- rownames(tinyRNA_C1epimutations_22G)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 
transitions_at_2 <- sum(tinyRNA_22G_UP_output_C1[,2])
transitions_at_4 <- sum(tinyRNA_22G_UP_output_C1[,3])
transitions_at_6 <- sum(tinyRNA_22G_UP_output_C1[,4])
transitions_at_8 <- sum(tinyRNA_22G_UP_output_C1[,5])
transitions_at_10 <- sum(tinyRNA_22G_UP_output_C1[,6])
transitions_at_12 <- sum(tinyRNA_22G_UP_output_C1[,7])
transitions_at_14 <- sum(tinyRNA_22G_UP_output_C1[,8])
transitions_at_16 <- sum(tinyRNA_22G_UP_output_C1[,9])
transitions_at_18 <- sum(tinyRNA_22G_UP_output_C1[,10])
transitions_at_20 <- sum(tinyRNA_22G_UP_output_C1[,11])

tinyRNA_22G_UP_C1_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                     transitions_at_4, 
                                                     transitions_at_6,
                                                     transitions_at_8, 
                                                     transitions_at_10, 
                                                     transitions_at_12, 
                                                     transitions_at_14, 
                                                     transitions_at_16, 
                                                     transitions_at_18, 
                                                     transitions_at_20)

#DOWN
tinyRNA_22G_DOWN_output_C1<- c()
for(i in 1:nrow(tinyRNA_C1epimutations_22G)){
  tinyRNA_22G_DOWN_output_C1 <-rbind(tinyRNA_22G_DOWN_output_C1, DOWN_transition_func(tinyRNA_C1epimutations_22G[i,], input_name=row.names(tinyRNA_C1epimutations_22G)[i]))}
colnames(tinyRNA_22G_DOWN_output_C1) <- c("0", "2", "4", "6", "8", "10", "12", "14", "16", "18", "20")
row.names(tinyRNA_22G_DOWN_output_C1) <- rownames(tinyRNA_C1epimutations_22G)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 
transitions_at_2 <- sum(tinyRNA_22G_DOWN_output_C1[,2])
transitions_at_4 <- sum(tinyRNA_22G_DOWN_output_C1[,3])
transitions_at_6 <- sum(tinyRNA_22G_DOWN_output_C1[,4])
transitions_at_8 <- sum(tinyRNA_22G_DOWN_output_C1[,5])
transitions_at_10 <- sum(tinyRNA_22G_DOWN_output_C1[,6])
transitions_at_12 <- sum(tinyRNA_22G_DOWN_output_C1[,7])
transitions_at_14 <- sum(tinyRNA_22G_DOWN_output_C1[,8])
transitions_at_16 <- sum(tinyRNA_22G_DOWN_output_C1[,9])
transitions_at_18 <- sum(tinyRNA_22G_DOWN_output_C1[,10])
transitions_at_20 <- sum(tinyRNA_22G_DOWN_output_C1[,11])

tinyRNA_22G_DOWN_C1_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                       transitions_at_4, 
                                                       transitions_at_6,
                                                       transitions_at_8, 
                                                       transitions_at_10, 
                                                       transitions_at_12, 
                                                       transitions_at_14, 
                                                       transitions_at_16, 
                                                       transitions_at_18, 
                                                       transitions_at_20)

#Control2
tinyRNA_C2epimutations_22G<-subset(tinyRNA_C2epimutations,tinyRNA_C2epimutations$Tag=="22G")

#UP
tinyRNA_22G_UP_output_C2<- c()
tinyRNA_ID <- paste(tinyRNA_C2epimutations_22G$Feature.ID, tinyRNA_C2epimutations_22G$Tag,tinyRNA_C2epimutations_22G$Feature.Name,tinyRNA_C2epimutations_22G$Feature.Class, sep="_")
row.names(tinyRNA_C2epimutations_22G)<-tinyRNA_ID
tinyRNA_C2epimutations_22G<-tinyRNA_C2epimutations_22G[,c(7:16)]
for(i in 1:nrow(tinyRNA_C2epimutations_22G)){
  tinyRNA_22G_UP_output_C2 <-rbind(tinyRNA_22G_UP_output_C2, UP_transition_func(tinyRNA_C2epimutations_22G[i,], input_name=row.names(tinyRNA_C2epimutations_22G)[i]))}
colnames(tinyRNA_22G_UP_output_C2) <- c("0", "2", "4", "8", "10", "12", "14", "16", "18", "20")
row.names(tinyRNA_22G_UP_output_C2) <- rownames(tinyRNA_C2epimutations_22G)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 
transitions_at_2 <- sum(tinyRNA_22G_UP_output_C2[,2])
transitions_at_4 <- sum(tinyRNA_22G_UP_output_C2[,3])
transitions_at_6 <- NA
transitions_at_8 <- sum(tinyRNA_22G_UP_output_C2[,4])
transitions_at_10 <- sum(tinyRNA_22G_UP_output_C2[,5])
transitions_at_12 <- sum(tinyRNA_22G_UP_output_C2[,6])
transitions_at_14 <- sum(tinyRNA_22G_UP_output_C2[,7])
transitions_at_16 <- sum(tinyRNA_22G_UP_output_C2[,8])
transitions_at_18 <- sum(tinyRNA_22G_UP_output_C2[,9])
transitions_at_20 <- sum(tinyRNA_22G_UP_output_C2[,10])

tinyRNA_22G_UP_C2_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                     transitions_at_4, 
                                                     transitions_at_6,
                                                     transitions_at_8, 
                                                     transitions_at_10, 
                                                     transitions_at_12, 
                                                     transitions_at_14, 
                                                     transitions_at_16, 
                                                     transitions_at_18, 
                                                     transitions_at_20)

#DOWN
tinyRNA_22G_DOWN_output_C2<- c()
for(i in 1:nrow(tinyRNA_C2epimutations_22G)){
  tinyRNA_22G_DOWN_output_C2 <-rbind(tinyRNA_22G_DOWN_output_C2, DOWN_transition_func(tinyRNA_C2epimutations_22G[i,], input_name=row.names(tinyRNA_C2epimutations_22G)[i]))}
colnames(tinyRNA_22G_DOWN_output_C2) <- c("0", "2", "4", "8", "10", "12", "14", "16", "18", "20")
row.names(tinyRNA_22G_DOWN_output_C2) <- rownames(tinyRNA_C2epimutations_22G)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 
transitions_at_2 <- sum(tinyRNA_22G_DOWN_output_C2[,2])
transitions_at_4 <- sum(tinyRNA_22G_DOWN_output_C2[,3])
transitions_at_6 <- NA
transitions_at_8 <- sum(tinyRNA_22G_DOWN_output_C2[,4])
transitions_at_10 <- sum(tinyRNA_22G_DOWN_output_C2[,5])
transitions_at_12 <- sum(tinyRNA_22G_DOWN_output_C2[,6])
transitions_at_14 <- sum(tinyRNA_22G_DOWN_output_C2[,7])
transitions_at_16 <- sum(tinyRNA_22G_DOWN_output_C2[,8])
transitions_at_18 <- sum(tinyRNA_22G_DOWN_output_C2[,9])
transitions_at_20 <- sum(tinyRNA_22G_DOWN_output_C2[,10])

tinyRNA_22G_DOWN_C2_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                       transitions_at_4, 
                                                       transitions_at_6,
                                                       transitions_at_8, 
                                                       transitions_at_10, 
                                                       transitions_at_12, 
                                                       transitions_at_14, 
                                                       transitions_at_16, 
                                                       transitions_at_18, 
                                                       transitions_at_20)


#Low1 new epi
tinyRNA_L1epimutations_22G<-subset(tinyRNA_L1epimutations,tinyRNA_L1epimutations$Tag=="22G")

#UP
tinyRNA_22G_UP_output_L1<- c()
tinyRNA_ID <- paste(tinyRNA_L1epimutations_22G$Feature.ID, tinyRNA_L1epimutations_22G$Tag,tinyRNA_L1epimutations_22G$Feature.Name,tinyRNA_L1epimutations_22G$Feature.Class, sep="_")
row.names(tinyRNA_L1epimutations_22G)<-tinyRNA_ID
tinyRNA_L1epimutations_22G<-tinyRNA_L1epimutations_22G[,c(7:16)]
for(i in 1:nrow(tinyRNA_L1epimutations_22G)){
  tinyRNA_22G_UP_output_L1 <-rbind(tinyRNA_22G_UP_output_L1, UP_transition_func(tinyRNA_L1epimutations_22G[i,], input_name=row.names(tinyRNA_L1epimutations_22G)[i]))}
colnames(tinyRNA_22G_UP_output_L1) <- c("0", "2", "4","6", "8", "12", "14", "16", "18", "20")
row.names(tinyRNA_22G_UP_output_L1) <- rownames(tinyRNA_L1epimutations_22G)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 
transitions_at_2 <- sum(tinyRNA_22G_UP_output_L1[,2])
transitions_at_4 <- sum(tinyRNA_22G_UP_output_L1[,3])
transitions_at_6 <- sum(tinyRNA_22G_UP_output_L1[,4])
transitions_at_8 <- sum(tinyRNA_22G_UP_output_L1[,5])
transitions_at_10 <- NA
transitions_at_12 <- sum(tinyRNA_22G_UP_output_L1[,6])
transitions_at_14 <- sum(tinyRNA_22G_UP_output_L1[,7])
transitions_at_16 <- sum(tinyRNA_22G_UP_output_L1[,8])
transitions_at_18 <- sum(tinyRNA_22G_UP_output_L1[,9])
transitions_at_20 <- sum(tinyRNA_22G_UP_output_L1[,10])

tinyRNA_22G_UP_L1_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                     transitions_at_4, 
                                                     transitions_at_6,
                                                     transitions_at_8, 
                                                     transitions_at_10, 
                                                     transitions_at_12, 
                                                     transitions_at_14, 
                                                     transitions_at_16, 
                                                     transitions_at_18, 
                                                     transitions_at_20)

#DOWN
tinyRNA_22G_DOWN_output_L1<- c()
for(i in 1:nrow(tinyRNA_L1epimutations_22G)){
  tinyRNA_22G_DOWN_output_L1 <-rbind(tinyRNA_22G_DOWN_output_L1, DOWN_transition_func(tinyRNA_L1epimutations_22G[i,], input_name=row.names(tinyRNA_L1epimutations_22G)[i]))}
colnames(tinyRNA_22G_DOWN_output_L1) <- c("0", "2", "4", "6","8", "12", "14", "16", "18", "20")
row.names(tinyRNA_22G_DOWN_output_L1) <- rownames(tinyRNA_L1epimutations_22G)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 
transitions_at_2 <- sum(tinyRNA_22G_DOWN_output_L1[,2])
transitions_at_4 <- sum(tinyRNA_22G_DOWN_output_L1[,3])
transitions_at_6 <- sum(tinyRNA_22G_DOWN_output_L1[,4])
transitions_at_8 <- sum(tinyRNA_22G_DOWN_output_L1[,5])
transitions_at_10 <- NA
transitions_at_12 <- sum(tinyRNA_22G_DOWN_output_L1[,6])
transitions_at_14 <- sum(tinyRNA_22G_DOWN_output_L1[,7])
transitions_at_16 <- sum(tinyRNA_22G_DOWN_output_L1[,8])
transitions_at_18 <- sum(tinyRNA_22G_DOWN_output_L1[,9])
transitions_at_20 <- sum(tinyRNA_22G_DOWN_output_L1[,10])

tinyRNA_22G_DOWN_L1_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                       transitions_at_4, 
                                                       transitions_at_6,
                                                       transitions_at_8, 
                                                       transitions_at_10, 
                                                       transitions_at_12, 
                                                       transitions_at_14, 
                                                       transitions_at_16, 
                                                       transitions_at_18, 
                                                       transitions_at_20)

#Low2 new epi
#UP
tinyRNA_22G_UP_output_L2<- c()
tinyRNA_ID <- paste(tinyRNA_L2epimutations_22G$Feature.ID, tinyRNA_L2epimutations_22G$Tag,tinyRNA_L2epimutations_22G$Feature.Name,tinyRNA_L2epimutations_22G$Feature.Class, sep="_")
row.names(tinyRNA_L2epimutations_22G)<-tinyRNA_ID
tinyRNA_L2epimutations_22G<-tinyRNA_L2epimutations_22G[,c(7:16)]
for(i in 1:nrow(tinyRNA_L2epimutations_22G)){
  tinyRNA_22G_UP_output_L2 <-rbind(tinyRNA_22G_UP_output_L2, UP_transition_func(tinyRNA_L2epimutations_22G[i,], input_name=row.names(tinyRNA_L2epimutations_22G)[i]))}
colnames(tinyRNA_22G_UP_output_L2) <- c("0", "2", "4", "8","10", "12", "14", "16", "18", "20")
row.names(tinyRNA_22G_UP_output_L2) <- rownames(tinyRNA_L2epimutations_22G)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 
transitions_at_2 <- sum(tinyRNA_22G_UP_output_L2[,2])
transitions_at_4 <- sum(tinyRNA_22G_UP_output_L2[,3])
transitions_at_6 <- NA
transitions_at_8 <- sum(tinyRNA_22G_UP_output_L2[,4])
transitions_at_10 <- sum(tinyRNA_22G_UP_output_L2[,5])
transitions_at_12 <- sum(tinyRNA_22G_UP_output_L2[,6])
transitions_at_14 <- sum(tinyRNA_22G_UP_output_L2[,7])
transitions_at_16 <- sum(tinyRNA_22G_UP_output_L2[,8])
transitions_at_18 <- sum(tinyRNA_22G_UP_output_L2[,9])
transitions_at_20 <- sum(tinyRNA_22G_UP_output_L2[,10])

tinyRNA_22G_UP_L2_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                     transitions_at_4, 
                                                     transitions_at_6,
                                                     transitions_at_8, 
                                                     transitions_at_10, 
                                                     transitions_at_12, 
                                                     transitions_at_14, 
                                                     transitions_at_16, 
                                                     transitions_at_18, 
                                                     transitions_at_20)


#DOWN
tinyRNA_22G_DOWN_output_L2<- c()
for(i in 1:nrow(tinyRNA_L2epimutations_22G)){
  tinyRNA_22G_DOWN_output_L2 <-rbind(tinyRNA_22G_DOWN_output_L2, DOWN_transition_func(tinyRNA_L2epimutations_22G[i,], input_name=row.names(tinyRNA_L2epimutations_22G)[i]))}
colnames(tinyRNA_22G_DOWN_output_L2) <- c("0", "2", "4", "6","8", "12", "14", "16", "18", "20")
row.names(tinyRNA_22G_DOWN_output_L2) <- rownames(tinyRNA_L2epimutations_22G)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 
transitions_at_2 <- sum(tinyRNA_22G_DOWN_output_L2[,2])
transitions_at_4 <- sum(tinyRNA_22G_DOWN_output_L2[,3])
transitions_at_6 <- NA
transitions_at_8 <- sum(tinyRNA_22G_DOWN_output_L2[,4])
transitions_at_10 <- sum(tinyRNA_22G_DOWN_output_L2[,5])
transitions_at_12 <- sum(tinyRNA_22G_DOWN_output_L2[,6])
transitions_at_14 <- sum(tinyRNA_22G_DOWN_output_L2[,7])
transitions_at_16 <- sum(tinyRNA_22G_DOWN_output_L2[,8])
transitions_at_18 <- sum(tinyRNA_22G_DOWN_output_L2[,9])
transitions_at_20 <- sum(tinyRNA_22G_DOWN_output_L2[,10])

tinyRNA_22G_DOWN_L2_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                       transitions_at_4, 
                                                       transitions_at_6,
                                                       transitions_at_8, 
                                                       transitions_at_10, 
                                                       transitions_at_12, 
                                                       transitions_at_14, 
                                                       transitions_at_16, 
                                                       transitions_at_18, 
                                                       transitions_at_20)

#High1 new epi
tinyRNA_H1epimutations_22G<-subset(tinyRNA_H1epimutations,tinyRNA_H1epimutations$Tag=="22G")

#UP
tinyRNA_22G_UP_output_H1<- c()
tinyRNA_ID <- paste(tinyRNA_H1epimutations_22G$Feature.ID, tinyRNA_H1epimutations_22G$Tag,tinyRNA_H1epimutations_22G$Feature.Name,tinyRNA_H1epimutations_22G$Feature.Class, sep="_")
row.names(tinyRNA_H1epimutations_22G)<-tinyRNA_ID
tinyRNA_H1epimutations_22G<-tinyRNA_H1epimutations_22G[,c(7:16)]
for(i in 1:nrow(tinyRNA_H1epimutations_22G)){
  tinyRNA_22G_UP_output_H1 <-rbind(tinyRNA_22G_UP_output_H1, UP_transition_func(tinyRNA_H1epimutations_22G[i,], input_name=row.names(tinyRNA_H1epimutations_22G)[i]))}
colnames(tinyRNA_22G_UP_output_H1) <- c("0", "2", "4", "8","10", "12", "14", "16", "18", "20")
row.names(tinyRNA_22G_UP_output_H1) <- rownames(tinyRNA_H1epimutations_22G)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 
transitions_at_2 <- sum(tinyRNA_22G_UP_output_H1[,2])
transitions_at_4 <- sum(tinyRNA_22G_UP_output_H1[,3])
transitions_at_6 <- NA
transitions_at_8 <- sum(tinyRNA_22G_UP_output_H1[,4])
transitions_at_10 <- sum(tinyRNA_22G_UP_output_H1[,5])
transitions_at_12 <- sum(tinyRNA_22G_UP_output_H1[,6])
transitions_at_14 <- sum(tinyRNA_22G_UP_output_H1[,7])
transitions_at_16 <- sum(tinyRNA_22G_UP_output_H1[,8])
transitions_at_18 <- sum(tinyRNA_22G_UP_output_H1[,9])
transitions_at_20 <- sum(tinyRNA_22G_UP_output_H1[,10])

tinyRNA_22G_UP_H1_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                     transitions_at_4, 
                                                     transitions_at_6,
                                                     transitions_at_8, 
                                                     transitions_at_10, 
                                                     transitions_at_12, 
                                                     transitions_at_14, 
                                                     transitions_at_16, 
                                                     transitions_at_18, 
                                                     transitions_at_20)


#DOWN
tinyRNA_22G_DOWN_output_H1<- c()
for(i in 1:nrow(tinyRNA_H1epimutations_22G)){
  tinyRNA_22G_DOWN_output_H1 <-rbind(tinyRNA_22G_DOWN_output_H1, DOWN_transition_func(tinyRNA_H1epimutations_22G[i,], input_name=row.names(tinyRNA_H1epimutations_22G)[i]))}
colnames(tinyRNA_22G_DOWN_output_H1) <- c("0", "2", "4", "8", "10", "12", "14", "16", "18", "20")
row.names(tinyRNA_22G_DOWN_output_H1) <- rownames(tinyRNA_H1epimutations_22G)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 
transitions_at_2 <- sum(tinyRNA_22G_DOWN_output_H1[,2])
transitions_at_4 <- sum(tinyRNA_22G_DOWN_output_H1[,3])
transitions_at_6 <- NA
transitions_at_8 <- sum(tinyRNA_22G_DOWN_output_H1[,4])
transitions_at_10 <- sum(tinyRNA_22G_DOWN_output_H1[,5])
transitions_at_12 <- sum(tinyRNA_22G_DOWN_output_H1[,6])
transitions_at_14 <- sum(tinyRNA_22G_DOWN_output_H1[,7])
transitions_at_16 <- sum(tinyRNA_22G_DOWN_output_H1[,8])
transitions_at_18 <- sum(tinyRNA_22G_DOWN_output_H1[,9])
transitions_at_20 <- sum(tinyRNA_22G_DOWN_output_H1[,10])

tinyRNA_22G_DOWN_H1_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                       transitions_at_4, 
                                                       transitions_at_6,
                                                       transitions_at_8, 
                                                       transitions_at_10, 
                                                       transitions_at_12, 
                                                       transitions_at_14, 
                                                       transitions_at_16, 
                                                       transitions_at_18, 
                                                       transitions_at_20)

#High2
tinyRNA_H2epimutations_22G<-subset(tinyRNA_H2epimutations,tinyRNA_H2epimutations$Tag=="22G")
#UP

tinyRNA_22G_UP_output_H2<- c()
tinyRNA_ID <- paste(tinyRNA_H2epimutations_22G$Feature.ID, tinyRNA_H2epimutations_22G$Tag,tinyRNA_H2epimutations_22G$Feature.Name,tinyRNA_H2epimutations_22G$Feature.Class, sep="_")
row.names(tinyRNA_H2epimutations_22G)<-tinyRNA_ID
tinyRNA_H2epimutations_22G<-tinyRNA_H2epimutations_22G[,c(7:17)]
for(i in 1:nrow(tinyRNA_H2epimutations_22G)){
  tinyRNA_22G_UP_output_H2 <-rbind(tinyRNA_22G_UP_output_H2, UP_transition_func(tinyRNA_H2epimutations_22G[i,], input_name=row.names(tinyRNA_H2epimutations_22G)[i]))}
colnames(tinyRNA_22G_UP_output_H2) <- c("0", "2", "4", "6", "8","10", "12", "14", "16", "18", "20")
row.names(tinyRNA_22G_UP_output_H2) <- rownames(tinyRNA_H2epimutations_22G)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tinyRNA_22G_UP_output_H2[,2])
transitions_at_4 <- sum(tinyRNA_22G_UP_output_H2[,3])
transitions_at_6 <- sum(tinyRNA_22G_UP_output_H2[,4])
transitions_at_8 <- sum(tinyRNA_22G_UP_output_H2[,5])
transitions_at_10 <- sum(tinyRNA_22G_UP_output_H2[,6])
transitions_at_12 <- sum(tinyRNA_22G_UP_output_H2[,7])
transitions_at_14 <- sum(tinyRNA_22G_UP_output_H2[,8])
transitions_at_16 <- sum(tinyRNA_22G_UP_output_H2[,9])
transitions_at_18 <- sum(tinyRNA_22G_UP_output_H2[,10])
transitions_at_20 <- sum(tinyRNA_22G_UP_output_H2[,11])

tinyRNA_22G_UP_H2_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                     transitions_at_4, 
                                                     transitions_at_6,
                                                     transitions_at_8, 
                                                     transitions_at_10, 
                                                     transitions_at_12, 
                                                     transitions_at_14, 
                                                     transitions_at_16, 
                                                     transitions_at_18, 
                                                     transitions_at_20)


#DOWN
tinyRNA_22G_DOWN_output_H2<- c()
for(i in 1:nrow(tinyRNA_H2epimutations_22G)){
  tinyRNA_22G_DOWN_output_H2 <-rbind(tinyRNA_22G_DOWN_output_H2, DOWN_transition_func(tinyRNA_H2epimutations_22G[i,], input_name=row.names(tinyRNA_H2epimutations_22G)[i]))}
colnames(tinyRNA_22G_DOWN_output_H2) <- c("0", "2", "4", "6", "8", "10", "12", "14", "16", "18", "20")
row.names(tinyRNA_22G_DOWN_output_H2) <- rownames(tinyRNA_H2epimutations_22G)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tinyRNA_22G_DOWN_output_H2[,2])
transitions_at_4 <- sum(tinyRNA_22G_DOWN_output_H2[,3])
transitions_at_6 <- sum(tinyRNA_22G_DOWN_output_H2[,4])
transitions_at_8 <- sum(tinyRNA_22G_DOWN_output_H2[,5])
transitions_at_10 <- sum(tinyRNA_22G_DOWN_output_H2[,6])
transitions_at_12 <- sum(tinyRNA_22G_DOWN_output_H2[,7])
transitions_at_14 <- sum(tinyRNA_22G_DOWN_output_H2[,8])
transitions_at_16 <- sum(tinyRNA_22G_DOWN_output_H2[,9])
transitions_at_18 <- sum(tinyRNA_22G_DOWN_output_H2[,10])
transitions_at_20 <- sum(tinyRNA_22G_DOWN_output_H2[,11])

tinyRNA_22G_DOWN_H2_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                       transitions_at_4, 
                                                       transitions_at_6,
                                                       transitions_at_8, 
                                                       transitions_at_10, 
                                                       transitions_at_12, 
                                                       transitions_at_14, 
                                                       transitions_at_16, 
                                                       transitions_at_18, 
                                                       transitions_at_20)

#Plot construction - Fig.4.A
names <- rownames(tinyRNA_22G_UP_C1_Table_of_new_epimutations)
rownames(tinyRNA_22G_UP_C1_Table_of_new_epimutations) <- NULL
tinyRNA_22G_UP_C1_Table_of_new_epimutations <- cbind(names,tinyRNA_22G_UP_C1_Table_of_new_epimutations)
colnames(tinyRNA_22G_UP_C1_Table_of_new_epimutations) <- c("Transitions","Up")
names <- rownames(tinyRNA_22G_DOWN_C1_Table_of_new_epimutations)
rownames(tinyRNA_22G_DOWN_C1_Table_of_new_epimutations) <- NULL
tinyRNA_22G_DOWN_C1_Table_of_new_epimutations <- cbind(names,tinyRNA_22G_DOWN_C1_Table_of_new_epimutations)
colnames(tinyRNA_22G_DOWN_C1_Table_of_new_epimutations) <- c("Transitions","Down")

tinyRNA_22G_new_epimut_C1 <- merge(tinyRNA_22G_UP_C1_Table_of_new_epimutations,tinyRNA_22G_DOWN_C1_Table_of_new_epimutations)
tinyRNA_22G_new_epimut_C1 <- tinyRNA_22G_new_epimut_C1 %>%
  # Creating an empty column:
  add_column(Lineage = "C1", .after="Transitions")
names <- rownames(tinyRNA_22G_UP_C2_Table_of_new_epimutations)
rownames(tinyRNA_22G_UP_C2_Table_of_new_epimutations) <- NULL
tinyRNA_22G_UP_C2_Table_of_new_epimutations <- cbind(names,tinyRNA_22G_UP_C2_Table_of_new_epimutations)
colnames(tinyRNA_22G_UP_C2_Table_of_new_epimutations) <- c("Transitions","Up")
names <- rownames(tinyRNA_22G_DOWN_C2_Table_of_new_epimutations)
rownames(tinyRNA_22G_DOWN_C2_Table_of_new_epimutations) <- NULL
tinyRNA_22G_DOWN_C2_Table_of_new_epimutations <- cbind(names,tinyRNA_22G_DOWN_C2_Table_of_new_epimutations)
colnames(tinyRNA_22G_DOWN_C2_Table_of_new_epimutations) <- c("Transitions","Down")

tinyRNA_22G_new_epimut_C2 <- merge(tinyRNA_22G_UP_C2_Table_of_new_epimutations,tinyRNA_22G_DOWN_C2_Table_of_new_epimutations)
tinyRNA_22G_new_epimut_C2 <- tinyRNA_22G_new_epimut_C2 %>%
  # Creating an empty column:
  add_column(Lineage = "C2", .after="Transitions")

tiny_22G_C <- rbind(tinyRNA_22G_new_epimut_C1, tinyRNA_22G_new_epimut_C2)
tiny_22G_C <- tiny_22G_C %>%
  # Creating an empty column:
  add_column(Condition = "Control", .after="Transitions")
names <- rownames(tinyRNA_22G_UP_L1_Table_of_new_epimutations)
rownames(tinyRNA_22G_UP_L1_Table_of_new_epimutations) <- NULL
tinyRNA_22G_UP_L1_Table_of_new_epimutations <- cbind(names,tinyRNA_22G_UP_L1_Table_of_new_epimutations)
colnames(tinyRNA_22G_UP_L1_Table_of_new_epimutations) <- c("Transitions","Up")
names <- rownames(tinyRNA_22G_DOWN_L1_Table_of_new_epimutations)
rownames(tinyRNA_22G_DOWN_L1_Table_of_new_epimutations) <- NULL
tinyRNA_22G_DOWN_L1_Table_of_new_epimutations <- cbind(names,tinyRNA_22G_DOWN_L1_Table_of_new_epimutations)
colnames(tinyRNA_22G_DOWN_L1_Table_of_new_epimutations) <- c("Transitions","Down")

tinyRNA_22G_new_epimut_L1 <- merge(tinyRNA_22G_UP_L1_Table_of_new_epimutations,tinyRNA_22G_DOWN_L1_Table_of_new_epimutations)
tinyRNA_22G_new_epimut_L1 <- tinyRNA_22G_new_epimut_L1 %>%
  # Creating an empty column:
  add_column(Lineage = "L1", .after="Transitions")
names <- rownames(tinyRNA_22G_UP_L2_Table_of_new_epimutations)
rownames(tinyRNA_22G_UP_L2_Table_of_new_epimutations) <- NULL
tinyRNA_22G_UP_L2_Table_of_new_epimutations <- cbind(names,tinyRNA_22G_UP_L2_Table_of_new_epimutations)
colnames(tinyRNA_22G_UP_L2_Table_of_new_epimutations) <- c("Transitions","Up")
names <- rownames(tinyRNA_22G_DOWN_L2_Table_of_new_epimutations)
rownames(tinyRNA_22G_DOWN_L2_Table_of_new_epimutations) <- NULL
tinyRNA_22G_DOWN_L2_Table_of_new_epimutations <- cbind(names,tinyRNA_22G_DOWN_L2_Table_of_new_epimutations)
colnames(tinyRNA_22G_DOWN_L2_Table_of_new_epimutations) <- c("Transitions","Down")

tinyRNA_22G_new_epimut_L2 <- merge(tinyRNA_22G_UP_L2_Table_of_new_epimutations,tinyRNA_22G_DOWN_L2_Table_of_new_epimutations)
tinyRNA_22G_new_epimut_L2 <- tinyRNA_22G_new_epimut_L2 %>%
  # Creating an empty column:
  add_column(Lineage = "L2", .after="Transitions")

tiny_22G_L <- rbind(tinyRNA_22G_new_epimut_L1, tinyRNA_22G_new_epimut_L2)
tiny_22G_L <- tiny_22G_L %>%
  # Creating an empty column:
  add_column(Condition = "Low dose", .after="Transitions")
names <- rownames(tinyRNA_22G_UP_H1_Table_of_new_epimutations)
rownames(tinyRNA_22G_UP_H1_Table_of_new_epimutations) <- NULL
tinyRNA_22G_UP_H1_Table_of_new_epimutations <- cbind(names,tinyRNA_22G_UP_H1_Table_of_new_epimutations)
colnames(tinyRNA_22G_UP_H1_Table_of_new_epimutations) <- c("Transitions","Up")
names <- rownames(tinyRNA_22G_DOWN_H1_Table_of_new_epimutations)
rownames(tinyRNA_22G_DOWN_H1_Table_of_new_epimutations) <- NULL
tinyRNA_22G_DOWN_H1_Table_of_new_epimutations <- cbind(names,tinyRNA_22G_DOWN_H1_Table_of_new_epimutations)
colnames(tinyRNA_22G_DOWN_H1_Table_of_new_epimutations) <- c("Transitions","Down")

tinyRNA_22G_new_epimut_H1 <- merge(tinyRNA_22G_UP_H1_Table_of_new_epimutations,tinyRNA_22G_DOWN_H1_Table_of_new_epimutations)
tinyRNA_22G_new_epimut_H1 <- tinyRNA_22G_new_epimut_H1 %>%
  # Creating an empty column:
  add_column(Lineage = "H1", .after="Transitions")
names <- rownames(tinyRNA_22G_UP_H2_Table_of_new_epimutations)
rownames(tinyRNA_22G_UP_H2_Table_of_new_epimutations) <- NULL
tinyRNA_22G_UP_H2_Table_of_new_epimutations <- cbind(names,tinyRNA_22G_UP_H2_Table_of_new_epimutations)
colnames(tinyRNA_22G_UP_H2_Table_of_new_epimutations) <- c("Transitions","Up")
names <- rownames(tinyRNA_22G_DOWN_H2_Table_of_new_epimutations)
rownames(tinyRNA_22G_DOWN_H2_Table_of_new_epimutations) <- NULL
tinyRNA_22G_DOWN_H2_Table_of_new_epimutations <- cbind(names,tinyRNA_22G_DOWN_H2_Table_of_new_epimutations)
colnames(tinyRNA_22G_DOWN_H2_Table_of_new_epimutations) <- c("Transitions","Down")

tinyRNA_22G_new_epimut_H2 <- merge(tinyRNA_22G_UP_H2_Table_of_new_epimutations,tinyRNA_22G_DOWN_H2_Table_of_new_epimutations)
tinyRNA_22G_new_epimut_H2 <- tinyRNA_22G_new_epimut_H2 %>%
  # Creating an empty column:
  add_column(Lineage = "H2", .after="Transitions")

tiny_22G_H <- rbind(tinyRNA_22G_new_epimut_H1, tinyRNA_22G_new_epimut_H2)
tiny_22G_H <- tiny_22G_H %>%
  # Creating an empty column:
  add_column(Condition = "High dose", .after="Transitions")

tiny_22G_newepi<-rbind(tiny_22G_C,tiny_22G_L,tiny_22G_H)
tiny_22G_newepi$Up <- as.numeric(tiny_22G_newepi$Up)
tiny_22G_newepi$Down <- as.numeric(tiny_22G_newepi$Down)
tiny_22G_newepi$Total=rowSums(cbind(tiny_22G_newepi$Up,tiny_22G_newepi$Down),na.rm=FALSE)

UP<-tiny_22G_newepi$Up
DOWN<-tiny_22G_newepi$Down

#Test for statistical significance
ggdensity(tiny_22G_newepi$Total, 
          main = "Density plot of epimutations_number",
          xlab = "Total number of epimutations_number")
shapiro.test(tiny_22G_newepi$Total)
kruskal.test(Total ~ Condition, data = tiny_22G_newepi)
dunnTest(Total ~ Condition, data = tiny_22G_newepi)

#Fig.4.A
tiny_22G_newepi_detailed$Generation <- fct_relevel(tiny_22G_newepi_detailed$Generation, c("2", "4", "6","8","10","12","14","16","18","20"))
tiny_22G_allLineage_rate <- ggplot(tiny_22G_newepi, aes(x=Condition, y=Total, color = Condition,group=Condition)) + 
  geom_boxplot(fatten = 1, lwd = 1, width=0.5)+
  scale_color_manual(values=c("cornflowerblue", "darkgreen", "red"))+
  labs(y = "Number of new epimutations", x = "\nCondition")+
  geom_dotplot(binaxis='y', binwidth = 0.25, stackdir='center',  dotsize=75, binpositions = "all", stackgroups = TRUE)+
  theme_bw()+
  theme_linedraw()+
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(color = "black", size = 25, angle = 0, face = "plain"),
        axis.text.y = element_text(color = "black", size = 25, angle = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 20, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, face = "bold"))+
  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))+  
  scale_fill_brewer(palette="Pastel2")+
  ggtitle(paste("tiny_22G_allLineage_rate"))

#Sup.Fig.5.A
Generation<-c("10","12","14","16","18","2","20","4","6","8",
              "10","12","14","16","18","2","20","4","6","8",
              "10","12","14","16","18","2","20","4","6","8",
              "10","12","14","16","18","2","20","4","6","8",
              "10","12","14","16","18","2","20","4","6","8",
              "10","12","14","16","18","2","20","4","6","8")
tiny_22G_newepi_detailed<-cbind(tiny_22G_newepi,Generation)
tiny_22G_newepi_detailed$Generation <- fct_relevel(tiny_22G_newepi_detailed$Generation, c("2", "4", "6","8","10","12","14","16","18","20"))
tiny_22G_newepi_detailed$Condition <- fct_relevel(tiny_22G_newepi_detailed$Condition, c("Control", "Low dose", "High dose"))
tiny_22G_newepi_detailed_plot <- ggplot(tiny_22G_newepi_detailed, aes(x=Generation , y=Total, color = Condition)) + 
  geom_boxplot(fatten = 1, lwd = 1, width=0.5)+
  scale_color_manual(values=c("cornflowerblue", "darkgreen", "red"))+
  labs(y = "Number of new epimutations", x = "\nGeneration")+
  geom_jitter(aes(shape=Lineage, size=8),
              position=position_jitter(width = 0.2,
                                       height = 0.2))+
  theme_bw()+
  theme(axis.text.x = element_text(color = "black", size = 20, angle = 0, face = "plain"),
        axis.text.y = element_text(color = "black", size = 20, angle = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 16, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 16, angle = 90, face = "bold"))+
  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))+  
  scale_fill_brewer(palette="Pastel2")+
  ggtitle(paste(""))
tiny_22G_newepi_detailed_plot 

#Get raw data
write.xlsx(tiny_22G_newepi_detailed,"Data_Fig_4_A_&_Sup_Fig_5_A.xlsx")
#--------------------------
#####Fig.4.B & Sup.Fig.5.B - 22G-RNAs epimutation duration#####
#Calculation duration epimutations sncRNAs
#Data preparation
#Remplacement of missing columns with data following generation
row.names(tinyRNA_C1epimutations)<-paste(tinyRNA_C1epimutations[,1], tinyRNA_C1epimutations[,2], tinyRNA_C1epimutations[,3], tinyRNA_C1epimutations[,4], sep = ":")
tinyRNA_C1epimutations <- tinyRNA_C1epimutations[,7:17]

tinyRNA_C2epimutationsbis<-cbind(tinyRNA_C2epimutations, rep(tinyRNA_C2epimutations[9],1))
tinyRNA_C2epimutationsbis<-tinyRNA_C2epimutationsbis[,c(1,2,3,4,5,6,7,8,17,9,10,11,12,13,14,15,16)]
colnames(tinyRNA_C2epimutationsbis) <- c('Feature.ID','Tag','Feature.Name','Feature.Class','genes','condition','0','2','4','6','8','10','12','14','16','18','20')
row.names(tinyRNA_C2epimutationsbis)<-paste(tinyRNA_C2epimutationsbis[,1], tinyRNA_C2epimutationsbis[,2], tinyRNA_C2epimutationsbis[,3], tinyRNA_C2epimutationsbis[,4], sep = ":")
tinyRNA_C2epimutationsbis <- tinyRNA_C2epimutationsbis[,7:17]

tinyRNA_L1epimutationsbis<-cbind(tinyRNA_L1epimutations, rep(tinyRNA_L1epimutations[11],1))
tinyRNA_L1epimutationsbis<-tinyRNA_L1epimutationsbis[,c(1,2,3,4,5,6,7,8,9,10,17,11,12,13,14,15,16)]
colnames(tinyRNA_L1epimutationsbis) <- c('Feature.ID','Tag','Feature.Name','Feature.Class','genes','condition','0','2','4','6','8','10','12','14','16','18','20')
row.names(tinyRNA_L1epimutationsbis)<-paste(tinyRNA_L1epimutationsbis[,1], tinyRNA_L1epimutationsbis[,2], tinyRNA_L1epimutationsbis[,3], tinyRNA_L1epimutationsbis[,4], sep = ":")
tinyRNA_L1epimutationsbis <- tinyRNA_L1epimutationsbis[,7:17]

tinyRNA_L2epimutationsbis<-cbind(tinyRNA_L2epimutations, rep(tinyRNA_L2epimutations[9],1))
tinyRNA_L2epimutationsbis<-tinyRNA_L2epimutationsbis[,c(1,2,3,4,5,6,7,8,17,9,10,11,12,13,14,15,16)]
colnames(tinyRNA_L2epimutationsbis) <- c('Feature.ID','Tag','Feature.Name','Feature.Class','genes','condition','0','2','4','6','8','10','12','14','16','18','20')
row.names(tinyRNA_L2epimutationsbis)<-paste(tinyRNA_L2epimutationsbis[,1], tinyRNA_L2epimutationsbis[,2], tinyRNA_L2epimutationsbis[,3], tinyRNA_L2epimutationsbis[,4], sep = ":")
tinyRNA_L2epimutationsbis <- tinyRNA_L2epimutationsbis[,7:17]

tinyRNA_H1epimutationsbis<-cbind(tinyRNA_H1epimutations, rep(tinyRNA_H1epimutations[9],1))
tinyRNA_H1epimutationsbis<-tinyRNA_H1epimutationsbis[,c(1,2,3,4,5,6,7,8,17,9,10,11,12,13,14,15,16)]
colnames(tinyRNA_H1epimutationsbis) <- c('Feature.ID','Tag','Feature.Name','Feature.Class','genes','condition','0','2','4','6','8','10','12','14','16','18','20')
row.names(tinyRNA_H1epimutationsbis)<-paste(tinyRNA_H1epimutationsbis[,1], tinyRNA_H1epimutationsbis[,2], tinyRNA_H1epimutationsbis[,3], tinyRNA_H1epimutationsbis[,4], sep = ":")
tinyRNA_H1epimutationsbis <- tinyRNA_H1epimutationsbis[,7:17]

row.names(tinyRNA_H2epimutations)<-paste(tinyRNA_H2epimutations[,1], tinyRNA_H2epimutations[,2], tinyRNA_H2epimutations[,3], tinyRNA_H2epimutations[,4], sep = ":")
tinyRNA_H2epimutations <- tinyRNA_H2epimutations[,7:17]

tiny_all_epi<-list(tinyRNA_C1epimutations,tinyRNA_C2epimutationsbis,tinyRNA_L1epimutationsbis,tinyRNA_L2epimutationsbis,tinyRNA_H1epimutationsbis,tinyRNA_H2epimutations)
names(tiny_all_epi)<-c("C1","C2","L1","L2","H1","H2")

#Calculation
tinyRNA_All_output_list_with_missing_as_muts <- list()

for(e in 1:length(tiny_all_epi)){
  
  main_frame <- as.data.frame(tiny_all_epi[[e]])
  
  Overall_output <- c()
  
  for(t in 1:nrow(main_frame)){
    
    vector_in <- main_frame[t, ]
    gen_names <- as.numeric(colnames(vector_in))
    
    number_transitions <- 0
    onset_gen <- 0
    tempL <- 0
    output_frame <- c()
    
    is_up <- 0
    is_down <- 0
    complete <- 0
    
    # first deal with no epimutations in the lineage for a locus/gene
    
    if(sum(abs(vector_in)) == 0){
      
      length = 0
      is_up <- 0
      is_down <- 0
      gene_name <- rownames(vector_in)
      gene <- strsplit(gene_name, ":")[[1]][4]
      
      number_transitions <- 0 
      onset_gen <- 0
      
      Lineage <- names(tiny_all_epi[e])
      complete <- 0
      
      save_mut <- data.frame(gene_name, gene, number_transitions, onset_gen, length,  complete,  Lineage, is_up, is_down)   
      
      output_frame <- rbind(output_frame, save_mut)
      
    }
    
    else
      
      if(sum(abs(vector_in)) > 0){
        
        for(i in 2:length(vector_in)){
          
          
          # if it is an UP transition, turn off any down transitions and start new epimutation
          
          if(vector_in[i]==1&vector_in[i-1]== 0|vector_in[i]== 1&vector_in[i-1]== -1) {
            
            if(is_down == 1){ # a down epimutation that is now turning off
              
              gene_name <- rownames(vector_in)
              gene <- strsplit(gene_name, ":")[[1]][4]
              
              number_transitions <- number_transitions + 1 
              
              length <-tempL
              
              Lineage <- names(tiny_all_epi[e])
              
              complete <- 1
              
              save_mut <- data.frame(gene_name, gene, number_transitions, onset_gen, length,  complete, Lineage, is_up, is_down)   
              
              output_frame <- rbind(output_frame, save_mut)
              
            }
            
            # start the UP transition
            
            is_up <- 1
            is_down <- 0
            tempL <- 1
            
            number_transitions <- 1 # there is a transition to an UP epimutation at generation[i]
            onset_gen <- gen_names[i]
            
            onset_check <- i -1
          }
          
          # if it is a new DOWN transition, turn off any up transitions and start new epimutation
          
          if(vector_in[i]==-1&vector_in[i-1]== 0|vector_in[i]== -1&vector_in[i-1]== 1) {
            
            if(is_up == 1){ # an up epimutation that is now turning off
              
              gene_name <- rownames(vector_in)
              gene <- strsplit(gene_name, ":")[[1]][4]
              
              number_transitions <- number_transitions + 1
              
              length <- tempL
              
              Lineage <- names(tiny_all_epi[e])
              
              complete <- 1
              
              save_mut <- data.frame(gene_name, gene, number_transitions, onset_gen, length,  complete, Lineage, is_up, is_down)   
              
              output_frame <- rbind(output_frame, save_mut)
              
            }
            
            # start the down transition
            
            is_down <- 1
            is_up <- 0
            tempL <- 1
            number_transitions <- 1 # there is a transition to a DOWN epimutation at generation[i]
            onset_gen <- gen_names[i]
            onset_check <- i-1
          }
          
          # if it is an UP epimutation, does it continue, does it turn off, does it extend to last gen?
          
          if(is_up==1){
            
            
            if(vector_in[i]==1&vector_in[i-1]==1)
              
            {tempL<-tempL+(gen_names[i]-gen_names[i-1])}  # the epimutation continues
            
            if(vector_in[i]==0&vector_in[i-1]==1){
              
              gene_name <- rownames(vector_in)
              gene <- strsplit(gene_name, ":")[[1]][4]
              
              number_transitions <- number_transitions +1  # there is a transition to OFF from an UP epimutation at generation[i]
              
              length <- tempL
              
              Lineage <- names(tiny_all_epi[e])
              
              complete <- 1
              
              save_mut <- data.frame(gene_name, gene, number_transitions, onset_gen, length,  complete, Lineage, is_up, is_down)   
              
              output_frame <- rbind(output_frame, save_mut)
              
              is_up <- 0
              is_down <- 0
              tempL <- 0
              length <- 0
              complete <- 0
              
            }
            
            # deal with an up epimutation extending to end of lineage
            
            if(vector_in[i]==1&i==length(vector_in)){
              
              gene_name <- rownames(vector_in)
              
              gene <- strsplit(gene_name, ":")[[1]][4]
              
              number_transitions <- number_transitions + 0
              
              length <- tempL
              
              Lineage <- names(tiny_all_epi[e])
              
              complete <- 0
              
              save_mut <- data.frame(gene_name, gene, number_transitions, onset_gen, length,  complete, Lineage, is_up, is_down)   
              
              output_frame <- rbind(output_frame, save_mut)
              
              is_up <- 0
              is_down <- 0
              length <- 0
              tempL <- 0
              complete <- 0
              
            }
            
          }
          
          # if it is a DOWN epimutation, does it continue, does it turn off, does it extend to last gen?
          
          if(is_down==1){
            
            if(vector_in[i]==-1&vector_in[i-1]==-1)
              
            {tempL<-tempL+(gen_names[i]-gen_names[i-1])}
            
            
            if(vector_in[i]==0&vector_in[i-1]==-1){
              
              gene_name <- rownames(vector_in)
              
              gene <- strsplit(gene_name, ":")[[1]][4]
              
              number_transitions <- number_transitions +1  # there is a transition to OFF from a DOWN epimutation at generation[i]
              
              length <- tempL
              
              Lineage <- names(tiny_all_epi[e])
              
              complete <- 1
              
              save_mut <- data.frame(gene_name, gene, number_transitions, onset_gen, length,  complete, Lineage, is_up, is_down)   
              
              output_frame <- rbind(output_frame, save_mut)
              
              is_up <- 0
              is_down <- 0
              tempL <- 0
              complete <- 0
              
            }
            
            # deal with a down epimutation extending to end of lineage
            
            if(vector_in[i]==-1&i==length(vector_in)){
              
              gene_name <- rownames(vector_in)
              
              gene <- strsplit(gene_name, ":")[[1]][4]
              
              number_transitions <- number_transitions + 0
              
              length <- tempL
              
              Lineage <- names(tiny_all_epi[e])
              
              complete <- 0
              
              save_mut <- data.frame(gene_name, gene, number_transitions, onset_gen, length,  complete, Lineage, is_up, is_down)   
              output_frame <- rbind(output_frame, save_mut)
              
              is_up <- 0
              is_down <- 0
              tempL <- 0
              complete <- 0
              
            }
            
          }    
          
        }} 
    
    Overall_output <- rbind(Overall_output, output_frame)  
    
  }
  
  tinyRNA_All_output_list_with_missing_as_muts[[e]] <- Overall_output
}

tiny_Allepiduration<-c()
names(tinyRNA_All_output_list_with_missing_as_muts) <- names(tiny_all_epi)
tinyRNA_All_output_list_with_missing_as_muts$C1 <- tinyRNA_All_output_list_with_missing_as_muts$C1 %>%
  # Creating an empty column:
  add_column(Condition = "Control")
tiny_Allepiduration<-tinyRNA_All_output_list_with_missing_as_muts$C1
tinyRNA_All_output_list_with_missing_as_muts$C2 <- tinyRNA_All_output_list_with_missing_as_muts$C2 %>%
  # Creating an empty column:
  add_column(Condition = "Control")
tiny_Allepiduration<-rbind(tiny_Allepiduration, tinyRNA_All_output_list_with_missing_as_muts$C2)
tinyRNA_All_output_list_with_missing_as_muts$L1 <- tinyRNA_All_output_list_with_missing_as_muts$L1 %>%
  # Creating an empty column:
  add_column(Condition = "Low dose")
tiny_Allepiduration<-rbind(tiny_Allepiduration, tinyRNA_All_output_list_with_missing_as_muts$L1)
tinyRNA_All_output_list_with_missing_as_muts$L2 <- tinyRNA_All_output_list_with_missing_as_muts$L2 %>%
  # Creating an empty column:
  add_column(Condition = "Low dose")
tiny_Allepiduration<-rbind(tiny_Allepiduration, tinyRNA_All_output_list_with_missing_as_muts$L2)
tinyRNA_All_output_list_with_missing_as_muts$H1 <- tinyRNA_All_output_list_with_missing_as_muts$H1 %>%
  # Creating an empty column:
  add_column(Condition = "High dose")
tiny_Allepiduration<-rbind(tiny_Allepiduration, tinyRNA_All_output_list_with_missing_as_muts$H1)
tinyRNA_All_output_list_with_missing_as_muts$H2 <- tinyRNA_All_output_list_with_missing_as_muts$H2 %>%
  # Creating an empty column:
  add_column(Condition = "High dose")
tiny_Allepiduration<-rbind(tiny_Allepiduration, tinyRNA_All_output_list_with_missing_as_muts$H2)
head(tiny_Allepiduration)
tiny_all_epimut<-subset(tiny_Allepiduration,tiny_Allepiduration$length !=0)
save(tiny_all_epimut,file="tiny_all_epimut.Rdata")

#22G RNAs
tiny_22G_epimut<-tiny_all_epimut %>% filter(grepl('22G', gene_name))
CoxMod<-coxph(Surv(length,complete)~Condition,data=tiny_22G_epimut)
ggforest(CoxMod, data=tiny_22G_epimut)

fit <- survfit(Surv(length, complete) ~ Condition, data = tiny_22G_epimut)

ggsurvplot(survfit(CoxMod,data=tiny_22G_epimut), palette = "#2E9FDF",cof.int=TRUE,
           ggtheme = theme_minimal())
condition_df <- with(tiny_22G_epimut,
                     data.frame(Condition = c("Control","Low dose","High dose")
                     )
)
condition_df

#Fig.4.B
ggsurvplot <- ggsurvplot(fit,  tiny_22G_epimut, censor = T, break.time.by= 1, pval=TRUE,pval.coord = c(0, 0.03),
                         font.main = c(16, "bold"),
                         palette = c("cornflowerblue", "darkgreen","red"),
                         font.x = c(20,"bold"),
                         font.y = c(20,"bold"),
                         font.legend = 20,
                         legend.labs=c("Control","Low dose","High dose"),
                         font.tickslab = 18)+
  # surv.median.line = "hv")+
  ggtitle("Changes survival all conditions")+
  xlab("Time (generations)")
ggsurvplot$plot +theme(plot.title = element_text(hjust = 0.5))

#Sup Fig 5. B
CoxMod<-coxph(Surv(length,complete)~Lineage,data=tiny_22G_epimut)
ggforest(model=CoxMod,data=tiny_22G_epimut,fontsize = 0.8, noDigits = 2)
fit_2 <- survfit(Surv(length, complete) ~ Lineage, data = tiny_22G_epimut)
ggsurvplot <- ggsurvplot(fit_2,  tiny_22G_epimut, censor = T, break.time.by= 1, pval=TRUE,pval.coord = c(0, 0.03),                         font.main = c(16, "bold"),
                         font.x = c(20,"bold"),
                         palette = c("#3399FF","blue","red","2E9FDF", "green","#006600"),
                         font.y = c(20,"bold"),
                         font.legend = 20,
                         legend.labs=c("C1","C2","H1","H2","L1","L2"),
                         font.tickslab = 18)+
  ggtitle("Changes survival all conditions")+
  xlab("Time (generations)")
ggsurvplot$plot +theme(plot.title = element_text(hjust = 0.5)) 

#Get raw data
write.xlsx(RNA_epimut,"Data_Fig_4_B_&_Sup_Fig_5_B.xlsx")
#--------------------------
#####Sup.Fig.7.A - Overlap 22G-RNA epimutation & gene exp.change####
load("tinyRNA_C1epimutations_22G.Rdata")
load("tinyRNA_C2epimutations_22G.Rdata")
load("tinyRNA_L1epimutations_22G.Rdata")
load("tinyRNA_L2epimutations_22G.Rdata")
load("tinyRNA_H1epimutations_22G.Rdata")
load("tinyRNA_H2epimutations_22G.Rdata")
#Data prep.
tinyNames<-tinyRNA_C1epimutations_22G$genes
tinyRNA_C1epimutations_22G <- subset(tinyRNA_C1epimutations_22G,  select = -Feature.ID)
tinyRNA_C1epimutations_22G <- subset(tinyRNA_C1epimutations_22G,  select = -Tag)
tinyRNA_C1epimutations_22G <- subset(tinyRNA_C1epimutations_22G,  select = -Feature.Name)
tinyRNA_C1epimutations_22G <- subset(tinyRNA_C1epimutations_22G,  select = -Feature.Class)
tinyRNA_C1epimutations_22G <- subset(tinyRNA_C1epimutations_22G,  select = -condition)
tinyRNA_C1epimutations_22G_unique<-unique(tinyRNA_C1epimutations_22G)
rownames(tinyRNA_C1epimutations_22G_unique)<-tinyRNA_C1epimutations_22G_unique$genes
tinyRNA_C1epimutations_22G_unique <- subset(tinyRNA_C1epimutations_22G_unique,  select = -genes)

tinyRNA_C2epimutations_22G <- subset(tinyRNA_C2epimutations_22G,  select = -Feature.ID)
tinyRNA_C2epimutations_22G <- subset(tinyRNA_C2epimutations_22G,  select = -Tag)
tinyRNA_C2epimutations_22G <- subset(tinyRNA_C2epimutations_22G,  select = -Feature.Name)
tinyRNA_C2epimutations_22G <- subset(tinyRNA_C2epimutations_22G,  select = -Feature.Class)
tinyRNA_C2epimutations_22G <- subset(tinyRNA_C2epimutations_22G,  select = -condition)
tinyRNA_C2epimutations_22G_unique<-unique(tinyRNA_C2epimutations_22G)
rownames(tinyRNA_C2epimutations_22G_unique)<-tinyRNA_C2epimutations_22G_unique$genes
tinyRNA_C2epimutations_22G_unique <- subset(tinyRNA_C2epimutations_22G_unique,  select = -genes)

tinyRNA_L1epimutations_22G <- subset(tinyRNA_L1epimutations_22G,  select = -Feature.ID)
tinyRNA_L1epimutations_22G <- subset(tinyRNA_L1epimutations_22G,  select = -Tag)
tinyRNA_L1epimutations_22G <- subset(tinyRNA_L1epimutations_22G,  select = -Feature.Name)
tinyRNA_L1epimutations_22G <- subset(tinyRNA_L1epimutations_22G,  select = -Feature.Class)
tinyRNA_L1epimutations_22G <- subset(tinyRNA_L1epimutations_22G,  select = -condition)
tinyRNA_L1epimutations_22G_unique<-unique(tinyRNA_L1epimutations_22G)
rownames(tinyRNA_L1epimutations_22G_unique)<-tinyRNA_L1epimutations_22G_unique$genes
tinyRNA_L1epimutations_22G_unique <- subset(tinyRNA_L1epimutations_22G_unique,  select = -genes)

tinyRNA_L2epimutations_22G <- subset(tinyRNA_L2epimutations_22G,  select = -Feature.ID)
tinyRNA_L2epimutations_22G <- subset(tinyRNA_L2epimutations_22G,  select = -Tag)
tinyRNA_L2epimutations_22G <- subset(tinyRNA_L2epimutations_22G,  select = -Feature.Name)
tinyRNA_L2epimutations_22G <- subset(tinyRNA_L2epimutations_22G,  select = -Feature.Class)
tinyRNA_L2epimutations_22G <- subset(tinyRNA_L2epimutations_22G,  select = -condition)
tinyRNA_L2epimutations_22G_unique<-unique(tinyRNA_L2epimutations_22G)
rownames(tinyRNA_L2epimutations_22G_unique)<-tinyRNA_L2epimutations_22G_unique$genes
tinyRNA_L2epimutations_22G_unique <- subset(tinyRNA_L2epimutations_22G_unique,  select = -genes)

tinyRNA_H1epimutations_22G <- subset(tinyRNA_H1epimutations_22G,  select = -Feature.ID)
tinyRNA_H1epimutations_22G <- subset(tinyRNA_H1epimutations_22G,  select = -Tag)
tinyRNA_H1epimutations_22G <- subset(tinyRNA_H1epimutations_22G,  select = -Feature.Name)
tinyRNA_H1epimutations_22G <- subset(tinyRNA_H1epimutations_22G,  select = -Feature.Class)
tinyRNA_H1epimutations_22G <- subset(tinyRNA_H1epimutations_22G,  select = -condition)
tinyRNA_H1epimutations_22G_unique<-unique(tinyRNA_H1epimutations_22G)
rownames(tinyRNA_H1epimutations_22G_unique)<-tinyRNA_H1epimutations_22G_unique$genes
tinyRNA_H1epimutations_22G_unique <- subset(tinyRNA_H1epimutations_22G_unique,  select = -genes)

tinyRNA_H2epimutations_22G <- subset(tinyRNA_H2epimutations_22G,  select = -Feature.ID)
tinyRNA_H2epimutations_22G <- subset(tinyRNA_H2epimutations_22G,  select = -Tag)
tinyRNA_H2epimutations_22G <- subset(tinyRNA_H2epimutations_22G,  select = -Feature.Name)
tinyRNA_H2epimutations_22G <- subset(tinyRNA_H2epimutations_22G,  select = -Feature.Class)
tinyRNA_H2epimutations_22G <- subset(tinyRNA_H2epimutations_22G,  select = -condition)
tinyRNA_H2epimutations_22G_unique<-unique(tinyRNA_H2epimutations_22G)
rownames(tinyRNA_H2epimutations_22G_unique)<-tinyRNA_H2epimutations_22G_unique$genes
tinyRNA_H2epimutations_22G_unique <- subset(tinyRNA_H2epimutations_22G_unique,  select = -genes)

Gene_Features<-read.xlsx("Gene_Features.xlsx")
piRNA_csr_hrde <- Gene_Features[, c(1, 8:10)]

piRNA_target_genes <- unique(piRNA_csr_hrde[piRNA_csr_hrde$piRNA_targets==TRUE, 1])

csr_genes <- unique(piRNA_csr_hrde[piRNA_csr_hrde$csr1_targets==TRUE,  1])

hrde_genes <- unique(piRNA_csr_hrde[piRNA_csr_hrde$hrde1_targets ==TRUE,  1])

#Creation of an RNA centric table with RNA genes, whether they map to small RNA, RNA expression change and small RNA epimutation status
#Control
Control_RNA_list <- list(RNA_C1epimutations, RNA_C2epimutationsbis)
Control_small_RNA_list <- list(tinyRNA_C1epimutations_22G_unique, tinyRNA_C2epimutations_22G_unique)
Control_smallRNA_genes <- rownames(tinyRNA_C1epimutations_22G_unique)

lin <- c("C1", "C2")
Control_RNA_and_smallRNA_integrated_table <- c()

for(x in 1:length(Control_RNA_list)){
  RNA_bin <- Control_RNA_list[[x]]
  smallRNA_bin <- Control_small_RNA_list[[x]]
  Lineage <- lin[[x]]
  integrated_table <- c()
  ep_table <- c()
  
  for(e in 1:length(rownames(RNA_bin))){
    coord <- rownames(RNA_bin)[e]
    gene <- strsplit(coord, ":")[[1]][4]
    RNA_UP <- 0
    RNA_DOWN <- 0
    RNA_MixedUPDOWN <- 0
    
    RNA_mut <- 0
    
    if(sum(abs(RNA_bin[coord, ])) >0){
      RNA_mut <- 1
      # Direction of RNA expression changes  
      # Determine the original direction of the RNA events    
      events <- sum(abs(RNA_bin[coord, ]))
      sum <- sum(RNA_bin[coord, ])
      
      if(sum == events){
        RNA_UP <- 1
      }     
      
      if(sum == -1*(events)){
        RNA_DOWN <- 1
      }     
      
      if(!abs(sum) == events){
        RNA_MixedUPDOWN <- 1
      }  
    }
    RNA_gens <- "0"
    RNA_inherited <- 0
    
    if(RNA_mut ==1){
      RNA_gens <-colnames(RNA_bin[(which(abs(RNA_bin[coord,])>0))])
      
      if(length(RNA_gens) > 1){
        RNA_gens <- paste(RNA_gens, collapse = "_")
        
        if(2 %in% diff(as.numeric(strsplit(RNA_gens, "_")[[1]]))){
          RNA_inherited <- 1
        } 
      }
    }
    # The RNA coding locus has to a) be present in the small RNA genes so as to get a small RNA signal data point
    # b) be present in the Ahringer data set so as to get associated annotations on chromatin domain distribution
    # First assume the coding locus with RNAseq data is not present in small RNA 
    maps_to_Ahr <- 0
    maps_to_smallRNA <- 0
    smallRNA_mut <- 0
    smallRNA_inherited <- 0
    smallRNA_ep_gens <- 0
    multiple_smallRNA_ep_gens  <- 0
    smallRNA_UP <- 0
    smallRNA_DOWN <- 0
    smallRNA_mixed_dir <- 0
    piRNA_cluster <- 0
    
    if(gene %in% piRNA_cluster_genes){
      piRNA_cluster <- 1
    }
    
    if(gene %in% Ahringer_single_gene_ref_table$Gene){
      maps_to_Ahr <- 1} 
    
    if(gene %in% Control_smallRNA_genes == T){
      maps_to_smallRNA <- 1
      
      if(sum(abs(smallRNA_bin[gene, ]))>0){
        smallRNA_mut <- 1
        target_row <- smallRNA_bin[gene, ]
        events <- sum(abs(smallRNA_bin[gene, ]))
        # Direction
        
        if(!-1 %in% target_row){
          smallRNA_UP <- 1
        }
        
        if(!1 %in% target_row){
          smallRNA_DOWN <- 1
        }
        
        if(abs(sum(target_row)) < events){
          smallRNA_mixed_dir <- 1
        }
      }
      if(smallRNA_mut ==1){
        # Generations
        smallRNA_ep_gens <- colnames(smallRNA_bin[(which(abs(smallRNA_bin[gene, ]) > 0))])
        
        if(length(smallRNA_ep_gens)> 1){
          smallRNA_ep_gens <- paste(smallRNA_ep_gens, collapse = "_")
          multiple_smallRNA_ep_gens <- 1
        }
        # Inherited
        inherited <- c()
        inherit <- 0 
        each_row <- colnames(smallRNA_bin[(which(abs(smallRNA_bin[gene, ]) > 0))])
        
        if(2 %in% diff(as.numeric(each_row))){
          inherit <- 1}
        inherited <- c(inherited, inherit) 
        
        if(sum(inherited > 0)){
          smallRNA_inherited <- 1
        }
      }}
    time_matched_to_RNA <- 0
    smallRNA_target_is_Active <- 0
    smallRNA_target_is_Regulated <- 0
    smallRNA_target_is_Border <- 0
    smallRNA_target_is_ChrX <- 0
    Mixed_Domain_smallRNA_targets <- 0
    # Determine Domain of target gene
    # Domain can only be determined if the gene also maps to Ahringer
    # Similarly, regulatory element chromatin epimutation status can be determined
    
    if(maps_to_Ahr == 1){
      Domain <-  Ahringer_single_gene_ref_table[Ahringer_single_gene_ref_table$Gene %in% gene , 4]
      
      if("A" %in% unique(Domain)){ 
        smallRNA_target_is_Active <- 1
      }  
      
      if("R" %in% unique(Domain)){ 
        smallRNA_target_is_Regulated <- 1
      }  
      if("border" %in% unique(Domain)){ 
        smallRNA_target_is_Border <- 1
      }
      if("." %in% unique(Domain)){ 
        smallRNA_target_is_ChrX <- 1
      }  
      if(length(unique(Domain)) > 1){
        Mixed_Domain_smallRNA_targets <- 1
      }
    }
    # Determine Direction matching
    Direction_match_to_RNA <- 0
    
    if(smallRNA_UP==1&RNA_UP==1){
      Direction_match_to_RNA <- 1
    }
    
    if(smallRNA_DOWN==1&RNA_DOWN==1){
      Direction_match_to_RNA <- 1
    } 
    # Determine time matching of small RNA to RNA
    
    if(RNA_mut==1&smallRNA_mut==1){
      smallRNA_gens <- unlist(str_split(smallRNA_ep_gens, "_"))
      rna_gens <- unlist(str_split(RNA_gens, "_"))
      matched_gens <- smallRNA_gens[smallRNA_gens %in% rna_gens]
      matched <-  length(which(smallRNA_gens %in% rna_gens))
      
      if(matched > 0){
        time_matched_to_RNA <- 1
      }
    }
    save <- data.frame(Lineage, gene, coord,RNA_gens, maps_to_Ahr, maps_to_smallRNA, RNA_mut, RNA_inherited, 
                       smallRNA_mut, smallRNA_target_is_Active, smallRNA_target_is_Regulated, 
                       smallRNA_target_is_Border, smallRNA_target_is_ChrX, Mixed_Domain_smallRNA_targets, 
                       smallRNA_inherited, time_matched_to_RNA,smallRNA_ep_gens,
                       smallRNA_UP, smallRNA_DOWN, smallRNA_mixed_dir, 
                       RNA_UP, RNA_DOWN, RNA_MixedUPDOWN, Direction_match_to_RNA, 
                       piRNA_cluster)
    integrated_table <- rbind(integrated_table, save)
  }
  Control_RNA_and_smallRNA_integrated_table <- rbind(Control_RNA_and_smallRNA_integrated_table, integrated_table)
}
save(Control_RNA_and_smallRNA_integrated_table,file="Control_RNA_and_smallRNA_integrated_table.Rdata")
#Low dose
Low_dose_RNA_list <- list(RNA_L1epimutations, RNA_L2epimutations)
Low_dose_small_RNA_list <- list(tinyRNA_L1epimutations_22G_unique, tinyRNA_L2epimutations_22G_unique)
Low_dose_smallRNA_genes <- rownames(tinyRNA_L1epimutations_22G_unique)

lin <- c("L1", "L2")
Low_dose_RNA_and_smallRNA_integrated_table <- c()

for(x in 1:length(Low_dose_RNA_list)){
  RNA_bin <- Low_dose_RNA_list[[x]]
  smallRNA_bin <- Low_dose_small_RNA_list[[x]]
  Lineage <- lin[[x]]
  integrated_table <- c()
  ep_table <- c()
  
  for(e in 1:length(rownames(RNA_bin))){
    coord <- rownames(RNA_bin)[e]
    gene <- strsplit(coord, ":")[[1]][4]
    RNA_UP <- 0
    RNA_DOWN <- 0
    RNA_MixedUPDOWN <- 0
    RNA_mut <- 0
    
    if(sum(abs(RNA_bin[coord, ])) >0){
      RNA_mut <- 1
      # Direction of RNA expression changes  
      # Determine the original direction of the RNA events    
      events <- sum(abs(RNA_bin[coord, ]))
      sum <- sum(RNA_bin[coord, ])
      
      if(sum == events){
        RNA_UP <- 1
      }     
      if(sum == -1*(events)){
        RNA_DOWN <- 1
      }     
      if(!abs(sum) == events){
        RNA_MixedUPDOWN <- 1
      }  
    }
    RNA_gens <- "0"
    RNA_inherited <- 0
    
    if(RNA_mut ==1){
      RNA_gens <-colnames(RNA_bin[(which(abs(RNA_bin[coord,])>0))])
      
      if(length(RNA_gens) > 1){
        RNA_gens <- paste(RNA_gens, collapse = "_")
        
        if(2 %in% diff(as.numeric(strsplit(RNA_gens, "_")[[1]]))){
          RNA_inherited <- 1
        } 
      }
    }
    # The RNA coding locus has to a) be present in the small RNA genes so as to get a small RNA signal data point
    # b) be present in the Ahringer data set so as to get associated annotations on chromatin domain distribution
    # First assume the coding locus with RNAseq data is not present in small RNA 
    maps_to_Ahr <- 0
    maps_to_smallRNA <- 0
    smallRNA_mut <- 0
    smallRNA_inherited <- 0
    smallRNA_ep_gens <- 0
    multiple_smallRNA_ep_gens  <- 0
    smallRNA_UP <- 0
    smallRNA_DOWN <- 0
    smallRNA_mixed_dir <- 0
    piRNA_cluster <- 0
    
    if(gene %in% piRNA_cluster_genes){
      piRNA_cluster <- 1
    }
    
    if(gene %in% Ahringer_single_gene_ref_table$Gene){
      maps_to_Ahr <- 1} 
    
    if(gene %in% Low_dose_smallRNA_genes == T){
      maps_to_smallRNA <- 1
      
      if(sum(abs(smallRNA_bin[gene, ]))>0){
        smallRNA_mut <- 1
        target_row <- smallRNA_bin[gene, ]
        events <- sum(abs(smallRNA_bin[gene, ]))
        # Direction
        
        if(!-1 %in% target_row){
          smallRNA_UP <- 1
        }
        if(!1 %in% target_row){
          smallRNA_DOWN <- 1
        }
        if(abs(sum(target_row)) < events){
          smallRNA_mixed_dir <- 1
        }
      }
      if(smallRNA_mut ==1){
        # Generations
        smallRNA_ep_gens <- colnames(smallRNA_bin[(which(abs(smallRNA_bin[gene, ]) > 0))]) 
        
        if(length(smallRNA_ep_gens)> 1){
          smallRNA_ep_gens <- paste(smallRNA_ep_gens, collapse = "_")
          multiple_smallRNA_ep_gens <- 1
        }
        # Inherited
        inherited <- c()
        inherit <- 0 
        each_row <- colnames(smallRNA_bin[(which(abs(smallRNA_bin[gene, ]) > 0))])
        
        if(2 %in% diff(as.numeric(each_row))){
          inherit <- 1}
        inherited <- c(inherited, inherit) 
        
        if(sum(inherited > 0)){
          smallRNA_inherited <- 1
        }
      }}
    time_matched_to_RNA <- 0
    smallRNA_target_is_Active <- 0
    smallRNA_target_is_Regulated <- 0
    smallRNA_target_is_Border <- 0
    smallRNA_target_is_ChrX <- 0
    Mixed_Domain_smallRNA_targets <- 0
    # Determine Domain of target gene
    # Domain can only be determined if the gene also maps to Ahringer
    # Similarly, regulatory element chromatin epimutation status can be determined
    
    if(maps_to_Ahr == 1){
      Domain <-  Ahringer_single_gene_ref_table[Ahringer_single_gene_ref_table$Gene %in% gene , 4]
      
      if("A" %in% unique(Domain)){ 
        smallRNA_target_is_Active <- 1
      }  
      if("R" %in% unique(Domain)){ 
        smallRNA_target_is_Regulated <- 1
      }  
      if("border" %in% unique(Domain)){ 
        smallRNA_target_is_Border <- 1
      }
      if("." %in% unique(Domain)){ 
        smallRNA_target_is_ChrX <- 1
      }  
      if(length(unique(Domain)) > 1){
        Mixed_Domain_smallRNA_targets <- 1
      }
    }
    # Determine Direction matching
    Direction_match_to_RNA <- 0
    
    if(smallRNA_UP==1&RNA_UP==1){
      Direction_match_to_RNA <- 1
    }
    
    if(smallRNA_DOWN==1&RNA_DOWN==1){
      Direction_match_to_RNA <- 1
    } 
    # Determine time matching of small RNA to RNA
    
    if(RNA_mut==1&smallRNA_mut==1){
      smallRNA_gens <- unlist(str_split(smallRNA_ep_gens, "_"))
      rna_gens <- unlist(str_split(RNA_gens, "_"))
      matched_gens <- smallRNA_gens[smallRNA_gens %in% rna_gens]
      matched <-  length(which(smallRNA_gens %in% rna_gens))
      
      if(matched > 0){
        time_matched_to_RNA <- 1
      }
    }
    save <- data.frame(Lineage, gene, coord,RNA_gens, maps_to_Ahr, maps_to_smallRNA, RNA_mut, RNA_inherited, 
                       smallRNA_mut, smallRNA_target_is_Active, smallRNA_target_is_Regulated, 
                       smallRNA_target_is_Border, smallRNA_target_is_ChrX, Mixed_Domain_smallRNA_targets, 
                       smallRNA_inherited, time_matched_to_RNA,
                       smallRNA_UP, smallRNA_DOWN, smallRNA_mixed_dir, 
                       RNA_UP, RNA_DOWN, RNA_MixedUPDOWN, Direction_match_to_RNA, 
                       piRNA_cluster)
    integrated_table <- rbind(integrated_table, save)
  }
  Low_dose_RNA_and_smallRNA_integrated_table <- rbind(Low_dose_RNA_and_smallRNA_integrated_table, integrated_table)
}
save(Low_dose_RNA_and_smallRNA_integrated_table,file="Low_dose_RNA_and_smallRNA_integrated_table.Rdata")
#High dose
High_dose_RNA_list <- list(RNA_H1epimutationsbis, RNA_H2epimutations)
High_dose_small_RNA_list <- list(tinyRNA_H1epimutations_22G_unique, tinyRNA_H2epimutations_22G_unique)
High_dose_smallRNA_genes <- rownames(tinyRNA_H1epimutations_22G_unique)

lin <- c("H1", "H2")
High_dose_RNA_and_smallRNA_integrated_table <- c()

for(x in 1:length(High_dose_RNA_list)){
  RNA_bin <- High_dose_RNA_list[[x]]
  smallRNA_bin <- High_dose_small_RNA_list[[x]]
  Lineage <- lin[[x]]
  integrated_table <- c()
  ep_table <- c()
  
  for(e in 1:length(rownames(RNA_bin))){
    coord <- rownames(RNA_bin)[e]
    gene <- strsplit(coord, ":")[[1]][4]
    RNA_UP <- 0
    RNA_DOWN <- 0
    RNA_MixedUPDOWN <- 0
    RNA_mut <- 0
    
    if(sum(abs(RNA_bin[coord, ])) >0){
      RNA_mut <- 1
      # Direction of RNA expression changes  
      # Determine the original direction of the RNA events    
      events <- sum(abs(RNA_bin[coord, ]))
      sum <- sum(RNA_bin[coord, ])
      
      if(sum == events){
        RNA_UP <- 1
      }     
      if(sum == -1*(events)){
        RNA_DOWN <- 1
      }     
      if(!abs(sum) == events){
        RNA_MixedUPDOWN <- 1
      }  
    }
    RNA_gens <- "0"
    RNA_inherited <- 0
    
    if(RNA_mut ==1){
      RNA_gens <- colnames(RNA_bin[(which(abs(RNA_bin[coord, ]) > 0))])
      
      if(length(RNA_gens) > 1){
        RNA_gens <- paste(RNA_gens, collapse = "_")
        
        if(2 %in% diff(as.numeric(strsplit(RNA_gens, "_")[[1]]))){
          RNA_inherited <- 1
        } 
      }
    }
    # The RNA coding locus has to a) be present in the small RNA genes so as to get a small RNA signal data point
    # b) be present in the Ahringer data set so as to get associated annotations on chromatin domain distribution
    # First assume the coding locus with RNAseq data is not present in small RNA 
    maps_to_Ahr <- 0
    maps_to_smallRNA <- 0
    smallRNA_mut <- 0
    smallRNA_inherited <- 0
    smallRNA_ep_gens <- 0
    multiple_smallRNA_ep_gens  <- 0
    smallRNA_UP <- 0
    smallRNA_DOWN <- 0
    smallRNA_mixed_dir <- 0
    piRNA_cluster <- 0
    
    if(gene %in% piRNA_cluster_genes){
      piRNA_cluster <- 1
    }
    if(gene %in% Ahringer_single_gene_ref_table$Gene){
      maps_to_Ahr <- 1} 
    
    if(gene %in% High_dose_smallRNA_genes == T){
      maps_to_smallRNA <- 1
      
      if(sum(abs(smallRNA_bin[gene, ]))>0){
        smallRNA_mut <- 1
        target_row <- smallRNA_bin[gene, ]
        events <- sum(abs(smallRNA_bin[gene, ]))
        # Direction
        
        if(!-1 %in% target_row){
          smallRNA_UP <- 1
        }
        if(!1 %in% target_row){
          smallRNA_DOWN <- 1
        }
        if(abs(sum(target_row)) < events){
          smallRNA_mixed_dir <- 1
        }
      }
      if(smallRNA_mut ==1){
        # Generations
        smallRNA_ep_gens <- colnames(smallRNA_bin[(which(abs(smallRNA_bin[gene, ]) > 0))]) 
        
        if(length(smallRNA_ep_gens)> 1){
          smallRNA_ep_gens <- paste(smallRNA_ep_gens, collapse = "_")
          multiple_smallRNA_ep_gens <- 1
        }
        # Inherited
        inherited <- c()
        inherit <- 0 
        each_row <- colnames(smallRNA_bin[(which(abs(smallRNA_bin[gene, ]) > 0))])
        
        if(2 %in% diff(as.numeric(each_row))){
          inherit <- 1}
        inherited <- c(inherited, inherit) 
        
        if(sum(inherited > 0)){
          smallRNA_inherited <- 1
        }
      }}
    time_matched_to_RNA <- 0
    smallRNA_target_is_Active <- 0
    smallRNA_target_is_Regulated <- 0
    smallRNA_target_is_Border <- 0
    smallRNA_target_is_ChrX <- 0
    Mixed_Domain_smallRNA_targets <- 0
    # Determine Domain of target gene
    # Domain can only be determined if the gene also maps to Ahringer
    # Similarly, regulatory element chromatin epimutation status can be determined
    
    if(maps_to_Ahr == 1){
      Domain <-  Ahringer_single_gene_ref_table[Ahringer_single_gene_ref_table$Gene %in% gene , 4]
      
      if("A" %in% unique(Domain)){ 
        smallRNA_target_is_Active <- 1
      }  
      if("R" %in% unique(Domain)){ 
        smallRNA_target_is_Regulated <- 1
      }  
      if("border" %in% unique(Domain)){ 
        smallRNA_target_is_Border <- 1
      }
      if("." %in% unique(Domain)){ 
        smallRNA_target_is_ChrX <- 1
      }  
      if(length(unique(Domain)) > 1){
        Mixed_Domain_smallRNA_targets <- 1
      }
    }
    # Determine Direction matching
    Direction_match_to_RNA <- 0
    
    if(smallRNA_UP==1&RNA_UP==1){
      Direction_match_to_RNA <- 1
    }
    
    if(smallRNA_DOWN==1&RNA_DOWN==1){
      Direction_match_to_RNA <- 1
    } 
    # Determine time matching of small RNA to RNA
    
    if(RNA_mut==1&smallRNA_mut==1){
      smallRNA_gens <- unlist(str_split(smallRNA_ep_gens, "_"))
      rna_gens <- unlist(str_split(RNA_gens, "_"))
      matched_gens <- smallRNA_gens[smallRNA_gens %in% rna_gens]
      matched <-  length(which(smallRNA_gens %in% rna_gens))
      
      if(matched > 0){
        time_matched_to_RNA <- 1
      }
    }
    save <- data.frame(Lineage, gene, coord, maps_to_Ahr, maps_to_smallRNA, RNA_mut, RNA_inherited, 
                       smallRNA_mut, smallRNA_target_is_Active, smallRNA_target_is_Regulated, 
                       smallRNA_target_is_Border, smallRNA_target_is_ChrX, Mixed_Domain_smallRNA_targets, 
                       smallRNA_inherited, time_matched_to_RNA,
                       smallRNA_UP, smallRNA_DOWN, smallRNA_mixed_dir, 
                       RNA_UP, RNA_DOWN, RNA_MixedUPDOWN, Direction_match_to_RNA, 
                       piRNA_cluster)
    integrated_table <- rbind(integrated_table, save)
  }
  High_dose_RNA_and_smallRNA_integrated_table <- rbind(High_dose_RNA_and_smallRNA_integrated_table, integrated_table)
}
save(High_dose_RNA_and_smallRNA_integrated_table,file="High_dose_RNA_and_smallRNA_integrated_table.Rdata")

#Overlaps calculation Control
#Inherited RNA seq changes
#Nb genes with RNA inherited and 22G epimut
Control_genes_RNA_muts_inherited_and_smallRNA_muts <- nrow(Control_RNA_and_smallRNA_integrated_table[Control_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Control_RNA_and_smallRNA_integrated_table$RNA_mut==1&Control_RNA_and_smallRNA_integrated_table$RNA_inherited==1&Control_RNA_and_smallRNA_integrated_table$smallRNA_mut==1, ])
#Nb genes with both RNA and 22G inherited epimut at matching time point
Control_RNA_inherit_smallRNA_mut_tm_inherit <- nrow(Control_RNA_and_smallRNA_integrated_table[Control_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1& Control_RNA_and_smallRNA_integrated_table$RNA_mut==1&Control_RNA_and_smallRNA_integrated_table$RNA_inherited==1&Control_RNA_and_smallRNA_integrated_table$smallRNA_mut==1&Control_RNA_and_smallRNA_integrated_table$smallRNA_inherited==1&Control_RNA_and_smallRNA_integrated_table$time_matched_to_RNA==1, ])/Control_genes_RNA_muts_inherited_and_smallRNA_muts*100
# genes with RNA inherit and 22G non inherited epimut at matching time point
Control_RNA_inherit_smallRNA_mut_tm_non_inherit <- nrow(Control_RNA_and_smallRNA_integrated_table[Control_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Control_RNA_and_smallRNA_integrated_table$RNA_mut==1&Control_RNA_and_smallRNA_integrated_table$RNA_inherited==1&Control_RNA_and_smallRNA_integrated_table$smallRNA_mut==1&Control_RNA_and_smallRNA_integrated_table$smallRNA_inherited==0&Control_RNA_and_smallRNA_integrated_table$time_matched_to_RNA==1, ])/Control_genes_RNA_muts_inherited_and_smallRNA_muts*100
#% genes with RNA inherit and 22G non inherited epimut at non matching time point
Control_RNA_inherit_smallRNA_mut_non_tm <-  nrow(Control_RNA_and_smallRNA_integrated_table[Control_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Control_RNA_and_smallRNA_integrated_table$RNA_mut==1&Control_RNA_and_smallRNA_integrated_table$RNA_inherited==1&Control_RNA_and_smallRNA_integrated_table$smallRNA_mut==1&Control_RNA_and_smallRNA_integrated_table$time_matched_to_RNA==0, ])/Control_genes_RNA_muts_inherited_and_smallRNA_muts*100

Control_RNA_Inherit_proportions <-  data.frame(c(Control_RNA_inherit_smallRNA_mut_non_tm, 
                                                 Control_RNA_inherit_smallRNA_mut_tm_non_inherit,
                                                 Control_RNA_inherit_smallRNA_mut_tm_inherit))
Names <- c(
  "Has small RNA epimutations - not time matched", 
  "Has time matched small RNA epimutations - not inherited", 
  "Has time matched small RNA epimutations - inherited")                                           
Species <- rep("Genes with inherited \nRNA expression changes\ntargeted by 22G epimutations", 3)

Control_RNA_Inherit <- cbind(Control_RNA_Inherit_proportions, Names, Species)
colnames(Control_RNA_Inherit) <- c("Percentages", "Epigenetic_feature_of_gene", "Species")

#Non Inherited RNA seq changes

Control_genes_RNA_muts_n_inherited_smallRNA_muts <- nrow(Control_RNA_and_smallRNA_integrated_table[Control_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Control_RNA_and_smallRNA_integrated_table$RNA_mut==1&Control_RNA_and_smallRNA_integrated_table$RNA_inherited==0&Control_RNA_and_smallRNA_integrated_table$smallRNA_mut==1, ])
Control_RNA_n_inherit_smallRNA_mut_tm_inherit <- nrow(Control_RNA_and_smallRNA_integrated_table[Control_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Control_RNA_and_smallRNA_integrated_table$RNA_mut==1&Control_RNA_and_smallRNA_integrated_table$RNA_inherited==0&Control_RNA_and_smallRNA_integrated_table$smallRNA_mut==1&Control_RNA_and_smallRNA_integrated_table$smallRNA_inherited==1&Control_RNA_and_smallRNA_integrated_table$time_matched_to_RNA==1, ])/Control_genes_RNA_muts_n_inherited_smallRNA_muts*100

Control_RNA_n_inherit_smallRNA_mut_tm_non_inherit <- nrow(Control_RNA_and_smallRNA_integrated_table[Control_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Control_RNA_and_smallRNA_integrated_table$RNA_mut==1&Control_RNA_and_smallRNA_integrated_table$RNA_inherited==0&Control_RNA_and_smallRNA_integrated_table$smallRNA_mut==1&Control_RNA_and_smallRNA_integrated_table$smallRNA_inherited==0&Control_RNA_and_smallRNA_integrated_table$time_matched_to_RNA==1, ])/Control_genes_RNA_muts_n_inherited_smallRNA_muts*100
Control_RNA_n_inherit_smallRNA_mut_non_tm <-  nrow(Control_RNA_and_smallRNA_integrated_table[Control_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Control_RNA_and_smallRNA_integrated_table$RNA_mut==1&Control_RNA_and_smallRNA_integrated_table$RNA_inherited==0&Control_RNA_and_smallRNA_integrated_table$smallRNA_mut==1&Control_RNA_and_smallRNA_integrated_table$time_matched_to_RNA==0, ])/Control_genes_RNA_muts_n_inherited_smallRNA_muts*100

Control_RNA_n_Inherit_proportions <-  data.frame(c(Control_RNA_n_inherit_smallRNA_mut_non_tm, 
                                                   Control_RNA_n_inherit_smallRNA_mut_tm_non_inherit,
                                                   Control_RNA_n_inherit_smallRNA_mut_tm_inherit))


Names <- c(
  "Has small RNA epimutations - not time matched", 
  "Has time matched small RNA epimutations - not inherited", 
  "Has time matched small RNA epimutations - inherited")                                             

Species <- rep("Genes with non inherited \nRNA expression changes\ntargeted by 22G epimutations", 3)

Control_RNA_N_Inherit <- cbind(Control_RNA_n_Inherit_proportions, Names, Species)

colnames(Control_RNA_N_Inherit) <- c("Percentages", "Epigenetic_feature_of_gene", "Species")

#Plot for RNA seq changes

Control_RNA_Inherited_features <- rbind(Control_RNA_N_Inherit, Control_RNA_Inherit)
Control_RNA_Inherited_features <- Control_RNA_Inherited_features[order(nrow(Control_RNA_Inherited_features):1), ] 

Control_RNA_Inherited_features$Epigenetic_feature_of_gene <- factor(Control_RNA_Inherited_features$Epigenetic_feature_of_gene, levels = unique(Control_RNA_Inherited_features$Epigenetic_feature_of_gene))
positions <- unique(Control_RNA_Inherited_features$Species)

# Proportions of genes with inherited and non inherited RNAseq changes that have small RNA changes

Control_RNA_w_smallRNA_proportion_plot <-
  
  ggplot(Control_RNA_Inherited_features, aes(fill=Epigenetic_feature_of_gene, y=Percentages, x=Species)) + 
  
  geom_bar(position="stack", stat="identity", width = 0.4, color="black")+
  scale_x_discrete(limits = positions)+
  scale_fill_manual(values=c("cornflowerblue", 
                             "cornsilk3", 
                             "coral",
                             "cornflowerblue", 
                             "cornsilk3", 
                             "coral"))+
  theme_bw()+
  labs(x= " ", 
       y = "Percentage of genes (%)\n")+
  theme(plot.margin = unit(c(2,1,2,2), "cm"))+
  theme(plot.title = element_text(size = 16))+
  labs(title = "Control data")+
  theme(text = element_text(size=14))+
  theme(axis.text.x=element_text(colour="black", size = 12))+
  guides(fill=guide_legend(ncol=1))+
  scale_fill_manual(name = "RNA expression change status",
                    values=c("cornflowerblue", 
                             "cornsilk3", 
                             "coral"),
                    labels=c("Simultaneous inherited 22G epimutation", 
                             "Simultaneous non inherited 22G epimutation", 
                             "Non simultaneous 22G epimutation"))

Control_RNA_w_smallRNA_proportion_plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(fill=guide_legend(title="Associated 22G epimutation status"))+
  theme(legend.text = element_text(size=14))+
  theme(legend.title = element_text(size=14))

#Get raw data 
write.xlsx(Control_RNA_Inherited_features,"Data_Sup_Fig_6.xlsx")

#Overlaps Low dose
#Inherited RNA seq changes
Low_dose_genes_RNA_muts_inherited_smallRNA_muts <- nrow(Low_dose_RNA_and_smallRNA_integrated_table[Low_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Low_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1&Low_dose_RNA_and_smallRNA_integrated_table$RNA_inherited==1&Low_dose_RNA_and_smallRNA_integrated_table$smallRNA_mut==1, ])
Low_dose_RNA_inherit_smallRNA_mut_tm_inherit <- nrow(Low_dose_RNA_and_smallRNA_integrated_table[Low_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1& Low_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1&Low_dose_RNA_and_smallRNA_integrated_table$RNA_inherited==1&Low_dose_RNA_and_smallRNA_integrated_table$smallRNA_mut==1&Low_dose_RNA_and_smallRNA_integrated_table$smallRNA_inherited==1&Low_dose_RNA_and_smallRNA_integrated_table$time_matched_to_RNA==1, ])/Low_dose_genes_RNA_muts_inherited_smallRNA_muts*100

Low_dose_RNA_inherit_smallRNA_mut_tm_non_inherit <- nrow(Low_dose_RNA_and_smallRNA_integrated_table[Low_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Low_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1&Low_dose_RNA_and_smallRNA_integrated_table$RNA_inherited==1&Low_dose_RNA_and_smallRNA_integrated_table$smallRNA_mut==1&Low_dose_RNA_and_smallRNA_integrated_table$smallRNA_inherited==0&Low_dose_RNA_and_smallRNA_integrated_table$time_matched_to_RNA==1, ])/Low_dose_genes_RNA_muts_inherited_smallRNA_muts*100
Low_dose_RNA_inherit_smallRNA_mut_non_tm <-  nrow(Low_dose_RNA_and_smallRNA_integrated_table[Low_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Low_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1&Low_dose_RNA_and_smallRNA_integrated_table$RNA_inherited==1&Low_dose_RNA_and_smallRNA_integrated_table$smallRNA_mut==1&Low_dose_RNA_and_smallRNA_integrated_table$time_matched_to_RNA==0, ])/Low_dose_genes_RNA_muts_inherited_smallRNA_muts*100

Low_dose_RNA_Inherit_proportions <-  data.frame(c(Low_dose_RNA_inherit_smallRNA_mut_non_tm, 
                                                  Low_dose_RNA_inherit_smallRNA_mut_tm_non_inherit,
                                                  Low_dose_RNA_inherit_smallRNA_mut_tm_inherit))
Names <- c(
  "Has small RNA epimutations - not time matched", 
  "Has time matched small RNA epimutations - not inherited", 
  "Has time matched small RNA epimutations - inherited")                                           
Species <- rep("Genes with inherited \nRNA expression changes\ntargeted by 22G RNA epimutations", 3)

Low_dose_RNA_Inherit <- cbind(Low_dose_RNA_Inherit_proportions, Names, Species)
colnames(Low_dose_RNA_Inherit) <- c("Percentages", "Epigenetic_feature_of_gene", "Species")

#Non Inherited RNA seq changes

Low_dose_genes_RNA_muts_n_inherited_smallRNA_muts <- nrow(Low_dose_RNA_and_smallRNA_integrated_table[Low_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Low_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1&Low_dose_RNA_and_smallRNA_integrated_table$RNA_inherited==0&Low_dose_RNA_and_smallRNA_integrated_table$smallRNA_mut==1, ])
Low_dose_RNA_n_inherit_smallRNA_mut_tm_inherit <- nrow(Low_dose_RNA_and_smallRNA_integrated_table[Low_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Low_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1&Low_dose_RNA_and_smallRNA_integrated_table$RNA_inherited==0&Low_dose_RNA_and_smallRNA_integrated_table$smallRNA_mut==1&Low_dose_RNA_and_smallRNA_integrated_table$smallRNA_inherited==1&Low_dose_RNA_and_smallRNA_integrated_table$time_matched_to_RNA==1, ])/Low_dose_genes_RNA_muts_n_inherited_smallRNA_muts*100

Low_dose_RNA_n_inherit_smallRNA_mut_tm_non_inherit <- nrow(Low_dose_RNA_and_smallRNA_integrated_table[Low_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Low_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1&Low_dose_RNA_and_smallRNA_integrated_table$RNA_inherited==0&Low_dose_RNA_and_smallRNA_integrated_table$smallRNA_mut==1&Low_dose_RNA_and_smallRNA_integrated_table$smallRNA_inherited==0&Low_dose_RNA_and_smallRNA_integrated_table$time_matched_to_RNA==1, ])/Low_dose_genes_RNA_muts_n_inherited_smallRNA_muts*100
Low_dose_RNA_n_inherit_smallRNA_mut_non_tm <-  nrow(Low_dose_RNA_and_smallRNA_integrated_table[Low_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Low_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1&Low_dose_RNA_and_smallRNA_integrated_table$RNA_inherited==0&Low_dose_RNA_and_smallRNA_integrated_table$smallRNA_mut==1&Low_dose_RNA_and_smallRNA_integrated_table$time_matched_to_RNA==0, ])/Low_dose_genes_RNA_muts_n_inherited_smallRNA_muts*100

Low_dose_RNA_n_Inherit_proportions <-  data.frame(c(Low_dose_RNA_n_inherit_smallRNA_mut_non_tm, 
                                                    Low_dose_RNA_n_inherit_smallRNA_mut_tm_non_inherit,
                                                    Low_dose_RNA_n_inherit_smallRNA_mut_tm_inherit))

Names <- c(
  "Has small RNA epimutations - not time matched", 
  "Has time matched small RNA epimutations - not inherited", 
  "Has time matched small RNA epimutations - inherited")                                             

Species <- rep("Genes with non inherited \nRNA expression changes\ntargeted by 22G RNA epimutations", 3)

Low_dose_RNA_N_Inherit <- cbind(Low_dose_RNA_n_Inherit_proportions, Names, Species)

colnames(Low_dose_RNA_N_Inherit) <- c("Percentages", "Epigenetic_feature_of_gene", "Species")

#Plot for RNA seq changes

Low_dose_RNA_Inherited_features <- rbind(Low_dose_RNA_N_Inherit, Low_dose_RNA_Inherit)
Low_dose_RNA_Inherited_features <- Low_dose_RNA_Inherited_features[order(nrow(Low_dose_RNA_Inherited_features):1), ] 

Low_dose_RNA_Inherited_features$Epigenetic_feature_of_gene <- factor(Low_dose_RNA_Inherited_features$Epigenetic_feature_of_gene, levels = unique(Low_dose_RNA_Inherited_features$Epigenetic_feature_of_gene))
positions <- unique(Low_dose_RNA_Inherited_features$Species)

# Proportions of genes with inherited and non inherited RNAseq changes that have small RNA changes

Low_dose_RNA_w_smallRNA_proportion_plot <-
  ggplot(Low_dose_RNA_Inherited_features, aes(fill=Epigenetic_feature_of_gene, y=Percentages, x=Species)) + 
  geom_bar(position="stack", stat="identity", width = 0.4, color="black")+
  scale_x_discrete(limits = positions)+
  scale_fill_manual(values=c("cornflowerblue", 
                             "cornsilk3", 
                             "coral",
                             "cornflowerblue", 
                             "cornsilk3", 
                             "coral"))+
  theme_bw()+
  labs(x= " ", 
       y = "Percentage of genes (%)\n")+
  theme(plot.margin = unit(c(2,1,2,2), "cm"))+
  theme(plot.title = element_text(size = 16))+
  labs(title = "Low dose data")+
  theme(text = element_text(size=14))+
  theme(axis.text.x=element_text(colour="black", size = 12))+
  guides(fill=guide_legend(ncol=1))+
  scale_fill_manual(name = "RNA expression change status",
                    values=c("cornflowerblue", 
                             "cornsilk3", 
                             "coral"),
                    labels=c("Simultaneous inherited 22G RNA epimutation", 
                             "Simultaneous non inherited 22G RNA epimutation", 
                             "Non simultaneous 22G RNA epimutation"))

Low_dose_RNA_w_smallRNA_proportion_plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(fill=guide_legend(title="Associated 22G RNA epimutation status"))+
  theme(legend.text = element_text(size=14))+
  theme(legend.title = element_text(size=14))

#Get raw data
write.xlsx(Low_dose_RNA_Inherited_features,"Data_Sup_Fig_6_LD.xlsx")

#Overlaps High dose
#Inherited RNA seq changes
High_dose_genes_RNA_muts_inherited_smallRNA_muts <- nrow(High_dose_RNA_and_smallRNA_integrated_table[High_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&High_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1&High_dose_RNA_and_smallRNA_integrated_table$RNA_inherited==1&High_dose_RNA_and_smallRNA_integrated_table$smallRNA_mut==1, ])
High_dose_RNA_inherit_smallRNA_mut_tm_inherit <- nrow(High_dose_RNA_and_smallRNA_integrated_table[High_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1& High_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1&High_dose_RNA_and_smallRNA_integrated_table$RNA_inherited==1&High_dose_RNA_and_smallRNA_integrated_table$smallRNA_mut==1&High_dose_RNA_and_smallRNA_integrated_table$smallRNA_inherited==1&High_dose_RNA_and_smallRNA_integrated_table$time_matched_to_RNA==1, ])/High_dose_genes_RNA_muts_inherited_smallRNA_muts*100

High_dose_RNA_inherit_smallRNA_mut_tm_non_inherit <- nrow(High_dose_RNA_and_smallRNA_integrated_table[High_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&High_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1&High_dose_RNA_and_smallRNA_integrated_table$RNA_inherited==1&High_dose_RNA_and_smallRNA_integrated_table$smallRNA_mut==1&High_dose_RNA_and_smallRNA_integrated_table$smallRNA_inherited==0&High_dose_RNA_and_smallRNA_integrated_table$time_matched_to_RNA==1, ])/High_dose_genes_RNA_muts_inherited_smallRNA_muts*100
High_dose_RNA_inherit_smallRNA_mut_non_tm <-  nrow(High_dose_RNA_and_smallRNA_integrated_table[High_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&High_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1&High_dose_RNA_and_smallRNA_integrated_table$RNA_inherited==1&High_dose_RNA_and_smallRNA_integrated_table$smallRNA_mut==1&High_dose_RNA_and_smallRNA_integrated_table$time_matched_to_RNA==0, ])/High_dose_genes_RNA_muts_inherited_smallRNA_muts*100

High_dose_RNA_Inherit_proportions <-  data.frame(c(High_dose_RNA_inherit_smallRNA_mut_non_tm, 
                                                   High_dose_RNA_inherit_smallRNA_mut_tm_non_inherit,
                                                   High_dose_RNA_inherit_smallRNA_mut_tm_inherit))
Names <- c(
  "Has small RNA epimutations - not time matched", 
  "Has time matched small RNA epimutations - not inherited", 
  "Has time matched small RNA epimutations - inherited")                                           
Species <- rep("Genes with inherited \nRNA expression changes\ntargeted by 22G RNA epimutations", 3)

High_dose_RNA_Inherit <- cbind(High_dose_RNA_Inherit_proportions, Names, Species)
colnames(High_dose_RNA_Inherit) <- c("Percentages", "Epigenetic_feature_of_gene", "Species")

#Non Inherited RNA seq changes
High_dose_genes_RNA_muts_n_inherited_smallRNA_muts <- nrow(High_dose_RNA_and_smallRNA_integrated_table[High_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&High_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1&High_dose_RNA_and_smallRNA_integrated_table$RNA_inherited==0&High_dose_RNA_and_smallRNA_integrated_table$smallRNA_mut==1, ])
High_dose_RNA_n_inherit_smallRNA_mut_tm_inherit <- nrow(High_dose_RNA_and_smallRNA_integrated_table[High_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&High_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1&High_dose_RNA_and_smallRNA_integrated_table$RNA_inherited==0&High_dose_RNA_and_smallRNA_integrated_table$smallRNA_mut==1&High_dose_RNA_and_smallRNA_integrated_table$smallRNA_inherited==1&High_dose_RNA_and_smallRNA_integrated_table$time_matched_to_RNA==1, ])/High_dose_genes_RNA_muts_n_inherited_smallRNA_muts*100

High_dose_RNA_n_inherit_smallRNA_mut_tm_non_inherit <- nrow(High_dose_RNA_and_smallRNA_integrated_table[High_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&High_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1&High_dose_RNA_and_smallRNA_integrated_table$RNA_inherited==0&High_dose_RNA_and_smallRNA_integrated_table$smallRNA_mut==1&High_dose_RNA_and_smallRNA_integrated_table$smallRNA_inherited==0&High_dose_RNA_and_smallRNA_integrated_table$time_matched_to_RNA==1, ])/High_dose_genes_RNA_muts_n_inherited_smallRNA_muts*100
High_dose_RNA_n_inherit_smallRNA_mut_non_tm <-  nrow(High_dose_RNA_and_smallRNA_integrated_table[High_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&High_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1&High_dose_RNA_and_smallRNA_integrated_table$RNA_inherited==0&High_dose_RNA_and_smallRNA_integrated_table$smallRNA_mut==1&High_dose_RNA_and_smallRNA_integrated_table$time_matched_to_RNA==0, ])/High_dose_genes_RNA_muts_n_inherited_smallRNA_muts*100

High_dose_RNA_n_Inherit_proportions <-  data.frame(c(High_dose_RNA_n_inherit_smallRNA_mut_non_tm, 
                                                     High_dose_RNA_n_inherit_smallRNA_mut_tm_non_inherit,
                                                     High_dose_RNA_n_inherit_smallRNA_mut_tm_inherit))

Names <- c(
  "Has small RNA epimutations - not time matched", 
  "Has time matched small RNA epimutations - not inherited", 
  "Has time matched small RNA epimutations - inherited")                                             

Species <- rep("Genes with non inherited \nRNA expression changes\ntargeted by 22G RNA epimutations", 3)

High_dose_RNA_N_Inherit <- cbind(High_dose_RNA_n_Inherit_proportions, Names, Species)

colnames(High_dose_RNA_N_Inherit) <- c("Percentages", "Epigenetic_feature_of_gene", "Species")

#Plot for RNA seq changes

High_dose_RNA_Inherited_features <- rbind(High_dose_RNA_N_Inherit, High_dose_RNA_Inherit)
High_dose_RNA_Inherited_features <- High_dose_RNA_Inherited_features[order(nrow(High_dose_RNA_Inherited_features):1), ] 

High_dose_RNA_Inherited_features$Epigenetic_feature_of_gene <- factor(High_dose_RNA_Inherited_features$Epigenetic_feature_of_gene, levels = unique(High_dose_RNA_Inherited_features$Epigenetic_feature_of_gene))
positions <- unique(High_dose_RNA_Inherited_features$Species)

# Proportions of genes with inherited and non inherited RNAseq changes that have small RNA changes

High_dose_RNA_w_smallRNA_proportion_plot <-
  ggplot(High_dose_RNA_Inherited_features, aes(fill=Epigenetic_feature_of_gene, y=Percentages, x=Species)) + 
  geom_bar(position="stack", stat="identity", width = 0.4, color="black")+
  scale_x_discrete(limits = positions)+
  scale_fill_manual(values=c("cornflowerblue", 
                             "cornsilk3", 
                             "coral",
                             "cornflowerblue", 
                             "cornsilk3", 
                             "coral"))+
  theme_bw()+
  labs(x= " ", 
       y = "Percentage of genes (%)\n")+
  theme(plot.margin = unit(c(2,1,2,2), "cm"))+
  theme(plot.title = element_text(size = 16))+
  labs(title = "High_dose data")+
  theme(text = element_text(size=14))+
  theme(axis.text.x=element_text(colour="black", size = 12))+
  guides(fill=guide_legend(ncol=1))+
  scale_fill_manual(name = "RNA expression change status",
                    values=c("cornflowerblue", 
                             "cornsilk3", 
                             "coral"),
                    labels=c("Simultaneous inherited 22G RNA epimutation", 
                             "Simultaneous non inherited 22G RNA epimutation", 
                             "Non simultaneous 22G RNA epimutation"))

High_dose_RNA_w_smallRNA_proportion_plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(fill=guide_legend(title="Associated 22G RNA epimutation status"))+
  theme(legend.text = element_text(size=14))+
  theme(legend.title = element_text(size=14))

#Get raw data
write.xlsx(High_dose_RNA_Inherited_features,"Data_Sup_Fig_6_HD.xlsx")
#--------------------------
####Sup.Table.6#####
#Odds calculation Control

# 1. What are the odds that genes have both 22G RNA epimutations and RNA changes ?

# background is all genes which have 22G RNA signal, so have to map to small RNA
background <- nrow(Control_RNA_and_smallRNA_integrated_table[Control_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1, ])

# of which have RNA changes
with_RNA_changes <- nrow(Control_RNA_and_smallRNA_integrated_table[Control_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Control_RNA_and_smallRNA_integrated_table$RNA_mut==1, ])

# of which have small RNA changes 
with_smallRNA_changes <- nrow(Control_RNA_and_smallRNA_integrated_table[Control_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Control_RNA_and_smallRNA_integrated_table$smallRNA_mut==1, ])

# of which have both 
with_both <- nrow(Control_RNA_and_smallRNA_integrated_table[Control_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Control_RNA_and_smallRNA_integrated_table$smallRNA_mut==1&Control_RNA_and_smallRNA_integrated_table$RNA_mut==1, ])
smallRNA_RNA <-  with_both
RNA_not_smallRNA <- with_RNA_changes - smallRNA_RNA
smallRNA_not_RNA <- with_smallRNA_changes  -  smallRNA_RNA
not_smallRNA_not_RNA <- background - (with_RNA_changes + with_smallRNA_changes -  smallRNA_RNA)

contingency_table <-
  rbind(c(smallRNA_RNA, RNA_not_smallRNA),
        c(smallRNA_not_RNA, not_smallRNA_not_RNA))
FT_out <- fisher.test(contingency_table)
pval_1 <- FT_out$p.value
OddsRatio_1 <- FT_out$estimate

# 2. What are the odds that genes with RNA expression changes have both INHERITED RNA expression changes and small RNA epimutations?

# background is all genes which map to small RNA and have AND RNA expression changes 
background <- nrow(Control_RNA_and_smallRNA_integrated_table[Control_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Control_RNA_and_smallRNA_integrated_table$RNA_mut==1, ])

# of which have inherited RNA changes
with_RNA_changes <- nrow(Control_RNA_and_smallRNA_integrated_table[Control_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Control_RNA_and_smallRNA_integrated_table$RNA_mut==1&Control_RNA_and_smallRNA_integrated_table$RNA_inherited==1, ])

# of which have small RNA changes of any kind
with_smallRNA_changes <- nrow(Control_RNA_and_smallRNA_integrated_table[Control_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Control_RNA_and_smallRNA_integrated_table$RNA_mut==1&Control_RNA_and_smallRNA_integrated_table$smallRNA_mut==1, ])

# of which have both inherited RNA changes and small RNA changes
with_both <- nrow(Control_RNA_and_smallRNA_integrated_table[Control_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Control_RNA_and_smallRNA_integrated_table$smallRNA_mut==1&Control_RNA_and_smallRNA_integrated_table$RNA_mut==1&Control_RNA_and_smallRNA_integrated_table$RNA_inherited==1, ])

smallRNA_RNA <-  with_both
RNA_not_smallRNA <- with_RNA_changes - smallRNA_RNA
smallRNA_not_RNA <- with_smallRNA_changes  -  smallRNA_RNA
not_smallRNA_not_RNA <- background - (with_RNA_changes + with_smallRNA_changes -  smallRNA_RNA)

contingency_table <-
  rbind(c(smallRNA_RNA, RNA_not_smallRNA),
        c(smallRNA_not_RNA, not_smallRNA_not_RNA))
FT_out <- fisher.test(contingency_table)
pval_2 <- FT_out$p.value
OddsRatio_2 <- FT_out$estimate

# Time matching
# 3. What are the odds that Genes with RNA expression changes have both inherited RNA expression changes and time matched small RNA epimutations?

# background is all genes which map to small RNA and have RNA expression changes
background <- nrow(Control_RNA_and_smallRNA_integrated_table[Control_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Control_RNA_and_smallRNA_integrated_table$RNA_mut==1, ])

# of which have inherited RNA changes
with_RNA_changes <- nrow(Control_RNA_and_smallRNA_integrated_table[Control_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Control_RNA_and_smallRNA_integrated_table$RNA_mut==1&Control_RNA_and_smallRNA_integrated_table$RNA_inherited==1, ])

# of which have time matched small RNA changes 
with_smallRNA_changes <- nrow(Control_RNA_and_smallRNA_integrated_table[Control_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Control_RNA_and_smallRNA_integrated_table$RNA_mut==1&Control_RNA_and_smallRNA_integrated_table$smallRNA_mut==1&Control_RNA_and_smallRNA_integrated_table$time_matched_to_RNA ==1, ])

# of which have both inherited RNA changes and time matched small RNA changes
with_both <- nrow(Control_RNA_and_smallRNA_integrated_table[Control_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Control_RNA_and_smallRNA_integrated_table$RNA_mut==1&Control_RNA_and_smallRNA_integrated_table$RNA_inherited==1&Control_RNA_and_smallRNA_integrated_table$smallRNA_mut==1&Control_RNA_and_smallRNA_integrated_table$time_matched_to_RNA==1, ])

smallRNA__RNA <-  with_both
RNA_not_smallRNA <- with_RNA_changes - smallRNA_RNA
smallRNA_not_RNA <- with_smallRNA_changes  -  smallRNA_RNA
not_smallRNA_not_RNA <- background - (with_RNA_changes + with_smallRNA_changes -  smallRNA_RNA)

contingency_table <-
  rbind(c(smallRNA_RNA, RNA_not_smallRNA),
        c(smallRNA_not_RNA, not_smallRNA_not_RNA))
FT_out <- fisher.test(contingency_table)
pval_3 <- FT_out$p.value
OddsRatio_3 <- FT_out$estimate

# 4. What are the odds that in Genes with RNA expression changes and small RNA epimutations, the small RNA epimutations are time matched and the RNA expression changes are inherited?

# background is all genes which have small RNA data points AND RNA expression changes and small RNA epimutations
background <- nrow(Control_RNA_and_smallRNA_integrated_table[Control_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Control_RNA_and_smallRNA_integrated_table$RNA_mut==1&Control_RNA_and_smallRNA_integrated_table$smallRNA_mut==1, ])

# of which have inherited RNA changes
with_RNA_changes <- nrow(Control_RNA_and_smallRNA_integrated_table[Control_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Control_RNA_and_smallRNA_integrated_table$smallRNA_mut==1&Control_RNA_and_smallRNA_integrated_table$RNA_mut==1&Control_RNA_and_smallRNA_integrated_table$RNA_inherited==1, ])

# of which have time matched small RNA changes
with_smallRNA_changes <- nrow(Control_RNA_and_smallRNA_integrated_table[Control_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Control_RNA_and_smallRNA_integrated_table$smallRNA_mut==1&Control_RNA_and_smallRNA_integrated_table$RNA_mut==1&Control_RNA_and_smallRNA_integrated_table$time_matched_to_RNA==1, ])

# of which have both inherited RNA changes and time matched smallRNA changes
with_both <- nrow(Control_RNA_and_smallRNA_integrated_table[Control_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Control_RNA_and_smallRNA_integrated_table$smallRNA_mut==1&Control_RNA_and_smallRNA_integrated_table$RNA_mut==1&Control_RNA_and_smallRNA_integrated_table$RNA_inherited==1&Control_RNA_and_smallRNA_integrated_table$time_matched_to_RNA==1, ])

smallRNA_RNA <-  with_both
RNA_not_smallRNA <- with_RNA_changes - smallRNA_RNA
smallRNA_not_RNA <- with_smallRNA_changes  -  smallRNA_RNA
not_smallRNA_not_RNA <- background - (with_RNA_changes + with_smallRNA_changes -  smallRNA_RNA)

contingency_table <-
  rbind(c(smallRNA_RNA, RNA_not_smallRNA),
        c(smallRNA_not_RNA, not_smallRNA_not_RNA))
FT_out <- fisher.test(contingency_table)
pval_4 <- FT_out$p.value
OddsRatio_4 <- FT_out$estimate

#5. What are the odds that genes have both inherited small RNA epimutations and inherited RNA changes ?

# background is all genes which have small RNA data points AND RNA expression changes and small RNA epimutations
background <- nrow(Control_RNA_and_smallRNA_integrated_table[Control_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Control_RNA_and_smallRNA_integrated_table$RNA_mut==1&Control_RNA_and_smallRNA_integrated_table$smallRNA_mut==1, ])

# of which have inherited RNA changes
with_RNA_changes <- nrow(Control_RNA_and_smallRNA_integrated_table[Control_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Control_RNA_and_smallRNA_integrated_table$smallRNA_mut==1&Control_RNA_and_smallRNA_integrated_table$RNA_mut==1&Control_RNA_and_smallRNA_integrated_table$RNA_inherited==1, ])

# of which have inherited small RNA changes
with_smallRNA_changes <- nrow(Control_RNA_and_smallRNA_integrated_table[Control_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Control_RNA_and_smallRNA_integrated_table$smallRNA_mut==1&Control_RNA_and_smallRNA_integrated_table$RNA_mut==1&Control_RNA_and_smallRNA_integrated_table$smallRNA_inherited==1, ])

# of which have both inherited RNA changes and time matched smallRNA changes
with_both <- nrow(Control_RNA_and_smallRNA_integrated_table[Control_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Control_RNA_and_smallRNA_integrated_table$smallRNA_mut==1&Control_RNA_and_smallRNA_integrated_table$RNA_mut==1&Control_RNA_and_smallRNA_integrated_table$RNA_inherited==1&Control_RNA_and_smallRNA_integrated_table$time_matched_to_RNA==1, ])

smallRNA_RNA <-  with_both
RNA_not_smallRNA <- with_RNA_changes - smallRNA_RNA
smallRNA_not_RNA <- with_smallRNA_changes  -  smallRNA_RNA
not_smallRNA_not_RNA <- background - (with_RNA_changes + with_smallRNA_changes -  smallRNA_RNA)

contingency_table <-
  rbind(c(smallRNA_RNA, RNA_not_smallRNA),
        c(smallRNA_not_RNA, not_smallRNA_not_RNA))
FT_out <- fisher.test(contingency_table)
pval_4 <- FT_out$p.value
OddsRatio_4 <- FT_out$estimate

# Table 2

# Create a table to show all the above results 1-4:
OR_frame <- as.data.frame(c(OddsRatio_1, OddsRatio_2, OddsRatio_3, OddsRatio_4))
Pval_frame <- as.data.frame(c(pval_1, pval_2, pval_3, pval_4))  
names <- c("Out of all genes, \nlikelihood of genes having both \nRNA epression changes and small RNA epimutations", 
           "Out of genes with RNAseq changes, \nlikelihood of genes having both \ninherited RNA expression changes and small RNA epimutations", 
           "Out of genes with RNAseq changes, \nlikelihood of genes having both \ninherited RNA expression changes and simultaneous small RNA epimutations",
           "Out of genes with RNAseq changes and small RNA changes, \nlikelihood of genes having both \n inherited RNA expression changes and simultaneous small RNA epimutations")


Table_2 <- cbind(names, OR_frame, Pval_frame)
colnames(Table_2) <- c("Association_tested", "Odds Ratio", "p-value")


#Odds calculation Low dose

# 1. What are the odds that genes have both small RNA epimutations and RNA changes ?

# background is all genes which have small RNA signal, so have to map to small RNA
background <- nrow(Low_dose_RNA_and_smallRNA_integrated_table[Low_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1, ])

# of which have RNA changes
with_RNA_changes <- nrow(Low_dose_RNA_and_smallRNA_integrated_table[Low_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Low_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1, ])

# of which have small RNA changes 
with_smallRNA_changes <- nrow(Low_dose_RNA_and_smallRNA_integrated_table[Low_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Low_dose_RNA_and_smallRNA_integrated_table$smallRNA_mut==1, ])

# of which have both 
with_both <- nrow(Low_dose_RNA_and_smallRNA_integrated_table[Low_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Low_dose_RNA_and_smallRNA_integrated_table$smallRNA_mut==1&Low_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1, ])

smallRNA_RNA <-  with_both
RNA_not_smallRNA <- with_RNA_changes - smallRNA_RNA
smallRNA_not_RNA <- with_smallRNA_changes  -  smallRNA_RNA
not_smallRNA_not_RNA <- background - (with_RNA_changes + with_smallRNA_changes -  smallRNA_RNA)

contingency_table <-
  rbind(c(smallRNA_RNA, RNA_not_smallRNA),
        c(smallRNA_not_RNA, not_smallRNA_not_RNA))
FT_out <- fisher.test(contingency_table)
pval_1 <- FT_out$p.value
OddsRatio_1 <- FT_out$estimate

# 2. What are the odds that genes with RNA expression changes have both INHERITED RNA expression changes and small RNA epimutations?

# background is all genes which map to small RNA and have AND RNA expression changes 
background <- nrow(Low_dose_RNA_and_smallRNA_integrated_table[Low_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Low_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1, ])

# of which have inherited RNA changes
with_RNA_changes <- nrow(Low_dose_RNA_and_smallRNA_integrated_table[Low_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Low_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1&Low_dose_RNA_and_smallRNA_integrated_table$RNA_inherited==1, ])

# of which have small RNA changes of any kind
with_smallRNA_changes <- nrow(Low_dose_RNA_and_smallRNA_integrated_table[Low_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Low_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1&Low_dose_RNA_and_smallRNA_integrated_table$smallRNA_mut==1, ])

# of which have both inherited RNA changes and small RNA changes
with_both <- nrow(Low_dose_RNA_and_smallRNA_integrated_table[Low_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Low_dose_RNA_and_smallRNA_integrated_table$smallRNA_mut==1&Low_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1&Low_dose_RNA_and_smallRNA_integrated_table$RNA_inherited==1, ])

smallRNA_RNA <-  with_both
RNA_not_smallRNA <- with_RNA_changes - smallRNA_RNA
smallRNA_not_RNA <- with_smallRNA_changes  -  smallRNA_RNA
not_smallRNA_not_RNA <- background - (with_RNA_changes + with_smallRNA_changes -  smallRNA_RNA)

contingency_table <-
  rbind(c(smallRNA_RNA, RNA_not_smallRNA),
        c(smallRNA_not_RNA, not_smallRNA_not_RNA))
FT_out <- fisher.test(contingency_table)
pval_2 <- FT_out$p.value
OddsRatio_2 <- FT_out$estimate

# Time matching
# 3. What are the odds that Genes with RNA expression changes have both inherited RNA expression changes and time matched small RNA epimutations?

# background is all genes which map to small RNA and have RNA expression changes
background <- nrow(Low_dose_RNA_and_smallRNA_integrated_table[Low_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Low_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1, ])

# of which have inherited RNA changes
with_RNA_changes <- nrow(Low_dose_RNA_and_smallRNA_integrated_table[Low_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Low_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1&Low_dose_RNA_and_smallRNA_integrated_table$RNA_inherited==1, ])

# of which have time matched small RNA changes 
with_smallRNA_changes <- nrow(Low_dose_RNA_and_smallRNA_integrated_table[Low_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Low_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1&Low_dose_RNA_and_smallRNA_integrated_table$smallRNA_mut==1&Low_dose_RNA_and_smallRNA_integrated_table$time_matched_to_RNA ==1, ])

# of which have both inherited RNA changes and time matched small RNA changes
with_both <- nrow(Low_dose_RNA_and_smallRNA_integrated_table[Low_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Low_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1&Low_dose_RNA_and_smallRNA_integrated_table$RNA_inherited==1&Low_dose_RNA_and_smallRNA_integrated_table$smallRNA_mut==1&Low_dose_RNA_and_smallRNA_integrated_table$time_matched_to_RNA==1, ])

smallRNA__RNA <-  with_both
RNA_not_smallRNA <- with_RNA_changes - smallRNA_RNA
smallRNA_not_RNA <- with_smallRNA_changes  -  smallRNA_RNA
not_smallRNA_not_RNA <- background - (with_RNA_changes + with_smallRNA_changes -  smallRNA_RNA)

contingency_table <-
  rbind(c(smallRNA_RNA, RNA_not_smallRNA),
        c(smallRNA_not_RNA, not_smallRNA_not_RNA))
FT_out <- fisher.test(contingency_table)
pval_3 <- FT_out$p.value
OddsRatio_3 <- FT_out$estimate

# 4. What are the odds that in Genes with RNA expression changes and small RNA epimutations, the small RNA epimutations are time matched and the RNA expression changes are inherited?

# background is all genes which have small RNA data points AND RNA expression changes and small RNA epimutations
background <- nrow(Low_dose_RNA_and_smallRNA_integrated_table[Low_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Low_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1&Low_dose_RNA_and_smallRNA_integrated_table$smallRNA_mut==1, ])

# of which have inherited RNA changes
with_RNA_changes <- nrow(Low_dose_RNA_and_smallRNA_integrated_table[Low_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Low_dose_RNA_and_smallRNA_integrated_table$smallRNA_mut==1&Low_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1&Low_dose_RNA_and_smallRNA_integrated_table$RNA_inherited==1, ])

# of which have time matched small RNA changes
with_smallRNA_changes <- nrow(Low_dose_RNA_and_smallRNA_integrated_table[Low_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Low_dose_RNA_and_smallRNA_integrated_table$smallRNA_mut==1&Low_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1&Low_dose_RNA_and_smallRNA_integrated_table$time_matched_to_RNA==1, ])

# of which have both inherited RNA changes and time matched smallRNA changes
with_both <- nrow(Low_dose_RNA_and_smallRNA_integrated_table[Low_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Low_dose_RNA_and_smallRNA_integrated_table$smallRNA_mut==1&Low_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1&Low_dose_RNA_and_smallRNA_integrated_table$RNA_inherited==1&Low_dose_RNA_and_smallRNA_integrated_table$time_matched_to_RNA==1, ])

smallRNA_RNA <-  with_both
RNA_not_smallRNA <- with_RNA_changes - smallRNA_RNA
smallRNA_not_RNA <- with_smallRNA_changes  -  smallRNA_RNA
not_smallRNA_not_RNA <- background - (with_RNA_changes + with_smallRNA_changes -  smallRNA_RNA)

contingency_table <-
  rbind(c(smallRNA_RNA, RNA_not_smallRNA),
        c(smallRNA_not_RNA, not_smallRNA_not_RNA))
FT_out <- fisher.test(contingency_table)
pval_4 <- FT_out$p.value
OddsRatio_4 <- FT_out$estimate

# Table 2

# Create a table to show all the above results 1-4:
OR_frame <- as.data.frame(c(OddsRatio_1, OddsRatio_2, OddsRatio_3, OddsRatio_4))
Pval_frame <- as.data.frame(c(pval_1, pval_2, pval_3, pval_4))  

names <- c("Out of all genes, \nlikelihood of genes having both \nRNA epression changes and small RNA epimutations", 
           "Out of genes with RNAseq changes, \nlikelihood of genes having both \ninherited RNA expression changes and small RNA epimutations", 
           "Out of genes with RNAseq changes, \nlikelihood of genes having both \ninherited RNA expression changes and simultaneous small RNA epimutations",
           "Out of genes with RNAseq changes and small RNA changes, \nlikelihood of genes having both \n inherited RNA expression changes and simultaneous small RNA epimutations")


Table_2 <- cbind(names, OR_frame, Pval_frame)
colnames(Table_2) <- c("Association_tested", "Odds Ratio", "p-value")


#Odds calculation High dose

# 1. What are the odds that genes have both small RNA epimutations and RNA changes ?

# background is all genes which have small RNA signal, so have to map to small RNA
background <- nrow(High_dose_RNA_and_smallRNA_integrated_table[High_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1, ])

# of which have RNA changes
with_RNA_changes <- nrow(High_dose_RNA_and_smallRNA_integrated_table[High_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&High_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1, ])

# of which have small RNA changes 
with_smallRNA_changes <- nrow(High_dose_RNA_and_smallRNA_integrated_table[High_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&High_dose_RNA_and_smallRNA_integrated_table$smallRNA_mut==1, ])

# of which have both 
with_both <- nrow(High_dose_RNA_and_smallRNA_integrated_table[High_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&High_dose_RNA_and_smallRNA_integrated_table$smallRNA_mut==1&High_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1, ])

smallRNA_RNA <-  with_both
RNA_not_smallRNA <- with_RNA_changes - smallRNA_RNA
smallRNA_not_RNA <- with_smallRNA_changes  -  smallRNA_RNA
not_smallRNA_not_RNA <- background - (with_RNA_changes + with_smallRNA_changes -  smallRNA_RNA)

contingency_table <-
  rbind(c(smallRNA_RNA, RNA_not_smallRNA),
        c(smallRNA_not_RNA, not_smallRNA_not_RNA))
FT_out <- fisher.test(contingency_table)
pval_1 <- FT_out$p.value
OddsRatio_1 <- FT_out$estimate

# 2. What are the odds that genes with RNA expression changes have both INHERITED RNA expression changes and small RNA epimutations?

# background is all genes which map to small RNA and have AND RNA expression changes 
background <- nrow(High_dose_RNA_and_smallRNA_integrated_table[High_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&High_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1, ])

# of which have inherited RNA changes
with_RNA_changes <- nrow(High_dose_RNA_and_smallRNA_integrated_table[High_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&High_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1&High_dose_RNA_and_smallRNA_integrated_table$RNA_inherited==1, ])

# of which have small RNA changes of any kind
with_smallRNA_changes <- nrow(High_dose_RNA_and_smallRNA_integrated_table[High_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&High_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1&High_dose_RNA_and_smallRNA_integrated_table$smallRNA_mut==1, ])

# of which have both inherited RNA changes and small RNA changes
with_both <- nrow(High_dose_RNA_and_smallRNA_integrated_table[High_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&High_dose_RNA_and_smallRNA_integrated_table$smallRNA_mut==1&High_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1&High_dose_RNA_and_smallRNA_integrated_table$RNA_inherited==1, ])

smallRNA_RNA <-  with_both
RNA_not_smallRNA <- with_RNA_changes - smallRNA_RNA
smallRNA_not_RNA <- with_smallRNA_changes  -  smallRNA_RNA
not_smallRNA_not_RNA <- background - (with_RNA_changes + with_smallRNA_changes -  smallRNA_RNA)

contingency_table <-
  rbind(c(smallRNA_RNA, RNA_not_smallRNA),
        c(smallRNA_not_RNA, not_smallRNA_not_RNA))
FT_out <- fisher.test(contingency_table)
pval_2 <- FT_out$p.value
OddsRatio_2 <- FT_out$estimate

# Time matching
# 3. What are the odds that Genes with RNA expression changes have both inherited RNA expression changes and time matched small RNA epimutations?

# background is all genes which map to small RNA and have RNA expression changes
background <- nrow(High_dose_RNA_and_smallRNA_integrated_table[High_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&High_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1, ])

# of which have inherited RNA changes
with_RNA_changes <- nrow(High_dose_RNA_and_smallRNA_integrated_table[High_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&High_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1&High_dose_RNA_and_smallRNA_integrated_table$RNA_inherited==1, ])

# of which have time matched small RNA changes 
with_smallRNA_changes <- nrow(High_dose_RNA_and_smallRNA_integrated_table[High_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&High_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1&High_dose_RNA_and_smallRNA_integrated_table$smallRNA_mut==1&High_dose_RNA_and_smallRNA_integrated_table$time_matched_to_RNA ==1, ])

# of which have both inherited RNA changes and time matched small RNA changes
with_both <- nrow(High_dose_RNA_and_smallRNA_integrated_table[High_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&High_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1&High_dose_RNA_and_smallRNA_integrated_table$RNA_inherited==1&High_dose_RNA_and_smallRNA_integrated_table$smallRNA_mut==1&High_dose_RNA_and_smallRNA_integrated_table$time_matched_to_RNA==1, ])

smallRNA__RNA <-  with_both
RNA_not_smallRNA <- with_RNA_changes - smallRNA_RNA
smallRNA_not_RNA <- with_smallRNA_changes  -  smallRNA_RNA
not_smallRNA_not_RNA <- background - (with_RNA_changes + with_smallRNA_changes -  smallRNA_RNA)

contingency_table <-
  rbind(c(smallRNA_RNA, RNA_not_smallRNA),
        c(smallRNA_not_RNA, not_smallRNA_not_RNA))
FT_out <- fisher.test(contingency_table)
pval_3 <- FT_out$p.value
OddsRatio_3 <- FT_out$estimate

# 4. What are the odds that in Genes with RNA expression changes and small RNA epimutations, the small RNA epimutations are time matched and the RNA expression changes are inherited?

# background is all genes which have small RNA data points AND RNA expression changes and small RNA epimutations
background <- nrow(High_dose_RNA_and_smallRNA_integrated_table[High_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&High_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1&High_dose_RNA_and_smallRNA_integrated_table$smallRNA_mut==1, ])

# of which have inherited RNA changes
with_RNA_changes <- nrow(High_dose_RNA_and_smallRNA_integrated_table[High_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&High_dose_RNA_and_smallRNA_integrated_table$smallRNA_mut==1&High_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1&High_dose_RNA_and_smallRNA_integrated_table$RNA_inherited==1, ])

# of which have time matched small RNA changes
with_smallRNA_changes <- nrow(Low_dose_RNA_and_smallRNA_integrated_table[Low_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Low_dose_RNA_and_smallRNA_integrated_table$smallRNA_mut==1&Low_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1&Low_dose_RNA_and_smallRNA_integrated_table$time_matched_to_RNA==1, ])

# of which have both inherited RNA changes and time matched smallRNA changes
with_both <- nrow(High_dose_RNA_and_smallRNA_integrated_table[High_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&High_dose_RNA_and_smallRNA_integrated_table$smallRNA_mut==1&High_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1&High_dose_RNA_and_smallRNA_integrated_table$RNA_inherited==1&High_dose_RNA_and_smallRNA_integrated_table$time_matched_to_RNA==1, ])

smallRNA_RNA <-  with_both
RNA_not_smallRNA <- with_RNA_changes - smallRNA_RNA
smallRNA_not_RNA <- with_smallRNA_changes  -  smallRNA_RNA
not_smallRNA_not_RNA <- background - (with_RNA_changes + with_smallRNA_changes -  smallRNA_RNA)

contingency_table <-
  rbind(c(smallRNA_RNA, RNA_not_smallRNA),
        c(smallRNA_not_RNA, not_smallRNA_not_RNA))
FT_out <- fisher.test(contingency_table)
pval_4 <- FT_out$p.value
OddsRatio_4 <- FT_out$estimate

# Create a table to show all the above results 1-4:
OR_frame <- as.data.frame(c(OddsRatio_1, OddsRatio_2, OddsRatio_3, OddsRatio_4))
Pval_frame <- as.data.frame(c(pval_1, pval_2, pval_3, pval_4))  

names <- c("Out of all genes, \nlikelihood of genes having both \nRNA epression changes and small RNA epimutations", 
           "Out of genes with RNAseq changes, \nlikelihood of genes having both \ninherited RNA expression changes and small RNA epimutations", 
           "Out of genes with RNAseq changes, \nlikelihood of genes having both \ninherited RNA expression changes and simultaneous small RNA epimutations",
           "Out of genes with RNAseq changes and small RNA changes, \nlikelihood of genes having both \n inherited RNA expression changes and simultaneous small RNA epimutations")

Table_2 <- cbind(names, OR_frame, Pval_frame)
colnames(Table_2) <- c("Association_tested", "Odds Ratio", "p-value")
#--------------------------
#####Fig.4.C & Sup.Fig.6.C - piRNAs epimutations####
load("tinyRNA_C1epimutations.Rdata")
load("tinyRNA_C2epimutations.Rdata")
load("tinyRNA_L1epimutations.Rdata")
load("tinyRNA_L2epimutations.Rdata")
load("tinyRNA_H1epimutations.Rdata")
load("tinyRNA_H2epimutations.Rdata")
#Control1 new epi
tinyRNA_C1epimutations_piRNAs<-subset(tinyRNA_C1epimutations,tinyRNA_C1epimutations$Feature.Class=="piRNA")
#UP
tinyRNA_piRNAs_UP_output_C1<- c()
tinyRNA_ID <- paste(tinyRNA_C1epimutations_piRNAs$Feature.ID, tinyRNA_C1epimutations_piRNAs$Tag,tinyRNA_C1epimutations_piRNAs$Feature.Name,tinyRNA_C1epimutations_piRNAs$Feature.Class, sep="_")
row.names(tinyRNA_C1epimutations_piRNAs)<-tinyRNA_ID
tinyRNA_C1epimutations_piRNAs<-tinyRNA_C1epimutations_piRNAs[,c(7:17)]
for(i in 1:nrow(tinyRNA_C1epimutations_piRNAs)){
  tinyRNA_piRNAs_UP_output_C1 <-rbind(tinyRNA_piRNAs_UP_output_C1, UP_transition_func(tinyRNA_C1epimutations_piRNAs[i,], input_name=row.names(tinyRNA_C1epimutations_piRNAs)[i]))}
colnames(tinyRNA_piRNAs_UP_output_C1) <- c("0", "2", "4", "6", "8", "10", "12", "14", "16", "18", "20")
row.names(tinyRNA_piRNAs_UP_output_C1) <- rownames(tinyRNA_C1epimutations_piRNAs)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tinyRNA_piRNAs_UP_output_C1[,2])
transitions_at_4 <- sum(tinyRNA_piRNAs_UP_output_C1[,3])
transitions_at_6 <- sum(tinyRNA_piRNAs_UP_output_C1[,4])
transitions_at_8 <- sum(tinyRNA_piRNAs_UP_output_C1[,5])
transitions_at_10 <- sum(tinyRNA_piRNAs_UP_output_C1[,6])
transitions_at_12 <- sum(tinyRNA_piRNAs_UP_output_C1[,7])
transitions_at_14 <- sum(tinyRNA_piRNAs_UP_output_C1[,8])
transitions_at_16 <- sum(tinyRNA_piRNAs_UP_output_C1[,9])
transitions_at_18 <- sum(tinyRNA_piRNAs_UP_output_C1[,10])
transitions_at_20 <- sum(tinyRNA_piRNAs_UP_output_C1[,11])

tinyRNA_piRNAs_UP_output_C1_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                               transitions_at_4, 
                                                               transitions_at_6,
                                                               transitions_at_8, 
                                                               transitions_at_10, 
                                                               transitions_at_12, 
                                                               transitions_at_14, 
                                                               transitions_at_16, 
                                                               transitions_at_18, 
                                                               transitions_at_20)

#DOWN
tinyRNA_piRNAs_DOWN_output_C1<- c()
for(i in 1:nrow(tinyRNA_C1epimutations_piRNAs)){
  tinyRNA_piRNAs_DOWN_output_C1 <-rbind(tinyRNA_piRNAs_DOWN_output_C1, DOWN_transition_func(tinyRNA_C1epimutations_piRNAs[i,], input_name=row.names(tinyRNA_C1epimutations_piRNAs)[i]))}
colnames(tinyRNA_piRNAs_DOWN_output_C1) <- c("0", "2", "4", "6", "8", "10", "12", "14", "16", "18", "20")
row.names(tinyRNA_piRNAs_DOWN_output_C1) <- rownames(tinyRNA_C1epimutations_piRNAs)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tinyRNA_piRNAs_DOWN_output_C1[,2])
transitions_at_4 <- sum(tinyRNA_piRNAs_DOWN_output_C1[,3])
transitions_at_6 <- sum(tinyRNA_piRNAs_DOWN_output_C1[,4])
transitions_at_8 <- sum(tinyRNA_piRNAs_DOWN_output_C1[,5])
transitions_at_10 <- sum(tinyRNA_piRNAs_DOWN_output_C1[,6])
transitions_at_12 <- sum(tinyRNA_piRNAs_DOWN_output_C1[,7])
transitions_at_14 <- sum(tinyRNA_piRNAs_DOWN_output_C1[,8])
transitions_at_16 <- sum(tinyRNA_piRNAs_DOWN_output_C1[,9])
transitions_at_18 <- sum(tinyRNA_piRNAs_DOWN_output_C1[,10])
transitions_at_20 <- sum(tinyRNA_piRNAs_DOWN_output_C1[,11])

tinyRNA_piRNAs_DOWN_output_C1_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                                 transitions_at_4, 
                                                                 transitions_at_6,
                                                                 transitions_at_8, 
                                                                 transitions_at_10, 
                                                                 transitions_at_12, 
                                                                 transitions_at_14, 
                                                                 transitions_at_16, 
                                                                 transitions_at_18, 
                                                                 transitions_at_20)

#Control2 new epi
tinyRNA_C2epimutations_piRNAs<-subset(tinyRNA_C2epimutations,tinyRNA_C2epimutations$Feature.Class=="piRNA")
#UP
tinyRNA_piRNAs_UP_output_C2<- c()
tinyRNA_ID <- paste(tinyRNA_C2epimutations_piRNAs$Feature.ID, tinyRNA_C2epimutations_piRNAs$Tag,tinyRNA_C2epimutations_piRNAs$Feature.Name,tinyRNA_C2epimutations_piRNAs$Feature.Class, sep="_")
row.names(tinyRNA_C2epimutations_piRNAs)<-tinyRNA_ID
tinyRNA_C2epimutations_piRNAs<-tinyRNA_C2epimutations_piRNAs[,c(7:16)]
for(i in 1:nrow(tinyRNA_C2epimutations_piRNAs)){
  tinyRNA_piRNAs_UP_output_C2 <-rbind(tinyRNA_piRNAs_UP_output_C2, UP_transition_func(tinyRNA_C2epimutations_piRNAs[i,], input_name=row.names(tinyRNA_C2epimutations_piRNAs)[i]))}
colnames(tinyRNA_piRNAs_UP_output_C2) <- c("0", "2", "4","8", "10", "12", "14", "16", "18", "20")
row.names(tinyRNA_piRNAs_UP_output_C2) <- rownames(tinyRNA_C2epimutations_piRNAs)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tinyRNA_piRNAs_UP_output_C2[,2])
transitions_at_4 <- sum(tinyRNA_piRNAs_UP_output_C2[,3])
transitions_at_6 <- NA
transitions_at_8 <- sum(tinyRNA_piRNAs_UP_output_C2[,4])
transitions_at_10 <- sum(tinyRNA_piRNAs_UP_output_C2[,5])
transitions_at_12 <- sum(tinyRNA_piRNAs_UP_output_C2[,6])
transitions_at_14 <- sum(tinyRNA_piRNAs_UP_output_C2[,7])
transitions_at_16 <- sum(tinyRNA_piRNAs_UP_output_C2[,8])
transitions_at_18 <- sum(tinyRNA_piRNAs_UP_output_C2[,9])
transitions_at_20 <- sum(tinyRNA_piRNAs_UP_output_C2[,10])

tinyRNA_piRNAs_UP_output_C2_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                               transitions_at_4, 
                                                               transitions_at_6,
                                                               transitions_at_8, 
                                                               transitions_at_10, 
                                                               transitions_at_12, 
                                                               transitions_at_14, 
                                                               transitions_at_16, 
                                                               transitions_at_18, 
                                                               transitions_at_20)

#DOWN
tinyRNA_piRNAs_DOWN_output_C2<- c()
for(i in 1:nrow(tinyRNA_C2epimutations_piRNAs)){
  tinyRNA_piRNAs_DOWN_output_C2 <-rbind(tinyRNA_piRNAs_DOWN_output_C2, DOWN_transition_func(tinyRNA_C2epimutations_piRNAs[i,], input_name=row.names(tinyRNA_C2epimutations_piRNAs)[i]))}
colnames(tinyRNA_piRNAs_DOWN_output_C2) <- c("0", "2", "4", "8", "10", "12", "14", "16", "18", "20")
row.names(tinyRNA_piRNAs_DOWN_output_C2) <- rownames(tinyRNA_C2epimutations_piRNAs)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tinyRNA_piRNAs_DOWN_output_C2[,2])
transitions_at_4 <- sum(tinyRNA_piRNAs_DOWN_output_C2[,3])
transitions_at_6 <- NA
transitions_at_8 <- sum(tinyRNA_piRNAs_DOWN_output_C2[,4])
transitions_at_10 <- sum(tinyRNA_piRNAs_DOWN_output_C2[,5])
transitions_at_12 <- sum(tinyRNA_piRNAs_DOWN_output_C2[,6])
transitions_at_14 <- sum(tinyRNA_piRNAs_DOWN_output_C2[,7])
transitions_at_16 <- sum(tinyRNA_piRNAs_DOWN_output_C2[,8])
transitions_at_18 <- sum(tinyRNA_piRNAs_DOWN_output_C2[,9])
transitions_at_20 <- sum(tinyRNA_piRNAs_DOWN_output_C2[,10])

tinyRNA_piRNAs_DOWN_output_C2_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                                 transitions_at_4, 
                                                                 transitions_at_6,
                                                                 transitions_at_8, 
                                                                 transitions_at_10, 
                                                                 transitions_at_12, 
                                                                 transitions_at_14, 
                                                                 transitions_at_16, 
                                                                 transitions_at_18, 
                                                                 transitions_at_20)

#Low1 new epi
tinyRNA_L1epimutations_piRNAs<-subset(tinyRNA_L1epimutations,tinyRNA_L1epimutations$Feature.Class=="piRNA")
#UP

tinyRNA_piRNAs_UP_output_L1<- c()
tinyRNA_ID <- paste(tinyRNA_L1epimutations_piRNAs$Feature.ID, tinyRNA_L1epimutations_piRNAs$Tag,tinyRNA_L1epimutations_piRNAs$Feature.Name,tinyRNA_L1epimutations_piRNAs$Feature.Class, sep="_")
row.names(tinyRNA_L1epimutations_piRNAs)<-tinyRNA_ID
tinyRNA_L1epimutations_piRNAs<-tinyRNA_L1epimutations_piRNAs[,c(7:16)]
for(i in 1:nrow(tinyRNA_L1epimutations_piRNAs)){
  tinyRNA_piRNAs_UP_output_L1 <-rbind(tinyRNA_piRNAs_UP_output_L1, UP_transition_func(tinyRNA_L1epimutations_piRNAs[i,], input_name=row.names(tinyRNA_L1epimutations_piRNAs)[i]))}
colnames(tinyRNA_piRNAs_UP_output_L1) <- c("0", "2", "4","6", "8", "12", "14", "16", "18", "20")
row.names(tinyRNA_piRNAs_UP_output_L1) <- rownames(tinyRNA_L1epimutations_piRNAs)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tinyRNA_piRNAs_UP_output_L1[,2])
transitions_at_4 <- sum(tinyRNA_piRNAs_UP_output_L1[,3])
transitions_at_6 <- sum(tinyRNA_piRNAs_UP_output_L1[,4])
transitions_at_8 <- sum(tinyRNA_piRNAs_UP_output_L1[,5])
transitions_at_10 <- NA
transitions_at_12 <- sum(tinyRNA_piRNAs_UP_output_L1[,6])
transitions_at_14 <- sum(tinyRNA_piRNAs_UP_output_L1[,7])
transitions_at_16 <- sum(tinyRNA_piRNAs_UP_output_L1[,8])
transitions_at_18 <- sum(tinyRNA_piRNAs_UP_output_L1[,9])
transitions_at_20 <- sum(tinyRNA_piRNAs_UP_output_L1[,10])

tinyRNA_piRNAs_UP_L1_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                        transitions_at_4, 
                                                        transitions_at_6,
                                                        transitions_at_8, 
                                                        transitions_at_10, 
                                                        transitions_at_12, 
                                                        transitions_at_14, 
                                                        transitions_at_16, 
                                                        transitions_at_18, 
                                                        transitions_at_20)

#DOWN
tinyRNA_piRNAs_DOWN_output_L1<- c()
for(i in 1:nrow(tinyRNA_L1epimutations_piRNAs)){
  tinyRNA_piRNAs_DOWN_output_L1 <-rbind(tinyRNA_piRNAs_DOWN_output_L1, DOWN_transition_func(tinyRNA_L1epimutations_piRNAs[i,], input_name=row.names(tinyRNA_L1epimutations_piRNAs)[i]))}
colnames(tinyRNA_piRNAs_DOWN_output_L1) <- c("0", "2", "4", "6","8", "12", "14", "16", "18", "20")
row.names(tinyRNA_piRNAs_DOWN_output_L1) <- rownames(tinyRNA_L1epimutations_piRNAs)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tinyRNA_piRNAs_DOWN_output_L1[,2])
transitions_at_4 <- sum(tinyRNA_piRNAs_DOWN_output_L1[,3])
transitions_at_6 <- sum(tinyRNA_piRNAs_DOWN_output_L1[,4])
transitions_at_8 <- sum(tinyRNA_piRNAs_DOWN_output_L1[,5])
transitions_at_10 <- NA
transitions_at_12 <- sum(tinyRNA_piRNAs_DOWN_output_L1[,6])
transitions_at_14 <- sum(tinyRNA_piRNAs_DOWN_output_L1[,7])
transitions_at_16 <- sum(tinyRNA_piRNAs_DOWN_output_L1[,8])
transitions_at_18 <- sum(tinyRNA_piRNAs_DOWN_output_L1[,9])
transitions_at_20 <- sum(tinyRNA_piRNAs_DOWN_output_L1[,10])

tinyRNA_piRNAs_DOWN_L1_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                          transitions_at_4, 
                                                          transitions_at_6,
                                                          transitions_at_8, 
                                                          transitions_at_10, 
                                                          transitions_at_12, 
                                                          transitions_at_14, 
                                                          transitions_at_16, 
                                                          transitions_at_18, 
                                                          transitions_at_20)

#Low2 new epi
tinyRNA_L2epimutations_piRNAs<-subset(tinyRNA_L2epimutations,tinyRNA_L2epimutations$Feature.Class=="piRNA")
#UP

tinyRNA_piRNAs_UP_output_L2<- c()
tinyRNA_ID <- paste(tinyRNA_L2epimutations_piRNAs$Feature.ID, tinyRNA_L2epimutations_piRNAs$Tag,tinyRNA_L2epimutations_piRNAs$Feature.Name,tinyRNA_L2epimutations_piRNAs$Feature.Class, sep="_")
row.names(tinyRNA_L2epimutations_piRNAs)<-tinyRNA_ID
tinyRNA_L2epimutations_piRNAs<-tinyRNA_L2epimutations_piRNAs[,c(7:16)]
for(i in 1:nrow(tinyRNA_L2epimutations_piRNAs)){
  tinyRNA_piRNAs_UP_output_L2 <-rbind(tinyRNA_piRNAs_UP_output_L2, UP_transition_func(tinyRNA_L2epimutations_piRNAs[i,], input_name=row.names(tinyRNA_L2epimutations_piRNAs)[i]))}
colnames(tinyRNA_piRNAs_UP_output_L2) <- c("0", "2", "4", "8","10", "12", "14", "16", "18", "20")
row.names(tinyRNA_piRNAs_UP_output_L2) <- rownames(tinyRNA_L2epimutations_piRNAs)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tinyRNA_piRNAs_UP_output_L2[,2])
transitions_at_4 <- sum(tinyRNA_piRNAs_UP_output_L2[,3])
transitions_at_6 <- NA
transitions_at_8 <- sum(tinyRNA_piRNAs_UP_output_L2[,4])
transitions_at_10 <- sum(tinyRNA_piRNAs_UP_output_L2[,5])
transitions_at_12 <- sum(tinyRNA_piRNAs_UP_output_L2[,6])
transitions_at_14 <- sum(tinyRNA_piRNAs_UP_output_L2[,7])
transitions_at_16 <- sum(tinyRNA_piRNAs_UP_output_L2[,8])
transitions_at_18 <- sum(tinyRNA_piRNAs_UP_output_L2[,9])
transitions_at_20 <- sum(tinyRNA_piRNAs_UP_output_L2[,10])

tinyRNA_piRNAs_UP_L2_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                        transitions_at_4, 
                                                        transitions_at_6,
                                                        transitions_at_8, 
                                                        transitions_at_10, 
                                                        transitions_at_12, 
                                                        transitions_at_14, 
                                                        transitions_at_16, 
                                                        transitions_at_18, 
                                                        transitions_at_20)


#DOWN
tinyRNA_piRNAs_DOWN_output_L2<- c()
for(i in 1:nrow(tinyRNA_L2epimutations_piRNAs)){
  tinyRNA_piRNAs_DOWN_output_L2 <-rbind(tinyRNA_piRNAs_DOWN_output_L2, DOWN_transition_func(tinyRNA_L2epimutations_piRNAs[i,], input_name=row.names(tinyRNA_L2epimutations_piRNAs)[i]))}
colnames(tinyRNA_piRNAs_DOWN_output_L2) <- c("0", "2", "4", "6","8", "12", "14", "16", "18", "20")
row.names(tinyRNA_piRNAs_DOWN_output_L2) <- rownames(tinyRNA_L2epimutations_piRNAs)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tinyRNA_piRNAs_DOWN_output_L2[,2])
transitions_at_4 <- sum(tinyRNA_piRNAs_DOWN_output_L2[,3])
transitions_at_6 <- NA
transitions_at_8 <- sum(tinyRNA_piRNAs_DOWN_output_L2[,4])
transitions_at_10 <- sum(tinyRNA_piRNAs_DOWN_output_L2[,5])
transitions_at_12 <- sum(tinyRNA_piRNAs_DOWN_output_L2[,6])
transitions_at_14 <- sum(tinyRNA_piRNAs_DOWN_output_L2[,7])
transitions_at_16 <- sum(tinyRNA_piRNAs_DOWN_output_L2[,8])
transitions_at_18 <- sum(tinyRNA_piRNAs_DOWN_output_L2[,9])
transitions_at_20 <- sum(tinyRNA_piRNAs_DOWN_output_L2[,10])

tinyRNA_piRNAs_DOWN_L2_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                          transitions_at_4, 
                                                          transitions_at_6,
                                                          transitions_at_8, 
                                                          transitions_at_10, 
                                                          transitions_at_12, 
                                                          transitions_at_14, 
                                                          transitions_at_16, 
                                                          transitions_at_18, 
                                                          transitions_at_20)

#High1 new epi
tinyRNA_H1epimutations_piRNAs<-subset(tinyRNA_H1epimutations,tinyRNA_H1epimutations$Feature.Class=="piRNA")
#UP

tinyRNA_piRNAs_UP_output_H1<- c()
tinyRNA_ID <- paste(tinyRNA_H1epimutations_piRNAs$Feature.ID, tinyRNA_H1epimutations_piRNAs$Tag,tinyRNA_H1epimutations_piRNAs$Feature.Name,tinyRNA_H1epimutations_piRNAs$Feature.Class, sep="_")
row.names(tinyRNA_H1epimutations_piRNAs)<-tinyRNA_ID
tinyRNA_H1epimutations_piRNAs<-tinyRNA_H1epimutations_piRNAs[,c(7:16)]
for(i in 1:nrow(tinyRNA_H1epimutations_piRNAs)){
  tinyRNA_piRNAs_UP_output_H1 <-rbind(tinyRNA_piRNAs_UP_output_H1, UP_transition_func(tinyRNA_H1epimutations_piRNAs[i,], input_name=row.names(tinyRNA_H1epimutations_piRNAs)[i]))}
colnames(tinyRNA_piRNAs_UP_output_H1) <- c("0", "2", "4", "8","10", "12", "14", "16", "18", "20")
row.names(tinyRNA_piRNAs_UP_output_H1) <- rownames(tinyRNA_H1epimutations_piRNAs)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tinyRNA_piRNAs_UP_output_H1[,2])
transitions_at_4 <- sum(tinyRNA_piRNAs_UP_output_H1[,3])
transitions_at_6 <- NA
transitions_at_8 <- sum(tinyRNA_piRNAs_UP_output_H1[,4])
transitions_at_10 <- sum(tinyRNA_piRNAs_UP_output_H1[,5])
transitions_at_12 <- sum(tinyRNA_piRNAs_UP_output_H1[,6])
transitions_at_14 <- sum(tinyRNA_piRNAs_UP_output_H1[,7])
transitions_at_16 <- sum(tinyRNA_piRNAs_UP_output_H1[,8])
transitions_at_18 <- sum(tinyRNA_piRNAs_UP_output_H1[,9])
transitions_at_20 <- sum(tinyRNA_piRNAs_UP_output_H1[,10])

tinyRNA_piRNAs_UP_H1_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                        transitions_at_4, 
                                                        transitions_at_6,
                                                        transitions_at_8, 
                                                        transitions_at_10, 
                                                        transitions_at_12, 
                                                        transitions_at_14, 
                                                        transitions_at_16, 
                                                        transitions_at_18, 
                                                        transitions_at_20)


#DOWN
tinyRNA_piRNAs_DOWN_output_H1<- c()
for(i in 1:nrow(tinyRNA_H1epimutations_piRNAs)){
  tinyRNA_piRNAs_DOWN_output_H1 <-rbind(tinyRNA_piRNAs_DOWN_output_H1, DOWN_transition_func(tinyRNA_H1epimutations_piRNAs[i,], input_name=row.names(tinyRNA_H1epimutations_piRNAs)[i]))}
colnames(tinyRNA_piRNAs_DOWN_output_H1) <- c("0", "2", "4", "8", "10", "12", "14", "16", "18", "20")
row.names(tinyRNA_piRNAs_DOWN_output_H1) <- rownames(tinyRNA_H1epimutations_piRNAs)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tinyRNA_piRNAs_DOWN_output_H1[,2])
transitions_at_4 <- sum(tinyRNA_piRNAs_DOWN_output_H1[,3])
transitions_at_6 <- NA
transitions_at_8 <- sum(tinyRNA_piRNAs_DOWN_output_H1[,4])
transitions_at_10 <- sum(tinyRNA_piRNAs_DOWN_output_H1[,5])
transitions_at_12 <- sum(tinyRNA_piRNAs_DOWN_output_H1[,6])
transitions_at_14 <- sum(tinyRNA_piRNAs_DOWN_output_H1[,7])
transitions_at_16 <- sum(tinyRNA_piRNAs_DOWN_output_H1[,8])
transitions_at_18 <- sum(tinyRNA_piRNAs_DOWN_output_H1[,9])
transitions_at_20 <- sum(tinyRNA_piRNAs_DOWN_output_H1[,10])

tinyRNA_piRNAs_DOWN_H1_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                          transitions_at_4, 
                                                          transitions_at_6,
                                                          transitions_at_8, 
                                                          transitions_at_10, 
                                                          transitions_at_12, 
                                                          transitions_at_14, 
                                                          transitions_at_16, 
                                                          transitions_at_18, 
                                                          transitions_at_20)

#High2 new epi
tinyRNA_H2epimutations_piRNAs<-subset(tinyRNA_H2epimutations,tinyRNA_H2epimutations$Feature.Class=="piRNA")
#UP

tinyRNA_piRNAs_UP_output_H2<- c()
tinyRNA_ID <- paste(tinyRNA_H2epimutations_piRNAs$Feature.ID, tinyRNA_H2epimutations_piRNAs$Tag,tinyRNA_H2epimutations_piRNAs$Feature.Name,tinyRNA_H2epimutations_piRNAs$Feature.Class, sep="_")
row.names(tinyRNA_H2epimutations_piRNAs)<-tinyRNA_ID
tinyRNA_H2epimutations_piRNAs<-tinyRNA_H2epimutations_piRNAs[,c(7:17)]
for(i in 1:nrow(tinyRNA_H2epimutations_piRNAs)){
  tinyRNA_piRNAs_UP_output_H2 <-rbind(tinyRNA_piRNAs_UP_output_H2, UP_transition_func(tinyRNA_H2epimutations_piRNAs[i,], input_name=row.names(tinyRNA_H2epimutations_piRNAs)[i]))}
colnames(tinyRNA_piRNAs_UP_output_H2) <- c("0", "2", "4", "6", "8","10", "12", "14", "16", "18", "20")
row.names(tinyRNA_piRNAs_UP_output_H2) <- rownames(tinyRNA_H2epimutations_piRNAs)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tinyRNA_piRNAs_UP_output_H2[,2])
transitions_at_4 <- sum(tinyRNA_piRNAs_UP_output_H2[,3])
transitions_at_6 <- sum(tinyRNA_piRNAs_UP_output_H2[,4])
transitions_at_8 <- sum(tinyRNA_piRNAs_UP_output_H2[,5])
transitions_at_10 <- sum(tinyRNA_piRNAs_UP_output_H2[,6])
transitions_at_12 <- sum(tinyRNA_piRNAs_UP_output_H2[,7])
transitions_at_14 <- sum(tinyRNA_piRNAs_UP_output_H2[,8])
transitions_at_16 <- sum(tinyRNA_piRNAs_UP_output_H2[,9])
transitions_at_18 <- sum(tinyRNA_piRNAs_UP_output_H2[,10])
transitions_at_20 <- sum(tinyRNA_piRNAs_UP_output_H2[,11])

tinyRNA_piRNAs_UP_H2_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                        transitions_at_4, 
                                                        transitions_at_6,
                                                        transitions_at_8, 
                                                        transitions_at_10, 
                                                        transitions_at_12, 
                                                        transitions_at_14, 
                                                        transitions_at_16, 
                                                        transitions_at_18, 
                                                        transitions_at_20)


#DOWN
tinyRNA_piRNAs_DOWN_output_H2<- c()
for(i in 1:nrow(tinyRNA_H2epimutations_piRNAs)){
  tinyRNA_piRNAs_DOWN_output_H2 <-rbind(tinyRNA_piRNAs_DOWN_output_H2, DOWN_transition_func(tinyRNA_H2epimutations_piRNAs[i,], input_name=row.names(tinyRNA_H2epimutations_piRNAs)[i]))}
colnames(tinyRNA_piRNAs_DOWN_output_H2) <- c("0", "2", "4", "6", "8", "10", "12", "14", "16", "18", "20")
row.names(tinyRNA_piRNAs_DOWN_output_H2) <- rownames(tinyRNA_H2epimutations_piRNAs)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tinyRNA_piRNAs_DOWN_output_H2[,2])
transitions_at_4 <- sum(tinyRNA_piRNAs_DOWN_output_H2[,3])
transitions_at_6 <- sum(tinyRNA_piRNAs_DOWN_output_H2[,4])
transitions_at_8 <- sum(tinyRNA_piRNAs_DOWN_output_H2[,5])
transitions_at_10 <- sum(tinyRNA_piRNAs_DOWN_output_H2[,6])
transitions_at_12 <- sum(tinyRNA_piRNAs_DOWN_output_H2[,7])
transitions_at_14 <- sum(tinyRNA_piRNAs_DOWN_output_H2[,8])
transitions_at_16 <- sum(tinyRNA_piRNAs_DOWN_output_H2[,9])
transitions_at_18 <- sum(tinyRNA_piRNAs_DOWN_output_H2[,10])
transitions_at_20 <- sum(tinyRNA_piRNAs_DOWN_output_H2[,11])

tinyRNA_piRNAs_DOWN_H2_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                          transitions_at_4, 
                                                          transitions_at_6,
                                                          transitions_at_8, 
                                                          transitions_at_10, 
                                                          transitions_at_12, 
                                                          transitions_at_14, 
                                                          transitions_at_16, 
                                                          transitions_at_18, 
                                                          transitions_at_20)

#Fig.4.C construction
names <- rownames(tinyRNA_piRNAs_UP_output_C1_Table_of_new_epimutations)
rownames(tinyRNA_piRNAs_UP_output_C1_Table_of_new_epimutations) <- NULL
tinyRNA_piRNAs_UP_output_C1_Table_of_new_epimutations <- cbind(names,tinyRNA_piRNAs_UP_output_C1_Table_of_new_epimutations)
colnames(tinyRNA_piRNAs_UP_output_C1_Table_of_new_epimutations) <- c("Transitions","Up")
names <- rownames(tinyRNA_piRNAs_DOWN_output_C1_Table_of_new_epimutations)
rownames(tinyRNA_piRNAs_DOWN_output_C1_Table_of_new_epimutations) <- NULL
tinyRNA_piRNAs_DOWN_output_C1_Table_of_new_epimutations <- cbind(names,tinyRNA_piRNAs_DOWN_output_C1_Table_of_new_epimutations)
colnames(tinyRNA_piRNAs_DOWN_output_C1_Table_of_new_epimutations) <- c("Transitions","Down")

tinyRNA_piRNAs_new_epimut_C1 <- merge(tinyRNA_piRNAs_UP_output_C1_Table_of_new_epimutations,tinyRNA_piRNAs_DOWN_output_C1_Table_of_new_epimutations)
tinyRNA_piRNAs_new_epimut_C1 <- tinyRNA_piRNAs_new_epimut_C1 %>%
  # Creating an empty column:
  add_column(Lineage = "C1", .after="Transitions")

names <- rownames(tinyRNA_piRNAs_UP_output_C2_Table_of_new_epimutations)
rownames(tinyRNA_piRNAs_UP_output_C2_Table_of_new_epimutations) <- NULL
tinyRNA_piRNAs_UP_output_C2_Table_of_new_epimutations <- cbind(names,tinyRNA_piRNAs_UP_output_C2_Table_of_new_epimutations)
colnames(tinyRNA_piRNAs_UP_output_C2_Table_of_new_epimutations) <- c("Transitions","Up")
names <- rownames(tinyRNA_piRNAs_DOWN_output_C2_Table_of_new_epimutations)
rownames(tinyRNA_piRNAs_DOWN_output_C2_Table_of_new_epimutations) <- NULL
tinyRNA_piRNAs_DOWN_output_C2_Table_of_new_epimutations <- cbind(names,tinyRNA_piRNAs_DOWN_output_C2_Table_of_new_epimutations)
colnames(tinyRNA_piRNAs_DOWN_output_C2_Table_of_new_epimutations) <- c("Transitions","Down")

tinyRNA_piRNAs_new_epimut_C2 <- merge(tinyRNA_piRNAs_UP_output_C2_Table_of_new_epimutations,tinyRNA_piRNAs_DOWN_output_C2_Table_of_new_epimutations)
tinyRNA_piRNAs_new_epimut_C2 <- tinyRNA_piRNAs_new_epimut_C2 %>%
  # Creating an empty column:
  add_column(Lineage = "C2", .after="Transitions")

tiny_piRNAs_C <- rbind(tinyRNA_piRNAs_new_epimut_C1, tinyRNA_piRNAs_new_epimut_C2)
tiny_piRNAs_C <- tiny_piRNAs_C %>%
  # Creating an empty column:
  add_column(Condition = "Control", .after="Transitions")

names <- rownames(tinyRNA_piRNAs_UP_L1_Table_of_new_epimutations)
rownames(tinyRNA_piRNAs_UP_L1_Table_of_new_epimutations) <- NULL
tinyRNA_piRNAs_UP_L1_Table_of_new_epimutations <- cbind(names,tinyRNA_piRNAs_UP_L1_Table_of_new_epimutations)
colnames(tinyRNA_piRNAs_UP_L1_Table_of_new_epimutations) <- c("Transitions","Up")
names <- rownames(tinyRNA_piRNAs_DOWN_L1_Table_of_new_epimutations)
rownames(tinyRNA_piRNAs_DOWN_L1_Table_of_new_epimutations) <- NULL
tinyRNA_piRNAs_DOWN_L1_Table_of_new_epimutations <- cbind(names,tinyRNA_piRNAs_DOWN_L1_Table_of_new_epimutations)
colnames(tinyRNA_piRNAs_DOWN_L1_Table_of_new_epimutations) <- c("Transitions","Down")

tinyRNA_piRNAs_new_epimut_L1 <- merge(tinyRNA_piRNAs_UP_L1_Table_of_new_epimutations,tinyRNA_piRNAs_DOWN_L1_Table_of_new_epimutations)
tinyRNA_piRNAs_new_epimut_L1 <- tinyRNA_piRNAs_new_epimut_L1 %>%
  # Creating an empty column:
  add_column(Lineage = "L1", .after="Transitions")
names <- rownames(tinyRNA_piRNAs_UP_L2_Table_of_new_epimutations)
rownames(tinyRNA_piRNAs_UP_L2_Table_of_new_epimutations) <- NULL
tinyRNA_piRNAs_UP_L2_Table_of_new_epimutations <- cbind(names,tinyRNA_piRNAs_UP_L2_Table_of_new_epimutations)
colnames(tinyRNA_piRNAs_UP_L2_Table_of_new_epimutations) <- c("Transitions","Up")
names <- rownames(tinyRNA_piRNAs_DOWN_L2_Table_of_new_epimutations)
rownames(tinyRNA_piRNAs_DOWN_L2_Table_of_new_epimutations) <- NULL
tinyRNA_piRNAs_DOWN_L2_Table_of_new_epimutations <- cbind(names,tinyRNA_piRNAs_DOWN_L2_Table_of_new_epimutations)
colnames(tinyRNA_piRNAs_DOWN_L2_Table_of_new_epimutations) <- c("Transitions","Down")

tinyRNA_piRNAs_new_epimut_L2 <- merge(tinyRNA_piRNAs_UP_L2_Table_of_new_epimutations,tinyRNA_piRNAs_DOWN_L2_Table_of_new_epimutations)
tinyRNA_piRNAs_new_epimut_L2 <- tinyRNA_piRNAs_new_epimut_L2 %>%
  # Creating an empty column:
  add_column(Lineage = "L2", .after="Transitions")

tiny_piRNAs_L <- rbind(tinyRNA_piRNAs_new_epimut_L1, tinyRNA_piRNAs_new_epimut_L2)
tiny_piRNAs_L <- tiny_piRNAs_L %>%
  # Creating an empty column:
  add_column(Condition = "Low dose", .after="Transitions")
names <- rownames(tinyRNA_piRNAs_UP_H1_Table_of_new_epimutations)
rownames(tinyRNA_piRNAs_UP_H1_Table_of_new_epimutations) <- NULL
tinyRNA_piRNAs_UP_H1_Table_of_new_epimutations <- cbind(names,tinyRNA_piRNAs_UP_H1_Table_of_new_epimutations)
colnames(tinyRNA_piRNAs_UP_H1_Table_of_new_epimutations) <- c("Transitions","Up")
names <- rownames(tinyRNA_piRNAs_DOWN_H1_Table_of_new_epimutations)
rownames(tinyRNA_piRNAs_DOWN_H1_Table_of_new_epimutations) <- NULL
tinyRNA_piRNAs_DOWN_H1_Table_of_new_epimutations <- cbind(names,tinyRNA_piRNAs_DOWN_H1_Table_of_new_epimutations)
colnames(tinyRNA_piRNAs_DOWN_H1_Table_of_new_epimutations) <- c("Transitions","Down")

tinyRNA_piRNAs_new_epimut_H1 <- merge(tinyRNA_piRNAs_UP_H1_Table_of_new_epimutations,tinyRNA_piRNAs_DOWN_H1_Table_of_new_epimutations)
tinyRNA_piRNAs_new_epimut_H1 <- tinyRNA_piRNAs_new_epimut_H1 %>%
  # Creating an empty column:
  add_column(Lineage = "H1", .after="Transitions")
names <- rownames(tinyRNA_piRNAs_UP_H2_Table_of_new_epimutations)
rownames(tinyRNA_piRNAs_UP_H2_Table_of_new_epimutations) <- NULL
tinyRNA_piRNAs_UP_H2_Table_of_new_epimutations <- cbind(names,tinyRNA_piRNAs_UP_H2_Table_of_new_epimutations)
colnames(tinyRNA_piRNAs_UP_H2_Table_of_new_epimutations) <- c("Transitions","Up")
names <- rownames(tinyRNA_piRNAs_DOWN_H2_Table_of_new_epimutations)
rownames(tinyRNA_piRNAs_DOWN_H2_Table_of_new_epimutations) <- NULL
tinyRNA_piRNAs_DOWN_H2_Table_of_new_epimutations <- cbind(names,tinyRNA_piRNAs_DOWN_H2_Table_of_new_epimutations)
colnames(tinyRNA_piRNAs_DOWN_H2_Table_of_new_epimutations) <- c("Transitions","Down")

tinyRNA_piRNAs_new_epimut_H2 <- merge(tinyRNA_piRNAs_UP_H2_Table_of_new_epimutations,tinyRNA_piRNAs_DOWN_H2_Table_of_new_epimutations)
tinyRNA_piRNAs_new_epimut_H2 <- tinyRNA_piRNAs_new_epimut_H2 %>%
  # Creating an empty column:
  add_column(Lineage = "H2", .after="Transitions")

tiny_piRNAs_H <- rbind(tinyRNA_piRNAs_new_epimut_H1, tinyRNA_piRNAs_new_epimut_H2)
tiny_piRNAs_H <- tiny_piRNAs_H %>%
  # Creating an empty column:
  add_column(Condition = "High dose", .after="Transitions")

tiny_piRNAs_newepi<-rbind(tiny_piRNAs_C,tiny_piRNAs_L,tiny_piRNAs_H)
tiny_piRNAs_newepi$Up <- as.numeric(tiny_piRNAs_newepi$Up)
tiny_piRNAs_newepi$Down <- as.numeric(tiny_piRNAs_newepi$Down)
tiny_piRNAs_newepi$Total=rowSums(cbind(tiny_piRNAs_newepi$Up,tiny_piRNAs_newepi$Down),na.rm=FALSE)

UP<-tiny_piRNAs_newepi$Up
DOWN<-tiny_piRNAs_newepi$Down

#Test for statistical significance
ggdensity(tiny_piRNAs_newepi$Total, 
          main = "Density plot of epimutations_number",
          xlab = "Total number of epimutations_number")
shapiro.test(tiny_piRNAs_newepi$Total)
kruskal.test(Total ~ Condition, data = tiny_piRNAs_newepi)
dunnTest(Total ~ Condition, data = tiny_piRNAs_newepi)

#Fig.4.C
tiny_piRNAs_newepi$Condition <- fct_relevel(tiny_piRNAs_newepi$Condition, c("Control", "Low dose","High dose"))
tiny_piRNAs_all_rate <- ggplot(tiny_piRNAs_newepi, aes(x=Condition, y=Total, color = Condition,group=Condition)) + 
  geom_boxplot(fatten = 1, lwd = 1, width=0.5)+
  scale_color_manual(values=c("cornflowerblue", "darkgreen", "red"))+
  labs(y = "Number of new epimutations", x = "\nCondition")+
  geom_dotplot(binaxis='y', binwidth = 0.25, stackdir='center',  dotsize=195, binpositions = "all", stackgroups = TRUE)+
  theme_bw()+
  theme_linedraw()+
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(color = "black", size = 25, angle = 0, face = "plain"),
        axis.text.y = element_text(color = "black", size = 25, angle = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 20, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, face = "bold"))+
  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))+  
  scale_fill_brewer(palette="Pastel2")+
  ggtitle(paste("tiny_piRNAs_allLineage_rate"))
tiny_piRNAs_all_rate

#Sup.Fig.5.C
Generation<-c("10","12","14","16","18","2","20","4","6","8",
              "10","12","14","16","18","2","20","4","6","8",
              "10","12","14","16","18","2","20","4","6","8",
              "10","12","14","16","18","2","20","4","6","8",
              "10","12","14","16","18","2","20","4","6","8",
              "10","12","14","16","18","2","20","4","6","8")
tiny_piRNAs_newepi_detailed<-cbind(tiny_piRNAs_newepi,Generation)
tiny_piRNAs_newepi_detailed$Generation <- fct_relevel(tiny_piRNAs_newepi_detailed$Generation, c("2", "4", "6","8","10","12","14","16","18","20"))
tiny_piRNAs_newepi_detailed$Condition <- fct_relevel(tiny_piRNAs_newepi_detailed$Condition, c("Control", "Low dose", "High dose"))
tiny_piRNAs_newepi_detailed_plot <- ggplot(tiny_piRNAs_newepi_detailed, aes(x=Generation , y=Total, color = Condition)) + 
  geom_boxplot(fatten = 1, lwd = 1, width=0.5)+
  scale_color_manual(values=c("cornflowerblue", "darkgreen", "red"))+
  labs(y = "Number of new epimutations", x = "\nGeneration")+
  geom_jitter(aes(shape=Lineage, size=8),
              position=position_jitter(width = 0.2,
                                       height = 0.2))+
  theme_bw()+
  theme(axis.text.x = element_text(color = "black", size = 20, angle = 0, face = "plain"),
        axis.text.y = element_text(color = "black", size = 20, angle = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 16, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 16, angle = 90, face = "bold"))+
  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))+  
  scale_fill_brewer(palette="Pastel2")+
  ggtitle(paste(""))
tiny_piRNAs_newepi_detailed_plot 

#Get raw data
write.xlsx(tiny_piRNAs_newepi_detailed,"Data_Fig_4_C_&_Sup_Fig_5_C.xlsx") 
#--------------------------
#####Fig.4.D & Sup.Fig.6.D - piRNAs epimutations duration#####
tiny_piRNAs_epimut<-tiny_all_epimut %>% filter(grepl('piRNA', gene_name))
#Fig.4.D
CoxMod<-coxph(Surv(length,complete)~Condition,data=tiny_piRNAs_epimut)
ggforest(CoxMod, data=tiny_piRNAs_epimut)

ggsurvplot(survfit(CoxMod,data=tiny_piRNAs_epimut), palette = "#2E9FDF",cof.int=TRUE,
           ggtheme = theme_minimal())
condition_df <- with(tiny_piRNAs_epimut,
                     data.frame(Condition = c("Control","Low dose","High dose")
                     )
)
condition_df

fit <- survfit(CoxMod, newdata = condition_df)
fit <- survfit(Surv(length, complete) ~ Condition, data = tiny_piRNAs_epimut)
ggsurvplot <- ggsurvplot(fit,  tiny_piRNAs_epimut, censor = T, break.time.by= 1, pval=TRUE,pval.coord = c(0, 0.03),
                         font.main = c(16, "bold"),
                         palette = c("cornflowerblue", "darkgreen","red"),
                         font.x = c(20,"bold"),
                         font.y = c(20,"bold"),
                         font.legend = 20,
                         legend.labs=c("Control","Low dose","High dose"),
                         font.tickslab = 18)+
  # surv.median.line = "hv")+
  ggtitle("Changes survival all conditions")+
  xlab("Time (generations)")
ggsurvplot$plot +theme(plot.title = element_text(hjust = 0.5))

#Sup.Fig.5.D
CoxMod<-coxph(Surv(length,complete)~Lineage,data=tiny_piRNAs_epimut)
ggforest(model=CoxMod,data=tiny_piRNAs_epimut,fontsize = 0.8, noDigits = 2)
fit_2 <- survfit(Surv(length, complete) ~ Lineage, data = tiny_piRNAs_epimut)
ggsurvplot <- ggsurvplot(fit_2,  tiny_piRNAs_epimut, censor = T, break.time.by= 1, pval=TRUE,pval.coord = c(0, 0.03),                         font.main = c(16, "bold"),
                         font.x = c(20,"bold"),
                         palette = c("#3399FF","blue","red","2E9FDF", "green","#006600"),
                         font.y = c(20,"bold"),
                         font.legend = 20,
                         legend.labs=c("C1","C2","H1","H2","L1","L2"),
                         font.tickslab = 18)+
  ggtitle("Changes survival all conditions")+
  xlab("Time (generations)")
ggsurvplot$plot +theme(plot.title = element_text(hjust = 0.5)) 

#Get raw data
write.xlsx(tiny_piRNAs_epimut,"Data_Fig_4_D_&_Sup_Fig_5_D.xlsx")
#--------------------------
#####Fig.4.E & Sup.Fig.6.E - miRNAs epimutations####
#Control1 new epi
tinyRNA_C1epimutations_miRNAs<-subset(tinyRNA_C1epimutations,tinyRNA_C1epimutations$Feature.Class=="miRNA")
#UP
tinyRNA_miRNAs_UP_output_C1<- c()
tinyRNA_ID <- paste(tinyRNA_C1epimutations_miRNAs$Feature.ID, tinyRNA_C1epimutations_miRNAs$Tag,tinyRNA_C1epimutations_miRNAs$Feature.Name,tinyRNA_C1epimutations_miRNAs$Feature.Class, sep="_")
row.names(tinyRNA_C1epimutations_miRNAs)<-tinyRNA_ID
tinyRNA_C1epimutations_miRNAs<-tinyRNA_C1epimutations_miRNAs[,c(7:17)]
for(i in 1:nrow(tinyRNA_C1epimutations_miRNAs)){
  tinyRNA_miRNAs_UP_output_C1 <-rbind(tinyRNA_miRNAs_UP_output_C1, UP_transition_func(tinyRNA_C1epimutations_miRNAs[i,], input_name=row.names(tinyRNA_C1epimutations_miRNAs)[i]))}
colnames(tinyRNA_miRNAs_UP_output_C1) <- c("0", "2", "4", "6", "8", "10", "12", "14", "16", "18", "20")
row.names(tinyRNA_miRNAs_UP_output_C1) <- rownames(tinyRNA_C1epimutations_miRNAs)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tinyRNA_miRNAs_UP_output_C1[,2])
transitions_at_4 <- sum(tinyRNA_miRNAs_UP_output_C1[,3])
transitions_at_6 <- sum(tinyRNA_miRNAs_UP_output_C1[,4])
transitions_at_8 <- sum(tinyRNA_miRNAs_UP_output_C1[,5])
transitions_at_10 <- sum(tinyRNA_miRNAs_UP_output_C1[,6])
transitions_at_12 <- sum(tinyRNA_miRNAs_UP_output_C1[,7])
transitions_at_14 <- sum(tinyRNA_miRNAs_UP_output_C1[,8])
transitions_at_16 <- sum(tinyRNA_miRNAs_UP_output_C1[,9])
transitions_at_18 <- sum(tinyRNA_miRNAs_UP_output_C1[,10])
transitions_at_20 <- sum(tinyRNA_miRNAs_UP_output_C1[,11])

tinyRNA_miRNAs_UP_output_C1_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                               transitions_at_4, 
                                                               transitions_at_6,
                                                               transitions_at_8, 
                                                               transitions_at_10, 
                                                               transitions_at_12, 
                                                               transitions_at_14, 
                                                               transitions_at_16, 
                                                               transitions_at_18, 
                                                               transitions_at_20)

#DOWN
tinyRNA_miRNAs_DOWN_output_C1<- c()
for(i in 1:nrow(tinyRNA_C1epimutations_miRNAs)){
  tinyRNA_miRNAs_DOWN_output_C1 <-rbind(tinyRNA_miRNAs_DOWN_output_C1, DOWN_transition_func(tinyRNA_C1epimutations_miRNAs[i,], input_name=row.names(tinyRNA_C1epimutations_miRNAs)[i]))}
colnames(tinyRNA_miRNAs_DOWN_output_C1) <- c("0", "2", "4", "6", "8", "10", "12", "14", "16", "18", "20")
row.names(tinyRNA_miRNAs_DOWN_output_C1) <- rownames(tinyRNA_C1epimutations_miRNAs)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tinyRNA_miRNAs_DOWN_output_C1[,2])
transitions_at_4 <- sum(tinyRNA_miRNAs_DOWN_output_C1[,3])
transitions_at_6 <- sum(tinyRNA_miRNAs_DOWN_output_C1[,4])
transitions_at_8 <- sum(tinyRNA_miRNAs_DOWN_output_C1[,5])
transitions_at_10 <- sum(tinyRNA_miRNAs_DOWN_output_C1[,6])
transitions_at_12 <- sum(tinyRNA_miRNAs_DOWN_output_C1[,7])
transitions_at_14 <- sum(tinyRNA_miRNAs_DOWN_output_C1[,8])
transitions_at_16 <- sum(tinyRNA_miRNAs_DOWN_output_C1[,9])
transitions_at_18 <- sum(tinyRNA_miRNAs_DOWN_output_C1[,10])
transitions_at_20 <- sum(tinyRNA_miRNAs_DOWN_output_C1[,11])

tinyRNA_miRNAs_DOWN_output_C1_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                                 transitions_at_4, 
                                                                 transitions_at_6,
                                                                 transitions_at_8, 
                                                                 transitions_at_10, 
                                                                 transitions_at_12, 
                                                                 transitions_at_14, 
                                                                 transitions_at_16, 
                                                                 transitions_at_18, 
                                                                 transitions_at_20)

#Control2 new epi
tinyRNA_C2epimutations_miRNAs<-subset(tinyRNA_C2epimutations,tinyRNA_C2epimutations$Feature.Class=="miRNA")
#UP
tinyRNA_miRNAs_UP_output_C2<- c()
tinyRNA_ID <- paste(tinyRNA_C2epimutations_miRNAs$Feature.ID, tinyRNA_C2epimutations_miRNAs$Tag,tinyRNA_C2epimutations_miRNAs$Feature.Name,tinyRNA_C2epimutations_miRNAs$Feature.Class, sep="_")
row.names(tinyRNA_C2epimutations_miRNAs)<-tinyRNA_ID
tinyRNA_C2epimutations_miRNAs<-tinyRNA_C2epimutations_miRNAs[,c(7:16)]
for(i in 1:nrow(tinyRNA_C2epimutations_miRNAs)){
  tinyRNA_miRNAs_UP_output_C2 <-rbind(tinyRNA_miRNAs_UP_output_C2, UP_transition_func(tinyRNA_C2epimutations_miRNAs[i,], input_name=row.names(tinyRNA_C2epimutations_miRNAs)[i]))}
colnames(tinyRNA_miRNAs_UP_output_C2) <- c("0", "2", "4","8", "10", "12", "14", "16", "18", "20")
row.names(tinyRNA_miRNAs_UP_output_C2) <- rownames(tinyRNA_C2epimutations_miRNAs)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tinyRNA_miRNAs_UP_output_C2[,2])
transitions_at_4 <- sum(tinyRNA_miRNAs_UP_output_C2[,3])
transitions_at_6 <- NA
transitions_at_8 <- sum(tinyRNA_miRNAs_UP_output_C2[,4])
transitions_at_10 <- sum(tinyRNA_miRNAs_UP_output_C2[,5])
transitions_at_12 <- sum(tinyRNA_miRNAs_UP_output_C2[,6])
transitions_at_14 <- sum(tinyRNA_miRNAs_UP_output_C2[,7])
transitions_at_16 <- sum(tinyRNA_miRNAs_UP_output_C2[,8])
transitions_at_18 <- sum(tinyRNA_miRNAs_UP_output_C2[,9])
transitions_at_20 <- sum(tinyRNA_miRNAs_UP_output_C2[,10])

tinyRNA_miRNAs_UP_output_C2_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                               transitions_at_4, 
                                                               transitions_at_6,
                                                               transitions_at_8, 
                                                               transitions_at_10, 
                                                               transitions_at_12, 
                                                               transitions_at_14, 
                                                               transitions_at_16, 
                                                               transitions_at_18, 
                                                               transitions_at_20)

#DOWN
tinyRNA_miRNAs_DOWN_output_C2<- c()
for(i in 1:nrow(tinyRNA_C2epimutations_miRNAs)){
  tinyRNA_miRNAs_DOWN_output_C2 <-rbind(tinyRNA_miRNAs_DOWN_output_C2, DOWN_transition_func(tinyRNA_C2epimutations_miRNAs[i,], input_name=row.names(tinyRNA_C2epimutations_miRNAs)[i]))}
colnames(tinyRNA_miRNAs_DOWN_output_C2) <- c("0", "2", "4", "8", "10", "12", "14", "16", "18", "20")
row.names(tinyRNA_miRNAs_DOWN_output_C2) <- rownames(tinyRNA_C2epimutations_miRNAs)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tinyRNA_miRNAs_DOWN_output_C2[,2])
transitions_at_4 <- sum(tinyRNA_miRNAs_DOWN_output_C2[,3])
transitions_at_6 <- NA
transitions_at_8 <- sum(tinyRNA_miRNAs_DOWN_output_C2[,4])
transitions_at_10 <- sum(tinyRNA_miRNAs_DOWN_output_C2[5])
transitions_at_12 <- sum(tinyRNA_miRNAs_DOWN_output_C2[,6])
transitions_at_14 <- sum(tinyRNA_miRNAs_DOWN_output_C2[,7])
transitions_at_16 <- sum(tinyRNA_miRNAs_DOWN_output_C2[,8])
transitions_at_18 <- sum(tinyRNA_miRNAs_DOWN_output_C2[,9])
transitions_at_20 <- sum(tinyRNA_miRNAs_DOWN_output_C2[,10])
#transitions_at_10 <- sum(tinyRNA_miRNAs_DOWN_output_C2[5]) =0
transitions_at_10 <- 0

tinyRNA_miRNAs_DOWN_output_C2_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                                 transitions_at_4, 
                                                                 transitions_at_6,
                                                                 transitions_at_8, 
                                                                 transitions_at_10, 
                                                                 transitions_at_12, 
                                                                 transitions_at_14, 
                                                                 transitions_at_16, 
                                                                 transitions_at_18, 
                                                                 transitions_at_20)

#Low1 new epi
tinyRNA_L1epimutations_miRNAs<-subset(tinyRNA_L1epimutations,tinyRNA_L1epimutations$Feature.Class=="miRNA")
#UP

tinyRNA_miRNAs_UP_output_L1<- c()
tinyRNA_ID <- paste(tinyRNA_L1epimutations_miRNAs$Feature.ID, tinyRNA_L1epimutations_miRNAs$Tag,tinyRNA_L1epimutations_miRNAs$Feature.Name,tinyRNA_L1epimutations_miRNAs$Feature.Class, sep="_")
row.names(tinyRNA_L1epimutations_miRNAs)<-tinyRNA_ID
tinyRNA_L1epimutations_miRNAs<-tinyRNA_L1epimutations_miRNAs[,c(7:16)]
for(i in 1:nrow(tinyRNA_L1epimutations_miRNAs)){
  tinyRNA_miRNAs_UP_output_L1 <-rbind(tinyRNA_miRNAs_UP_output_L1, UP_transition_func(tinyRNA_L1epimutations_miRNAs[i,], input_name=row.names(tinyRNA_L1epimutations_miRNAs)[i]))}
colnames(tinyRNA_miRNAs_UP_output_L1) <- c("0", "2", "4","6", "8", "12", "14", "16", "18", "20")
row.names(tinyRNA_miRNAs_UP_output_L1) <- rownames(tinyRNA_L1epimutations_miRNAs)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tinyRNA_miRNAs_UP_output_L1[,2])
transitions_at_4 <- sum(tinyRNA_miRNAs_UP_output_L1[,3])
transitions_at_6 <- sum(tinyRNA_miRNAs_UP_output_L1[,4])
transitions_at_8 <- sum(tinyRNA_miRNAs_UP_output_L1[,5])
transitions_at_10 <- NA
transitions_at_12 <- sum(tinyRNA_miRNAs_UP_output_L1[,6])
transitions_at_14 <- sum(tinyRNA_miRNAs_UP_output_L1[,7])
transitions_at_16 <- sum(tinyRNA_miRNAs_UP_output_L1[,8])
transitions_at_18 <- sum(tinyRNA_miRNAs_UP_output_L1[,9])
transitions_at_20 <- sum(tinyRNA_miRNAs_UP_output_L1[,10])

tinyRNA_miRNAs_UP_L1_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                        transitions_at_4, 
                                                        transitions_at_6,
                                                        transitions_at_8, 
                                                        transitions_at_10, 
                                                        transitions_at_12, 
                                                        transitions_at_14, 
                                                        transitions_at_16, 
                                                        transitions_at_18, 
                                                        transitions_at_20)

#DOWN
tinyRNA_miRNAs_DOWN_output_L1<- c()
for(i in 1:nrow(tinyRNA_L1epimutations_miRNAs)){
  tinyRNA_miRNAs_DOWN_output_L1 <-rbind(tinyRNA_miRNAs_DOWN_output_L1, DOWN_transition_func(tinyRNA_L1epimutations_miRNAs[i,], input_name=row.names(tinyRNA_L1epimutations_miRNAs)[i]))}
colnames(tinyRNA_miRNAs_DOWN_output_L1) <- c("0", "2", "4", "6","8", "12", "14", "16", "18", "20")
row.names(tinyRNA_miRNAs_DOWN_output_L1) <- rownames(tinyRNA_L1epimutations_miRNAs)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tinyRNA_miRNAs_DOWN_output_L1[,2])
transitions_at_4 <- sum(tinyRNA_miRNAs_DOWN_output_L1[,3])
transitions_at_6 <- sum(tinyRNA_miRNAs_DOWN_output_L1[,4])
transitions_at_8 <- sum(tinyRNA_miRNAs_DOWN_output_L1[,5])
transitions_at_10 <- NA
transitions_at_12 <- sum(tinyRNA_miRNAs_DOWN_output_L1[,6])
transitions_at_14 <- sum(tinyRNA_miRNAs_DOWN_output_L1[,7])
transitions_at_16 <- sum(tinyRNA_miRNAs_DOWN_output_L1[,8])
transitions_at_18 <- sum(tinyRNA_miRNAs_DOWN_output_L1[,9])
transitions_at_20 <- sum(tinyRNA_miRNAs_DOWN_output_L1[,10])

tinyRNA_miRNAs_DOWN_L1_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                          transitions_at_4, 
                                                          transitions_at_6,
                                                          transitions_at_8, 
                                                          transitions_at_10, 
                                                          transitions_at_12, 
                                                          transitions_at_14, 
                                                          transitions_at_16, 
                                                          transitions_at_18, 
                                                          transitions_at_20)

#Low2 new epi
tinyRNA_L2epimutations_miRNAs<-subset(tinyRNA_L2epimutations,tinyRNA_L2epimutations$Feature.Class=="miRNA")
#UP

tinyRNA_miRNAs_UP_output_L2<- c()
tinyRNA_ID <- paste(tinyRNA_L2epimutations_miRNAs$Feature.ID, tinyRNA_L2epimutations_miRNAs$Tag,tinyRNA_L2epimutations_miRNAs$Feature.Name,tinyRNA_L2epimutations_miRNAs$Feature.Class, sep="_")
row.names(tinyRNA_L2epimutations_miRNAs)<-tinyRNA_ID
tinyRNA_L2epimutations_miRNAs<-tinyRNA_L2epimutations_miRNAs[,c(7:16)]
for(i in 1:nrow(tinyRNA_L2epimutations_miRNAs)){
  tinyRNA_miRNAs_UP_output_L2 <-rbind(tinyRNA_miRNAs_UP_output_L2, UP_transition_func(tinyRNA_L2epimutations_miRNAs[i,], input_name=row.names(tinyRNA_L2epimutations_miRNAs)[i]))}
colnames(tinyRNA_miRNAs_UP_output_L2) <- c("0", "2", "4", "8","10", "12", "14", "16", "18", "20")
row.names(tinyRNA_miRNAs_UP_output_L2) <- rownames(tinyRNA_L2epimutations_miRNAs)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tinyRNA_miRNAs_UP_output_L2[,2])
transitions_at_4 <- sum(tinyRNA_miRNAs_UP_output_L2[,3])
transitions_at_6 <- NA
transitions_at_8 <- sum(tinyRNA_miRNAs_UP_output_L2[,4])
transitions_at_10 <- sum(tinyRNA_miRNAs_UP_output_L2[,5])
transitions_at_12 <- sum(tinyRNA_miRNAs_UP_output_L2[,6])
transitions_at_14 <- sum(tinyRNA_miRNAs_UP_output_L2[,7])
transitions_at_16 <- sum(tinyRNA_miRNAs_UP_output_L2[,8])
transitions_at_18 <- sum(tinyRNA_miRNAs_UP_output_L2[,9])
transitions_at_20 <- sum(tinyRNA_miRNAs_UP_output_L2[,10])

tinyRNA_miRNAs_UP_L2_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                        transitions_at_4, 
                                                        transitions_at_6,
                                                        transitions_at_8, 
                                                        transitions_at_10, 
                                                        transitions_at_12, 
                                                        transitions_at_14, 
                                                        transitions_at_16, 
                                                        transitions_at_18, 
                                                        transitions_at_20)


#DOWN
tinyRNA_miRNAs_DOWN_output_L2<- c()
for(i in 1:nrow(tinyRNA_L2epimutations_miRNAs)){
  tinyRNA_miRNAs_DOWN_output_L2 <-rbind(tinyRNA_miRNAs_DOWN_output_L2, DOWN_transition_func(tinyRNA_L2epimutations_miRNAs[i,], input_name=row.names(tinyRNA_L2epimutations_miRNAs)[i]))}
colnames(tinyRNA_miRNAs_DOWN_output_L2) <- c("0", "2", "4", "6","8", "12", "14", "16", "18", "20")
row.names(tinyRNA_miRNAs_DOWN_output_L2) <- rownames(tinyRNA_L2epimutations_miRNAs)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tinyRNA_miRNAs_DOWN_output_L2[,2])
transitions_at_4 <- sum(tinyRNA_miRNAs_DOWN_output_L2[,3])
transitions_at_6 <- NA
transitions_at_8 <- sum(tinyRNA_miRNAs_DOWN_output_L2[,4])
transitions_at_10 <- sum(tinyRNA_miRNAs_DOWN_output_L2[,5])
transitions_at_12 <- sum(tinyRNA_miRNAs_DOWN_output_L2[,6])
transitions_at_14 <- sum(tinyRNA_miRNAs_DOWN_output_L2[,7])
transitions_at_16 <- sum(tinyRNA_miRNAs_DOWN_output_L2[,8])
transitions_at_18 <- sum(tinyRNA_miRNAs_DOWN_output_L2[,9])
transitions_at_20 <- sum(tinyRNA_miRNAs_DOWN_output_L2[,10])

tinyRNA_miRNAs_DOWN_L2_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                          transitions_at_4, 
                                                          transitions_at_6,
                                                          transitions_at_8, 
                                                          transitions_at_10, 
                                                          transitions_at_12, 
                                                          transitions_at_14, 
                                                          transitions_at_16, 
                                                          transitions_at_18, 
                                                          transitions_at_20)

#High1 new epi
tinyRNA_H1epimutations_miRNAs<-subset(tinyRNA_H1epimutations,tinyRNA_H1epimutations$Feature.Class=="miRNA")
#UP

tinyRNA_miRNAs_UP_output_H1<- c()
tinyRNA_ID <- paste(tinyRNA_H1epimutations_miRNAs$Feature.ID, tinyRNA_H1epimutations_miRNAs$Tag,tinyRNA_H1epimutations_miRNAs$Feature.Name,tinyRNA_H1epimutations_miRNAs$Feature.Class, sep="_")
row.names(tinyRNA_H1epimutations_miRNAs)<-tinyRNA_ID
tinyRNA_H1epimutations_miRNAs<-tinyRNA_H1epimutations_miRNAs[,c(7:16)]
for(i in 1:nrow(tinyRNA_H1epimutations_miRNAs)){
  tinyRNA_miRNAs_UP_output_H1 <-rbind(tinyRNA_miRNAs_UP_output_H1, UP_transition_func(tinyRNA_H1epimutations_miRNAs[i,], input_name=row.names(tinyRNA_H1epimutations_miRNAs)[i]))}
colnames(tinyRNA_miRNAs_UP_output_H1) <- c("0", "2", "4", "8","10", "12", "14", "16", "18", "20")
row.names(tinyRNA_miRNAs_UP_output_H1) <- rownames(tinyRNA_H1epimutations_miRNAs)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tinyRNA_miRNAs_UP_output_H1[,2])
transitions_at_4 <- sum(tinyRNA_miRNAs_UP_output_H1[,3])
transitions_at_6 <- NA
transitions_at_8 <- sum(tinyRNA_miRNAs_UP_output_H1[,4])
transitions_at_10 <- sum(tinyRNA_miRNAs_UP_output_H1[,5])
transitions_at_12 <- sum(tinyRNA_miRNAs_UP_output_H1[,6])
transitions_at_14 <- sum(tinyRNA_miRNAs_UP_output_H1[,7])
transitions_at_16 <- sum(tinyRNA_miRNAs_UP_output_H1[,8])
transitions_at_18 <- sum(tinyRNA_miRNAs_UP_output_H1[,9])
transitions_at_20 <- sum(tinyRNA_miRNAs_UP_output_H1[,10])

tinyRNA_miRNAs_UP_H1_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                        transitions_at_4, 
                                                        transitions_at_6,
                                                        transitions_at_8, 
                                                        transitions_at_10, 
                                                        transitions_at_12, 
                                                        transitions_at_14, 
                                                        transitions_at_16, 
                                                        transitions_at_18, 
                                                        transitions_at_20)


#DOWN
tinyRNA_miRNAs_DOWN_output_H1<- c()
for(i in 1:nrow(tinyRNA_H1epimutations_miRNAs)){
  tinyRNA_miRNAs_DOWN_output_H1 <-rbind(tinyRNA_miRNAs_DOWN_output_H1, DOWN_transition_func(tinyRNA_H1epimutations_miRNAs[i,], input_name=row.names(tinyRNA_H1epimutations_miRNAs)[i]))}
colnames(tinyRNA_miRNAs_DOWN_output_H1) <- c("0", "2", "4", "8", "10", "12", "14", "16", "18", "20")
row.names(tinyRNA_miRNAs_DOWN_output_H1) <- rownames(tinyRNA_H1epimutations_miRNAs)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tinyRNA_miRNAs_DOWN_output_H1[,2])
transitions_at_4 <- sum(tinyRNA_miRNAs_DOWN_output_H1[,3])
transitions_at_6 <- NA
transitions_at_8 <- sum(tinyRNA_miRNAs_DOWN_output_H1[,4])
transitions_at_10 <- sum(tinyRNA_miRNAs_DOWN_output_H1[,5])
transitions_at_12 <- sum(tinyRNA_miRNAs_DOWN_output_H1[,6])
transitions_at_14 <- sum(tinyRNA_miRNAs_DOWN_output_H1[,7])
transitions_at_16 <- sum(tinyRNA_miRNAs_DOWN_output_H1[,8])
transitions_at_18 <- sum(tinyRNA_miRNAs_DOWN_output_H1[,9])
transitions_at_20 <- sum(tinyRNA_miRNAs_DOWN_output_H1[,10])

tinyRNA_miRNAs_DOWN_H1_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                          transitions_at_4, 
                                                          transitions_at_6,
                                                          transitions_at_8, 
                                                          transitions_at_10, 
                                                          transitions_at_12, 
                                                          transitions_at_14, 
                                                          transitions_at_16, 
                                                          transitions_at_18, 
                                                          transitions_at_20)

#High2 new epi
tinyRNA_H2epimutations_miRNAs<-subset(tinyRNA_H2epimutations,tinyRNA_H2epimutations$Feature.Class=="miRNA")
#UP

tinyRNA_miRNAs_UP_output_H2<- c()
tinyRNA_ID <- paste(tinyRNA_H2epimutations_miRNAs$Feature.ID, tinyRNA_H2epimutations_miRNAs$Tag,tinyRNA_H2epimutations_miRNAs$Feature.Name,tinyRNA_H2epimutations_miRNAs$Feature.Class, sep="_")
row.names(tinyRNA_H2epimutations_miRNAs)<-tinyRNA_ID
tinyRNA_H2epimutations_miRNAs<-tinyRNA_H2epimutations_miRNAs[,c(7:17)]
for(i in 1:nrow(tinyRNA_H2epimutations_miRNAs)){
  tinyRNA_miRNAs_UP_output_H2 <-rbind(tinyRNA_miRNAs_UP_output_H2, UP_transition_func(tinyRNA_H2epimutations_miRNAs[i,], input_name=row.names(tinyRNA_H2epimutations_miRNAs)[i]))}
colnames(tinyRNA_miRNAs_UP_output_H2) <- c("0", "2", "4", "6", "8","10", "12", "14", "16", "18", "20")
row.names(tinyRNA_miRNAs_UP_output_H2) <- rownames(tinyRNA_H2epimutations_miRNAs)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tinyRNA_miRNAs_UP_output_H2[,2])
transitions_at_4 <- sum(tinyRNA_miRNAs_UP_output_H2[,3])
transitions_at_6 <- sum(tinyRNA_miRNAs_UP_output_H2[,4])
transitions_at_8 <- sum(tinyRNA_miRNAs_UP_output_H2[,5])
transitions_at_10 <- sum(tinyRNA_miRNAs_UP_output_H2[,6])
transitions_at_12 <- sum(tinyRNA_miRNAs_UP_output_H2[,7])
transitions_at_14 <- sum(tinyRNA_miRNAs_UP_output_H2[,8])
transitions_at_16 <- sum(tinyRNA_miRNAs_UP_output_H2[,9])
transitions_at_18 <- sum(tinyRNA_miRNAs_UP_output_H2[,10])
transitions_at_20 <- sum(tinyRNA_miRNAs_UP_output_H2[,11])

tinyRNA_miRNAs_UP_H2_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                        transitions_at_4, 
                                                        transitions_at_6,
                                                        transitions_at_8, 
                                                        transitions_at_10, 
                                                        transitions_at_12, 
                                                        transitions_at_14, 
                                                        transitions_at_16, 
                                                        transitions_at_18, 
                                                        transitions_at_20)


#DOWN
tinyRNA_miRNAs_DOWN_output_H2<- c()
for(i in 1:nrow(tinyRNA_H2epimutations_miRNAs)){
  tinyRNA_miRNAs_DOWN_output_H2 <-rbind(tinyRNA_miRNAs_DOWN_output_H2, DOWN_transition_func(tinyRNA_H2epimutations_miRNAs[i,], input_name=row.names(tinyRNA_H2epimutations_miRNAs)[i]))}
colnames(tinyRNA_miRNAs_DOWN_output_H2) <- c("0", "2", "4", "6", "8", "10", "12", "14", "16", "18", "20")
row.names(tinyRNA_miRNAs_DOWN_output_H2) <- rownames(tinyRNA_H2epimutations_miRNAs)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tinyRNA_miRNAs_DOWN_output_H2[,2])
transitions_at_4 <- sum(tinyRNA_miRNAs_DOWN_output_H2[,3])
transitions_at_6 <- sum(tinyRNA_miRNAs_DOWN_output_H2[,4])
transitions_at_8 <- sum(tinyRNA_miRNAs_DOWN_output_H2[,5])
transitions_at_10 <- sum(tinyRNA_miRNAs_DOWN_output_H2[,6])
transitions_at_12 <- sum(tinyRNA_miRNAs_DOWN_output_H2[,7])
transitions_at_14 <- sum(tinyRNA_miRNAs_DOWN_output_H2[,8])
transitions_at_16 <- sum(tinyRNA_miRNAs_DOWN_output_H2[,9])
transitions_at_18 <- sum(tinyRNA_miRNAs_DOWN_output_H2[,10])
transitions_at_20 <- sum(tinyRNA_miRNAs_DOWN_output_H2[,11])

tinyRNA_miRNAs_DOWN_H2_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                          transitions_at_4, 
                                                          transitions_at_6,
                                                          transitions_at_8, 
                                                          transitions_at_10, 
                                                          transitions_at_12, 
                                                          transitions_at_14, 
                                                          transitions_at_16, 
                                                          transitions_at_18, 
                                                          transitions_at_20)

#Plot construction
names <- rownames(tinyRNA_miRNAs_UP_output_C1_Table_of_new_epimutations)
rownames(tinyRNA_miRNAs_UP_output_C1_Table_of_new_epimutations) <- NULL
tinyRNA_miRNAs_UP_output_C1_Table_of_new_epimutations <- cbind(names,tinyRNA_miRNAs_UP_output_C1_Table_of_new_epimutations)
colnames(tinyRNA_miRNAs_UP_output_C1_Table_of_new_epimutations) <- c("Transitions","Up")
names <- rownames(tinyRNA_miRNAs_DOWN_output_C1_Table_of_new_epimutations)
rownames(tinyRNA_miRNAs_DOWN_output_C1_Table_of_new_epimutations) <- NULL
tinyRNA_miRNAs_DOWN_output_C1_Table_of_new_epimutations <- cbind(names,tinyRNA_miRNAs_DOWN_output_C1_Table_of_new_epimutations)
colnames(tinyRNA_miRNAs_DOWN_output_C1_Table_of_new_epimutations) <- c("Transitions","Down")

tinyRNA_miRNAs_new_epimut_C1 <- merge(tinyRNA_miRNAs_UP_output_C1_Table_of_new_epimutations,tinyRNA_miRNAs_DOWN_output_C1_Table_of_new_epimutations)
tinyRNA_miRNAs_new_epimut_C1 <- tinyRNA_miRNAs_new_epimut_C1 %>%
  # Creating an empty column:
  add_column(Lineage = "C1", .after="Transitions")

names <- rownames(tinyRNA_miRNAs_UP_output_C2_Table_of_new_epimutations)
rownames(tinyRNA_miRNAs_UP_output_C2_Table_of_new_epimutations) <- NULL
tinyRNA_miRNAs_UP_output_C2_Table_of_new_epimutations <- cbind(names,tinyRNA_miRNAs_UP_output_C2_Table_of_new_epimutations)
colnames(tinyRNA_miRNAs_UP_output_C2_Table_of_new_epimutations) <- c("Transitions","Up")
names <- rownames(tinyRNA_miRNAs_DOWN_output_C2_Table_of_new_epimutations)
rownames(tinyRNA_miRNAs_DOWN_output_C2_Table_of_new_epimutations) <- NULL
tinyRNA_miRNAs_DOWN_output_C2_Table_of_new_epimutations <- cbind(names,tinyRNA_miRNAs_DOWN_output_C2_Table_of_new_epimutations)
colnames(tinyRNA_miRNAs_DOWN_output_C2_Table_of_new_epimutations) <- c("Transitions","Down")

tinyRNA_miRNAs_new_epimut_C2 <- merge(tinyRNA_miRNAs_UP_output_C2_Table_of_new_epimutations,tinyRNA_miRNAs_DOWN_output_C2_Table_of_new_epimutations)
tinyRNA_miRNAs_new_epimut_C2 <- tinyRNA_miRNAs_new_epimut_C2 %>%
  # Creating an empty column:
  add_column(Lineage = "C2", .after="Transitions")

tiny_miRNAs_C <- rbind(tinyRNA_miRNAs_new_epimut_C1, tinyRNA_miRNAs_new_epimut_C2)
tiny_miRNAs_C <- tiny_miRNAs_C %>%
  # Creating an empty column:
  add_column(Condition = "Control", .after="Transitions")

names <- rownames(tinyRNA_miRNAs_UP_L1_Table_of_new_epimutations)
rownames(tinyRNA_miRNAs_UP_L1_Table_of_new_epimutations) <- NULL
tinyRNA_miRNAs_UP_L1_Table_of_new_epimutations <- cbind(names,tinyRNA_miRNAs_UP_L1_Table_of_new_epimutations)
colnames(tinyRNA_miRNAs_UP_L1_Table_of_new_epimutations) <- c("Transitions","Up")
names <- rownames(tinyRNA_miRNAs_DOWN_L1_Table_of_new_epimutations)
rownames(tinyRNA_miRNAs_DOWN_L1_Table_of_new_epimutations) <- NULL
tinyRNA_miRNAs_DOWN_L1_Table_of_new_epimutations <- cbind(names,tinyRNA_miRNAs_DOWN_L1_Table_of_new_epimutations)
colnames(tinyRNA_miRNAs_DOWN_L1_Table_of_new_epimutations) <- c("Transitions","Down")

tinyRNA_miRNAs_new_epimut_L1 <- merge(tinyRNA_miRNAs_UP_L1_Table_of_new_epimutations,tinyRNA_miRNAs_DOWN_L1_Table_of_new_epimutations)
tinyRNA_miRNAs_new_epimut_L1 <- tinyRNA_miRNAs_new_epimut_L1 %>%
  # Creating an empty column:
  add_column(Lineage = "L1", .after="Transitions")
names <- rownames(tinyRNA_miRNAs_UP_L2_Table_of_new_epimutations)
rownames(tinyRNA_miRNAs_UP_L2_Table_of_new_epimutations) <- NULL
tinyRNA_miRNAs_UP_L2_Table_of_new_epimutations <- cbind(names,tinyRNA_miRNAs_UP_L2_Table_of_new_epimutations)
colnames(tinyRNA_miRNAs_UP_L2_Table_of_new_epimutations) <- c("Transitions","Up")
names <- rownames(tinyRNA_miRNAs_DOWN_L2_Table_of_new_epimutations)
rownames(tinyRNA_miRNAs_DOWN_L2_Table_of_new_epimutations) <- NULL
tinyRNA_miRNAs_DOWN_L2_Table_of_new_epimutations <- cbind(names,tinyRNA_miRNAs_DOWN_L2_Table_of_new_epimutations)
colnames(tinyRNA_miRNAs_DOWN_L2_Table_of_new_epimutations) <- c("Transitions","Down")

tinyRNA_miRNAs_new_epimut_L2 <- merge(tinyRNA_miRNAs_UP_L2_Table_of_new_epimutations,tinyRNA_miRNAs_DOWN_L2_Table_of_new_epimutations)
tinyRNA_miRNAs_new_epimut_L2 <- tinyRNA_miRNAs_new_epimut_L2 %>%
  # Creating an empty column:
  add_column(Lineage = "L2", .after="Transitions")

tiny_miRNAs_L <- rbind(tinyRNA_miRNAs_new_epimut_L1, tinyRNA_miRNAs_new_epimut_L2)
tiny_miRNAs_L <- tiny_miRNAs_L %>%
  # Creating an empty column:
  add_column(Condition = "Low dose", .after="Transitions")
names <- rownames(tinyRNA_miRNAs_UP_H1_Table_of_new_epimutations)
rownames(tinyRNA_miRNAs_UP_H1_Table_of_new_epimutations) <- NULL
tinyRNA_miRNAs_UP_H1_Table_of_new_epimutations <- cbind(names,tinyRNA_miRNAs_UP_H1_Table_of_new_epimutations)
colnames(tinyRNA_miRNAs_UP_H1_Table_of_new_epimutations) <- c("Transitions","Up")
names <- rownames(tinyRNA_miRNAs_DOWN_H1_Table_of_new_epimutations)
rownames(tinyRNA_miRNAs_DOWN_H1_Table_of_new_epimutations) <- NULL
tinyRNA_miRNAs_DOWN_H1_Table_of_new_epimutations <- cbind(names,tinyRNA_miRNAs_DOWN_H1_Table_of_new_epimutations)
colnames(tinyRNA_miRNAs_DOWN_H1_Table_of_new_epimutations) <- c("Transitions","Down")

tinyRNA_miRNAs_new_epimut_H1 <- merge(tinyRNA_miRNAs_UP_H1_Table_of_new_epimutations,tinyRNA_miRNAs_DOWN_H1_Table_of_new_epimutations)
tinyRNA_miRNAs_new_epimut_H1 <- tinyRNA_miRNAs_new_epimut_H1 %>%
  # Creating an empty column:
  add_column(Lineage = "H1", .after="Transitions")
names <- rownames(tinyRNA_miRNAs_UP_H2_Table_of_new_epimutations)
rownames(tinyRNA_miRNAs_UP_H2_Table_of_new_epimutations) <- NULL
tinyRNA_miRNAs_UP_H2_Table_of_new_epimutations <- cbind(names,tinyRNA_miRNAs_UP_H2_Table_of_new_epimutations)
colnames(tinyRNA_miRNAs_UP_H2_Table_of_new_epimutations) <- c("Transitions","Up")
names <- rownames(tinyRNA_miRNAs_DOWN_H2_Table_of_new_epimutations)
rownames(tinyRNA_miRNAs_DOWN_H2_Table_of_new_epimutations) <- NULL
tinyRNA_miRNAs_DOWN_H2_Table_of_new_epimutations <- cbind(names,tinyRNA_miRNAs_DOWN_H2_Table_of_new_epimutations)
colnames(tinyRNA_miRNAs_DOWN_H2_Table_of_new_epimutations) <- c("Transitions","Down")

tinyRNA_miRNAs_new_epimut_H2 <- merge(tinyRNA_miRNAs_UP_H2_Table_of_new_epimutations,tinyRNA_miRNAs_DOWN_H2_Table_of_new_epimutations)
tinyRNA_miRNAs_new_epimut_H2 <- tinyRNA_miRNAs_new_epimut_H2 %>%
  # Creating an empty column:
  add_column(Lineage = "H2", .after="Transitions")

tiny_miRNAs_H <- rbind(tinyRNA_miRNAs_new_epimut_H1, tinyRNA_miRNAs_new_epimut_H2)
tiny_miRNAs_H <- tiny_miRNAs_H %>%
  # Creating an empty column:
  add_column(Condition = "High dose", .after="Transitions")

tiny_miRNAs_newepi<-rbind(tiny_miRNAs_C,tiny_miRNAs_L,tiny_miRNAs_H)
tiny_miRNAs_newepi$Up <- as.numeric(tiny_miRNAs_newepi$Up)
tiny_miRNAs_newepi$Down <- as.numeric(tiny_miRNAs_newepi$Down)
tiny_miRNAs_newepi$Total=rowSums(cbind(tiny_miRNAs_newepi$Up,tiny_miRNAs_newepi$Down),na.rm=FALSE)

#Test for statistical significance
ggdensity(tiny_miRNAs_newepi$Total, 
          main = "Density plot of epimutations_number",
          xlab = "Total number of epimutations_number")
shapiro.test(tiny_miRNAs_newepi$Total)
kruskal.test(Total ~ Condition, data = tiny_miRNAs_newepi)
dunnTest(Total ~ Condition, data = tiny_miRNAs_newepi)

#Fig.4.E
tiny_miRNAs_newepi$Condition<-fct_relevel(tiny_miRNAs_newepi$Condition,c("Control","Low dose","High dose"))
tiny_miRNAs_allLineage_rate <- ggplot(tiny_miRNAs_newepi, aes(x=Condition, y=Total, color = Condition,group=Condition)) + 
  geom_boxplot(fatten = 1, lwd = 1, width=0.5)+
  scale_color_manual(values=c("cornflowerblue", "darkgreen", "red"))+
  labs(y = "Number of new epimutations", x = "\nCondition")+
  geom_dotplot(binaxis='y', binwidth = 0.25, stackdir='center',  dotsize=30, binpositions = "all", stackgroups = TRUE)+
  theme_bw()+
  theme_linedraw()+
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(color = "black", size = 25, angle = 0, face = "plain"),
        axis.text.y = element_text(color = "black", size = 25, angle = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 20, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, face = "bold"))+
  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))+  
  scale_fill_brewer(palette="Pastel2")+
  ggtitle(paste(""))
tiny_miRNAs_allLineage_rate

#Sup.Fig.5.E
Generation<-c("10","12","14","16","18","2","20","4","6","8",
              "10","12","14","16","18","2","20","4","6","8",
              "10","12","14","16","18","2","20","4","6","8",
              "10","12","14","16","18","2","20","4","6","8",
              "10","12","14","16","18","2","20","4","6","8",
              "10","12","14","16","18","2","20","4","6","8")
tiny_miRNAs_newepi_detailed<-cbind(tiny_miRNAs_newepi,Generation)
tiny_miRNAs_newepi_detailed$Generation <- fct_relevel(tiny_miRNAs_newepi_detailed$Generation, c("2", "4", "6","8","10","12","14","16","18","20"))
tiny_miRNAs_newepi_detailed$Condition <- fct_relevel(tiny_miRNAs_newepi_detailed$Condition, c("Control", "Low dose", "High dose"))
tiny_miRNAs_newepi_detailed_plot <- ggplot(tiny_miRNAs_newepi_detailed, aes(x=Generation , y=Total, color = Condition)) + 
  geom_boxplot(fatten = 1, lwd = 1, width=0.5)+
  scale_color_manual(values=c("cornflowerblue", "darkgreen", "red"))+
  labs(y = "Number of new epimutations", x = "\nGeneration")+
  geom_jitter(aes(shape=Lineage, size=8),
              position=position_jitter(width = 0.2,
                                       height = 0.2))+
  theme_bw()+
  theme(axis.text.x = element_text(color = "black", size = 20, angle = 0, face = "plain"),
        axis.text.y = element_text(color = "black", size = 20, angle = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 16, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 16, angle = 90, face = "bold"))+
  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))+  
  scale_fill_brewer(palette="Pastel2")+
  ggtitle(paste(""))
tiny_miRNAs_newepi_detailed_plot 

#Get raw data
write.xlsx(tiny_miRNAs_newepi_detailed,"Data_Fig_4_E_&_Sup_Fig_5_E.xlsx")
#--------------------------
#####Fig.4.F & Sup.Fig.6.F - miRNAs epimutations duration####
tiny_miRNAs_epimut<-tiny_all_epimut %>% filter(grepl('miRNA', gene_name))
#Fig.4.F
CoxMod<-coxph(Surv(length,complete)~Condition,data=tiny_miRNAs_epimut)
ggforest(CoxMod, data=tiny_miRNAs_epimut)

ggsurvplot(survfit(CoxMod,data=tiny_miRNAs_epimut), palette = "#2E9FDF",cof.int=TRUE,
           ggtheme = theme_minimal())
condition_df <- with(tiny_miRNAs_epimut,
                     data.frame(Condition = c("Control","Low dose","High dose")
                     )
)
condition_df
fit <- survfit(CoxMod, newdata = condition_df)
fit <- survfit(Surv(length, complete) ~ Condition, data = tiny_miRNAs_epimut)
ggsurvplot <- ggsurvplot(fit,  tiny_miRNAs_epimut, censor = T, break.time.by= 1, pval=TRUE,pval.coord = c(0, 0.03),
                         font.main = c(16, "bold"),
                         palette = c("cornflowerblue", "darkgreen","red"),
                         font.x = c(20,"bold"),
                         font.y = c(20,"bold"),
                         font.legend = 20,
                         legend.labs=c("Control","Low dose","High dose"),
                         font.tickslab = 18)+
  # surv.median.line = "hv")+
  ggtitle("Changes survival all conditions")+
  xlab("Time (generations)")
ggsurvplot$plot +theme(plot.title = element_text(hjust = 0.5))

#Sup.Fig.5.F
CoxMod<-coxph(Surv(length,complete)~Lineage,data=tiny_miRNAs_epimut)
ggforest(model=CoxMod,data=tiny_miRNAs_epimut,fontsize = 0.8, noDigits = 2)
fit_2 <- survfit(Surv(length, complete) ~ Lineage, data = tiny_miRNAs_epimut)
ggsurvplot <- ggsurvplot(fit_2,  tiny_miRNAs_epimut, censor = T, break.time.by= 1, pval=TRUE,pval.coord = c(0, 0.03),                         font.main = c(16, "bold"),
                         font.x = c(20,"bold"),
                         palette = c("#3399FF","blue","red","2E9FDF", "green","#006600"),
                         font.y = c(20,"bold"),
                         font.legend = 20,
                         legend.labs=c("C1","C2","H1","H2","L1","L2"),
                         font.tickslab = 18)+
  ggtitle("Changes survival all conditions")+
  xlab("Time (generations)")
ggsurvplot$plot +theme(plot.title = element_text(hjust = 0.5)) 

#Get raw data
write.xlsx(tiny_miRNAs_epimut,"Data_Fig_4_F_&_Sup_Fig_5_F.xlsx")
#--------------------------
#####Fig.4.G & Sup.Fig.6.G - 26G-RNAs epimutations####
#Control1 new epi
tinyRNA_C1epimutations_26G<-subset(tinyRNA_C1epimutations,tinyRNA_C1epimutations$Tag=="26G")
#UP
tinyRNA_26G_UP_output_C1<- c()
tinyRNA_ID <- paste(tinyRNA_C1epimutations_26G$Feature.ID, tinyRNA_C1epimutations_26G$Tag,tinyRNA_C1epimutations_26G$Feature.Name,tinyRNA_C1epimutations_26G$Feature.Class, sep="_")
row.names(tinyRNA_C1epimutations_26G)<-tinyRNA_ID
tinyRNA_C1epimutations_26G<-tinyRNA_C1epimutations_26G[,c(7:17)]
for(i in 1:nrow(tinyRNA_C1epimutations_26G)){
  tinyRNA_26G_UP_output_C1 <-rbind(tinyRNA_26G_UP_output_C1, UP_transition_func(tinyRNA_C1epimutations_26G[i,], input_name=row.names(tinyRNA_C1epimutations_26G)[i]))}
colnames(tinyRNA_26G_UP_output_C1) <- c("0", "2", "4", "6", "8", "10", "12", "14", "16", "18", "20")
row.names(tinyRNA_26G_UP_output_C1) <- rownames(tinyRNA_C1epimutations_26G)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tinyRNA_26G_UP_output_C1[,2])
transitions_at_4 <- sum(tinyRNA_26G_UP_output_C1[,3])
transitions_at_6 <- sum(tinyRNA_26G_UP_output_C1[,4])
transitions_at_8 <- sum(tinyRNA_26G_UP_output_C1[,5])
transitions_at_10 <- sum(tinyRNA_26G_UP_output_C1[,6])
transitions_at_12 <- sum(tinyRNA_26G_UP_output_C1[,7])
transitions_at_14 <- sum(tinyRNA_26G_UP_output_C1[,8])
transitions_at_16 <- sum(tinyRNA_26G_UP_output_C1[,9])
transitions_at_18 <- sum(tinyRNA_26G_UP_output_C1[,10])
transitions_at_20 <- sum(tinyRNA_26G_UP_output_C1[,11])

tinyRNA_26G_UP_output_C1_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                            transitions_at_4, 
                                                            transitions_at_6,
                                                            transitions_at_8, 
                                                            transitions_at_10, 
                                                            transitions_at_12, 
                                                            transitions_at_14, 
                                                            transitions_at_16, 
                                                            transitions_at_18, 
                                                            transitions_at_20)

#DOWN
tinyRNA_26G_DOWN_output_C1<- c()
for(i in 1:nrow(tinyRNA_C1epimutations_26G)){
  tinyRNA_26G_DOWN_output_C1 <-rbind(tinyRNA_26G_DOWN_output_C1, DOWN_transition_func(tinyRNA_C1epimutations_26G[i,], input_name=row.names(tinyRNA_C1epimutations_26G)[i]))}
colnames(tinyRNA_26G_DOWN_output_C1) <- c("0", "2", "4", "6", "8", "10", "12", "14", "16", "18", "20")
row.names(tinyRNA_26G_DOWN_output_C1) <- rownames(tinyRNA_C1epimutations_26G)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tinyRNA_26G_DOWN_output_C1[,2])
transitions_at_4 <- sum(tinyRNA_26G_DOWN_output_C1[,3])
transitions_at_6 <- sum(tinyRNA_26G_DOWN_output_C1[,4])
transitions_at_8 <- sum(tinyRNA_26G_DOWN_output_C1[,5])
transitions_at_10 <- sum(tinyRNA_26G_DOWN_output_C1[,6])
transitions_at_12 <- sum(tinyRNA_26G_DOWN_output_C1[,7])
transitions_at_14 <- sum(tinyRNA_26G_DOWN_output_C1[,8])
transitions_at_16 <- sum(tinyRNA_26G_DOWN_output_C1[,9])
transitions_at_18 <- sum(tinyRNA_26G_DOWN_output_C1[,10])
transitions_at_20 <- sum(tinyRNA_26G_DOWN_output_C1[,11])

tinyRNA_26G_DOWN_output_C1_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                              transitions_at_4, 
                                                              transitions_at_6,
                                                              transitions_at_8, 
                                                              transitions_at_10, 
                                                              transitions_at_12, 
                                                              transitions_at_14, 
                                                              transitions_at_16, 
                                                              transitions_at_18, 
                                                              transitions_at_20)

#Control2 new epi
tinyRNA_C2epimutations_26G<-subset(tinyRNA_C2epimutations,tinyRNA_C2epimutations$Tag=="26G")
#UP
tinyRNA_26G_UP_output_C2<- c()
tinyRNA_ID <- paste(tinyRNA_C2epimutations_26G$Feature.ID, tinyRNA_C2epimutations_26G$Tag,tinyRNA_C2epimutations_26G$Feature.Name,tinyRNA_C2epimutations_26G$Feature.Class, sep="_")
row.names(tinyRNA_C2epimutations_26G)<-tinyRNA_ID
tinyRNA_C2epimutations_26G<-tinyRNA_C2epimutations_26G[,c(7:16)]
for(i in 1:nrow(tinyRNA_C2epimutations_26G)){
  tinyRNA_26G_UP_output_C2 <-rbind(tinyRNA_26G_UP_output_C2, UP_transition_func(tinyRNA_C2epimutations_26G[i,], input_name=row.names(tinyRNA_C2epimutations_26G)[i]))}
colnames(tinyRNA_26G_UP_output_C2) <- c("0", "2", "4","8", "10", "12", "14", "16", "18", "20")
row.names(tinyRNA_26G_UP_output_C2) <- rownames(tinyRNA_C2epimutations_26G)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tinyRNA_26G_UP_output_C2[,2])
transitions_at_4 <- sum(tinyRNA_26G_UP_output_C2[,3])
transitions_at_6 <- NA
transitions_at_8 <- sum(tinyRNA_26G_UP_output_C2[,4])
transitions_at_10 <- sum(tinyRNA_26G_UP_output_C2[,5])
transitions_at_12 <- sum(tinyRNA_26G_UP_output_C2[,6])
transitions_at_14 <- sum(tinyRNA_26G_UP_output_C2[,7])
transitions_at_16 <- sum(tinyRNA_26G_UP_output_C2[,8])
transitions_at_18 <- sum(tinyRNA_26G_UP_output_C2[,9])
transitions_at_20 <- sum(tinyRNA_26G_UP_output_C2[,10])

tinyRNA_26G_UP_output_C2_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                            transitions_at_4, 
                                                            transitions_at_6,
                                                            transitions_at_8, 
                                                            transitions_at_10, 
                                                            transitions_at_12, 
                                                            transitions_at_14, 
                                                            transitions_at_16, 
                                                            transitions_at_18, 
                                                            transitions_at_20)

#DOWN
tinyRNA_26G_DOWN_output_C2<- c()
for(i in 1:nrow(tinyRNA_C2epimutations_26G)){
  tinyRNA_26G_DOWN_output_C2 <-rbind(tinyRNA_26G_DOWN_output_C2, DOWN_transition_func(tinyRNA_C2epimutations_26G[i,], input_name=row.names(tinyRNA_C2epimutations_26G)[i]))}
colnames(tinyRNA_26G_DOWN_output_C2) <- c("0", "2", "4", "8", "10", "12", "14", "16", "18", "20")
row.names(tinyRNA_26G_DOWN_output_C2) <- rownames(tinyRNA_C2epimutations_26G)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tinyRNA_26G_DOWN_output_C2[,2])
transitions_at_4 <- sum(tinyRNA_26G_DOWN_output_C2[,3])
transitions_at_6 <- NA
transitions_at_8 <- sum(tinyRNA_26G_DOWN_output_C2[,4])
transitions_at_10 <- sum(tinyRNA_26G_DOWN_output_C2[,5])
transitions_at_12 <- sum(tinyRNA_26G_DOWN_output_C2[,6])
transitions_at_14 <- sum(tinyRNA_26G_DOWN_output_C2[,7])
transitions_at_16 <- sum(tinyRNA_26G_DOWN_output_C2[,8])
transitions_at_18 <- sum(tinyRNA_26G_DOWN_output_C2[,9])
transitions_at_20 <- sum(tinyRNA_26G_DOWN_output_C2[,10])

tinyRNA_26G_DOWN_output_C2_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                              transitions_at_4, 
                                                              transitions_at_6,
                                                              transitions_at_8, 
                                                              transitions_at_10, 
                                                              transitions_at_12, 
                                                              transitions_at_14, 
                                                              transitions_at_16, 
                                                              transitions_at_18, 
                                                              transitions_at_20)

#Low1 new epi
tinyRNA_L1epimutations_26G<-subset(tinyRNA_L1epimutations,tinyRNA_L1epimutations$Tag=="26G")
#UP
tinyRNA_26G_UP_output_L1<- c()
tinyRNA_ID <- paste(tinyRNA_L1epimutations_26G$Feature.ID, tinyRNA_L1epimutations_26G$Tag,tinyRNA_L1epimutations_26G$Feature.Name,tinyRNA_L1epimutations_26G$Feature.Class, sep="_")
row.names(tinyRNA_L1epimutations_26G)<-tinyRNA_ID
tinyRNA_L1epimutations_26G<-tinyRNA_L1epimutations_26G[,c(7:16)]
for(i in 1:nrow(tinyRNA_L1epimutations_26G)){
  tinyRNA_26G_UP_output_L1 <-rbind(tinyRNA_26G_UP_output_L1, UP_transition_func(tinyRNA_L1epimutations_26G[i,], input_name=row.names(tinyRNA_L1epimutations_26G)[i]))}
colnames(tinyRNA_26G_UP_output_L1) <- c("0", "2", "4","6", "8", "12", "14", "16", "18", "20")
row.names(tinyRNA_26G_UP_output_L1) <- rownames(tinyRNA_L1epimutations_26G)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tinyRNA_26G_UP_output_L1[,2])
transitions_at_4 <- sum(tinyRNA_26G_UP_output_L1[,3])
transitions_at_6 <- sum(tinyRNA_26G_UP_output_L1[,4])
transitions_at_8 <- sum(tinyRNA_26G_UP_output_L1[,5])
transitions_at_10 <- NA
transitions_at_12 <- sum(tinyRNA_26G_UP_output_L1[,6])
transitions_at_14 <- sum(tinyRNA_26G_UP_output_L1[,7])
transitions_at_16 <- sum(tinyRNA_26G_UP_output_L1[,8])
transitions_at_18 <- sum(tinyRNA_26G_UP_output_L1[,9])
transitions_at_20 <- sum(tinyRNA_26G_UP_output_L1[,10])

tinyRNA_26G_UP_L1_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                     transitions_at_4, 
                                                     transitions_at_6,
                                                     transitions_at_8, 
                                                     transitions_at_10, 
                                                     transitions_at_12, 
                                                     transitions_at_14, 
                                                     transitions_at_16, 
                                                     transitions_at_18, 
                                                     transitions_at_20)

#DOWN
tinyRNA_26G_DOWN_output_L1<- c()
for(i in 1:nrow(tinyRNA_L1epimutations_26G)){
  tinyRNA_26G_DOWN_output_L1 <-rbind(tinyRNA_26G_DOWN_output_L1, DOWN_transition_func(tinyRNA_L1epimutations_26G[i,], input_name=row.names(tinyRNA_L1epimutations_26G)[i]))}
colnames(tinyRNA_26G_DOWN_output_L1) <- c("0", "2", "4", "6","8", "12", "14", "16", "18", "20")
row.names(tinyRNA_26G_DOWN_output_L1) <- rownames(tinyRNA_L1epimutations_26G)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tinyRNA_26G_DOWN_output_L1[,2])
transitions_at_4 <- sum(tinyRNA_26G_DOWN_output_L1[,3])
transitions_at_6 <- sum(tinyRNA_26G_DOWN_output_L1[,4])
transitions_at_8 <- sum(tinyRNA_26G_DOWN_output_L1[,5])
transitions_at_10 <- NA
transitions_at_12 <- sum(tinyRNA_26G_DOWN_output_L1[,6])
transitions_at_14 <- sum(tinyRNA_26G_DOWN_output_L1[,7])
transitions_at_16 <- sum(tinyRNA_26G_DOWN_output_L1[,8])
transitions_at_18 <- sum(tinyRNA_26G_DOWN_output_L1[,9])
transitions_at_20 <- sum(tinyRNA_26G_DOWN_output_L1[,10])

tinyRNA_26G_DOWN_L1_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                       transitions_at_4, 
                                                       transitions_at_6,
                                                       transitions_at_8, 
                                                       transitions_at_10, 
                                                       transitions_at_12, 
                                                       transitions_at_14, 
                                                       transitions_at_16, 
                                                       transitions_at_18, 
                                                       transitions_at_20)

#Low2 new epi
tinyRNA_L2epimutations_26G<-subset(tinyRNA_L2epimutations,tinyRNA_L2epimutations$Tag=="26G")
#UP

tinyRNA_26G_UP_output_L2<- c()
tinyRNA_ID <- paste(tinyRNA_L2epimutations_26G$Feature.ID, tinyRNA_L2epimutations_26G$Tag,tinyRNA_L2epimutations_26G$Feature.Name,tinyRNA_L2epimutations_26G$Feature.Class, sep="_")
row.names(tinyRNA_L2epimutations_26G)<-tinyRNA_ID
tinyRNA_L2epimutations_26G<-tinyRNA_L2epimutations_26G[,c(7:16)]
for(i in 1:nrow(tinyRNA_L2epimutations_26G)){
  tinyRNA_26G_UP_output_L2 <-rbind(tinyRNA_26G_UP_output_L2, UP_transition_func(tinyRNA_L2epimutations_26G[i,], input_name=row.names(tinyRNA_L2epimutations_26G)[i]))}
colnames(tinyRNA_26G_UP_output_L2) <- c("0", "2", "4", "8","10", "12", "14", "16", "18", "20")
row.names(tinyRNA_26G_UP_output_L2) <- rownames(tinyRNA_L2epimutations_26G)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tinyRNA_26G_UP_output_L2[,2])
transitions_at_4 <- sum(tinyRNA_26G_UP_output_L2[,3])
transitions_at_6 <- NA
transitions_at_8 <- sum(tinyRNA_26G_UP_output_L2[,4])
transitions_at_10 <- sum(tinyRNA_26G_UP_output_L2[,5])
transitions_at_12 <- sum(tinyRNA_26G_UP_output_L2[,6])
transitions_at_14 <- sum(tinyRNA_26G_UP_output_L2[,7])
transitions_at_16 <- sum(tinyRNA_26G_UP_output_L2[,8])
transitions_at_18 <- sum(tinyRNA_26G_UP_output_L2[,9])
transitions_at_20 <- sum(tinyRNA_26G_UP_output_L2[,10])

tinyRNA_26G_UP_L2_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                     transitions_at_4, 
                                                     transitions_at_6,
                                                     transitions_at_8, 
                                                     transitions_at_10, 
                                                     transitions_at_12, 
                                                     transitions_at_14, 
                                                     transitions_at_16, 
                                                     transitions_at_18, 
                                                     transitions_at_20)


#DOWN
tinyRNA_26G_DOWN_output_L2<- c()
for(i in 1:nrow(tinyRNA_L2epimutations_26G)){
  tinyRNA_26G_DOWN_output_L2 <-rbind(tinyRNA_26G_DOWN_output_L2, DOWN_transition_func(tinyRNA_L2epimutations_26G[i,], input_name=row.names(tinyRNA_L2epimutations_26G)[i]))}
colnames(tinyRNA_26G_DOWN_output_L2) <- c("0", "2", "4", "6","8", "12", "14", "16", "18", "20")
row.names(tinyRNA_26G_DOWN_output_L2) <- rownames(tinyRNA_L2epimutations_26G)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tinyRNA_26G_DOWN_output_L2[,2])
transitions_at_4 <- sum(tinyRNA_26G_DOWN_output_L2[,3])
transitions_at_6 <- NA
transitions_at_8 <- sum(tinyRNA_26G_DOWN_output_L2[,4])
transitions_at_10 <- sum(tinyRNA_26G_DOWN_output_L2[,5])
transitions_at_12 <- sum(tinyRNA_26G_DOWN_output_L2[,6])
transitions_at_14 <- sum(tinyRNA_26G_DOWN_output_L2[,7])
transitions_at_16 <- sum(tinyRNA_26G_DOWN_output_L2[,8])
transitions_at_18 <- sum(tinyRNA_26G_DOWN_output_L2[,9])
transitions_at_20 <- sum(tinyRNA_26G_DOWN_output_L2[,10])

tinyRNA_26G_DOWN_L2_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                       transitions_at_4, 
                                                       transitions_at_6,
                                                       transitions_at_8, 
                                                       transitions_at_10, 
                                                       transitions_at_12, 
                                                       transitions_at_14, 
                                                       transitions_at_16, 
                                                       transitions_at_18, 
                                                       transitions_at_20)

#High1 new epi
tinyRNA_H1epimutations_26G<-subset(tinyRNA_H1epimutations,tinyRNA_H1epimutations$Tag=="26G")

tinyRNA_26G_UP_output_H1<- c()
tinyRNA_ID <- paste(tinyRNA_H1epimutations_26G$Feature.ID, tinyRNA_H1epimutations_26G$Tag,tinyRNA_H1epimutations_26G$Feature.Name,tinyRNA_H1epimutations_26G$Feature.Class, sep="_")
row.names(tinyRNA_H1epimutations_26G)<-tinyRNA_ID
tinyRNA_H1epimutations_26G<-tinyRNA_H1epimutations_26G[,c(7:16)]
for(i in 1:nrow(tinyRNA_H1epimutations_26G)){
  tinyRNA_26G_UP_output_H1 <-rbind(tinyRNA_26G_UP_output_H1, UP_transition_func(tinyRNA_H1epimutations_26G[i,], input_name=row.names(tinyRNA_H1epimutations_26G)[i]))}
colnames(tinyRNA_26G_UP_output_H1) <- c("0", "2", "4", "8","10", "12", "14", "16", "18", "20")
row.names(tinyRNA_26G_UP_output_H1) <- rownames(tinyRNA_H1epimutations_26G)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tinyRNA_26G_UP_output_H1[,2])
transitions_at_4 <- sum(tinyRNA_26G_UP_output_H1[,3])
transitions_at_6 <- NA
transitions_at_8 <- sum(tinyRNA_26G_UP_output_H1[,4])
transitions_at_10 <- sum(tinyRNA_26G_UP_output_H1[,5])
transitions_at_12 <- sum(tinyRNA_26G_UP_output_H1[,6])
transitions_at_14 <- sum(tinyRNA_26G_UP_output_H1[,7])
transitions_at_16 <- sum(tinyRNA_26G_UP_output_H1[,8])
transitions_at_18 <- sum(tinyRNA_26G_UP_output_H1[,9])
transitions_at_20 <- sum(tinyRNA_26G_UP_output_H1[,10])

tinyRNA_26G_UP_H1_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                     transitions_at_4, 
                                                     transitions_at_6,
                                                     transitions_at_8, 
                                                     transitions_at_10, 
                                                     transitions_at_12, 
                                                     transitions_at_14, 
                                                     transitions_at_16, 
                                                     transitions_at_18, 
                                                     transitions_at_20)


#DOWN
tinyRNA_26G_DOWN_output_H1<- c()
for(i in 1:nrow(tinyRNA_H1epimutations_26G)){
  tinyRNA_26G_DOWN_output_H1 <-rbind(tinyRNA_26G_DOWN_output_H1, DOWN_transition_func(tinyRNA_H1epimutations_26G[i,], input_name=row.names(tinyRNA_H1epimutations_26G)[i]))}
colnames(tinyRNA_26G_DOWN_output_H1) <- c("0", "2", "4", "8", "10", "12", "14", "16", "18", "20")
row.names(tinyRNA_26G_DOWN_output_H1) <- rownames(tinyRNA_H1epimutations_26G)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tinyRNA_26G_DOWN_output_H1[,2])
transitions_at_4 <- sum(tinyRNA_26G_DOWN_output_H1[,3])
transitions_at_6 <- NA
transitions_at_8 <- sum(tinyRNA_26G_DOWN_output_H1[,4])
transitions_at_10 <- sum(tinyRNA_26G_DOWN_output_H1[,5])
transitions_at_12 <- sum(tinyRNA_26G_DOWN_output_H1[,6])
transitions_at_14 <- sum(tinyRNA_26G_DOWN_output_H1[,7])
transitions_at_16 <- sum(tinyRNA_26G_DOWN_output_H1[,8])
transitions_at_18 <- sum(tinyRNA_26G_DOWN_output_H1[,9])
transitions_at_20 <- sum(tinyRNA_26G_DOWN_output_H1[,10])

tinyRNA_26G_DOWN_H1_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                       transitions_at_4, 
                                                       transitions_at_6,
                                                       transitions_at_8, 
                                                       transitions_at_10, 
                                                       transitions_at_12, 
                                                       transitions_at_14, 
                                                       transitions_at_16, 
                                                       transitions_at_18, 
                                                       transitions_at_20)

#High2 new epi
tinyRNA_H2epimutations_26G<-subset(tinyRNA_H2epimutations,tinyRNA_H2epimutations$Tag=="26G")
#UP
tinyRNA_26G_UP_output_H2<- c()
tinyRNA_ID <- paste(tinyRNA_H2epimutations_26G$Feature.ID, tinyRNA_H2epimutations_26G$Tag,tinyRNA_H2epimutations_26G$Feature.Name,tinyRNA_H2epimutations_26G$Feature.Class, sep="_")
row.names(tinyRNA_H2epimutations_26G)<-tinyRNA_ID
tinyRNA_H2epimutations_26G<-tinyRNA_H2epimutations_26G[,c(7:17)]
for(i in 1:nrow(tinyRNA_H2epimutations_26G)){
  tinyRNA_26G_UP_output_H2 <-rbind(tinyRNA_26G_UP_output_H2, UP_transition_func(tinyRNA_H2epimutations_26G[i,], input_name=row.names(tinyRNA_H2epimutations_26G)[i]))}
colnames(tinyRNA_26G_UP_output_H2) <- c("0", "2", "4", "6", "8","10", "12", "14", "16", "18", "20")
row.names(tinyRNA_26G_UP_output_H2) <- rownames(tinyRNA_H2epimutations_26G)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tinyRNA_26G_UP_output_H2[,2])
transitions_at_4 <- sum(tinyRNA_26G_UP_output_H2[,3])
transitions_at_6 <- sum(tinyRNA_26G_UP_output_H2[,4])
transitions_at_8 <- sum(tinyRNA_26G_UP_output_H2[,5])
transitions_at_10 <- sum(tinyRNA_26G_UP_output_H2[,6])
transitions_at_12 <- sum(tinyRNA_26G_UP_output_H2[,7])
transitions_at_14 <- sum(tinyRNA_26G_UP_output_H2[,8])
transitions_at_16 <- sum(tinyRNA_26G_UP_output_H2[,9])
transitions_at_18 <- sum(tinyRNA_26G_UP_output_H2[,10])
transitions_at_20 <- sum(tinyRNA_26G_UP_output_H2[,11])

tinyRNA_26G_UP_H2_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                     transitions_at_4, 
                                                     transitions_at_6,
                                                     transitions_at_8, 
                                                     transitions_at_10, 
                                                     transitions_at_12, 
                                                     transitions_at_14, 
                                                     transitions_at_16, 
                                                     transitions_at_18, 
                                                     transitions_at_20)


#DOWN
tinyRNA_26G_DOWN_output_H2<- c()
for(i in 1:nrow(tinyRNA_H2epimutations_26G)){
  tinyRNA_26G_DOWN_output_H2 <-rbind(tinyRNA_26G_DOWN_output_H2, DOWN_transition_func(tinyRNA_H2epimutations_26G[i,], input_name=row.names(tinyRNA_H2epimutations_26G)[i]))}
colnames(tinyRNA_26G_DOWN_output_H2) <- c("0", "2", "4", "6", "8", "10", "12", "14", "16", "18", "20")
row.names(tinyRNA_26G_DOWN_output_H2) <- rownames(tinyRNA_H2epimutations_26G)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tinyRNA_26G_DOWN_output_H2[,2])
transitions_at_4 <- sum(tinyRNA_26G_DOWN_output_H2[,3])
transitions_at_6 <- sum(tinyRNA_26G_DOWN_output_H2[,4])
transitions_at_8 <- sum(tinyRNA_26G_DOWN_output_H2[,5])
transitions_at_10 <- sum(tinyRNA_26G_DOWN_output_H2[,6])
transitions_at_12 <- sum(tinyRNA_26G_DOWN_output_H2[,7])
transitions_at_14 <- sum(tinyRNA_26G_DOWN_output_H2[,8])
transitions_at_16 <- sum(tinyRNA_26G_DOWN_output_H2[,9])
transitions_at_18 <- sum(tinyRNA_26G_DOWN_output_H2[,10])
transitions_at_20 <- sum(tinyRNA_26G_DOWN_output_H2[,11])

tinyRNA_26G_DOWN_H2_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                       transitions_at_4, 
                                                       transitions_at_6,
                                                       transitions_at_8, 
                                                       transitions_at_10, 
                                                       transitions_at_12, 
                                                       transitions_at_14, 
                                                       transitions_at_16, 
                                                       transitions_at_18, 
                                                       transitions_at_20)

#Plot construction
names <- rownames(tinyRNA_26G_UP_output_C1_Table_of_new_epimutations)
rownames(tinyRNA_26G_UP_output_C1_Table_of_new_epimutations) <- NULL
tinyRNA_26G_UP_output_C1_Table_of_new_epimutations <- cbind(names,tinyRNA_26G_UP_output_C1_Table_of_new_epimutations)
colnames(tinyRNA_26G_UP_output_C1_Table_of_new_epimutations) <- c("Transitions","Up")
names <- rownames(tinyRNA_26G_DOWN_output_C1_Table_of_new_epimutations)
rownames(tinyRNA_26G_DOWN_output_C1_Table_of_new_epimutations) <- NULL
tinyRNA_26G_DOWN_output_C1_Table_of_new_epimutations <- cbind(names,tinyRNA_26G_DOWN_output_C1_Table_of_new_epimutations)
colnames(tinyRNA_26G_DOWN_output_C1_Table_of_new_epimutations) <- c("Transitions","Down")

tinyRNA_26G_new_epimut_C1 <- merge(tinyRNA_26G_UP_output_C1_Table_of_new_epimutations,tinyRNA_26G_DOWN_output_C1_Table_of_new_epimutations)
tinyRNA_26G_new_epimut_C1 <- tinyRNA_26G_new_epimut_C1 %>%
  # Creating an empty column:
  add_column(Lineage = "C1", .after="Transitions")

names <- rownames(tinyRNA_26G_UP_output_C2_Table_of_new_epimutations)
rownames(tinyRNA_26G_UP_output_C2_Table_of_new_epimutations) <- NULL
tinyRNA_26G_UP_output_C2_Table_of_new_epimutations <- cbind(names,tinyRNA_26G_UP_output_C2_Table_of_new_epimutations)
colnames(tinyRNA_26G_UP_output_C2_Table_of_new_epimutations) <- c("Transitions","Up")
names <- rownames(tinyRNA_26G_DOWN_output_C2_Table_of_new_epimutations)
rownames(tinyRNA_26G_DOWN_output_C2_Table_of_new_epimutations) <- NULL
tinyRNA_26G_DOWN_output_C2_Table_of_new_epimutations <- cbind(names,tinyRNA_26G_DOWN_output_C2_Table_of_new_epimutations)
colnames(tinyRNA_26G_DOWN_output_C2_Table_of_new_epimutations) <- c("Transitions","Down")

tinyRNA_26G_new_epimut_C2 <- merge(tinyRNA_26G_UP_output_C2_Table_of_new_epimutations,tinyRNA_26G_DOWN_output_C2_Table_of_new_epimutations)
tinyRNA_26G_new_epimut_C2 <- tinyRNA_26G_new_epimut_C2 %>%
  # Creating an empty column:
  add_column(Lineage = "C2", .after="Transitions")

tiny_26G_C <- rbind(tinyRNA_26G_new_epimut_C1, tinyRNA_26G_new_epimut_C2)
tiny_26G_C <- tiny_26G_C %>%
  # Creating an empty column:
  add_column(Condition = "Control", .after="Transitions")

names <- rownames(tinyRNA_26G_UP_L1_Table_of_new_epimutations)
rownames(tinyRNA_26G_UP_L1_Table_of_new_epimutations) <- NULL
tinyRNA_26G_UP_L1_Table_of_new_epimutations <- cbind(names,tinyRNA_26G_UP_L1_Table_of_new_epimutations)
colnames(tinyRNA_26G_UP_L1_Table_of_new_epimutations) <- c("Transitions","Up")
names <- rownames(tinyRNA_26G_DOWN_L1_Table_of_new_epimutations)
rownames(tinyRNA_26G_DOWN_L1_Table_of_new_epimutations) <- NULL
tinyRNA_26G_DOWN_L1_Table_of_new_epimutations <- cbind(names,tinyRNA_26G_DOWN_L1_Table_of_new_epimutations)
colnames(tinyRNA_26G_DOWN_L1_Table_of_new_epimutations) <- c("Transitions","Down")

tinyRNA_26G_new_epimut_L1 <- merge(tinyRNA_26G_UP_L1_Table_of_new_epimutations,tinyRNA_26G_DOWN_L1_Table_of_new_epimutations)
tinyRNA_26G_new_epimut_L1 <- tinyRNA_26G_new_epimut_L1 %>%
  # Creating an empty column:
  add_column(Lineage = "L1", .after="Transitions")
names <- rownames(tinyRNA_26G_UP_L2_Table_of_new_epimutations)
rownames(tinyRNA_26G_UP_L2_Table_of_new_epimutations) <- NULL
tinyRNA_26G_UP_L2_Table_of_new_epimutations <- cbind(names,tinyRNA_26G_UP_L2_Table_of_new_epimutations)
colnames(tinyRNA_26G_UP_L2_Table_of_new_epimutations) <- c("Transitions","Up")
names <- rownames(tinyRNA_26G_DOWN_L2_Table_of_new_epimutations)
rownames(tinyRNA_26G_DOWN_L2_Table_of_new_epimutations) <- NULL
tinyRNA_26G_DOWN_L2_Table_of_new_epimutations <- cbind(names,tinyRNA_26G_DOWN_L2_Table_of_new_epimutations)
colnames(tinyRNA_26G_DOWN_L2_Table_of_new_epimutations) <- c("Transitions","Down")

tinyRNA_26G_new_epimut_L2 <- merge(tinyRNA_26G_UP_L2_Table_of_new_epimutations,tinyRNA_26G_DOWN_L2_Table_of_new_epimutations)
tinyRNA_26G_new_epimut_L2 <- tinyRNA_26G_new_epimut_L2 %>%
  # Creating an empty column:
  add_column(Lineage = "L2", .after="Transitions")

tiny_26G_L <- rbind(tinyRNA_26G_new_epimut_L1, tinyRNA_26G_new_epimut_L2)
tiny_26G_L <- tiny_26G_L %>%
  # Creating an empty column:
  add_column(Condition = "Low dose", .after="Transitions")
names <- rownames(tinyRNA_26G_UP_H1_Table_of_new_epimutations)
rownames(tinyRNA_26G_UP_H1_Table_of_new_epimutations) <- NULL
tinyRNA_26G_UP_H1_Table_of_new_epimutations <- cbind(names,tinyRNA_26G_UP_H1_Table_of_new_epimutations)
colnames(tinyRNA_26G_UP_H1_Table_of_new_epimutations) <- c("Transitions","Up")
names <- rownames(tinyRNA_26G_DOWN_H1_Table_of_new_epimutations)
rownames(tinyRNA_26G_DOWN_H1_Table_of_new_epimutations) <- NULL
tinyRNA_26G_DOWN_H1_Table_of_new_epimutations <- cbind(names,tinyRNA_26G_DOWN_H1_Table_of_new_epimutations)
colnames(tinyRNA_26G_DOWN_H1_Table_of_new_epimutations) <- c("Transitions","Down")

tinyRNA_26G_new_epimut_H1 <- merge(tinyRNA_26G_UP_H1_Table_of_new_epimutations,tinyRNA_26G_DOWN_H1_Table_of_new_epimutations)
tinyRNA_26G_new_epimut_H1 <- tinyRNA_26G_new_epimut_H1 %>%
  # Creating an empty column:
  add_column(Lineage = "H1", .after="Transitions")
names <- rownames(tinyRNA_26G_UP_H2_Table_of_new_epimutations)
rownames(tinyRNA_26G_UP_H2_Table_of_new_epimutations) <- NULL
tinyRNA_26G_UP_H2_Table_of_new_epimutations <- cbind(names,tinyRNA_26G_UP_H2_Table_of_new_epimutations)
colnames(tinyRNA_26G_UP_H2_Table_of_new_epimutations) <- c("Transitions","Up")
names <- rownames(tinyRNA_26G_DOWN_H2_Table_of_new_epimutations)
rownames(tinyRNA_26G_DOWN_H2_Table_of_new_epimutations) <- NULL
tinyRNA_26G_DOWN_H2_Table_of_new_epimutations <- cbind(names,tinyRNA_26G_DOWN_H2_Table_of_new_epimutations)
colnames(tinyRNA_26G_DOWN_H2_Table_of_new_epimutations) <- c("Transitions","Down")

tinyRNA_26G_new_epimut_H2 <- merge(tinyRNA_26G_UP_H2_Table_of_new_epimutations,tinyRNA_26G_DOWN_H2_Table_of_new_epimutations)
tinyRNA_26G_new_epimut_H2 <- tinyRNA_26G_new_epimut_H2 %>%
  # Creating an empty column:
  add_column(Lineage = "H2", .after="Transitions")

tiny_26G_H <- rbind(tinyRNA_26G_new_epimut_H1, tinyRNA_26G_new_epimut_H2)
tiny_26G_H <- tiny_26G_H %>%
  # Creating an empty column:
  add_column(Condition = "High dose", .after="Transitions")

tiny_26G_newepi<-rbind(tiny_26G_C,tiny_26G_L,tiny_26G_H)
tiny_26G_newepi$Up <- as.numeric(tiny_26G_newepi$Up)
tiny_26G_newepi$Down <- as.numeric(tiny_26G_newepi$Down)
tiny_26G_newepi$Total=rowSums(cbind(tiny_26G_newepi$Up,tiny_26G_newepi$Down),na.rm=FALSE)

#Test for statistical significance
ggdensity(tiny_26G_newepi$Total, 
          main = "Density plot of epimutations_number",
          xlab = "Total number of epimutations_number")
shapiro.test(tiny_26G_newepi$Total)
kruskal.test(Total ~ Condition, data = tiny_26G_newepi)
dunnTest(Total ~ Condition, data = tiny_26G_newepi)

#Fig.4.G
tiny_26G_newepi$Condition<-fct_relevel(tiny_26G_newepi$Condition,c("Control","Low dose","High dose"))
tiny_26G_allLineage_rate <- ggplot(tiny_26G_newepi, aes(x=Condition, y=Total, color = Condition,group=Condition)) + 
  geom_boxplot(fatten = 1, lwd = 1, width=0.5)+
  scale_color_manual(values=c("cornflowerblue", "darkgreen", "red"))+
  labs(y = "Number of new epimutations", x = "\nCondition")+
  geom_dotplot(binaxis='y', binwidth = 0.25, stackdir='center',  dotsize=15, binpositions = "all", stackgroups = TRUE)+
  theme_bw()+
  theme_linedraw()+
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(color = "black", size = 25, angle = 0, face = "plain"),
        axis.text.y = element_text(color = "black", size = 25, angle = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 20, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, face = "bold"))+
  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))+  
  scale_fill_brewer(palette="Pastel2")+
  ggtitle(paste(""))
tiny_26G_allLineage_rate

#Sup.Fig.5.G
Generation<-c("10","12","14","16","18","2","20","4","6","8",
              "10","12","14","16","18","2","20","4","6","8",
              "10","12","14","16","18","2","20","4","6","8",
              "10","12","14","16","18","2","20","4","6","8",
              "10","12","14","16","18","2","20","4","6","8",
              "10","12","14","16","18","2","20","4","6","8")
tiny_26G_newepi_detailed<-cbind(tiny_26G_newepi,Generation)
tiny_26G_newepi_detailed$Generation <- fct_relevel(tiny_26G_newepi_detailed$Generation, c("2", "4", "6","8","10","12","14","16","18","20"))
tiny_26G_newepi_detailed$Condition <- fct_relevel(tiny_26G_newepi_detailed$Condition, c("Control", "Low dose", "High dose"))
tiny_26G_newepi_detailed_plot <- ggplot(tiny_26G_newepi_detailed, aes(x=Generation , y=Total, color = Condition)) + 
  geom_boxplot(fatten = 1, lwd = 1, width=0.5)+
  scale_color_manual(values=c("cornflowerblue", "darkgreen", "red"))+
  labs(y = "Number of new epimutations", x = "\nGeneration")+
  geom_jitter(aes(shape=Lineage, size=8),
              position=position_jitter(width = 0.2,
                                       height = 0.2))+
  theme_bw()+
  theme(axis.text.x = element_text(color = "black", size = 20, angle = 0, face = "plain"),
        axis.text.y = element_text(color = "black", size = 20, angle = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 16, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 16, angle = 90, face = "bold"))+
  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))+  
  scale_fill_brewer(palette="Pastel2")+
  ggtitle(paste(""))
tiny_26G_newepi_detailed_plot 

#Get raw data
write.xlsx(tiny_26G_newepi_detailed,"Data_Fig_4_G_&_Sup_Fig_5_G.xlsx")
#--------------------------
#####Fig.4.H & Sup.Fig.6.H - 26G-RNAs epimutations duration####
tiny_26G_epimut<-tiny_all_epimut %>% filter(grepl('26G', gene_name))
#Fig.4.H
CoxMod<-coxph(Surv(length,complete)~Condition,data=tiny_26G_epimut)
ggforest(CoxMod, data=tiny_26G_epimut)

ggsurvplot(survfit(CoxMod,data=tiny_26G_epimut), palette = "#2E9FDF",cof.int=TRUE,
           ggtheme = theme_minimal())
condition_df <- with(tiny_26G_epimut,
                     data.frame(Condition = c("Control","Low dose","High dose")
                     )
)
condition_df
fit <- survfit(CoxMod, newdata = condition_df)
fit <- survfit(Surv(length, complete) ~ Condition, data = tiny_26G_epimut)
ggsurvplot <- ggsurvplot(fit,  tiny_26G_epimut, censor = T, break.time.by= 1, pval=TRUE,pval.coord = c(0, 0.03),
                         font.main = c(16, "bold"),
                         palette = c("cornflowerblue", "darkgreen","red"),
                         font.x = c(20,"bold"),
                         font.y = c(20,"bold"),
                         font.legend = 20,
                         legend.labs=c("Control","Low dose","High dose"),
                         font.tickslab = 18)+
  # surv.median.line = "hv")+
  ggtitle("Changes survival all conditions")+
  xlab("Time (generations)")
ggsurvplot$plot +theme(plot.title = element_text(hjust = 0.5))

#Sup.Fig.5.H
CoxMod<-coxph(Surv(length,complete)~Lineage,data=tiny_26G_epimut)
ggforest(model=CoxMod,data=tiny_26G_epimut,fontsize = 0.8, noDigits = 2)
fit_2 <- survfit(Surv(length, complete) ~ Lineage, data = tiny_26G_epimut)
ggsurvplot <- ggsurvplot(fit_2,  tiny_26G_epimut, censor = T, break.time.by= 1, pval=TRUE,pval.coord = c(0, 0.03),                         font.main = c(16, "bold"),
                         font.x = c(20,"bold"),
                         palette = c("#3399FF","blue","red","2E9FDF", "green","#006600"),
                         font.y = c(20,"bold"),
                         font.legend = 20,
                         legend.labs=c("C1","C2","H1","H2","L1","L2"),
                         font.tickslab = 18)+
  ggtitle("Changes survival all conditions")+
  xlab("Time (generations)")
ggsurvplot$plot +theme(plot.title = element_text(hjust = 0.5)) 

#Get the raw data
write.xlsx(tiny_26G_epimut,"Data_Fig_4_H_&_Sup_Fig_5_H.xlsx")
#--------------------------
#####Fig.5.B & C - work on tRNAs-derived fragments####
#getcounts
AllFiles<-list.files(pattern="tRNA.bed")
A0<-read.table(AllFiles[1], sep="\t")
getCounts<-function(file_in){
  file_in[which(file_in[,16]==0),13]<-paste("NA_count=0")
  file_in<-file_in[-which(duplicated(file_in[,13])),]
  Input<-read.table(text=file_in[,13], sep="=")
  Input<-as.numeric(Input[,2])
  Output<-cbind(file_in,Input)
  return(Output)
}

SumStartPos<-function(file_in,rangeStart,rangeEnd){
  rangeNums<-seq(rangeStart,rangeEnd,length=rangeEnd-rangeStart+1)
  if(!(ncol(file_in)==17)){stop("no counts detected")}
  else{
    outputN<-c()
    outputS<-c()
    
    for(i in 1:length(rangeNums)){
      sel<-which(file_in[,11]-file_in[,4]==rangeNums[i]&file_in[,7]=="+"&file_in[,15]=="+")
      outputN<-c(outputN, length(sel))
      if(length(sel)>0){
        outputS<-c(outputS, sum(file_in[sel,17]))
        }else{outputS<-c(outputS, 0)}
      }
  
    outputTot<-data.frame(rangeNums,outputN,outputS)
    colnames(outputTot)<-c("startDist","totalSeqs","totalReads")
    return(outputTot)
  }
}

A0Counts<-getCounts(A0)
A0StartsPlus<-SumStartPos(A0Counts[which(A0Counts[,7]=="+"&A0Counts[,15]=="+"),],0,70)

SumStartNeg<-function(file_in,rangeStart,rangeEnd){
  rangeNums<-seq(rangeStart,rangeEnd,length=rangeEnd-rangeStart+1)
  if(!(ncol(file_in)==17)){stop("no counts detected")}
  else{
    outputN<-c()
    outputS<-c()
    
    for(i in 1:length(rangeNums)){
      sel<-which(file_in[,5]-file_in[,12]==rangeNums[i]&file_in[,7]=="-"&file_in[,15]=="-")
      outputN<-c(outputN, length(sel))
      if(length(sel)>0){
        outputS<-c(outputS, sum(file_in[sel,17]))
      }else{outputS<-c(outputS, 0)}
   }
    outputTot<-data.frame(rangeNums,outputN,outputS)
    colnames(outputTot)<-c("startDist","totalSeqs","totalReads")
    return(outputTot)
  }
}
A0StartsNeg<-SumStartPos(A0Counts[which(A0Counts[,7]=="-"&A0Counts[,15]=="-"),],0,70)
StartsPosAll<-list()
StartsNegAll<-list()
for(i in 1:length(AllFiles)){
  Temp<-read.table(AllFiles[i], sep="\t")
  TempCounts<-getCounts(Temp)
  StartsPosAll[[i]]<-SumStartPos(TempCounts[which(TempCounts[,7]=="+"&TempCounts[,15]=="+"),],0,70)
  StartsNegAll[[i]]<-SumStartNeg(TempCounts[which(TempCounts[,7]=="-"&TempCounts[,15]=="-"),],0,70)
  
}

load("Allfiles_SF.Rdata")
     NormPosAll<-list()
     for(i in 1:length(StartsPosAll)){
       NormPosAll[[i]]<-StartsPosAll[[i]][,3]/AllFilesSF[i,2]
     }
     NormNegAll<-list()
     for(i in 1:length(StartsNegAll)){
       NormNegAll[[i]]<-StartsNegAll[[i]][,3]/AllFilesSF[i,2]
     }
     
     names(NormPosAll)<-AllFilesSF[,1]

     ColTabCont<-grep("^[AB]", names(NormPosAll))
     ColTabLow<-grep("^[CD]", names(NormPosAll))
     ColTabHigh<-grep("^[EF]", names(NormPosAll))
     
     NormContPlus<-c()
     for(i in 1:length(ColTabCont)){
       NormContPlus<-cbind(NormContPlus,NormPosAll[[ColTabCont[i]]])
     }
     
     NormLowPlus<-c()
     
     for(i in 1:length(ColTabLow)){
       NormLowPlus<-cbind(NormLowPlus, NormPosAll[[ColTabLow[i]]])
       
     }
     NormHighPlus<-c()
     for(i in 1:length(ColTabHigh)){
       NormHighPlus<-cbind(NormHighPlus, NormPosAll[[ColTabHigh[i]]])
       
     }
     
     NormContNeg<-c()
     for(i in 1:length(ColTabCont)){
       NormContNeg<-cbind(NormContNeg, NormNegAll[[ColTabCont[i]]])
       
     }
     
     NormLowNeg<-c()
     for(i in 1:length(ColTabLow)){
       
       NormLowNeg<-cbind(NormLowNeg, NormNegAll[[ColTabLow[i]]])
     }
     
     NormHighNeg<-c()
     for(i in 1:length(ColTabHigh)){
       
       NormHighNeg<-cbind(NormHighNeg,NormNegAll[[ColTabHigh[i]]])
     }

#Fig.5.B          
plot(apply(NormHighPlus,1,median),col="cornflowerblue",lwd=2,type="l",ylab="total normalized reads",xlab="5_end relative  to tRNA start")
points(apply(NormContPlus,1,median),col="red", lwd=2, type="l")
points(apply(NormLowPlus,1,median),col="darkgreen",lwd=2, type="l")
legend("topright", c("median control","median low Cisplatin","median high Cisplatin"),col=c("cornflowerblue","darkgreen","red"),lwd=2,lty=1)
dev.copy(pdf, "tRNA_cleavage_points.pdf")
dev.off()

#Fig.5.C
vioplot(NormContPlus[41,],NormLowPlus[41,],NormHighPlus[41,],col=c("cornflowerblue","darkgreen","red"),names=c("control","low dose","high dose"), ylab="total normalized tRA 3' ends")

#Get raw data
NormHighPlus<-as.data.frame(NormHighPlus)
NormContPlus<-as.data.frame(NormContPlus)
NormLowPlus<-as.data.frame(NormLowPlus)
write.xlsx(NormHighPlus,"Data_Fig_5_B_&_C_HD.xlsx")
write.xlsx(NormContPlus,"Data_Fig_5_B_&_C_Ctr.xlsx")
write.xlsx(NormLowPlus,"Data_Fig_5_B_&_C_LD.xlsx")


#Fig.4.E
tRNA_volcanoData<-tRNA_
plot(log2(tRNA_volcanoData[,2])-log2(tRNA_volcanoData[,1]),-1*log10((tRNA_volcanoData[,3])))
plot(log2(tRNA_volcanoData[,2])-log2(tRNA_volcanoData[,1]),-1*log10((tRNA_volcanoData[,3])),ylab="pvalue",xlab="log2(high_dose/control)",pch=18,col="red",xlim=c(-1.5,1.5))
abline(h=1.3)
text(-0.5,1.35,"p<0.05")

save(tRNA_normAbun,file="tRNA_normalizedCounts_filteredAbundant.Rdata")
row.names(tRNA_volcanoData)<-row.names(tRNA_normAbun)
colnames(tRNA_volcanoData)<-c("median_cont","median_high","pvalue")
write.table(tRNA_volcanoData, "tRNA_volcanodata.txt", sep="\t")
tRNA_volcanoData[which(tRNA_volcanoData[,3]<0.05),]
text(1.3,2,"GLN_TTG")
text(1.33,1.55,"VAL_AAC")
text(1.3,1.4,"SER_CGA")
dev.copy(pdf, "tRNA_volcanoplot_highvlow.pdf")
dev.off()
#--------------------------
#####Sup.Fig.8.B & C - work on tRNAs-derived fragments####
#Sup.Fig.7.B
ColTabCont_C1<-grep("^[A]", names(NormPosAll))
ColTabCont_C2<-grep("^[B]", names(NormPosAll))
ColTabLow_L1<-grep("^[C]", names(NormPosAll))
ColTabLow_L2<-grep("^[D]", names(NormPosAll))
ColTabHigh_H1<-grep("^[E]", names(NormPosAll))
ColTabHigh_H2<-grep("^[F]", names(NormPosAll))

NormContPlus_C1<-c()
for(i in 1:length(ColTabCont_C1)){
  NormContPlus_C1<-cbind(NormContPlus_C1,NormPosAll[[ColTabCont_C1[i]]])
}

NormContPlus_C2<-c()
for(i in 1:length(ColTabCont_C2)){
  NormContPlus_C2<-cbind(NormContPlus_C2,NormPosAll[[ColTabCont_C2[i]]])
}

NormLowPlus_L1<-c()
for(i in 1:length(ColTabLow_L1)){
  NormLowPlus_L1<-cbind(NormLowPlus_L1,NormPosAll[[ColTabLow_L1[i]]])
}

NormLowPlus_L2<-c()
for(i in 1:length(ColTabLow_L2)){
  NormLowPlus_L2<-cbind(NormLowPlus_L2,NormPosAll[[ColTabLow_L2[i]]])
}

NormHighPlus_H1<-c()
for(i in 1:length(ColTabHigh_H1)){
  NormHighPlus_H1<-cbind(NormHighPlus_H1,NormPosAll[[ColTabHigh_H1[i]]])
}

NormHighPlus_H2<-c()
for(i in 1:length(ColTabHigh_H2)){
  NormHighPlus_H2<-cbind(NormHighPlus_H2,NormPosAll[[ColTabHigh_H2[i]]])
}

NormContNeg_C1<-c()
for(i in 1:length(ColTabCont_C1)){
  NormContNeg_C1<-cbind(NormContNeg_C1,NormPosAll[[ColTabCont_C1[i]]])
}

NormContNeg_C2<-c()
for(i in 1:length(ColTabCont_C2)){
  NormContNeg_C2<-cbind(NormContNeg_C2,NormPosAll[[ColTabCont_C2[i]]])
}

NormLowNeg_L1<-c()
for(i in 1:length(ColTabLow_L1)){
  NormLowNeg_L1<-cbind(NormLowNeg_L1,NormPosAll[[ColTabLow_L1[i]]])
}

NormLowNeg_L2<-c()
for(i in 1:length(ColTabLow_L2)){
  NormLowNeg_L2<-cbind(NormLowNeg_L2,NormPosAll[[ColTabLow_L2[i]]])
}

NormHighNeg_H1<-c()
for(i in 1:length(ColTabHigh_H1)){
  NormHighNeg_H1<-cbind(NormHighNeg_H1,NormPosAll[[ColTabHigh_H1[i]]])
}

NormHighNeg_H2<-c()
for(i in 1:length(ColTabHigh_H2)){
  NormHighNeg_H2<-cbind(NormHighNeg_H2,NormPosAll[[ColTabHigh_H2[i]]])
}

#Sup.Fig.7.B          
plot(apply(NormHighPlus_H2,1,median),col="darkred",lwd=2,type="l",ylab="total normalized reads",xlab="5_end relative  to tRNA start")
points(apply(NormContPlus_C2,1,median),col="blue", lwd=2, type="l")
points(apply(NormLowPlus_L1,1,median),col="green",lwd=2, type="l")
points(apply(NormLowPlus_L2,1,median),col="darkgreen",lwd=2, type="l")
points(apply(NormHighPlus_H1,1,median),col="red",lwd=2, type="l")
points(apply(NormContPlus_C1,1,median),col="cornflowerblue",lwd=2, type="l")
legend("topright", c("median C1","median C2","median L1","median L2", "median H1", "median H2"),col=c("cornflowerblue","blue","green","darkgreen","red","darkred"),lwd=2,lty=1)

#Sup.Fig.7.C
vioplot(NormContPlus_C1[41,],NormContPlus_C2[41,],NormLowPlus_L1[41,],NormLowPlus_L2[41,],NormHighPlus_H1[41,],NormHighPlus_H2[41,], col=c("cornflowerblue","blue","green","darkgreen","red","darkred"),names=c("C1","C2","L1","L2","H1","H2"), ylab="total normalized tRA 3' ends")

#Get raw data
NormContPlus_C1<-as.data.frame(NormContPlus_C1)
NormContPlus_C2<-as.data.frame(NormContPlus_C2)
NormLowPlus_L1<-as.data.frame(NormLowPlus_L1)
NormLowPlus_L2<-as.data.frame(NormLowPlus_L2)
NormHighPlus_H1<-as.data.frame(NormHighPlus_H1)
NormHighPlus_H2<-as.data.frame(NormHighPlus_H2)
write.xlsx(NormContPlus_C1,"NormContPlus_C1.xlsx")
write.xlsx(NormContPlus_C2,"NormContPlus_C2.xlsx")
write.xlsx(NormLowPlus_L1,"NormLowPlus_L1.xlsx")
write.xlsx(NormLowPlus_L2,"NormLowPlus_L2.xlsx")
write.xlsx(NormHighPlus_H1,"NormHighPlus_H1.xlsx")
write.xlsx(NormHighPlus_H2,"NormHighPlus_H2.xlsx")
#--------------------------
#####Fig.5.A & E - work on tRNAs-derived fragments#####
#Ranalysis1.Rdata
#get unnormalized counts
load("annotatedtRNAIDs.Rdata")
AllFiles<-list.files(pattern=".bed")

CountsAll<-matrix(0,nrow=nrow(tRNAids),ncol=length(AllFiles))

for(i in 1:length(AllFiles)){
  Tmp<-read.table(AllFiles[i], sep="\t", stringsAsFactors=F)
  Tmp<-getCounts(Tmp)
  for(j in 1:nrow(CountsAll)){
    k<-which(Tmp[,9]==tRNAids[j,1])
    if(length(k>0)){
      CountsAll[j,i]<-sum(Tmp[k,ncol(Tmp)])}
    
  }
}
colnames(CountsAll)<-c("A_0","A_10","A_12","A_14","A_16","A_18","A_2","A_20","A_4","A_6","A_8",
                        "B_0","B_10","B_12","B_14","B_16","B_18","B_2","B_20","B_4","B_8",
                          "C_0","C_12","C_14","C_16","C_18","C_2","C_20","C_4","C_6","C_8",
                          "D_0","D_10","D_12","D_14","D_16","D_18","D_2","D_20","D_4","D_8",
                          "E_0","E_10","E_12","E_14","E_16","E_18","E_2","E_20","E_4","E_8",
                          "F_0","F_10","F_12","F_14","F_16","F_18","F_2","F_20","F_4","F_6","F_8")
ColTabCont<-grep("^[AB]", colnames(CountsAll))
ColTabLow<-grep("^[CD]", colnames(CountsAll))
ColTabHigh<-grep("^[EF]", colnames(CountsAll))

CountsAllSens<-matrix(0,nrow=nrow(tRNAids),ncol=length(AllFiles))

for(i in 1:length(AllFiles)){
  Tmp<-read.table(AllFiles[i], sep="\t", stringsAsFactors=F)
  Tmp<-Tmp[which(Tmp[,7]==Tmp[,15]),]
  Tmp<-getCounts(Tmp)
  for(j in 1:nrow(CountsAllSens)){
    k<-which(Tmp[,9]==tRNAids[j,1])
    if(length(k>0)){
      CountsAllSens[j,i]<-sum(Tmp[k,ncol(Tmp)])}
    
  }
}
unnorm<-read.table("sRNA_align.csv", sep=",", header=T)
unnorm<-read.xlsx("sRNA_align.xlsx")

unnormDESeq<-unnorm[,4:ncol(unnorm)]

unnormDESeq<-round(unnormDESeq, digits=0)

colInfo<-colnames(unnormDESeq)
colInfo<-cbind(colInfo, rep("Control", length(colInfo)))
colInfo[grep("[CD]", colInfo[,1]),2]<-"Low"
colInfo[grep("[EF]", colInfo[,1]),2]<-"High"
colInfo<-as.data.frame(colInfo)
colInfo[,2]<-as.factor(colInfo[,2])
colnames(colInfo)<-c("ID", "CisplatinDose")
DESeq_allsRNA<-DESeqDataSetFromMatrix(unnormDESeq,colData=colInfo, design=~CisplatinDose)
DESeq_allsRNA<-DESeq(DESeq_allsRNA)
sRNA_sizeFactor<-sizeFactors(DESeq_allsRNA)

namesFiles<-paste(namesFiles[,1],namesFiles[,2], sep="_")
namesFiles<-c("A_0","A_2","A_4","A_6","A_8","A_10","A_12","A_14","A_16","A_18","A_20",
              "B_0","B_2","B_4","B_8","B_10","B_12","B_14","B_16","B_18","B_20",
              "C_0","C_2","C_4","C_6","C_8","C_12","C_14","C_16","C_18","C_20",
              "D_0","D_2","D_4","D_8","D_10","D_12","D_14","D_16","D_18","D_20",
              "E_0","E_2","E_4","E_8","E_10","E_12","E_14","E_16","E_18","E_20",
              "F_0","F_2","F_4","F_6","F_8","F_10","F_12","F_14","F_16","F_18","F_20")
FileMatches<-match(namesFiles,names(sRNA_sizeFactor))
namesFiles<-data.frame(namesFiles,sRNA_sizeFactor[FileMatches])
colnames(namesFiles)<-c("ID","DESeq_sizeFactor")

tRNA_norm<-CountsAllSens
for(i in 1:ncol(CountsAllSens)){
  tRNA_norm[,i]<-CountsAllSens[,i]*namesFiles[i,2]
  
}

tRNA_normAbun<-tRNA_norm[which(apply(tRNA_norm,1,max)>500),]
#remove unannotated tRNAs
tRNA_normAbun<-tRNA_normAbun[-grep("NA", row.names(tRNA_normAbun)),]

#Fig.4.A
boxplot(apply(tRNA_normAbun[,ColTabCont],2,sum),apply(tRNA_normAbun[,ColTabLow],2,sum),apply(tRNA_normAbun[,ColTabHigh],2,sum)
        ,col=c("cornflowerblue","darkgreen","red"),names=c("control","low cisplatin","high cisplatin"),ylab="total norm tRNA reads")
dev.copy(pdf, "DESeq normalized TEcounts Total.pdf")
dev.off()

wilcox.test(apply(tRNA_normAbun[,ColTabCont],2,sum),apply(tRNA_normAbun[,ColTabHigh],2,sum))
#p=0.11 not significantly higher. 

var.test(apply(tRNA_normAbun[,ColTabCont],2,sum),apply(tRNA_normAbun[,ColTabHigh],2,sum))
#p=0.039 significantly higher in high dose

tRNA_volcanoData<-matrix(0, ncol=3,nrow=nrow(tRNA_normAbun))
for(i in 1:nrow(tRNA_normAbun)){
  
  WT<-wilcox.test(tRNA_normAbun[i,ColTabCont],tRNA_normAbun[i,ColTabHigh])
  tRNA_volcanoData[i,1]<-median(tRNA_normAbun[i,ColTabCont])
  tRNA_volcanoData[i,2]<-median(tRNA_normAbun[i,ColTabHigh])
  tRNA_volcanoData[i,3]<-WT$p.value
}
row.names(tRNA_volcanoData)<-row.names(tRNA_normAbun)
colnames(tRNA_volcanoData)<-c("median_cont","median_high","pvalue")

#--------------------------
#####Sup.Fig.8.A#####
#Ranalysis1.Rdata
#get unnormalized counts
load("annotatedtRNAIDs.Rdata")
AllFiles<-list.files(pattern=".bed")

CountsAll<-matrix(0,nrow=nrow(tRNAids),ncol=length(AllFiles))

for(i in 1:length(AllFiles)){
  Tmp<-read.table(AllFiles[i], sep="\t", stringsAsFactors=F)
  Tmp<-getCounts(Tmp)
  for(j in 1:nrow(CountsAll)){
    k<-which(Tmp[,9]==tRNAids[j,1])
    if(length(k>0)){
      CountsAll[j,i]<-sum(Tmp[k,ncol(Tmp)])}
    
  }
}
colnames(CountsAll)<-c("A_0","A_10","A_12","A_14","A_16","A_18","A_2","A_20","A_4","A_6","A_8",
                       "B_0","B_10","B_12","B_14","B_16","B_18","B_2","B_20","B_4","B_8",
                       "C_0","C_12","C_14","C_16","C_18","C_2","C_20","C_4","C_6","C_8",
                       "D_0","D_10","D_12","D_14","D_16","D_18","D_2","D_20","D_4","D_8",
                       "E_0","E_10","E_12","E_14","E_16","E_18","E_2","E_20","E_4","E_8",
                       "F_0","F_10","F_12","F_14","F_16","F_18","F_2","F_20","F_4","F_6","F_8")
ColTabCont_C1<-grep("^[A]", colnames(CountsAll))
ColTabCont_C2<-grep("^[B]", colnames(CountsAll))
ColTabLow_L1<-grep("^[C]", colnames(CountsAll))
ColTabLow_L2<-grep("^[D]", colnames(CountsAll))
ColTabHigh_H1<-grep("^[E]", colnames(CountsAll))
ColTabHigh_H2<-grep("^[F]", colnames(CountsAll))

#Sup.Fig.7.A
boxplot(apply(tRNA_normAbun[,ColTabCont_C1],2,sum),apply(tRNA_normAbun[,ColTabCont_C2],2,sum),apply(tRNA_normAbun[,ColTabLow_L1],2,sum),apply(tRNA_normAbun[,ColTabLow_L2],2,sum),apply(tRNA_normAbun[,ColTabHigh_H1],2,sum),apply(tRNA_normAbun[,ColTabHigh_H2],2,sum)
        ,col=c("cornflowerblue","blue","green","darkgreen","red","darkred"),names=c("C1","C2","L1","L2","H1","H2"),ylab="total norm tRNA reads")

#Get the raw data
write.xlsx(tRNA_normAbun,"Data_Fig_5_A &_E_&_Sup_Fig_7_A.xlsx")
#--------------------------
#####Fig.5.F - Association tRNAs fragments & AGOs####
#Ago IP analysis
GetCount<-function(fileIn){M<-read.table(text=fileIn[,10], sep="-");M<-as.numeric(M[,2]);tRNACount<-M;return(tRNACount)}
AllFiles<-list.files(pattern="tRNA.bed")
head(AllFiles)
AllSums<-c()
AlltRNA<-list()
for(i in 1:length(AllFiles)){AlltRNA[[i]]<-read.table(AllFiles[i], sep="\t", stringsAsFactors=F);AlltRNA[[i]]<-AlltRNA[[i]][which(AlltRNA[[i]][,13]>18&AlltRNA[[i]][,6]==AlltRNA[[i]][,12]),];AlltRNA[[i]]<-cbind(AlltRNA[[i]],GetCount(AlltRNA[[i]]));AllSums<-c(AllSums, sum(AlltRNA[[i]][,14]))}
head(AllSums)

#now we need the normalization factor for each file
NontRNA<-list.files(pattern=".fasta.bed$")
GetCountsAll<-function(file_in){temp<-read.table(file_in);temp<-read.table(text=temp[,4], sep="-");out<-as.numeric(temp[,2]); return(sum(out))}
Allreadcount<-c()
for(i in 1:length(NontRNA)){Allreadcount<-c(Allreadcount,GetCountsAll(NontRNA[i]))}

AllSumsNorm<-AllSums/Allreadcount

#fileinfo
fileInfo<-read.table("filelist.csv", sep=",", stringsAsFactors=F)

R1Files<-grep("_R1", fileInfo[,3])
Inputs<-grep("Input", fileInfo[,4])
Reals<-grep("Real", fileInfo[,4])
R1Reals<-intersect(R1Files,Reals)
R1Inputs<-intersect(R1Files,Inputs)
R1IPSums<-data.frame(fileInfo[R1Reals,],AllSumsNorm[R1Reals])
R1InputSums<-data.frame(fileInfo[R1Inputs,],AllSumsNorm[R1Inputs])

R2Files<-grep("_R2", fileInfo[,3])
R2Reals<-intersect(R2Files,Reals)
R2Inputs<-intersect(R2Files,Inputs)
R2IPSums<-data.frame(fileInfo[R2Reals,],AllSumsNorm[R2Reals])
R2InputSums<-data.frame(fileInfo[R2Inputs,],AllSumsNorm[R2Inputs])
r2M<-match(R2IPSums[,2],R2InputSums[,2])
barplot(R2IPSums[,5]/R2InputSums[r2M,5], names=R2IPSums[,2], las=2)

#--------------------------
#####Fig.5.G - tRNAs-derived fragments epimutations####
load("/Users/manonfallet/Downloads/tRNA_normalizedCounts_filteredAbundant.Rdata")
#Data preparation
tRNA_normAbun
colnames(tRNA_normAbun)<-c("C1_0","C1_10","C1_12","C1_14","C1_16","C1_18","C1_2","C1_20","C1_4","C1_6","C1_8",
                           "C2_0","C2_10","C2_12","C2_14","C2_16","C2_18","C2_2","C2_20","C2_4","C2_8",
                           "L1_0","L1_12","L1_14","L1_16","L1_18","L1_2","L1_20","L1_4","L1_6","L1_8",
                           "L2_0","L2_10","L2_12","L2_14","L2_16","L2_18","L2_2","L2_20","L2_4","L2_8",
                           "H1_0","H1_10","H1_12","H1_14","H1_16","H1_18","H1_2","H1_20","H1_4","H1_8",
                           "H2_0","H2_10","H2_12","H2_14","H2_16","H2_18","H2_2","H2_20","H2_4","H2_6","H2_8")
#Epimutations detection
#Create data frame per lineage
tRNA_normAbun<-as.data.frame(tRNA_normAbun)
tRNA_normAbun_Lineage_C1 <- cbind(tRNA_normAbun$C1_0,tRNA_normAbun$C1_2,tRNA_normAbun$C1_4,tRNA_normAbun$C1_6,tRNA_normAbun$C1_8,tRNA_normAbun$C1_10,tRNA_normAbun$C1_12,tRNA_normAbun$C1_14,tRNA_normAbun$C1_16,tRNA_normAbun$C1_18,tRNA_normAbun$C1_20)
colnames(tRNA_normAbun_Lineage_C1) <- c("C1_0", "C1_2", "C1_4", "C1_6", "C1_8", "C1_10","C1_12","C1_14", "C1_16", "C1_18","C1_20")
tRNA_normAbun_Lineage_C2 <- cbind(tRNA_normAbun$C2_0,tRNA_normAbun$C2_2,tRNA_normAbun$C2_4,tRNA_normAbun$C2_8,tRNA_normAbun$C2_10,tRNA_normAbun$C2_12,tRNA_normAbun$C2_14,tRNA_normAbun$C2_16,tRNA_normAbun$C2_18,tRNA_normAbun$C2_20)
colnames(tRNA_normAbun_Lineage_C2) <- c("C2_0", "C2_2", "C2_4", "C2_8", "C2_10","C2_12","C2_14", "C2_16", "C2_18","C2_20")
tRNA_normAbun_Lineage_L1 <- cbind(tRNA_normAbun$L1_0,tRNA_normAbun$L1_2,tRNA_normAbun$L1_4,tRNA_normAbun$L1_6,tRNA_normAbun$L1_8,tRNA_normAbun$L1_12,tRNA_normAbun$L1_14,tRNA_normAbun$L1_16,tRNA_normAbun$L1_18,tRNA_normAbun$L1_20)
colnames(tRNA_normAbun_Lineage_L1) <- c("L1_0", "L1_2", "L1_4", "L1_6", "L1_8", "L1_12", "L1_14", "L1_16", "L1_18","L1_20")
tRNA_normAbun_Lineage_L2 <- cbind(tRNA_normAbun$L2_0,tRNA_normAbun$L2_2,tRNA_normAbun$L2_4,tRNA_normAbun$L2_8,tRNA_normAbun$L2_10,tRNA_normAbun$L2_12,tRNA_normAbun$L2_14,tRNA_normAbun$L2_16,tRNA_normAbun$L2_18,tRNA_normAbun$L2_20)
colnames(tRNA_normAbun_Lineage_L2) <- c("L2_0", "L2_2", "L2_4", "L2_8", "L2_10", "L2_12", "L2_14", "L2_16", "L2_18","L2_20")
tRNA_normAbun_Lineage_H1 <- cbind(tRNA_normAbun$H1_0,tRNA_normAbun$H1_2,tRNA_normAbun$H1_4,tRNA_normAbun$H1_8,tRNA_normAbun$H1_10,tRNA_normAbun$H1_12,tRNA_normAbun$H1_14,tRNA_normAbun$H1_16,tRNA_normAbun$H1_18,tRNA_normAbun$H1_20)
colnames(tRNA_normAbun_Lineage_H1) <- c("H1_0", "H1_2", "H1_4", "H1_8", "H1_10", "H1_12", "H1_14", "H1_16", "H1_18", "H1_20")
tRNA_normAbun_Lineage_H2 <- cbind(tRNA_normAbun$H2_0,tRNA_normAbun$H2_2,tRNA_normAbun$H2_4,tRNA_normAbun$H2_6,tRNA_normAbun$H2_8,tRNA_normAbun$H2_10,tRNA_normAbun$H2_12,tRNA_normAbun$H2_14,tRNA_normAbun$H2_16,tRNA_normAbun$H2_18,tRNA_normAbun$H2_20)
colnames(tRNA_normAbun_Lineage_H2) <- c("H2_0", "H2_2", "H2_4", "H2_6", "H2_8", "H2_10", "H2_12", "H2_14", "H2_16", "H2_18","H2_20")

#C1
#Identify epimutations in tinyRNA using linear model in each lineage 
tRNA_Lineage_C1_2_to_20<-tRNA_normAbun_Lineage_C1[,2:11]
tRNA_z_scores_table_C1<-matrix(0,ncol=10, nrow = nrow(tRNA_normAbun_Lineage_C1))

for (i in 1:ncol(tRNA_Lineage_C1_2_to_20)) {
  tRNA_EpimutC1<-
    lm(tRNA_normAbun_Lineage_C1[,1]~tRNA_Lineage_C1_2_to_20[,i])
  tRNA_ResidC1<-tRNA_EpimutC1$residuals
  tRNA_z_scores_table_C1[,i]<-(tRNA_ResidC1-mean(tRNA_ResidC1))/sd(tRNA_ResidC1)
}

#keep only relevant epimutations and binarised data
tRNA_binarised_z_scores_table_C1<-matrix(0,ncol=ncol(tRNA_z_scores_table_C1),nrow=nrow(tRNA_z_scores_table_C1))

for (i in 1:ncol(tRNA_z_scores_table_C1)) {
  for(j in 1:nrow(tRNA_z_scores_table_C1)){
    
    if(tRNA_z_scores_table_C1[j,i]> (-2.25) & tRNA_z_scores_table_C1[j,i]< 2.25) {
      tRNA_binarised_z_scores_table_C1[j,i]<-0
    }
    if(tRNA_z_scores_table_C1[j,i] > 2.25){
      tRNA_binarised_z_scores_table_C1[j,i]<-1
    }
    if(tRNA_z_scores_table_C1[j,i] < (-2.25)){
      tRNA_binarised_z_scores_table_C1[j,i]<- (-1)
    }
  }  
}
tRNA_C1epimutations<-tRNA_binarised_z_scores_table_C1 %>% as.data.frame() 
colnames(tRNA_C1epimutations) <- c('F2','4','6','8','10','12','14','16','18','20')
tRNA_C1epimutations <- tRNA_C1epimutations %>%
  # Creating an empty column:
  add_column(F0 = 0, .before='F2')
tRNA_C1epimutations <- tRNA_C1epimutations %>%
  # Creating a condition column:
  add_column(condition = "C1", .before="F0")
tRNA_C1epimutations <- cbind(rownames(tRNA_normAbun),tRNA_C1epimutations)
colnames(tRNA_C1epimutations) <- c('tRNAs','condition','0','2','4','6','8','10','12','14','16','18','20')
save(tRNA_C1epimutations,file="tRNA_C1epimutations.Rdata")

#C2
#Identify epimutations in tinyRNA using linear model in each lineage 
tRNA_Lineage_C2_2_to_20<-tRNA_normAbun_Lineage_C2[,2:10]
tRNA_z_scores_table_C2<-matrix(0,ncol=9, nrow = nrow(tRNA_normAbun_Lineage_C2))

for (i in 1:ncol(tRNA_Lineage_C2_2_to_20)) {
  tRNA_EpimutC2<-
    lm(tRNA_normAbun_Lineage_C2[,1]~tRNA_Lineage_C2_2_to_20[,i])
  tRNA_ResidC2<-tRNA_EpimutC2$residuals
  tRNA_z_scores_table_C2[,i]<-(tRNA_ResidC2-mean(tRNA_ResidC2))/sd(tRNA_ResidC2)
}

#keep only relevant epimutations and binarised data
tRNA_binarised_z_scores_table_C2<-matrix(0,ncol=ncol(tRNA_z_scores_table_C2),nrow=nrow(tRNA_z_scores_table_C2))

for (i in 1:ncol(tRNA_z_scores_table_C2)) {
  for(j in 1:nrow(tRNA_z_scores_table_C2)){
    
    if(tRNA_z_scores_table_C2[j,i]> (-2.25) & tRNA_z_scores_table_C2[j,i]< 2.25) {
      tRNA_binarised_z_scores_table_C2[j,i]<-0
    }
    if(tRNA_z_scores_table_C2[j,i] > 2.25){
      tRNA_binarised_z_scores_table_C2[j,i]<-1
    }
    if(tRNA_z_scores_table_C2[j,i] < (-2.25)){
      tRNA_binarised_z_scores_table_C2[j,i]<- (-1)
    }
  }  
}
tRNA_C2epimutations<-tRNA_binarised_z_scores_table_C2 %>% as.data.frame() 
colnames(tRNA_C2epimutations) <- c('F2','4','8','10','12','14','16','18','20')
tRNA_C2epimutations <- tRNA_C2epimutations %>%
  # Creating an empty column:
  add_column(F0 = 0, .before='F2')
tRNA_C2epimutations <- tRNA_C2epimutations %>%
  # Creating a condition column:
  add_column(condition = "C2", .before="F0")
tRNA_C2epimutations <- cbind(rownames(tRNA_normAbun),tRNA_C2epimutations)
colnames(tRNA_C2epimutations) <- c('tRNAs','condition','0','2','4','8','10','12','14','16','18','20')

save(tRNA_C2epimutations,file="tRNA_C2epimutations.Rdata")

#L1
#Identify epimutations in tinyRNA using linear model in each lineage 
# Replace all values of 0 with 
tRNA_Lineage_L1_2_to_20<-tRNA_normAbun_Lineage_L1[,2:10]
tRNA_z_scores_table_L1<-matrix(0,ncol=9, nrow = nrow(tRNA_normAbun_Lineage_L1))

for (i in 1:ncol(tRNA_Lineage_L1_2_to_20)) {
  tRNA_EpimutL1<-
    lm(tRNA_normAbun_Lineage_L1[,1]~tRNA_Lineage_L1_2_to_20[,i])
  tRNA_ResidL1<-tRNA_EpimutL1$residuals
  tRNA_z_scores_table_L1[,i]<-(tRNA_ResidL1-mean(tRNA_ResidL1))/sd(tRNA_ResidL1)
}

#keep only relevant epimutations and binarised data
tRNA_binarised_z_scores_table_L1<-matrix(0,ncol=ncol(tRNA_z_scores_table_L1),nrow=nrow(tRNA_z_scores_table_L1))

for (i in 1:ncol(tRNA_z_scores_table_L1)) {
  for(j in 1:nrow(tRNA_z_scores_table_L1)){
    
    if(tRNA_z_scores_table_L1[j,i]> (-2.25) & tRNA_z_scores_table_L1[j,i]< 2.25) {
      tRNA_binarised_z_scores_table_L1[j,i]<-0
    }
    if(tRNA_z_scores_table_L1[j,i] > 2.25){
      tRNA_binarised_z_scores_table_L1[j,i]<-1
    }
    if(tRNA_z_scores_table_L1[j,i] < (-2.25)){
      tRNA_binarised_z_scores_table_L1[j,i]<- (-1)
    }
  }  
}
tRNA_L1epimutations<-tRNA_binarised_z_scores_table_L1 %>% as.data.frame() 
colnames(tRNA_L1epimutations) <- c('F2','4','6','8','12','14','16','18','20')
tRNA_L1epimutations <- tRNA_L1epimutations %>%
  # Creating an empty column:
  add_column(F0 = 0, .before='F2')
tRNA_L1epimutations <- tRNA_L1epimutations %>%
  # Creating a condition column:
  add_column(condition = "L1", .before="F0")
tRNA_L1epimutations <- cbind(rownames(tRNA_normAbun),tRNA_L1epimutations)
colnames(tRNA_L1epimutations) <- c('tRNAs', 'condition','0','2','4','6','8','12','14','16','18','20')
save(tRNA_L1epimutations,file="tRNA_L1epimutations.Rdata")

#L2
#Identify epimutations in tinyRNA using linear model in each lineage 
tRNA_Lineage_L2_2_to_20<-tRNA_normAbun_Lineage_L2[,2:10]
tRNA_z_scores_table_L2<-matrix(0,ncol=9, nrow = nrow(tRNA_normAbun_Lineage_L2))

for (i in 1:ncol(tRNA_Lineage_L2_2_to_20)) {
  tRNA_EpimutL2<-
    lm(tRNA_normAbun_Lineage_L2[,1]~tRNA_Lineage_L2_2_to_20[,i])
  tRNA_ResidL2<-tRNA_EpimutL2$residuals
  tRNA_z_scores_table_L2[,i]<-(tRNA_ResidL2-mean(tRNA_ResidL2))/sd(tRNA_ResidL2)
}

#keep only relevant epimutations and binarised data
tRNA_binarised_z_scores_table_L2<-matrix(0,ncol=ncol(tRNA_z_scores_table_L2),nrow=nrow(tRNA_z_scores_table_L2))

for (i in 1:ncol(tRNA_z_scores_table_L2)) {
  for(j in 1:nrow(tRNA_z_scores_table_L2)){
    
    if(tRNA_z_scores_table_L2[j,i]> (-2.25) & tRNA_z_scores_table_L2[j,i]< 2.25) {
      tRNA_binarised_z_scores_table_L2[j,i]<-0
    }
    if(tRNA_z_scores_table_L2[j,i] > 2.25){
      tRNA_binarised_z_scores_table_L2[j,i]<-1
    }
    if(tRNA_z_scores_table_L2[j,i] < (-2.25)){
      tRNA_binarised_z_scores_table_L2[j,i]<- (-1)
    }
  }  
}
tRNA_L2epimutations<-tRNA_binarised_z_scores_table_L2 %>% as.data.frame() 
colnames(tRNA_L2epimutations) <- c('F2','4','8','10','12','14','16','18','20')
tRNA_L2epimutations <- tRNA_L2epimutations %>%
  # Creating an empty column:
  add_column(F0 = 0, .before='F2')
tRNA_L2epimutations <- tRNA_L2epimutations %>%
  # Creating a condition column:
  add_column(condition = "L2", .before="F0")
tRNA_L2epimutations <- cbind(rownames(tRNA_normAbun),tRNA_L2epimutations)
colnames(tRNA_L2epimutations) <- c('tRNAs','condition','0','2','4','8','10','12','14','16','18','20')
save(tRNA_L2epimutations,file="tRNA_L2epimutations.Rdata")

#H1
#Identify epimutations in tinyRNA using linear model in each lineage 
tRNA_Lineage_H1_2_to_20<-tRNA_normAbun_Lineage_H1[,2:10]
tRNA_z_scores_table_H1<-matrix(0,ncol=9, nrow = nrow(tRNA_normAbun_Lineage_H1))

for (i in 1:ncol(tRNA_Lineage_H1_2_to_20)) {
  tRNA_EpimutH1<-
    lm(tRNA_normAbun_Lineage_H1[,1]~tRNA_Lineage_H1_2_to_20[,i])
  tRNA_ResidH1<-tRNA_EpimutH1$residuals
  tRNA_z_scores_table_H1[,i]<-(tRNA_ResidH1-mean(tRNA_ResidH1))/sd(tRNA_ResidH1)
}

#keep only relevant epimutations and binarised data
tRNA_binarised_z_scores_table_H1<-matrix(0,ncol=ncol(tRNA_z_scores_table_H1),nrow=nrow(tRNA_z_scores_table_H1))

for (i in 1:ncol(tRNA_z_scores_table_H1)) {
  for(j in 1:nrow(tRNA_z_scores_table_H1)){
    
    if(tRNA_z_scores_table_H1[j,i]> (-2.25) & tRNA_z_scores_table_H1[j,i]< 2.25) {
      tRNA_binarised_z_scores_table_H1[j,i]<-0
    }
    if(tRNA_z_scores_table_H1[j,i] > 2.25){
      tRNA_binarised_z_scores_table_H1[j,i]<-1
    }
    if(tRNA_z_scores_table_H1[j,i] < (-2.25)){
      tRNA_binarised_z_scores_table_H1[j,i]<- (-1)
    }
  }  
}
tRNA_H1epimutations<-tRNA_binarised_z_scores_table_H1 %>% as.data.frame() 
colnames(tRNA_H1epimutations) <- c('F2','4','8','10','12','14','16','18','20')
tRNA_H1epimutations <- tRNA_H1epimutations %>%
  # Creating an empty column:
  add_column(F0 = 0, .before='F2')
tRNA_H1epimutations <- tRNA_H1epimutations %>%
  # Creating a condition column:
  add_column(condition = "H1", .before="F0")
tRNA_H1epimutations <- cbind(rownames(tRNA_normAbun),tRNA_H1epimutations)
colnames(tRNA_H1epimutations) <- c('tRNAs','condition','0','2','4','8','10','12','14','16','18','20')
save(tRNA_H1epimutations,file="tRNA_H1epimutations.Rdata")

#H2
#Identify epimutations in tinyRNA using linear model in each lineage 
tRNA_Lineage_H2_2_to_20<-tRNA_normAbun_Lineage_H2[,2:11]
tRNA_z_scores_table_H2<-matrix(0,ncol=10, nrow = nrow(tRNA_normAbun_Lineage_H2))

for (i in 1:ncol(tRNA_Lineage_H2_2_to_20)) {
  tRNA_EpimutH2<-
    lm(tRNA_normAbun_Lineage_H2[,1]~tRNA_Lineage_H2_2_to_20[,i])
  tRNA_ResidH2<-tRNA_EpimutH2$residuals
  tRNA_z_scores_table_H2[,i]<-(tRNA_ResidH2-mean(tRNA_ResidH2))/sd(tRNA_ResidH2)
}

#keep only relevant epimutations and binarised data
tRNA_binarised_z_scores_table_H2<-matrix(0,ncol=ncol(tRNA_z_scores_table_H2),nrow=nrow(tRNA_z_scores_table_H2))

for (i in 1:ncol(tRNA_z_scores_table_H2)) {
  for(j in 1:nrow(tRNA_z_scores_table_H2)){
    
    if(tRNA_z_scores_table_H2[j,i]> (-2.25) & tRNA_z_scores_table_H2[j,i]< 2.25) {
      tRNA_binarised_z_scores_table_H2[j,i]<-0
    }
    if(tRNA_z_scores_table_H2[j,i] > 2.25){
      tRNA_binarised_z_scores_table_H2[j,i]<-1
    }
    if(tRNA_z_scores_table_H2[j,i] < (-2.25)){
      tRNA_binarised_z_scores_table_H2[j,i]<- (-1)
    }
  }  
}
tRNA_H2epimutations<-tRNA_binarised_z_scores_table_H2 %>% as.data.frame() 
colnames(tRNA_H2epimutations) <- c('F2','4','6','8','10','12','14','16','18','20')
tRNA_H2epimutations <- tRNA_H2epimutations %>%
  # Creating an empty column:
  add_column(F0 = 0, .before='F2')
tRNA_H2epimutations <- tRNA_H2epimutations %>%
  # Creating a condition column:
  add_column(condition = "H2", .before="F0")
tRNA_H2epimutations <- cbind(rownames(tRNA_normAbun),tRNA_H2epimutations)
colnames(tRNA_H2epimutations) <- c('tRNAs','condition','0','2','4','6','8','10','12','14','16','18','20')
save(tRNA_H2epimutations,file="tRNA_H2epimutations.Rdata")

#Calculate the new epimutations each generation
#Control1 new epi
#UP

tRNA_UP_output_C1<- c()
row.names(tRNA_C1epimutations)<-tRNA_C1epimutations$tRNAs
tRNA_C1epimutations<-tRNA_C1epimutations[,c(3:13)]
for(i in 1:nrow(tRNA_C1epimutations)){
  tRNA_UP_output_C1 <-rbind(tRNA_UP_output_C1, UP_transition_func(tRNA_C1epimutations[i,], input_name=row.names(tRNA_C1epimutations)[i]))}
colnames(tRNA_UP_output_C1) <- c("0", "2", "4", "6", "8", "10", "12", "14", "16", "18", "20")
row.names(tRNA_UP_output_C1) <- rownames(tRNA_C1epimutations)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tRNA_UP_output_C1[,2])
transitions_at_4 <- sum(tRNA_UP_output_C1[,3])
transitions_at_6 <- sum(tRNA_UP_output_C1[,4])
transitions_at_8 <- sum(tRNA_UP_output_C1[,5])
transitions_at_10 <- sum(tRNA_UP_output_C1[,6])
transitions_at_12 <- sum(tRNA_UP_output_C1[,7])
transitions_at_14 <- sum(tRNA_UP_output_C1[,8])
transitions_at_16 <- sum(tRNA_UP_output_C1[,9])
transitions_at_18 <- sum(tRNA_UP_output_C1[,10])
transitions_at_20 <- sum(tRNA_UP_output_C1[,11])

tRNA_UP_C1_Table_of_new_epimutations <- rbind(transitions_at_2,
                                              transitions_at_4, 
                                              transitions_at_6,
                                              transitions_at_8, 
                                              transitions_at_10, 
                                              transitions_at_12, 
                                              transitions_at_14, 
                                              transitions_at_16, 
                                              transitions_at_18, 
                                              transitions_at_20)
colnames(tRNA_UP_C1_Table_of_new_epimutations) <- c("Up")

#DOWN
tRNA_DOWN_output_C1<- c()
for(i in 1:nrow(tRNA_C1epimutations)){
  tRNA_DOWN_output_C1 <-rbind(tRNA_DOWN_output_C1, DOWN_transition_func(tRNA_C1epimutations[i,], input_name=row.names(tinyRNA_C1epimutations)[i]))}
colnames(tRNA_DOWN_output_C1) <- c("0", "2", "4", "6", "8", "10", "12", "14", "16", "18", "20")
row.names(tRNA_DOWN_output_C1) <- rownames(tRNA_C1epimutations)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tRNA_DOWN_output_C1[,2])
transitions_at_4 <- sum(tRNA_DOWN_output_C1[,3])
transitions_at_6 <- sum(tRNA_DOWN_output_C1[,4])
transitions_at_8 <- sum(tRNA_DOWN_output_C1[,5])
transitions_at_10 <- sum(tRNA_DOWN_output_C1[,6])
transitions_at_12 <- sum(tRNA_DOWN_output_C1[,7])
transitions_at_14 <- sum(tRNA_DOWN_output_C1[,8])
transitions_at_16 <- sum(tRNA_DOWN_output_C1[,9])
transitions_at_18 <- sum(tRNA_DOWN_output_C1[,10])
transitions_at_20 <- sum(tRNA_DOWN_output_C1[,11])

tRNA_DOWN_C1_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                transitions_at_4, 
                                                transitions_at_6,
                                                transitions_at_8, 
                                                transitions_at_10, 
                                                transitions_at_12, 
                                                transitions_at_14, 
                                                transitions_at_16, 
                                                transitions_at_18, 
                                                transitions_at_20)
colnames(tRNA_DOWN_C1_Table_of_new_epimutations) <- c("Down")
#Control2 new epi
#UP

tRNA_UP_output_C2<- c()
row.names(tRNA_C2epimutations)<-tRNA_C2epimutations$tRNAs
tRNA_C2epimutations<-tRNA_C2epimutations[,c(3:12)]
for(i in 1:nrow(tRNA_C2epimutations)){
  tRNA_UP_output_C2 <-rbind(tRNA_UP_output_C2, UP_transition_func(tRNA_C2epimutations[i,], input_name=row.names(tRNA_C2epimutations)[i]))}
colnames(tRNA_UP_output_C2) <- c("0", "2", "4", "8", "10", "12", "14", "16", "18", "20")
row.names(tRNA_UP_output_C2) <- rownames(tRNA_C2epimutations)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tRNA_UP_output_C2[,2])
transitions_at_4 <- sum(tRNA_UP_output_C2[,3])
transitions_at_6 <- NA
transitions_at_8 <- sum(tRNA_UP_output_C2[,4])
transitions_at_10 <- sum(tRNA_UP_output_C2[,5])
transitions_at_12 <- sum(tRNA_UP_output_C2[,6])
transitions_at_14 <- sum(tRNA_UP_output_C2[,7])
transitions_at_16 <- sum(tRNA_UP_output_C2[,8])
transitions_at_18 <- sum(tRNA_UP_output_C2[,9])
transitions_at_20 <- sum(tRNA_UP_output_C2[,10])

tRNA_UP_C2_Table_of_new_epimutations <- rbind(transitions_at_2,
                                              transitions_at_4, 
                                              transitions_at_6,
                                              transitions_at_8, 
                                              transitions_at_10, 
                                              transitions_at_12, 
                                              transitions_at_14, 
                                              transitions_at_16, 
                                              transitions_at_18, 
                                              transitions_at_20)
colnames(tRNA_UP_C2_Table_of_new_epimutations) <- c("Up")

#DOWN
tRNA_DOWN_output_C2<- c()
for(i in 1:nrow(tRNA_C2epimutations)){
  tRNA_DOWN_output_C2 <-rbind(tRNA_DOWN_output_C2, DOWN_transition_func(tRNA_C2epimutations[i,], input_name=row.names(tRNA_C2epimutations)[i]))}
colnames(tRNA_DOWN_output_C2) <- c("0", "2", "4", "8", "10", "12", "14", "16", "18", "20")
row.names(tRNA_DOWN_output_C2) <- rownames(tRNA_C2epimutations)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tRNA_DOWN_output_C2[,2])
transitions_at_4 <- sum(tRNA_DOWN_output_C2[,3])
transitions_at_6 <- NA
transitions_at_8 <- sum(tRNA_DOWN_output_C2[,4])
transitions_at_10 <- sum(tRNA_DOWN_output_C2[,5])
transitions_at_12 <- sum(tRNA_DOWN_output_C2[,6])
transitions_at_14 <- sum(tRNA_DOWN_output_C2[,7])
transitions_at_16 <- sum(tRNA_DOWN_output_C2[,8])
transitions_at_18 <- sum(tRNA_DOWN_output_C2[,9])
transitions_at_20 <- sum(tRNA_DOWN_output_C2[,10])

tRNA_DOWN_C2_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                transitions_at_4, 
                                                transitions_at_6,
                                                transitions_at_8, 
                                                transitions_at_10, 
                                                transitions_at_12, 
                                                transitions_at_14, 
                                                transitions_at_16, 
                                                transitions_at_18, 
                                                transitions_at_20)
colnames(tRNA_DOWN_C2_Table_of_new_epimutations) <- c("Down")
#Low1 new epi
#UP

tRNA_UP_output_L1<- c()
row.names(tRNA_L1epimutations)<-tRNA_L1epimutations$tRNAs
tRNA_L1epimutations<-tRNA_L1epimutations[,c(3:12)]
for(i in 1:nrow(tRNA_L1epimutations)){
  tRNA_UP_output_L1 <-rbind(tRNA_UP_output_L1, UP_transition_func(tRNA_L1epimutations[i,], input_name=row.names(tRNA_L1epimutations)[i]))}
colnames(tRNA_UP_output_L1) <- c("0", "2", "4","6", "8", "12", "14", "16", "18", "20")
row.names(tRNA_UP_output_L1) <- rownames(tRNA_L1epimutations)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tRNA_UP_output_L1[,2])
transitions_at_4 <- sum(tRNA_UP_output_L1[,3])
transitions_at_6 <- sum(tRNA_UP_output_L1[,4])
transitions_at_8 <- sum(tRNA_UP_output_L1[,5])
transitions_at_10 <- NA
transitions_at_12 <- sum(tRNA_UP_output_L1[,6])
transitions_at_14 <- sum(tRNA_UP_output_L1[,7])
transitions_at_16 <- sum(tRNA_UP_output_L1[,8])
transitions_at_18 <- sum(tRNA_UP_output_L1[,9])
transitions_at_20 <- sum(tRNA_UP_output_L1[,10])

tRNA_UP_L1_Table_of_new_epimutations <- rbind(transitions_at_2,
                                              transitions_at_4, 
                                              transitions_at_6,
                                              transitions_at_8, 
                                              transitions_at_10, 
                                              transitions_at_12, 
                                              transitions_at_14, 
                                              transitions_at_16, 
                                              transitions_at_18, 
                                              transitions_at_20)
colnames(tRNA_UP_L1_Table_of_new_epimutations) <- c("Up")

#DOWN
tRNA_DOWN_output_L1<- c()
for(i in 1:nrow(tRNA_L1epimutations)){
  tRNA_DOWN_output_L1 <-rbind(tRNA_DOWN_output_L1, DOWN_transition_func(tRNA_L1epimutations[i,], input_name=row.names(tRNA_L1epimutations)[i]))}
colnames(tRNA_DOWN_output_L1) <- c("0", "2", "4", "6","8", "12", "14", "16", "18", "20")
row.names(tRNA_DOWN_output_L1) <- rownames(tRNA_L1epimutations)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tRNA_DOWN_output_L1[,2])
transitions_at_4 <- sum(tRNA_DOWN_output_L1[,3])
transitions_at_6 <- sum(tRNA_DOWN_output_L1[,4])
transitions_at_8 <- sum(tRNA_DOWN_output_L1[,5])
transitions_at_10 <- NA
transitions_at_12 <- sum(tRNA_DOWN_output_L1[,6])
transitions_at_14 <- sum(tRNA_DOWN_output_L1[,7])
transitions_at_16 <- sum(tRNA_DOWN_output_L1[,8])
transitions_at_18 <- sum(tRNA_DOWN_output_L1[,9])
transitions_at_20 <- sum(tRNA_DOWN_output_L1[,10])

tRNA_DOWN_L1_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                transitions_at_4, 
                                                transitions_at_6,
                                                transitions_at_8, 
                                                transitions_at_10, 
                                                transitions_at_12, 
                                                transitions_at_14, 
                                                transitions_at_16, 
                                                transitions_at_18, 
                                                transitions_at_20)
colnames(tRNA_DOWN_L1_Table_of_new_epimutations) <- c("Down")
#Low2 new epi
#UP

tRNA_UP_output_L2<- c()
row.names(tRNA_L2epimutations)<-tRNA_L2epimutations$tRNAs
tRNA_L2epimutations<-tRNA_L2epimutations[,c(3:12)]
for(i in 1:nrow(tRNA_L2epimutations)){
  tRNA_UP_output_L2 <-rbind(tRNA_UP_output_L2, UP_transition_func(tRNA_L2epimutations[i,], input_name=row.names(tRNA_L2epimutations)[i]))}
colnames(tRNA_UP_output_L2) <- c("0", "2", "4", "8","10", "12", "14", "16", "18", "20")
row.names(tRNA_UP_output_L2) <- rownames(tRNA_L2epimutations)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tRNA_UP_output_L2[,2])
transitions_at_4 <- sum(tRNA_UP_output_L2[,3])
transitions_at_6 <- NA
transitions_at_8 <- sum(tRNA_UP_output_L2[,4])
transitions_at_10 <- sum(tRNA_UP_output_L2[,5])
transitions_at_12 <- sum(tRNA_UP_output_L2[,6])
transitions_at_14 <- sum(tRNA_UP_output_L2[,7])
transitions_at_16 <- sum(tRNA_UP_output_L2[,8])
transitions_at_18 <- sum(tRNA_UP_output_L2[,9])
transitions_at_20 <- sum(tRNA_UP_output_L2[,10])

tRNA_UP_L2_Table_of_new_epimutations <- rbind(transitions_at_2,
                                              transitions_at_4, 
                                              transitions_at_6,
                                              transitions_at_8, 
                                              transitions_at_10, 
                                              transitions_at_12, 
                                              transitions_at_14, 
                                              transitions_at_16, 
                                              transitions_at_18, 
                                              transitions_at_20)
colnames(tRNA_UP_L2_Table_of_new_epimutations) <- c("Up")

#DOWN
tRNA_DOWN_output_L2<- c()
for(i in 1:nrow(tRNA_L2epimutations)){
  tRNA_DOWN_output_L2 <-rbind(tRNA_DOWN_output_L2, DOWN_transition_func(tRNA_L2epimutations[i,], input_name=row.names(tRNA_L2epimutations)[i]))}
colnames(tRNA_DOWN_output_L2) <- c("0", "2", "4", "6","8", "12", "14", "16", "18", "20")
row.names(tRNA_DOWN_output_L2) <- rownames(tRNA_L2epimutations)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tRNA_DOWN_output_L2[,2])
transitions_at_4 <- sum(tRNA_DOWN_output_L2[,3])
transitions_at_6 <- NA
transitions_at_8 <- sum(tRNA_DOWN_output_L2[,4])
transitions_at_10 <- sum(tRNA_DOWN_output_L2[,5])
transitions_at_12 <- sum(tRNA_DOWN_output_L2[,6])
transitions_at_14 <- sum(tRNA_DOWN_output_L2[,7])
transitions_at_16 <- sum(tRNA_DOWN_output_L2[,8])
transitions_at_18 <- sum(tRNA_DOWN_output_L2[,9])
transitions_at_20 <- sum(tRNA_DOWN_output_L2[,10])

tRNA_DOWN_L2_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                transitions_at_4, 
                                                transitions_at_6,
                                                transitions_at_8, 
                                                transitions_at_10, 
                                                transitions_at_12, 
                                                transitions_at_14, 
                                                transitions_at_16, 
                                                transitions_at_18, 
                                                transitions_at_20)
colnames(tRNA_DOWN_L2_Table_of_new_epimutations) <- c("Down")
#High1 new epi
#UP
tRNA_UP_output_H1<- c()
row.names(tRNA_H1epimutations)<-tRNA_H1epimutations$tRNAs
tRNA_H1epimutations<-tRNA_H1epimutations[,c(3:12)]
for(i in 1:nrow(tRNA_H1epimutations)){
  tRNA_UP_output_H1 <-rbind(tRNA_UP_output_H1, UP_transition_func(tRNA_H1epimutations[i,], input_name=row.names(tRNA_H1epimutations)[i]))}
colnames(tRNA_UP_output_H1) <- c("0", "2", "4", "8","10", "12", "14", "16", "18", "20")
row.names(tRNA_UP_output_H1) <- rownames(tRNA_H1epimutations)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tRNA_UP_output_H1[,2])
transitions_at_4 <- sum(tRNA_UP_output_H1[,3])
transitions_at_6 <- NA
transitions_at_8 <- sum(tRNA_UP_output_H1[,4])
transitions_at_10 <- sum(tRNA_UP_output_H1[,5])
transitions_at_12 <- sum(tRNA_UP_output_H1[,6])
transitions_at_14 <- sum(tRNA_UP_output_H1[,7])
transitions_at_16 <- sum(tRNA_UP_output_H1[,8])
transitions_at_18 <- sum(tRNA_UP_output_H1[,9])
transitions_at_20 <- sum(tRNA_UP_output_H1[,10])

tRNA_UP_H1_Table_of_new_epimutations <- rbind(transitions_at_2,
                                              transitions_at_4, 
                                              transitions_at_6,
                                              transitions_at_8, 
                                              transitions_at_10, 
                                              transitions_at_12, 
                                              transitions_at_14, 
                                              transitions_at_16, 
                                              transitions_at_18, 
                                              transitions_at_20)
colnames(tRNA_UP_H1_Table_of_new_epimutations) <- c("Up")

#DOWN
tRNA_DOWN_output_H1<- c()
for(i in 1:nrow(tRNA_H1epimutations)){
  tRNA_DOWN_output_H1 <-rbind(tRNA_DOWN_output_H1, DOWN_transition_func(tRNA_H1epimutations[i,], input_name=row.names(tRNA_H1epimutations)[i]))}
colnames(tRNA_DOWN_output_H1) <- c("0", "2", "4", "8", "10", "12", "14", "16", "18", "20")
row.names(tRNA_DOWN_output_H1) <- rownames(tRNA_H1epimutations)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tRNA_DOWN_output_H1[,2])
transitions_at_4 <- sum(tRNA_DOWN_output_H1[,3])
transitions_at_6 <- NA
transitions_at_8 <- sum(tRNA_DOWN_output_H1[,4])
transitions_at_10 <- sum(tRNA_DOWN_output_H1[,5])
transitions_at_12 <- sum(tRNA_DOWN_output_H1[,6])
transitions_at_14 <- sum(tRNA_DOWN_output_H1[,7])
transitions_at_16 <- sum(tRNA_DOWN_output_H1[,8])
transitions_at_18 <- sum(tRNA_DOWN_output_H1[,9])
transitions_at_20 <- sum(tRNA_DOWN_output_H1[,10])

tRNA_DOWN_H1_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                transitions_at_4, 
                                                transitions_at_6,
                                                transitions_at_8, 
                                                transitions_at_10, 
                                                transitions_at_12, 
                                                transitions_at_14, 
                                                transitions_at_16, 
                                                transitions_at_18, 
                                                transitions_at_20)
colnames(tRNA_DOWN_H1_Table_of_new_epimutations) <- c("Down")

#High2 new epi
#UP

tRNA_UP_output_H2<- c()
row.names(tRNA_H2epimutations)<-tRNA_H2epimutations$tRNAs
tRNA_H2epimutations<-tRNA_H2epimutations[,c(3:13)]
for(i in 1:nrow(tRNA_H2epimutations)){
  tRNA_UP_output_H2 <-rbind(tRNA_UP_output_H2, UP_transition_func(tRNA_H2epimutations[i,], input_name=row.names(tRNA_H2epimutations)[i]))}
colnames(tRNA_UP_output_H2) <- c("0", "2", "4", "6", "8","10", "12", "14", "16", "18", "20")
row.names(tRNA_UP_output_H2) <- rownames(tRNA_H2epimutations)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tRNA_UP_output_H2[,2])
transitions_at_4 <- sum(tRNA_UP_output_H2[,3])
transitions_at_6 <- sum(tRNA_UP_output_H2[,4])
transitions_at_8 <- sum(tRNA_UP_output_H2[,5])
transitions_at_10 <- sum(tRNA_UP_output_H2[,6])
transitions_at_12 <- sum(tRNA_UP_output_H2[,7])
transitions_at_14 <- sum(tRNA_UP_output_H2[,8])
transitions_at_16 <- sum(tRNA_UP_output_H2[,9])
transitions_at_18 <- sum(tRNA_UP_output_H2[,10])
transitions_at_20 <- sum(tRNA_UP_output_H2[,11])

tRNA_UP_H2_Table_of_new_epimutations <- rbind(transitions_at_2,
                                              transitions_at_4, 
                                              transitions_at_6,
                                              transitions_at_8, 
                                              transitions_at_10, 
                                              transitions_at_12, 
                                              transitions_at_14, 
                                              transitions_at_16, 
                                              transitions_at_18, 
                                              transitions_at_20)
colnames(tRNA_UP_H2_Table_of_new_epimutations) <- c("Up")

#DOWN
tRNA_DOWN_output_H2<- c()
for(i in 1:nrow(tRNA_H2epimutations)){
  tRNA_DOWN_output_H2 <-rbind(tRNA_DOWN_output_H2, DOWN_transition_func(tRNA_H2epimutations[i,], input_name=row.names(tRNA_H2epimutations)[i]))}
colnames(tRNA_DOWN_output_H2) <- c("0", "2", "4", "6", "8", "10", "12", "14", "16", "18", "20")
row.names(tRNA_DOWN_output_H2) <- rownames(tRNA_H2epimutations)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(tRNA_DOWN_output_H2[,2])
transitions_at_4 <- sum(tRNA_DOWN_output_H2[,3])
transitions_at_6 <- sum(tRNA_DOWN_output_H2[,4])
transitions_at_8 <- sum(tRNA_DOWN_output_H2[,5])
transitions_at_10 <- sum(tRNA_DOWN_output_H2[,6])
transitions_at_12 <- sum(tRNA_DOWN_output_H2[,7])
transitions_at_14 <- sum(tRNA_DOWN_output_H2[,8])
transitions_at_16 <- sum(tRNA_DOWN_output_H2[,9])
transitions_at_18 <- sum(tRNA_DOWN_output_H2[,10])
transitions_at_20 <- sum(tRNA_DOWN_output_H2[,11])

tRNA_DOWN_H2_Table_of_new_epimutations <- rbind(transitions_at_2,
                                                transitions_at_4, 
                                                transitions_at_6,
                                                transitions_at_8, 
                                                transitions_at_10, 
                                                transitions_at_12, 
                                                transitions_at_14, 
                                                transitions_at_16, 
                                                transitions_at_18, 
                                                transitions_at_20)
colnames(tRNA_DOWN_H2_Table_of_new_epimutations) <- c("Down")
#Plot creation
Transitions<-row.names(tRNA_UP_C1_Table_of_new_epimutations)
colnames(tRNA_UP_C1_Table_of_new_epimutations)<-c("UP")
colnames(tRNA_DOWN_C1_Table_of_new_epimutations)<-c("DOWN")
tRNA_allepiC1 <- cbind(Transitions,tRNA_UP_C1_Table_of_new_epimutations,tRNA_DOWN_C1_Table_of_new_epimutations)
tRNA_allepiC1<-as.data.frame(tRNA_allepiC1)
tRNA_allepiC1 <- tRNA_allepiC1 %>%
  # Creating an empty column:
  add_column(Lineage = "C1", .before="UP")
colnames(tRNA_UP_C2_Table_of_new_epimutations)<-c("UP")
colnames(tRNA_DOWN_C2_Table_of_new_epimutations)<-c("DOWN")
tRNA_allepiC2 <- cbind(Transitions,tRNA_UP_C2_Table_of_new_epimutations,tRNA_DOWN_C2_Table_of_new_epimutations)
tRNA_allepiC2<-as.data.frame(tRNA_allepiC2)
tRNA_allepiC2 <- tRNA_allepiC2 %>%
  # Creating an empty column:
  add_column(Lineage = "C2", .before="UP")
tRNA_allC <- rbind(tRNA_allepiC1, tRNA_allepiC2)
tRNA_allC <- tRNA_allC %>%
  # Creating an empty column:
  add_column(Condition = "Control", .before="Lineage")

colnames(tRNA_UP_L1_Table_of_new_epimutations)<-c("UP")
colnames(tRNA_DOWN_L1_Table_of_new_epimutations)<-c("DOWN")
tRNA_allepiL1 <- cbind(Transitions,tRNA_UP_L1_Table_of_new_epimutations,tRNA_DOWN_L1_Table_of_new_epimutations)
tRNA_allepiL1<-as.data.frame(tRNA_allepiL1)
tRNA_allepiL1 <- tRNA_allepiL1 %>%
  # Creating an empty column:
  add_column(Lineage = "L1", .before="UP")
colnames(tRNA_UP_L2_Table_of_new_epimutations)<-c("UP")
colnames(tRNA_DOWN_L2_Table_of_new_epimutations)<-c("DOWN")
tRNA_allepiL2 <- cbind(Transitions,tRNA_UP_L2_Table_of_new_epimutations,tRNA_DOWN_L2_Table_of_new_epimutations)
tRNA_allepiL2<-as.data.frame(tRNA_allepiL2)
tRNA_allepiL2 <- tRNA_allepiL2 %>%
  # Creating an empty column:
  add_column(Lineage = "L2", .before="UP")
tRNA_allL <- rbind(tRNA_allepiL1, tRNA_allepiL2)
tRNA_allL <- tRNA_allL %>%
  # Creating an empty column:
  add_column(Condition = "Low dose", .before="Lineage")

colnames(tRNA_UP_H1_Table_of_new_epimutations)<-c("UP")
colnames(tRNA_DOWN_H1_Table_of_new_epimutations)<-c("DOWN")
tRNA_allepiH1 <- cbind(Transitions,tRNA_UP_H1_Table_of_new_epimutations,tRNA_DOWN_H1_Table_of_new_epimutations)
tRNA_allepiH1<-as.data.frame(tRNA_allepiH1)
tRNA_allepiH1 <- tRNA_allepiH1 %>%
  # Creating an empty column:
  add_column(Lineage = "H1", .before="UP")
colnames(tRNA_UP_H2_Table_of_new_epimutations)<-c("UP")
colnames(tRNA_DOWN_H2_Table_of_new_epimutations)<-c("DOWN")
tRNA_allepiH2 <- cbind(Transitions,tRNA_UP_H2_Table_of_new_epimutations,tRNA_DOWN_H2_Table_of_new_epimutations)
tRNA_allepiH2<-as.data.frame(tRNA_allepiH2)
tRNA_allepiH2 <- tRNA_allepiH2 %>%
  # Creating an empty column:
  add_column(Lineage = "H2", .before="UP")
tRNA_allH <- rbind(tRNA_allepiH1, tRNA_allepiH2)
tRNA_allH <- tRNA_allH %>%
  # Creating an empty column:
  add_column(Condition = "High dose", .before="UP")
tRNA_allnewepi<-rbind(tRNA_allC,tRNA_allL,tRNA_allH)
rownames(tRNA_allnewepi)<-NULL

tRNA_allnewepi$UP <- as.numeric(tRNA_allnewepi$UP)
tRNA_allnewepi$DOWN <- as.numeric(tRNA_allnewepi$DOWN)
tRNA_allnewepi$Total=rowSums(cbind(tRNA_allnewepi$UP,tRNA_allnewepi$DOWN),na.rm=FALSE)

#Fig.5.G
tRNA_allnewepi$Condition <- fct_relevel(tRNA_allnewepi$Condition, c("Control", "Low dose", "High dose"))
p1 <- ggplot(tRNA_allnewepi, aes(x=Condition, y=Total, color = Condition,group=Condition)) + 
  geom_boxplot(fatten = 1, lwd = 1, width=0.5)+
  scale_color_manual(values=c("cornflowerblue", "darkgreen", "red"))+
  labs(y = "Number of new epimutations", x = "\nCondition")+
  geom_dotplot(binaxis='y', binwidth = 0.25, stackdir='center',  dotsize=0.9, binpositions = "all", stackgroups = TRUE)+
  theme_bw()+
  theme_linedraw()+
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(color = "black", size = 20, angle = 0, face = "plain"),
        axis.text.y = element_text(color = "black", size = 20, angle = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 16, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 16, angle = 90, face = "bold"))+
  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))+  
  scale_fill_brewer(palette="Pastel2")+
  ggtitle(paste(""))

ggdensity(tRNA_allnewepi$Total, 
          main = "Density plot of epimutations_number",
          xlab = "Total number of epimutations_number")
shapiro.test(tRNA_allnewepi$Total)
kruskal.test(Total ~ Condition, data = tRNA_allnewepi)
dunnTest(Total ~ Condition, data = tRNA_allnewepi)
pairwise.wilcox.test(tRNA_allnewepi$Total, tRNA_allnewepi$Condition, p.adj='bonferroni', exact=F)

#Sup.Fig.6.G
Generation<-c("2","4","6","8","10","12","14","16","18","20",
              "2","4","6","8","10","12","14","16","18","20",
              "2","4","6","8","10","12","14","16","18","20",
              "2","4","6","8","10","12","14","16","18","20",
              "2","4","6","8","10","12","14","16","18","20",
              "2","4","6","8","10","12","14","16","18","20")
tRNA_allnewepi_detailed<-cbind(tRNA_allnewepi,Generation)
tRNA_allnewepi_detailed$Generation <- fct_relevel(tRNA_allnewepi_detailed$Generation, c("2", "4", "6","8","10","12","14","16","18","20"))
tRNA_allnewepi_detailed$Condition <- fct_relevel(tRNA_allnewepi_detailed$Condition, c("Control", "Low dose", "High dose"))
tRNA_allnewepi_detailed$Lineage  <- fct_relevel(tRNA_allnewepi_detailed$Lineage , c("C1", "C2", "L1","L2","H1","H2"))
tRNA_allnewepi_detailed_plot <- ggplot(tRNA_allnewepi_detailed, aes(x=Generation , y=Total, color = Condition)) + 
  geom_boxplot(fatten = 1, lwd = 1, width=0.5)+
  scale_color_manual(values=c("cornflowerblue", "darkgreen", "red"))+
  labs(y = "Number of new epimutations", x = "\nGeneration")+
  geom_jitter(aes(shape=Lineage, size=8),
              position=position_jitter(width = 0.2,
                                       height = 0.2))+
  theme_bw()+
  theme(axis.text.x = element_text(color = "black", size = 20, angle = 0, face = "plain"),
        axis.text.y = element_text(color = "black", size = 20, angle = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 16, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 16, angle = 90, face = "bold"))+
  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))+  
  scale_fill_brewer(palette="Pastel2")+
  ggtitle(paste(""))
tRNA_allnewepi_detailed_plot

#Get raw data
write.xlsx(tRNA_allnewepi_detailed,"Data_Fig_5_G_&_Sup_Fig_6_G.xlsx")
#--------------------------
#####Sup.Fig.10.C####
#Calc duration epimut
#Data preparation
#Remplacement of missing columns with data following generation
tRNA_C2epimutationsbis<-cbind(tRNA_C2epimutations, rep(tRNA_C2epimutations[4],1))
tRNA_C2epimutationsbis<-tRNA_C2epimutationsbis[,c(1,2,3,11,4,5,6,7,8,9,10)]
colnames(tRNA_C2epimutationsbis) <- c('0','2','4','6','8','10','12','14','16','18','20')

tRNA_L1epimutationsbis<-cbind(tRNA_L1epimutations, rep(tRNA_L1epimutations[6],1))
tRNA_L1epimutationsbis<-tRNA_L1epimutationsbis[,c(1,2,3,4,5,11,6,7,8,9,10)]
colnames(tRNA_L1epimutationsbis) <- c('0','2','4','6','8','10','12','14','16','18','20')
write.xlsx(tRNA_L1epimutationsbis,"tRNA_L1epimutationsbis.xlsx")

tRNA_L2epimutationsbis<-cbind(tRNA_L2epimutations, rep(tRNA_L2epimutations[4],1))
tRNA_L2epimutationsbis<-tRNA_L2epimutationsbis[,c(1,2,3,11,4,5,6,7,8,9,10)]
colnames(tRNA_L2epimutationsbis) <- c('0','2','4','6','8','10','12','14','16','18','20')
write.xlsx(tRNA_L2epimutationsbis,"tRNA_L2epimutationsbis.xlsx")

tRNA_H1epimutationsbis<-cbind(tRNA_H1epimutations, rep(tRNA_H1epimutations[4],1))
tRNA_H1epimutationsbis<-tRNA_H1epimutationsbis[,c(1,2,3,11,4,5,6,7,8,9,10)]
colnames(tRNA_H1epimutationsbis) <- c('0','2','4','6','8','10','12','14','16','18','20')
write.xlsx(tRNA_H1epimutationsbis,"tRNA_H1epimutationsbis.xlsx")

write.xlsx(tRNA_H2epimutations,"tRNA_H2epimutations.xlsx")

tRNAs_all_epi<-list(tRNA_C1epimutations,tRNA_C2epimutationsbis,tRNA_L1epimutationsbis,tRNA_L2epimutationsbis,tRNA_H1epimutationsbis,tRNA_H2epimutations)
names(tRNAs_all_epi)<-c("C1","C2","L1","L2","H1","H2")

#Calculation
tRNA_sorted_output_list_with_missing_as_muts <- list()

for(e in 1:length(tRNAs_all_epi)){
  
  main_frame <- as.data.frame(tRNAs_all_epi[[e]])
  
  Overall_output <- c()
  
  for(t in 1:nrow(main_frame)){
    
    vector_in <- main_frame[t, ]
    gen_names <- as.numeric(colnames(vector_in))
    
    number_transitions <- 0
    onset_gen <- 0
    tempL <- 0
    output_frame <- c()
    
    is_up <- 0
    is_down <- 0
    complete <- 0
    
    # first deal with no epimutations in the lineage for a locus/gene
    
    if(sum(abs(vector_in)) == 0){
      
      length = 0
      is_up <- 0
      is_down <- 0
      gene_name <- rownames(vector_in)
      gene <- strsplit(gene_name, ":")[[1]][4]
      
      number_transitions <- 0 
      onset_gen <- 0
      
      Lineage <- names(tRNAs_all_epi[e])
      complete <- 0
      
      save_mut <- data.frame(gene_name, gene, number_transitions, onset_gen, length,  complete,  Lineage, is_up, is_down)   
      
      output_frame <- rbind(output_frame, save_mut)
      
    }
    
    else
      
      if(sum(abs(vector_in)) > 0){
        
        for(i in 2:length(vector_in)){
          
          
          # if it is an UP transition, turn off any down transitions and start new epimutation
          
          if(vector_in[i]==1&vector_in[i-1]== 0|vector_in[i]== 1&vector_in[i-1]== -1) {
            
            if(is_down == 1){ # a down epimutation that is now turning off
              
              gene_name <- rownames(vector_in)
              gene <- strsplit(gene_name, ":")[[1]][4]
              
              number_transitions <- number_transitions + 1 
              
              length <-tempL
              
              Lineage <- names(tRNAs_all_epi[e])
              
              complete <- 1
              
              save_mut <- data.frame(gene_name, gene, number_transitions, onset_gen, length,  complete, Lineage, is_up, is_down)   
              
              output_frame <- rbind(output_frame, save_mut)
              
            }
            
            # start the UP transition
            
            is_up <- 1
            is_down <- 0
            tempL <- 1
            
            number_transitions <- 1 # there is a transition to an UP epimutation at generation[i]
            onset_gen <- gen_names[i]
            
            onset_check <- i -1
          }
          
          # if it is a new DOWN transition, turn off any up transitions and start new epimutation
          
          if(vector_in[i]==-1&vector_in[i-1]== 0|vector_in[i]== -1&vector_in[i-1]== 1) {
            
            if(is_up == 1){ # an up epimutation that is now turning off
              
              gene_name <- rownames(vector_in)
              gene <- strsplit(gene_name, ":")[[1]][4]
              
              number_transitions <- number_transitions + 1
              
              length <- tempL
              
              Lineage <- names(tRNAs_all_epi[e])
              
              complete <- 1
              
              save_mut <- data.frame(gene_name, gene, number_transitions, onset_gen, length,  complete, Lineage, is_up, is_down)   
              
              output_frame <- rbind(output_frame, save_mut)
              
            }
            
            # start the down transition
            
            is_down <- 1
            is_up <- 0
            tempL <- 1
            number_transitions <- 1 # there is a transition to a DOWN epimutation at generation[i]
            onset_gen <- gen_names[i]
            onset_check <- i-1
          }
          
          # if it is an UP epimutation, does it continue, does it turn off, does it extend to last gen?
          
          if(is_up==1){
            
            
            if(vector_in[i]==1&vector_in[i-1]==1)
              
            {tempL<-tempL+(gen_names[i]-gen_names[i-1])}  # the epimutation continues
            
            if(vector_in[i]==0&vector_in[i-1]==1){
              
              gene_name <- rownames(vector_in)
              gene <- strsplit(gene_name, ":")[[1]][4]
              
              number_transitions <- number_transitions +1  # there is a transition to OFF from an UP epimutation at generation[i]
              
              length <- tempL
              
              Lineage <- names(tRNAs_all_epi[e])
              
              complete <- 1
              
              save_mut <- data.frame(gene_name, gene, number_transitions, onset_gen, length,  complete, Lineage, is_up, is_down)   
              
              output_frame <- rbind(output_frame, save_mut)
              
              is_up <- 0
              is_down <- 0
              tempL <- 0
              length <- 0
              complete <- 0
              
            }
            
            # deal with an up epimutation extending to end of lineage
            
            if(vector_in[i]==1&i==length(vector_in)){
              
              gene_name <- rownames(vector_in)
              
              gene <- strsplit(gene_name, ":")[[1]][4]
              
              number_transitions <- number_transitions + 0
              
              length <- tempL
              
              Lineage <- names(tRNAs_all_epi[e])
              
              complete <- 0
              
              save_mut <- data.frame(gene_name, gene, number_transitions, onset_gen, length,  complete, Lineage, is_up, is_down)   
              
              output_frame <- rbind(output_frame, save_mut)
              
              is_up <- 0
              is_down <- 0
              length <- 0
              tempL <- 0
              complete <- 0
              
            }
            
          }
          
          # if it is a DOWN epimutation, does it continue, does it turn off, does it extend to last gen?
          
          if(is_down==1){
            
            if(vector_in[i]==-1&vector_in[i-1]==-1)
              
            {tempL<-tempL+(gen_names[i]-gen_names[i-1])}
            
            
            if(vector_in[i]==0&vector_in[i-1]==-1){
              
              gene_name <- rownames(vector_in)
              
              gene <- strsplit(gene_name, ":")[[1]][4]
              
              number_transitions <- number_transitions +1  # there is a transition to OFF from a DOWN epimutation at generation[i]
              
              length <- tempL
              
              Lineage <- names(tRNAs_all_epi[e])
              
              complete <- 1
              
              save_mut <- data.frame(gene_name, gene, number_transitions, onset_gen, length,  complete, Lineage, is_up, is_down)   
              
              output_frame <- rbind(output_frame, save_mut)
              
              is_up <- 0
              is_down <- 0
              tempL <- 0
              complete <- 0
              
            }
            
            # deal with a down epimutation extending to end of lineage
            
            if(vector_in[i]==-1&i==length(vector_in)){
              
              gene_name <- rownames(vector_in)
              
              gene <- strsplit(gene_name, ":")[[1]][4]
              
              number_transitions <- number_transitions + 0
              
              length <- tempL
              
              Lineage <- names(tRNAs_all_epi[e])
              
              complete <- 0
              
              save_mut <- data.frame(gene_name, gene, number_transitions, onset_gen, length,  complete, Lineage, is_up, is_down)   
              output_frame <- rbind(output_frame, save_mut)
              
              is_up <- 0
              is_down <- 0
              tempL <- 0
              complete <- 0
              
            }
            
          }    
          
        }} 
    
    Overall_output <- rbind(Overall_output, output_frame)  
    
  }
  
  tRNA_sorted_output_list_with_missing_as_muts[[e]] <- Overall_output
}

tRNA_sorted_epiduration<-c()
names(tRNA_sorted_output_list_with_missing_as_muts) <- names(tRNAs_all_epi)
tRNA_sorted_output_list_with_missing_as_muts$C1 <- tRNA_sorted_output_list_with_missing_as_muts$C1 %>%
  # Creating an empty column:
  add_column(Condition = "Control")
tRNA_sorted_epiduration<-tRNA_sorted_output_list_with_missing_as_muts$C1
tRNA_sorted_output_list_with_missing_as_muts$C2 <- tRNA_sorted_output_list_with_missing_as_muts$C2 %>%
  # Creating an empty column:
  add_column(Condition = "Control")
tRNA_sorted_epiduration<-rbind(tRNA_sorted_epiduration, tRNA_sorted_output_list_with_missing_as_muts$C2)
tRNA_sorted_output_list_with_missing_as_muts$L1 <- tRNA_sorted_output_list_with_missing_as_muts$L1 %>%
  # Creating an empty column:
  add_column(Condition = "Low dose")
tRNA_sorted_epiduration<-rbind(tRNA_sorted_epiduration, tRNA_sorted_output_list_with_missing_as_muts$L1)
tRNA_sorted_output_list_with_missing_as_muts$L2 <- tRNA_sorted_output_list_with_missing_as_muts$L2 %>%
  # Creating an empty column:
  add_column(Condition = "Low dose")
tRNA_sorted_epiduration<-rbind(tRNA_sorted_epiduration, tRNA_sorted_output_list_with_missing_as_muts$L2)
tRNA_sorted_output_list_with_missing_as_muts$H1 <- tRNA_sorted_output_list_with_missing_as_muts$H1 %>%
  # Creating an empty column:
  add_column(Condition = "High dose")
tRNA_sorted_epiduration<-rbind(tRNA_sorted_epiduration, tRNA_sorted_output_list_with_missing_as_muts$H1)
tRNA_sorted_output_list_with_missing_as_muts$H2 <- tRNA_sorted_output_list_with_missing_as_muts$H2 %>%
  # Creating an empty column:
  add_column(Condition = "High dose")
tRNA_sorted_epiduration<-rbind(tRNA_sorted_epiduration, tRNA_sorted_output_list_with_missing_as_muts$H2)
head(tRNA_sorted_epiduration)
tRNA_sorted_epimut_wt_0<-subset(tRNA_sorted_epiduration,tRNA_sorted_epiduration$length !=0)
save(tRNA_sorted_epimut_wt_0,file="tRNA_sorted_epimut_wt_0.Rdata")
#Survival analysis
tRNA_sorted_epimut_wt_0$Condition <- fct_relevel(tRNA_sorted_epimut_wt_0$Condition, c("Control", "Low dose","High dose"))
CoxMod<-coxph(Surv(length,complete)~Condition,data=tRNA_sorted_epimut_wt_0)
ggforest(model=CoxMod,data=tRNA_sorted_epimut_wt_0,fontsize = 0.8, noDigits = 2)

#Get raw data
write.xlsx(tRNA_sorted_epimut_wt_0,"Data_Sup_Fig_9_C.xlsx")
#--------------------------
#####Fig.5.H - Comparion 22G-RNAs vs. tRNAs fragments epimutations duration####
#Combined 22G and tRNAs
tiny_22G_epimut<-tiny_all_epimut %>% filter(grepl('22G', gene_name))
tRNA_sorted_epimut_wt_0

#Plot construction
colnames(tiny_22G_epimut)[1] <- "names"
tiny_22G_epimut <- cbind(rep("22G", length = nrow(tiny_22G_epimut)), tiny_22G_epimut)
colnames(tiny_22G_epimut)[1] <- "Data_type"

colnames(tRNA_sorted_epimut_wt_0)[1] <- "names"
tRNA_sorted_epimut_wt_0 <- cbind(rep("tRNAs", length = nrow(tRNA_sorted_epimut_wt_0)), tRNA_sorted_epimut_wt_0)
colnames(tRNA_sorted_epimut_wt_0)[1] <- "Data_type"
tRNA_sorted_epimut_wt_0 <- subset(tRNA_sorted_epimut_wt_0,  select = -12)

ALL_SETS_epimutations <- rbind(tiny_22G_epimut,tRNA_sorted_epimut_wt_0)
All_SETS_surv_object <- Surv(time = ALL_SETS_epimutations$length, event = ALL_SETS_epimutations$complete)
ALL_SETS_fit <- survfit(All_SETS_surv_object ~ 1, data = ALL_SETS_epimutations)
ALL_SETS_plot <- ggsurvplot(ALL_SETS_fit, ALL_SETS_epimutations, censor = T)+
  ggtitle("All data compare") +
  xlab("Time (generations)")

ALL_SETS_plot$plot <- ALL_SETS_plot$plot + theme(legend.text = element_text(size = 14))
ALL_SETS_fit_data_type <-  survfit(All_SETS_surv_object ~ Data_type, data = ALL_SETS_epimutations)
ALL_SETS_plot_data_type <- ggsurvplot(ALL_SETS_fit_data_type,  ALL_SETS_epimutations, censor = T, break.time.by= 1, pval = T,  pval.coord = c(0, 0.05), 
                                      font.main = c(16, "bold"),
                                      font.x = 14,
                                      font.y = 14,
                                      font.legend = 14,
                                      font.tickslab = 12)+
  # surv.median.line = "hv")+
  ggtitle("Changes survival all conditions")+
  xlab("Time (generations)")

ALL_SETS_plot_data_type$plot +theme(plot.title = element_text(hjust = 0.5))
# n.b. p value derived in ggsurv plot through log rank test
CoxMod<-coxph(Surv(length,complete)~Data_type+Condition,data=ALL_SETS_epimutations)
ggforest(CoxMod)

#Get raw data
write.xlsx(ALL_SETS_epimutations,"Data_Fig_5_H.xlsx")
#--------------------------
#####Sup.Fig.9.A & B####
GetCountsAll<-function(file_in){temp<-read.table(file_in);temp<-read.table(text=temp[,4], sep="-");out<-as.numeric(temp[,2]); return(sum(out))}

GetCount<-function(fileIn){temp<-read.table(fileIn, sep="\t", stringsAsFactors=F);sel<-which(temp[,6]==temp[,12]&temp[,ncol(temp)]>25); M<-read.table(text=temp[sel,10], sep="-", stringsAsFactors=F);tRNACount<-data.frame(temp[sel,4],as.numeric(M[,2]));return(tRNACount)}



tRNAMuts<-list.files(pattern="tRNA.bed")


tRNAMutsums<-c()
tRNAMutDist<-list()
for(i in 1:length(tRNAMuts)){
  CountDist<-GetCount(tRNAMuts[i])
  
  tRNAMutsums<-c(tRNAMutsums,sum(CountDist[,2]))
  tRNAMutDist[[i]]<-CountDist
  
  
}

ContFiles<-list.files(pattern="fasta.bed$")
allreads<-c()

for(i in 1:length(ContFiles)){allreads<-c(allreads,GetCountsAll(ContFiles[i]))}

NormtRNAAll<-data.frame(c(rep("ergo1",3),rep("N2",3),rep("rde1",3),rep("wago10",3)),tRNAMutsums/allreads)
NormtRNAAll<-NormtRNAAll[c(4,5,6,1,2,3,7:12),]
pointCol<-c("black","black","black","blue","blue","blue","orange","orange","orange","brown","brown","brown")
plot(c(1,1,1,2,2,2,3,3,3,4,4,4),NormtRNAAll[,2]*1e6,xaxt="null",xlab="mutant", ylab="total tRNAcounts rpm",pch=18,cex=2,col=pointCol,ylim=c(2000,7000))
axis(side=1, at=c(1,2,3,4),labels=unique(NormtRNAAll[,1]))
dev.copy(pdf, "tRNAoverallcounts.pdf")
dev.off()

#wago10 shows a small decrease, which is interesting. 
#to test whether this is real, we need to look systematically at all tRNAs


CollectType<-function(tRNADist){
  
  #sortbytRNA
  Label<-read.table(text=tRNADist[,1],sep="_")
  Label<-Label[,3]
  Label2<-read.table(text=Label,sep="-")
  Label2<-paste(Label2[,2],Label2[,3],sep=":")
  Types<-table(Label2)
  Out<-rep(0, length=length(Types))
  for(i in 1:length(Types)){
    sel<-which(Label2==names(Types[i]))
    Out[i]<-sum(tRNADist[sel,2])
  }
  
  return(cbind(Types,Out))
  
  
  
}

tRNAReadsAll<-list()
AllTypes<-c()
for(i in 1:length(tRNAMutDist)){
  
  tRNAReadsAll[[i]]<-CollectType(tRNAMutDist[[i]])
  AllTypes<-c(AllTypes, row.names(tRNAReadsAll[[i]]))
}

AllTypes<-unique(AllTypes)

selFiles<-c(4,5,6,10,11,12)
Wago10vN2<-c()
for(i in 1:length(selFiles)){
  A<-match(AllTypes,row.names(tRNAReadsAll[[selFiles[i]]]))
  Wago10vN2<-cbind(Wago10vN2,tRNAReadsAll[[selFiles[i]]][A,2])
  
}
colnames(Wago10vN2)<-c("N2_r1","N2_r2","N2_r3","Wago10_r1","Wago10_r2","Wago10_r3")
for(i in 1:6){Wago10vN2[,i]<-Wago10vN2[,i]/allreads[selFiles[i]]}
for(i in 1:6){
  L<-which(Wago10vN2[,i]>0)
  Wago10vN2[-L,i]<-0
  
}
Wago10vN2<-Wago10vN2[which(apply(Wago10vN2,1,max)>0),]

VolData<-c()
for(i in 1:nrow(Wago10vN2)){
  WT<-t.test(Wago10vN2[i,1:3],Wago10vN2[i,4:6])
  VolData<-c(VolData,WT$p.value)
}
VolData<-cbind(VolData,apply(Wago10vN2[,1:3],1,mean),apply(Wago10vN2[,4:6],1,mean))

plot(log2(VolData[,2]),log2(VolData[,3]),pch=18,col="red",xlab="log2mean_N2",ylab="log2mean_Wago10",xlim=c(-25,-9),ylim=c(-25,-9))
LM<-lm(log2(VolData[,3])~log2(VolData[,2]))
abline(LM)
text(-9,-11.2,"GluTTC",cex=0.5)
text(-9.5,-10.5,"GluCTC",cex=0.5)
dev.copy(pdf, "tRNA_types_wagovWT.pdf")
dev.off()

#are any tRNA types differentially expressed?
#group by tRNA
tRNAtypes<-read.table(text=row.names(Wago10vN2),sep=":")
tRNAtypeTab<-table(tRNAtypes[,1])
tRNAtypeTab<-cbind(tRNAtypeTab,matrix(0, ncol=3,nrow=length(tRNAtypeTab)))
for(i in 1:nrow(tRNAtypeTab)){
  L<-which(tRNAtypes[,1]==row.names(tRNAtypeTab)[i])
  if(length(L==1)){N2_read<-Wago10vN2[L,1:3]
  Wago_read<-Wago10vN2[L,4:6]
  }else{
    N2_read<-apply(Wago10vN2[L,1:3],2,sum)
    Wago_read<-apply(Wago10vN2[L,4:6],2,sum)}
  WT<-wilcox.test(N2_read,Wago_read)
  tRNAtypeTab[i,2:4]<-c(median(N2_read),median(Wago_read),WT$p.value)
}
row.names(tRNAtypeTab)[(which(p.adjust(tRNAtypeTab[,4])<0.05))]

#[1] "Glu" 

plot(c(rep(1,length=6),rep(2,length=6)),c(Wago10vN2[c(13:14),1:3],Wago10vN2[c(13:14),4:6])*1e6,col=c(rep("black",length=6),rep("orange",length=6)),pch=c(18,18,18,16,16,16,18,18,18,16,16,16),xaxt="null",xlim=c(0,3),ylab="reads/million",xlab="",cex=1.5)
axis(side=1, at=c(1,2),labels=c("N2","wago10"))
legend("topright", c("GluCTC","GluTCC"),pch=c(18,16))
text(2.5,500,paste("Adjusted p=",round(p.adjust(tRNAtypeTab[,4])["Glu"],digits=3)))
dev.copy(pdf, "GlutRNA_Wago10vN2.pdf")
dev.off()



save.image("tRNA_analysis_inMutantData.Rdata")


#--------------------------
####Sup.Fig.10.A & B#####
#Table of number of epimut per loci
table_loci<-read.xlsx("table_epimut_loci_conditions.xlsx")
colnames(table_loci)<-c("loci","nb_epimut_control","max_length_control","nb_epimut_low","max_length_low","nb_epimut_high","max_length_high")
table_nb_control<-table_loci[,1:2]
colnames(table_nb_control)<-c("loci","Total")
table_nb_control <- table_nb_control %>%
  # Creating an empty column:
  add_column(Condition = "Control", .before="Total")

table_nb_low<-cbind(table_loci[,1], table_loci[,4])
colnames(table_nb_low)<-c("loci","Total")
table_nb_low<-as.data.frame(table_nb_low)
table_nb_low <- table_nb_low %>%
  # Creating an empty column:
  add_column(Condition = "Low dose", .before="Total")

table_nb_high<-cbind(table_loci[,1], table_loci[,6])
colnames(table_nb_high)<-c("loci","Total")
table_nb_high<-as.data.frame(table_nb_high)
table_nb_high <- table_nb_high %>%
  # Creating an empty column:
  add_column(Condition = "High dose", .before="Total")

table_nb_epimut<-rbind(table_nb_control,table_nb_low,table_nb_high)
table_nb_epimut_wt_0<-subset(table_nb_epimut, table_nb_epimut$Total !=0)

table_nb_epimut_wt_0$Condition <- fct_relevel(table_nb_epimut_wt_0$Condition, c("Control", "Low dose", "High dose"))
p<-ggplot(table_nb_epimut_wt_0, aes(x=loci, y=Total, fill=Condition)) +
  geom_bar(stat="identity", position=position_dodge())
p + scale_y_discrete(limits=c("1", "2","3","4","5")) + coord_flip() + labs(y= "Total number of epimutations per condition", x = "loci")+scale_color_manual(values=c("cornflowerblue", "coral1", "red"))

#Get raw data
write.xlsx(table_nb_epimut_wt_0,"Data_Sup_Fig_9_A.xlsx")

#Table of max length of epimut per loci
table_loci<-read.xlsx("table_epimut_loci_conditions.xlsx")
colnames(table_loci)<-c("loci","nb_epimut_control","max_length_control","nb_epimut_low","max_length_low","nb_epimut_high","max_length_high")

table_lg_control<-cbind(table_loci[,1], table_loci[,3])
colnames(table_lg_control)<-c("loci","Max_Length")
table_lg_control<-as.data.frame(table_lg_control)
table_lg_control <- table_lg_control %>%
  # Creating an empty column:
  add_column(Condition = "Control", .before="Max_Length")

table_lg_low<-cbind(table_loci[,1], table_loci[,5])
colnames(table_lg_low)<-c("loci","Max_Length")
table_lg_low<-as.data.frame(table_lg_low)
table_lg_low <- table_lg_low %>%
  # Creating an empty column:
  add_column(Condition = "Low dose", .before="Max_Length")

table_lg_high<-cbind(table_loci[,1], table_loci[,7])
colnames(table_lg_high)<-c("loci","Max_Length")
table_lg_high<-as.data.frame(table_lg_high)
table_lg_high <- table_lg_high %>%
  # Creating an empty column:
  add_column(Condition = "High dose", .before="Max_Length")

table_lg_epimut<-rbind(table_lg_control,table_lg_low,table_lg_high)
table_lg_epimut_wt_0<-subset(table_lg_epimut, table_lg_epimut$Max_Length !=0)

table_lg_epimut_wt_0$Condition <- fct_relevel(table_lg_epimut_wt_0$Condition, c("Control", "Low dose", "High dose"))
p<-ggplot(table_lg_epimut_wt_0, aes(x=loci, y=Max_Length, fill=Condition)) +
  geom_bar(stat="identity", position=position_dodge())
p + scale_y_discrete(limits=c("1", "2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19")) + coord_flip() + labs(y= "Max length of epimutations per condition", x = "loci")

#Get raw data
write.xlsx(table_lg_epimut_wt_0,"Data_Sup_Fig_9_B.xlsx")
#--------------------------
#####Fig.6.A - GO term analysis for tRNAs associated genes####
#Data prep
setwd("~/Documents/Cisplatin project/Data_analysis/Analysis")
tRNA_target<-read.xlsx("tRNA_targetedGenes2.xlsx")
tRNA_target_names<-tRNA_target$WB

load("RNA_C1epimutations.Rdata")
RNA_C1epimutations$genes <-sapply(strsplit(RNA_C1epimutations$genes,":"), `[`, 5)
RNA_C1epimutations_target <- RNA_C1epimutations[RNA_C1epimutations$genes %in% tRNA_target_names, ]

load("RNA_C2epimutations.Rdata")
RNA_C2epimutations$genes <-sapply(strsplit(RNA_C2epimutations$genes,":"), `[`, 5)
RNA_C2epimutations_target <- RNA_C2epimutations[RNA_C2epimutations$genes %in% tRNA_target_names, ]

load("RNA_L1epimutations.Rdata")
RNA_L1epimutations$genes <-sapply(strsplit(RNA_L1epimutations$genes,":"), `[`, 5)
RNA_L1epimutations_target <- RNA_L1epimutations[RNA_L1epimutations$genes %in% tRNA_target_names, ] 

load("RNA_L2epimutations.Rdata")
RNA_L2epimutations$genes <-sapply(strsplit(RNA_L2epimutations$genes,":"), `[`, 5)
RNA_L2epimutations_target <- RNA_L2epimutations[RNA_L2epimutations$genes %in% tRNA_target_names, ] 

load("RNA_H1epimutations.Rdata")
RNA_H1epimutations$genes <-sapply(strsplit(RNA_H1epimutations$genes,":"), `[`, 5)
RNA_H1epimutations_target <- RNA_H1epimutations[RNA_H1epimutations$genes %in% tRNA_target_names, ] 

load("RNA_H2epimutations.Rdata")
RNA_H2epimutations$genes <-sapply(strsplit(RNA_H2epimutations$genes,":"), `[`, 5)
RNA_H2epimutations_target <- RNA_H2epimutations[RNA_H2epimutations$genes %in% tRNA_target_names, ] 

colnames(RNA_C2epimutations_target) <- c('genes','condition','0','2','4','6','8','10','F12','14','16','18','20')
RNA_C2epimutations_target <- subset(RNA_C2epimutations_target,  select = -F12)
RNA_C2epimutations_target<-cbind(RNA_C2epimutations_target, rep(RNA_C2epimutations_target[9],1))
RNA_C2epimutations_target<-RNA_C2epimutations_target[,c(1,2,3,4,5,6,7,8,13,9,10,11,12)]
colnames(RNA_C2epimutations_target) <- c('genes','condition','0','2','4','6','8','10','12','14','16','18','20')

rownames(RNA_C1epimutations_target)<-RNA_C1epimutations_target$genes
RNA_C1epimutations_target <- subset(RNA_C1epimutations_target,  select = -genes)
RNA_C1epimutations_target <- subset(RNA_C1epimutations_target,  select = -condition)

rownames(RNA_C2epimutations_target)<-RNA_C2epimutations_target$genes
RNA_C2epimutations_target <- subset(RNA_C2epimutations_target,  select = -genes)
RNA_C2epimutations_target <- subset(RNA_C2epimutations_target,  select = -condition)

rownames(RNA_L1epimutations_target)<-RNA_L1epimutations_target$genes
RNA_L1epimutations_target <- subset(RNA_L1epimutations_target,  select = -genes)
RNA_L1epimutations_target <- subset(RNA_L1epimutations_target,  select = -condition)
write.xlsx(RNA_L1epimutations_target,"RNA_L1epimutations_target.xlsx")

rownames(RNA_L2epimutations_target)<-RNA_L2epimutations_target$genes
RNA_L2epimutations_target <- subset(RNA_L2epimutations_target,  select = -genes)
RNA_L2epimutations_target <- subset(RNA_L2epimutations_target,  select = -condition)

rownames(RNA_H1epimutations_target)<-RNA_H1epimutations_target$genes
RNA_H1epimutations_target <- subset(RNA_H1epimutations_target,  select = -genes)
RNA_H1epimutations_target <- subset(RNA_H1epimutations_target,  select = -condition)
colnames(RNA_H1epimutations_target) <- c('0','2','4','6','8','10','12','14','16','18','F20')
RNA_H1epimutations_target <- subset(RNA_H1epimutations_target,  select = -F20)

rownames(RNA_H2epimutations_target)<-RNA_H2epimutations_target$genes
RNA_H2epimutations_target <- subset(RNA_H2epimutations_target,  select = -genes)
RNA_H2epimutations_target <- subset(RNA_H2epimutations_target,  select = -condition)

#EnrichR
library(enrichR)

# tell the enrichR package to use worms
setEnrichrSite("WormEnrichr")

# find gene set databases available
dbs <- listEnrichrDbs()

# We need to tell enrichR which databases (from the selection in dbs) we would like to query.
# We can start with KEGG and the three 2018 GO databases

chosendbs <- c("KEGG_2019",
               "GO_Cellular_Component_2018",
               "GO_Molecular_Function_2018",
               "GO_Biological_Process_2018", 
               "InterPro_Domains_2019")

#genes associated to tRNAs vs all genes
# The test gene list
Sim_test_list <- unique(tRNA_target$Name)

# The background list
#Import one sample table to extract genes names
setwd("~/Documents/Cisplatin project/Data_analysis/Count data")
C_elegans_gene_names <- read.xlsx("C_elegans_gene_names2.xlsx")
setwd("~/Documents/Cisplatin project/Data_analysis/Analysis")
Sim_bg_list <- unique(C_elegans_gene_names[,4])

# Send both lists to enrichr
Sim_epimutated_enriched_list  <- enrichr(Sim_test_list, chosendbs)
Sim_epimutated_bg_list  <- enrichr(Sim_bg_list, chosendbs)

Sim_test_bg_list <- list(Sim_epimutated_enriched_list, Sim_epimutated_bg_list)

# Annotate with name of libraries

Sim_test_bg_results <- lapply(Sim_test_bg_list, function(x){
  
  lapply(1:5, function(j){
    
    if(nrow(x[[j]]) > 0){
      x[[j]][1:nrow(x[[j]]), "library"] <- names(Sim_test_bg_list[[1]])[j]
      return(x[[j]])
    } 
    
  })
  
})

# then we combine them and put them in order of significance using the Old.Adjusted.P.value
Sim_test_bg_RESULTS <- lapply(Sim_test_bg_results, function(x){
  
  tempdf <- do.call(rbind, x)
  tempdf[order(tempdf$Old.Adjusted.P.value), ]
  
})

split_value_list <- list()

for(i in 1:length(Sim_test_bg_RESULTS)){
  
  split <- strsplit(Sim_test_bg_RESULTS[[i]][, 2], "/")
  
  split_value <- c()
  
  for(j in 1:length(split)){
    
    save <- split[[j]][1]
    
    split_value <- as.numeric(c(split_value, save))
  }
  
  Sim_test_bg_RESULTS[[i]]$Overlap <- split_value
  
}

# get the test_term, the test_not, the mean(sample_term) and the sample_not

# test_term is overlap
# test_not is length of gene set minus overlap

# then add 1 to all if either test_term OR sample_term = 0
# then only pursue if test_term &/OR sample_term > 4
# let the terms of interest be the terms for the enrichr result for the test list, i.e. with simultaneously epimutated genes

terms_of_interest <- Sim_test_bg_RESULTS[[1]]$Term

# and if the term is not present in the other gene list then assume there is no enrichment for that term in those lists
Test_List <- Sim_test_bg_RESULTS[[1]]
BG_List <- Sim_test_bg_RESULTS[[2]]

table_comp <- c()

for(i in 1:length(terms_of_interest)){
  
  Term_select <- terms_of_interest[[i]]
  
  sample_term <- 0
  
  if(Term_select %in% BG_List$Term == T){
    
    sample_term <- BG_List[which(BG_List$Term %in% Term_select), 2]}
  sample_not <- length(Sim_bg_list) - sample_term
  
  test_term <- Test_List[which(Test_List$Term %in% Term_select), 2]
  test_not <- length(Sim_test_list) - test_term
  
  # save the components in a row in a new table
  
  component_row <- data.frame(test_term, test_not, sample_term, sample_not)
  
  rownames(component_row) <- Term_select
  
  table_comp <- rbind(table_comp, component_row)
  colnames(table_comp) <- c("test_term", "test_not", "sample_term", "sample_not")
  
}

# Now we will adjust  all values by adding 1 if there are any zeros

zeros <- colSums(table_comp == 0)

if(sum(zeros) > 0){
  
  table_comp <-  table_comp + 1
}

# Now modify the table to remove any rows where test term or mean sample term are < 5
# Now apply the function only if at least one category has 5 or more genes
# It is only reasonable to assess for relative enrichment/depletion if either the test list or sample list has > 5 genes in that ontology category

#table_comp <- table_comp[(table_comp$test_term < 2 & table_comp$sample_term < 2)==F, ]

# Now we have found terms where there are at least 5 genes in 1 test or sample category
# we have adjusted the data for each gene set for any 0 values relative to the specific test list
# Now we can do the analysis but only on terms present in each table_comp list

Pval_col <- c()
OR_col <- c()

for(q in 1:nrow(table_comp)){
  
  test_term <- table_comp[q, 1] 
  test_not <- table_comp[q, 2]
  mean_Sample_term <- table_comp[q, 3]
  sample_not <- table_comp[q, 4]
  
  contingency_table <-
    rbind(c(test_term, test_not),
          c(mean_Sample_term, sample_not))
  FT_out <- fisher.test(contingency_table)
  pval <- FT_out$p.value
  OR <- FT_out$estimate
  
  Pval_col <- c(Pval_col, pval)
  OR_col<- c(OR_col, OR)
  
}    

OR_pval_frame <- data.frame(Pval_col, OR_col)

colnames(OR_pval_frame) <- c("Pval", "OR")
rownames(OR_pval_frame) <- rownames(table_comp)

# Bonferroni correction

OR_pval_frame[, 1] <- p.adjust(OR_pval_frame[, 1], method="bonferroni", n = nrow(OR_pval_frame)) 

# Make plots for this

# i) In each plot order them by the largest significant enrichment out of all of the comparisons 

OR_pval_frame_ordered <- OR_pval_frame[order(OR_pval_frame[, 2], decreasing = T), ]

# Make a bubble plot 

P_1 <- as.numeric(OR_pval_frame_ordered[, 1])

OR_1 <- log10(as.numeric(OR_pval_frame_ordered[, 2]))

neg_log_p_1 <- -log(P_1)

bubble_table_1 <- data.frame(rownames(OR_pval_frame_ordered), OR_1,  P_1, neg_log_p_1)

colnames(bubble_table_1) <- c("Term", "log_OR", "p_val", "neg_log_P")

bubble_table_1$Term <- factor(bubble_table_1$Term, levels = rev(bubble_table_1$Term))

find_size <- c()
find_alpha <- c()

for(i in 1:nrow(bubble_table_1)){
  
  alpha <- 0.3
  size_point <- 14
  
  if(as.numeric(bubble_table_1[i, 3]) < 0.1){
    alpha <- 1
    size_point <- bubble_table_1[i,4]
    
  }
  find_alpha <- c(find_alpha, alpha)  
  find_size <- c(find_size, size_point)
}

bubble_table <- cbind(bubble_table_1, find_alpha, find_size)

# just plot first 10

allsignifLD_DEGs<-subset(bubble_table,bubble_table$p_val<0.05)
write.xlsx(allsignifLD_DEGs,"allsignifLD_DEGs.xlsx")
bubble_table <- bubble_table[1:10, ]
trunc_terms <- as.character(bubble_table$Term)

# Trunc terms manually entered as Y axis labels in Adobe Illustrator

trunc_terms <- c(
  
  "positive regulation of organelle assembly",
  "positive regulation of cytoskeleton organization",                                               
  "endosomal transport",                
  "transmembrane transport",                                  
  "actin binding",                            
  "beta-Alanine metabolism",              
  "Mitophagy",                                       
  "ABC transporters",     
  "Propanoate metabolism" ,                 
  "positive regulation of protein complex assembly"
)

guides_merge <- function(gdefs) {
  gdefs <- lapply(gdefs, function(g) { g$hash <- paste(g$order, g$hash, sep = "z"); g})
  tapply(gdefs, sapply(gdefs, function(g)g$hash), function(gs)Reduce(guide_merge, gs))
}
environment(guides_merge) <- environment(ggplot)
assignInNamespace("guides_merge", guides_merge, pos = "package:ggplot2")

# Supplementary Figure 5

bubble_Sim_RNA_chrom <-
  
  ggplot(bubble_table, aes(y=Term, x=as.numeric(log_OR)))+
  geom_point(aes(color=Term, size=as.numeric(find_size), alpha=find_alpha))+
  scale_x_continuous(breaks = c(3.09,3.10,3.11,3.12), limits = c(3.09,3.12)) +
  geom_vline(xintercept = 0, colour = "grey")+
  theme_bw()+
  
  
  scale_size_continuous(range = c(2, 15), breaks = c(10, 5), 
                        limits = c(4, 51))+  
  
  scale_alpha(name = paste("Transparency indicates \nsignificance of enrichment"), range = c(0.3, 1),
              breaks = seq(1, 0.3, length = 2),
              limits = c(0.29, 1.05), 
              labels = c("p value significant < 0.05", "p value not significant > 0.05"))+
  
  ggtitle(paste("Gene ontology of tRNAs associated genes"))

bubble_Sim_RNA_chrom <- bubble_Sim_RNA_chrom + labs(y="Gene Ontology Terms\n", x = "log10(Odds Ratio for enrichment)")+
  theme(axis.title.x=element_text(face = "bold", size=14))+
  theme(axis.title.y=element_text(face = "bold", size=14))+
  theme(plot.title = element_text(face = "bold", size=16))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_discrete(labels= rev(trunc_terms))+
  # scale_y_discrete(labels= "")+
  theme(axis.text.x = element_text(color="#000000", size=10))+
  theme(axis.text.y = element_text(color="#000000", size = 14))+
  guides(size = guide_legend(order = 3), colour = "none", alpha = guide_legend(order = 2))+
  theme(legend.text=element_text(color="#000000", size=12))

bubble_Sim_RNA_chrom

#Get raw data
write.xlsx(bubble_table,"Data_Fig_6_A.xlsx")
#--------------------------
#######Fig.6.B - Overlap tRNAs fragments epimutation & gene exp.change ######### 
#Creation of integrated tables
# We will need piRNA cluster genes for this
# Coordinates for piRNA cluster genes (5 - 7.5 & 13 - 17 MB from Ruby et al. 2006) have been intersected with RNA coordinates 

RNA_piRNA_domains <- read.table("RNA_Coords_intersected_piRNA_c.bed", header = FALSE)
piRNA_domains <- RNA_piRNA_domains[RNA_piRNA_domains$V5==1, ]

# When we are considering the genes with chromatin domain annotations we can only consider the genes which map to Ahringer as this is where we get the domain annotations from. 
all_Ahr_RNA_genes <- unique(Ahringer_single_gene_ref_table$Gene)
piRNA_cluster_genes <- intersect(unique(piRNA_domains$V4), all_Ahr_RNA_genes)

#Control
RNA_list <- list(RNA_C1epimutations_target, RNA_C2epimutations_target)
lin <- c("C1", "C2")

target_RNA_Control_integrated_table <- c()

for(x in 1:length(RNA_list)){
  RNA_bin <- RNA_list[[x]]
  Lineage <- lin[[x]]
  integrated_table <- c()
  ep_table <- c()
  
  for(e in 1:length(rownames(RNA_bin))){
    coord <- rownames(RNA_bin)[e]
    gene <- strsplit(coord, ":")[[1]][4]
    AllUP <- 0
    AllDOWN <- 0
    MixedUPDOWN <- 0
    RNA_mut <- 0
    
    if(sum(abs(RNA_bin[coord, ])) >0){
      RNA_mut <- 1
      # Direction of RNA expression changes  
      # Determine the original direction of the RNA events    
      events <- sum(abs(RNA_bin[coord, ]))
      sum <- sum(RNA_bin[coord, ])
      
      if(sum == events){
        AllUP <- 1
      }     
      
      if(sum == -1*(events)){
        AllDOWN <- 1
      }     
      
      if(!abs(sum) == events){
        MixedUPDOWN <- 1
      }  
    }
    
    RNA_gens <- "0"
    RNA_inherited <- 0
    
    if(RNA_mut ==1){
      RNA_gens<-colnames(RNA_bin[(which(abs(RNA_bin[coord,])>0))])
      
      if(length(RNA_gens) > 1){
        RNA_gens <- paste(RNA_gens, collapse = "_")
        
        if(2 %in% diff(as.numeric(strsplit(RNA_gens, "_")[[1]]))){
          RNA_inherited <- 1
        } 
      }
    }
    save <- data.frame(Lineage, gene, coord, RNA_mut, RNA_gens, RNA_inherited, 
                       AllUP, AllDOWN, MixedUPDOWN)
    integrated_table <- rbind(integrated_table, save)
  }
  target_RNA_Control_integrated_table <- rbind(target_RNA_Control_integrated_table, integrated_table)
}
colnames(target_RNA_Control_integrated_table) <- c(
  "Lineage", "gene", "RNA_coord", "is_RNA_exp_change",
  "RNA_mut_gens", "is_RNA_inherited","RNA_All_UP", "RNA_All_DOWN", "RNA_Mixed_Muts"
)
save(target_RNA_Control_integrated_table,file="target_RNA_Control_integrated_table.Rdata")
write.xlsx(target_RNA_Control_integrated_table,"target_RNA_Control_integrated_table.xlsx")

#Low dose
RNA_list <- list(RNA_L1epimutations_target, RNA_L2epimutations_target)
lin <- c("L1", "L2")

target_RNA_Low_dose_integrated_table <- c()

for(x in 1:length(RNA_list)){
  RNA_bin <- RNA_list[[x]]
  Lineage <- lin[[x]]
  integrated_table <- c()
  ep_table <- c()
  
  for(e in 1:length(rownames(RNA_bin))){
    coord <- rownames(RNA_bin)[e]
    gene <- strsplit(coord, ":")[[1]][4]
    AllUP <- 0
    AllDOWN <- 0
    MixedUPDOWN <- 0
    RNA_mut <- 0
    
    if(sum(abs(RNA_bin[coord, ])) >0){
      RNA_mut <- 1
      # Direction of RNA expression changes  
      # Determine the original direction of the RNA events    
      events <- sum(abs(RNA_bin[coord, ]))
      sum <- sum(RNA_bin[coord, ])
      
      if(sum == events){
        AllUP <- 1
      }     
      
      if(sum == -1*(events)){
        AllDOWN <- 1
      }     
      
      if(!abs(sum) == events){
        MixedUPDOWN <- 1
      }  
    }
    
    RNA_gens <- "0"
    RNA_inherited <- 0
    
    if(RNA_mut ==1){
      RNA_gens<-colnames(RNA_bin[(which(abs(RNA_bin[coord,])>0))])
      
      if(length(RNA_gens) > 1){
        RNA_gens <- paste(RNA_gens, collapse = "_")
        
        if(2 %in% diff(as.numeric(strsplit(RNA_gens, "_")[[1]]))){
          RNA_inherited <- 1
        } 
      }
    }
    save <- data.frame(Lineage, gene, coord, RNA_mut, RNA_gens, RNA_inherited, 
                       AllUP, AllDOWN, MixedUPDOWN)
    integrated_table <- rbind(integrated_table, save)
  }
  target_RNA_Low_dose_integrated_table <- rbind(target_RNA_Low_dose_integrated_table, integrated_table)
}
colnames(target_RNA_Low_dose_integrated_table) <- c("Lineage", "gene", "RNA_coord", "is_RNA_exp_change","RNA_mut_gens", "is_RNA_inherited",  "RNA_All_UP", "RNA_All_DOWN", "RNA_Mixed_Muts")
save(target_RNA_Low_dose_integrated_table,file="target_RNA_Low_dose_integrated_table.Rdata")
write.xlsx(target_RNA_Low_dose_integrated_table,"target_RNA_Low_dose_integrated_table.xlsx")

#High dose
RNA_list <- list(RNA_H1epimutations_target, RNA_H2epimutations_target)
lin <- c("H1", "H2")

target_RNA_High_dose_integrated_table <- c()

for(x in 1:length(RNA_list)){
  RNA_bin <- RNA_list[[x]]
  Lineage <- lin[[x]]
  integrated_table <- c()
  ep_table <- c()
  
  for(e in 1:length(rownames(RNA_bin))){
    coord <- rownames(RNA_bin)[e]
    gene <- strsplit(coord, ":")[[1]][4]
    AllUP <- 0
    AllDOWN <- 0
    MixedUPDOWN <- 0
    RNA_mut <- 0
    
    if(sum(abs(RNA_bin[coord, ])) >0){
      RNA_mut <- 1
      # Direction of RNA expression changes  
      # Determine the original direction of the RNA events    
      events <- sum(abs(RNA_bin[coord, ]))
      sum <- sum(RNA_bin[coord, ])
      
      if(sum == events){
        AllUP <- 1
      }     
      
      if(sum == -1*(events)){
        AllDOWN <- 1
      }     
      
      if(!abs(sum) == events){
        MixedUPDOWN <- 1
      }  
    }
    
    RNA_gens <- "0"
    RNA_inherited <- 0
    
    if(RNA_mut ==1){
      RNA_gens<-colnames(RNA_bin[(which(abs(RNA_bin[coord,])>0))])
      if(length(RNA_gens) > 1){
        RNA_gens <- paste(RNA_gens, collapse = "_")
        
        if(2 %in% diff(as.numeric(strsplit(RNA_gens, "_")[[1]]))){
          RNA_inherited <- 1
        } 
      }
    }
    save <- data.frame(Lineage, gene, coord, RNA_mut, RNA_gens, RNA_inherited, 
                       AllUP, AllDOWN, MixedUPDOWN)
    integrated_table <- rbind(integrated_table, save)
  }
  target_RNA_High_dose_integrated_table <- rbind(target_RNA_High_dose_integrated_table, integrated_table)
}
colnames(target_RNA_High_dose_integrated_table) <- c("Lineage", "gene", "RNA_coord", "is_RNA_exp_change","RNA_mut_gens", "is_RNA_inherited",  "RNA_All_UP", "RNA_All_DOWN", "RNA_Mixed_Muts")
save(target_RNA_High_dose_integrated_table,file="target_RNA_High_dose_integrated_table.Rdata")
write.xlsx(target_RNA_High_dose_integrated_table,"target_RNA_High_dose_integrated_table.xlsx")

#Overlaps calculation

#Control
#We kept only the genes associated with tRNAs up to 2 mismatchs
target_RNA_tRNA_Control_integrated_table<-read.xlsx("target_RNA_&_tRNA_Control_integrated_table2.xlsx")
#Inherited RNA seq changes
#Nb genes with target RNA inherited and tRNA epimut
Control_genes_target_RNA_muts_inherited_and_tRNA_muts <- nrow(target_RNA_tRNA_Control_integrated_table[target_RNA_tRNA_Control_integrated_table$is_RNA_exp_change==1&target_RNA_tRNA_Control_integrated_table$is_RNA_inherited==1&target_RNA_tRNA_Control_integrated_table$tRNA_mut==1, ])
#Nb genes with both RNA and 22G inherited epimut at matching time point
Control_target_RNA_inherit_tRNA_mut_tm_inherit <- nrow(target_RNA_tRNA_Control_integrated_table[target_RNA_tRNA_Control_integrated_table$is_RNA_exp_change==1&target_RNA_tRNA_Control_integrated_table$is_RNA_inherited==1&target_RNA_tRNA_Control_integrated_table$tRNA_mut==1&target_RNA_tRNA_Control_integrated_table$tRNA_inherited==1&target_RNA_tRNA_Control_integrated_table$Time_match_to_tRNA==1, ])/Control_genes_target_RNA_muts_inherited_and_tRNA_muts*100
# genes with RNA inherit and 22G non inherited epimut at matching time point
Control_target_RNA_inherit_tRNA_mut_tm_non_inherit <- nrow(target_RNA_tRNA_Control_integrated_table[target_RNA_tRNA_Control_integrated_table$is_RNA_exp_change==1&target_RNA_tRNA_Control_integrated_table$is_RNA_inherited==1&target_RNA_tRNA_Control_integrated_table$tRNA_mut==1&target_RNA_tRNA_Control_integrated_table$tRNA_inherited==0&target_RNA_tRNA_Control_integrated_table$Time_match_to_tRNA==1, ])/Control_genes_target_RNA_muts_inherited_and_tRNA_muts*100
#% genes with RNA inherit and 22G non inherited epimut at non matching time point
Control_target_RNA_inherit_tRNA_mut_non_tm <-  nrow(target_RNA_tRNA_Control_integrated_table[target_RNA_tRNA_Control_integrated_table$is_RNA_exp_change==1&target_RNA_tRNA_Control_integrated_table$is_RNA_inherited==1&target_RNA_tRNA_Control_integrated_table$tRNA_mut==1&target_RNA_tRNA_Control_integrated_table$Time_match_to_tRNA==0, ])/Control_genes_target_RNA_muts_inherited_and_tRNA_muts*100

Control_target_RNA_Inherit_proportions <-  data.frame(c(Control_target_RNA_inherit_tRNA_mut_non_tm, 
                                                        Control_target_RNA_inherit_tRNA_mut_tm_non_inherit,
                                                        Control_target_RNA_inherit_tRNA_mut_tm_inherit))
Names <- c(
  "Has tRNA epimutations - not time matched", 
  "Has time matched tRNA epimutations - not inherited", 
  "Has time matched tRNA epimutations - inherited")                                           
Species <- rep("Genes with inherited \nRNA expression changes\ntargeted by tRNAs epimutations", 3)

Control_target_RNA_Inherit <- cbind(Control_target_RNA_Inherit_proportions, Names, Species)
colnames(Control_target_RNA_Inherit) <- c("Percentages", "Epigenetic_feature_of_gene", "Species")

#Non Inherited RNA seq changes

Control_genes_target_RNA_muts_n_inherited_tRNA_muts <- nrow(target_RNA_tRNA_Control_integrated_table[target_RNA_tRNA_Control_integrated_table$is_RNA_exp_change==1&target_RNA_tRNA_Control_integrated_table$is_RNA_inherited==0&target_RNA_tRNA_Control_integrated_table$tRNA_mut==1, ])
Control_target_RNA_n_inherit_tRNA_mut_tm_inherit <- nrow(target_RNA_tRNA_Control_integrated_table[target_RNA_tRNA_Control_integrated_table$is_RNA_exp_change==1&target_RNA_tRNA_Control_integrated_table$is_RNA_inherited==0&target_RNA_tRNA_Control_integrated_table$tRNA_mut==1&target_RNA_tRNA_Control_integrated_table$tRNA_inherited==1&target_RNA_tRNA_Control_integrated_table$Time_match_to_tRNA==1, ])/Control_genes_target_RNA_muts_n_inherited_tRNA_muts*100
Control_target_RNA_n_inherit_tRNA_mut_tm_non_inherit <- nrow(target_RNA_tRNA_Control_integrated_table[target_RNA_tRNA_Control_integrated_table$is_RNA_exp_change==1&target_RNA_tRNA_Control_integrated_table$is_RNA_inherited==0&target_RNA_tRNA_Control_integrated_table$tRNA_mut==1&target_RNA_tRNA_Control_integrated_table$tRNA_inherited==0&target_RNA_tRNA_Control_integrated_table$Time_match_to_tRNA==1, ])/Control_genes_target_RNA_muts_n_inherited_tRNA_muts*100
Control_target_RNA_n_inherit_tRNA_mut_non_tm <-  nrow(target_RNA_tRNA_Control_integrated_table[target_RNA_tRNA_Control_integrated_table$is_RNA_exp_change==1&target_RNA_tRNA_Control_integrated_table$is_RNA_inherited==0&target_RNA_tRNA_Control_integrated_table$tRNA_mut==1&target_RNA_tRNA_Control_integrated_table$Time_match_to_tRNA==0, ])/Control_genes_target_RNA_muts_n_inherited_tRNA_muts*100

Control_target_RNA_n_Inherit_proportions <-  data.frame(c(Control_target_RNA_n_inherit_tRNA_mut_non_tm, 
                                                          Control_target_RNA_n_inherit_tRNA_mut_tm_non_inherit,
                                                          Control_target_RNA_n_inherit_tRNA_mut_tm_inherit))


Names <- c(
  "Has tRNA epimutations - not time matched", 
  "Has time matched tRNA epimutations - not inherited", 
  "Has time matched tRNA epimutations - inherited")                                             

Species <- rep("Genes with non inherited \nRNA expression changes\ntargeted by tRNAs epimutations", 3)

Control_target_RNA_N_Inherit <- cbind(Control_target_RNA_n_Inherit_proportions, Names, Species)

colnames(Control_target_RNA_N_Inherit) <- c("Percentages", "Epigenetic_feature_of_gene", "Species")

#Plot for Control

Control_target_RNA_Inherited_features <- rbind(Control_target_RNA_N_Inherit, Control_target_RNA_Inherit)
Control_target_RNA_Inherited_features <- Control_target_RNA_Inherited_features[order(nrow(Control_target_RNA_Inherited_features):1), ] 

Control_target_RNA_Inherited_features$Epigenetic_feature_of_gene <- factor(Control_target_RNA_Inherited_features$Epigenetic_feature_of_gene, levels = unique(Control_target_RNA_Inherited_features$Epigenetic_feature_of_gene))
positions <- unique(Control_target_RNA_Inherited_features$Species)

# Proportions of genes with inherited and non inherited RNAseq changes that have small RNA changes

Control_target_RNA_w_tRNA_proportion_plot <-
  
  ggplot(Control_target_RNA_Inherited_features, aes(fill=Epigenetic_feature_of_gene, y=Percentages, x=Species)) + 
  
  geom_bar(position="stack", stat="identity", width = 0.4, color="black")+
  scale_x_discrete(limits = positions)+
  scale_fill_manual(values=c("cornflowerblue", 
                             "cornsilk3", 
                             "coral",
                             "cornflowerblue", 
                             "cornsilk3", 
                             "coral"))+
  theme_bw()+
  labs(x= " ", 
       y = "Percentage of genes (%)\n")+
  theme(plot.margin = unit(c(2,1,2,2), "cm"))+
  theme(plot.title = element_text(size = 16))+
  labs(title = "Control data")+
  theme(text = element_text(size=14))+
  theme(axis.text.x=element_text(colour="black", size = 2))+
  theme(axis.text.y=element_text(colour="black", size = 16))+
  guides(fill=guide_legend(ncol=1))+
  scale_fill_manual(name = "RNA expression change status",
                    values=c("cornflowerblue", 
                             "cornsilk3", 
                             "coral"),
                    labels=c("Simultaneous inherited tRNA epimutation", 
                             "Simultaneous non inherited tRNA epimutation", 
                             "Non simultaneous tRNA epimutation"))

Control_target_RNA_w_tRNA_proportion_plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(fill=guide_legend(title="Associated tRNAs epimutation status"))+
  theme(legend.text = element_text(size=16))+
  theme(legend.title = element_text(size=16))

#Get the raw data
write.xlsx(Control_target_RNA_Inherited_features,"Data_Fig_6_B_Control.xlsx")

#Low dose
target_RNA_tRNAs_Low_dose_integrated_table<-read.xlsx("target_RNA_&_tRNA_Low_dose_integrated_table2.xlsx")
#Inherited RNA seq changes
Low_dose_genes_target_RNA_muts_inherited_tRNA_muts <- nrow(target_RNA_tRNAs_Low_dose_integrated_table[target_RNA_tRNAs_Low_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_Low_dose_integrated_table$is_RNA_inherited==1&target_RNA_tRNAs_Low_dose_integrated_table$tRNA_mut==1, ])
Low_dose_target_RNA_inherit_tRNA_mut_tm_inherit <- nrow(target_RNA_tRNAs_Low_dose_integrated_table[target_RNA_tRNAs_Low_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_Low_dose_integrated_table$is_RNA_inherited==1&target_RNA_tRNAs_Low_dose_integrated_table$tRNA_mut==1&target_RNA_tRNAs_Low_dose_integrated_table$tRNA_inherited==1&target_RNA_tRNAs_Low_dose_integrated_table$Time_match_to_tRNA==1, ])/Low_dose_genes_target_RNA_muts_inherited_tRNA_muts*100

Low_dose_target_RNA_inherit_tRNA_mut_tm_non_inherit <- nrow(target_RNA_tRNAs_Low_dose_integrated_table[target_RNA_tRNAs_Low_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_Low_dose_integrated_table$is_RNA_inherited==1&target_RNA_tRNAs_Low_dose_integrated_table$tRNA_mut==1&target_RNA_tRNAs_Low_dose_integrated_table$tRNA_inherited==0&target_RNA_tRNAs_Low_dose_integrated_table$Time_match_to_tRNA==1, ])/Low_dose_genes_target_RNA_muts_inherited_tRNA_muts*100
Low_dose_target_RNA_inherit_tRNA_mut_non_tm <-  nrow(target_RNA_tRNAs_Low_dose_integrated_table[target_RNA_tRNAs_Low_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_Low_dose_integrated_table$is_RNA_inherited==1&target_RNA_tRNAs_Low_dose_integrated_table$tRNA_mut==1&target_RNA_tRNAs_Low_dose_integrated_table$Time_match_to_tRNA==0, ])/Low_dose_genes_target_RNA_muts_inherited_tRNA_muts*100

Low_dose_target_RNA_Inherit_proportions <-  data.frame(c(Low_dose_target_RNA_inherit_tRNA_mut_non_tm, 
                                                         Low_dose_target_RNA_inherit_tRNA_mut_tm_non_inherit,
                                                         Low_dose_target_RNA_inherit_tRNA_mut_tm_inherit))
Names <- c(
  "Has tRNA epimutations - not time matched", 
  "Has time matched tRNA epimutations - not inherited", 
  "Has time matched tRNA epimutations - inherited")                                           
Species <- rep("Genes with inherited \nRNA expression changes\ntargeted by tRNA epimutations", 3)

Low_dose_target_RNA_Inherit <- cbind(Low_dose_target_RNA_Inherit_proportions, Names, Species)
colnames(Low_dose_target_RNA_Inherit) <- c("Percentages", "Epigenetic_feature_of_gene", "Species")

#Non Inherited RNA seq changes
Low_dose_genes_target_RNA_muts_n_inherited_tRNA_muts <- nrow(target_RNA_tRNAs_Low_dose_integrated_table[target_RNA_tRNAs_Low_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_Low_dose_integrated_table$is_RNA_inherited==0&target_RNA_tRNAs_Low_dose_integrated_table$tRNA_mut==1, ])
Low_dose_target_RNA_n_inherit_tRNA_mut_tm_inherit <- nrow(target_RNA_tRNAs_Low_dose_integrated_table[target_RNA_tRNAs_Low_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_Low_dose_integrated_table$is_RNA_inherited==0&target_RNA_tRNAs_Low_dose_integrated_table$tRNA_mut==1&target_RNA_tRNAs_Low_dose_integrated_table$tRNA_inherited==1&target_RNA_tRNAs_Low_dose_integrated_table$Time_match_to_tRNA==1, ])/Low_dose_genes_target_RNA_muts_n_inherited_tRNA_muts*100

Low_dose_target_RNA_n_inherit_tRNA_mut_tm_non_inherit <- nrow(target_RNA_tRNAs_Low_dose_integrated_table[target_RNA_tRNAs_Low_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_Low_dose_integrated_table$is_RNA_inherited==0&target_RNA_tRNAs_Low_dose_integrated_table$tRNA_mut==1&target_RNA_tRNAs_Low_dose_integrated_table$tRNA_inherited==0&target_RNA_tRNAs_Low_dose_integrated_table$Time_match_to_tRNA==1, ])/Low_dose_genes_target_RNA_muts_n_inherited_tRNA_muts*100
Low_dose_target_RNA_n_inherit_tRNA_mut_non_tm <-  nrow(target_RNA_tRNAs_Low_dose_integrated_table[target_RNA_tRNAs_Low_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_Low_dose_integrated_table$is_RNA_inherited==0&target_RNA_tRNAs_Low_dose_integrated_table$tRNA_mut==1&target_RNA_tRNAs_Low_dose_integrated_table$Time_match_to_tRNA==0, ])/Low_dose_genes_target_RNA_muts_n_inherited_tRNA_muts*100

Low_dose_target_RNA_n_Inherit_proportions <-  data.frame(c(Low_dose_target_RNA_n_inherit_tRNA_mut_non_tm, 
                                                           Low_dose_target_RNA_n_inherit_tRNA_mut_tm_non_inherit,
                                                           Low_dose_target_RNA_n_inherit_tRNA_mut_tm_inherit))

Names <- c(
  "Has tRNA epimutations - not time matched", 
  "Has time matched tRNA epimutations - not inherited", 
  "Has time matched tRNA epimutations - inherited")                                             

Species <- rep("Genes with non inherited \nRNA expression changes\ntargeted by tRNA epimutations", 3)

Low_dose_target_RNA_N_Inherit <- cbind(Low_dose_target_RNA_n_Inherit_proportions, Names, Species)

colnames(Low_dose_target_RNA_N_Inherit) <- c("Percentages", "Epigenetic_feature_of_gene", "Species")

#Plot for Low dose

Low_dose_target_RNA_Inherited_features <- rbind(Low_dose_target_RNA_N_Inherit, Low_dose_target_RNA_Inherit)
Low_dose_target_RNA_Inherited_features <- Low_dose_target_RNA_Inherited_features[order(nrow(Low_dose_target_RNA_Inherited_features):1), ] 

Low_dose_target_RNA_Inherited_features$Epigenetic_feature_of_gene <- factor(Low_dose_target_RNA_Inherited_features$Epigenetic_feature_of_gene, levels = unique(Low_dose_target_RNA_Inherited_features$Epigenetic_feature_of_gene))
positions <- unique(Low_dose_target_RNA_Inherited_features$Species)

# Proportions of genes with inherited and non inherited RNAseq changes that have small RNA changes

Low_dose_target_RNA_w_smallRNA_proportion_plot <-
  
  ggplot(Low_dose_target_RNA_Inherited_features, aes(fill=Epigenetic_feature_of_gene, y=Percentages, x=Species)) + 
  
  geom_bar(position="stack", stat="identity", width = 0.4, color="black")+
  scale_x_discrete(limits = positions)+
  scale_fill_manual(values=c("cornflowerblue", 
                             "cornsilk3", 
                             "coral",
                             "cornflowerblue", 
                             "cornsilk3", 
                             "coral"))+
  theme_bw()+
  labs(x= " ", 
       y = "Percentage of genes (%)\n")+
  theme(plot.margin = unit(c(2,1,2,2), "cm"))+
  theme(plot.title = element_text(size = 16))+
  labs(title = "Low dose data")+
  theme(text = element_text(size=14))+
  theme(axis.text.x=element_text(colour="black", size = 2))+
  theme(axis.text.y=element_text(colour="black", size = 16))+
  guides(fill=guide_legend(ncol=1))+
  scale_fill_manual(name = "RNA expression change status",
                    values=c("cornflowerblue", 
                             "cornsilk3", 
                             "coral"),
                    labels=c("Simultaneous inherited tRNA epimutation", 
                             "Simultaneous non inherited tRNA epimutation", 
                             "Non simultaneous tRNA epimutation"))

Low_dose_target_RNA_w_smallRNA_proportion_plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(fill=guide_legend(title="Associated tRNA epimutation status"))+
  theme(legend.text = element_text(size=16))+
  theme(legend.title = element_text(size=16))

#Get raw data
write.xlsx(Low_dose_target_RNA_Inherited_features,"Data_Fig_6_B_Low_dose.xlsx")

#High dose
target_RNA_tRNAs_High_dose_integrated_table<-read.xlsx("target_RNA_&_tRNA_High_dose_integrated_table2.xlsx")
#Inherited RNA seq changes
High_dose_genes_RNA_muts_inherited_smallRNA_muts <- nrow(target_RNA_tRNAs_High_dose_integrated_table[target_RNA_tRNAs_High_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_High_dose_integrated_table$is_RNA_inherited==1&target_RNA_tRNAs_High_dose_integrated_table$tRNA_mut==1, ])
High_dose_RNA_inherit_smallRNA_mut_tm_inherit <- nrow(target_RNA_tRNAs_High_dose_integrated_table[target_RNA_tRNAs_High_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_High_dose_integrated_table$is_RNA_inherited==1&target_RNA_tRNAs_High_dose_integrated_table$tRNA_mut==1&target_RNA_tRNAs_High_dose_integrated_table$tRNA_inherited==1&target_RNA_tRNAs_High_dose_integrated_table$Time_match_to_tRNA==1, ])/High_dose_genes_RNA_muts_inherited_smallRNA_muts*100

High_dose_RNA_inherit_smallRNA_mut_tm_non_inherit <- nrow(target_RNA_tRNAs_High_dose_integrated_table[target_RNA_tRNAs_High_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_High_dose_integrated_table$is_RNA_inherited==1&target_RNA_tRNAs_High_dose_integrated_table$tRNA_mut==1&target_RNA_tRNAs_High_dose_integrated_table$tRNA_inherited==0&target_RNA_tRNAs_High_dose_integrated_table$Time_match_to_tRNA==1, ])/High_dose_genes_RNA_muts_inherited_smallRNA_muts*100
High_dose_RNA_inherit_smallRNA_mut_non_tm <-  nrow(target_RNA_tRNAs_High_dose_integrated_table[target_RNA_tRNAs_High_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_High_dose_integrated_table$is_RNA_inherited==1&target_RNA_tRNAs_High_dose_integrated_table$tRNA_mut==1&target_RNA_tRNAs_High_dose_integrated_table$Time_match_to_tRNA==0, ])/High_dose_genes_RNA_muts_inherited_smallRNA_muts*100

High_dose_RNA_Inherit_proportions <-  data.frame(c(High_dose_RNA_inherit_smallRNA_mut_non_tm, 
                                                   High_dose_RNA_inherit_smallRNA_mut_tm_non_inherit,
                                                   High_dose_RNA_inherit_smallRNA_mut_tm_inherit))
Names <- c(
  "Has small RNA epimutations - not time matched", 
  "Has time matched small RNA epimutations - not inherited", 
  "Has time matched small RNA epimutations - inherited")                                           
Species <- rep("Genes with inherited \nRNA expression changes\ntargeted by tRNA epimutations", 3)

High_dose_RNA_Inherit <- cbind(High_dose_RNA_Inherit_proportions, Names, Species)
colnames(High_dose_RNA_Inherit) <- c("Percentages", "Epigenetic_feature_of_gene", "Species")

#Non Inherited RNA seq changes

High_dose_genes_RNA_muts_n_inherited_smallRNA_muts <- nrow(target_RNA_tRNAs_High_dose_integrated_table[target_RNA_tRNAs_High_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_High_dose_integrated_table$is_RNA_inherited==0&target_RNA_tRNAs_High_dose_integrated_table$tRNA_mut==1, ])
High_dose_RNA_n_inherit_smallRNA_mut_tm_inherit <- nrow(target_RNA_tRNAs_High_dose_integrated_table[target_RNA_tRNAs_High_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_High_dose_integrated_table$is_RNA_inherited==0&target_RNA_tRNAs_High_dose_integrated_table$tRNA_mut==1&target_RNA_tRNAs_High_dose_integrated_table$tRNA_inherited==1&target_RNA_tRNAs_High_dose_integrated_table$Time_match_to_tRNA==1, ])/High_dose_genes_RNA_muts_n_inherited_smallRNA_muts*100

High_dose_RNA_n_inherit_smallRNA_mut_tm_non_inherit <- nrow(target_RNA_tRNAs_High_dose_integrated_table[target_RNA_tRNAs_High_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_High_dose_integrated_table$is_RNA_inherited==0&target_RNA_tRNAs_High_dose_integrated_table$tRNA_mut==1&target_RNA_tRNAs_High_dose_integrated_table$tRNA_inherited==0&target_RNA_tRNAs_High_dose_integrated_table$Time_match_to_tRNA==1, ])/High_dose_genes_RNA_muts_n_inherited_smallRNA_muts*100
High_dose_RNA_n_inherit_smallRNA_mut_non_tm <-  nrow(target_RNA_tRNAs_High_dose_integrated_table[target_RNA_tRNAs_High_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_High_dose_integrated_table$is_RNA_inherited==0&target_RNA_tRNAs_High_dose_integrated_table$tRNA_mut==1&target_RNA_tRNAs_High_dose_integrated_table$Time_match_to_tRNA==0, ])/High_dose_genes_RNA_muts_n_inherited_smallRNA_muts*100

High_dose_RNA_n_Inherit_proportions <-  data.frame(c(High_dose_RNA_n_inherit_smallRNA_mut_non_tm, 
                                                     High_dose_RNA_n_inherit_smallRNA_mut_tm_non_inherit,
                                                     High_dose_RNA_n_inherit_smallRNA_mut_tm_inherit))

Names <- c(
  "Has small RNA epimutations - not time matched", 
  "Has time matched small RNA epimutations - not inherited", 
  "Has time matched small RNA epimutations - inherited")                                             

Species <- rep("Genes with non inherited \nRNA expression changes\ntargeted by tRNA epimutations", 3)

High_dose_RNA_N_Inherit <- cbind(High_dose_RNA_n_Inherit_proportions, Names, Species)

colnames(High_dose_RNA_N_Inherit) <- c("Percentages", "Epigenetic_feature_of_gene", "Species")

#Plot for High dose

High_dose_RNA_Inherited_features <- rbind(High_dose_RNA_N_Inherit, High_dose_RNA_Inherit)
High_dose_RNA_Inherited_features <- High_dose_RNA_Inherited_features[order(nrow(High_dose_RNA_Inherited_features):1), ] 

High_dose_RNA_Inherited_features$Epigenetic_feature_of_gene <- factor(High_dose_RNA_Inherited_features$Epigenetic_feature_of_gene, levels = unique(High_dose_RNA_Inherited_features$Epigenetic_feature_of_gene))
positions <- unique(High_dose_RNA_Inherited_features$Species)

# Proportions of genes with inherited and non inherited RNAseq changes that have small RNA changes

High_dose_RNA_w_smallRNA_proportion_plot <-
  
  ggplot(High_dose_RNA_Inherited_features, aes(fill=Epigenetic_feature_of_gene, y=Percentages, x=Species)) + 
  
  geom_bar(position="stack", stat="identity", width = 0.4, color="black")+
  scale_x_discrete(limits = positions)+
  scale_fill_manual(values=c("cornflowerblue", 
                             "cornsilk3", 
                             "coral",
                             "cornflowerblue", 
                             "cornsilk3", 
                             "coral"))+
  theme_bw()+
  labs(x= " ", 
       y = "Percentage of genes (%)\n")+
  theme(plot.margin = unit(c(2,1,2,2), "cm"))+
  theme(plot.title = element_text(size = 16))+
  labs(title = "High dose data")+
  theme(text = element_text(size=14))+
  theme(axis.text.x=element_text(colour="black", size = 2))+
  theme(axis.text.y=element_text(colour="black", size = 16))+
  guides(fill=guide_legend(ncol=1))+
  scale_fill_manual(name = "RNA expression change status",
                    values=c("cornflowerblue", 
                             "cornsilk3", 
                             "coral"),
                    labels=c("Simultaneous inherited tRNA epimutation", 
                             "Simultaneous non inherited tRNA epimutation", 
                             "Non simultaneous tRNA epimutation"))

High_dose_RNA_w_smallRNA_proportion_plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(fill=guide_legend(title="Associated tRNA epimutation status"))+
  theme(legend.text = element_text(size=16))+
  theme(legend.title = element_text(size=16))

#Get raw data
write.xlsx(High_dose_RNA_Inherited_features,"Data_Fig_6_B_High_dose.xlsx")

#-------------
#######Table 3 - Odds of association tRNAs fragments epimutation & gene exp.change#####
#Odds calculation Control
# 1. What are the odds that genes have both tRNA epimutations and RNA changes ?

# background is all genes which have tRNA signal, so have to map to small RNA

Control_background <- nrow(target_RNA_tRNA_Control_integrated_table)

# of which have RNA changes

Control_with_target_RNA_changes <- nrow(target_RNA_tRNA_Control_integrated_table[target_RNA_tRNA_Control_integrated_table$is_RNA_exp_change==1, ])

# of which have tRNA changes 

Control_with_tRNA_changes <- nrow(target_RNA_tRNA_Control_integrated_table[target_RNA_tRNA_Control_integrated_table$tRNA_mut==1, ])

# of which have both 

Control_with_both <- nrow(target_RNA_tRNA_Control_integrated_table[target_RNA_tRNA_Control_integrated_table$tRNA_mut==1&target_RNA_tRNA_Control_integrated_table$is_RNA_exp_change==1, ])

Control_tRNA_RNA <-  Control_with_both

Control_RNA_not_tRNA <- Control_with_target_RNA_changes - Control_tRNA_RNA

Control_tRNA_not_RNA <- Control_with_tRNA_changes  -  Control_tRNA_RNA

Control_not_tRNA_not_RNA <- Control_background - (Control_with_target_RNA_changes + Control_with_tRNA_changes -  Control_tRNA_RNA)

contingency_table <-
  rbind(c(Control_tRNA_RNA, Control_RNA_not_tRNA),
        c(Control_tRNA_not_RNA, Control_not_tRNA_not_RNA))
FT_out <- fisher.test(contingency_table)
pval_1 <- FT_out$p.value
OddsRatio_1 <- FT_out$estimate

# 2. What are the odds that genes with RNA expression changes have both INHERITED RNA expression changes and small RNA epimutations?

# background is all genes which map to small RNA and have AND RNA expression changes 

background <- nrow(target_RNA_tRNA_Control_integrated_table[target_RNA_tRNA_Control_integrated_table$is_RNA_exp_change==1, ])

# of which have inherited RNA changes

with_inh_RNA_changes <- nrow(target_RNA_tRNA_Control_integrated_table[target_RNA_tRNA_Control_integrated_table$is_RNA_exp_change==1&target_RNA_tRNA_Control_integrated_table$is_RNA_inherited==1, ])

# of which have tRNA changes of any kind

with_tRNA_changes <- nrow(target_RNA_tRNA_Control_integrated_table[target_RNA_tRNA_Control_integrated_table$is_RNA_exp_change==1&target_RNA_tRNA_Control_integrated_table$tRNA_mut==1, ])

# of which have both inherited RNA changes and small RNA changes

with_both <- nrow(target_RNA_tRNA_Control_integrated_table[target_RNA_tRNA_Control_integrated_table$tRNA_mut==1&target_RNA_tRNA_Control_integrated_table$is_RNA_exp_change==1&target_RNA_tRNA_Control_integrated_table$is_RNA_inherited==1, ])
with_both <- nrow(target_RNA_tRNA_Control_integrated_table[target_RNA_tRNA_Control_integrated_table$tRNA_mut==1&target_RNA_tRNA_Control_integrated_table$is_RNA_exp_change==1&target_RNA_tRNA_Control_integrated_table$Time_match_to_tRNA==1, ])

tRNA_RNA <-  with_both

RNA_not_tRNA <- with_inh_RNA_changes - tRNA_RNA

tRNA_not_RNA <- with_tRNA_changes  -  tRNA_RNA

not_tRNA_not_RNA <- background - (with_inh_RNA_changes + with_tRNA_changes -  tRNA_RNA)


contingency_table <-
  rbind(c(tRNA_RNA, RNA_not_tRNA),
        c(tRNA_not_RNA, not_tRNA_not_RNA))
FT_out <- fisher.test(contingency_table)
pval_2 <- FT_out$p.value
OddsRatio_2 <- FT_out$estimate

contingency_table <-
  rbind(c(75, 912),
        c(6, 25))
FT_out <- fisher.test(contingency_table)
pval_2 <- FT_out$p.value
OddsRatio_2 <- FT_out$estimate
# Time matching

# 3. What are the odds that Genes with RNA expression changes have both inherited RNA expression changes and time matched small RNA epimutations?


# background is all genes which map to small RNA and have RNA expression changes

background <- nrow(target_RNA_tRNA_Control_integrated_table[target_RNA_tRNA_Control_integrated_table$is_RNA_exp_change==1, ])

# of which have inherited RNA changes

with_RNA_changes <- nrow(target_RNA_tRNA_Control_integrated_table[target_RNA_tRNA_Control_integrated_table$is_RNA_exp_change==1&target_RNA_tRNA_Control_integrated_table$is_RNA_inherited==1, ])

# of which have time matched small RNA changes 

with_tRNA_changes <- nrow(target_RNA_tRNA_Control_integrated_table[target_RNA_tRNA_Control_integrated_table$is_RNA_exp_change==1&target_RNA_tRNA_Control_integrated_table$tRNA_mut==1&target_RNA_tRNA_Control_integrated_table$Time_match_to_tRNA ==1, ])


# of which have both inherited RNA changes and time matched small RNA changes

with_both <- nrow(target_RNA_tRNA_Control_integrated_table[target_RNA_tRNA_Control_integrated_table$is_RNA_exp_change==1&target_RNA_tRNA_Control_integrated_table$is_RNA_inherited==1&target_RNA_tRNA_Control_integrated_table$tRNA_mut==1&target_RNA_tRNA_Control_integrated_table$Time_match_to_tRNA==1, ])
tRNA__RNA <-  with_both

RNA_not_tRNA <- with_RNA_changes - tRNA_RNA

tRNA_not_RNA <- with_tRNA_changes  -  tRNA_RNA

not_tRNA_not_RNA <- background - (with_RNA_changes + with_tRNA_changes -  tRNA_RNA)


contingency_table <-
  rbind(c(tRNA__RNA, RNA_not_tRNA),
        c(tRNA_not_RNA, not_tRNA_not_RNA))
FT_out <- fisher.test(contingency_table)
pval_3 <- FT_out$p.value
OddsRatio_3 <- FT_out$estimate


# 4. What are the odds that in Genes with RNA expression changes and small RNA epimutations, the small RNA epimutations are time matched and the RNA expression changes are inherited?

# background is all genes which have small RNA data points AND RNA expression changes and small RNA epimutations

background <- nrow(target_RNA_tRNA_Control_integrated_table[target_RNA_tRNA_Control_integrated_table$is_RNA_exp_change==1&target_RNA_tRNA_Control_integrated_table$tRNA_mut==1, ])

# of which have inherited RNA changes

with_RNA_changes <- nrow(target_RNA_tRNA_Control_integrated_table[target_RNA_tRNA_Control_integrated_table$tRNA_mut==1&target_RNA_tRNA_Control_integrated_table$is_RNA_exp_change==1&target_RNA_tRNA_Control_integrated_table$is_RNA_inherited==1, ])

# of which have time matched small RNA changes

with_tRNA_changes <- nrow(target_RNA_tRNA_Control_integrated_table[target_RNA_tRNA_Control_integrated_table$tRNA_mut==1&target_RNA_tRNA_Control_integrated_table$is_RNA_exp_change==1&target_RNA_tRNA_Control_integrated_table$Time_match_to_tRNA==1, ])

# of which have both inherited RNA changes and time matched smallRNA changes

with_both <- nrow(target_RNA_tRNA_Control_integrated_table[target_RNA_tRNA_Control_integrated_table$tRNA_mut==1&target_RNA_tRNA_Control_integrated_table$is_RNA_exp_change==1&target_RNA_tRNA_Control_integrated_table$is_RNA_inherited==1&target_RNA_tRNA_Control_integrated_table$Time_match_to_tRNA==1, ])

tRNA_RNA <-  with_both

RNA_not_tRNA <- with_RNA_changes - tRNA_RNA

tRNA_not_RNA <- with_tRNA_changes  -  tRNA_RNA

not_tRNA_not_RNA <- background - (with_RNA_changes + with_tRNA_changes -  tRNA_RNA)


contingency_table <-
  rbind(c(tRNA_RNA, RNA_not_tRNA),
        c(tRNA_not_RNA, not_tRNA_not_RNA))
FT_out <- fisher.test(contingency_table)
pval_4 <- FT_out$p.value
OddsRatio_4 <- FT_out$estimate

#5. What are the odds that genes have both inherited small RNA epimutations and inherited RNA changes ?
# background is all genes which have small RNA data points AND RNA expression changes and small RNA epimutations

background <- nrow(target_RNA_tRNA_Control_integrated_table[target_RNA_tRNA_Control_integrated_table$is_RNA_exp_change==1&target_RNA_tRNA_Control_integrated_table$tRNA_mut==1, ])

# of which have inherited RNA changes

with_RNA_changes <- nrow(target_RNA_tRNA_Control_integrated_table[target_RNA_tRNA_Control_integrated_table$tRNA_mut==1&target_RNA_tRNA_Control_integrated_table$is_RNA_exp_change==1&target_RNA_tRNA_Control_integrated_table$is_RNA_inherited==1, ])

# of which have inherited small RNA changes

with_tRNA_changes <- nrow(target_RNA_tRNA_Control_integrated_table[target_RNA_tRNA_Control_integrated_table$tRNA_mut==1&target_RNA_tRNA_Control_integrated_table$is_RNA_exp_change==1&target_RNA_tRNA_Control_integrated_table$tRNA_inherited==1, ])


# of which have both inherited RNA changes and time matched smallRNA changes

with_both <- nrow(target_RNA_tRNA_Control_integrated_table[target_RNA_tRNA_Control_integrated_table$tRNA_mut==1&target_RNA_tRNA_Control_integrated_table$is_RNA_exp_change==1&target_RNA_tRNA_Control_integrated_table$is_RNA_inherited==1&target_RNA_tRNA_Control_integrated_table$Time_match_to_tRNA==1, ])

tRNA_RNA <-  with_both

RNA_not_tRNA <- with_RNA_changes - tRNA_RNA

tRNA_not_RNA <- with_tRNA_changes  -  tRNA_RNA

not_tRNA_not_RNA <- background - (with_RNA_changes + with_tRNA_changes -  tRNA_RNA)


contingency_table <-
  rbind(c(tRNA_RNA, RNA_not_tRNA),
        c(tRNA_not_RNA, not_tRNA_not_RNA))
FT_out <- fisher.test(contingency_table)
pval_4 <- FT_out$p.value
OddsRatio_4 <- FT_out$estimate

# Table 2

# Create a table to show all the above results 1-4:

OR_frame <- as.data.frame(c(OddsRatio_1, OddsRatio_2, OddsRatio_3, OddsRatio_4))
Pval_frame <- as.data.frame(c(pval_1, pval_2, pval_3, pval_4))  

names <- c("Out of all genes, \nlikelihood of genes having both \nRNA epression changes and small RNA epimutations", 
           "Out of genes with RNAseq changes, \nlikelihood of genes having both \ninherited RNA expression changes and small RNA epimutations", 
           "Out of genes with RNAseq changes, \nlikelihood of genes having both \ninherited RNA expression changes and simultaneous small RNA epimutations",
           "Out of genes with RNAseq changes and small RNA changes, \nlikelihood of genes having both \n inherited RNA expression changes and simultaneous small RNA epimutations")


Table_2 <- cbind(names, OR_frame, Pval_frame)

colnames(Table_2) <- c("Association_tested", "Odds Ratio", "p-value")

#Odds calculation Low dose
# 1. What are the odds that genes have both small RNA epimutations and RNA changes ?

# background is all genes which have small RNA signal, so have to map to small RNA

Low_dose_background <- nrow(target_RNA_tRNAs_Low_dose_integrated_table)

# of which have RNA changes

Low_dose_with_RNA_changes <- nrow(target_RNA_tRNAs_Low_dose_integrated_table[target_RNA_tRNAs_Low_dose_integrated_table$is_RNA_exp_change==1, ])

# of which have tRNA changes 

Low_dose_with_tRNA_changes <- nrow(target_RNA_tRNAs_Low_dose_integrated_table[target_RNA_tRNAs_Low_dose_integrated_table$tRNA_mut==1, ])

# of which have both 

Low_dose_with_both <- nrow(target_RNA_tRNAs_Low_dose_integrated_table[target_RNA_tRNAs_Low_dose_integrated_table$tRNA_mut==1&target_RNA_tRNAs_Low_dose_integrated_table$is_RNA_exp_change==1, ])

Low_dose_tRNA_RNA <-  Low_dose_with_both

Low_dose_RNA_not_tRNA <- Low_dose_with_RNA_changes - Low_dose_tRNA_RNA

Low_dose_tRNA_not_RNA <- Low_dose_with_tRNA_changes  -  Low_dose_tRNA_RNA

Low_dose_not_tRNA_not_RNA <- Low_dose_background - (Low_dose_with_RNA_changes + Low_dose_with_tRNA_changes -  Low_dose_tRNA_RNA)

contingency_table <-
  rbind(c(Low_dose_tRNA_RNA, Low_dose_RNA_not_tRNA),
        c(Low_dose_tRNA_not_RNA, Low_dose_not_tRNA_not_RNA))
FT_out <- fisher.test(contingency_table)
pval_1 <- FT_out$p.value
OddsRatio_1 <- FT_out$estimate

# 2. What are the odds that genes with RNA expression changes have both INHERITED RNA expression changes and small RNA epimutations?

# background is all genes which map to small RNA and have AND RNA expression changes 

background <- nrow(target_RNA_tRNAs_Low_dose_integrated_table[target_RNA_tRNAs_Low_dose_integrated_table$is_RNA_exp_change==1, ])

# of which have inherited RNA changes

with_RNA_changes <- nrow(target_RNA_tRNAs_Low_dose_integrated_table[target_RNA_tRNAs_Low_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_Low_dose_integrated_table$is_RNA_inherited==1, ])

# of which have small RNA changes of any kind

with_tRNA_changes <- nrow(target_RNA_tRNAs_Low_dose_integrated_table[target_RNA_tRNAs_Low_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_Low_dose_integrated_table$tRNA_mut==1, ])

# of which have both inherited RNA changes and small RNA changes

with_both <- nrow(target_RNA_tRNAs_Low_dose_integrated_table[target_RNA_tRNAs_Low_dose_integrated_table$tRNA_mut==1&target_RNA_tRNAs_Low_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_Low_dose_integrated_table$is_RNA_inherited==1, ])

tRNA_RNA <-  with_both

RNA_not_tRNA <- with_RNA_changes - tRNA_RNA

tRNA_not_RNA <- with_tRNA_changes  -  tRNA_RNA

not_tRNA_not_RNA <- background - (with_RNA_changes + with_tRNA_changes -  tRNA_RNA)

contingency_table <-
  rbind(c(tRNA_RNA, RNA_not_tRNA),
        c(tRNA_not_RNA, not_tRNA_not_RNA))
FT_out <- fisher.test(contingency_table)
pval_2 <- FT_out$p.value
OddsRatio_2 <- FT_out$estimate

# Time matching

# 3. What are the odds that Genes with RNA expression changes have both inherited RNA expression changes and time matched small RNA epimutations?


# background is all genes which map to small RNA and have RNA expression changes

background <- nrow(target_RNA_tRNAs_Low_dose_integrated_table[target_RNA_tRNAs_Low_dose_integrated_table$is_RNA_exp_change==1, ])

# of which have inherited RNA changes

with_RNA_changes <- nrow(target_RNA_tRNAs_Low_dose_integrated_table[target_RNA_tRNAs_Low_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_Low_dose_integrated_table$is_RNA_inherited==1, ])

# of which have time matched small RNA changes 

with_tRNA_changes <- nrow(target_RNA_tRNAs_Low_dose_integrated_table[target_RNA_tRNAs_Low_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_Low_dose_integrated_table$tRNA_mut==1&target_RNA_tRNAs_Low_dose_integrated_table$Time_match_to_tRNA ==1, ])

# of which have both inherited RNA changes and time matched small RNA changes

with_both <- nrow(target_RNA_tRNAs_Low_dose_integrated_table[target_RNA_tRNAs_Low_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_Low_dose_integrated_table$is_RNA_inherited==1&target_RNA_tRNAs_Low_dose_integrated_table$tRNA_mut==1&target_RNA_tRNAs_Low_dose_integrated_table$Time_match_to_tRNA==1, ])

tRNA_RNA <-  with_both

RNA_not_tRNA <- with_RNA_changes - tRNA_RNA

tRNA_not_RNA <- with_tRNA_changes  -  tRNA_RNA

not_tRNA_not_RNA <- background - (with_RNA_changes + with_tRNA_changes -  tRNA_RNA)


contingency_table <-
  rbind(c(tRNA_RNA, RNA_not_tRNA),
        c(tRNA_not_RNA, not_tRNA_not_RNA))
FT_out <- fisher.test(contingency_table)
pval_3 <- FT_out$p.value
OddsRatio_3 <- FT_out$estimate

contingency_table <-
  rbind(c(73, 936),
        c(6, 18))
FT_out <- fisher.test(contingency_table)
pval_3 <- FT_out$p.value
OddsRatio_3 <- FT_out$estimate
# 4. What are the odds that in Genes with RNA expression changes and small RNA epimutations, the small RNA epimutations are time matched and the RNA expression changes are inherited?

# background is all genes which have small RNA data points AND RNA expression changes and small RNA epimutations

background <- nrow(target_RNA_tRNAs_Low_dose_integrated_table[target_RNA_tRNAs_Low_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_Low_dose_integrated_table$tRNA_mut==1, ])

# of which have inherited RNA changes

with_RNA_changes <- nrow(target_RNA_tRNAs_Low_dose_integrated_table[target_RNA_tRNAs_Low_dose_integrated_table$tRNA_mut==1&target_RNA_tRNAs_Low_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_Low_dose_integrated_table$is_RNA_inherited==1, ])

# of which have time matched small RNA changes

with_tRNA_changes <- nrow(target_RNA_tRNAs_Low_dose_integrated_table[target_RNA_tRNAs_Low_dose_integrated_table$tRNA_mut==1&target_RNA_tRNAs_Low_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_Low_dose_integrated_table$Time_match_to_tRNA==1, ])

# of which have both inherited RNA changes and time matched smallRNA changes

with_both <- nrow(target_RNA_tRNAs_Low_dose_integrated_table[target_RNA_tRNAs_Low_dose_integrated_table$tRNA_mut==1&target_RNA_tRNAs_Low_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_Low_dose_integrated_table$is_RNA_inherited==1&target_RNA_tRNAs_Low_dose_integrated_table$Time_match_to_tRNA==1, ])

tRNA_RNA <-  with_both

RNA_not_tRNA <- with_RNA_changes - tRNA_RNA

tRNA_not_RNA <- with_tRNA_changes  -  tRNA_RNA

not_tRNA_not_RNA <- background - (with_RNA_changes + with_tRNA_changes -  tRNA_RNA)

contingency_table <-
  rbind(c(tRNA_RNA, RNA_not_tRNA),
        c(tRNA_not_RNA, not_tRNA_not_RNA))
FT_out <- fisher.test(contingency_table)
pval_4 <- FT_out$p.value
OddsRatio_4 <- FT_out$estimate

# Table 2

# Create a table to show all the above results 1-4:

OR_frame <- as.data.frame(c(OddsRatio_1, OddsRatio_2, OddsRatio_3, OddsRatio_4))
Pval_frame <- as.data.frame(c(pval_1, pval_2, pval_3, pval_4))  

names <- c("Out of all genes, \nlikelihood of genes having both \nRNA epression changes and small RNA epimutations", 
           "Out of genes with RNAseq changes, \nlikelihood of genes having both \ninherited RNA expression changes and small RNA epimutations", 
           "Out of genes with RNAseq changes, \nlikelihood of genes having both \ninherited RNA expression changes and simultaneous small RNA epimutations",
           "Out of genes with RNAseq changes and small RNA changes, \nlikelihood of genes having both \n inherited RNA expression changes and simultaneous small RNA epimutations")


Table_2 <- cbind(names, OR_frame, Pval_frame)

colnames(Table_2) <- c("Association_tested", "Odds Ratio", "p-value")

#Odds calculation High dose
# 1. What are the odds that genes have both small RNA epimutations and RNA changes ?

# background is all genes which have small RNA signal, so have to map to small RNA

High_dose_background <- nrow(target_RNA_tRNAs_High_dose_integrated_table)

# of which have RNA changes

High_dose_with_RNA_changes <- nrow(target_RNA_tRNAs_High_dose_integrated_table[target_RNA_tRNAs_High_dose_integrated_table$is_RNA_exp_change==1, ])

# of which have small RNA changes 

High_dose_with_tRNA_changes <- nrow(target_RNA_tRNAs_High_dose_integrated_table[target_RNA_tRNAs_High_dose_integrated_table$tRNA_mut==1, ])

# of which have both 

High_dose_with_both <- nrow(target_RNA_tRNAs_High_dose_integrated_table[target_RNA_tRNAs_High_dose_integrated_table$tRNA_mut==1&target_RNA_tRNAs_High_dose_integrated_table$is_RNA_exp_change==1, ])

High_dose_tRNA_RNA <-  High_dose_with_both

High_dose_RNA_not_tRNA <- High_dose_with_RNA_changes - High_dose_tRNA_RNA

High_dose_tRNA_not_RNA <- High_dose_with_tRNA_changes  -  High_dose_tRNA_RNA

High_dose_not_tRNA_not_RNA <- High_dose_background - (High_dose_with_RNA_changes + High_dose_with_tRNA_changes -  tRNA_RNA)

contingency_table <-
  rbind(c(High_dose_tRNA_RNA, High_dose_RNA_not_tRNA),
        c(High_dose_tRNA_not_RNA, High_dose_not_tRNA_not_RNA))
FT_out <- fisher.test(contingency_table)
pval_1 <- FT_out$p.value
OddsRatio_1 <- FT_out$estimate

# 2. What are the odds that genes with RNA expression changes have both INHERITED RNA expression changes and small RNA epimutations?

# background is all genes which map to small RNA and have AND RNA expression changes 

background <- nrow(target_RNA_tRNAs_High_dose_integrated_table[target_RNA_tRNAs_High_dose_integrated_table$is_RNA_exp_change==1, ])

# of which have inherited RNA changes

with_RNA_changes <- nrow(target_RNA_tRNAs_High_dose_integrated_table[target_RNA_tRNAs_High_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_High_dose_integrated_table$is_RNA_inherited==1, ])

# of which have small RNA changes of any kind

with_tRNA_changes <- nrow(target_RNA_tRNAs_High_dose_integrated_table[target_RNA_tRNAs_High_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_High_dose_integrated_table$tRNA_mut==1, ])

# of which have both inherited RNA changes and small RNA changes

with_both <- nrow(target_RNA_tRNAs_High_dose_integrated_table[target_RNA_tRNAs_High_dose_integrated_table$tRNA_mut==1&target_RNA_tRNAs_High_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_High_dose_integrated_table$is_RNA_inherited==1, ])

tRNA_RNA <-  with_both

RNA_not_tRNA <- with_RNA_changes - tRNA_RNA

tRNA_not_RNA <- with_tRNA_changes  -  tRNA_RNA

not_tRNA_not_RNA <- background - (with_RNA_changes + with_tRNA_changes -  tRNA_RNA)

contingency_table <-
  rbind(c(tRNA_RNA, RNA_not_tRNA),
        c(tRNA_not_RNA, not_tRNA_not_RNA))
FT_out <- fisher.test(contingency_table)
pval_2 <- FT_out$p.value
OddsRatio_2 <- FT_out$estimate

# Time matching

# 3. What are the odds that Genes with RNA expression changes have both inherited RNA expression changes and time matched small RNA epimutations?


# background is all genes which map to small RNA and have RNA expression changes

background <- nrow(target_RNA_tRNAs_High_dose_integrated_table[target_RNA_tRNAs_High_dose_integrated_table$is_RNA_exp_change==1, ])

# of which have inherited RNA changes

with_RNA_changes <- nrow(target_RNA_tRNAs_High_dose_integrated_table[target_RNA_tRNAs_High_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_High_dose_integrated_table$is_RNA_inherited==1, ])

# of which have time matched small RNA changes 

with_tRNA_changes <- nrow(target_RNA_tRNAs_High_dose_integrated_table[target_RNA_tRNAs_High_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_High_dose_integrated_table$tRNA_mut==1&target_RNA_tRNAs_High_dose_integrated_table$Time_match_to_tRNA ==1, ])

# of which have both inherited RNA changes and time matched small RNA changes

with_both <- nrow(target_RNA_tRNAs_High_dose_integrated_table[target_RNA_tRNAs_High_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_High_dose_integrated_table$is_RNA_inherited==1&target_RNA_tRNAs_High_dose_integrated_table$tRNA_mut==1&target_RNA_tRNAs_High_dose_integrated_table$Time_match_to_tRNA==1, ])

tRNA_RNA <-  with_both

RNA_not_tRNA <- with_RNA_changes - tRNA_RNA

tRNA_not_RNA <- with_tRNA_changes  -  tRNA_RNA

not_tRNA_not_RNA <- background - (with_RNA_changes + with_tRNA_changes -  tRNA_RNA)

contingency_table <-
  rbind(c(tRNA_RNA, RNA_not_tRNA),
        c(tRNA_not_RNA, not_tRNA_not_RNA))
FT_out <- fisher.test(contingency_table)
pval_3 <- FT_out$p.value
OddsRatio_3 <- FT_out$estimate

# 4. What are the odds that in Genes with RNA expression changes and small RNA epimutations, the small RNA epimutations are time matched and the RNA expression changes are inherited?

# background is all genes which have small RNA data points AND RNA expression changes and small RNA epimutations

background <- nrow(target_RNA_tRNAs_High_dose_integrated_table[target_RNA_tRNAs_High_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_High_dose_integrated_table$tRNA_mut==1, ])

# of which have inherited RNA changes

with_RNA_changes <- nrow(target_RNA_tRNAs_High_dose_integrated_table[target_RNA_tRNAs_High_dose_integrated_table$tRNA_mut==1&target_RNA_tRNAs_High_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_High_dose_integrated_table$is_RNA_inherited==1, ])

# of which have time matched small RNA changes

with_tRNA_changes <- nrow(target_RNA_tRNAs_High_dose_integrated_table[target_RNA_tRNAs_High_dose_integrated_table$tRNA_mut==1&target_RNA_tRNAs_High_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_High_dose_integrated_table$Time_match_to_tRNA==1, ])

# of which have both inherited RNA changes and time matched smallRNA changes

with_both <- nrow(target_RNA_tRNAs_High_dose_integrated_table[target_RNA_tRNAs_High_dose_integrated_table$tRNA_mut==1&target_RNA_tRNAs_High_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_High_dose_integrated_table$is_RNA_inherited==1&target_RNA_tRNAs_High_dose_integrated_table$Time_match_to_tRNA==1, ])

tRNA_RNA <-  with_both

RNA_not_tRNA <- with_RNA_changes - tRNA_RNA

tRNA_not_RNA <- with_tRNA_changes  -  tRNA_RNA

not_tRNA_not_RNA <- background - (with_RNA_changes + with_tRNA_changes -  tRNA_RNA)

contingency_table <-
  rbind(c(tRNA_RNA, RNA_not_tRNA),
        c(tRNA_not_RNA, not_tRNA_not_RNA))
FT_out <- fisher.test(contingency_table)
pval_4 <- FT_out$p.value
OddsRatio_4 <- FT_out$estimate

# Create a table to show all the above results 1-4:

OR_frame <- as.data.frame(c(OddsRatio_1, OddsRatio_2, OddsRatio_3, OddsRatio_4))
Pval_frame <- as.data.frame(c(pval_1, pval_2, pval_3, pval_4))  

names <- c("Out of all genes, \nlikelihood of genes having both \nRNA epression changes and small RNA epimutations", 
           "Out of genes with RNAseq changes, \nlikelihood of genes having both \ninherited RNA expression changes and small RNA epimutations", 
           "Out of genes with RNAseq changes, \nlikelihood of genes having both \ninherited RNA expression changes and simultaneous small RNA epimutations",
           "Out of genes with RNAseq changes and small RNA changes, \nlikelihood of genes having both \n inherited RNA expression changes and simultaneous small RNA epimutations")

Table_2 <- cbind(names, OR_frame, Pval_frame)

colnames(Table_2) <- c("Association_tested", "Odds Ratio", "p-value")


#-------------
#######Fig.6.C - Concordance tRNAs fragments epimutation & gene exp.change####
#Calculation for Control
# RNA exp change  with tRNA concordancy

# Now look within the category of tRNA epimutated genes only, i.e. remove genes which may have inherited RNA changes but no tRNA epimutation
# Do we see any bias for directionality in terms of enrichment for inherited RNA changes?   

# Direction matched 
# background is all genes 
background <- nrow(target_RNA_tRNA_Control_integrated_table)

# of which is RNA expression changes only
RNA_only <- nrow(target_RNA_tRNA_Control_integrated_table[target_RNA_tRNA_Control_integrated_table$is_RNA_exp_change==1&target_RNA_tRNA_Control_integrated_table$tRNA_mut==0, ])

# of which is tRNA epimutations only
tRNA_only <-  nrow(target_RNA_tRNA_Control_integrated_table[target_RNA_tRNA_Control_integrated_table$is_RNA_exp_change==0&target_RNA_tRNA_Control_integrated_table$tRNA_mut==1,])

# of which is both time-matched concordant RNA and tRNA epimutations
concordant <- nrow(target_RNA_tRNA_Control_integrated_table[target_RNA_tRNA_Control_integrated_table$is_RNA_exp_change==1&target_RNA_tRNA_Control_integrated_table$tRNA_mut==1&target_RNA_tRNA_Control_integrated_table$Time_match_to_tRNA==1&target_RNA_tRNA_Control_integrated_table$Direction_match==1,])


contingency_table <- rbind(c(concordant, tRNA_only), 
                           c(RNA_only, background - (tRNA_only + concordant + RNA_only)))

FT_out <- fisher.test(contingency_table)
tRNA_concord_pval <- FT_out$p.value
tRNA_concord_OddsRatio <- FT_out$estimate

# Direction not matched 
# of which is both inherited RNA and discordant tRNA  epimutations
discordant <- nrow(target_RNA_tRNA_Control_integrated_table[target_RNA_tRNA_Control_integrated_table$is_RNA_exp_change==1&target_RNA_tRNA_Control_integrated_table$tRNA_mut==1&target_RNA_tRNA_Control_integrated_table$Time_match_to_tRNA==1&target_RNA_tRNA_Control_integrated_table$Direction_match==0,])

contingency_table <- rbind(c(discordant, tRNA_only), 
                           c(RNA_only, background - (tRNA_only + discordant + RNA_only)))

FT_out <- fisher.test(contingency_table)
tRNA_discord_pval <- FT_out$p.value
tRNA_discord_OddsRatio <- FT_out$estimate

bar_plot <-  as.data.frame(c(log2(tRNA_concord_OddsRatio), log2(tRNA_discord_OddsRatio)))

names <- c("All genes concordant epimutations", 
           "All genes discordant epimutations")

Pvalues <- c(tRNA_concord_pval, tRNA_discord_pval)
Histo_Bars_1 <- cbind(bar_plot, names, Pvalues)
colnames(Histo_Bars_1) <- c("Log_OR", "Names", "Pvalues")

# When we split into UP and DOWN and look within the category of epimutated genes only, i.e. remove genes which may have inherited RNA changes but no tRNA epimutation
# Do we see any bias for directionality in terms of enrichment for inherited RNA changes?   

# UP
# Direction matched 
both_UP <- nrow(target_RNA_tRNA_Control_integrated_table[target_RNA_tRNA_Control_integrated_table$is_RNA_exp_change==1&target_RNA_tRNA_Control_integrated_table$tRNA_mut==1&target_RNA_tRNA_Control_integrated_table$Time_match_to_tRNA==1&target_RNA_tRNA_Control_integrated_table$Direction_match==1&target_RNA_tRNA_Control_integrated_table$RNA_All_UP==1,])

contingency_table <- rbind(c(both_UP, tRNA_only), 
                           c(RNA_only, background - (tRNA_only + both_UP + RNA_only)))

FT_out <- fisher.test(contingency_table)
tRNA_UP_concord_pval <- FT_out$p.value
tRNA_UP_concord_OddsRatio <- FT_out$estimate

# Direction not matched 
UP_discordant <- nrow(target_RNA_tRNA_Control_integrated_table[target_RNA_tRNA_Control_integrated_table$is_RNA_exp_change==1&target_RNA_tRNA_Control_integrated_table$tRNA_mut==1&target_RNA_tRNA_Control_integrated_table$Time_match_to_tRNA==1&target_RNA_tRNA_Control_integrated_table$Direction_match==0&target_RNA_tRNA_Control_integrated_table$RNA_All_UP==1,])

contingency_table <- rbind(c(UP_discordant, tRNA_only), 
                           c(RNA_only, background - (tRNA_only + UP_discordant + RNA_only)))

FT_out <- fisher.test(contingency_table)
tRNA_UP_discord_pval <- FT_out$p.value
tRNA_UP_discord_OddsRatio <- FT_out$estimate

# DOWN
# Direction matched 

both_DOWN <- nrow(target_RNA_tRNA_Control_integrated_table[target_RNA_tRNA_Control_integrated_table$is_RNA_exp_change==1&target_RNA_tRNA_Control_integrated_table$tRNA_mut==1&target_RNA_tRNA_Control_integrated_table$Time_match_to_tRNA==1&target_RNA_tRNA_Control_integrated_table$Direction_match==1&target_RNA_tRNA_Control_integrated_table$RNA_All_DOWN==1,])

contingency_table <- rbind(c(both_DOWN, tRNA_only), 
                           c(RNA_only, background - (tRNA_only + both_DOWN + RNA_only)))

FT_out <- fisher.test(contingency_table)
tRNA_DOWN_concord_pval <- FT_out$p.value
tRNA_DOWN_concord_OddsRatio <- FT_out$estimate

# Direction not matched 

DOWN_discordant <- nrow(target_RNA_tRNA_Control_integrated_table[target_RNA_tRNA_Control_integrated_table$is_RNA_exp_change==1&target_RNA_tRNA_Control_integrated_table$tRNA_mut==1&target_RNA_tRNA_Control_integrated_table$Time_match_to_tRNA==1&target_RNA_tRNA_Control_integrated_table$Direction_match==0&target_RNA_tRNA_Control_integrated_table$RNA_All_DOWN==1,])

contingency_table <- rbind(c(DOWN_discordant, tRNA_only), 
                           c(RNA_only, background - (tRNA_only + DOWN_discordant + RNA_only)))

FT_out <- fisher.test(contingency_table)
tRNA_DOWN_discord_pval <- FT_out$p.value
tRNA_DOWN_discord_OddsRatio <- FT_out$estimate

bar_plot <-  as.data.frame(c(log2(tRNA_UP_concord_OddsRatio), log2(tRNA_DOWN_concord_OddsRatio), 
                             log2(tRNA_UP_discord_OddsRatio), log2(tRNA_DOWN_discord_OddsRatio)))

names <- c("Concordant UP",
           "Concordant DOWN",
           "Discordant UP", 
           "Discordant DOWN")

Pvalues <- c(tRNA_UP_concord_pval, tRNA_DOWN_concord_pval,
             tRNA_UP_discord_pval, tRNA_DOWN_discord_pval)

Histo_Bars_2 <- cbind(bar_plot, names, Pvalues)
colnames(Histo_Bars_2) <- c("Log_OR", "Names", "Pvalues")
E_Histo_Bars <- rbind(Histo_Bars_1, Histo_Bars_2)
E_Histo_Bars$Names <- factor(E_Histo_Bars$Names, levels = E_Histo_Bars$Names)
neg_log_P <- -log2(E_Histo_Bars$Pvalues)
E_Histo_Bars <- cbind(E_Histo_Bars, neg_log_P)

find_size <- c()
find_alpha <- c()

for(i in 1:nrow(E_Histo_Bars)){
  
  alpha <- 0.3
  size_point <- E_Histo_Bars[i,4]
  
  if(as.numeric(E_Histo_Bars[i, 3]) < 0.05){
    alpha <- 1
    size_point <- E_Histo_Bars[i,4]
    
  }
  find_alpha <- c(find_alpha, alpha)  
  find_size <- c(find_size, size_point)
}

E_Histo_Bars_tRNA_C <- cbind(E_Histo_Bars, find_alpha, find_size)
colnames(E_Histo_Bars_tRNA_C) <- c("log_OR", "Names", "Pvalues", "neg_log_P", "find_alpha", "find_size")
E_Histo_Bars_tRNA_C$Names <- factor(E_Histo_Bars_tRNA_C$Names, levels = rev(E_Histo_Bars_tRNA_C$Names))

guides_merge <- function(gdefs) {
  gdefs <- lapply(gdefs, function(g) { g$hash <- paste(g$order, g$hash, sep = "z"); g})
  tapply(gdefs, sapply(gdefs, function(g)g$hash), function(gs)Reduce(guide_merge, gs))
}
environment(guides_merge) <- environment(ggplot)
assignInNamespace("guides_merge", guides_merge, pos = "package:ggplot2")

RNA_tRNA_UP_DOWN <-
  
  ggplot(E_Histo_Bars_tRNA_C, aes(y=Names, x=as.numeric(log_OR)))+
  geom_point(aes(color=Names, size=find_size, alpha=find_alpha))+
  expand_limits(x=c(-3, 2))+
  geom_vline(xintercept = 0, colour = "grey")+
  theme_bw()+
  
  scale_size_continuous(range = c(1, 15), breaks = c(1, 3, 5), 
                        limits = c(1, 10))+
  
  scale_alpha(name = paste("Transparency indicates \nsignificance of enrichment"), range = c(0.3, 1),
              breaks = seq(1, 0.3, length = 2),
              limits = c(0.29, 1.05), 
              labels = c("p value significant < 0.05", "p value not significant > 0.05"))+
  
  ggtitle(paste(""))

RNA_tRNA_UP_DOWN <- RNA_tRNA_UP_DOWN + labs(y="Concordance subset\n", x=bquote(paste('log2(Odds Ratio) for enrichment for inherited RNAseq changes')))+
  theme(axis.title.x=element_text(face = "bold", size=14))+
  theme(axis.title.y=element_text(face = "bold", size=14))+
  theme(plot.title = element_text(face = "bold", size=16))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_discrete(labels= rev(E_Histo_Bars$Names))+
  theme(axis.text.x = element_text(color="#000000", size=12))+
  theme(axis.text.y = element_text(color="#000000", size = 16))+
  theme(legend.text=element_text(color="#000000", size=12))

RNA_tRNA_UP_DOWN 

#Calculation for Low dose
# RNA exp change  with tRNA concordancy
# Direction matched 
# background is all genes 
background <- nrow(target_RNA_tRNAs_Low_dose_integrated_table)

# of which is RNA expression changes only
RNA_only <- nrow(target_RNA_tRNAs_Low_dose_integrated_table[target_RNA_tRNAs_Low_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_Low_dose_integrated_table$tRNA_mut==0, ])

# of which is tRNA epimutations only
tRNA_only <-  nrow(target_RNA_tRNAs_Low_dose_integrated_table[target_RNA_tRNAs_Low_dose_integrated_table$is_RNA_exp_change==0&target_RNA_tRNAs_Low_dose_integrated_table$tRNA_mut==1,])

# of which is both time-matched concordant RNA and small RNA epimutations
concordant <- nrow(target_RNA_tRNAs_Low_dose_integrated_table[target_RNA_tRNAs_Low_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_Low_dose_integrated_table$tRNA_mut==1&target_RNA_tRNAs_Low_dose_integrated_table$Time_match_to_tRNA==1&target_RNA_tRNAs_Low_dose_integrated_table$Direction_match==1,])


contingency_table <- rbind(c(concordant, tRNA_only), 
                           c(RNA_only, background - (tRNA_only + concordant + RNA_only)))

FT_out <- fisher.test(contingency_table)
tRNA_concord_pval <- FT_out$p.value
tRNA_concord_OddsRatio <- FT_out$estimate

# Direction not matched 
# of which is both inherited RNA and discordant tRNA  epimutations
discordant <- nrow(target_RNA_tRNAs_Low_dose_integrated_table[target_RNA_tRNAs_Low_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_Low_dose_integrated_table$tRNA_mut==1&target_RNA_tRNAs_Low_dose_integrated_table$Time_match_to_tRNA==1&target_RNA_tRNAs_Low_dose_integrated_table$Direction_match==0,])

contingency_table <- rbind(c(discordant, tRNA_only), 
                           c(RNA_only, background - (tRNA_only + discordant + RNA_only)))

FT_out <- fisher.test(contingency_table)
tRNA_discord_pval <- FT_out$p.value
tRNA_discord_OddsRatio <- FT_out$estimate

bar_plot <-  as.data.frame(c(log2(tRNA_concord_OddsRatio), log2(tRNA_discord_OddsRatio)))

names <- c("All genes concordant epimutations", 
           "All genes discordant epimutations")

Pvalues <- c(tRNA_concord_pval, tRNA_discord_pval)
Histo_Bars_1 <- cbind(bar_plot, names, Pvalues)
colnames(Histo_Bars_1) <- c("Log_OR", "Names", "Pvalues")

# When we split into UP and DOWN and look within the category of epimutated genes only, i.e. remove genes which may have inherited RNA changes but no tRNA epimutation
# Do we see any bias for directionality in terms of enrichment for inherited RNA changes?   

# UP
# Direction matched 
both_UP <- nrow(target_RNA_tRNAs_Low_dose_integrated_table[target_RNA_tRNAs_Low_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_Low_dose_integrated_table$tRNA_mut==1&target_RNA_tRNAs_Low_dose_integrated_table$Time_match_to_tRNA==1&target_RNA_tRNAs_Low_dose_integrated_table$Direction_match==1&target_RNA_tRNAs_Low_dose_integrated_table$RNA_All_UP==1,])

contingency_table <- rbind(c(both_UP, tRNA_only), 
                           c(RNA_only, background - (tRNA_only + both_UP + RNA_only)))

FT_out <- fisher.test(contingency_table)
tRNA_UP_concord_pval <- FT_out$p.value
tRNA_UP_concord_OddsRatio <- FT_out$estimate

# Direction not matched 
UP_discordant <- nrow(target_RNA_tRNAs_Low_dose_integrated_table[target_RNA_tRNAs_Low_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_Low_dose_integrated_table$tRNA_mut==1&target_RNA_tRNAs_Low_dose_integrated_table$Time_match_to_tRNA==1&target_RNA_tRNAs_Low_dose_integrated_table$Direction_match==0&target_RNA_tRNAs_Low_dose_integrated_table$RNA_All_UP==1,])

contingency_table <- rbind(c(UP_discordant, tRNA_only), 
                           c(RNA_only, background - (tRNA_only + UP_discordant + RNA_only)))

FT_out <- fisher.test(contingency_table)
tRNA_UP_discord_pval <- FT_out$p.value
tRNA_UP_discord_OddsRatio <- FT_out$estimate

# DOWN
# Direction matched 

both_DOWN <- nrow(target_RNA_tRNAs_Low_dose_integrated_table[target_RNA_tRNAs_Low_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_Low_dose_integrated_table$tRNA_mut==1&target_RNA_tRNAs_Low_dose_integrated_table$Time_match_to_tRNA==1&target_RNA_tRNAs_Low_dose_integrated_table$Direction_match==1&target_RNA_tRNAs_Low_dose_integrated_table$RNA_All_DOWN==1,])

contingency_table <- rbind(c(both_DOWN, tRNA_only), 
                           c(RNA_only, background - (tRNA_only + both_DOWN + RNA_only)))

FT_out <- fisher.test(contingency_table)
tRNA_DOWN_concord_pval <- FT_out$p.value
tRNA_DOWN_concord_OddsRatio <- FT_out$estimate

# Direction not matched 

DOWN_discordant <- nrow(target_RNA_tRNAs_Low_dose_integrated_table[target_RNA_tRNAs_Low_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_Low_dose_integrated_table$tRNA_mut==1&target_RNA_tRNAs_Low_dose_integrated_table$Time_match_to_tRNA==1&target_RNA_tRNAs_Low_dose_integrated_table$Direction_match==0&target_RNA_tRNAs_Low_dose_integrated_table$RNA_All_DOWN==1,])

contingency_table <- rbind(c(DOWN_discordant, tRNA_only), 
                           c(RNA_only, background - (tRNA_only + DOWN_discordant + RNA_only)))

FT_out <- fisher.test(contingency_table)
tRNA_DOWN_discord_pval <- FT_out$p.value
tRNA_DOWN_discord_OddsRatio <- FT_out$estimate

bar_plot <-  as.data.frame(c(log2(tRNA_UP_concord_OddsRatio), log2(tRNA_DOWN_concord_OddsRatio), 
                             log2(tRNA_UP_discord_OddsRatio), log2(tRNA_DOWN_discord_OddsRatio)))

names <- c("Concordant UP",
           "Concordant DOWN",
           "Discordant UP", 
           "Discordant DOWN")

Pvalues <- c(tRNA_UP_concord_pval, tRNA_DOWN_concord_pval,
             tRNA_UP_discord_pval, tRNA_DOWN_discord_pval)

Histo_Bars_2 <- cbind(bar_plot, names, Pvalues)
colnames(Histo_Bars_2) <- c("Log_OR", "Names", "Pvalues")
E_Histo_Bars <- rbind(Histo_Bars_1, Histo_Bars_2)
E_Histo_Bars$Names <- factor(E_Histo_Bars$Names, levels = E_Histo_Bars$Names)
neg_log_P <- -log2(E_Histo_Bars$Pvalues)
E_Histo_Bars <- cbind(E_Histo_Bars, neg_log_P)

find_size <- c()
find_alpha <- c()

for(i in 1:nrow(E_Histo_Bars)){
  
  alpha <- 0.3
  size_point <- E_Histo_Bars[i,4]
  
  if(as.numeric(E_Histo_Bars[i, 3]) < 0.05){
    alpha <- 1
    size_point <- E_Histo_Bars[i,4]
    
  }
  find_alpha <- c(find_alpha, alpha)  
  find_size <- c(find_size, size_point)
}

E_Histo_Bars_tRNA_LD <- cbind(E_Histo_Bars, find_alpha, find_size)
colnames(E_Histo_Bars_tRNA_LD) <- c("log_OR", "Names", "Pvalues", "neg_log_P", "find_alpha", "find_size")
E_Histo_Bars_tRNA_LD$Names <- factor(E_Histo_Bars_tRNA_LD$Names, levels = rev(E_Histo_Bars_tRNA_LD$Names))

guides_merge <- function(gdefs) {
  gdefs <- lapply(gdefs, function(g) { g$hash <- paste(g$order, g$hash, sep = "z"); g})
  tapply(gdefs, sapply(gdefs, function(g)g$hash), function(gs)Reduce(guide_merge, gs))
}
environment(guides_merge) <- environment(ggplot)
assignInNamespace("guides_merge", guides_merge, pos = "package:ggplot2")

RNA_tRNA_UP_DOWN <-
  
  ggplot(E_Histo_Bars_tRNA_LD, aes(y=Names, x=as.numeric(log_OR)))+
  geom_point(aes(color=Names, size=find_size, alpha=find_alpha))+
  expand_limits(x=c(-3, 2))+
  geom_vline(xintercept = 0, colour = "grey")+
  theme_bw()+
  
  scale_size_continuous(range = c(1, 15), breaks = c(1, 3, 5), 
                        limits = c(1, 10))+
  
  scale_alpha(name = paste("Transparency indicates \nsignificance of enrichment"), range = c(0.3, 1),
              breaks = seq(1, 0.3, length = 2),
              limits = c(0.29, 1.05), 
              labels = c("p value significant < 0.05", "p value not significant > 0.05"))+
  
  ggtitle(paste(""))

RNA_tRNA_UP_DOWN <- RNA_tRNA_UP_DOWN + labs(y="Concordance subset\n", x=bquote(paste('log2(Odds Ratio) for enrichment for inherited RNAseq changes')))+
  theme(axis.title.x=element_text(face = "bold", size=14))+
  theme(axis.title.y=element_text(face = "bold", size=14))+
  theme(plot.title = element_text(face = "bold", size=16))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_discrete(labels= rev(E_Histo_Bars$Names))+
  theme(axis.text.x = element_text(color="#000000", size=12))+
  theme(axis.text.y = element_text(color="#000000", size = 16))+
  theme(legend.text=element_text(color="#000000", size=12))

RNA_tRNA_UP_DOWN 

#Calculation for High dose
# RNA exp change  with tRNA concordancy
# Direction matched 
# background is all genes 
background <- nrow(target_RNA_tRNAs_High_dose_integrated_table)

# of which is RNA expression changes only
RNA_only <- nrow(target_RNA_tRNAs_High_dose_integrated_table[target_RNA_tRNAs_High_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_High_dose_integrated_table$tRNA_mut==0, ])

# of which is tRNA epimutations only
tRNA_only <-  nrow(target_RNA_tRNAs_High_dose_integrated_table[target_RNA_tRNAs_High_dose_integrated_table$is_RNA_exp_change==0&target_RNA_tRNAs_High_dose_integrated_table$tRNA_mut==1,])

# of which is both time-matched concordant RNA and small RNA epimutations
concordant <- nrow(target_RNA_tRNAs_High_dose_integrated_table[target_RNA_tRNAs_High_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_High_dose_integrated_table$tRNA_mut==1&target_RNA_tRNAs_High_dose_integrated_table$Time_match_to_tRNA==1&target_RNA_tRNAs_High_dose_integrated_table$Direction_match==1,])


contingency_table <- rbind(c(concordant, tRNA_only), 
                           c(RNA_only, background - (tRNA_only + concordant + RNA_only)))

FT_out <- fisher.test(contingency_table)
tRNA_concord_pval <- FT_out$p.value
tRNA_concord_OddsRatio <- FT_out$estimate

# Direction not matched 
# of which is both inherited RNA and discordant tRNA  epimutations
discordant <- nrow(target_RNA_tRNAs_High_dose_integrated_table[target_RNA_tRNAs_High_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_High_dose_integrated_table$tRNA_mut==1&target_RNA_tRNAs_High_dose_integrated_table$Time_match_to_tRNA==1&target_RNA_tRNAs_High_dose_integrated_table$Direction_match==0,])

contingency_table <- rbind(c(discordant, tRNA_only), 
                           c(RNA_only, background - (tRNA_only + discordant + RNA_only)))

FT_out <- fisher.test(contingency_table)
tRNA_discord_pval <- FT_out$p.value
tRNA_discord_OddsRatio <- FT_out$estimate

bar_plot <-  as.data.frame(c(log2(tRNA_concord_OddsRatio), log2(tRNA_discord_OddsRatio)))

names <- c("All genes concordant epimutations", 
           "All genes discordant epimutations")

Pvalues <- c(tRNA_concord_pval, tRNA_discord_pval)
Histo_Bars_1 <- cbind(bar_plot, names, Pvalues)
colnames(Histo_Bars_1) <- c("Log_OR", "Names", "Pvalues")

# When we split into UP and DOWN and look within the category of epimutated genes only, i.e. remove genes which may have inherited RNA changes but no tRNA epimutation
# Do we see any bias for directionality in terms of enrichment for inherited RNA changes?   

# UP
# Direction matched 
both_UP <- nrow(target_RNA_tRNAs_High_dose_integrated_table[target_RNA_tRNAs_High_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_High_dose_integrated_table$tRNA_mut==1&target_RNA_tRNAs_High_dose_integrated_table$Time_match_to_tRNA==1&target_RNA_tRNAs_High_dose_integrated_table$Direction_match==1&target_RNA_tRNAs_High_dose_integrated_table$RNA_All_UP==1,])

contingency_table <- rbind(c(both_UP, tRNA_only), 
                           c(RNA_only, background - (tRNA_only + both_UP + RNA_only)))

FT_out <- fisher.test(contingency_table)
tRNA_UP_concord_pval <- FT_out$p.value
tRNA_UP_concord_OddsRatio <- FT_out$estimate

# Direction not matched 
UP_discordant <- nrow(target_RNA_tRNAs_High_dose_integrated_table[target_RNA_tRNAs_High_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_High_dose_integrated_table$tRNA_mut==1&target_RNA_tRNAs_High_dose_integrated_table$Time_match_to_tRNA==1&target_RNA_tRNAs_High_dose_integrated_table$Direction_match==0&target_RNA_tRNAs_High_dose_integrated_table$RNA_All_UP==1,])

contingency_table <- rbind(c(UP_discordant, tRNA_only), 
                           c(RNA_only, background - (tRNA_only + UP_discordant + RNA_only)))

FT_out <- fisher.test(contingency_table)
tRNA_UP_discord_pval <- FT_out$p.value
tRNA_UP_discord_OddsRatio <- FT_out$estimate

# DOWN
# Direction matched 

both_DOWN <- nrow(target_RNA_tRNAs_High_dose_integrated_table[target_RNA_tRNAs_High_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_High_dose_integrated_table$tRNA_mut==1&target_RNA_tRNAs_High_dose_integrated_table$Time_match_to_tRNA==1&target_RNA_tRNAs_High_dose_integrated_table$Direction_match==1&target_RNA_tRNAs_High_dose_integrated_table$RNA_All_DOWN==1,])

contingency_table <- rbind(c(both_DOWN, tRNA_only), 
                           c(RNA_only, background - (tRNA_only + both_DOWN + RNA_only)))

FT_out <- fisher.test(contingency_table)
tRNA_DOWN_concord_pval <- FT_out$p.value
tRNA_DOWN_concord_OddsRatio <- FT_out$estimate

# Direction not matched 

DOWN_discordant <- nrow(target_RNA_tRNAs_High_dose_integrated_table[target_RNA_tRNAs_High_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_High_dose_integrated_table$tRNA_mut==1&target_RNA_tRNAs_High_dose_integrated_table$Time_match_to_tRNA==1&target_RNA_tRNAs_High_dose_integrated_table$Direction_match==0&target_RNA_tRNAs_High_dose_integrated_table$RNA_All_DOWN==1,])

contingency_table <- rbind(c(DOWN_discordant, tRNA_only), 
                           c(RNA_only, background - (tRNA_only + DOWN_discordant + RNA_only)))

FT_out <- fisher.test(contingency_table)
tRNA_DOWN_discord_pval <- FT_out$p.value
tRNA_DOWN_discord_OddsRatio <- FT_out$estimate

bar_plot <-  as.data.frame(c(log2(tRNA_UP_concord_OddsRatio), log2(tRNA_DOWN_concord_OddsRatio), 
                             log2(tRNA_UP_discord_OddsRatio), log2(tRNA_DOWN_discord_OddsRatio)))

names <- c("Concordant UP",
           "Concordant DOWN",
           "Discordant UP", 
           "Discordant DOWN")

Pvalues <- c(tRNA_UP_concord_pval, tRNA_DOWN_concord_pval,
             tRNA_UP_discord_pval, tRNA_DOWN_discord_pval)

Histo_Bars_2 <- cbind(bar_plot, names, Pvalues)
colnames(Histo_Bars_2) <- c("Log_OR", "Names", "Pvalues")
E_Histo_Bars <- rbind(Histo_Bars_1, Histo_Bars_2)
E_Histo_Bars$Names <- factor(E_Histo_Bars$Names, levels = E_Histo_Bars$Names)
neg_log_P <- -log2(E_Histo_Bars$Pvalues)
E_Histo_Bars <- cbind(E_Histo_Bars, neg_log_P)

find_size <- c()
find_alpha <- c()

for(i in 1:nrow(E_Histo_Bars)){
  
  alpha <- 0.3
  size_point <- E_Histo_Bars[i,4]
  
  if(as.numeric(E_Histo_Bars[i, 3]) < 0.05){
    alpha <- 1
    size_point <- E_Histo_Bars[i,4]
    
  }
  find_alpha <- c(find_alpha, alpha)  
  find_size <- c(find_size, size_point)
}

E_Histo_Bars_tRNA_HD <- cbind(E_Histo_Bars, find_alpha, find_size)
colnames(E_Histo_Bars_tRNA_HD) <- c("log_OR", "Names", "Pvalues", "neg_log_P", "find_alpha", "find_size")
E_Histo_Bars_tRNA_HD$Names <- factor(E_Histo_Bars_tRNA_HD$Names, levels = rev(E_Histo_Bars_tRNA_HD$Names))

guides_merge <- function(gdefs) {
  gdefs <- lapply(gdefs, function(g) { g$hash <- paste(g$order, g$hash, sep = "z"); g})
  tapply(gdefs, sapply(gdefs, function(g)g$hash), function(gs)Reduce(guide_merge, gs))
}
environment(guides_merge) <- environment(ggplot)
assignInNamespace("guides_merge", guides_merge, pos = "package:ggplot2")

RNA_tRNA_UP_DOWN <-
  
  ggplot(E_Histo_Bars_tRNA_HD, aes(y=Names, x=as.numeric(log_OR)))+
  geom_point(aes(color=Names, size=find_size, alpha=find_alpha))+
  expand_limits(x=c(-3, 2))+
  geom_vline(xintercept = 0, colour = "grey")+
  theme_bw()+
  
  scale_size_continuous(range = c(1, 15), breaks = c(1, 3, 5), 
                        limits = c(1, 10))+
  
  scale_alpha(name = paste("Transparency indicates \nsignificance of enrichment"), range = c(0.3, 1),
              breaks = seq(1, 0.3, length = 2),
              limits = c(0.29, 1.05), 
              labels = c("p value significant < 0.05", "p value not significant > 0.05"))+
  
  ggtitle(paste(""))

RNA_tRNA_UP_DOWN <- RNA_tRNA_UP_DOWN + labs(y="Concordance subset\n", x=bquote(paste('log2(Odds Ratio) for enrichment for inherited RNAseq changes')))+
  theme(axis.title.x=element_text(face = "bold", size=14))+
  theme(axis.title.y=element_text(face = "bold", size=14))+
  theme(plot.title = element_text(face = "bold", size=16))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_discrete(labels= rev(E_Histo_Bars$Names))+
  theme(axis.text.x = element_text(color="#000000", size=12))+
  theme(axis.text.y = element_text(color="#000000", size = 16))+
  theme(legend.text=element_text(color="#000000", size=12))

RNA_tRNA_UP_DOWN 

#Bubble plot all conditions together
Histo_all<-rbind(E_Histo_Bars_tRNA_C,E_Histo_Bars_tRNA_LD,E_Histo_Bars_tRNA_HD)
conditions<-rep(c("Control","Low dose","High dose"),each=6)
Histo_all<-cbind(Histo_all,conditions)
Histo_all$conditions<-fct_relevel(Histo_all$conditions,c("Control","Low dose","High dose"))

Histo_all_plot <-
  
  ggplot(Histo_all, aes(y=Names, x=as.numeric(log_OR)))+
  geom_point(aes(color=conditions, size=(find_size+5), alpha=find_alpha))+
  expand_limits(x=c(-1, 4))+
  geom_vline(xintercept = 0, colour = "grey")+
  theme_bw()+
  
  scale_size_continuous(range = c(0, 15), breaks = c(1, 5,10,15), 
                        limits = c(0, 15))+
  
  scale_alpha(name = paste("Transparency indicates \nsignificance of enrichment"), range = c(0.3, 1),
              breaks = seq(1, 0.3, length = 2),
              limits = c(0.29, 1.05), 
              labels = c("p value significant < 0.05", "p value not significant > 0.05"))+
  
  ggtitle(paste(""))

#Fig.6.C

Histo_all_plot <- Histo_all_plot + labs(y="Concordance subset\n", x=bquote(paste('log2(Odds Ratio) for enrichment for RNAseq changes')))+
  theme(axis.title.x=element_text(face = "bold", size=14))+
  theme(axis.title.y=element_text(face = "bold", size=14))+
  theme(plot.title = element_text(face = "bold", size=16))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_discrete(labels= rev(E_Histo_Bars$Names))+
  scale_color_manual(values= c("cornflowerblue","darkgreen","red"))+
  theme(axis.text.x = element_text(color="#000000", size=12))+
  theme(axis.text.y = element_text(color="#000000", size = 16))+
  theme(legend.text=element_text(color="#000000", size=12))

Histo_all_plot 

#Get raw data
write.xlsx(Histo_all,"Data_Fig_6_C.xlsx")
#-------------
#######Table 4 - Comparison between tRNAs and 22G epimutations odds of association to gene exp. change####
##Control
#Nb genes with RNA epimutation and simultaneous 22G epimutation 
Control_genes_RNA_muts_and_smallRNA_muts <- nrow(Control_RNA_and_smallRNA_integrated_table[Control_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Control_RNA_and_smallRNA_integrated_table$RNA_mut==1&Control_RNA_and_smallRNA_integrated_table$smallRNA_mut==1&Control_RNA_and_smallRNA_integrated_table$time_matched_to_RNA==1, ])
#Nb genes with RNA epimutation 
Control_genes_RNA_muts_asso22G <- nrow(Control_RNA_and_smallRNA_integrated_table[Control_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Control_RNA_and_smallRNA_integrated_table$RNA_mut==1,])
#Nb genes with RNA epimutation and simultaneous tRNA epimutation in control condition
Control_genes_RNA_muts_and_tRNA_muts <- nrow(target_RNA_tRNA_Control_integrated_table[target_RNA_tRNA_Control_integrated_table$is_RNA_exp_change==1&target_RNA_tRNA_Control_integrated_table$tRNA_mut==1&target_RNA_tRNA_Control_integrated_table$Time_match_to_tRNA==1,])
#Nb genes with RNA epimutation 
Control_genes_RNA_muts_assotRNAs <- nrow(target_RNA_tRNA_Control_integrated_table[target_RNA_tRNA_Control_integrated_table$is_RNA_exp_change==1,])

contingency_table <-
  rbind(c(Control_genes_RNA_muts_and_smallRNA_muts, Control_genes_RNA_muts_and_tRNA_muts),
        c(Control_genes_RNA_muts_asso22G-Control_genes_RNA_muts_and_smallRNA_muts, Control_genes_RNA_muts_assotRNAs-Control_genes_RNA_muts_and_tRNA_muts))
FT_out <- fisher.test(contingency_table)
pval_1 <- FT_out$p.value
OddsRatio_1 <- FT_out$estimate

##Low dose
#Nb genes with RNA epimutation and simultaneous 22G epimutation 
Low_dose_genes_RNA_muts_and_smallRNA_muts <- nrow(Low_dose_RNA_and_smallRNA_integrated_table[Low_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Low_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1&Low_dose_RNA_and_smallRNA_integrated_table$smallRNA_mut==1&Low_dose_RNA_and_smallRNA_integrated_table$time_matched_to_RNA==1, ])
#Nb genes with RNA epimutation 
Low_dose_genes_RNA_muts_asso22G <- nrow(Low_dose_RNA_and_smallRNA_integrated_table[Low_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&Low_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1,])
#Nb genes with RNA epimutation and simultaneous tRNA epimutation in control condition
Low_dose_genes_RNA_muts_and_tRNA_muts <- nrow(target_RNA_tRNAs_Low_dose_integrated_table[target_RNA_tRNAs_Low_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_Low_dose_integrated_table$tRNA_mut==1&target_RNA_tRNAs_Low_dose_integrated_table$Time_match_to_tRNA==1,])
#Nb genes with RNA epimutation 
Low_dose_genes_RNA_muts_assotRNAs <- nrow(target_RNA_tRNAs_Low_dose_integrated_table[target_RNA_tRNAs_Low_dose_integrated_table$is_RNA_exp_change==1,])

contingency_table <-
  rbind(c(Low_dose_genes_RNA_muts_and_smallRNA_muts, Low_dose_genes_RNA_muts_and_tRNA_muts),
        c(Low_dose_genes_RNA_muts_asso22G-Low_dose_genes_RNA_muts_and_smallRNA_muts, Low_dose_genes_RNA_muts_assotRNAs-Low_dose_genes_RNA_muts_and_tRNA_muts))
FT_out <- fisher.test(contingency_table)
pval_1 <- FT_out$p.value
OddsRatio_1 <- FT_out$estimate

##High dose
#Nb genes with RNA epimutation and simultaneous 22G epimutation 
High_dose_genes_RNA_muts_and_smallRNA_muts <- nrow(High_dose_RNA_and_smallRNA_integrated_table[High_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&High_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1&High_dose_RNA_and_smallRNA_integrated_table$smallRNA_mut==1&High_dose_RNA_and_smallRNA_integrated_table$time_matched_to_RNA==1, ])
#Nb genes with RNA epimutation 
High_dose_genes_RNA_muts_asso22G <- nrow(High_dose_RNA_and_smallRNA_integrated_table[High_dose_RNA_and_smallRNA_integrated_table$maps_to_smallRNA==1&High_dose_RNA_and_smallRNA_integrated_table$RNA_mut==1,])
#Nb genes with RNA epimutation and simultaneous tRNA epimutation in control condition
High_dose_genes_RNA_muts_and_tRNA_muts <- nrow(target_RNA_tRNAs_High_dose_integrated_table[target_RNA_tRNAs_High_dose_integrated_table$is_RNA_exp_change==1&target_RNA_tRNAs_High_dose_integrated_table$tRNA_mut==1&target_RNA_tRNAs_High_dose_integrated_table$Time_match_to_tRNA==1,])
#Nb genes with RNA epimutation 
High_dose_genes_RNA_muts_assotRNAs <- nrow(target_RNA_tRNAs_High_dose_integrated_table[target_RNA_tRNAs_High_dose_integrated_table$is_RNA_exp_change==1,])

contingency_table <-
  rbind(c(High_dose_genes_RNA_muts_and_smallRNA_muts, High_dose_genes_RNA_muts_and_tRNA_muts),
        c(High_dose_genes_RNA_muts_asso22G-High_dose_genes_RNA_muts_and_smallRNA_muts, High_dose_genes_RNA_muts_assotRNAs-High_dose_genes_RNA_muts_and_tRNA_muts))
FT_out <- fisher.test(contingency_table)
pval_1 <- FT_out$p.value
OddsRatio_1 <- FT_out$estimate


#-------------
#######Fig.7.A & B####
#Construct Z-scores plot for each gene associated with a tRNA
#Data preparation
tRNA_z_scores_table_C1<-tRNA_z_scores_table_C1 %>% as.data.frame() 
colnames(tRNA_z_scores_table_C1) <- c('F2','4','6','8','10','12','14','16','18','20')
tRNA_z_scores_table_C1 <- tRNA_z_scores_table_C1 %>%
  # Creating an empty column:
  add_column(F0 = 0, .before='F2')
tRNA_z_scores_table_C1 <- tRNA_z_scores_table_C1 %>%
  # Creating a condition column:
  add_column(condition = "C1", .before="F0")
tRNA_z_scores_table_C1 <- cbind(rownames(tRNA_normAbun),tRNA_z_scores_table_C1)
colnames(tRNA_z_scores_table_C1) <- c('target','condition','0','2','4','6','8','10','12','14','16','18','20')
save(tRNA_z_scores_table_C1,file="tRNA_z_scores_table_C1.Rdata")

tRNA_z_scores_table_C2<-tRNA_z_scores_table_C2 %>% as.data.frame() 
colnames(tRNA_z_scores_table_C2) <- c('F2','4','8','10','12','14','16','18','20')
tRNA_z_scores_table_C2 <- tRNA_z_scores_table_C2 %>%
  # Creating an empty column:
  add_column(F0 = 0, .before='F2')
tRNA_z_scores_table_C2 <- tRNA_z_scores_table_C2 %>%
  # Creating a condition column:
  add_column(condition = "C2", .before="F0")
tRNA_z_scores_table_C2 <- cbind(rownames(tRNA_normAbun),tRNA_z_scores_table_C2)
colnames(tRNA_z_scores_table_C2) <- c('target','condition','0','2','4','8','10','12','14','16','18','20')
save(tRNA_z_scores_table_C2,file="tRNA_z_scores_table_C2.Rdata")

tRNA_z_scores_table_L1<-tRNA_z_scores_table_L1 %>% as.data.frame() 
colnames(tRNA_z_scores_table_L1) <- c('F2','4','6','8','12','14','16','18','20')
tRNA_z_scores_table_L1 <- tRNA_z_scores_table_L1 %>%
  # Creating an empty column:
  add_column(F0 = 0, .before='F2')
tRNA_z_scores_table_L1 <- tRNA_z_scores_table_L1 %>%
  # Creating a condition column:
  add_column(condition = "L1", .before="F0")
tRNA_z_scores_table_L1 <- cbind(rownames(tRNA_normAbun),tRNA_z_scores_table_L1)
colnames(tRNA_z_scores_table_L1) <- c('target', 'condition','0','2','4','6','8','12','14','16','18','20')
save(tRNA_z_scores_table_L1,file="tRNA_z_scores_table_L1.Rdata")

tRNA_z_scores_table_L2<-tRNA_z_scores_table_L2 %>% as.data.frame() 
colnames(tRNA_z_scores_table_L2) <- c('F2','4','8','10','12','14','16','18','20')
tRNA_z_scores_table_L2 <- tRNA_z_scores_table_L2 %>%
  # Creating an empty column:
  add_column(F0 = 0, .before='F2')
tRNA_z_scores_table_L2 <- tRNA_z_scores_table_L2 %>%
  # Creating a condition column:
  add_column(condition = "L2", .before="F0")
tRNA_z_scores_table_L2 <- cbind(rownames(tRNA_normAbun),tRNA_z_scores_table_L2)
colnames(tRNA_z_scores_table_L2) <- c('target','condition','0','2','4','8','10','12','14','16','18','20')
save(tRNA_z_scores_table_L2,file="tRNA_z_scores_table_L2.Rdata")

tRNA_z_scores_table_H1<-tRNA_z_scores_table_H1 %>% as.data.frame() 
colnames(tRNA_z_scores_table_H1) <- c('F2','4','8','10','12','14','16','18','20')
tRNA_z_scores_table_H1 <- tRNA_z_scores_table_H1 %>%
  # Creating an empty column:
  add_column(F0 = 0, .before='F2')
tRNA_z_scores_table_H1 <- tRNA_z_scores_table_H1 %>%
  # Creating a condition column:
  add_column(condition = "H1", .before="F0")
tRNA_z_scores_table_H1 <- cbind(rownames(tRNA_normAbun),tRNA_z_scores_table_H1)
colnames(tRNA_z_scores_table_H1) <- c('target','condition','0','2','4','8','10','12','14','16','18','20')
save(tRNA_z_scores_table_H1,file="tRNA_z_scores_table_H1.Rdata")

tRNA_z_scores_table_H2<-tRNA_z_scores_table_H2 %>% as.data.frame() 
colnames(tRNA_z_scores_table_H2) <- c('F2','4','6','8','10','12','14','16','18','20')
tRNA_z_scores_table_H2 <- tRNA_z_scores_table_H2 %>%
  # Creating an empty column:
  add_column(F0 = 0, .before='F2')
tRNA_z_scores_table_H2 <- tRNA_z_scores_table_H2 %>%
  # Creating a condition column:
  add_column(condition = "H2", .before="F0")
tRNA_z_scores_table_H2 <- cbind(rownames(tRNA_normAbun),tRNA_z_scores_table_H2)
colnames(tRNA_z_scores_table_H2) <- c('target','condition','0','2','4','6','8','10','12','14','16','18','20')
save(tRNA_z_scores_table_H2,file="tRNA_z_scores_table_H2.Rdata")

#Correction missing generations
tRNA_z_scores_table_C2bis<-cbind(tRNA_z_scores_table_C2, rep(tRNA_z_scores_table_C2[6],1))
tRNA_z_scores_table_C2bis<-tRNA_z_scores_table_C2bis[,c(1,2,3,4,5,13,6,7,8,9,10,11,12)]
colnames(tRNA_z_scores_table_C2bis) <- c("target","condition",'0','2','4','6','8','10','12','14','16','18','20')
write.xlsx(tRNA_z_scores_table_C1,"tRNA_z_scores_table_C1.xlsx")

tRNA_z_scores_table_L1bis<-cbind(tRNA_z_scores_table_L1, rep(tRNA_z_scores_table_L1[8],1))
tRNA_z_scores_table_L1bis<-tRNA_z_scores_table_L1bis[,c(1,2,3,4,5,6,7,13,8,9,10,11,12)]
colnames(tRNA_z_scores_table_L1bis) <- c("target","condition",'0','2','4','6','8','10','12','14','16','18','20')

tRNA_z_scores_table_L2bis<-cbind(tRNA_z_scores_table_L2, rep(tRNA_z_scores_table_L2[6],1))
tRNA_z_scores_table_L2bis<-tRNA_z_scores_table_L2bis[,c(1,2,3,4,5,13,6,7,8,9,10,11,12)]
colnames(tRNA_z_scores_table_L2bis) <- c("target","condition",'0','2','4','6','8','10','12','14','16','18','20')

tRNA_z_scores_table_H1bis<-cbind(tRNA_z_scores_table_H1, rep(tRNA_z_scores_table_H1[6],1))
tRNA_z_scores_table_H1bis<-tRNA_z_scores_table_H1bis[,c(1,2,3,4,5,13,6,7,8,9,10,11,12)]
colnames(tRNA_z_scores_table_H1bis) <- c("target","condition",'0','2','4','6','8','10','12','14','16','18','20')

#Association with WB names
setwd("~/Documents/Cisplatin project/Data_analysis/Analysis")
tRNA_target<-read.xlsx("tRNA_targetedGenes.xlsx")
tRNA_target_names<-tRNA_target$WB
tRNA_target$target<-paste(tRNA_target$AA,tRNA_target$codon,sep=".")

tRNA_z_scores_table_C1$target<-c("Trp.CCA","Gly.GCC","Pro.AGG","Leu.AAG","Tyr.GTA","Lys.CTT","Cys.GCA","Glu.CTC","Arg.TCG","Gly.TCC","Glu.TTC","Gly.GCC",
                                 "Asp.GTC","Gly.GCC","Cys.GCA","Ser.CGA","Val.AAC","Met.CAT","Ala.AGC","Gly.GCC","His.GTG","Val.TAC","Gly.GCC","Ala.AGC",
                                 "Gln.TTG","Pro.TGG","Leu.AAG","Leu.CAA","Pro.CGG","Gln.CTG","Tyr.GTA","Leu.AAG","Pro.AGG","Trp.CCA","Gly.GCC")
tRNA_z_scores_table_C1_WB <- merge(tRNA_z_scores_table_C1, tRNA_target, by = "target")
tRNA_z_scores_table_C1_WB<-cbind(tRNA_z_scores_table_C1_WB$WB,tRNA_z_scores_table_C1_WB[,1:13])
tRNA_z_scores_table_C1_WB<-tRNA_z_scores_table_C1_WB %>% distinct()

tRNA_z_scores_table_C2bis$target<-c("Trp.CCA","Gly.GCC","Pro.AGG","Leu.AAG","Tyr.GTA","Lys.CTT","Cys.GCA","Glu.CTC","Arg.TCG","Gly.TCC","Glu.TTC","Gly.GCC",
                                    "Asp.GTC","Gly.GCC","Cys.GCA","Ser.CGA","Val.AAC","Met.CAT","Ala.AGC","Gly.GCC","His.GTG","Val.TAC","Gly.GCC","Ala.AGC",
                                    "Gln.TTG","Pro.TGG","Leu.AAG","Leu.CAA","Pro.CGG","Gln.CTG","Tyr.GTA","Leu.AAG","Pro.AGG","Trp.CCA","Gly.GCC")
tRNA_z_scores_table_C2bis_WB <- merge(tRNA_z_scores_table_C2bis, tRNA_target, by = "target")
tRNA_z_scores_table_C2bis_WB<-cbind(tRNA_z_scores_table_C2bis_WB$WB,tRNA_z_scores_table_C2bis_WB[,1:13])
tRNA_z_scores_table_C2bis_WB<-tRNA_z_scores_table_C2bis_WB %>% distinct()

tRNA_z_scores_table_L1bis$target<-c("Trp.CCA","Gly.GCC","Pro.AGG","Leu.AAG","Tyr.GTA","Lys.CTT","Cys.GCA","Glu.CTC","Arg.TCG","Gly.TCC","Glu.TTC","Gly.GCC",
                                    "Asp.GTC","Gly.GCC","Cys.GCA","Ser.CGA","Val.AAC","Met.CAT","Ala.AGC","Gly.GCC","His.GTG","Val.TAC","Gly.GCC","Ala.AGC",
                                    "Gln.TTG","Pro.TGG","Leu.AAG","Leu.CAA","Pro.CGG","Gln.CTG","Tyr.GTA","Leu.AAG","Pro.AGG","Trp.CCA","Gly.GCC")
tRNA_z_scores_table_L1bis_WB <- merge(tRNA_z_scores_table_L1bis, tRNA_target, by = "target")
tRNA_z_scores_table_L1bis_WB<-cbind(tRNA_z_scores_table_L1bis_WB$WB,tRNA_z_scores_table_L1bis_WB[,1:13])
tRNA_z_scores_table_L1bis_WB<-tRNA_z_scores_table_L1bis_WB %>% distinct()

tRNA_z_scores_table_L2bis$target<-c("Trp.CCA","Gly.GCC","Pro.AGG","Leu.AAG","Tyr.GTA","Lys.CTT","Cys.GCA","Glu.CTC","Arg.TCG","Gly.TCC","Glu.TTC","Gly.GCC",
                                    "Asp.GTC","Gly.GCC","Cys.GCA","Ser.CGA","Val.AAC","Met.CAT","Ala.AGC","Gly.GCC","His.GTG","Val.TAC","Gly.GCC","Ala.AGC",
                                    "Gln.TTG","Pro.TGG","Leu.AAG","Leu.CAA","Pro.CGG","Gln.CTG","Tyr.GTA","Leu.AAG","Pro.AGG","Trp.CCA","Gly.GCC")
tRNA_z_scores_table_L2bis_WB <- merge(tRNA_z_scores_table_L2bis, tRNA_target, by = "target")
tRNA_z_scores_table_L2bis_WB<-cbind(tRNA_z_scores_table_L2bis_WB$WB,tRNA_z_scores_table_L2bis_WB[,1:13])
tRNA_z_scores_table_L2bis_WB<-tRNA_z_scores_table_L2bis_WB %>% distinct()

tRNA_z_scores_table_H1bis$target<-c("Trp.CCA","Gly.GCC","Pro.AGG","Leu.AAG","Tyr.GTA","Lys.CTT","Cys.GCA","Glu.CTC","Arg.TCG","Gly.TCC","Glu.TTC","Gly.GCC",
                                    "Asp.GTC","Gly.GCC","Cys.GCA","Ser.CGA","Val.AAC","Met.CAT","Ala.AGC","Gly.GCC","His.GTG","Val.TAC","Gly.GCC","Ala.AGC",
                                    "Gln.TTG","Pro.TGG","Leu.AAG","Leu.CAA","Pro.CGG","Gln.CTG","Tyr.GTA","Leu.AAG","Pro.AGG","Trp.CCA","Gly.GCC")
tRNA_z_scores_table_H1bis_WB <- merge(tRNA_z_scores_table_H1bis, tRNA_target, by = "target")
tRNA_z_scores_table_H1bis_WB<-cbind(tRNA_z_scores_table_H1bis_WB$WB,tRNA_z_scores_table_H1bis_WB[,1:13])
tRNA_z_scores_table_H1bis_WB<-tRNA_z_scores_table_H1bis_WB %>% distinct()

tRNA_z_scores_table_H2$target<-c("Trp.CCA","Gly.GCC","Pro.AGG","Leu.AAG","Tyr.GTA","Lys.CTT","Cys.GCA","Glu.CTC","Arg.TCG","Gly.TCC","Glu.TTC","Gly.GCC",
                                 "Asp.GTC","Gly.GCC","Cys.GCA","Ser.CGA","Val.AAC","Met.CAT","Ala.AGC","Gly.GCC","His.GTG","Val.TAC","Gly.GCC","Ala.AGC",
                                 "Gln.TTG","Pro.TGG","Leu.AAG","Leu.CAA","Pro.CGG","Gln.CTG","Tyr.GTA","Leu.AAG","Pro.AGG","Trp.CCA","Gly.GCC")
tRNA_z_scores_table_H2_WB <- merge(tRNA_z_scores_table_H2, tRNA_target, by = "target")
tRNA_z_scores_table_H2_WB<-cbind(tRNA_z_scores_table_H2_WB$WB,tRNA_z_scores_table_H2_WB[,1:13])
tRNA_z_scores_table_H2_WB<-tRNA_z_scores_table_H2_WB %>% distinct()

#Import RNA epimutations data
RNA_z_scores_table_C1<-RNA_z_scores_table_C1 %>% as.data.frame() 
colnames(RNA_z_scores_table_C1) <- c('F2','4','6','8','10','12','14','16','18','20')
RNA_z_scores_table_C1 <- RNA_z_scores_table_C1 %>%
  # Creating an empty column:
  add_column(F0 = 0, .before='F2')
RNA_z_scores_table_C1 <- RNA_z_scores_table_C1 %>%
  # Creating a condition column:
  add_column(condition = "C1", .before="F0")
RNA_z_scores_table_C1 <- RNA_z_scores_table_C1 %>%
  # Creating a genes column:
  add_column(genes = row.names(RNA_Lineage_C1), .before="condition")
colnames(RNA_z_scores_table_C1) <- c('genes','condition','0','2','4','6','8','10','12','14','16','18','20')
save(RNA_z_scores_table_C1,file="RNA_z_scores_table_C1.Rdata")

RNA_z_scores_table_C2<-RNA_z_scores_table_C2 %>% as.data.frame() 
colnames(RNA_z_scores_table_C2) <- c('F2','4','6','8','10','14','16','18','20')
RNA_z_scores_table_C2 <- RNA_z_scores_table_C2 %>%
  # Creating an empty column:
  add_column(F0 = 0, .before='F2')
RNA_z_scores_table_C2 <- RNA_z_scores_table_C2 %>%
  # Creating a condition column:
  add_column(condition = "C2", .before="F0")
RNA_z_scores_table_C2 <- RNA_z_scores_table_C2 %>%
  # Creating a genes column:
  add_column(genes = row.names(RNA_Lineage_C2), .before="condition")
RNA_z_scores_table_C2 <- RNA_z_scores_table_C2 %>%
  # Creating a NA column:
  add_column(F12 = "NA", .before="14")
colnames(RNA_z_scores_table_C2) <- c('genes','condition','0','2','4','6','8','10','12','14','16','18','20')
save(RNA_z_scores_table_C2,file="RNA_z_scores_table_C2.Rdata")

RNA_z_scores_table_L1 <- as.data.frame(RNA_z_scores_table_L1)
colnames(RNA_z_scores_table_L1) <- c('F2','4','6','8','10','12','14','16','18','20')
RNA_z_scores_table_L1 <- RNA_z_scores_table_L1 %>%
  # Creating an empty column:
  add_column(F0 = 0, .before="F2")
RNA_z_scores_table_L1 <- RNA_z_scores_table_L1 %>%
  # Creating an empty column:
  add_column(condition = 'L1', .before="F0")
RNA_z_scores_table_L1 <- RNA_z_scores_table_L1 %>%
  # Creating a genes column:
  add_column(genes = row.names(RNA_Lineage_L1), .before="condition")
colnames(RNA_z_scores_table_L1) <- c('genes','condition','0','2','4','6','8','10','12','14','16','18','20')
save(RNA_z_scores_table_L1,file="RNA_z_scores_table_L1.Rdata")

RNA_z_scores_table_L2 <- as.data.frame(RNA_z_scores_table_L2)
colnames(RNA_z_scores_table_L2) <- c('F2','4','6','8','10','12','14','16','18','20')
RNA_z_scores_table_L2 <- RNA_z_scores_table_L2 %>%
  # Creating an empty column:
  add_column(F0 = 0, .before="F2")
RNA_z_scores_table_L2 <- RNA_z_scores_table_L2 %>%
  # Creating an empty column:
  add_column(condition = 'L2', .before="F0")
RNA_z_scores_table_L2 <- RNA_z_scores_table_L2 %>%
  # Creating a genes column:
  add_column(genes = row.names(RNA_Lineage_L2), .before="condition")
colnames(RNA_z_scores_table_L2) <- c('genes','condition','0','2','4','6','8','10','12','14','16','18','20')
save(RNA_z_scores_table_L2,file="RNA_z_scores_table_L2.Rdata")

RNA_z_scores_table_H1<-RNA_z_scores_table_H1 %>% as.data.frame() 
colnames(RNA_z_scores_table_H1) <- c('F2','4','6','8','10','12','14','16','18')
RNA_z_scores_table_H1 <- RNA_z_scores_table_H1 %>%
  # Creating an empty column:
  add_column(F0 = 0, .before='F2')
RNA_z_scores_table_H1 <- RNA_z_scores_table_H1 %>%
  # Creating a condition column:
  add_column(condition = "H1", .before="F0")
RNA_z_scores_table_H1 <- RNA_z_scores_table_H1 %>%
  # Creating a genes column:
  add_column(genes = row.names(RNA_Lineage_H1), .before="condition")
RNA_z_scores_table_H1 <- RNA_z_scores_table_H1 %>%
  # Creating a NA column for the missing generation:
  add_column(F20 = "NA", .after="18")
RNA_z_scores_table_H1['F20'] <- RNA_z_scores_table_H1['18']
colnames(RNA_z_scores_table_H1) <- c('genes','condition','0','2','4','6','8','10','12','14','16','18','20')
save(RNA_z_scores_table_H1,file="RNA_z_scores_table_H1.Rdata")

RNA_z_scores_table_H2 <- as.data.frame(RNA_z_scores_table_H2)
colnames(RNA_z_scores_table_H2) <- c('F2','4','6','8','10','12','14','16','18','20')
RNA_z_scores_table_H2 <- RNA_z_scores_table_H2 %>%
  # Creating an empty column:
  add_column(F0 = 0, .before="F2")
RNA_z_scores_table_H2 <- RNA_z_scores_table_H2 %>%
  # Creating an empty column:
  add_column(condition = 'H2', .before="F0")
RNA_z_scores_table_H2 <- RNA_z_scores_table_H2 %>%
  # Creating a genes column:
  add_column(genes = row.names(RNA_Lineage_H2), .before="condition")
colnames(RNA_z_scores_table_H2) <- c('genes','condition','0','2','4','6','8','10','12','14','16','18','20')
save(RNA_z_scores_table_H2,file="RNA_z_scores_table_H2.Rdata")

#Import 22G epimutations data
z_scores_22G_C1<-tinyRNA_z_scores_table_C1 %>% as.data.frame() 
colnames(z_scores_22G_C1) <- c('F2','4','6','8','10','12','14','16','18','20')
names<- cbind(alltiny[,1:4],alltiny[,67])
z_scores_22G_C1<-cbind(names,z_scores_22G_C1)
z_scores_22G_C1<-subset(z_scores_22G_C1,z_scores_22G_C1$Tag=="22G")
z_scores_table_22G_C1<-cbind(z_scores_22G_C1$Feature.ID,z_scores_22G_C1[,6:15])
z_scores_table_22G_C1 <- z_scores_table_22G_C1 %>%
  # Creating an empty column:
  add_column(F0 = 0, .before='F2')
z_scores_table_22G_C1 <- z_scores_table_22G_C1 %>%
  # Creating a condition column:
  add_column(condition = "C1", .before="F0")
colnames(z_scores_table_22G_C1) <- c('genes','condition','0','2','4','6','8','10','12','14','16','18','20')
save(z_scores_table_22G_C1,file="22G-RNAs_z_scores_table_C1.Rdata")

z_scores_22G_H1<-tinyRNA_z_scores_table_H1 %>% as.data.frame() 
colnames(z_scores_22G_H1) <- c('F2','4','8','10','12','14','16','18','20')
z_scores_22G_H1bis<-cbind(z_scores_22G_H1, rep(z_scores_22G_H1[3],1))
z_scores_22G_H1bis<-z_scores_22G_H1bis[,c(1,2,10,3,4,5,6,7,8,9)]
names<- cbind(alltiny[,1:4],alltiny[,67])
z_scores_22G_H1bis<-cbind(names,z_scores_22G_H1bis)
z_scores_22G_H1bis<-subset(z_scores_22G_H1bis,z_scores_22G_H1bis$Tag=="22G")
z_scores_table_22G_H1<-cbind(z_scores_22G_H1bis$Feature.ID,z_scores_22G_H1bis[,6:15])
z_scores_table_22G_H1 <- z_scores_table_22G_H1 %>%
  # Creating an empty column:
  add_column(F0 = 0, .before='F2')
z_scores_table_22G_H1 <- z_scores_table_22G_H1 %>%
  # Creating a condition column:
  add_column(condition = "H1", .before="F0")
colnames(z_scores_table_22G_H1) <- c('genes','condition','0','2','4','6','8','10','12','14','16','18','20')
save(z_scores_table_22G_H1,file="22G-RNAs_z_scores_table_H1.Rdata")
#High dose
#gene 1
WBGene00001909_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[1,]
WBGene00001909_tRNA_H1<-WBGene00001909_tRNA_H1[,5:14]
WBGene00001909_tRNA_H1<-t(WBGene00001909_tRNA_H1)
generations<-rownames(WBGene00001909_tRNA_H1)
WBGene00001909_tRNA_H1<-cbind(generations,WBGene00001909_tRNA_H1)
colnames(WBGene00001909_tRNA_H1)<-c("Generation","Zscore")
WBGene00001909_tRNA_H1<-as.data.frame(WBGene00001909_tRNA_H1)
WBGene00001909_tRNA_H1$Generation<-as.numeric(WBGene00001909_tRNA_H1$Generation)
rownames(WBGene00001909_tRNA_H1)<-c()

WBGene00001909_tRNA_H2<-tRNA_z_scores_table_H2_WB[1,]
WBGene00001909_tRNA_H2<-WBGene00001909_tRNA_H2[,5:14]
WBGene00001909_tRNA_H2<-t(WBGene00001909_tRNA_H2)
generations<-rownames(WBGene00001909_tRNA_H2)
WBGene00001909_tRNA_H2<-cbind(generations,WBGene00001909_tRNA_H2)
colnames(WBGene00001909_tRNA_H2)<-c("Generation","Zscore")
WBGene00001909_tRNA_H2<-as.data.frame(WBGene00001909_tRNA_H2)
WBGene00001909_tRNA_H2$Generation<-as.numeric(WBGene00001909_tRNA_H2$Generation)
rownames(WBGene00001909_tRNA_H2)<-c()

RNA_z_scores_table_H1$genes<-sapply(strsplit(RNA_z_scores_table_H1$genes,":"), `[`, 5)
WBGene00001909_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00001909")
WBGene00001909_RNA_H1<-WBGene00001909_RNA_H1[,4:13]
WBGene00001909_RNA_H1<-t(WBGene00001909_RNA_H1)
generations<-rownames(WBGene00001909_RNA_H1)
WBGene00001909_RNA_H1<-cbind(generations,WBGene00001909_RNA_H1)
rownames(WBGene00001909_RNA_H1)<-c()
colnames(WBGene00001909_RNA_H1)<-c("Generation","Zscore")
WBGene00001909_RNA_H1<-as.data.frame(WBGene00001909_RNA_H1)

RNA_z_scores_table_H2$genes<-sapply(strsplit(RNA_z_scores_table_H2$genes,":"), `[`, 5)
WBGene00001909_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00001909")
WBGene00001909_RNA_H2<-WBGene00001909_RNA_H2[,4:13]
WBGene00001909_RNA_H2<-t(WBGene00001909_RNA_H2)
generations<-rownames(WBGene00001909_RNA_H2)
WBGene00001909_RNA_H2<-cbind(generations,WBGene00001909_RNA_H2)
rownames(WBGene00001909_RNA_H2)<-c()
colnames(WBGene00001909_RNA_H2)<-c("Generation","Zscore")
WBGene00001909_RNA_H2<-as.data.frame(WBGene00001909_RNA_H2)

plot(WBGene00001909_tRNA_H1$Generation,WBGene00001909_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="his-35 (WBGene00001909)/Ala.AGC ")
lines(WBGene00001909_tRNA_H2$Generation,WBGene00001909_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00001909_RNA_H1$Generation,WBGene00001909_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00001909_RNA_H2$Generation,WBGene00001909_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)

#gene 2
WBGene00001909_bis_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[2,]
WBGene00001909_bis_tRNA_H1<-WBGene00001909_bis_tRNA_H1[,5:14]
WBGene00001909_bis_tRNA_H1<-t(WBGene00001909_bis_tRNA_H1)
generations<-rownames(WBGene00001909_bis_tRNA_H1)
WBGene00001909_bis_tRNA_H1<-cbind(generations,WBGene00001909_bis_tRNA_H1)
colnames(WBGene00001909_bis_tRNA_H1)<-c("Generation","Zscore")
WBGene00001909_bis_tRNA_H1<-as.data.frame(WBGene00001909_bis_tRNA_H1)
WBGene00001909_bis_tRNA_H1$Generation<-as.numeric(WBGene00001909_bis_tRNA_H1$Generation)
rownames(WBGene00001909_bis_tRNA_H1)<-c()

WBGene00001909_bis_tRNA_H2<-tRNA_z_scores_table_H2_WB[2,]
WBGene00001909_bis_tRNA_H2<-WBGene00001909_bis_tRNA_H2[,5:14]
WBGene00001909_bis_tRNA_H2<-t(WBGene00001909_bis_tRNA_H2)
generations<-rownames(WBGene00001909_bis_tRNA_H2)
WBGene00001909_bis_tRNA_H2<-cbind(generations,WBGene00001909_bis_tRNA_H2)
colnames(WBGene00001909_bis_tRNA_H2)<-c("Generation","Zscore")
WBGene00001909_bis_tRNA_H2<-as.data.frame(WBGene00001909_bis_tRNA_H2)
WBGene00001909_bis_tRNA_H2$Generation<-as.numeric(WBGene00001909_bis_tRNA_H2$Generation)
rownames(WBGene00001909_bis_tRNA_H2)<-c()

plot(WBGene00001909_bis_tRNA_H1$Generation,WBGene00001909_bis_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="his-35 (WBGene00001909)/Ala.AGC.2")
lines(WBGene00001909_bis_tRNA_H2$Generation,WBGene00001909_bis_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00001909_RNA_H1$Generation,WBGene00001909_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00001909_RNA_H2$Generation,WBGene00001909_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)
#gene 3
WBGene00010843_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[3,]
WBGene00010843_tRNA_H1<-WBGene00010843_tRNA_H1[,5:14]
WBGene00010843_tRNA_H1<-t(WBGene00010843_tRNA_H1)
generations<-rownames(WBGene00010843_tRNA_H1)
WBGene00010843_tRNA_H1<-cbind(generations,WBGene00010843_tRNA_H1)
colnames(WBGene00010843_tRNA_H1)<-c("Generation","Zscore")
WBGene00010843_tRNA_H1<-as.data.frame(WBGene00010843_tRNA_H1)
WBGene00010843_tRNA_H1$Generation<-as.numeric(WBGene00010843_tRNA_H1$Generation)
rownames(WBGene00010843_tRNA_H1)<-c()

WBGene00010843_tRNA_H2<-tRNA_z_scores_table_H2_WB[3,]
WBGene00010843_tRNA_H2<-WBGene00010843_tRNA_H2[,5:14]
WBGene00010843_tRNA_H2<-t(WBGene00010843_tRNA_H2)
generations<-rownames(WBGene00010843_tRNA_H2)
WBGene00010843_tRNA_H2<-cbind(generations,WBGene00010843_tRNA_H2)
colnames(WBGene00010843_tRNA_H2)<-c("Generation","Zscore")
WBGene00010843_tRNA_H2<-as.data.frame(WBGene00010843_tRNA_H2)
WBGene00010843_tRNA_H2$Generation<-as.numeric(WBGene00010843_tRNA_H2$Generation)
rownames(WBGene00010843_tRNA_H2)<-c()

WBGene00010843_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00010843")
WBGene00010843_RNA_H1<-WBGene00010843_RNA_H1[,4:13]
WBGene00010843_RNA_H1<-t(WBGene00010843_RNA_H1)
generations<-rownames(WBGene00010843_RNA_H1)
WBGene00010843_RNA_H1<-cbind(generations,WBGene00010843_RNA_H1)
rownames(WBGene00010843_RNA_H1)<-c()
colnames(WBGene00010843_RNA_H1)<-c("Generation","Zscore")
WBGene00010843_RNA_H1<-as.data.frame(WBGene00010843_RNA_H1)

WBGene00010843_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00010843")
WBGene00010843_RNA_H2<-WBGene00010843_RNA_H2[,4:13]
WBGene00010843_RNA_H2<-t(WBGene00010843_RNA_H2)
generations<-rownames(WBGene00010843_RNA_H2)
WBGene00010843_RNA_H2<-cbind(generations,WBGene00010843_RNA_H2)
rownames(WBGene00010843_RNA_H2)<-c()
colnames(WBGene00010843_RNA_H2)<-c("Generation","Zscore")
WBGene00010843_RNA_H2<-as.data.frame(WBGene00010843_RNA_H2)
WBGene00010843_RNA_H2$Zscore[WBGene00010843_RNA_H2$Zscore == 'NA']<--1.097255

plot(WBGene00010843_tRNA_H1$Generation,WBGene00010843_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="M03H11.6 (WBGene00010843)/Arg.TCG")
lines(WBGene00010843_tRNA_H2$Generation,WBGene00010843_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00010843_RNA_H1$Generation,WBGene00010843_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00010843_RNA_H2$Generation,WBGene00010843_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)
#gene 4
WBGene00006792_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[4,]
WBGene00006792_tRNA_H1<-WBGene00006792_tRNA_H1[,5:14]
WBGene00006792_tRNA_H1<-t(WBGene00006792_tRNA_H1)
generations<-rownames(WBGene00006792_tRNA_H1)
WBGene00006792_tRNA_H1<-cbind(generations,WBGene00006792_tRNA_H1)
colnames(WBGene00006792_tRNA_H1)<-c("Generation","Zscore")
WBGene00006792_tRNA_H1<-as.data.frame(WBGene00006792_tRNA_H1)
WBGene00006792_tRNA_H1$Generation<-as.numeric(WBGene00006792_tRNA_H1$Generation)
rownames(WBGene00006792_tRNA_H1)<-c()

WBGene00006792_tRNA_H2<-tRNA_z_scores_table_H2_WB[4,]
WBGene00006792_tRNA_H2<-WBGene00006792_tRNA_H2[,5:14]
WBGene00006792_tRNA_H2<-t(WBGene00006792_tRNA_H2)
generations<-rownames(WBGene00006792_tRNA_H2)
WBGene00006792_tRNA_H2<-cbind(generations,WBGene00006792_tRNA_H2)
colnames(WBGene00006792_tRNA_H2)<-c("Generation","Zscore")
WBGene00006792_tRNA_H2<-as.data.frame(WBGene00006792_tRNA_H2)
WBGene00006792_tRNA_H2$Generation<-as.numeric(WBGene00006792_tRNA_H2$Generation)
rownames(WBGene00006792_tRNA_H2)<-c()

WBGene00006792_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00006792")
WBGene00006792_RNA_H1<-WBGene00006792_RNA_H1[,4:13]
WBGene00006792_RNA_H1<-t(WBGene00006792_RNA_H1)
generations<-rownames(WBGene00006792_RNA_H1)
WBGene00006792_RNA_H1<-cbind(generations,WBGene00006792_RNA_H1)
rownames(WBGene00006792_RNA_H1)<-c()
colnames(WBGene00006792_RNA_H1)<-c("Generation","Zscore")
WBGene00006792_RNA_H1<-as.data.frame(WBGene00006792_RNA_H1)

WBGene00006792_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00006792")
WBGene00006792_RNA_H2<-WBGene00006792_RNA_H2[,4:13]
WBGene00006792_RNA_H2<-t(WBGene00006792_RNA_H2)
generations<-rownames(WBGene00006792_RNA_H2)
WBGene00006792_RNA_H2<-cbind(generations,WBGene00006792_RNA_H2)
rownames(WBGene00006792_RNA_H2)<-c()
colnames(WBGene00006792_RNA_H2)<-c("Generation","Zscore")
WBGene00006792_RNA_H2<-as.data.frame(WBGene00006792_RNA_H2)
WBGene00006792_RNA_H2$Zscore[WBGene00006792_RNA_H2$Zscore == 'NA']<- -0.4514275

plot(WBGene00006792_tRNA_H1$Generation,WBGene00006792_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="unc-58 (WBGene00006792)/Arg.TCG")
lines(WBGene00006792_tRNA_H2$Generation,WBGene00006792_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00006792_RNA_H1$Generation,WBGene00006792_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00006792_RNA_H2$Generation,WBGene00006792_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)

#gene 5
WBGene00015483_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[5,]
WBGene00015483_tRNA_H1<-WBGene00015483_tRNA_H1[,5:14]
WBGene00015483_tRNA_H1<-t(WBGene00015483_tRNA_H1)
generations<-rownames(WBGene00015483_tRNA_H1)
WBGene00015483_tRNA_H1<-cbind(generations,WBGene00015483_tRNA_H1)
colnames(WBGene00015483_tRNA_H1)<-c("Generation","Zscore")
WBGene00015483_tRNA_H1<-as.data.frame(WBGene00015483_tRNA_H1)
WBGene00015483_tRNA_H1$Generation<-as.numeric(WBGene00015483_tRNA_H1$Generation)
rownames(WBGene00015483_tRNA_H1)<-c()

WBGene00015483_tRNA_H2<-tRNA_z_scores_table_H2_WB[5,]
WBGene00015483_tRNA_H2<-WBGene00015483_tRNA_H2[,5:14]
WBGene00015483_tRNA_H2<-t(WBGene00015483_tRNA_H2)
generations<-rownames(WBGene00015483_tRNA_H2)
WBGene00015483_tRNA_H2<-cbind(generations,WBGene00015483_tRNA_H2)
colnames(WBGene00015483_tRNA_H2)<-c("Generation","Zscore")
WBGene00015483_tRNA_H2<-as.data.frame(WBGene00015483_tRNA_H2)
WBGene00015483_tRNA_H2$Generation<-as.numeric(WBGene00015483_tRNA_H2$Generation)
rownames(WBGene00015483_tRNA_H2)<-c()

WBGene00015483_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00015483")
WBGene00015483_RNA_H1<-WBGene00015483_RNA_H1[,4:13]
WBGene00015483_RNA_H1<-t(WBGene00015483_RNA_H1)
generations<-rownames(WBGene00015483_RNA_H1)
WBGene00015483_RNA_H1<-cbind(generations,WBGene00015483_RNA_H1)
rownames(WBGene00015483_RNA_H1)<-c()
colnames(WBGene00015483_RNA_H1)<-c("Generation","Zscore")
WBGene00015483_RNA_H1<-as.data.frame(WBGene00015483_RNA_H1)

WBGene00015483_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00015483")
WBGene00015483_RNA_H2<-WBGene00015483_RNA_H2[,4:13]
WBGene00015483_RNA_H2<-t(WBGene00015483_RNA_H2)
generations<-rownames(WBGene00015483_RNA_H2)
WBGene00015483_RNA_H2<-cbind(generations,WBGene00015483_RNA_H2)
rownames(WBGene00015483_RNA_H2)<-c()
colnames(WBGene00015483_RNA_H2)<-c("Generation","Zscore")
WBGene00015483_RNA_H2<-as.data.frame(WBGene00015483_RNA_H2)
WBGene00015483_RNA_H2$Zscore[WBGene00015483_RNA_H2$Zscore == 'NA']<- -0.151748

plot(WBGene00015483_tRNA_H1$Generation,WBGene00015483_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="C05D11.5 (WBGene00015483)/Asp.GTC")
lines(WBGene00015483_tRNA_H2$Generation,WBGene00015483_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00015483_RNA_H1$Generation,WBGene00015483_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00015483_RNA_H2$Generation,WBGene00015483_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)
#gene 6
WBGene00018672_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[6,]
WBGene00018672_tRNA_H1<-WBGene00018672_tRNA_H1[,5:14]
WBGene00018672_tRNA_H1<-t(WBGene00018672_tRNA_H1)
generations<-rownames(WBGene00018672_tRNA_H1)
WBGene00018672_tRNA_H1<-cbind(generations,WBGene00018672_tRNA_H1)
colnames(WBGene00018672_tRNA_H1)<-c("Generation","Zscore")
WBGene00018672_tRNA_H1<-as.data.frame(WBGene00018672_tRNA_H1)
WBGene00018672_tRNA_H1$Generation<-as.numeric(WBGene00018672_tRNA_H1$Generation)
rownames(WBGene00018672_tRNA_H1)<-c()

WBGene00018672_tRNA_H2<-tRNA_z_scores_table_H2_WB[6,]
WBGene00018672_tRNA_H2<-WBGene00018672_tRNA_H2[,5:14]
WBGene00018672_tRNA_H2<-t(WBGene00018672_tRNA_H2)
generations<-rownames(WBGene00018672_tRNA_H2)
WBGene00018672_tRNA_H2<-cbind(generations,WBGene00018672_tRNA_H2)
colnames(WBGene00018672_tRNA_H2)<-c("Generation","Zscore")
WBGene00018672_tRNA_H2<-as.data.frame(WBGene00018672_tRNA_H2)
WBGene00018672_tRNA_H2$Generation<-as.numeric(WBGene00018672_tRNA_H2$Generation)
rownames(WBGene00018672_tRNA_H2)<-c()

WBGene00018672_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00018672")
WBGene00018672_RNA_H1<-WBGene00018672_RNA_H1[,4:13]
WBGene00018672_RNA_H1<-t(WBGene00018672_RNA_H1)
generations<-rownames(WBGene00018672_RNA_H1)
WBGene00018672_RNA_H1<-cbind(generations,WBGene00018672_RNA_H1)
rownames(WBGene00018672_RNA_H1)<-c()
colnames(WBGene00018672_RNA_H1)<-c("Generation","Zscore")
WBGene00018672_RNA_H1<-as.data.frame(WBGene00018672_RNA_H1)

WBGene00018672_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00018672")
WBGene00018672_RNA_H2<-WBGene00018672_RNA_H2[,4:13]
WBGene00018672_RNA_H2<-t(WBGene00018672_RNA_H2)
generations<-rownames(WBGene00018672_RNA_H2)
WBGene00018672_RNA_H2<-cbind(generations,WBGene00018672_RNA_H2)
rownames(WBGene00018672_RNA_H2)<-c()
colnames(WBGene00018672_RNA_H2)<-c("Generation","Zscore")
WBGene00018672_RNA_H2<-as.data.frame(WBGene00018672_RNA_H2)
WBGene00018672_RNA_H2$Zscore[WBGene00018672_RNA_H2$Zscore == 'NA']<- -1.000521

plot(WBGene00018672_tRNA_H1$Generation,WBGene00018672_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="sorf-2 (WBGene00018672)/Cys.GCA")
lines(WBGene00018672_tRNA_H2$Generation,WBGene00018672_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00018672_RNA_H1$Generation,WBGene00018672_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00018672_RNA_H2$Generation,WBGene00018672_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)
#gene 7
WBGene00000218_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[7,]
WBGene00000218_tRNA_H1<-WBGene00000218_tRNA_H1[,5:14]
WBGene00000218_tRNA_H1<-t(WBGene00000218_tRNA_H1)
generations<-rownames(WBGene00000218_tRNA_H1)
WBGene00000218_tRNA_H1<-cbind(generations,WBGene00000218_tRNA_H1)
colnames(WBGene00000218_tRNA_H1)<-c("Generation","Zscore")
WBGene00000218_tRNA_H1<-as.data.frame(WBGene00000218_tRNA_H1)
WBGene00000218_tRNA_H1$Generation<-as.numeric(WBGene00000218_tRNA_H1$Generation)
rownames(WBGene00000218_tRNA_H1)<-c()

WBGene00000218_tRNA_H2<-tRNA_z_scores_table_H2_WB[7,]
WBGene00000218_tRNA_H2<-WBGene00000218_tRNA_H2[,5:14]
WBGene00000218_tRNA_H2<-t(WBGene00000218_tRNA_H2)
generations<-rownames(WBGene00000218_tRNA_H2)
WBGene00000218_tRNA_H2<-cbind(generations,WBGene00000218_tRNA_H2)
colnames(WBGene00000218_tRNA_H2)<-c("Generation","Zscore")
WBGene00000218_tRNA_H2<-as.data.frame(WBGene00000218_tRNA_H2)
WBGene00000218_tRNA_H2$Generation<-as.numeric(WBGene00000218_tRNA_H2$Generation)
rownames(WBGene00000218_tRNA_H2)<-c()

WBGene00000218_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00000218")
WBGene00000218_RNA_H1<-WBGene00000218_RNA_H1[,4:13]
WBGene00000218_RNA_H1<-t(WBGene00000218_RNA_H1)
generations<-rownames(WBGene00000218_RNA_H1)
WBGene00000218_RNA_H1<-cbind(generations,WBGene00000218_RNA_H1)
rownames(WBGene00000218_RNA_H1)<-c()
colnames(WBGene00000218_RNA_H1)<-c("Generation","Zscore")
WBGene00000218_RNA_H1<-as.data.frame(WBGene00000218_RNA_H1)

WBGene00000218_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00000218")
WBGene00000218_RNA_H2<-WBGene00000218_RNA_H2[,4:13]
WBGene00000218_RNA_H2<-t(WBGene00000218_RNA_H2)
generations<-rownames(WBGene00000218_RNA_H2)
WBGene00000218_RNA_H2<-cbind(generations,WBGene00000218_RNA_H2)
rownames(WBGene00000218_RNA_H2)<-c()
colnames(WBGene00000218_RNA_H2)<-c("Generation","Zscore")
WBGene00000218_RNA_H2<-as.data.frame(WBGene00000218_RNA_H2)
WBGene00000218_RNA_H2$Zscore[WBGene00000218_RNA_H2$Zscore == 'NA']<- -1.153808

plot(WBGene00000218_tRNA_H1$Generation,WBGene00000218_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="asp-5 (WBGene00000218)/Cys.GCA")
lines(WBGene00000218_tRNA_H2$Generation,WBGene00000218_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00000218_RNA_H1$Generation,WBGene00000218_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00000218_RNA_H2$Generation,WBGene00000218_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)
#gene 8
WBGene00018672_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[8,]
WBGene00018672_tRNA_H1<-WBGene00018672_tRNA_H1[,5:14]
WBGene00018672_tRNA_H1<-t(WBGene00018672_tRNA_H1)
generations<-rownames(WBGene00018672_tRNA_H1)
WBGene00018672_tRNA_H1<-cbind(generations,WBGene00018672_tRNA_H1)
colnames(WBGene00018672_tRNA_H1)<-c("Generation","Zscore")
WBGene00018672_tRNA_H1<-as.data.frame(WBGene00018672_tRNA_H1)
WBGene00018672_tRNA_H1$Generation<-as.numeric(WBGene00018672_tRNA_H1$Generation)
rownames(WBGene00018672_tRNA_H1)<-c()

WBGene00018672_tRNA_H2<-tRNA_z_scores_table_H2_WB[8,]
WBGene00018672_tRNA_H2<-WBGene00018672_tRNA_H2[,5:14]
WBGene00018672_tRNA_H2<-t(WBGene00018672_tRNA_H2)
generations<-rownames(WBGene00018672_tRNA_H2)
WBGene00018672_tRNA_H2<-cbind(generations,WBGene00018672_tRNA_H2)
colnames(WBGene00018672_tRNA_H2)<-c("Generation","Zscore")
WBGene00018672_tRNA_H2<-as.data.frame(WBGene00018672_tRNA_H2)
WBGene00018672_tRNA_H2$Generation<-as.numeric(WBGene00018672_tRNA_H2$Generation)
rownames(WBGene00018672_tRNA_H2)<-c()

WBGene00018672_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00018672")
WBGene00018672_RNA_H1<-WBGene00018672_RNA_H1[,4:13]
WBGene00018672_RNA_H1<-t(WBGene00018672_RNA_H1)
generations<-rownames(WBGene00018672_RNA_H1)
WBGene00018672_RNA_H1<-cbind(generations,WBGene00018672_RNA_H1)
rownames(WBGene00018672_RNA_H1)<-c()
colnames(WBGene00018672_RNA_H1)<-c("Generation","Zscore")
WBGene00018672_RNA_H1<-as.data.frame(WBGene00018672_RNA_H1)

WBGene00018672_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00018672")
WBGene00018672_RNA_H2<-WBGene00018672_RNA_H2[,4:13]
WBGene00018672_RNA_H2<-t(WBGene00018672_RNA_H2)
generations<-rownames(WBGene00018672_RNA_H2)
WBGene00018672_RNA_H2<-cbind(generations,WBGene00018672_RNA_H2)
rownames(WBGene00018672_RNA_H2)<-c()
colnames(WBGene00018672_RNA_H2)<-c("Generation","Zscore")
WBGene00018672_RNA_H2<-as.data.frame(WBGene00018672_RNA_H2)
WBGene00018672_RNA_H2$Zscore[WBGene00018672_RNA_H2$Zscore == 'NA']<- 1.000521

plot(WBGene00018672_tRNA_H1$Generation,WBGene00018672_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="sorf-2 (WBGene00018672)/Cys.GCA")
lines(WBGene00018672_tRNA_H2$Generation,WBGene00018672_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00018672_RNA_H1$Generation,WBGene00018672_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00018672_RNA_H2$Generation,WBGene00018672_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)
#gene 9
WBGene00000218_bis_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[9,]
WBGene00000218_bis_tRNA_H1<-WBGene00000218_bis_tRNA_H1[,5:14]
WBGene00000218_bis_tRNA_H1<-t(WBGene00000218_bis_tRNA_H1)
generations<-rownames(WBGene00000218_bis_tRNA_H1)
WBGene00000218_bis_tRNA_H1<-cbind(generations,WBGene00000218_bis_tRNA_H1)
colnames(WBGene00000218_bis_tRNA_H1)<-c("Generation","Zscore")
WBGene00000218_bis_tRNA_H1<-as.data.frame(WBGene00000218_bis_tRNA_H1)
WBGene00000218_bis_tRNA_H1$Generation<-as.numeric(WBGene00000218_bis_tRNA_H1$Generation)
rownames(WBGene00000218_bis_tRNA_H1)<-c()

WBGene00000218_tRNA_H2<-tRNA_z_scores_table_H2_WB[9,]
WBGene00000218_tRNA_H2<-WBGene00000218_tRNA_H2[,5:14]
WBGene00000218_tRNA_H2<-t(WBGene00000218_tRNA_H2)
generations<-rownames(WBGene00000218_tRNA_H2)
WBGene00000218_tRNA_H2<-cbind(generations,WBGene00000218_tRNA_H2)
colnames(WBGene00000218_tRNA_H2)<-c("Generation","Zscore")
WBGene00000218_tRNA_H2<-as.data.frame(WBGene00000218_tRNA_H2)
WBGene00000218_tRNA_H2$Generation<-as.numeric(WBGene00000218_tRNA_H2$Generation)
rownames(WBGene00000218_tRNA_H2)<-c()

WBGene00000218_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00000218")
WBGene00000218_RNA_H1<-WBGene00000218_RNA_H1[,4:13]
WBGene00000218_RNA_H1<-t(WBGene00000218_RNA_H1)
generations<-rownames(WBGene00000218_RNA_H1)
WBGene00000218_RNA_H1<-cbind(generations,WBGene00000218_RNA_H1)
rownames(WBGene00000218_RNA_H1)<-c()
colnames(WBGene00000218_RNA_H1)<-c("Generation","Zscore")
WBGene00000218_RNA_H1<-as.data.frame(WBGene00000218_RNA_H1)

WBGene00000218_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00000218")
WBGene00000218_RNA_H2<-WBGene00000218_RNA_H2[,4:13]
WBGene00000218_RNA_H2<-t(WBGene00000218_RNA_H2)
generations<-rownames(WBGene00000218_RNA_H2)
WBGene00000218_RNA_H2<-cbind(generations,WBGene00000218_RNA_H2)
rownames(WBGene00000218_RNA_H2)<-c()
colnames(WBGene00000218_RNA_H2)<-c("Generation","Zscore")
WBGene00000218_RNA_H2<-as.data.frame(WBGene00000218_RNA_H2)
WBGene00000218_RNA_H2$Zscore[WBGene00000218_RNA_H2$Zscore == 'NA']<- -1.153808

plot(WBGene00000218_tRNA_H1$Generation,WBGene00000218_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="asp-5 (WBGene00000218)/Cys.GCA.2")
lines(WBGene00000218_tRNA_H2$Generation,WBGene00000218_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00000218_RNA_H1$Generation,WBGene00000218_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00000218_RNA_H2$Generation,WBGene00000218_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)
#gene 10
WBGene00019472_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[10,]
WBGene00019472_tRNA_H1<-WBGene00019472_tRNA_H1[,5:14]
WBGene00019472_tRNA_H1<-t(WBGene00019472_tRNA_H1)
generations<-rownames(WBGene00019472_tRNA_H1)
WBGene00019472_tRNA_H1<-cbind(generations,WBGene00019472_tRNA_H1)
colnames(WBGene00019472_tRNA_H1)<-c("Generation","Zscore")
WBGene00019472_tRNA_H1<-as.data.frame(WBGene00019472_tRNA_H1)
WBGene00019472_tRNA_H1$Generation<-as.numeric(WBGene00019472_tRNA_H1$Generation)
rownames(WBGene00019472_tRNA_H1)<-c()

WBGene00019472_tRNA_H2<-tRNA_z_scores_table_H2_WB[10,]
WBGene00019472_tRNA_H2<-WBGene00019472_tRNA_H2[,5:14]
WBGene00019472_tRNA_H2<-t(WBGene00019472_tRNA_H2)
generations<-rownames(WBGene00019472_tRNA_H2)
WBGene00019472_tRNA_H2<-cbind(generations,WBGene00019472_tRNA_H2)
colnames(WBGene00019472_tRNA_H2)<-c("Generation","Zscore")
WBGene00019472_tRNA_H2<-as.data.frame(WBGene00019472_tRNA_H2)
WBGene00019472_tRNA_H2$Generation<-as.numeric(WBGene00019472_tRNA_H2$Generation)
rownames(WBGene00019472_tRNA_H2)<-c()

WBGene00019472_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00019472")
WBGene00019472_RNA_H1<-WBGene00019472_RNA_H1[,4:13]
WBGene00019472_RNA_H1<-t(WBGene00019472_RNA_H1)
generations<-rownames(WBGene00019472_RNA_H1)
WBGene00019472_RNA_H1<-cbind(generations,WBGene00019472_RNA_H1)
rownames(WBGene00019472_RNA_H1)<-c()
colnames(WBGene00019472_RNA_H1)<-c("Generation","Zscore")
WBGene00019472_RNA_H1<-as.data.frame(WBGene00019472_RNA_H1)

WBGene00019472_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00019472")
WBGene00019472_RNA_H2<-WBGene00019472_RNA_H2[,4:13]
WBGene00019472_RNA_H2<-t(WBGene00019472_RNA_H2)
generations<-rownames(WBGene00019472_RNA_H2)
WBGene00019472_RNA_H2<-cbind(generations,WBGene00019472_RNA_H2)
rownames(WBGene00019472_RNA_H2)<-c()
colnames(WBGene00019472_RNA_H2)<-c("Generation","Zscore")
WBGene00019472_RNA_H2<-as.data.frame(WBGene00019472_RNA_H2)
WBGene00019472_RNA_H2$Zscore[WBGene00019472_RNA_H2$Zscore == 'NA']<- -0.531126

#no associated 22G

plot(WBGene00019472_tRNA_H1$Generation,WBGene00019472_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Generation",ylab="Zscore",main="cyp-35B1 (WBGene00019472)")
lines(WBGene00019472_RNA_H1$Generation,WBGene00019472_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00019472_tRNA_H1_bis$Generation,WBGene00019472_tRNA_H1_bis$Zscore, col="cadetblue1",type="l")
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(16, -3, legend=c("RNA", "Gln.CTG_tRNA","Gln.TTG_tRNA"),
       col=c("red", "blue","cadetblue1"), lty=c(1,1,1), cex=0.8)

#Get raw data
write.xlsx(WBGene00019472_tRNA_H1,"Data_Gln.CTG_tRNA.xlsx")
write.xlsx(WBGene00019472_RNA_H1,"RNA.xlsx")
write.xlsx(WBGene00019472_tRNA_H1_bis,"Data_Gln.TTG_tRNA.xlsx")

#gene 11
WBGene00009926_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[11,]
WBGene00009926_tRNA_H1<-WBGene00009926_tRNA_H1[,5:14]
WBGene00009926_tRNA_H1<-t(WBGene00009926_tRNA_H1)
generations<-rownames(WBGene00009926_tRNA_H1)
WBGene00009926_tRNA_H1<-cbind(generations,WBGene00009926_tRNA_H1)
colnames(WBGene00009926_tRNA_H1)<-c("Generation","Zscore")
WBGene00009926_tRNA_H1<-as.data.frame(WBGene00009926_tRNA_H1)
WBGene00009926_tRNA_H1$Generation<-as.numeric(WBGene00009926_tRNA_H1$Generation)
rownames(WBGene00009926_tRNA_H1)<-c()

WBGene00009926_tRNA_H2<-tRNA_z_scores_table_H2_WB[11,]
WBGene00009926_tRNA_H2<-WBGene00009926_tRNA_H2[,5:14]
WBGene00009926_tRNA_H2<-t(WBGene00009926_tRNA_H2)
generations<-rownames(WBGene00009926_tRNA_H2)
WBGene00009926_tRNA_H2<-cbind(generations,WBGene00009926_tRNA_H2)
colnames(WBGene00009926_tRNA_H2)<-c("Generation","Zscore")
WBGene00009926_tRNA_H2<-as.data.frame(WBGene00009926_tRNA_H2)
WBGene00009926_tRNA_H2$Generation<-as.numeric(WBGene00009926_tRNA_H2$Generation)
rownames(WBGene00009926_tRNA_H2)<-c()

WBGene00009926_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00009926")
WBGene00009926_RNA_H1<-WBGene00009926_RNA_H1[,4:13]
WBGene00009926_RNA_H1<-t(WBGene00009926_RNA_H1)
generations<-rownames(WBGene00009926_RNA_H1)
WBGene00009926_RNA_H1<-cbind(generations,WBGene00009926_RNA_H1)
rownames(WBGene00009926_RNA_H1)<-c()
colnames(WBGene00009926_RNA_H1)<-c("Generation","Zscore")
WBGene00009926_RNA_H1<-as.data.frame(WBGene00009926_RNA_H1)

WBGene00009926_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00009926")
WBGene00009926_RNA_H2<-WBGene00009926_RNA_H2[,4:13]
WBGene00009926_RNA_H2<-t(WBGene00009926_RNA_H2)
generations<-rownames(WBGene00009926_RNA_H2)
WBGene00009926_RNA_H2<-cbind(generations,WBGene00009926_RNA_H2)
rownames(WBGene00009926_RNA_H2)<-c()
colnames(WBGene00009926_RNA_H2)<-c("Generation","Zscore")
WBGene00009926_RNA_H2<-as.data.frame(WBGene00009926_RNA_H2)
WBGene00009926_RNA_H2$Zscore[WBGene00009926_RNA_H2$Zscore == 'NA']<- 0.5026323

plot(WBGene00009926_tRNA_H1$Generation,WBGene00009926_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="noah-2 (WBGene00009926)/Gln.CTG")
lines(WBGene00009926_tRNA_H2$Generation,WBGene00009926_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00009926_RNA_H1$Generation,WBGene00009926_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00009926_RNA_H2$Generation,WBGene00009926_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)
#gene 12
WBGene00019472_tRNA_H1_bis<-tRNA_z_scores_table_H1bis_WB[12,]
WBGene00019472_tRNA_H1_bis<-WBGene00019472_tRNA_H1_bis[,5:14]
WBGene00019472_tRNA_H1_bis<-t(WBGene00019472_tRNA_H1_bis)
generations<-rownames(WBGene00019472_tRNA_H1_bis)
WBGene00019472_tRNA_H1_bis<-cbind(generations,WBGene00019472_tRNA_H1_bis)
colnames(WBGene00019472_tRNA_H1_bis)<-c("Generation","Zscore")
WBGene00019472_tRNA_H1_bis<-as.data.frame(WBGene00019472_tRNA_H1_bis)
WBGene00019472_tRNA_H1_bis$Generation<-as.numeric(WBGene00019472_tRNA_H1_bis$Generation)
rownames(WBGene00019472_tRNA_H1_bis)<-c()

WBGene00019472_tRNA_H2_bis<-tRNA_z_scores_table_H2_WB[12,]
WBGene00019472_tRNA_H2_bis<-WBGene00019472_tRNA_H2_bis[,5:14]
WBGene00019472_tRNA_H2_bis<-t(WBGene00019472_tRNA_H2_bis)
generations<-rownames(WBGene00019472_tRNA_H2_bis)
WBGene00019472_tRNA_H2_bis<-cbind(generations,WBGene00019472_tRNA_H2_bis)
colnames(WBGene00019472_tRNA_H2_bis)<-c("Generation","Zscore")
WBGene00019472_tRNA_H2_bis<-as.data.frame(WBGene00019472_tRNA_H2_bis)
WBGene00019472_tRNA_H2_bis$Generation<-as.numeric(WBGene00019472_tRNA_H2_bis$Generation)
rownames(WBGene00019472_tRNA_H2_bis)<-c()

WBGene00019472_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00019472")
WBGene00019472_RNA_H1<-WBGene00019472_RNA_H1[,4:13]
WBGene00019472_RNA_H1<-t(WBGene00019472_RNA_H1)
generations<-rownames(WBGene00019472_RNA_H1)
WBGene00019472_RNA_H1<-cbind(generations,WBGene00019472_RNA_H1)
rownames(WBGene00019472_RNA_H1)<-c()
colnames(WBGene00019472_RNA_H1)<-c("Generation","Zscore")
WBGene00019472_RNA_H1<-as.data.frame(WBGene00019472_RNA_H1)

WBGene00019472_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00019472")
WBGene00019472_RNA_H2<-WBGene00019472_RNA_H2[,4:13]
WBGene00019472_RNA_H2<-t(WBGene00019472_RNA_H2)
generations<-rownames(WBGene00019472_RNA_H2)
WBGene00019472_RNA_H2<-cbind(generations,WBGene00019472_RNA_H2)
rownames(WBGene00019472_RNA_H2)<-c()
colnames(WBGene00019472_RNA_H2)<-c("Generation","Zscore")
WBGene00019472_RNA_H2<-as.data.frame(WBGene00019472_RNA_H2)
WBGene00019472_RNA_H2$Zscore[WBGene00019472_RNA_H2$Zscore == 'NA']<- -0.531126

plot(WBGene00019472_tRNA_H1$Generation,WBGene00019472_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="cyp-35B1 (WBGene00019472)/Gln.TTG")
lines(WBGene00019472_tRNA_H2$Generation,WBGene00019472_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00019472_RNA_H1$Generation,WBGene00019472_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00019472_RNA_H2$Generation,WBGene00019472_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)
#gene 13
WBGene00007660_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[13,]
WBGene00007660_tRNA_H1<-WBGene00007660_tRNA_H1[,5:14]
WBGene00007660_tRNA_H1<-t(WBGene00007660_tRNA_H1)
generations<-rownames(WBGene00007660_tRNA_H1)
WBGene00007660_tRNA_H1<-cbind(generations,WBGene00007660_tRNA_H1)
colnames(WBGene00007660_tRNA_H1)<-c("Generation","Zscore")
WBGene00007660_tRNA_H1<-as.data.frame(WBGene00007660_tRNA_H1)
WBGene00007660_tRNA_H1$Generation<-as.numeric(WBGene00007660_tRNA_H1$Generation)
rownames(WBGene00007660_tRNA_H1)<-c()

WBGene00007660_tRNA_H2<-tRNA_z_scores_table_H2_WB[13,]
WBGene00007660_tRNA_H2<-WBGene00007660_tRNA_H2[,5:14]
WBGene00007660_tRNA_H2<-t(WBGene00007660_tRNA_H2)
generations<-rownames(WBGene00007660_tRNA_H2)
WBGene00007660_tRNA_H2<-cbind(generations,WBGene00007660_tRNA_H2)
colnames(WBGene00007660_tRNA_H2)<-c("Generation","Zscore")
WBGene00007660_tRNA_H2<-as.data.frame(WBGene00007660_tRNA_H2)
WBGene00007660_tRNA_H2$Generation<-as.numeric(WBGene00007660_tRNA_H2$Generation)
rownames(WBGene00007660_tRNA_H2)<-c()

WBGene00007660_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00007660")
WBGene00007660_RNA_H1<-WBGene00007660_RNA_H1[,4:13]
WBGene00007660_RNA_H1<-t(WBGene00007660_RNA_H1)
generations<-rownames(WBGene00007660_RNA_H1)
WBGene00007660_RNA_H1<-cbind(generations,WBGene00007660_RNA_H1)
rownames(WBGene00007660_RNA_H1)<-c()
colnames(WBGene00007660_RNA_H1)<-c("Generation","Zscore")
WBGene00007660_RNA_H1<-as.data.frame(WBGene00007660_RNA_H1)

WBGene00007660_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00007660")
WBGene00007660_RNA_H2<-WBGene00007660_RNA_H2[,4:13]
WBGene00007660_RNA_H2<-t(WBGene00007660_RNA_H2)
generations<-rownames(WBGene00007660_RNA_H2)
WBGene00007660_RNA_H2<-cbind(generations,WBGene00007660_RNA_H2)
rownames(WBGene00007660_RNA_H2)<-c()
colnames(WBGene00007660_RNA_H2)<-c("Generation","Zscore")
WBGene00007660_RNA_H2<-as.data.frame(WBGene00007660_RNA_H2)
WBGene00007660_RNA_H2$Zscore[WBGene00007660_RNA_H2$Zscore == 'NA']<- -2.699975

plot(WBGene00007660_tRNA_H1$Generation,WBGene00007660_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="pals-6 (WBGene00007660)/Gln.TTG")
lines(WBGene00007660_tRNA_H2$Generation,WBGene00007660_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00007660_RNA_H1$Generation,WBGene00007660_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00007660_RNA_H2$Generation,WBGene00007660_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)

#gene 14
WBGene00020699_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[14,]
WBGene00020699_tRNA_H1<-WBGene00020699_tRNA_H1[,5:14]
WBGene00020699_tRNA_H1<-t(WBGene00020699_tRNA_H1)
generations<-rownames(WBGene00020699_tRNA_H1)
WBGene00020699_tRNA_H1<-cbind(generations,WBGene00020699_tRNA_H1)
colnames(WBGene00020699_tRNA_H1)<-c("Generation","Zscore")
WBGene00020699_tRNA_H1<-as.data.frame(WBGene00020699_tRNA_H1)
WBGene00020699_tRNA_H1$Generation<-as.numeric(WBGene00020699_tRNA_H1$Generation)
rownames(WBGene00020699_tRNA_H1)<-c()

WBGene00020699_tRNA_H2<-tRNA_z_scores_table_H2_WB[14,]
WBGene00020699_tRNA_H2<-WBGene00020699_tRNA_H2[,5:14]
WBGene00020699_tRNA_H2<-t(WBGene00020699_tRNA_H2)
generations<-rownames(WBGene00020699_tRNA_H2)
WBGene00020699_tRNA_H2<-cbind(generations,WBGene00020699_tRNA_H2)
colnames(WBGene00020699_tRNA_H2)<-c("Generation","Zscore")
WBGene00020699_tRNA_H2<-as.data.frame(WBGene00020699_tRNA_H2)
WBGene00020699_tRNA_H2$Generation<-as.numeric(WBGene00020699_tRNA_H2$Generation)
rownames(WBGene00020699_tRNA_H2)<-c()

WBGene00020699_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00020699")
WBGene00020699_RNA_H1<-WBGene00020699_RNA_H1[,4:13]
WBGene00020699_RNA_H1<-t(WBGene00020699_RNA_H1)
generations<-rownames(WBGene00020699_RNA_H1)
WBGene00020699_RNA_H1<-cbind(generations,WBGene00020699_RNA_H1)
rownames(WBGene00020699_RNA_H1)<-c()
colnames(WBGene00020699_RNA_H1)<-c("Generation","Zscore")
WBGene00020699_RNA_H1<-as.data.frame(WBGene00020699_RNA_H1)

WBGene00020699_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00020699")
WBGene00020699_RNA_H2<-WBGene00020699_RNA_H2[,4:13]
WBGene00020699_RNA_H2<-t(WBGene00020699_RNA_H2)
generations<-rownames(WBGene00020699_RNA_H2)
WBGene00020699_RNA_H2<-cbind(generations,WBGene00020699_RNA_H2)
rownames(WBGene00020699_RNA_H2)<-c()
colnames(WBGene00020699_RNA_H2)<-c("Generation","Zscore")
WBGene00020699_RNA_H2<-as.data.frame(WBGene00020699_RNA_H2)
WBGene00020699_RNA_H2$Zscore[WBGene00020699_RNA_H2$Zscore == 'NA']<- 0.2111528

plot(WBGene00020699_tRNA_H1$Generation,WBGene00020699_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="T22F3.10 (WBGene00020699)/Glu.CTC")
lines(WBGene00020699_tRNA_H2$Generation,WBGene00020699_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00020699_RNA_H1$Generation,WBGene00020699_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00020699_RNA_H2$Generation,WBGene00020699_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)

#gene 15
WBGene00006759_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[15,]
WBGene00006759_tRNA_H1<-WBGene00006759_tRNA_H1[,5:14]
WBGene00006759_tRNA_H1<-t(WBGene00006759_tRNA_H1)
generations<-rownames(WBGene00006759_tRNA_H1)
WBGene00006759_tRNA_H1<-cbind(generations,WBGene00006759_tRNA_H1)
colnames(WBGene00006759_tRNA_H1)<-c("Generation","Zscore")
WBGene00006759_tRNA_H1<-as.data.frame(WBGene00006759_tRNA_H1)
WBGene00006759_tRNA_H1$Generation<-as.numeric(WBGene00006759_tRNA_H1$Generation)
rownames(WBGene00006759_tRNA_H1)<-c()

WBGene00006759_tRNA_H2<-tRNA_z_scores_table_H2_WB[15,]
WBGene00006759_tRNA_H2<-WBGene00006759_tRNA_H2[,5:14]
WBGene00006759_tRNA_H2<-t(WBGene00006759_tRNA_H2)
generations<-rownames(WBGene00006759_tRNA_H2)
WBGene00006759_tRNA_H2<-cbind(generations,WBGene00006759_tRNA_H2)
colnames(WBGene00006759_tRNA_H2)<-c("Generation","Zscore")
WBGene00006759_tRNA_H2<-as.data.frame(WBGene00006759_tRNA_H2)
WBGene00006759_tRNA_H2$Generation<-as.numeric(WBGene00006759_tRNA_H2$Generation)
rownames(WBGene00006759_tRNA_H2)<-c()

WBGene00006759_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00006759")
WBGene00006759_RNA_H1<-WBGene00006759_RNA_H1[,4:13]
WBGene00006759_RNA_H1<-t(WBGene00006759_RNA_H1)
generations<-rownames(WBGene00006759_RNA_H1)
WBGene00006759_RNA_H1<-cbind(generations,WBGene00006759_RNA_H1)
rownames(WBGene00006759_RNA_H1)<-c()
colnames(WBGene00006759_RNA_H1)<-c("Generation","Zscore")
WBGene00006759_RNA_H1<-as.data.frame(WBGene00006759_RNA_H1)

WBGene00006759_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00006759")
WBGene00006759_RNA_H2<-WBGene00006759_RNA_H2[,4:13]
WBGene00006759_RNA_H2<-t(WBGene00006759_RNA_H2)
generations<-rownames(WBGene00006759_RNA_H2)
WBGene00006759_RNA_H2<-cbind(generations,WBGene00006759_RNA_H2)
rownames(WBGene00006759_RNA_H2)<-c()
colnames(WBGene00006759_RNA_H2)<-c("Generation","Zscore")
WBGene00006759_RNA_H2<-as.data.frame(WBGene00006759_RNA_H2)
WBGene00006759_RNA_H2$Zscore[WBGene00006759_RNA_H2$Zscore == 'NA']<- 0.2111528

plot(WBGene00006759_tRNA_H1$Generation,WBGene00006759_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Generation",ylab="Zscore",main="unc-22 (WBGene00006759)")
#lines(WBGene00006759_tRNA_H2$Generation,WBGene00006759_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00006759_RNA_H1$Generation,WBGene00006759_RNA_H1$Zscore, col="red",type="l")
#lines(WBGene00006759_RNA_H2$Generation,WBGene00006759_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(16, 5, legend=c("RNA", "Glu.CTC_tRNA"),
       col=c("red", "blue","black","black"), lty=c(1,1), cex=0.8)

#Get raw data
write.xlsx(WBGene00006759_tRNA_H1,"Data_Glu.CTC_tRNA.xlsx")
write.xlsx(WBGene00006759_RNA_H1,"Data_RNA.xlsx")

#gene 16
WBGene00014171_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[16,]
WBGene00014171_tRNA_H1<-WBGene00014171_tRNA_H1[,5:14]
WBGene00014171_tRNA_H1<-t(WBGene00014171_tRNA_H1)
generations<-rownames(WBGene00014171_tRNA_H1)
WBGene00014171_tRNA_H1<-cbind(generations,WBGene00014171_tRNA_H1)
colnames(WBGene00014171_tRNA_H1)<-c("Generation","Zscore")
WBGene00014171_tRNA_H1<-as.data.frame(WBGene00014171_tRNA_H1)
WBGene00014171_tRNA_H1$Generation<-as.numeric(WBGene00014171_tRNA_H1$Generation)
rownames(WBGene00014171_tRNA_H1)<-c()

WBGene00014171_tRNA_H2<-tRNA_z_scores_table_H2bis_WB[16,]
WBGene00014171_tRNA_H2<-WBGene00014171_tRNA_H2[,5:14]
WBGene00014171_tRNA_H2<-t(WBGene00014171_tRNA_H2)
generations<-rownames(WBGene00014171_tRNA_H2)
WBGene00014171_tRNA_H2<-cbind(generations,WBGene00014171_tRNA_H2)
colnames(WBGene00014171_tRNA_H2)<-c("Generation","Zscore")
WBGene00014171_tRNA_H2<-as.data.frame(WBGene00014171_tRNA_H2)
WBGene00014171_tRNA_H2$Generation<-as.numeric(WBGene00014171_tRNA_H2$Generation)
rownames(WBGene00014171_tRNA_H2)<-c()

WBGene00014171_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00014171")
WBGene00014171_RNA_H1<-WBGene00014171_RNA_H1[,4:13]
WBGene00014171_RNA_H1<-t(WBGene00014171_RNA_H1)
generations<-rownames(WBGene00014171_RNA_H1)
WBGene00014171_RNA_H1<-cbind(generations,WBGene00014171_RNA_H1)
rownames(WBGene00014171_RNA_H1)<-c()
colnames(WBGene00014171_RNA_H1)<-c("Generation","Zscore")
WBGene00014171_RNA_H1<-as.data.frame(WBGene00014171_RNA_H1)

WBGene00014171_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00014171")
WBGene00014171_RNA_H2<-WBGene00014171_RNA_H2[,4:13]
WBGene00014171_RNA_H2<-t(WBGene00014171_RNA_H2)
generations<-rownames(WBGene00014171_RNA_H2)
WBGene00014171_RNA_H2<-cbind(generations,WBGene00014171_RNA_H2)
rownames(WBGene00014171_RNA_H2)<-c()
colnames(WBGene00014171_RNA_H2)<-c("Generation","Zscore")
WBGene00014171_RNA_H2<-as.data.frame(WBGene00014171_RNA_H2)
WBGene00014171_RNA_H2$Zscore[WBGene00014171_RNA_H2$Zscore == 'NA']<- -0.6608825

plot(WBGene00014171_tRNA_H1$Generation,WBGene00014171_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="nep-260 (WBGene00014171)/Glu.TTC")
lines(WBGene00014171_tRNA_H2$Generation,WBGene00014171_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00014171_RNA_H1$Generation,WBGene00014171_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00014171_RNA_H2$Generation,WBGene00014171_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)

#gene 17
WBGene00011522_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[17,]
WBGene00011522_tRNA_H1<-WBGene00011522_tRNA_H1[,5:14]
WBGene00011522_tRNA_H1<-t(WBGene00011522_tRNA_H1)
generations<-rownames(WBGene00011522_tRNA_H1)
WBGene00011522_tRNA_H1<-cbind(generations,WBGene00011522_tRNA_H1)
colnames(WBGene00011522_tRNA_H1)<-c("Generation","Zscore")
WBGene00011522_tRNA_H1<-as.data.frame(WBGene00011522_tRNA_H1)
WBGene00011522_tRNA_H1$Generation<-as.numeric(WBGene00011522_tRNA_H1$Generation)
rownames(WBGene00011522_tRNA_H1)<-c()

WBGene00011522_tRNA_H2<-tRNA_z_scores_table_H2_WB[17,]
WBGene00011522_tRNA_H2<-WBGene00011522_tRNA_H2[,5:14]
WBGene00011522_tRNA_H2<-t(WBGene00011522_tRNA_H2)
generations<-rownames(WBGene00011522_tRNA_H2)
WBGene00011522_tRNA_H2<-cbind(generations,WBGene00011522_tRNA_H2)
colnames(WBGene00011522_tRNA_H2)<-c("Generation","Zscore")
WBGene00011522_tRNA_H2<-as.data.frame(WBGene00011522_tRNA_H2)
WBGene00011522_tRNA_H2$Generation<-as.numeric(WBGene00011522_tRNA_H2$Generation)
rownames(WBGene00011522_tRNA_H2)<-c()

WBGene00011522_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00011522")
WBGene00011522_RNA_H1<-WBGene00011522_RNA_H1[,4:13]
WBGene00011522_RNA_H1<-t(WBGene00011522_RNA_H1)
generations<-rownames(WBGene00011522_RNA_H1)
WBGene00011522_RNA_H1<-cbind(generations,WBGene00011522_RNA_H1)
rownames(WBGene00011522_RNA_H1)<-c()
colnames(WBGene00011522_RNA_H1)<-c("Generation","Zscore")
WBGene00011522_RNA_H1<-as.data.frame(WBGene00011522_RNA_H1)

WBGene00011522_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00011522")
WBGene00011522_RNA_H2<-WBGene00011522_RNA_H2[,4:13]
WBGene00011522_RNA_H2<-t(WBGene00011522_RNA_H2)
generations<-rownames(WBGene00011522_RNA_H2)
WBGene00011522_RNA_H2<-cbind(generations,WBGene00011522_RNA_H2)
rownames(WBGene00011522_RNA_H2)<-c()
colnames(WBGene00011522_RNA_H2)<-c("Generation","Zscore")
WBGene00011522_RNA_H2<-as.data.frame(WBGene00011522_RNA_H2)
WBGene00011522_RNA_H2$Zscore[WBGene00011522_RNA_H2$Zscore == 'NA']<- 0.7278732

plot(WBGene00011522_tRNA_H1$Generation,WBGene00011522_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="srap-1 (WBGene00011522)/Gly.GCC")
lines(WBGene00011522_tRNA_H2$Generation,WBGene00011522_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00011522_RNA_H1$Generation,WBGene00011522_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00011522_RNA_H2$Generation,WBGene00011522_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)

#gene 18
WBGene00011522_bis_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[18,]
WBGene00011522_bis_tRNA_H1<-WBGene00011522_bis_tRNA_H1[,5:14]
WBGene00011522_bis_tRNA_H1<-t(WBGene00011522_bis_tRNA_H1)
generations<-rownames(WBGene00011522_bis_tRNA_H1)
WBGene00011522_bis_tRNA_H1<-cbind(generations,WBGene00011522_bis_tRNA_H1)
colnames(WBGene00011522_bis_tRNA_H1)<-c("Generation","Zscore")
WBGene00011522_bis_tRNA_H1<-as.data.frame(WBGene00011522_bis_tRNA_H1)
WBGene00011522_bis_tRNA_H1$Generation<-as.numeric(WBGene00011522_bis_tRNA_H1$Generation)
rownames(WBGene00011522_bis_tRNA_H1)<-c()

WBGene00011522_bis_tRNA_H2<-tRNA_z_scores_table_H2_WB[18,]
WBGene00011522_bis_tRNA_H2<-WBGene00011522_bis_tRNA_H2[,5:14]
WBGene00011522_bis_tRNA_H2<-t(WBGene00011522_bis_tRNA_H2)
generations<-rownames(WBGene00011522_bis_tRNA_H2)
WBGene00011522_bis_tRNA_H2<-cbind(generations,WBGene00011522_bis_tRNA_H2)
colnames(WBGene00011522_bis_tRNA_H2)<-c("Generation","Zscore")
WBGene00011522_bis_tRNA_H2<-as.data.frame(WBGene00011522_bis_tRNA_H2)
WBGene00011522_bis_tRNA_H2$Generation<-as.numeric(WBGene00011522_bis_tRNA_H2$Generation)
rownames(WBGene00011522_bis_tRNA_H2)<-c()

WBGene00011522_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00011522")
WBGene00011522_RNA_H1<-WBGene00011522_RNA_H1[,4:13]
WBGene00011522_RNA_H1<-t(WBGene00011522_RNA_H1)
generations<-rownames(WBGene00011522_RNA_H1)
WBGene00011522_RNA_H1<-cbind(generations,WBGene00011522_RNA_H1)
rownames(WBGene00011522_RNA_H1)<-c()
colnames(WBGene00011522_RNA_H1)<-c("Generation","Zscore")
WBGene00011522_RNA_H1<-as.data.frame(WBGene00011522_RNA_H1)

WBGene00011522_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00011522")
WBGene00011522_RNA_H2<-WBGene00011522_RNA_H2[,4:13]
WBGene00011522_RNA_H2<-t(WBGene00011522_RNA_H2)
generations<-rownames(WBGene00011522_RNA_H2)
WBGene00011522_RNA_H2<-cbind(generations,WBGene00011522_RNA_H2)
rownames(WBGene00011522_RNA_H2)<-c()
colnames(WBGene00011522_RNA_H2)<-c("Generation","Zscore")
WBGene00011522_RNA_H2<-as.data.frame(WBGene00011522_RNA_H2)
WBGene00011522_RNA_H2$Zscore[WBGene00011522_RNA_H2$Zscore == 'NA']<- 0.7278732

plot(WBGene00011522_bis_tRNA_H1$Generation,WBGene00011522_bis_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="srap-1 (WBGene00011522)/Gly.GCC.2")
lines(WBGene00011522_bis_tRNA_H2$Generation,WBGene00011522_bis_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00011522_RNA_H1$Generation,WBGene00011522_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00011522_RNA_H2$Generation,WBGene00011522_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)

#gene 19
WBGene00011522_ter_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[19,]
WBGene00011522_ter_tRNA_H1<-WBGene00011522_ter_tRNA_H1[,5:14]
WBGene00011522_ter_tRNA_H1<-t(WBGene00011522_ter_tRNA_H1)
generations<-rownames(WBGene00011522_ter_tRNA_H1)
WBGene00011522_ter_tRNA_H1<-cbind(generations,WBGene00011522_ter_tRNA_H1)
colnames(WBGene00011522_ter_tRNA_H1)<-c("Generation","Zscore")
WBGene00011522_ter_tRNA_H1<-as.data.frame(WBGene00011522_ter_tRNA_H1)
WBGene00011522_ter_tRNA_H1$Generation<-as.numeric(WBGene00011522_ter_tRNA_H1$Generation)
rownames(WBGene00011522_ter_tRNA_H1)<-c()

WBGene00011522_ter_tRNA_H2<-tRNA_z_scores_table_H2_WB[19,]
WBGene00011522_ter_tRNA_H2<-WBGene00011522_ter_tRNA_H2[,5:14]
WBGene00011522_ter_tRNA_H2<-t(WBGene00011522_ter_tRNA_H2)
generations<-rownames(WBGene00011522_ter_tRNA_H2)
WBGene00011522_ter_tRNA_H2<-cbind(generations,WBGene00011522_ter_tRNA_H2)
colnames(WBGene00011522_ter_tRNA_H2)<-c("Generation","Zscore")
WBGene00011522_ter_tRNA_H2<-as.data.frame(WBGene00011522_ter_tRNA_H2)
WBGene00011522_ter_tRNA_H2$Generation<-as.numeric(WBGene00011522_ter_tRNA_H2$Generation)
rownames(WBGene00011522_ter_tRNA_H2)<-c()

WBGene00011522_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00011522")
WBGene00011522_RNA_H1<-WBGene00011522_RNA_H1[,4:13]
WBGene00011522_RNA_H1<-t(WBGene00011522_RNA_H1)
generations<-rownames(WBGene00011522_RNA_H1)
WBGene00011522_RNA_H1<-cbind(generations,WBGene00011522_RNA_H1)
rownames(WBGene00011522_RNA_H1)<-c()
colnames(WBGene00011522_RNA_H1)<-c("Generation","Zscore")
WBGene00011522_RNA_H1<-as.data.frame(WBGene00011522_RNA_H1)

WBGene00011522_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00011522")
WBGene00011522_RNA_H2<-WBGene00011522_RNA_H2[,4:13]
WBGene00011522_RNA_H2<-t(WBGene00011522_RNA_H2)
generations<-rownames(WBGene00011522_RNA_H2)
WBGene00011522_RNA_H2<-cbind(generations,WBGene00011522_RNA_H2)
rownames(WBGene00011522_RNA_H2)<-c()
colnames(WBGene00011522_RNA_H2)<-c("Generation","Zscore")
WBGene00011522_RNA_H2<-as.data.frame(WBGene00011522_RNA_H2)
WBGene00011522_RNA_H2$Zscore[WBGene00011522_RNA_H2$Zscore == 'NA']<- 0.7278732

plot(WBGene00011522_ter_tRNA_H1$Generation,WBGene00011522_ter_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="srap-1 (WBGene00011522)/Gly.GCC.3")
lines(WBGene00011522_ter_tRNA_H2$Generation,WBGene00011522_ter_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00011522_RNA_H1$Generation,WBGene00011522_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00011522_RNA_H2$Generation,WBGene00011522_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)

#gene 20
WBGene00011522_4_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[20,]
WBGene00011522_4_tRNA_H1<-WBGene00011522_4_tRNA_H1[,5:14]
WBGene00011522_4_tRNA_H1<-t(WBGene00011522_4_tRNA_H1)
generations<-rownames(WBGene00011522_4_tRNA_H1)
WBGene00011522_4_tRNA_H1<-cbind(generations,WBGene00011522_4_tRNA_H1)
colnames(WBGene00011522_4_tRNA_H1)<-c("Generation","Zscore")
WBGene00011522_4_tRNA_H1<-as.data.frame(WBGene00011522_4_tRNA_H1)
WBGene00011522_4_tRNA_H1$Generation<-as.numeric(WBGene00011522_4_tRNA_H1$Generation)
rownames(WBGene00011522_4_tRNA_H1)<-c()

WBGene00011522_4_tRNA_H2<-tRNA_z_scores_table_H2_WB[20,]
WBGene00011522_4_tRNA_H2<-WBGene00011522_4_tRNA_H2[,5:14]
WBGene00011522_4_tRNA_H2<-t(WBGene00011522_4_tRNA_H2)
generations<-rownames(WBGene00011522_4_tRNA_H2)
WBGene00011522_4_tRNA_H2<-cbind(generations,WBGene00011522_4_tRNA_H2)
colnames(WBGene00011522_4_tRNA_H2)<-c("Generation","Zscore")
WBGene00011522_4_tRNA_H2<-as.data.frame(WBGene00011522_4_tRNA_H2)
WBGene00011522_4_tRNA_H2$Generation<-as.numeric(WBGene00011522_4_tRNA_H2$Generation)
rownames(WBGene00011522_4_tRNA_H2)<-c()

WBGene00011522_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00011522")
WBGene00011522_RNA_H1<-WBGene00011522_RNA_H1[,4:13]
WBGene00011522_RNA_H1<-t(WBGene00011522_RNA_H1)
generations<-rownames(WBGene00011522_RNA_H1)
WBGene00011522_RNA_H1<-cbind(generations,WBGene00011522_RNA_H1)
rownames(WBGene00011522_RNA_H1)<-c()
colnames(WBGene00011522_RNA_H1)<-c("Generation","Zscore")
WBGene00011522_RNA_H1<-as.data.frame(WBGene00011522_RNA_H1)

WBGene00011522_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00011522")
WBGene00011522_RNA_H2<-WBGene00011522_RNA_H2[,4:13]
WBGene00011522_RNA_H2<-t(WBGene00011522_RNA_H2)
generations<-rownames(WBGene00011522_RNA_H2)
WBGene00011522_RNA_H2<-cbind(generations,WBGene00011522_RNA_H2)
rownames(WBGene00011522_RNA_H2)<-c()
colnames(WBGene00011522_RNA_H2)<-c("Generation","Zscore")
WBGene00011522_RNA_H2<-as.data.frame(WBGene00011522_RNA_H2)
WBGene00011522_RNA_H2$Zscore[WBGene00011522_RNA_H2$Zscore == 'NA']<- 0.7278732

plot(WBGene00011522_4_tRNA_H1$Generation,WBGene00011522_4_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="srap-1 (WBGene00011522)/Gly.GCC.4")
lines(WBGene00011522_4_tRNA_H2$Generation,WBGene00011522_4_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00011522_RNA_H1$Generation,WBGene00011522_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00011522_RNA_H2$Generation,WBGene00011522_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)

#gene 21
WBGene00011522_5_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[21,]
WBGene00011522_5_tRNA_H1<-WBGene00011522_5_tRNA_H1[,5:14]
WBGene00011522_5_tRNA_H1<-t(WBGene00011522_5_tRNA_H1)
generations<-rownames(WBGene00011522_5_tRNA_H1)
WBGene00011522_5_tRNA_H1<-cbind(generations,WBGene00011522_5_tRNA_H1)
colnames(WBGene00011522_5_tRNA_H1)<-c("Generation","Zscore")
WBGene00011522_5_tRNA_H1<-as.data.frame(WBGene00011522_5_tRNA_H1)
WBGene00011522_5_tRNA_H1$Generation<-as.numeric(WBGene00011522_5_tRNA_H1$Generation)
rownames(WBGene00011522_5_tRNA_H1)<-c()

WBGene00011522_5_tRNA_H2<-tRNA_z_scores_table_H2_WB[21,]
WBGene00011522_5_tRNA_H2<-WBGene00011522_5_tRNA_H2[,5:14]
WBGene00011522_5_tRNA_H2<-t(WBGene00011522_5_tRNA_H2)
generations<-rownames(WBGene00011522_5_tRNA_H2)
WBGene00011522_5_tRNA_H2<-cbind(generations,WBGene00011522_5_tRNA_H2)
colnames(WBGene00011522_5_tRNA_H2)<-c("Generation","Zscore")
WBGene00011522_5_tRNA_H2<-as.data.frame(WBGene00011522_5_tRNA_H2)
WBGene00011522_5_tRNA_H2$Generation<-as.numeric(WBGene00011522_5_tRNA_H2$Generation)
rownames(WBGene00011522_5_tRNA_H2)<-c()

WBGene00011522_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00011522")
WBGene00011522_RNA_H1<-WBGene00011522_RNA_H1[,4:13]
WBGene00011522_RNA_H1<-t(WBGene00011522_RNA_H1)
generations<-rownames(WBGene00011522_RNA_H1)
WBGene00011522_RNA_H1<-cbind(generations,WBGene00011522_RNA_H1)
rownames(WBGene00011522_RNA_H1)<-c()
colnames(WBGene00011522_RNA_H1)<-c("Generation","Zscore")
WBGene00011522_RNA_H1<-as.data.frame(WBGene00011522_RNA_H1)

WBGene00011522_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00011522")
WBGene00011522_RNA_H2<-WBGene00011522_RNA_H2[,4:13]
WBGene00011522_RNA_H2<-t(WBGene00011522_RNA_H2)
generations<-rownames(WBGene00011522_RNA_H2)
WBGene00011522_RNA_H2<-cbind(generations,WBGene00011522_RNA_H2)
rownames(WBGene00011522_RNA_H2)<-c()
colnames(WBGene00011522_RNA_H2)<-c("Generation","Zscore")
WBGene00011522_RNA_H2<-as.data.frame(WBGene00011522_RNA_H2)
WBGene00011522_RNA_H2$Zscore[WBGene00011522_RNA_H2$Zscore == 'NA']<- 0.7278732

plot(WBGene00011522_5_tRNA_H1$Generation,WBGene00011522_5_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="srap-1 (WBGene00011522)/Gly.GCC.5")
lines(WBGene00011522_5_tRNA_H2$Generation,WBGene00011522_5_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00011522_RNA_H1$Generation,WBGene00011522_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00011522_RNA_H2$Generation,WBGene00011522_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)

#gene 22
WBGene00011522_6_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[22,]
WBGene00011522_6_tRNA_H1<-WBGene00011522_6_tRNA_H1[,5:14]
WBGene00011522_6_tRNA_H1<-t(WBGene00011522_6_tRNA_H1)
generations<-rownames(WBGene00011522_6_tRNA_H1)
WBGene00011522_6_tRNA_H1<-cbind(generations,WBGene00011522_6_tRNA_H1)
colnames(WBGene00011522_6_tRNA_H1)<-c("Generation","Zscore")
WBGene00011522_6_tRNA_H1<-as.data.frame(WBGene00011522_6_tRNA_H1)
WBGene00011522_6_tRNA_H1$Generation<-as.numeric(WBGene00011522_6_tRNA_H1$Generation)
rownames(WBGene00011522_6_tRNA_H1)<-c()

WBGene00011522_6_tRNA_H2<-tRNA_z_scores_table_H2_WB[22,]
WBGene00011522_6_tRNA_H2<-WBGene00011522_6_tRNA_H2[,5:14]
WBGene00011522_6_tRNA_H2<-t(WBGene00011522_6_tRNA_H2)
generations<-rownames(WBGene00011522_6_tRNA_H2)
WBGene00011522_6_tRNA_H2<-cbind(generations,WBGene00011522_6_tRNA_H2)
colnames(WBGene00011522_6_tRNA_H2)<-c("Generation","Zscore")
WBGene00011522_6_tRNA_H2<-as.data.frame(WBGene00011522_6_tRNA_H2)
WBGene00011522_6_tRNA_H2$Generation<-as.numeric(WBGene00011522_6_tRNA_H2$Generation)
rownames(WBGene00011522_6_tRNA_H2)<-c()

WBGene00011522_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00011522")
WBGene00011522_RNA_H1<-WBGene00011522_RNA_H1[,4:13]
WBGene00011522_RNA_H1<-t(WBGene00011522_RNA_H1)
generations<-rownames(WBGene00011522_RNA_H1)
WBGene00011522_RNA_H1<-cbind(generations,WBGene00011522_RNA_H1)
rownames(WBGene00011522_RNA_H1)<-c()
colnames(WBGene00011522_RNA_H1)<-c("Generation","Zscore")
WBGene00011522_RNA_H1<-as.data.frame(WBGene00011522_RNA_H1)

WBGene00011522_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00011522")
WBGene00011522_RNA_H2<-WBGene00011522_RNA_H2[,4:13]
WBGene00011522_RNA_H2<-t(WBGene00011522_RNA_H2)
generations<-rownames(WBGene00011522_RNA_H2)
WBGene00011522_RNA_H2<-cbind(generations,WBGene00011522_RNA_H2)
rownames(WBGene00011522_RNA_H2)<-c()
colnames(WBGene00011522_RNA_H2)<-c("Generation","Zscore")
WBGene00011522_RNA_H2<-as.data.frame(WBGene00011522_RNA_H2)
WBGene00011522_RNA_H2$Zscore[WBGene00011522_RNA_H2$Zscore == 'NA']<- 0.7278732

plot(WBGene00011522_6_tRNA_H1$Generation,WBGene00011522_6_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="srap-1 (WBGene00011522)/Gly.GCC.6")
lines(WBGene00011522_6_tRNA_H2$Generation,WBGene00011522_6_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00011522_RNA_H1$Generation,WBGene00011522_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00011522_RNA_H2$Generation,WBGene00011522_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)


#gene 23
WBGene00016984_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[23,]
WBGene00016984_tRNA_H1<-WBGene00016984_tRNA_H1[,5:14]
WBGene00016984_tRNA_H1<-t(WBGene00016984_tRNA_H1)
generations<-rownames(WBGene00016984_tRNA_H1)
WBGene00016984_tRNA_H1<-cbind(generations,WBGene00016984_tRNA_H1)
colnames(WBGene00016984_tRNA_H1)<-c("Generation","Zscore")
WBGene00016984_tRNA_H1<-as.data.frame(WBGene00016984_tRNA_H1)
WBGene00016984_tRNA_H1$Generation<-as.numeric(WBGene00016984_tRNA_H1$Generation)
rownames(WBGene00016984_tRNA_H1)<-c()

WBGene00016984_tRNA_H2<-tRNA_z_scores_table_H2_WB[23,]
WBGene00016984_tRNA_H2<-WBGene00016984_tRNA_H2[,5:14]
WBGene00016984_tRNA_H2<-t(WBGene00016984_tRNA_H2)
generations<-rownames(WBGene00016984_tRNA_H2)
WBGene00016984_tRNA_H2<-cbind(generations,WBGene00016984_tRNA_H2)
colnames(WBGene00016984_tRNA_H2)<-c("Generation","Zscore")
WBGene00016984_tRNA_H2<-as.data.frame(WBGene00016984_tRNA_H2)
WBGene00016984_tRNA_H2$Generation<-as.numeric(WBGene00016984_tRNA_H2$Generation)
rownames(WBGene00016984_tRNA_H2)<-c()

WBGene00016984_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00016984")
WBGene00016984_RNA_H1<-WBGene00016984_RNA_H1[,4:13]
WBGene00016984_RNA_H1<-t(WBGene00016984_RNA_H1)
generations<-rownames(WBGene00016984_RNA_H1)
WBGene00016984_RNA_H1<-cbind(generations,WBGene00016984_RNA_H1)
rownames(WBGene00016984_RNA_H1)<-c()
colnames(WBGene00016984_RNA_H1)<-c("Generation","Zscore")
WBGene00016984_RNA_H1<-as.data.frame(WBGene00016984_RNA_H1)

WBGene00016984_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00016984")
WBGene00016984_RNA_H2<-WBGene00016984_RNA_H2[,4:13]
WBGene00016984_RNA_H2<-t(WBGene00016984_RNA_H2)
generations<-rownames(WBGene00016984_RNA_H2)
WBGene00016984_RNA_H2<-cbind(generations,WBGene00016984_RNA_H2)
rownames(WBGene00016984_RNA_H2)<-c()
colnames(WBGene00016984_RNA_H2)<-c("Generation","Zscore")
WBGene00016984_RNA_H2<-as.data.frame(WBGene00016984_RNA_H2)
WBGene00016984_RNA_H2$Zscore[WBGene00016984_RNA_H2$Zscore == 'NA']<- 0.08864671

plot(WBGene00016984_tRNA_H1$Generation,WBGene00016984_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="npr-8 (WBGene00016984)/Gly.TCC")
lines(WBGene00016984_tRNA_H2$Generation,WBGene00016984_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00016984_RNA_H1$Generation,WBGene00016984_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00016984_RNA_H2$Generation,WBGene00016984_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)


#gene 24
WBGene00003728_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[24,]
WBGene00003728_tRNA_H1<-WBGene00003728_tRNA_H1[,5:14]
WBGene00003728_tRNA_H1<-t(WBGene00003728_tRNA_H1)
generations<-rownames(WBGene00003728_tRNA_H1)
WBGene00003728_tRNA_H1<-cbind(generations,WBGene00003728_tRNA_H1)
colnames(WBGene00003728_tRNA_H1)<-c("Generation","Zscore")
WBGene00003728_tRNA_H1<-as.data.frame(WBGene00003728_tRNA_H1)
WBGene00003728_tRNA_H1$Generation<-as.numeric(WBGene00003728_tRNA_H1$Generation)
rownames(WBGene00003728_tRNA_H1)<-c()

WBGene00003728_tRNA_H2<-tRNA_z_scores_table_H2_WB[24,]
WBGene00003728_tRNA_H2<-WBGene00003728_tRNA_H2[,5:14]
WBGene00003728_tRNA_H2<-t(WBGene00003728_tRNA_H2)
generations<-rownames(WBGene00003728_tRNA_H2)
WBGene00003728_tRNA_H2<-cbind(generations,WBGene00003728_tRNA_H2)
colnames(WBGene00003728_tRNA_H2)<-c("Generation","Zscore")
WBGene00003728_tRNA_H2<-as.data.frame(WBGene00003728_tRNA_H2)
WBGene00003728_tRNA_H2$Generation<-as.numeric(WBGene00003728_tRNA_H2$Generation)
rownames(WBGene00003728_tRNA_H2)<-c()

WBGene00003728_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00003728")
WBGene00003728_RNA_H1<-WBGene00003728_RNA_H1[,4:13]
WBGene00003728_RNA_H1<-t(WBGene00003728_RNA_H1)
generations<-rownames(WBGene00003728_RNA_H1)
WBGene00003728_RNA_H1<-cbind(generations,WBGene00003728_RNA_H1)
rownames(WBGene00003728_RNA_H1)<-c()
colnames(WBGene00003728_RNA_H1)<-c("Generation","Zscore")
WBGene00003728_RNA_H1<-as.data.frame(WBGene00003728_RNA_H1)

WBGene00003728_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00003728")
WBGene00003728_RNA_H2<-WBGene00003728_RNA_H2[,4:13]
WBGene00003728_RNA_H2<-t(WBGene00003728_RNA_H2)
generations<-rownames(WBGene00003728_RNA_H2)
WBGene00003728_RNA_H2<-cbind(generations,WBGene00003728_RNA_H2)
rownames(WBGene00003728_RNA_H2)<-c()
colnames(WBGene00003728_RNA_H2)<-c("Generation","Zscore")
WBGene00003728_RNA_H2<-as.data.frame(WBGene00003728_RNA_H2)
WBGene00003728_RNA_H2$Zscore[WBGene00003728_RNA_H2$Zscore == 'NA']<- 0.2665558

plot(WBGene00003728_tRNA_H1$Generation,WBGene00003728_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="nhr-138 (WBGene00003728)/Gly.TCC")
lines(WBGene00003728_tRNA_H2$Generation,WBGene00003728_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00003728_RNA_H1$Generation,WBGene00003728_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00003728_RNA_H2$Generation,WBGene00003728_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)


#gene 25
WBGene00005502_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[25,]
WBGene00005502_tRNA_H1<-WBGene00005502_tRNA_H1[,5:14]
WBGene00005502_tRNA_H1<-t(WBGene00005502_tRNA_H1)
generations<-rownames(WBGene00005502_tRNA_H1)
WBGene00005502_tRNA_H1<-cbind(generations,WBGene00005502_tRNA_H1)
colnames(WBGene00005502_tRNA_H1)<-c("Generation","Zscore")
WBGene00005502_tRNA_H1<-as.data.frame(WBGene00005502_tRNA_H1)
WBGene00005502_tRNA_H1$Generation<-as.numeric(WBGene00005502_tRNA_H1$Generation)
rownames(WBGene00005502_tRNA_H1)<-c()

WBGene00005502_tRNA_H2<-tRNA_z_scores_table_H2_WB[25,]
WBGene00005502_tRNA_H2<-WBGene00005502_tRNA_H2[,5:14]
WBGene00005502_tRNA_H2<-t(WBGene00005502_tRNA_H2)
generations<-rownames(WBGene00005502_tRNA_H2)
WBGene00005502_tRNA_H2<-cbind(generations,WBGene00005502_tRNA_H2)
colnames(WBGene00005502_tRNA_H2)<-c("Generation","Zscore")
WBGene00005502_tRNA_H2<-as.data.frame(WBGene00005502_tRNA_H2)
WBGene00005502_tRNA_H2$Generation<-as.numeric(WBGene00005502_tRNA_H2$Generation)
rownames(WBGene00005502_tRNA_H2)<-c()

WBGene00005502_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00005502")
WBGene00005502_RNA_H1<-WBGene00005502_RNA_H1[,4:13]
WBGene00005502_RNA_H1<-t(WBGene00005502_RNA_H1)
generations<-rownames(WBGene00005502_RNA_H1)
WBGene00005502_RNA_H1<-cbind(generations,WBGene00005502_RNA_H1)
rownames(WBGene00005502_RNA_H1)<-c()
colnames(WBGene00005502_RNA_H1)<-c("Generation","Zscore")
WBGene00005502_RNA_H1<-as.data.frame(WBGene00005502_RNA_H1)

WBGene00005502_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00005502")
WBGene00005502_RNA_H2<-WBGene00005502_RNA_H2[,4:13]
WBGene00005502_RNA_H2<-t(WBGene00005502_RNA_H2)
generations<-rownames(WBGene00005502_RNA_H2)
WBGene00005502_RNA_H2<-cbind(generations,WBGene00005502_RNA_H2)
rownames(WBGene00005502_RNA_H2)<-c()
colnames(WBGene00005502_RNA_H2)<-c("Generation","Zscore")
WBGene00005502_RNA_H2<-as.data.frame(WBGene00005502_RNA_H2)
WBGene00005502_RNA_H2$Zscore[WBGene00005502_RNA_H2$Zscore == 'NA']<- -0.531126

plot(WBGene00005502_tRNA_H1$Generation,WBGene00005502_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="srh-298 (WBGene00005502)/His.GTG")
lines(WBGene00005502_tRNA_H2$Generation,WBGene00005502_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00005502_RNA_H1$Generation,WBGene00005502_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00005502_RNA_H2$Generation,WBGene00005502_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)



#gene 26
WBGene00022056_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[26,]
WBGene00022056_tRNA_H1<-WBGene00022056_tRNA_H1[,5:14]
WBGene00022056_tRNA_H1<-t(WBGene00022056_tRNA_H1)
generations<-rownames(WBGene00022056_tRNA_H1)
WBGene00022056_tRNA_H1<-cbind(generations,WBGene00022056_tRNA_H1)
colnames(WBGene00022056_tRNA_H1)<-c("Generation","Zscore")
WBGene00022056_tRNA_H1<-as.data.frame(WBGene00022056_tRNA_H1)
WBGene00022056_tRNA_H1$Generation<-as.numeric(WBGene00022056_tRNA_H1$Generation)
rownames(WBGene00022056_tRNA_H1)<-c()

WBGene00022056_tRNA_H2<-tRNA_z_scores_table_H2_WB[26,]
WBGene00022056_tRNA_H2<-WBGene00022056_tRNA_H2[,5:14]
WBGene00022056_tRNA_H2<-t(WBGene00022056_tRNA_H2)
generations<-rownames(WBGene00022056_tRNA_H2)
WBGene00022056_tRNA_H2<-cbind(generations,WBGene00022056_tRNA_H2)
colnames(WBGene00022056_tRNA_H2)<-c("Generation","Zscore")
WBGene00022056_tRNA_H2<-as.data.frame(WBGene00022056_tRNA_H2)
WBGene00022056_tRNA_H2$Generation<-as.numeric(WBGene00022056_tRNA_H2$Generation)
rownames(WBGene00022056_tRNA_H2)<-c()

WBGene00022056_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00022056")
WBGene00022056_RNA_H1<-WBGene00022056_RNA_H1[,4:13]
WBGene00022056_RNA_H1<-t(WBGene00022056_RNA_H1)
generations<-rownames(WBGene00022056_RNA_H1)
WBGene00022056_RNA_H1<-cbind(generations,WBGene00022056_RNA_H1)
rownames(WBGene00022056_RNA_H1)<-c()
colnames(WBGene00022056_RNA_H1)<-c("Generation","Zscore")
WBGene00022056_RNA_H1<-as.data.frame(WBGene00022056_RNA_H1)

WBGene00022056_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00022056")
WBGene00022056_RNA_H2<-WBGene00022056_RNA_H2[,4:13]
WBGene00022056_RNA_H2<-t(WBGene00022056_RNA_H2)
generations<-rownames(WBGene00022056_RNA_H2)
WBGene00022056_RNA_H2<-cbind(generations,WBGene00022056_RNA_H2)
rownames(WBGene00022056_RNA_H2)<-c()
colnames(WBGene00022056_RNA_H2)<-c("Generation","Zscore")
WBGene00022056_RNA_H2<-as.data.frame(WBGene00022056_RNA_H2)
WBGene00022056_RNA_H2$Zscore[WBGene00022056_RNA_H2$Zscore == 'NA']<- 0.6498078

plot(WBGene00022056_tRNA_H1$Generation,WBGene00022056_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="ddx-35 (WBGene00022056)/His.GTG")
lines(WBGene00022056_tRNA_H2$Generation,WBGene00022056_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00022056_RNA_H1$Generation,WBGene00022056_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00022056_RNA_H2$Generation,WBGene00022056_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)



#gene 27
WBGene00022654_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[27,]
WBGene00022654_tRNA_H1<-WBGene00022654_tRNA_H1[,5:14]
WBGene00022654_tRNA_H1<-t(WBGene00022654_tRNA_H1)
generations<-rownames(WBGene00022654_tRNA_H1)
WBGene00022654_tRNA_H1<-cbind(generations,WBGene00022654_tRNA_H1)
colnames(WBGene00022654_tRNA_H1)<-c("Generation","Zscore")
WBGene00022654_tRNA_H1<-as.data.frame(WBGene00022654_tRNA_H1)
WBGene00022654_tRNA_H1$Generation<-as.numeric(WBGene00022654_tRNA_H1$Generation)
rownames(WBGene00022654_tRNA_H1)<-c()

WBGene00022654_tRNA_H2<-tRNA_z_scores_table_H2_WB[27,]
WBGene00022654_tRNA_H2<-WBGene00022654_tRNA_H2[,5:14]
WBGene00022654_tRNA_H2<-t(WBGene00022654_tRNA_H2)
generations<-rownames(WBGene00022654_tRNA_H2)
WBGene00022654_tRNA_H2<-cbind(generations,WBGene00022654_tRNA_H2)
colnames(WBGene00022654_tRNA_H2)<-c("Generation","Zscore")
WBGene00022654_tRNA_H2<-as.data.frame(WBGene00022654_tRNA_H2)
WBGene00022654_tRNA_H2$Generation<-as.numeric(WBGene00022654_tRNA_H2$Generation)
rownames(WBGene00022654_tRNA_H2)<-c()

WBGene00022654_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00022654")
WBGene00022654_RNA_H1<-WBGene00022654_RNA_H1[,4:13]
WBGene00022654_RNA_H1<-t(WBGene00022654_RNA_H1)
generations<-rownames(WBGene00022654_RNA_H1)
WBGene00022654_RNA_H1<-cbind(generations,WBGene00022654_RNA_H1)
rownames(WBGene00022654_RNA_H1)<-c()
colnames(WBGene00022654_RNA_H1)<-c("Generation","Zscore")
WBGene00022654_RNA_H1<-as.data.frame(WBGene00022654_RNA_H1)

WBGene00022654_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00022654")
WBGene00022654_RNA_H2<-WBGene00022654_RNA_H2[,4:13]
WBGene00022654_RNA_H2<-t(WBGene00022654_RNA_H2)
generations<-rownames(WBGene00022654_RNA_H2)
WBGene00022654_RNA_H2<-cbind(generations,WBGene00022654_RNA_H2)
rownames(WBGene00022654_RNA_H2)<-c()
colnames(WBGene00022654_RNA_H2)<-c("Generation","Zscore")
WBGene00022654_RNA_H2<-as.data.frame(WBGene00022654_RNA_H2)
WBGene00022654_RNA_H2$Zscore[WBGene00022654_RNA_H2$Zscore == 'NA']<- -0.3026966

plot(WBGene00022654_tRNA_H2$Generation,WBGene00022654_tRNA_H2$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Generation",ylab="Zscore",main="ZK105.3 (WBGene00022654)")
lines(WBGene00022654_RNA_H2$Generation,WBGene00022654_RNA_H2$Zscore, col="red",type="l",lty=1)
lines(WBGene00022654_bis_tRNA_H2$Generation,WBGene00022654_bis_tRNA_H2$Zscore, col="cadetblue",type="l",lty=1)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(15, 4.7, legend=c("RNA", "Leu.AAG.1_tRNA","Leu.AAG.2_tRNA"),
       col=c("red", "blue","cadetblue"), lty=c(1,1,1), cex=0.8)

#Get the raw data
write.xlsx(WBGene00022654_tRNA_H2,"Data_Leu.AAG.1_tRNA.xlsx")
write.xlsx(WBGene00022654_RNA_H2,"Data_RNA.xlsx")
write.xlsx(WBGene00022654_bis_tRNA_H2,"Data_Leu.AAG.2_tRNA.xlsx")

#gene 28
WBGene00019600_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[28,]
WBGene00019600_tRNA_H1<-WBGene00019600_tRNA_H1[,5:14]
WBGene00019600_tRNA_H1<-t(WBGene00019600_tRNA_H1)
generations<-rownames(WBGene00019600_tRNA_H1)
WBGene00019600_tRNA_H1<-cbind(generations,WBGene00019600_tRNA_H1)
colnames(WBGene00019600_tRNA_H1)<-c("Generation","Zscore")
WBGene00019600_tRNA_H1<-as.data.frame(WBGene00019600_tRNA_H1)
WBGene00019600_tRNA_H1$Generation<-as.numeric(WBGene00019600_tRNA_H1$Generation)
rownames(WBGene00019600_tRNA_H1)<-c()

WBGene00019600_tRNA_H2<-tRNA_z_scores_table_H2_WB[28,]
WBGene00019600_tRNA_H2<-WBGene00019600_tRNA_H2[,5:14]
WBGene00019600_tRNA_H2<-t(WBGene00019600_tRNA_H2)
generations<-rownames(WBGene00019600_tRNA_H2)
WBGene00019600_tRNA_H2<-cbind(generations,WBGene00019600_tRNA_H2)
colnames(WBGene00019600_tRNA_H2)<-c("Generation","Zscore")
WBGene00019600_tRNA_H2<-as.data.frame(WBGene00019600_tRNA_H2)
WBGene00019600_tRNA_H2$Generation<-as.numeric(WBGene00019600_tRNA_H2$Generation)
rownames(WBGene00019600_tRNA_H2)<-c()

WBGene00019600_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00019600")
WBGene00019600_RNA_H1<-WBGene00019600_RNA_H1[,4:13]
WBGene00019600_RNA_H1<-t(WBGene00019600_RNA_H1)
generations<-rownames(WBGene00019600_RNA_H1)
WBGene00019600_RNA_H1<-cbind(generations,WBGene00019600_RNA_H1)
rownames(WBGene00019600_RNA_H1)<-c()
colnames(WBGene00019600_RNA_H1)<-c("Generation","Zscore")
WBGene00019600_RNA_H1<-as.data.frame(WBGene00019600_RNA_H1)

WBGene00019600_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00019600")
WBGene00019600_RNA_H2<-WBGene00019600_RNA_H2[,4:13]
WBGene00019600_RNA_H2<-t(WBGene00019600_RNA_H2)
generations<-rownames(WBGene00019600_RNA_H2)
WBGene00019600_RNA_H2<-cbind(generations,WBGene00019600_RNA_H2)
rownames(WBGene00019600_RNA_H2)<-c()
colnames(WBGene00019600_RNA_H2)<-c("Generation","Zscore")
WBGene00019600_RNA_H2<-as.data.frame(WBGene00019600_RNA_H2)
WBGene00019600_RNA_H2$Zscore[WBGene00019600_RNA_H2$Zscore == 'NA']<- 1.719649

plot(WBGene00019600_tRNA_H1$Generation,WBGene00019600_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="rga-3 (WBGene00019600)/Leu.AAG")
lines(WBGene00019600_tRNA_H2$Generation,WBGene00019600_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00019600_RNA_H1$Generation,WBGene00019600_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00019600_RNA_H2$Generation,WBGene00019600_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)



#gene 29
WBGene00022654_bis_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[29,]
WBGene00022654_bis_tRNA_H1<-WBGene00022654_bis_tRNA_H1[,5:14]
WBGene00022654_bis_tRNA_H1<-t(WBGene00022654_bis_tRNA_H1)
generations<-rownames(WBGene00022654_bis_tRNA_H1)
WBGene00022654_bis_tRNA_H1<-cbind(generations,WBGene00022654_bis_tRNA_H1)
colnames(WBGene00022654_bis_tRNA_H1)<-c("Generation","Zscore")
WBGene00022654_bis_tRNA_H1<-as.data.frame(WBGene00022654_bis_tRNA_H1)
WBGene00022654_bis_tRNA_H1$Generation<-as.numeric(WBGene00022654_bis_tRNA_H1$Generation)
rownames(WBGene00022654_bis_tRNA_H1)<-c()

WBGene00022654_bis_tRNA_H2<-tRNA_z_scores_table_H2_WB[29,]
WBGene00022654_bis_tRNA_H2<-WBGene00022654_bis_tRNA_H2[,5:14]
WBGene00022654_bis_tRNA_H2<-t(WBGene00022654_bis_tRNA_H2)
generations<-rownames(WBGene00022654_bis_tRNA_H2)
WBGene00022654_bis_tRNA_H2<-cbind(generations,WBGene00022654_bis_tRNA_H2)
colnames(WBGene00022654_bis_tRNA_H2)<-c("Generation","Zscore")
WBGene00022654_bis_tRNA_H2<-as.data.frame(WBGene00022654_bis_tRNA_H2)
WBGene00022654_bis_tRNA_H2$Generation<-as.numeric(WBGene00022654_bis_tRNA_H2$Generation)
rownames(WBGene00022654_bis_tRNA_H2)<-c()

WBGene00022654_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00022654")
WBGene00022654_RNA_H1<-WBGene00022654_RNA_H1[,4:13]
WBGene00022654_RNA_H1<-t(WBGene00022654_RNA_H1)
generations<-rownames(WBGene00022654_RNA_H1)
WBGene00022654_RNA_H1<-cbind(generations,WBGene00022654_RNA_H1)
rownames(WBGene00022654_RNA_H1)<-c()
colnames(WBGene00022654_RNA_H1)<-c("Generation","Zscore")
WBGene00022654_RNA_H1<-as.data.frame(WBGene00022654_RNA_H1)

WBGene00022654_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00022654")
WBGene00022654_RNA_H2<-WBGene00022654_RNA_H2[,4:13]
WBGene00022654_RNA_H2<-t(WBGene00022654_RNA_H2)
generations<-rownames(WBGene00022654_RNA_H2)
WBGene00022654_RNA_H2<-cbind(generations,WBGene00022654_RNA_H2)
rownames(WBGene00022654_RNA_H2)<-c()
colnames(WBGene00022654_RNA_H2)<-c("Generation","Zscore")
WBGene00022654_RNA_H2<-as.data.frame(WBGene00022654_RNA_H2)
WBGene00022654_RNA_H2$Zscore[WBGene00022654_RNA_H2$Zscore == 'NA']<- -0.3026966

plot(WBGene00022654_bis_tRNA_H1$Generation,WBGene00022654_bis_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="ZK105.3 (WBGene00022654)/Leu.AAG.2")
lines(WBGene00022654_bis_tRNA_H2$Generation,WBGene00022654_bis_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00022654_RNA_H1$Generation,WBGene00022654_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00022654_RNA_H2$Generation,WBGene00022654_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)



#gene 30
WBGene00019600_bis_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[30,]
WBGene00019600_bis_tRNA_H1<-WBGene00019600_bis_tRNA_H1[,5:14]
WBGene00019600_bis_tRNA_H1<-t(WBGene00019600_bis_tRNA_H1)
generations<-rownames(WBGene00019600_bis_tRNA_H1)
WBGene00019600_bis_tRNA_H1<-cbind(generations,WBGene00019600_bis_tRNA_H1)
colnames(WBGene00019600_bis_tRNA_H1)<-c("Generation","Zscore")
WBGene00019600_bis_tRNA_H1<-as.data.frame(WBGene00019600_bis_tRNA_H1)
WBGene00019600_bis_tRNA_H1$Generation<-as.numeric(WBGene00019600_bis_tRNA_H1$Generation)
rownames(WBGene00019600_bis_tRNA_H1)<-c()

WBGene00019600_bis_tRNA_H2<-tRNA_z_scores_table_H2_WB[30,]
WBGene00019600_bis_tRNA_H2<-WBGene00019600_bis_tRNA_H2[,5:14]
WBGene00019600_bis_tRNA_H2<-t(WBGene00019600_bis_tRNA_H2)
generations<-rownames(WBGene00019600_bis_tRNA_H2)
WBGene00019600_bis_tRNA_H2<-cbind(generations,WBGene00019600_bis_tRNA_H2)
colnames(WBGene00019600_bis_tRNA_H2)<-c("Generation","Zscore")
WBGene00019600_bis_tRNA_H2<-as.data.frame(WBGene00019600_bis_tRNA_H2)
WBGene00019600_bis_tRNA_H2$Generation<-as.numeric(WBGene00019600_bis_tRNA_H2$Generation)
rownames(WBGene00019600_bis_tRNA_H2)<-c()

WBGene00019600_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00019600")
WBGene00019600_RNA_H1<-WBGene00019600_RNA_H1[,4:13]
WBGene00019600_RNA_H1<-t(WBGene00019600_RNA_H1)
generations<-rownames(WBGene00019600_RNA_H1)
WBGene00019600_RNA_H1<-cbind(generations,WBGene00019600_RNA_H1)
rownames(WBGene00019600_RNA_H1)<-c()
colnames(WBGene00019600_RNA_H1)<-c("Generation","Zscore")
WBGene00019600_RNA_H1<-as.data.frame(WBGene00019600_RNA_H1)

WBGene00019600_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00019600")
WBGene00019600_RNA_H2<-WBGene00019600_RNA_H2[,4:13]
WBGene00019600_RNA_H2<-t(WBGene00019600_RNA_H2)
generations<-rownames(WBGene00019600_RNA_H2)
WBGene00019600_RNA_H2<-cbind(generations,WBGene00019600_RNA_H2)
rownames(WBGene00019600_RNA_H2)<-c()
colnames(WBGene00019600_RNA_H2)<-c("Generation","Zscore")
WBGene00019600_RNA_H2<-as.data.frame(WBGene00019600_RNA_H2)
WBGene00019600_RNA_H2$Zscore[WBGene00019600_RNA_H2$Zscore == 'NA']<- 1.719649

plot(WBGene00019600_bis_tRNA_H1$Generation,WBGene00019600_bis_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="rga-3 (WBGene00019600)/Leu.AAG.2")
lines(WBGene00019600_bis_tRNA_H2$Generation,WBGene00019600_bis_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00019600_RNA_H1$Generation,WBGene00019600_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00019600_RNA_H2$Generation,WBGene00019600_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)



#gene 31
WBGene00022654_ter_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[31,]
WBGene00022654_ter_tRNA_H1<-WBGene00022654_ter_tRNA_H1[,5:14]
WBGene00022654_ter_tRNA_H1<-t(WBGene00022654_ter_tRNA_H1)
generations<-rownames(WBGene00022654_ter_tRNA_H1)
WBGene00022654_ter_tRNA_H1<-cbind(generations,WBGene00022654_ter_tRNA_H1)
colnames(WBGene00022654_ter_tRNA_H1)<-c("Generation","Zscore")
WBGene00022654_ter_tRNA_H1<-as.data.frame(WBGene00022654_ter_tRNA_H1)
WBGene00022654_ter_tRNA_H1$Generation<-as.numeric(WBGene00022654_ter_tRNA_H1$Generation)
rownames(WBGene00022654_ter_tRNA_H1)<-c()

WBGene00022654_ter_tRNA_H2<-tRNA_z_scores_table_H2_WB[31,]
WBGene00022654_ter_tRNA_H2<-WBGene00022654_ter_tRNA_H2[,5:14]
WBGene00022654_ter_tRNA_H2<-t(WBGene00022654_ter_tRNA_H2)
generations<-rownames(WBGene00022654_ter_tRNA_H2)
WBGene00022654_ter_tRNA_H2<-cbind(generations,WBGene00022654_ter_tRNA_H2)
colnames(WBGene00022654_ter_tRNA_H2)<-c("Generation","Zscore")
WBGene00022654_ter_tRNA_H2<-as.data.frame(WBGene00022654_ter_tRNA_H2)
WBGene00022654_ter_tRNA_H2$Generation<-as.numeric(WBGene00022654_ter_tRNA_H2$Generation)
rownames(WBGene00022654_ter_tRNA_H2)<-c()

WBGene00022654_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00022654")
WBGene00022654_RNA_H1<-WBGene00022654_RNA_H1[,4:13]
WBGene00022654_RNA_H1<-t(WBGene00022654_RNA_H1)
generations<-rownames(WBGene00022654_RNA_H1)
WBGene00022654_RNA_H1<-cbind(generations,WBGene00022654_RNA_H1)
rownames(WBGene00022654_RNA_H1)<-c()
colnames(WBGene00022654_RNA_H1)<-c("Generation","Zscore")
WBGene00022654_RNA_H1<-as.data.frame(WBGene00022654_RNA_H1)

WBGene00022654_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00022654")
WBGene00022654_RNA_H2<-WBGene00022654_RNA_H2[,4:13]
WBGene00022654_RNA_H2<-t(WBGene00022654_RNA_H2)
generations<-rownames(WBGene00022654_RNA_H2)
WBGene00022654_RNA_H2<-cbind(generations,WBGene00022654_RNA_H2)
rownames(WBGene00022654_RNA_H2)<-c()
colnames(WBGene00022654_RNA_H2)<-c("Generation","Zscore")
WBGene00022654_RNA_H2<-as.data.frame(WBGene00022654_RNA_H2)
WBGene00022654_RNA_H2$Zscore[WBGene00022654_RNA_H2$Zscore == 'NA']<- -0.3026966

plot(WBGene00022654_ter_tRNA_H1$Generation,WBGene00022654_ter_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="ZK105.3 (WBGene00022654)/Leu.AAG.3")
lines(WBGene00022654_ter_tRNA_H2$Generation,WBGene00022654_ter_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00022654_RNA_H1$Generation,WBGene00022654_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00022654_RNA_H2$Generation,WBGene00022654_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)



#gene 32
WBGene00019600_ter_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[32,]
WBGene00019600_ter_tRNA_H1<-WBGene00019600_ter_tRNA_H1[,5:14]
WBGene00019600_ter_tRNA_H1<-t(WBGene00019600_ter_tRNA_H1)
generations<-rownames(WBGene00019600_ter_tRNA_H1)
WBGene00019600_ter_tRNA_H1<-cbind(generations,WBGene00019600_ter_tRNA_H1)
colnames(WBGene00019600_ter_tRNA_H1)<-c("Generation","Zscore")
WBGene00019600_ter_tRNA_H1<-as.data.frame(WBGene00019600_ter_tRNA_H1)
WBGene00019600_ter_tRNA_H1$Generation<-as.numeric(WBGene00019600_ter_tRNA_H1$Generation)
rownames(WBGene00019600_ter_tRNA_H1)<-c()

WBGene00019600_ter_tRNA_H2<-tRNA_z_scores_table_H2_WB[32,]
WBGene00019600_ter_tRNA_H2<-WBGene00019600_ter_tRNA_H2[,5:14]
WBGene00019600_ter_tRNA_H2<-t(WBGene00019600_ter_tRNA_H2)
generations<-rownames(WBGene00019600_ter_tRNA_H2)
WBGene00019600_ter_tRNA_H2<-cbind(generations,WBGene00019600_ter_tRNA_H2)
colnames(WBGene00019600_ter_tRNA_H2)<-c("Generation","Zscore")
WBGene00019600_ter_tRNA_H2<-as.data.frame(WBGene00019600_ter_tRNA_H2)
WBGene00019600_ter_tRNA_H2$Generation<-as.numeric(WBGene00019600_ter_tRNA_H2$Generation)
rownames(WBGene00019600_ter_tRNA_H2)<-c()

WBGene00019600_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00019600")
WBGene00019600_RNA_H1<-WBGene00019600_RNA_H1[,4:13]
WBGene00019600_RNA_H1<-t(WBGene00019600_RNA_H1)
generations<-rownames(WBGene00019600_RNA_H1)
WBGene00019600_RNA_H1<-cbind(generations,WBGene00019600_RNA_H1)
rownames(WBGene00019600_RNA_H1)<-c()
colnames(WBGene00019600_RNA_H1)<-c("Generation","Zscore")
WBGene00019600_RNA_H1<-as.data.frame(WBGene00019600_RNA_H1)

WBGene00019600_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00019600")
WBGene00019600_RNA_H2<-WBGene00019600_RNA_H2[,4:13]
WBGene00019600_RNA_H2<-t(WBGene00019600_RNA_H2)
generations<-rownames(WBGene00019600_RNA_H2)
WBGene00019600_RNA_H2<-cbind(generations,WBGene00019600_RNA_H2)
rownames(WBGene00019600_RNA_H2)<-c()
colnames(WBGene00019600_RNA_H2)<-c("Generation","Zscore")
WBGene00019600_RNA_H2<-as.data.frame(WBGene00019600_RNA_H2)
WBGene00019600_RNA_H2$Zscore[WBGene00019600_RNA_H2$Zscore == 'NA']<- 1.719649

plot(WBGene00019600_ter_tRNA_H1$Generation,WBGene00019600_ter_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="rga-3 (WBGene00019600)/Leu.AAG.3")
lines(WBGene00019600_ter_tRNA_H2$Generation,WBGene00019600_ter_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00019600_RNA_H1$Generation,WBGene00019600_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00019600_RNA_H2$Generation,WBGene00019600_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)



#gene 33
WBGene00004416_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[33,]
WBGene00004416_tRNA_H1<-WBGene00004416_tRNA_H1[,5:14]
WBGene00004416_tRNA_H1<-t(WBGene00004416_tRNA_H1)
generations<-rownames(WBGene00004416_tRNA_H1)
WBGene00004416_tRNA_H1<-cbind(generations,WBGene00004416_tRNA_H1)
colnames(WBGene00004416_tRNA_H1)<-c("Generation","Zscore")
WBGene00004416_tRNA_H1<-as.data.frame(WBGene00004416_tRNA_H1)
WBGene00004416_tRNA_H1$Generation<-as.numeric(WBGene00004416_tRNA_H1$Generation)
rownames(WBGene00004416_tRNA_H1)<-c()

WBGene00004416_tRNA_H2<-tRNA_z_scores_table_H2_WB[33,]
WBGene00004416_tRNA_H2<-WBGene00004416_tRNA_H2[,5:14]
WBGene00004416_tRNA_H2<-t(WBGene00004416_tRNA_H2)
generations<-rownames(WBGene00004416_tRNA_H2)
WBGene00004416_tRNA_H2<-cbind(generations,WBGene00004416_tRNA_H2)
colnames(WBGene00004416_tRNA_H2)<-c("Generation","Zscore")
WBGene00004416_tRNA_H2<-as.data.frame(WBGene00004416_tRNA_H2)
WBGene00004416_tRNA_H2$Generation<-as.numeric(WBGene00004416_tRNA_H2$Generation)
rownames(WBGene00004416_tRNA_H2)<-c()

WBGene00004416_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00004416")
WBGene00004416_RNA_H1<-WBGene00004416_RNA_H1[,4:13]
WBGene00004416_RNA_H1<-t(WBGene00004416_RNA_H1)
generations<-rownames(WBGene00004416_RNA_H1)
WBGene00004416_RNA_H1<-cbind(generations,WBGene00004416_RNA_H1)
rownames(WBGene00004416_RNA_H1)<-c()
colnames(WBGene00004416_RNA_H1)<-c("Generation","Zscore")
WBGene00004416_RNA_H1<-as.data.frame(WBGene00004416_RNA_H1)

WBGene00004416_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00004416")
WBGene00004416_RNA_H2<-WBGene00004416_RNA_H2[,4:13]
WBGene00004416_RNA_H2<-t(WBGene00004416_RNA_H2)
generations<-rownames(WBGene00004416_RNA_H2)
WBGene00004416_RNA_H2<-cbind(generations,WBGene00004416_RNA_H2)
rownames(WBGene00004416_RNA_H2)<-c()
colnames(WBGene00004416_RNA_H2)<-c("Generation","Zscore")
WBGene00004416_RNA_H2<-as.data.frame(WBGene00004416_RNA_H2)
WBGene00004416_RNA_H2$Zscore[WBGene00004416_RNA_H2$Zscore == 'NA']<- -0.4744733

plot(WBGene00004416_tRNA_H1$Generation,WBGene00004416_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="rpl-5 (WBGene00004416)/Leu.CAA")
lines(WBGene00004416_tRNA_H2$Generation,WBGene00004416_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00004416_RNA_H1$Generation,WBGene00004416_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00004416_RNA_H2$Generation,WBGene00004416_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)



#gene 34
WBGene00008215_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[34,]
WBGene00008215_tRNA_H1<-WBGene00008215_tRNA_H1[,5:14]
WBGene00008215_tRNA_H1<-t(WBGene00008215_tRNA_H1)
generations<-rownames(WBGene00008215_tRNA_H1)
WBGene00008215_tRNA_H1<-cbind(generations,WBGene00008215_tRNA_H1)
colnames(WBGene00008215_tRNA_H1)<-c("Generation","Zscore")
WBGene00008215_tRNA_H1<-as.data.frame(WBGene00008215_tRNA_H1)
WBGene00008215_tRNA_H1$Generation<-as.numeric(WBGene00008215_tRNA_H1$Generation)
rownames(WBGene00008215_tRNA_H1)<-c()

WBGene00008215_tRNA_H2<-tRNA_z_scores_table_H2_WB[34,]
WBGene00008215_tRNA_H2<-WBGene00008215_tRNA_H2[,5:14]
WBGene00008215_tRNA_H2<-t(WBGene00008215_tRNA_H2)
generations<-rownames(WBGene00008215_tRNA_H2)
WBGene00008215_tRNA_H2<-cbind(generations,WBGene00008215_tRNA_H2)
colnames(WBGene00008215_tRNA_H2)<-c("Generation","Zscore")
WBGene00008215_tRNA_H2<-as.data.frame(WBGene00008215_tRNA_H2)
WBGene00008215_tRNA_H2$Generation<-as.numeric(WBGene00008215_tRNA_H2$Generation)
rownames(WBGene00008215_tRNA_H2)<-c()

WBGene00008215_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00008215")
WBGene00008215_RNA_H1<-WBGene00008215_RNA_H1[,4:13]
WBGene00008215_RNA_H1<-t(WBGene00008215_RNA_H1)
generations<-rownames(WBGene00008215_RNA_H1)
WBGene00008215_RNA_H1<-cbind(generations,WBGene00008215_RNA_H1)
rownames(WBGene00008215_RNA_H1)<-c()
colnames(WBGene00008215_RNA_H1)<-c("Generation","Zscore")
WBGene00008215_RNA_H1<-as.data.frame(WBGene00008215_RNA_H1)

WBGene00008215_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00008215")
WBGene00008215_RNA_H2<-WBGene00008215_RNA_H2[,4:13]
WBGene00008215_RNA_H2<-t(WBGene00008215_RNA_H2)
generations<-rownames(WBGene00008215_RNA_H2)
WBGene00008215_RNA_H2<-cbind(generations,WBGene00008215_RNA_H2)
rownames(WBGene00008215_RNA_H2)<-c()
colnames(WBGene00008215_RNA_H2)<-c("Generation","Zscore")
WBGene00008215_RNA_H2<-as.data.frame(WBGene00008215_RNA_H2)
WBGene00008215_RNA_H2$Zscore[WBGene00008215_RNA_H2$Zscore == 'NA']<- -1.137299

plot(WBGene00008215_tRNA_H1$Generation,WBGene00008215_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="C49F8.3 (WBGene00008215)/Lys.CTT")
lines(WBGene00008215_tRNA_H2$Generation,WBGene00008215_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00008215_RNA_H1$Generation,WBGene00008215_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00008215_RNA_H2$Generation,WBGene00008215_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)



#gene 35
WBGene00009285_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[35,]
WBGene00009285_tRNA_H1<-WBGene00009285_tRNA_H1[,5:14]
WBGene00009285_tRNA_H1<-t(WBGene00009285_tRNA_H1)
generations<-rownames(WBGene00009285_tRNA_H1)
WBGene00009285_tRNA_H1<-cbind(generations,WBGene00009285_tRNA_H1)
colnames(WBGene00009285_tRNA_H1)<-c("Generation","Zscore")
WBGene00009285_tRNA_H1<-as.data.frame(WBGene00009285_tRNA_H1)
WBGene00009285_tRNA_H1$Generation<-as.numeric(WBGene00009285_tRNA_H1$Generation)
rownames(WBGene00009285_tRNA_H1)<-c()

WBGene00009285_tRNA_H2<-tRNA_z_scores_table_H2_WB[35,]
WBGene00009285_tRNA_H2<-WBGene00009285_tRNA_H2[,5:14]
WBGene00009285_tRNA_H2<-t(WBGene00009285_tRNA_H2)
generations<-rownames(WBGene00009285_tRNA_H2)
WBGene00009285_tRNA_H2<-cbind(generations,WBGene00009285_tRNA_H2)
colnames(WBGene00009285_tRNA_H2)<-c("Generation","Zscore")
WBGene00009285_tRNA_H2<-as.data.frame(WBGene00009285_tRNA_H2)
WBGene00009285_tRNA_H2$Generation<-as.numeric(WBGene00009285_tRNA_H2$Generation)
rownames(WBGene00009285_tRNA_H2)<-c()

WBGene00009285_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00009285")
WBGene00009285_RNA_H1<-WBGene00009285_RNA_H1[,4:13]
WBGene00009285_RNA_H1<-t(WBGene00009285_RNA_H1)
generations<-rownames(WBGene00009285_RNA_H1)
WBGene00009285_RNA_H1<-cbind(generations,WBGene00009285_RNA_H1)
rownames(WBGene00009285_RNA_H1)<-c()
colnames(WBGene00009285_RNA_H1)<-c("Generation","Zscore")
WBGene00009285_RNA_H1<-as.data.frame(WBGene00009285_RNA_H1)

WBGene00009285_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00009285")
WBGene00009285_RNA_H2<-WBGene00009285_RNA_H2[,4:13]
WBGene00009285_RNA_H2<-t(WBGene00009285_RNA_H2)
generations<-rownames(WBGene00009285_RNA_H2)
WBGene00009285_RNA_H2<-cbind(generations,WBGene00009285_RNA_H2)
rownames(WBGene00009285_RNA_H2)<-c()
colnames(WBGene00009285_RNA_H2)<-c("Generation","Zscore")
WBGene00009285_RNA_H2<-as.data.frame(WBGene00009285_RNA_H2)
WBGene00009285_RNA_H2$Zscore[WBGene00009285_RNA_H2$Zscore == 'NA']<- 2.090692

plot(WBGene00009285_tRNA_H1$Generation,WBGene00009285_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="F31C3.3 (WBGene00009285)/Lys.CTT")
lines(WBGene00009285_tRNA_H2$Generation,WBGene00009285_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00009285_RNA_H1$Generation,WBGene00009285_RNA_H2$Zscore, col="red",type="l")
lines(WBGene00009285_RNA_H2$Generation,WBGene00009285_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)



#gene 36
WBGene00017301_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[36,]
WBGene00017301_tRNA_H1<-WBGene00017301_tRNA_H1[,5:14]
WBGene00017301_tRNA_H1<-t(WBGene00017301_tRNA_H1)
generations<-rownames(WBGene00017301_tRNA_H1)
WBGene00017301_tRNA_H1<-cbind(generations,WBGene00017301_tRNA_H1)
colnames(WBGene00017301_tRNA_H1)<-c("Generation","Zscore")
WBGene00017301_tRNA_H1<-as.data.frame(WBGene00017301_tRNA_H1)
WBGene00017301_tRNA_H1$Generation<-as.numeric(WBGene00017301_tRNA_H1$Generation)
rownames(WBGene00017301_tRNA_H1)<-c()

WBGene00017301_tRNA_H2<-tRNA_z_scores_table_H2_WB[36,]
WBGene00017301_tRNA_H2<-WBGene00017301_tRNA_H2[,5:14]
WBGene00017301_tRNA_H2<-t(WBGene00017301_tRNA_H2)
generations<-rownames(WBGene00017301_tRNA_H2)
WBGene00017301_tRNA_H2<-cbind(generations,WBGene00017301_tRNA_H2)
colnames(WBGene00017301_tRNA_H2)<-c("Generation","Zscore")
WBGene00017301_tRNA_H2<-as.data.frame(WBGene00017301_tRNA_H2)
WBGene00017301_tRNA_H2$Generation<-as.numeric(WBGene00017301_tRNA_H2$Generation)
rownames(WBGene00017301_tRNA_H2)<-c()

WBGene00017301_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00017301")
WBGene00017301_RNA_H1<-WBGene00017301_RNA_H1[,4:13]
WBGene00017301_RNA_H1<-t(WBGene00017301_RNA_H1)
generations<-rownames(WBGene00017301_RNA_H1)
WBGene00017301_RNA_H1<-cbind(generations,WBGene00017301_RNA_H1)
rownames(WBGene00017301_RNA_H1)<-c()
colnames(WBGene00017301_RNA_H1)<-c("Generation","Zscore")
WBGene00017301_RNA_H1<-as.data.frame(WBGene00017301_RNA_H1)

WBGene00017301_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00017301")
WBGene00017301_RNA_H2<-WBGene00017301_RNA_H2[,4:13]
WBGene00017301_RNA_H2<-t(WBGene00017301_RNA_H2)
generations<-rownames(WBGene00017301_RNA_H2)
WBGene00017301_RNA_H2<-cbind(generations,WBGene00017301_RNA_H2)
rownames(WBGene00017301_RNA_H2)<-c()
colnames(WBGene00017301_RNA_H2)<-c("Generation","Zscore")
WBGene00017301_RNA_H2<-as.data.frame(WBGene00017301_RNA_H2)
WBGene00017301_RNA_H2$Zscore[WBGene00017301_RNA_H2$Zscore == 'NA']<- 0.7343569

plot(WBGene00017301_tRNA_H1$Generation,WBGene00017301_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="hach-1 (WBGene00017301)/Met.CAT")
lines(WBGene00017301_tRNA_H2$Generation,WBGene00017301_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00017301_RNA_H1$Generation,WBGene00017301_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00017301_RNA_H2$Generation,WBGene00017301_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)



#gene 37
WBGene00005087_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[37,]
WBGene00005087_tRNA_H1<-WBGene00005087_tRNA_H1[,5:14]
WBGene00005087_tRNA_H1<-t(WBGene00005087_tRNA_H1)
generations<-rownames(WBGene00005087_tRNA_H1)
WBGene00005087_tRNA_H1<-cbind(generations,WBGene00005087_tRNA_H1)
colnames(WBGene00005087_tRNA_H1)<-c("Generation","Zscore")
WBGene00005087_tRNA_H1<-as.data.frame(WBGene00005087_tRNA_H1)
WBGene00005087_tRNA_H1$Generation<-as.numeric(WBGene00005087_tRNA_H1$Generation)
rownames(WBGene00005087_tRNA_H1)<-c()

WBGene00005087_tRNA_H2<-tRNA_z_scores_table_H2_WB[37,]
WBGene00005087_tRNA_H2<-WBGene00005087_tRNA_H2[,5:14]
WBGene00005087_tRNA_H2<-t(WBGene00005087_tRNA_H2)
generations<-rownames(WBGene00005087_tRNA_H2)
WBGene00005087_tRNA_H2<-cbind(generations,WBGene00005087_tRNA_H2)
colnames(WBGene00005087_tRNA_H2)<-c("Generation","Zscore")
WBGene00005087_tRNA_H2<-as.data.frame(WBGene00005087_tRNA_H2)
WBGene00005087_tRNA_H2$Generation<-as.numeric(WBGene00005087_tRNA_H2$Generation)
rownames(WBGene00005087_tRNA_H2)<-c()

WBGene00005087_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00005087")
WBGene00005087_RNA_H1<-WBGene00005087_RNA_H1[,4:13]
WBGene00005087_RNA_H1<-t(WBGene00005087_RNA_H1)
generations<-rownames(WBGene00005087_RNA_H1)
WBGene00005087_RNA_H1<-cbind(generations,WBGene00005087_RNA_H1)
rownames(WBGene00005087_RNA_H1)<-c()
colnames(WBGene00005087_RNA_H1)<-c("Generation","Zscore")
WBGene00005087_RNA_H1<-as.data.frame(WBGene00005087_RNA_H1)

WBGene00005087_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00005087")
WBGene00005087_RNA_H2<-WBGene00005087_RNA_H2[,4:13]
WBGene00005087_RNA_H2<-t(WBGene00005087_RNA_H2)
generations<-rownames(WBGene00005087_RNA_H2)
WBGene00005087_RNA_H2<-cbind(generations,WBGene00005087_RNA_H2)
rownames(WBGene00005087_RNA_H2)<-c()
colnames(WBGene00005087_RNA_H2)<-c("Generation","Zscore")
WBGene00005087_RNA_H2<-as.data.frame(WBGene00005087_RNA_H2)
WBGene00005087_RNA_H2$Zscore[WBGene00005087_RNA_H2$Zscore == 'NA']<- 1.433058

plot(WBGene00005087_tRNA_H1$Generation,WBGene00005087_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="hach-1 (WBGene00005087)/Pro.AGG")
lines(WBGene00005087_tRNA_H2$Generation,WBGene00005087_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00005087_RNA_H1$Generation,WBGene00005087_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00005087_RNA_H2$Generation,WBGene00005087_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)



#gene 38
WBGene00018900_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[38,]
WBGene00018900_tRNA_H1<-WBGene00018900_tRNA_H1[,5:14]
WBGene00018900_tRNA_H1<-t(WBGene00018900_tRNA_H1)
generations<-rownames(WBGene00018900_tRNA_H1)
WBGene00018900_tRNA_H1<-cbind(generations,WBGene00018900_tRNA_H1)
colnames(WBGene00018900_tRNA_H1)<-c("Generation","Zscore")
WBGene00018900_tRNA_H1<-as.data.frame(WBGene00018900_tRNA_H1)
WBGene00018900_tRNA_H1$Generation<-as.numeric(WBGene00018900_tRNA_H1$Generation)
rownames(WBGene00018900_tRNA_H1)<-c()

WBGene00018900_tRNA_H2<-tRNA_z_scores_table_H2_WB[38,]
WBGene00018900_tRNA_H2<-WBGene00018900_tRNA_H2[,5:14]
WBGene00018900_tRNA_H2<-t(WBGene00018900_tRNA_H2)
generations<-rownames(WBGene00018900_tRNA_H2)
WBGene00018900_tRNA_H2<-cbind(generations,WBGene00018900_tRNA_H2)
colnames(WBGene00018900_tRNA_H2)<-c("Generation","Zscore")
WBGene00018900_tRNA_H2<-as.data.frame(WBGene00018900_tRNA_H2)
WBGene00018900_tRNA_H2$Generation<-as.numeric(WBGene00018900_tRNA_H2$Generation)
rownames(WBGene00018900_tRNA_H2)<-c()

WBGene00018900_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00018900")
WBGene00018900_RNA_H1<-WBGene00018900_RNA_H1[,4:13]
WBGene00018900_RNA_H1<-t(WBGene00018900_RNA_H1)
generations<-rownames(WBGene00018900_RNA_H1)
WBGene00018900_RNA_H1<-cbind(generations,WBGene00018900_RNA_H1)
rownames(WBGene00018900_RNA_H1)<-c()
colnames(WBGene00018900_RNA_H1)<-c("Generation","Zscore")
WBGene00018900_RNA_H1<-as.data.frame(WBGene00018900_RNA_H1)

WBGene00018900_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00018900")
WBGene00018900_RNA_H2<-WBGene00018900_RNA_H2[,4:13]
WBGene00018900_RNA_H2<-t(WBGene00018900_RNA_H2)
generations<-rownames(WBGene00018900_RNA_H2)
WBGene00018900_RNA_H2<-cbind(generations,WBGene00018900_RNA_H2)
rownames(WBGene00018900_RNA_H2)<-c()
colnames(WBGene00018900_RNA_H2)<-c("Generation","Zscore")
WBGene00018900_RNA_H2<-as.data.frame(WBGene00018900_RNA_H2)
WBGene00018900_RNA_H2$Zscore[WBGene00018900_RNA_H2$Zscore == 'NA']<- 1.90804

plot(WBGene00018900_tRNA_H1$Generation,WBGene00018900_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="rod-1 (WBGene00018900)/Pro.AGG")
lines(WBGene00018900_tRNA_H2$Generation,WBGene00018900_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00018900_RNA_H1$Generation,WBGene00018900_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00018900_RNA_H2$Generation,WBGene00018900_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)



#gene 39
WBGene00005087_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[39,]
WBGene00005087_tRNA_H1<-WBGene00005087_tRNA_H1[,5:14]
WBGene00005087_tRNA_H1<-t(WBGene00005087_tRNA_H1)
generations<-rownames(WBGene00005087_tRNA_H1)
WBGene00005087_tRNA_H1<-cbind(generations,WBGene00005087_tRNA_H1)
colnames(WBGene00005087_tRNA_H1)<-c("Generation","Zscore")
WBGene00005087_tRNA_H1<-as.data.frame(WBGene00005087_tRNA_H1)
WBGene00005087_tRNA_H1$Generation<-as.numeric(WBGene00005087_tRNA_H1$Generation)
rownames(WBGene00005087_tRNA_H1)<-c()

WBGene00005087_tRNA_H2<-tRNA_z_scores_table_H2_WB[39,]
WBGene00005087_tRNA_H2<-WBGene00005087_tRNA_H2[,5:14]
WBGene00005087_tRNA_H2<-t(WBGene00005087_tRNA_H2)
generations<-rownames(WBGene00005087_tRNA_H2)
WBGene00005087_tRNA_H2<-cbind(generations,WBGene00005087_tRNA_H2)
colnames(WBGene00005087_tRNA_H2)<-c("Generation","Zscore")
WBGene00005087_tRNA_H2<-as.data.frame(WBGene00005087_tRNA_H2)
WBGene00005087_tRNA_H2$Generation<-as.numeric(WBGene00005087_tRNA_H2$Generation)
rownames(WBGene00005087_tRNA_H2)<-c()

WBGene00005087_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00005087")
WBGene00005087_RNA_H1<-WBGene00005087_RNA_H1[,4:13]
WBGene00005087_RNA_H1<-t(WBGene00005087_RNA_H1)
generations<-rownames(WBGene00005087_RNA_H1)
WBGene00005087_RNA_H1<-cbind(generations,WBGene00005087_RNA_H1)
rownames(WBGene00005087_RNA_H1)<-c()
colnames(WBGene00005087_RNA_H1)<-c("Generation","Zscore")
WBGene00005087_RNA_H1<-as.data.frame(WBGene00005087_RNA_H1)

WBGene00005087_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00005087")
WBGene00005087_RNA_H2<-WBGene00005087_RNA_H2[,4:13]
WBGene00005087_RNA_H2<-t(WBGene00005087_RNA_H2)
generations<-rownames(WBGene00005087_RNA_H2)
WBGene00005087_RNA_H2<-cbind(generations,WBGene00005087_RNA_H2)
rownames(WBGene00005087_RNA_H2)<-c()
colnames(WBGene00005087_RNA_H2)<-c("Generation","Zscore")
WBGene00005087_RNA_H2<-as.data.frame(WBGene00005087_RNA_H2)
WBGene00005087_RNA_H2$Zscore[WBGene00005087_RNA_H2$Zscore == 'NA']<- 1.90804

plot(WBGene00005087_tRNA_H1$Generation,WBGene00005087_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="srd-9 (WBGene00005087)/Pro.AGG")
lines(WBGene00005087_tRNA_H2$Generation,WBGene00005087_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00005087_RNA_H1$Generation,WBGene00005087_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00005087_RNA_H2$Generation,WBGene00005087_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)



#gene 40
WBGene00018900_bis_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[40,]
WBGene00018900_bis_tRNA_H1<-WBGene00018900_bis_tRNA_H1[,5:14]
WBGene00018900_bis_tRNA_H1<-t(WBGene00018900_bis_tRNA_H1)
generations<-rownames(WBGene00018900_bis_tRNA_H1)
WBGene00018900_bis_tRNA_H1<-cbind(generations,WBGene00018900_bis_tRNA_H1)
colnames(WBGene00018900_bis_tRNA_H1)<-c("Generation","Zscore")
WBGene00018900_bis_tRNA_H1<-as.data.frame(WBGene00018900_bis_tRNA_H1)
WBGene00018900_bis_tRNA_H1$Generation<-as.numeric(WBGene00018900_bis_tRNA_H1$Generation)
rownames(WBGene00018900_bis_tRNA_H1)<-c()

WBGene00018900_bis_tRNA_H2<-tRNA_z_scores_table_H2_WB[40,]
WBGene00018900_bis_tRNA_H2<-WBGene00018900_bis_tRNA_H2[,5:14]
WBGene00018900_bis_tRNA_H2<-t(WBGene00018900_bis_tRNA_H2)
generations<-rownames(WBGene00018900_bis_tRNA_H2)
WBGene00018900_bis_tRNA_H2<-cbind(generations,WBGene00018900_bis_tRNA_H2)
colnames(WBGene00018900_bis_tRNA_H2)<-c("Generation","Zscore")
WBGene00018900_bis_tRNA_H2<-as.data.frame(WBGene00018900_bis_tRNA_H2)
WBGene00018900_bis_tRNA_H2$Generation<-as.numeric(WBGene00018900_bis_tRNA_H2$Generation)
rownames(WBGene00018900_bis_tRNA_H2)<-c()

WBGene00018900_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00018900")
WBGene00018900_RNA_H1<-WBGene00018900_RNA_H1[,4:13]
WBGene00018900_RNA_H1<-t(WBGene00018900_RNA_H1)
generations<-rownames(WBGene00018900_RNA_H1)
WBGene00018900_RNA_H1<-cbind(generations,WBGene00018900_RNA_H1)
rownames(WBGene00018900_RNA_H1)<-c()
colnames(WBGene00018900_RNA_H1)<-c("Generation","Zscore")
WBGene00018900_RNA_H1<-as.data.frame(WBGene00018900_RNA_H1)

WBGene00018900_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00018900")
WBGene00018900_RNA_H2<-WBGene00018900_RNA_H2[,4:13]
WBGene00018900_RNA_H2<-t(WBGene00018900_RNA_H2)
generations<-rownames(WBGene00018900_RNA_H2)
WBGene00018900_RNA_H2<-cbind(generations,WBGene00018900_RNA_H2)
rownames(WBGene00018900_RNA_H2)<-c()
colnames(WBGene00018900_RNA_H2)<-c("Generation","Zscore")
WBGene00018900_RNA_H2<-as.data.frame(WBGene00018900_RNA_H2)
WBGene00018900_RNA_H2$Zscore[WBGene00018900_RNA_H2$Zscore == 'NA']<- 1.90804

plot(WBGene00018900_bis_tRNA_H1$Generation,WBGene00018900_bis_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="rod-1 (WBGene00018900)/Pro.AGG.2")
lines(WBGene00018900_bis_tRNA_H2$Generation,WBGene00018900_bis_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00018900_RNA_H1$Generation,WBGene00018900_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00018900_RNA_H2$Generation,WBGene00018900_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)



#gene 41
WBGene00012097_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[41,]
WBGene00012097_tRNA_H1<-WBGene00012097_tRNA_H1[,5:14]
WBGene00012097_tRNA_H1<-t(WBGene00012097_tRNA_H1)
generations<-rownames(WBGene00012097_tRNA_H1)
WBGene00012097_tRNA_H1<-cbind(generations,WBGene00012097_tRNA_H1)
colnames(WBGene00012097_tRNA_H1)<-c("Generation","Zscore")
WBGene00012097_tRNA_H1<-as.data.frame(WBGene00012097_tRNA_H1)
WBGene00012097_tRNA_H1$Generation<-as.numeric(WBGene00012097_tRNA_H1$Generation)
rownames(WBGene00012097_tRNA_H1)<-c()

WBGene00012097_tRNA_H2<-tRNA_z_scores_table_H2_WB[41,]
WBGene00012097_tRNA_H2<-WBGene00012097_tRNA_H2[,5:14]
WBGene00012097_tRNA_H2<-t(WBGene00012097_tRNA_H2)
generations<-rownames(WBGene00012097_tRNA_H2)
WBGene00012097_tRNA_H2<-cbind(generations,WBGene00012097_tRNA_H2)
colnames(WBGene00012097_tRNA_H2)<-c("Generation","Zscore")
WBGene00012097_tRNA_H2<-as.data.frame(WBGene00012097_tRNA_H2)
WBGene00012097_tRNA_H2$Generation<-as.numeric(WBGene00012097_tRNA_H2$Generation)
rownames(WBGene00012097_tRNA_H2)<-c()

WBGene00012097_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00012097")
WBGene00012097_RNA_H1<-WBGene00012097_RNA_H1[,4:13]
WBGene00012097_RNA_H1<-t(WBGene00012097_RNA_H1)
generations<-rownames(WBGene00012097_RNA_H1)
WBGene00012097_RNA_H1<-cbind(generations,WBGene00012097_RNA_H1)
rownames(WBGene00012097_RNA_H1)<-c()
colnames(WBGene00012097_RNA_H1)<-c("Generation","Zscore")
WBGene00012097_RNA_H1<-as.data.frame(WBGene00012097_RNA_H1)

WBGene00012097_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00012097")
WBGene00012097_RNA_H2<-WBGene00012097_RNA_H2[,4:13]
WBGene00012097_RNA_H2<-t(WBGene00012097_RNA_H2)
generations<-rownames(WBGene00012097_RNA_H2)
WBGene00012097_RNA_H2<-cbind(generations,WBGene00012097_RNA_H2)
rownames(WBGene00012097_RNA_H2)<-c()
colnames(WBGene00012097_RNA_H2)<-c("Generation","Zscore")
WBGene00012097_RNA_H2<-as.data.frame(WBGene00012097_RNA_H2)
WBGene00012097_RNA_H2$Zscore[WBGene00012097_RNA_H2$Zscore == 'NA']<- 0.2103809

plot(WBGene00012097_tRNA_H1$Generation,WBGene00012097_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="abcf-2 (WBGene00012097)/Pro.CGG")
lines(WBGene00012097_tRNA_H2$Generation,WBGene00012097_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00012097_RNA_H1$Generation,WBGene00012097_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00012097_RNA_H2$Generation,WBGene00012097_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)



#gene 42
WBGene00008906_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[42,]
WBGene00008906_tRNA_H1<-WBGene00008906_tRNA_H1[,5:14]
WBGene00008906_tRNA_H1<-t(WBGene00008906_tRNA_H1)
generations<-rownames(WBGene00008906_tRNA_H1)
WBGene00008906_tRNA_H1<-cbind(generations,WBGene00008906_tRNA_H1)
colnames(WBGene00008906_tRNA_H1)<-c("Generation","Zscore")
WBGene00008906_tRNA_H1<-as.data.frame(WBGene00008906_tRNA_H1)
WBGene00008906_tRNA_H1$Generation<-as.numeric(WBGene00008906_tRNA_H1$Generation)
rownames(WBGene00008906_tRNA_H1)<-c()

WBGene00008906_tRNA_H2<-tRNA_z_scores_table_H2_WB[42,]
WBGene00008906_tRNA_H2<-WBGene00008906_tRNA_H2[,5:14]
WBGene00008906_tRNA_H2<-t(WBGene00008906_tRNA_H2)
generations<-rownames(WBGene00008906_tRNA_H2)
WBGene00008906_tRNA_H2<-cbind(generations,WBGene00008906_tRNA_H2)
colnames(WBGene00008906_tRNA_H2)<-c("Generation","Zscore")
WBGene00008906_tRNA_H2<-as.data.frame(WBGene00008906_tRNA_H2)
WBGene00008906_tRNA_H2$Generation<-as.numeric(WBGene00008906_tRNA_H2$Generation)
rownames(WBGene00008906_tRNA_H2)<-c()

WBGene00008906_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00008906")
WBGene00008906_RNA_H1<-WBGene00008906_RNA_H1[,4:13]
WBGene00008906_RNA_H1<-t(WBGene00008906_RNA_H1)
generations<-rownames(WBGene00008906_RNA_H1)
WBGene00008906_RNA_H1<-cbind(generations,WBGene00008906_RNA_H1)
rownames(WBGene00008906_RNA_H1)<-c()
colnames(WBGene00008906_RNA_H1)<-c("Generation","Zscore")
WBGene00008906_RNA_H1<-as.data.frame(WBGene00008906_RNA_H1)

WBGene00008906_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00008906")
WBGene00008906_RNA_H2<-WBGene00008906_RNA_H2[,4:13]
WBGene00008906_RNA_H2<-t(WBGene00008906_RNA_H2)
generations<-rownames(WBGene00008906_RNA_H2)
WBGene00008906_RNA_H2<-cbind(generations,WBGene00008906_RNA_H2)
rownames(WBGene00008906_RNA_H2)<-c()
colnames(WBGene00008906_RNA_H2)<-c("Generation","Zscore")
WBGene00008906_RNA_H2<-as.data.frame(WBGene00008906_RNA_H2)
WBGene00008906_RNA_H2$Zscore[WBGene00008906_RNA_H2$Zscore == 'NA']<- 0.5653183

plot(WBGene00008906_tRNA_H1$Generation,WBGene00008906_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="oac-17 (WBGene00008906)/Pro.TGG")
lines(WBGene00008906_tRNA_H2$Generation,WBGene00008906_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00008906_RNA_H1$Generation,WBGene00008906_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00008906_RNA_H2$Generation,WBGene00008906_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)



#gene 43
WBGene00018900_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[43,]
WBGene00018900_tRNA_H1<-WBGene00018900_tRNA_H1[,5:14]
WBGene00018900_tRNA_H1<-t(WBGene00018900_tRNA_H1)
generations<-rownames(WBGene00018900_tRNA_H1)
WBGene00018900_tRNA_H1<-cbind(generations,WBGene00018900_tRNA_H1)
colnames(WBGene00018900_tRNA_H1)<-c("Generation","Zscore")
WBGene00018900_tRNA_H1<-as.data.frame(WBGene00018900_tRNA_H1)
WBGene00018900_tRNA_H1$Generation<-as.numeric(WBGene00018900_tRNA_H1$Generation)
rownames(WBGene00018900_tRNA_H1)<-c()

WBGene00018900_tRNA_H2<-tRNA_z_scores_table_H2_WB[43,]
WBGene00018900_tRNA_H2<-WBGene00018900_tRNA_H2[,5:14]
WBGene00018900_tRNA_H2<-t(WBGene00018900_tRNA_H2)
generations<-rownames(WBGene00018900_tRNA_H2)
WBGene00018900_tRNA_H2<-cbind(generations,WBGene00018900_tRNA_H2)
colnames(WBGene00018900_tRNA_H2)<-c("Generation","Zscore")
WBGene00018900_tRNA_H2<-as.data.frame(WBGene00018900_tRNA_H2)
WBGene00018900_tRNA_H2$Generation<-as.numeric(WBGene00018900_tRNA_H2$Generation)
rownames(WBGene00018900_tRNA_H2)<-c()

WBGene00018900_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00018900")
WBGene00018900_RNA_H1<-WBGene00018900_RNA_H1[,4:13]
WBGene00018900_RNA_H1<-t(WBGene00018900_RNA_H1)
generations<-rownames(WBGene00018900_RNA_H1)
WBGene00018900_RNA_H1<-cbind(generations,WBGene00018900_RNA_H1)
rownames(WBGene00018900_RNA_H1)<-c()
colnames(WBGene00018900_RNA_H1)<-c("Generation","Zscore")
WBGene00018900_RNA_H1<-as.data.frame(WBGene00018900_RNA_H1)

WBGene00018900_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00018900")
WBGene00018900_RNA_H2<-WBGene00018900_RNA_H2[,4:13]
WBGene00018900_RNA_H2<-t(WBGene00018900_RNA_H2)
generations<-rownames(WBGene00018900_RNA_H2)
WBGene00018900_RNA_H2<-cbind(generations,WBGene00018900_RNA_H2)
rownames(WBGene00018900_RNA_H2)<-c()
colnames(WBGene00018900_RNA_H2)<-c("Generation","Zscore")
WBGene00018900_RNA_H2<-as.data.frame(WBGene00018900_RNA_H2)
WBGene00018900_RNA_H2$Zscore[WBGene00018900_RNA_H2$Zscore == 'NA']<- 1.90804

plot(WBGene00018900_tRNA_H1$Generation,WBGene00018900_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="rod-1 (WBGene00018900)/Pro.TGG")
lines(WBGene00018900_tRNA_H2$Generation,WBGene00008906_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00018900_RNA_H1$Generation,WBGene00018900_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00018900_RNA_H2$Generation,WBGene00018900_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)



#gene 44
WBGene00000199_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[44,]
WBGene00000199_tRNA_H1<-WBGene00000199_tRNA_H1[,5:14]
WBGene00000199_tRNA_H1<-t(WBGene00000199_tRNA_H1)
generations<-rownames(WBGene00000199_tRNA_H1)
WBGene00000199_tRNA_H1<-cbind(generations,WBGene00000199_tRNA_H1)
colnames(WBGene00000199_tRNA_H1)<-c("Generation","Zscore")
WBGene00000199_tRNA_H1<-as.data.frame(WBGene00000199_tRNA_H1)
WBGene00000199_tRNA_H1$Generation<-as.numeric(WBGene00000199_tRNA_H1$Generation)
rownames(WBGene00000199_tRNA_H1)<-c()

WBGene00000199_tRNA_H2<-tRNA_z_scores_table_H2_WB[44,]
WBGene00000199_tRNA_H2<-WBGene00000199_tRNA_H2[,5:14]
WBGene00000199_tRNA_H2<-t(WBGene00000199_tRNA_H2)
generations<-rownames(WBGene00000199_tRNA_H2)
WBGene00000199_tRNA_H2<-cbind(generations,WBGene00000199_tRNA_H2)
colnames(WBGene00000199_tRNA_H2)<-c("Generation","Zscore")
WBGene00000199_tRNA_H2<-as.data.frame(WBGene00000199_tRNA_H2)
WBGene00000199_tRNA_H2$Generation<-as.numeric(WBGene00000199_tRNA_H2$Generation)
rownames(WBGene00000199_tRNA_H2)<-c()

WBGene00000199_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00000199")
WBGene00000199_RNA_H1<-WBGene00000199_RNA_H1[,4:13]
WBGene00000199_RNA_H1<-t(WBGene00000199_RNA_H1)
generations<-rownames(WBGene00000199_RNA_H1)
WBGene00000199_RNA_H1<-cbind(generations,WBGene00000199_RNA_H1)
rownames(WBGene00000199_RNA_H1)<-c()
colnames(WBGene00000199_RNA_H1)<-c("Generation","Zscore")
WBGene00000199_RNA_H1<-as.data.frame(WBGene00000199_RNA_H1)

WBGene00000199_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00000199")
WBGene00000199_RNA_H2<-WBGene00000199_RNA_H2[,4:13]
WBGene00000199_RNA_H2<-t(WBGene00000199_RNA_H2)
generations<-rownames(WBGene00000199_RNA_H2)
WBGene00000199_RNA_H2<-cbind(generations,WBGene00000199_RNA_H2)
rownames(WBGene00000199_RNA_H2)<-c()
colnames(WBGene00000199_RNA_H2)<-c("Generation","Zscore")
WBGene00000199_RNA_H2<-as.data.frame(WBGene00000199_RNA_H2)
WBGene00000199_RNA_H2$Zscore[WBGene00000199_RNA_H2$Zscore == 'NA']<- -0.02877572

plot(WBGene00000199_tRNA_H1$Generation,WBGene00000199_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="arx-1 (WBGene00000199)/Ser.CGA")
lines(WBGene00000199_tRNA_H2$Generation,WBGene00000199_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00000199_RNA_H1$Generation,WBGene00000199_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00000199_RNA_H2$Generation,WBGene00000199_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)



#gene 45
WBGene00012407_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[45,]
WBGene00012407_tRNA_H1<-WBGene00012407_tRNA_H1[,5:14]
WBGene00012407_tRNA_H1<-t(WBGene00012407_tRNA_H1)
generations<-rownames(WBGene00012407_tRNA_H1)
WBGene00012407_tRNA_H1<-cbind(generations,WBGene00012407_tRNA_H1)
colnames(WBGene00012407_tRNA_H1)<-c("Generation","Zscore")
WBGene00012407_tRNA_H1<-as.data.frame(WBGene00012407_tRNA_H1)
WBGene00012407_tRNA_H1$Generation<-as.numeric(WBGene00012407_tRNA_H1$Generation)
rownames(WBGene00012407_tRNA_H1)<-c()

WBGene00012407_tRNA_H2<-tRNA_z_scores_table_H2_WB[45,]
WBGene00012407_tRNA_H2<-WBGene00012407_tRNA_H2[,5:14]
WBGene00012407_tRNA_H2<-t(WBGene00012407_tRNA_H2)
generations<-rownames(WBGene00012407_tRNA_H2)
WBGene00012407_tRNA_H2<-cbind(generations,WBGene00012407_tRNA_H2)
colnames(WBGene00012407_tRNA_H2)<-c("Generation","Zscore")
WBGene00012407_tRNA_H2<-as.data.frame(WBGene00012407_tRNA_H2)
WBGene00012407_tRNA_H2$Generation<-as.numeric(WBGene00012407_tRNA_H2$Generation)
rownames(WBGene00012407_tRNA_H2)<-c()

WBGene00012407_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00012407")
WBGene00012407_RNA_H1<-WBGene00012407_RNA_H1[,4:13]
WBGene00012407_RNA_H1<-t(WBGene00012407_RNA_H1)
generations<-rownames(WBGene00012407_RNA_H1)
WBGene00012407_RNA_H1<-cbind(generations,WBGene00012407_RNA_H1)
rownames(WBGene00012407_RNA_H1)<-c()
colnames(WBGene00012407_RNA_H1)<-c("Generation","Zscore")
WBGene00012407_RNA_H1<-as.data.frame(WBGene00012407_RNA_H1)

WBGene00012407_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00012407")
WBGene00012407_RNA_H2<-WBGene00012407_RNA_H2[,4:13]
WBGene00012407_RNA_H2<-t(WBGene00012407_RNA_H2)
generations<-rownames(WBGene00012407_RNA_H2)
WBGene00012407_RNA_H2<-cbind(generations,WBGene00012407_RNA_H2)
rownames(WBGene00012407_RNA_H2)<-c()
colnames(WBGene00012407_RNA_H2)<-c("Generation","Zscore")
WBGene00012407_RNA_H2<-as.data.frame(WBGene00012407_RNA_H2)
WBGene00012407_RNA_H2$Zscore[WBGene00012407_RNA_H2$Zscore == 'NA']<- -0.3562984

plot(WBGene00012407_tRNA_H1$Generation,WBGene00012407_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="Y7A5A.1 (WBGene00012407)/Ser.CGA")
lines(WBGene00012407_tRNA_H2$Generation,WBGene00012407_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00012407_RNA_H1$Generation,WBGene00012407_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00012407_RNA_H2$Generation,WBGene00012407_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)



#gene 46
WBGene00018547_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[46,]
WBGene00018547_tRNA_H1<-WBGene00018547_tRNA_H1[,5:14]
WBGene00018547_tRNA_H1<-t(WBGene00018547_tRNA_H1)
generations<-rownames(WBGene00018547_tRNA_H1)
WBGene00018547_tRNA_H1<-cbind(generations,WBGene00018547_tRNA_H1)
colnames(WBGene00018547_tRNA_H1)<-c("Generation","Zscore")
WBGene00018547_tRNA_H1<-as.data.frame(WBGene00018547_tRNA_H1)
WBGene00018547_tRNA_H1$Generation<-as.numeric(WBGene00018547_tRNA_H1$Generation)
rownames(WBGene00018547_tRNA_H1)<-c()

WBGene00018547_tRNA_H2<-tRNA_z_scores_table_H2_WB[46,]
WBGene00018547_tRNA_H2<-WBGene00018547_tRNA_H2[,5:14]
WBGene00018547_tRNA_H2<-t(WBGene00018547_tRNA_H2)
generations<-rownames(WBGene00018547_tRNA_H2)
WBGene00018547_tRNA_H2<-cbind(generations,WBGene00018547_tRNA_H2)
colnames(WBGene00018547_tRNA_H2)<-c("Generation","Zscore")
WBGene00018547_tRNA_H2<-as.data.frame(WBGene00018547_tRNA_H2)
WBGene00018547_tRNA_H2$Generation<-as.numeric(WBGene00018547_tRNA_H2$Generation)
rownames(WBGene00018547_tRNA_H2)<-c()

WBGene00018547_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00018547")
WBGene00018547_RNA_H1<-WBGene00018547_RNA_H1[,4:13]
WBGene00018547_RNA_H1<-t(WBGene00018547_RNA_H1)
generations<-rownames(WBGene00018547_RNA_H1)
WBGene00018547_RNA_H1<-cbind(generations,WBGene00018547_RNA_H1)
rownames(WBGene00018547_RNA_H1)<-c()
colnames(WBGene00018547_RNA_H1)<-c("Generation","Zscore")
WBGene00018547_RNA_H1<-as.data.frame(WBGene00018547_RNA_H1)

WBGene00018547_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00018547")
WBGene00018547_RNA_H2<-WBGene00018547_RNA_H2[,4:13]
WBGene00018547_RNA_H2<-t(WBGene00018547_RNA_H2)
generations<-rownames(WBGene00018547_RNA_H2)
WBGene00018547_RNA_H2<-cbind(generations,WBGene00018547_RNA_H2)
rownames(WBGene00018547_RNA_H2)<-c()
colnames(WBGene00018547_RNA_H2)<-c("Generation","Zscore")
WBGene00018547_RNA_H2<-as.data.frame(WBGene00018547_RNA_H2)
WBGene00018547_RNA_H2$Zscore[WBGene00018547_RNA_H2$Zscore == 'NA']<- 1.811032

plot(WBGene00018547_tRNA_H1$Generation,WBGene00018547_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="clec-78 (WBGene00018547)/Trp.CCA")
lines(WBGene00018547_tRNA_H2$Generation,WBGene00018547_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00018547_RNA_H1$Generation,WBGene00018547_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00018547_RNA_H2$Generation,WBGene00018547_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)



#gene 47
WBGene00018547_bis_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[47,]
WBGene00018547_bis_tRNA_H1<-WBGene00018547_bis_tRNA_H1[,5:14]
WBGene00018547_bis_tRNA_H1<-t(WBGene00018547_bis_tRNA_H1)
generations<-rownames(WBGene00018547_bis_tRNA_H1)
WBGene00018547_bis_tRNA_H1<-cbind(generations,WBGene00018547_bis_tRNA_H1)
colnames(WBGene00018547_bis_tRNA_H1)<-c("Generation","Zscore")
WBGene00018547_bis_tRNA_H1<-as.data.frame(WBGene00018547_bis_tRNA_H1)
WBGene00018547_bis_tRNA_H1$Generation<-as.numeric(WBGene00018547_bis_tRNA_H1$Generation)
rownames(WBGene00018547_bis_tRNA_H1)<-c()

WBGene00018547_bis_tRNA_H2<-tRNA_z_scores_table_H2_WB[47,]
WBGene00018547_bis_tRNA_H2<-WBGene00018547_bis_tRNA_H2[,5:14]
WBGene00018547_bis_tRNA_H2<-t(WBGene00018547_bis_tRNA_H2)
generations<-rownames(WBGene00018547_bis_tRNA_H2)
WBGene00018547_bis_tRNA_H2<-cbind(generations,WBGene00018547_bis_tRNA_H2)
colnames(WBGene00018547_bis_tRNA_H2)<-c("Generation","Zscore")
WBGene00018547_bis_tRNA_H2<-as.data.frame(WBGene00018547_bis_tRNA_H2)
WBGene00018547_bis_tRNA_H2$Generation<-as.numeric(WBGene00018547_bis_tRNA_H2$Generation)
rownames(WBGene00018547_bis_tRNA_H2)<-c()

WBGene00018547_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00018547")
WBGene00018547_RNA_H1<-WBGene00018547_RNA_H1[,4:13]
WBGene00018547_RNA_H1<-t(WBGene00018547_RNA_H1)
generations<-rownames(WBGene00018547_RNA_H1)
WBGene00018547_RNA_H1<-cbind(generations,WBGene00018547_RNA_H1)
rownames(WBGene00018547_RNA_H1)<-c()
colnames(WBGene00018547_RNA_H1)<-c("Generation","Zscore")
WBGene00018547_RNA_H1<-as.data.frame(WBGene00018547_RNA_H1)

WBGene00018547_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00018547")
WBGene00018547_RNA_H2<-WBGene00018547_RNA_H2[,4:13]
WBGene00018547_RNA_H2<-t(WBGene00018547_RNA_H2)
generations<-rownames(WBGene00018547_RNA_H2)
WBGene00018547_RNA_H2<-cbind(generations,WBGene00018547_RNA_H2)
rownames(WBGene00018547_RNA_H2)<-c()
colnames(WBGene00018547_RNA_H2)<-c("Generation","Zscore")
WBGene00018547_RNA_H2<-as.data.frame(WBGene00018547_RNA_H2)
WBGene00018547_RNA_H2$Zscore[WBGene00018547_RNA_H2$Zscore == 'NA']<- 1.811032

plot(WBGene00018547_bis_tRNA_H1$Generation,WBGene00018547_bis_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="clec-78 (WBGene00018547)/Trp.CCA.2")
lines(WBGene00018547_bis_tRNA_H2$Generation,WBGene00018547_bis_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00018547_RNA_H1$Generation,WBGene00018547_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00018547_RNA_H2$Generation,WBGene00018547_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)



#gene 48
WBGene00009069_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[48,]
WBGene00009069_tRNA_H1<-WBGene00009069_tRNA_H1[,5:14]
WBGene00009069_tRNA_H1<-t(WBGene00009069_tRNA_H1)
generations<-rownames(WBGene00009069_tRNA_H1)
WBGene00009069_tRNA_H1<-cbind(generations,WBGene00009069_tRNA_H1)
colnames(WBGene00009069_tRNA_H1)<-c("Generation","Zscore")
WBGene00009069_tRNA_H1<-as.data.frame(WBGene00009069_tRNA_H1)
WBGene00009069_tRNA_H1$Generation<-as.numeric(WBGene00009069_tRNA_H1$Generation)
rownames(WBGene00009069_tRNA_H1)<-c()

WBGene00009069_tRNA_H2<-tRNA_z_scores_table_H2_WB[48,]
WBGene00009069_tRNA_H2<-WBGene00009069_tRNA_H2[,5:14]
WBGene00009069_tRNA_H2<-t(WBGene00009069_tRNA_H2)
generations<-rownames(WBGene00009069_tRNA_H2)
WBGene00009069_tRNA_H2<-cbind(generations,WBGene00009069_tRNA_H2)
colnames(WBGene00009069_tRNA_H2)<-c("Generation","Zscore")
WBGene00009069_tRNA_H2<-as.data.frame(WBGene00009069_tRNA_H2)
WBGene00009069_tRNA_H2$Generation<-as.numeric(WBGene00009069_tRNA_H2$Generation)
rownames(WBGene00009069_tRNA_H2)<-c()

WBGene00009069_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00009069")
WBGene00009069_RNA_H1<-WBGene00009069_RNA_H1[,4:13]
WBGene00009069_RNA_H1<-t(WBGene00009069_RNA_H1)
generations<-rownames(WBGene00009069_RNA_H1)
WBGene00009069_RNA_H1<-cbind(generations,WBGene00009069_RNA_H1)
rownames(WBGene00009069_RNA_H1)<-c()
colnames(WBGene00009069_RNA_H1)<-c("Generation","Zscore")
WBGene00009069_RNA_H1<-as.data.frame(WBGene00009069_RNA_H1)

WBGene00009069_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00009069")
WBGene00009069_RNA_H2<-WBGene00009069_RNA_H2[,4:13]
WBGene00009069_RNA_H2<-t(WBGene00009069_RNA_H2)
generations<-rownames(WBGene00009069_RNA_H2)
WBGene00009069_RNA_H2<-cbind(generations,WBGene00009069_RNA_H2)
rownames(WBGene00009069_RNA_H2)<-c()
colnames(WBGene00009069_RNA_H2)<-c("Generation","Zscore")
WBGene00009069_RNA_H2<-as.data.frame(WBGene00009069_RNA_H2)
WBGene00009069_RNA_H2$Zscore[WBGene00009069_RNA_H2$Zscore == 'NA']<- -0.5111332

plot(WBGene00009069_tRNA_H1$Generation,WBGene00009069_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="F23A7.4 (WBGene00009069)/Tyr.GTA")
lines(WBGene00009069_tRNA_H2$Generation,WBGene00009069_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00009069_RNA_H1$Generation,WBGene00009069_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00009069_RNA_H2$Generation,WBGene00009069_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)



#gene 49
WBGene00044638_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[49,]
WBGene00044638_tRNA_H1<-WBGene00044638_tRNA_H1[,5:14]
WBGene00044638_tRNA_H1<-t(WBGene00044638_tRNA_H1)
generations<-rownames(WBGene00044638_tRNA_H1)
WBGene00044638_tRNA_H1<-cbind(generations,WBGene00044638_tRNA_H1)
colnames(WBGene00044638_tRNA_H1)<-c("Generation","Zscore")
WBGene00044638_tRNA_H1<-as.data.frame(WBGene00044638_tRNA_H1)
WBGene00044638_tRNA_H1$Generation<-as.numeric(WBGene00044638_tRNA_H1$Generation)
rownames(WBGene00044638_tRNA_H1)<-c()

WBGene00044638_tRNA_H2<-tRNA_z_scores_table_H2_WB[49,]
WBGene00044638_tRNA_H2<-WBGene00044638_tRNA_H2[,5:14]
WBGene00044638_tRNA_H2<-t(WBGene00044638_tRNA_H2)
generations<-rownames(WBGene00044638_tRNA_H2)
WBGene00044638_tRNA_H2<-cbind(generations,WBGene00044638_tRNA_H2)
colnames(WBGene00044638_tRNA_H2)<-c("Generation","Zscore")
WBGene00044638_tRNA_H2<-as.data.frame(WBGene00044638_tRNA_H2)
WBGene00044638_tRNA_H2$Generation<-as.numeric(WBGene00044638_tRNA_H2$Generation)
rownames(WBGene00044638_tRNA_H2)<-c()

WBGene00044638_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00044638")
WBGene00044638_RNA_H1<-WBGene00044638_RNA_H1[,4:13]
WBGene00044638_RNA_H1<-t(WBGene00044638_RNA_H1)
generations<-rownames(WBGene00044638_RNA_H1)
WBGene00044638_RNA_H1<-cbind(generations,WBGene00044638_RNA_H1)
rownames(WBGene00044638_RNA_H1)<-c()
colnames(WBGene00044638_RNA_H1)<-c("Generation","Zscore")
WBGene00044638_RNA_H1<-as.data.frame(WBGene00044638_RNA_H1)

WBGene00044638_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00044638")
WBGene00044638_RNA_H2<-WBGene00044638_RNA_H2[,4:13]
WBGene00044638_RNA_H2<-t(WBGene00044638_RNA_H2)
generations<-rownames(WBGene00044638_RNA_H2)
WBGene00044638_RNA_H2<-cbind(generations,WBGene00044638_RNA_H2)
rownames(WBGene00044638_RNA_H2)<-c()
colnames(WBGene00044638_RNA_H2)<-c("Generation","Zscore")
WBGene00044638_RNA_H2<-as.data.frame(WBGene00044638_RNA_H2)
WBGene00044638_RNA_H2$Zscore[WBGene00044638_RNA_H2$Zscore == 'NA']<- -0.2009013

plot(WBGene00044638_tRNA_H1$Generation,WBGene00044638_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="F23A7.8 (WBGene00044638)/Tyr.GTA")
lines(WBGene00044638_tRNA_H2$Generation,WBGene00044638_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00044638_RNA_H1$Generation,WBGene00044638_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00044638_RNA_H2$Generation,WBGene00044638_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)



#gene 50
WBGene00009069_bis_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[50,]
WBGene00009069_bis_tRNA_H1<-WBGene00009069_bis_tRNA_H1[,5:14]
WBGene00009069_bis_tRNA_H1<-t(WBGene00009069_bis_tRNA_H1)
generations<-rownames(WBGene00009069_bis_tRNA_H1)
WBGene00009069_bis_tRNA_H1<-cbind(generations,WBGene00009069_bis_tRNA_H1)
colnames(WBGene00009069_bis_tRNA_H1)<-c("Generation","Zscore")
WBGene00009069_bis_tRNA_H1<-as.data.frame(WBGene00009069_bis_tRNA_H1)
WBGene00009069_bis_tRNA_H1$Generation<-as.numeric(WBGene00009069_bis_tRNA_H1$Generation)
rownames(WBGene00009069_bis_tRNA_H1)<-c()

WBGene00009069_bis_tRNA_H2<-tRNA_z_scores_table_H2_WB[50,]
WBGene00009069_bis_tRNA_H2<-WBGene00009069_bis_tRNA_H2[,5:14]
WBGene00009069_bis_tRNA_H2<-t(WBGene00009069_bis_tRNA_H2)
generations<-rownames(WBGene00009069_bis_tRNA_H2)
WBGene00009069_bis_tRNA_H2<-cbind(generations,WBGene00009069_bis_tRNA_H2)
colnames(WBGene00009069_bis_tRNA_H2)<-c("Generation","Zscore")
WBGene00009069_bis_tRNA_H2<-as.data.frame(WBGene00009069_bis_tRNA_H2)
WBGene00009069_bis_tRNA_H2$Generation<-as.numeric(WBGene00009069_bis_tRNA_H2$Generation)
rownames(WBGene00009069_bis_tRNA_H2)<-c()

WBGene00009069_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00009069")
WBGene00009069_RNA_H1<-WBGene00009069_RNA_H1[,4:13]
WBGene00009069_RNA_H1<-t(WBGene00009069_RNA_H1)
generations<-rownames(WBGene00009069_RNA_H1)
WBGene00009069_RNA_H1<-cbind(generations,WBGene00009069_RNA_H1)
rownames(WBGene00009069_RNA_H1)<-c()
colnames(WBGene00009069_RNA_H1)<-c("Generation","Zscore")
WBGene00009069_RNA_H1<-as.data.frame(WBGene00009069_RNA_H1)

WBGene00009069_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00009069")
WBGene00009069_RNA_H2<-WBGene00009069_RNA_H2[,4:13]
WBGene00009069_RNA_H2<-t(WBGene00009069_RNA_H2)
generations<-rownames(WBGene00009069_RNA_H2)
WBGene00009069_RNA_H2<-cbind(generations,WBGene00009069_RNA_H2)
rownames(WBGene00009069_RNA_H2)<-c()
colnames(WBGene00009069_RNA_H2)<-c("Generation","Zscore")
WBGene00009069_RNA_H2<-as.data.frame(WBGene00009069_RNA_H2)
WBGene00009069_RNA_H2$Zscore[WBGene00009069_RNA_H2$Zscore == 'NA']<- -0.5111332

plot(WBGene00009069_bis_tRNA_H1$Generation,WBGene00009069_bis_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="F23A7.4 (WBGene00009069)/Tyr.GTA.2")
lines(WBGene00009069_bis_tRNA_H2$Generation,WBGene00009069_bis_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00009069_RNA_H1$Generation,WBGene00009069_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00009069_RNA_H2$Generation,WBGene00009069_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)



#gene 51
WBGene00044638_bis_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[51,]
WBGene00044638_bis_tRNA_H1<-WBGene00044638_bis_tRNA_H1[,5:14]
WBGene00044638_bis_tRNA_H1<-t(WBGene00044638_bis_tRNA_H1)
generations<-rownames(WBGene00044638_bis_tRNA_H1)
WBGene00044638_bis_tRNA_H1<-cbind(generations,WBGene00044638_bis_tRNA_H1)
colnames(WBGene00044638_bis_tRNA_H1)<-c("Generation","Zscore")
WBGene00044638_bis_tRNA_H1<-as.data.frame(WBGene00044638_bis_tRNA_H1)
WBGene00044638_bis_tRNA_H1$Generation<-as.numeric(WBGene00044638_bis_tRNA_H1$Generation)
rownames(WBGene00044638_bis_tRNA_H1)<-c()

WBGene00044638_bis_tRNA_H2<-tRNA_z_scores_table_H2_WB[51,]
WBGene00044638_bis_tRNA_H2<-WBGene00044638_bis_tRNA_H2[,5:14]
WBGene00044638_bis_tRNA_H2<-t(WBGene00044638_bis_tRNA_H2)
generations<-rownames(WBGene00044638_bis_tRNA_H2)
WBGene00044638_bis_tRNA_H2<-cbind(generations,WBGene00044638_bis_tRNA_H2)
colnames(WBGene00044638_bis_tRNA_H2)<-c("Generation","Zscore")
WBGene00044638_bis_tRNA_H2<-as.data.frame(WBGene00044638_bis_tRNA_H2)
WBGene00044638_bis_tRNA_H2$Generation<-as.numeric(WBGene00044638_bis_tRNA_H2$Generation)
rownames(WBGene00044638_bis_tRNA_H2)<-c()

WBGene00044638_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00044638")
WBGene00044638_RNA_H1<-WBGene00044638_RNA_H1[,4:13]
WBGene00044638_RNA_H1<-t(WBGene00044638_RNA_H1)
generations<-rownames(WBGene00044638_RNA_H1)
WBGene00044638_RNA_H1<-cbind(generations,WBGene00044638_RNA_H1)
rownames(WBGene00044638_RNA_H1)<-c()
colnames(WBGene00044638_RNA_H1)<-c("Generation","Zscore")
WBGene00044638_RNA_H1<-as.data.frame(WBGene00044638_RNA_H1)

WBGene00044638_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00044638")
WBGene00044638_RNA_H2<-WBGene00044638_RNA_H2[,4:13]
WBGene00044638_RNA_H2<-t(WBGene00044638_RNA_H2)
generations<-rownames(WBGene00044638_RNA_H2)
WBGene00044638_RNA_H2<-cbind(generations,WBGene00044638_RNA_H2)
rownames(WBGene00044638_RNA_H2)<-c()
colnames(WBGene00044638_RNA_H2)<-c("Generation","Zscore")
WBGene00044638_RNA_H2<-as.data.frame(WBGene00044638_RNA_H2)
WBGene00044638_RNA_H2$Zscore[WBGene00044638_RNA_H2$Zscore == 'NA']<- -0.2009013

plot(WBGene00044638_bis_tRNA_H1$Generation,WBGene00044638_bis_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="F23A7.8 (WBGene00044638)/Tyr.GTA.2")
lines(WBGene00044638_bis_tRNA_H2$Generation,WBGene00044638_bis_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00044638_RNA_H1$Generation,WBGene00044638_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00044638_RNA_H2$Generation,WBGene00044638_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)



#gene 52
WBGene00010923_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[52,]
WBGene00010923_tRNA_H1<-WBGene00010923_tRNA_H1[,5:14]
WBGene00010923_tRNA_H1<-t(WBGene00010923_tRNA_H1)
generations<-rownames(WBGene00010923_tRNA_H1)
WBGene00010923_tRNA_H1<-cbind(generations,WBGene00010923_tRNA_H1)
colnames(WBGene00010923_tRNA_H1)<-c("Generation","Zscore")
WBGene00010923_tRNA_H1<-as.data.frame(WBGene00010923_tRNA_H1)
WBGene00010923_tRNA_H1$Generation<-as.numeric(WBGene00010923_tRNA_H1$Generation)
rownames(WBGene00010923_tRNA_H1)<-c()

WBGene00010923_tRNA_H2<-tRNA_z_scores_table_H2_WB[52,]
WBGene00010923_tRNA_H2<-WBGene00010923_tRNA_H2[,5:14]
WBGene00010923_tRNA_H2<-t(WBGene00010923_tRNA_H2)
generations<-rownames(WBGene00010923_tRNA_H2)
WBGene00010923_tRNA_H2<-cbind(generations,WBGene00010923_tRNA_H2)
colnames(WBGene00010923_tRNA_H2)<-c("Generation","Zscore")
WBGene00010923_tRNA_H2<-as.data.frame(WBGene00010923_tRNA_H2)
WBGene00010923_tRNA_H2$Generation<-as.numeric(WBGene00010923_tRNA_H2$Generation)
rownames(WBGene00010923_tRNA_H2)<-c()

WBGene00010923_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00010923")
WBGene00010923_RNA_H1<-WBGene00010923_RNA_H1[,4:13]
WBGene00010923_RNA_H1<-t(WBGene00010923_RNA_H1)
generations<-rownames(WBGene00010923_RNA_H1)
WBGene00010923_RNA_H1<-cbind(generations,WBGene00010923_RNA_H1)
rownames(WBGene00010923_RNA_H1)<-c()
colnames(WBGene00010923_RNA_H1)<-c("Generation","Zscore")
WBGene00010923_RNA_H1<-as.data.frame(WBGene00010923_RNA_H1)

WBGene00010923_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00010923")
WBGene00010923_RNA_H2<-WBGene00010923_RNA_H2[,4:13]
WBGene00010923_RNA_H2<-t(WBGene00010923_RNA_H2)
generations<-rownames(WBGene00010923_RNA_H2)
WBGene00010923_RNA_H2<-cbind(generations,WBGene00010923_RNA_H2)
rownames(WBGene00010923_RNA_H2)<-c()
colnames(WBGene00010923_RNA_H2)<-c("Generation","Zscore")
WBGene00010923_RNA_H2<-as.data.frame(WBGene00010923_RNA_H2)
WBGene00010923_RNA_H2$Zscore[WBGene00010923_RNA_H2$Zscore == 'NA']<- 1.163906

plot(WBGene00010923_tRNA_H1$Generation,WBGene00010923_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="rle-1 (WBGene00010923)/Val.AAC")
lines(WBGene00010923_tRNA_H2$Generation,WBGene00010923_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00010923_RNA_H1$Generation,WBGene00010923_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00010923_RNA_H2$Generation,WBGene00010923_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)



#gene 53
WBGene00003411_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[53,]
WBGene00003411_tRNA_H1<-WBGene00003411_tRNA_H1[,5:14]
WBGene00003411_tRNA_H1<-t(WBGene00003411_tRNA_H1)
generations<-rownames(WBGene00003411_tRNA_H1)
WBGene00003411_tRNA_H1<-cbind(generations,WBGene00003411_tRNA_H1)
colnames(WBGene00003411_tRNA_H1)<-c("Generation","Zscore")
WBGene00003411_tRNA_H1<-as.data.frame(WBGene00003411_tRNA_H1)
WBGene00003411_tRNA_H1$Generation<-as.numeric(WBGene00003411_tRNA_H1$Generation)
rownames(WBGene00003411_tRNA_H1)<-c()

WBGene00003411_tRNA_H2<-tRNA_z_scores_table_H2_WB[53,]
WBGene00003411_tRNA_H2<-WBGene00003411_tRNA_H2[,5:14]
WBGene00003411_tRNA_H2<-t(WBGene00003411_tRNA_H2)
generations<-rownames(WBGene00003411_tRNA_H2)
WBGene00003411_tRNA_H2<-cbind(generations,WBGene00003411_tRNA_H2)
colnames(WBGene00003411_tRNA_H2)<-c("Generation","Zscore")
WBGene00003411_tRNA_H2<-as.data.frame(WBGene00003411_tRNA_H2)
WBGene00003411_tRNA_H2$Generation<-as.numeric(WBGene00003411_tRNA_H2$Generation)
rownames(WBGene00003411_tRNA_H2)<-c()

WBGene00003411_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00003411")
WBGene00003411_RNA_H1<-WBGene00003411_RNA_H1[,4:13]
WBGene00003411_RNA_H1<-t(WBGene00003411_RNA_H1)
generations<-rownames(WBGene00003411_RNA_H1)
WBGene00003411_RNA_H1<-cbind(generations,WBGene00003411_RNA_H1)
rownames(WBGene00003411_RNA_H1)<-c()
colnames(WBGene00003411_RNA_H1)<-c("Generation","Zscore")
WBGene00003411_RNA_H1<-as.data.frame(WBGene00003411_RNA_H1)

WBGene00003411_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00003411")
WBGene00003411_RNA_H2<-WBGene00003411_RNA_H2[,4:13]
WBGene00003411_RNA_H2<-t(WBGene00003411_RNA_H2)
generations<-rownames(WBGene00003411_RNA_H2)
WBGene00003411_RNA_H2<-cbind(generations,WBGene00003411_RNA_H2)
rownames(WBGene00003411_RNA_H2)<-c()
colnames(WBGene00003411_RNA_H2)<-c("Generation","Zscore")
WBGene00003411_RNA_H2<-as.data.frame(WBGene00003411_RNA_H2)
WBGene00003411_RNA_H2$Zscore[WBGene00003411_RNA_H2$Zscore == 'NA']<- 1.097042

plot(WBGene00003411_tRNA_H1$Generation,WBGene00003411_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="mrp-5 (WBGene00003411)/Val.AAC")
lines(WBGene00003411_tRNA_H2$Generation,WBGene00003411_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00003411_RNA_H1$Generation,WBGene00003411_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00003411_RNA_H2$Generation,WBGene00003411_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)



#gene 54
WBGene00000818_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[54,]
WBGene00000818_tRNA_H1<-WBGene00000818_tRNA_H1[,5:14]
WBGene00000818_tRNA_H1<-t(WBGene00000818_tRNA_H1)
generations<-rownames(WBGene00000818_tRNA_H1)
WBGene00000818_tRNA_H1<-cbind(generations,WBGene00000818_tRNA_H1)
colnames(WBGene00000818_tRNA_H1)<-c("Generation","Zscore")
WBGene00000818_tRNA_H1<-as.data.frame(WBGene00000818_tRNA_H1)
WBGene00000818_tRNA_H1$Generation<-as.numeric(WBGene00000818_tRNA_H1$Generation)
rownames(WBGene00000818_tRNA_H1)<-c()

WBGene00000818_tRNA_H2<-tRNA_z_scores_table_H2_WB[54,]
WBGene00000818_tRNA_H2<-WBGene00000818_tRNA_H2[,5:14]
WBGene00000818_tRNA_H2<-t(WBGene00000818_tRNA_H2)
generations<-rownames(WBGene00000818_tRNA_H2)
WBGene00000818_tRNA_H2<-cbind(generations,WBGene00000818_tRNA_H2)
colnames(WBGene00000818_tRNA_H2)<-c("Generation","Zscore")
WBGene00000818_tRNA_H2<-as.data.frame(WBGene00000818_tRNA_H2)
WBGene00000818_tRNA_H2$Generation<-as.numeric(WBGene00000818_tRNA_H2$Generation)
rownames(WBGene00000818_tRNA_H2)<-c()

WBGene00000818_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00000818")
WBGene00000818_RNA_H1<-WBGene00000818_RNA_H1[,4:13]
WBGene00000818_RNA_H1<-t(WBGene00000818_RNA_H1)
generations<-rownames(WBGene00000818_RNA_H1)
WBGene00000818_RNA_H1<-cbind(generations,WBGene00000818_RNA_H1)
rownames(WBGene00000818_RNA_H1)<-c()
colnames(WBGene00000818_RNA_H1)<-c("Generation","Zscore")
WBGene00000818_RNA_H1<-as.data.frame(WBGene00000818_RNA_H1)

WBGene00000818_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00000818")
WBGene00000818_RNA_H2<-WBGene00000818_RNA_H2[,4:13]
WBGene00000818_RNA_H2<-t(WBGene00000818_RNA_H2)
generations<-rownames(WBGene00000818_RNA_H2)
WBGene00000818_RNA_H2<-cbind(generations,WBGene00000818_RNA_H2)
rownames(WBGene00000818_RNA_H2)<-c()
colnames(WBGene00000818_RNA_H2)<-c("Generation","Zscore")
WBGene00000818_RNA_H2<-as.data.frame(WBGene00000818_RNA_H2)
WBGene00000818_RNA_H2$Zscore[WBGene00000818_RNA_H2$Zscore == 'NA']<- 0.3818171

plot(WBGene00000818_tRNA_H1$Generation,WBGene00000818_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="csn-6 (WBGene00000818)/Val.TAC")
lines(WBGene00000818_tRNA_H2$Generation,WBGene00000818_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00000818_RNA_H1$Generation,WBGene00000818_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00000818_RNA_H2$Generation,WBGene00000818_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)



#gene 55
WBGene00021562_tRNA_H1<-tRNA_z_scores_table_H1bis_WB[55,]
WBGene00021562_tRNA_H1<-WBGene00021562_tRNA_H1[,5:14]
WBGene00021562_tRNA_H1<-t(WBGene00021562_tRNA_H1)
generations<-rownames(WBGene00021562_tRNA_H1)
WBGene00021562_tRNA_H1<-cbind(generations,WBGene00021562_tRNA_H1)
colnames(WBGene00021562_tRNA_H1)<-c("Generation","Zscore")
WBGene00021562_tRNA_H1<-as.data.frame(WBGene00021562_tRNA_H1)
WBGene00021562_tRNA_H1$Generation<-as.numeric(WBGene00021562_tRNA_H1$Generation)
rownames(WBGene00021562_tRNA_H1)<-c()

WBGene00021562_tRNA_H2<-tRNA_z_scores_table_H2_WB[55,]
WBGene00021562_tRNA_H2<-WBGene00021562_tRNA_H2[,5:14]
WBGene00021562_tRNA_H2<-t(WBGene00021562_tRNA_H2)
generations<-rownames(WBGene00021562_tRNA_H2)
WBGene00021562_tRNA_H2<-cbind(generations,WBGene00021562_tRNA_H2)
colnames(WBGene00021562_tRNA_H2)<-c("Generation","Zscore")
WBGene00021562_tRNA_H2<-as.data.frame(WBGene00021562_tRNA_H2)
WBGene00021562_tRNA_H2$Generation<-as.numeric(WBGene00021562_tRNA_H2$Generation)
rownames(WBGene00021562_tRNA_H2)<-c()

WBGene00021562_RNA_H1<-RNA_z_scores_table_H1 %>% filter(RNA_z_scores_table_H1$genes %in% "WBGene00021562")
WBGene00021562_RNA_H1<-WBGene00021562_RNA_H1[,4:13]
WBGene00021562_RNA_H1<-t(WBGene00021562_RNA_H1)
generations<-rownames(WBGene00021562_RNA_H1)
WBGene00021562_RNA_H1<-cbind(generations,WBGene00021562_RNA_H1)
rownames(WBGene00021562_RNA_H1)<-c()
colnames(WBGene00021562_RNA_H1)<-c("Generation","Zscore")
WBGene00021562_RNA_H1<-as.data.frame(WBGene00021562_RNA_H1)

WBGene00021562_RNA_H2<-RNA_z_scores_table_H2 %>% filter(RNA_z_scores_table_H2$genes %in% "WBGene00021562")
WBGene00021562_RNA_H2<-WBGene00021562_RNA_H2[,4:13]
WBGene00021562_RNA_H2<-t(WBGene00021562_RNA_H2)
generations<-rownames(WBGene00021562_RNA_H2)
WBGene00021562_RNA_H2<-cbind(generations,WBGene00021562_RNA_H2)
rownames(WBGene00021562_RNA_H2)<-c()
colnames(WBGene00021562_RNA_H2)<-c("Generation","Zscore")
WBGene00021562_RNA_H2<-as.data.frame(WBGene00021562_RNA_H2)
WBGene00021562_RNA_H2$Zscore[WBGene00021562_RNA_H2$Zscore == 'NA']<- 0.1673972

plot(WBGene00021562_tRNA_H1$Generation,WBGene00021562_tRNA_H1$Zscore, col="blue",type="l",ylim=c(-5,5),xlab="Zscore",ylab="Generation",main="nuo-5 (WBGene00021562)/Val.TAC")
lines(WBGene00021562_tRNA_H2$Generation,WBGene00021562_tRNA_H2$Zscore, col="blue",type="l",lty=3)
lines(WBGene00021562_RNA_H1$Generation,WBGene00021562_RNA_H1$Zscore, col="red",type="l")
lines(WBGene00021562_RNA_H2$Generation,WBGene00021562_RNA_H2$Zscore, col="red",type="l",lty=3)
abline(h=2.5,lty=2)
abline(h=-2.5,lty=2)
legend(17, 5, legend=c("RNA", "tRNA","H1","H2"),
       col=c("red", "blue","black","black"), lty=c(1,1,1,3), cex=0.8)














#-------------
#######Sup.Fig.5 - Work on  genes involved in apoptosis####
apo<-read.xlsx("apoptosis_genes.xlsx")
apo<-apo$gene

#Gene expression epimutations on apoptotic genes
load("RNA_Control_integrated_table.Rdata")

RNA_Control_integrated_table_C1<-subset(RNA_Control_integrated_table,RNA_Control_integrated_table$Lineage=="C1")
RNA_Control_integrated_table_C2<-subset(RNA_Control_integrated_table,RNA_Control_integrated_table$Lineage=="C2")

C1_RNA_epi_apo <- filter(RNA_Control_integrated_table_C1, gene %in% apo)
C2_RNA_epi_apo <- filter(RNA_Control_integrated_table_C2, gene %in% apo)

load("RNA_Low_dose_integrated_table.Rdata")

RNA_Low_dose_integrated_table_L1<-subset(RNA_Low_dose_integrated_table,RNA_Low_dose_integrated_table$Lineage=="L1")
RNA_Low_dose_integrated_table_L2<-subset(RNA_Low_dose_integrated_table,RNA_Low_dose_integrated_table$Lineage=="L2")

L1_RNA_epi_apo <- filter(RNA_Low_dose_integrated_table_L1, gene %in% apo)
L2_RNA_epi_apo <- filter(RNA_Low_dose_integrated_table_L2, gene %in% apo)

load("RNA_High_dose_integrated_table.Rdata")

RNA_High_dose_integrated_table_H1<-subset(RNA_High_dose_integrated_table,RNA_High_dose_integrated_table$Lineage=="H1")
RNA_High_dose_integrated_table_H2<-subset(RNA_High_dose_integrated_table,RNA_High_dose_integrated_table$Lineage=="H2")

H1_RNA_epi_apo <- filter(RNA_High_dose_integrated_table_H1, gene %in% apo)
H2_RNA_epi_apo <- filter(RNA_High_dose_integrated_table_H2, gene %in% apo)

C1_RNA_epi_apo<-subset(C1_RNA_epi_apo,C1_RNA_epi_apo$is_RNA_exp_change==1)
C2_RNA_epi_apo<-subset(C2_RNA_epi_apo,C2_RNA_epi_apo$is_RNA_exp_change==1)
L1_RNA_epi_apo<-subset(L1_RNA_epi_apo,L1_RNA_epi_apo$is_RNA_exp_change==1)
L2_RNA_epi_apo<-subset(L2_RNA_epi_apo,L2_RNA_epi_apo$is_RNA_exp_change==1)
H1_RNA_epi_apo<-subset(H1_RNA_epi_apo,H1_RNA_epi_apo$is_RNA_exp_change==1)
H2_RNA_epi_apo<-subset(H2_RNA_epi_apo,H2_RNA_epi_apo$is_RNA_exp_change==1)

All_apo_epimut<-rbind(C1_RNA_epi_apo,C2_RNA_epi_apo,L1_RNA_epi_apo,L2_RNA_epi_apo,H1_RNA_epi_apo,H2_RNA_epi_apo)
write.xlsx(All_apo_epimut,"All_apo_epimut.xlsx")

ap<-read.xlsx("epimut_apo.xlsx")

#Total nb of epimutations
ap$Condition <- fct_relevel(ap$Condition, c("Control","Low dose","High dose"))
ap_plot <- ggplot(ap, aes(x=Condition  , y=Nbr.RNA.epimut, color = Condition)) + 
  geom_boxplot(fatten = 1, lwd = 1, width=0.5)+
  scale_color_manual(values=c("cornflowerblue", "darkgreen", "red"))+
  labs(y = "Number of new epimutations", x = "\nGeneration")+
  geom_jitter(aes(shape=Lineage, size=8),
              position=position_jitter(width = 0.2,
                                       height = 0.2))+
  theme_bw()+
  theme(axis.text.x = element_text(color = "black", size = 20, angle = 0, face = "plain"),
        axis.text.y = element_text(color = "black", size = 20, angle = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 16, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 16, angle = 90, face = "bold"))+
  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))+  
  scale_fill_brewer(palette="Pastel2")+
  ggtitle(paste(""))
ap_plot

ggdensity(ap$Nbr.RNA.epimut, 
          main = "Density plot of epimutations_number",
          xlab = "Total number of epimutations_number")
shapiro.test(ap$Nbr.RNA.epimut)
kruskal.test(Nbr.RNA.epimut ~ Condition, data = ap)
dunnTest(Nbr.RNA.epimut ~ Condition, data = ap)

#####Epimutation rate
load("RNA_C1epimutations.Rdata")
load("RNA_C2epimutations.Rdata")
load("RNA_L1epimutations.Rdata")
load("RNA_L2epimutations.Rdata")
load("RNA_H1epimutations.Rdata")
load("RNA_H2epimutations.Rdata")

####Select only genes involved in apoptosis
genes_names <- RNA_C1epimutations$genes %>%
  strsplit( ":" ) %>%
  sapply( "[", 4 )
RNA_C1epimutations_apo <- filter(RNA_C1epimutations, genes_names %in% apo)
RNA_C2epimutations_apo <- filter(RNA_C2epimutations, genes_names %in% apo)
RNA_L1epimutations_apo <- filter(RNA_L1epimutations, genes_names %in% apo)
RNA_L2epimutations_apo <- filter(RNA_L2epimutations, genes_names %in% apo)
RNA_H1epimutations_apo <- filter(RNA_H1epimutations, genes_names %in% apo)
RNA_H2epimutations_apo <- filter(RNA_H2epimutations, genes_names %in% apo)

#Calculate the new epimutations each generation
#Function creation
# Defining the number of transitions UP/DOWN per generation for each data type
# This code counts the number of transitions in chromatin/gene expression state in each generation

# UP transition function
#define function
UP_transition_func<-function(vector_in,input_name){
  Gen_ON <- vector()
  for(i in 2:length(vector_in)){
    if(vector_in[i]==1&vector_in[i-1]== 0){
      Gen_ON[[i]] <-1 # there is an UP transition to ON at generation[i]
    }
    
    if(vector_in[i]== 1&vector_in[i-1]== -1){
      Gen_ON[[i]] <-1 # there is an UP transition to ON at generation[i]
    }
    
    if(vector_in[i]== 1&vector_in[i-1]=="NA"){
      Gen_ON[[i]] <-1 # there is an UP transition to ON at generation[i]
    }
    
    if(vector_in[i]== -1&vector_in[i-1]== 0){
      Gen_ON[[i]] <- 0 # this is a DOWN transition
    }
    
    if(vector_in[i]== -1&vector_in[i-1]== 1){
      Gen_ON[[i]] <- 0 # this is a DOWN epimutation
    }
    
    if(vector_in[i]== -1&vector_in[i-1]== "NA"){
      Gen_ON[[i]] <- 0 # this is a DOWN epimutation
    }
    
    if(vector_in[i]==1&vector_in[i-1]==1){
      Gen_ON[[i]] <- 0 # inside a run of 1 1 epimutations at generation[i]
    }
    
    if(vector_in[i]==-1&vector_in[i-1]==-1){
      Gen_ON[[i]] <- 0 # inside a run of -1 -1 epimutations at generation[i]
    }
    
    if(vector_in[i]==0&vector_in[i-1]==0){
      Gen_ON[[i]] <- 0 # between epimutations at generation[i]
    }
    
    if(vector_in[i]==0&vector_in[i-1]==1){
      Gen_ON[[i]] <- 0 # there is a transition to OFF at generation[i]
    }
    
    if(vector_in[i]==0&vector_in[i-1]==-1){
      Gen_ON[[i]] <- 0 # there is a transition to OFF at generation[i]
    }
    
    if(vector_in[i]==0&vector_in[i-1]=="NA"){
      Gen_ON[[i]] <- 0 # there is a transition to OFF at generation[i]
    } 
    output <- Gen_ON
  }
  
  return(output)
}

# DOWN transition function

#define function
DOWN_transition_func<-function(vector_in,input_name){
  Gen_ON <- vector()
  for(i in 2:length(vector_in)){
    if(vector_in[i]== -1&vector_in[i-1]== 0){
      Gen_ON[[i]] <-1 # there is a DOWN transition to ON at generation[i]
    }
    
    if(vector_in[i]== -1&vector_in[i-1]== 1){
      Gen_ON[[i]] <-1 # there is a DOWN transition to ON at generation[i]
    }
    
    if(vector_in[i]== -1&vector_in[i-1]== "NA"){
      Gen_ON[[i]] <-1 # there is a DOWN transition to ON at generation[i]
    }
    
    if(vector_in[i]== 1&vector_in[i-1]== 0){
      Gen_ON[[i]] <- 0 # this is an UP transition
    }
    
    if(vector_in[i]== 1&vector_in[i-1]== -1){
      Gen_ON[[i]] <- 0 # this is an UP epimutation
    }
    
    if(vector_in[i]==1&vector_in[i-1]==1){
      Gen_ON[[i]] <- 0 # inside a run of 1 1 epimutations at generation[i]
    }
    
    if(vector_in[i]==1&vector_in[i-1]=="NA"){
      Gen_ON[[i]] <- 0 # this is an UP epimutation
    }
    
    if(vector_in[i]==-1&vector_in[i-1]==-1){
      Gen_ON[[i]] <- 0 # inside a run of -1 -1 epimutations at generation[i]
    }
    
    if(vector_in[i]==0&vector_in[i-1]==0){
      Gen_ON[[i]] <- 0 # between epimutations at generation[i]
    }
    
    if(vector_in[i]==0&vector_in[i-1]==1){
      Gen_ON[[i]] <- 0 # there is a transition to OFF at generation[i]
    }
    
    if(vector_in[i]==0&vector_in[i-1]==-1){
      Gen_ON[[i]] <- 0 # there is a transition to OFF at generation[i]
    }
    
    if(vector_in[i]==0&vector_in[i-1]=="NA"){
      Gen_ON[[i]] <- 0 # there is a transition to OFF at generation[i]
    }
    
    output <- Gen_ON
  }
  
  return(output)
}
#Control1 new epi
#UP
UP_output_RNA_C1<- c()
row.names(RNA_C1epimutations_apo)<-RNA_C1epimutations_apo$genes
RNA_C1epimutations_apo<-RNA_C1epimutations_apo[,c(3:13)]
for(i in 1:nrow(RNA_C1epimutations_apo)){
  UP_output_RNA_C1 <-rbind(UP_output_RNA_C1, UP_transition_func(RNA_C1epimutations_apo[i,], input_name=row.names(RNA_C1epimutations_apo)[i]))}
colnames(UP_output_RNA_C1) <- c("0", "2", "4", "6", "8", "10", "12", "14", "16", "18", "20")
row.names(UP_output_RNA_C1) <- rownames(RNA_C1epimutations_apo)
# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 
transitions_at_2 <- sum(UP_output_RNA_C1[,2])
transitions_at_4 <- sum(UP_output_RNA_C1[,3])
transitions_at_6 <- sum(UP_output_RNA_C1[,4])
transitions_at_8 <- sum(UP_output_RNA_C1[,5])
transitions_at_10 <- sum(UP_output_RNA_C1[,6])
transitions_at_12 <- sum(UP_output_RNA_C1[,7])
transitions_at_14 <- sum(UP_output_RNA_C1[,8])
transitions_at_16 <- sum(UP_output_RNA_C1[,9])
transitions_at_18 <- sum(UP_output_RNA_C1[,10])
transitions_at_20 <- sum(UP_output_RNA_C1[,11])

UP_RNA_C1_Table_of_new_epimutations <- rbind(transitions_at_2,
                                             transitions_at_4, 
                                             transitions_at_6,
                                             transitions_at_8, 
                                             transitions_at_10, 
                                             transitions_at_12, 
                                             transitions_at_14, 
                                             transitions_at_16, 
                                             transitions_at_18, 
                                             transitions_at_20)

#DOWN
DOWN_output_RNA_C1<- c()
for(i in 1:nrow(RNA_C1epimutations_apo)){
  DOWN_output_RNA_C1 <-rbind(DOWN_output_RNA_C1, DOWN_transition_func(RNA_C1epimutations_apo[i,], input_name=row.names(RNA_C1epimutations_apo)[i]))}
colnames(DOWN_output_RNA_C1) <- c("0", "2", "4", "6", "8", "10", "12", "14", "16", "18", "20")
row.names(DOWN_output_RNA_C1) <- rownames(RNA_C1epimutations_apo)
# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 
transitions_at_2 <- sum(DOWN_output_RNA_C1[,2])
transitions_at_4 <- sum(DOWN_output_RNA_C1[,3])
transitions_at_6 <- sum(DOWN_output_RNA_C1[,4])
transitions_at_8 <- sum(DOWN_output_RNA_C1[,5])
transitions_at_10 <- sum(DOWN_output_RNA_C1[,6])
transitions_at_12 <- sum(DOWN_output_RNA_C1[,7])
transitions_at_14 <- sum(DOWN_output_RNA_C1[,8])
transitions_at_16 <- sum(DOWN_output_RNA_C1[,9])
transitions_at_18 <- sum(DOWN_output_RNA_C1[,10])
transitions_at_20 <- sum(DOWN_output_RNA_C1[,11])
DOWN_RNA_C1_Table_of_new_epimutations <- rbind(transitions_at_2,
                                               transitions_at_4, 
                                               transitions_at_6,
                                               transitions_at_8, 
                                               transitions_at_10, 
                                               transitions_at_12, 
                                               transitions_at_14, 
                                               transitions_at_16, 
                                               transitions_at_18, 
                                               transitions_at_20)

#Control2 new epi
#UP
UP_output_RNA_C2<- c()
row.names(RNA_C2epimutations_apo)<-RNA_C2epimutations_apo$genes
RNA_C2epimutations_apo<-RNA_C2epimutations_apo[,c(3:13)]
for(i in 1:nrow(RNA_C2epimutations_apo)){
  UP_output_RNA_C2 <-rbind(UP_output_RNA_C2, UP_transition_func(RNA_C2epimutations_apo[i,], input_name=row.names(RNA_C2epimutations_apo)[i]))}
colnames(UP_output_RNA_C2) <- c("0", "2", "4", "6", "8", "10", "12","14", "16", "18", "20")
row.names(UP_output_RNA_C2) <- rownames(RNA_C2epimutations_apo)
# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 
transitions_at_2 <- sum(UP_output_RNA_C2[,2])
transitions_at_4 <- sum(UP_output_RNA_C2[,3])
transitions_at_6 <- sum(UP_output_RNA_C2[,4])
transitions_at_8 <- sum(UP_output_RNA_C2[,5])
transitions_at_10 <- sum(UP_output_RNA_C2[,6])
transitions_at_12 <- NA
transitions_at_14 <- sum(UP_output_RNA_C2[,8])
transitions_at_16 <- sum(UP_output_RNA_C2[,9])
transitions_at_18 <- sum(UP_output_RNA_C2[,10])
transitions_at_20 <- sum(UP_output_RNA_C2[,11])
UP_RNA_C2_Table_of_new_epimutations <- rbind(transitions_at_2,
                                             transitions_at_4, 
                                             transitions_at_6,
                                             transitions_at_8, 
                                             transitions_at_10, 
                                             transitions_at_12, 
                                             transitions_at_14, 
                                             transitions_at_16, 
                                             transitions_at_18, 
                                             transitions_at_20)

#DOWN
DOWN_output_RNA_C2<- c()
for(i in 1:nrow(RNA_C2epimutations_apo)){
  DOWN_output_RNA_C2 <-rbind(DOWN_output_RNA_C2, DOWN_transition_func(RNA_C2epimutations_apo[i,], input_name=row.names(RNA_C2epimutations_apo)[i]))}
colnames(DOWN_output_RNA_C2) <- c("0", "2", "4", "6", "8", "10", "12", "14", "16", "18", "20")
row.names(DOWN_output_RNA_C2) <- rownames(RNA_C2epimutations_apo)
# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 
transitions_at_2 <- sum(DOWN_output_RNA_C2[,2])
transitions_at_4 <- sum(DOWN_output_RNA_C2[,3])
transitions_at_6 <- sum(DOWN_output_RNA_C2[,4])
transitions_at_8 <- sum(DOWN_output_RNA_C2[,5])
transitions_at_10 <- sum(DOWN_output_RNA_C2[,6])
transitions_at_12 <- NA
transitions_at_14 <- sum(DOWN_output_RNA_C2[,8])
transitions_at_16 <- sum(DOWN_output_RNA_C2[,9])
transitions_at_18 <- sum(DOWN_output_RNA_C2[,10])
transitions_at_20 <- sum(DOWN_output_RNA_C2[,11])
DOWN_RNA_C2_Table_of_new_epimutations <- rbind(transitions_at_2,
                                               transitions_at_4, 
                                               transitions_at_6,
                                               transitions_at_8, 
                                               transitions_at_10, 
                                               transitions_at_12, 
                                               transitions_at_14, 
                                               transitions_at_16, 
                                               transitions_at_18, 
                                               transitions_at_20)

#Low1 new epi
#UP
UP_output_RNA_L1<- c()
row.names(RNA_L1epimutations_apo)<-RNA_L1epimutations_apo$genes
RNA_L1epimutations_apo<-RNA_L1epimutations_apo[,c(3:13)]
for(i in 1:nrow(RNA_L1epimutations_apo)){
  UP_output_RNA_L1 <-rbind(UP_output_RNA_L1, UP_transition_func(RNA_L1epimutations_apo[i,], input_name=row.names(RNA_L1epimutations_apo)[i]))}
colnames(UP_output_RNA_L1) <- c("0", "2", "4", "6", "8", "10", "12","14", "16", "18", "20")
row.names(UP_output_RNA_L1) <- rownames(RNA_L1epimutations_apo)
# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 
transitions_at_2 <- sum(UP_output_RNA_L1[,2])
transitions_at_4 <- sum(UP_output_RNA_L1[,3])
transitions_at_6 <- sum(UP_output_RNA_L1[,4])
transitions_at_8 <- sum(UP_output_RNA_L1[,5])
transitions_at_10 <- sum(UP_output_RNA_L1[,6])
transitions_at_12 <- sum(UP_output_RNA_L1[,7])
transitions_at_14 <- sum(UP_output_RNA_L1[,8])
transitions_at_16 <- sum(UP_output_RNA_L1[,9])
transitions_at_18 <- sum(UP_output_RNA_L1[,10])
transitions_at_20 <- sum(UP_output_RNA_L1[,11])
UP_RNA_L1_Table_of_new_epimutations <- rbind(transitions_at_2,
                                             transitions_at_4, 
                                             transitions_at_6,
                                             transitions_at_8, 
                                             transitions_at_10, 
                                             transitions_at_12, 
                                             transitions_at_14, 
                                             transitions_at_16, 
                                             transitions_at_18, 
                                             transitions_at_20)

#DOWN
DOWN_output_RNA_L1<- c()
for(i in 1:nrow(RNA_L1epimutations_apo)){
  DOWN_output_RNA_L1 <-rbind(DOWN_output_RNA_L1, DOWN_transition_func(RNA_L1epimutations_apo[i,], input_name=row.names(RNA_L1epimutations_apo)[i]))}
colnames(DOWN_output_RNA_L1) <- c("0", "2", "4", "6", "8", "10", "12", "14", "16", "18", "20")
row.names(DOWN_output_RNA_L1) <- rownames(RNA_L1epimutations_apo)
# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 
transitions_at_2 <- sum(DOWN_output_RNA_L1[,2])
transitions_at_4 <- sum(DOWN_output_RNA_L1[,3])
transitions_at_6 <- sum(DOWN_output_RNA_L1[,4])
transitions_at_8 <- sum(DOWN_output_RNA_L1[,5])
transitions_at_10 <- sum(DOWN_output_RNA_L1[,6])
transitions_at_12 <- sum(DOWN_output_RNA_L1[,7])
transitions_at_14 <- sum(DOWN_output_RNA_L1[,8])
transitions_at_16 <- sum(DOWN_output_RNA_L1[,9])
transitions_at_18 <- sum(DOWN_output_RNA_L1[,10])
transitions_at_20 <- sum(DOWN_output_RNA_L1[,11])
DOWN_RNA_L1_Table_of_new_epimutations <- rbind(transitions_at_2,
                                               transitions_at_4, 
                                               transitions_at_6,
                                               transitions_at_8, 
                                               transitions_at_10, 
                                               transitions_at_12, 
                                               transitions_at_14, 
                                               transitions_at_16, 
                                               transitions_at_18, 
                                               transitions_at_20)

#Low2 new epi
#UP
UP_output_RNA_L2<- c()
row.names(RNA_L2epimutations_apo)<-RNA_L2epimutations_apo$genes
RNA_L2epimutations_apo<-RNA_L2epimutations_apo[,c(3:13)]
for(i in 1:nrow(RNA_L2epimutations_apo)){
  UP_output_RNA_L2 <-rbind(UP_output_RNA_L2, UP_transition_func(RNA_L2epimutations_apo[i,], input_name=row.names(RNA_L2epimutations_apo)[i]))}
colnames(UP_output_RNA_L2) <- c("0", "2", "4", "6", "8", "10", "12","14", "16", "18", "20")
row.names(UP_output_RNA_L2) <- rownames(RNA_L2epimutations_apo)
# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 
transitions_at_2 <- sum(UP_output_RNA_L2[,2])
transitions_at_4 <- sum(UP_output_RNA_L2[,3])
transitions_at_6 <- sum(UP_output_RNA_L2[,4])
transitions_at_8 <- sum(UP_output_RNA_L2[,5])
transitions_at_10 <- sum(UP_output_RNA_L2[,6])
transitions_at_12 <- sum(UP_output_RNA_L2[,7])
transitions_at_14 <- sum(UP_output_RNA_L2[,8])
transitions_at_16 <- sum(UP_output_RNA_L2[,9])
transitions_at_18 <- sum(UP_output_RNA_L2[,10])
transitions_at_20 <- sum(UP_output_RNA_L2[,11])
UP_RNA_L2_Table_of_new_epimutations <- rbind(transitions_at_2,
                                             transitions_at_4, 
                                             transitions_at_6,
                                             transitions_at_8, 
                                             transitions_at_10, 
                                             transitions_at_12, 
                                             transitions_at_14, 
                                             transitions_at_16, 
                                             transitions_at_18, 
                                             transitions_at_20)

#DOWN
DOWN_output_RNA_L2<- c()
for(i in 1:nrow(RNA_L2epimutations_apo)){
  DOWN_output_RNA_L2 <-rbind(DOWN_output_RNA_L2, DOWN_transition_func(RNA_L2epimutations_apo[i,], input_name=row.names(RNA_L2epimutations_apo)[i]))}
colnames(DOWN_output_RNA_L2) <- c("0", "2", "4", "6", "8", "10", "12", "14", "16", "18", "20")
row.names(DOWN_output_RNA_L2) <- rownames(RNA_L2epimutations_apo)
# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 
transitions_at_2 <- sum(DOWN_output_RNA_L2[,2])
transitions_at_4 <- sum(DOWN_output_RNA_L2[,3])
transitions_at_6 <- sum(DOWN_output_RNA_L2[,4])
transitions_at_8 <- sum(DOWN_output_RNA_L2[,5])
transitions_at_10 <- sum(DOWN_output_RNA_L2[,6])
transitions_at_12 <- sum(DOWN_output_RNA_L2[,7])
transitions_at_14 <- sum(DOWN_output_RNA_L2[,8])
transitions_at_16 <- sum(DOWN_output_RNA_L2[,9])
transitions_at_18 <- sum(DOWN_output_RNA_L2[,10])
transitions_at_20 <- sum(DOWN_output_RNA_L2[,11])
DOWN_RNA_L2_Table_of_new_epimutations <- rbind(transitions_at_2,
                                               transitions_at_4, 
                                               transitions_at_6,
                                               transitions_at_8, 
                                               transitions_at_10, 
                                               transitions_at_12, 
                                               transitions_at_14, 
                                               transitions_at_16, 
                                               transitions_at_18, 
                                               transitions_at_20)

#High1 new epi
#UP
UP_output_RNA_H1<- c()
row.names(RNA_H1epimutations_apo)<-RNA_H1epimutations_apo$genes
RNA_H1epimutations_apo<-RNA_H1epimutations_apo[,c(3:13)]
for(i in 1:nrow(RNA_H1epimutations_apo)){
  UP_output_RNA_H1 <-rbind(UP_output_RNA_H1, UP_transition_func(RNA_H1epimutations_apo[i,], input_name=row.names(RNA_H1epimutations_apo)[i]))}
colnames(UP_output_RNA_H1) <- c("0", "2", "4", "6", "8", "10", "12","14", "16", "18")
row.names(UP_output_RNA_H1) <- rownames(RNA_H1epimutations_apo)
# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 
transitions_at_2 <- sum(UP_output_RNA_H1[,2])
transitions_at_4 <- sum(UP_output_RNA_H1[,3])
transitions_at_6 <- sum(UP_output_RNA_H1[,4])
transitions_at_8 <- sum(UP_output_RNA_H1[,5])
transitions_at_10 <- sum(UP_output_RNA_H1[,6])
transitions_at_12 <- sum(UP_output_RNA_H1[,7])
transitions_at_14 <- sum(UP_output_RNA_H1[,8])
transitions_at_16 <- sum(UP_output_RNA_H1[,9])
transitions_at_18 <- sum(UP_output_RNA_H1[,10])
UP_RNA_H1_Table_of_new_epimutations <- rbind(transitions_at_2,
                                             transitions_at_4, 
                                             transitions_at_6,
                                             transitions_at_8, 
                                             transitions_at_10, 
                                             transitions_at_12, 
                                             transitions_at_14, 
                                             transitions_at_16, 
                                             transitions_at_18)

#DOWN
DOWN_output_RNA_H1<- c()
for(i in 1:nrow(RNA_H1epimutations_apo)){
  DOWN_output_RNA_H1 <-rbind(DOWN_output_RNA_H1, DOWN_transition_func(RNA_H1epimutations_apo[i,], input_name=row.names(RNA_H1epimutations_apo)[i]))}
colnames(DOWN_output_RNA_H1) <- c("0", "2", "4", "6", "8", "10", "12", "14", "16", "18")
row.names(DOWN_output_RNA_H1) <- rownames(RNA_H1epimutations_apo)
# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 
transitions_at_2 <- sum(DOWN_output_RNA_H1[,2])
transitions_at_4 <- sum(DOWN_output_RNA_H1[,3])
transitions_at_6 <- sum(DOWN_output_RNA_H1[,4])
transitions_at_8 <- sum(DOWN_output_RNA_H1[,5])
transitions_at_10 <- sum(DOWN_output_RNA_H1[,6])
transitions_at_12 <- sum(DOWN_output_RNA_H1[,7])
transitions_at_14 <- sum(DOWN_output_RNA_H1[,8])
transitions_at_16 <- sum(DOWN_output_RNA_H1[,9])
transitions_at_18 <- sum(DOWN_output_RNA_H1[,10])
DOWN_RNA_H1_Table_of_new_epimutations <- rbind(transitions_at_2,
                                               transitions_at_4, 
                                               transitions_at_6,
                                               transitions_at_8, 
                                               transitions_at_10, 
                                               transitions_at_12, 
                                               transitions_at_14, 
                                               transitions_at_16, 
                                               transitions_at_18)
#High2 new epi
#UP
UP_output_RNA_H2<- c()
row.names(RNA_H2epimutations_apo)<-RNA_H2epimutations_apo$genes
RNA_H2epimutations_apo<-RNA_H2epimutations_apo[,c(3:13)]
for(i in 1:nrow(RNA_H2epimutations_apo)){
  UP_output_RNA_H2 <-rbind(UP_output_RNA_H2, UP_transition_func(RNA_H2epimutations_apo[i,], input_name=row.names(RNA_H2epimutations_apo)[i]))}
colnames(UP_output_RNA_H2) <- c("0", "2", "4", "6", "8", "10", "12","14", "16", "18","20")
row.names(UP_output_RNA_H2) <- rownames(RNA_H2epimutations_apo)
# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 
transitions_at_2 <- sum(UP_output_RNA_H2[,2])
transitions_at_4 <- sum(UP_output_RNA_H2[,3])
transitions_at_6 <- sum(UP_output_RNA_H2[,4])
transitions_at_8 <- sum(UP_output_RNA_H2[,5])
transitions_at_10 <- sum(UP_output_RNA_H2[,6])
transitions_at_12 <- sum(UP_output_RNA_H2[,7])
transitions_at_14 <- sum(UP_output_RNA_H2[,8])
transitions_at_16 <- sum(UP_output_RNA_H2[,9])
transitions_at_18 <- sum(UP_output_RNA_H2[,10])
transitions_at_20 <- sum(UP_output_RNA_H2[,11])
UP_RNA_H2_Table_of_new_epimutations <- rbind(transitions_at_2,
                                             transitions_at_4, 
                                             transitions_at_6,
                                             transitions_at_8, 
                                             transitions_at_10, 
                                             transitions_at_12, 
                                             transitions_at_14, 
                                             transitions_at_16, 
                                             transitions_at_18,
                                             transitions_at_20)

#DOWN
DOWN_output_RNA_H2<- c()
for(i in 1:nrow(RNA_H2epimutations_apo)){
  DOWN_output_RNA_H2 <-rbind(DOWN_output_RNA_H2, DOWN_transition_func(RNA_H2epimutations_apo[i,], input_name=row.names(RNA_H2epimutations_apo)[i]))}
colnames(DOWN_output_RNA_H2) <- c("0", "2", "4", "6", "8", "10", "12", "14", "16", "18","20")
row.names(DOWN_output_RNA_H2) <- rownames(RNA_H2epimutations_apo)
# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 
transitions_at_2 <- sum(DOWN_output_RNA_H2[,2])
transitions_at_4 <- sum(DOWN_output_RNA_H2[,3])
transitions_at_6 <- sum(DOWN_output_RNA_H2[,4])
transitions_at_8 <- sum(DOWN_output_RNA_H2[,5])
transitions_at_10 <- sum(DOWN_output_RNA_H2[,6])
transitions_at_12 <- sum(DOWN_output_RNA_H2[,7])
transitions_at_14 <- sum(DOWN_output_RNA_H2[,8])
transitions_at_16 <- sum(DOWN_output_RNA_H2[,9])
transitions_at_18 <- sum(DOWN_output_RNA_H2[,10])
transitions_at_20 <- sum(DOWN_output_RNA_H2[,11])
DOWN_RNA_H2_Table_of_new_epimutations <- rbind(transitions_at_2,
                                               transitions_at_4, 
                                               transitions_at_6,
                                               transitions_at_8, 
                                               transitions_at_10, 
                                               transitions_at_12, 
                                               transitions_at_14, 
                                               transitions_at_16, 
                                               transitions_at_18,
                                               transitions_at_20)
#Allnewepi
Transitions<-row.names(UP_RNA_C1_Table_of_new_epimutations)
colnames(UP_RNA_C1_Table_of_new_epimutations)<-c("UP")
colnames(DOWN_RNA_C1_Table_of_new_epimutations)<-c("DOWN")
RNA_apo_epiC1 <- cbind(Transitions,UP_RNA_C1_Table_of_new_epimutations,DOWN_RNA_C1_Table_of_new_epimutations)
RNA_apo_epiC1<-as.data.frame(RNA_apo_epiC1)
RNA_apo_epiC1 <- RNA_apo_epiC1 %>%
  # Creating an empty column:
  add_column(Lineage = "C1", .before="UP")
colnames(UP_RNA_C2_Table_of_new_epimutations)<-c("UP")
colnames(DOWN_RNA_C2_Table_of_new_epimutations)<-c("DOWN")
RNA_apo_epiC2 <- cbind(Transitions,UP_RNA_C2_Table_of_new_epimutations,DOWN_RNA_C2_Table_of_new_epimutations)
RNA_apo_epiC2<-as.data.frame(RNA_apo_epiC2)
RNA_apo_epiC2 <- RNA_apo_epiC2 %>%
  # Creating an empty column:
  add_column(Lineage = "C2", .before="UP")
RNA_apo_epiC <- rbind(RNA_apo_epiC1, RNA_apo_epiC2)
RNA_apo_epiC <- RNA_apo_epi_C %>%
  # Creating an empty column:
  add_column(Condition = "Control", .before="Lineage")

colnames(UP_RNA_L1_Table_of_new_epimutations)<-c("UP")
colnames(DOWN_RNA_L1_Table_of_new_epimutations)<-c("DOWN")
RNA_apo_epiL1 <- cbind(Transitions,UP_RNA_L1_Table_of_new_epimutations,DOWN_RNA_L1_Table_of_new_epimutations)
RNA_apo_epiL1<-as.data.frame(RNA_apo_epiL1)
RNA_apo_epiL1 <- RNA_apo_epiL1 %>%
  # Creating an empty column:
  add_column(Lineage = "L1", .before="UP")
colnames(UP_RNA_L2_Table_of_new_epimutations)<-c("UP")
colnames(DOWN_RNA_L2_Table_of_new_epimutations)<-c("DOWN")
RNA_apo_epiL2 <- cbind(Transitions,UP_RNA_L2_Table_of_new_epimutations,DOWN_RNA_L2_Table_of_new_epimutations)
RNA_apo_epiL2<-as.data.frame(RNA_apo_epiL2)
RNA_apo_epiL2 <- RNA_apo_epiL2 %>%
  # Creating an empty column:
  add_column(Lineage = "L2", .before="UP")
RNA_apo_epiL <- rbind(RNA_apo_epiL1, RNA_apo_epiL2)
RNA_apo_epiL <- RNA_apo_epiL %>%
  # Creating an empty column:
  add_column(Condition = "Low dose", .before="Lineage")

colnames(UP_RNA_H1_Table_of_new_epimutations)<-c("UP")
colnames(DOWN_RNA_H1_Table_of_new_epimutations)<-c("DOWN")
RNA_apo_epiH1 <- cbind(Transitions,UP_RNA_H1_Table_of_new_epimutations,DOWN_RNA_H1_Table_of_new_epimutations)
RNA_apo_epiH1<-as.data.frame(RNA_apo_epiH1)
RNA_apo_epiH1 <- RNA_apo_epiH1 %>%
  # Creating an empty column:
  add_column(Lineage = "H1", .before="UP")
colnames(UP_RNA_H2_Table_of_new_epimutations)<-c("UP")
colnames(DOWN_RNA_H2_Table_of_new_epimutations)<-c("DOWN")
RNA_apo_epiH2 <- cbind(Transitions,UP_RNA_H2_Table_of_new_epimutations,DOWN_RNA_H2_Table_of_new_epimutations)
RNA_apo_epiH2<-as.data.frame(RNA_apo_epiH2)
RNA_apo_epiH2 <- RNA_apo_epiH2 %>%
  # Creating an empty column:
  add_column(Lineage = "H2", .before="UP")
RNA_apo_epiH <- rbind(RNA_apo_epiH1, RNA_apo_epiH2)
RNA_apo_epiH <- RNA_apo_epiH %>%
  # Creating an empty column:
  add_column(Condition = "High dose", .before="UP")
RNA_allnewepi_apo<-rbind(RNA_apo_epiC,RNA_apo_epiL,RNA_apo_epiH)
rownames(RNA_allnewepi_apo)<-NULL

RNA_allnewepi_apo$UP <- as.numeric(RNA_allnewepi_apo$UP)
RNA_allnewepi_apo$DOWN <- as.numeric(RNA_allnewepi_apo$DOWN)
RNA_allnewepi_apo$Total=rowSums(cbind(RNA_allnewepi_apo$UP,RNA_allnewepi_apo$DOWN),na.rm=FALSE)

RNA_allnewepi_apo$Condition <- fct_relevel(RNA_allnewepi_apo$Condition, c("Control", "Low dose", "High dose"))

#Sup.Fig.5
ggdensity(RNA_allnewepi_apo$Total, 
          main = "Density plot of epimutations_number",
          xlab = "Total number of epimutations_number")
shapiro.test(RNA_allnewepi_apo$Total)
kruskal.test(Total ~ Condition, data = RNA_allnewepi_apo)
dunnTest(Total ~ Condition, data = RNA_allnewepi_apo)

RNA_allLineage_rate_apo <- ggplot(RNA_allnewepi_apo, aes(x=Condition, y=Total, color = Condition,group=Condition)) + 
  geom_boxplot(fatten = 1, lwd = 1, width=0.5)+
  scale_color_manual(values=c("cornflowerblue", "darkgreen", "red"))+
  labs(y = "Number of new epimutations", x = "\nCondition")+
  geom_dotplot(binaxis='y', binwidth = 0.25, stackdir='center',  dotsize=1.5, binpositions = "all", stackgroups = TRUE)+
  theme_bw()+
  theme_linedraw()+
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(color = "black", size = 20, angle = 0, face = "plain"),
        axis.text.y = element_text(color = "black", size = 20, angle = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 16, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 16, angle = 90, face = "bold"))+
  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))+  
  scale_fill_brewer(palette="Pastel2")
RNA_allLineage_rate_apo 

#Get the raw data
write.xlsx(RNA_allnewepi_apo,"Data_Sup_Fig_5.xlsx")

#tRNAs epimutations associated to apoptotic genes
tRNAgenes<-read.xlsx("tRNA_targetedGenes2.xlsx")
tRNAgenes_apo <- filter(tRNAgenes, Name %in% apo)
