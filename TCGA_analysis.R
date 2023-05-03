### Date:2023.5.1
### Author by zhangjk

### load packages
library(dplyr,quietly = T,verbose = F)
library(ggplot2,quietly = T,verbose = F)
library(survival,quietly = T,verbose = F)

rm(list=ls())
setwd("~")

### load data
TCGA_fpkm_exp <- "GDC-PANCAN.htseq_fpkm-uq.tsv" ### obtained from UCSC Xena 
TCGA_sur <- "Survival_SupplementalTable_S1_20171025_xena_sp" ### obtained from UCSC Xena 
TCGA_cBioportal <- "combined_study_clinical_data.tsv" ### obtained from cBioportal
TCGA_cnv <- "GDC-PANCAN.masked_cnv.tsv"

TCGA_fpkm_exp <- read.table(file = TCGA_fpkm_exp,header = T,sep = "\t",check.names = F,row.names = 1)
TCGA_sur <- read.csv(file = TCGA_sur,header = T,sep = "\t",check.names = F)
TCGA_cBioportal <- read.csv(file = TCGA_cBioportal,header = T,sep = "\t",check.names = F)

### Here,we select three genes for downstream analysis.
tem <- data.frame(id=(colnames(TCGA_fpkm_exp)),t(TCGA_fpkm_exp[c("ENSG00000149218",
                                                                 "ENSG00000012048",
                                                                 "ENSG00000152457"),]))
colnames(tem)[2:4] <- c("ENDOD1","BRCA1","DCLRE1C")
tem$id <- substring(tem$id,1,15)
tem <- tem %>% group_by(id) %>%
  summarise_each(funs(mean))

tem <- merge(tem,TCGA_cBioportal,by.x = "id",by.y = "Sample ID")
tem <- merge(tem,TCGA_sur[,c("sample","age_at_initial_pathologic_diagnosis",
                             "OS","OS.time","gender","race")],by.x = "id",by.y = "sample")
#tem <- tem[!is.na(tem$`Fraction Genome Altered`),]



### pan-cancer correlation plot 
pan_ana <- function(data,index) {
  indicator <- switch(index,
                      TMB = "TMB (nonsynonymous)",
                      MSI_sensor = "MSIsensor Score",
                      MSI_MANTIS = "MSI MANTIS Score",
                      CNV = "Aneuploidy Score",
                      CNA_burden = "Fraction Genome Altered")
  cor <- list()
  pvalue <- list()
  n = 1
  options(warn = -1)
  for (i in unique(data$`TCGA PanCanAtlas Cancer Type Acronym`)) {
    #print(i)
    tt <- subset(data,data$`TCGA PanCanAtlas Cancer Type Acronym` == i)
    tt <- tt[!is.na(tt[,indicator]),]
    if (dim(tt)[1] > 10) {
      test <- cor.test(tt$ENDOD1,tt[,indicator],method = "spearman")
      cor[[n]] <- as.vector(test[[4]])
      pvalue[[n]] <- as.vector(test[[3]])
      n = n + 1
    } else  {
      cor[[n]] <- NA
      pvalue[[n]] <- NA
    }
    n = n + 1
  }
  
  tem1 <- data.frame(type=unique(data$`TCGA PanCanAtlas Cancer Type Acronym`),cor=unlist(cor),pvalue=unlist(pvalue))
  tem1[,2:3] <- lapply(tem1[,2:3], function(x) as.numeric(x))
  tem1 <- na.omit(tem1)
  tem1 <- tem1 %>%
    mutate(p_color = case_when(
      (pvalue <= 0.05 & cor > 0) ~  "Red",
      (pvalue <= 0.05 & cor < 0) ~ "Blue",
      pvalue > 0.05 ~ "Grey"))
  
  tem1 <- tem1[order(tem1$cor,decreasing = T),]
  pan <- as.data.frame(t(c("PanCan \n Average",mean(tem1[tem1$pvalue <= 0.05,]$cor),"0.05","Yellow")))
  colnames(pan) <- colnames(tem1)
  tem1 <- rbind(tem1,pan)
  tem1$cor <- as.numeric(tem1$cor)
  tem1$p_color <- factor(tem1$p_color,levels = c("Red","Grey","Blue","Yellow"),labels = c("Red","Grey","Blue","Yellow"))
  tem1$type <- factor(tem1$type, levels=tem1$type, ordered=TRUE)
  ggplot(data = tem1,mapping = aes(x = type, y = cor,fill = p_color,color = "black")) +
    geom_bar(stat="identity", color="black", position=position_dodge(), width = 0.65) +
    #coord_cartesian(ylim = c(1500,2000)) +
    theme_bw() +
    theme(axis.line = element_line(linetype = "solid"),
          axis.title = element_text(size = 17),
          axis.ticks = element_line(size = 2),
          axis.text.x = element_text(family = "Times",angle = 90,
                                     size = 15),
          axis.text.y = element_text(size = 12),
          legend.position = c(0.7,0.8)) +
    scale_fill_manual(values=c("red","grey","blue","yellow"),
                      labels = c("Red","Grey","Blue","Yellow")) +
    labs(fill = "Correlation P-value", size = 13,x = NULL, y = paste("Spearman\'s Rho \n ( ENDOD1", " vs  ",index ,")",sep = " "))
}

pan_ana(tem,"TMB")
pan_ana(tem,"MSI_sensor")
pan_ana(tem,"MSI_MANTIS")
pan_ana(tem,"CNV")
pan_ana(tem,"CNA_burden")


### CNV plot
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x <= (qnt[1] - H)] <- NA
  y[x >= (qnt[2] + H)] <- NA
  y
}

give.n <- function(x){
  return(c(y = median(x)*1.1, label = length(x)))
  # experiment with the multiplier to find the perfect position
}


CNV <- read.table(TCGA_cnv,header = T,sep = "\t",check.names = F)
cancer_type <- "HNSC"

CNV$Chrom <- paste("chr",CNV$Chrom,sep = "")
CNV <- CNV[,c(2:4,1,5)]
CNV$sample <- substring(CNV$sample,1,15)

inter <- CNV[CNV$sample %in% 
                 intersect(CNV$sample,tem[tem$`TCGA PanCanAtlas Cancer Type Acronym` == cancer_type,]$id) ,]

inter$id <- paste(inter$Chrom,inter$Start,inter$End,sep = "_")
inter_cnv <- inter[,c(4,6)]

#inter_cnv <- subset(inter,inter$V6 != ".")
#inter_cnv <- inter_cnv[,c(4,10)]

inter_cnv <- unique(inter_cnv)
a <- aggregate(inter_cnv$id,by = list(inter_cnv$sample),FUN = length)
tt <- merge(tem,a,by.x = "id",by.y = "Group.1")

tt$group <- ifelse(tt$ENDOD1 >= median(tt$ENDOD1),ifelse(tt$BRCA1 >= median(tt$BRCA1),"ENDOD1 High \n BRCA1 High","High \n Low"),
                   ifelse(tt$BRCA1 >= median(tt$BRCA1),"Low \n High","Low \n Low"))
table(tt$group)

tt <- tt %>%
  group_by(group) %>%
  mutate(x = remove_outliers(x))
head(tt)

ggplot(tt,aes(x = group,y = log2(x),group = group)) +
  geom_violin(aes(fill = group), color = "white",how.legend = FALSE) +
  geom_boxplot(outlier.shape = NA,notch = T,notchwidth = 0.7,width=.2, size=1.05) + xlab(label = NULL) + ylab(bquote(paste("Log"[2],"(Masked CNV counts)" ,sep = ""))) +
  ggtitle(NULL) +
  stat_summary(fun.data = give.n, geom = "text", fun.y = median, color = "red",face="bold",size=6,
               position = position_dodge(width = 0.75)) +
  scale_fill_brewer(palette="RdPu") + theme_grey() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),axis.title = element_text(size = 18),strip.text = element_text(size = 16),
        legend.position = "none") +
  geom_signif(comparisons = list(c("Low \n High","ENDOD1 High \n BRCA1 High","High \n Low"),
                                 c("ENDOD1 High \n BRCA1 High","High \n Low"),
                                 c("Low \n High","Low \n Low")),
              map_signif_level=function(p) paste("P = ",format(p,scientific=TRUE,digits = 2),sep = ""),textsize=4,test=wilcox.test,step_increase=0.15, fontface = 'italic') +
  facet_wrap(~ type,scales = "free")