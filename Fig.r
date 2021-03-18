#!/usr/bin/R
#   _  ___     _____  _____ ____    ____  _   _    _    
#  | |/ / |   |__  / |_   _/ ___|  |  _ \| \ | |  / \   
#  | ' /| |     / /    | | \___ \  | |_) |  \| | / _ \  
#  | . \| |___ / /_    | |  ___) | |  _ <| |\  |/ ___ \ 
#  |_|\_\_____/____|   |_| |____/  |_| \_\_| \_/_/   \_\
#                                                       
#                     _           _     
#    __ _ _ __   __ _| |_   _ ___(_)___ 
#   / _` | '_ \ / _` | | | | / __| / __|
#  | (_| | | | | (_| | | |_| \__ \ \__ \
#   \__,_|_| |_|\__,_|_|\__, |___/_|___/
#                       |___/           
#    
# Library packages and loading data
setwd("$YOUR_DIR")
pacman::p_load(tidyverse, reshape2,scales, ggpubr, factoextra, FactoMineR, psych, vegan, patchwork, ggsci)
metadata <- read.delim("metadata.tsv", row.names = 1, stringsAsFactors = F)
species <- read.delim("species.tsv", row.names = 1, stringsAsFactors = F)
genus <- read.delim("genus.tsv", row.names = 1, stringsAsFactors = F)
### 
map <- metadata
map_A <- map[!duplicated(map$Individual_ID), ] %>% filter(sam2adm <= 3)
map_D <- map[!duplicated(map$Individual_ID, fromLast = T), ] %>% filter(dis2sam <= 7)

## Fig3D
####### Fig 3D #########
### Survive      ROC  AUC
##################
library(survminer)
library(survival)
#strep_list = c("Streptococcus parasanguinis","Streptococcus sp. oral taxon 431" "Streptococcus pneumoniae","Streptococcus salivarius","Streptococcus sanguinis","Streptococcus oralis","Streptococcus koreensis")
### Select gamlss  DE genus
df1 <- species[strep_list[1:7],rownames(map_A)] %>% as.matrix() # %>% .[unique(c(ts_t1_gam_list)),] %>% as.data.frame()
# KM surv
map_A$dis2adm <-  as.numeric(as.Date(map_A$Admission)-as.Date(map_A$Dischrge)) %>% abs
OS_m <- dplyr::select(map_A, Death, dis2adm) # dis2sym os_time dis2sam os_a2d
OS_m$Death=ifelse(OS_m$Death=="Deceased",1,0)
OS_m$os_time=OS_m$dis2adm
OS_m$os_time[OS_m$os_time>30]=30
mySurv <- with(OS_m, Surv(os_time, Death))
log_rank_p <- apply(df1, 1, function(gene) {
  # gene=exprSet[1,]
  OS_m$group <- ifelse(gene > median(gene), "high", "low")
  data.survdiff <- survdiff(mySurv ~ group, data = OS_m)
  p.val <- 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  return(p.val)
})
# log_rank_p=p.adjust(log_rank_p,method = "bonferroni")
table(log_rank_p < 0.05)
rf_list1 <- log_rank_p[log_rank_p < 0.05] %>% names()
rf_list1
df3 <- cbind(t(df1), OS_m)
i="Streptococcus parasanguinis"
  tmp <- dplyr::select(df3, Death, os_time, i)
  tmp$group <- ifelse(tmp[, i] > median(tmp[, i]), "High", "Low")
  mySurv_OS <- with(tmp, Surv(os_time, Death))
  sfit <- survfit(mySurv_OS ~ group, data = tmp)
  p <- ggsurvplot(sfit, pval = TRUE, conf.int = TRUE, surv.median.line = "hv",risk.table = TRUE,risk.table.col = "strata", ggtheme = theme_bw(), palette = c("#E7B800", "#2E9FDF")) + ggtitle(i)


### Figure S3A
df=genus[1,c(map_A_share$Sample_ID2,map_D_share$Sample_ID2)] %>% as.data.frame() %>% rownames_to_column("Sample_ID2") %>% right_join(metadata,.)
colnames(df)[ncol(df)]="value"
df$TT=ifelse(df$Sample_ID2 %in% map_A_share$Sample_ID2,"A","D")

df$value[df$value < 1e-6] <- 1e-6
ggpaired(df,
  x = "TT", y = "value", facet.by = "Group",
  line.color = "gray",line.size = 0.5,id = "Individual_ID", color = "TT", palette = "jama",
) +
  stat_compare_means(paired = TRUE)+  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  )


###Fig 3A
genus[1,filter(metadata,Death %in% c(0,1)) %>% rownames()] %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("Sample_ID2") %>%
  right_join(metadata, .) %>%
  mutate(Survival = as.factor(Death)) %>%
  mutate(Group = ifelse(Death == 1, "Deceased", "Survival")) %>%
  ggboxplot(., "sam2adm_correct", "Streptococcus", fill = "Survival", palette = Death_color) +
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) + stat_compare_means(aes(group = Survival), label = "p.signif") + ylab("Streptococcus relative abundance")+xlab("Days after admission")

## Fig 3C
df=species[strep_list[1:7], c(map_A %>% rownames(), filter(metadata,Death %in% c("health")) %>% rownames()] %>%as.data.frame() %>% 
  rownames_to_column("tax") %>%
  melt(variable.name = "Sample_ID2") %>%
  mutate(rela = value + 1e-5) %>%
  right_join(metadata, .)
##Summary different streptococcus species fold change 
df %>% group_by(tax,Group) %>% summarise(median=median(rela)) %>% dcast(tax~Group) %>% mutate(Healthy_Deceased= Healthy/Deceased,Healthy_Recovered= Healthy/Recovered,Recovered_Deceased= Recovered/Deceased)
df$Death[df$Death=="health"]="Healthy"
df$Death[df$Death=="1"]="Deceased"
df$Death[df$Death=="0"]="Recovered"
  df=mutate(df,Survival = as.factor(Death))
  df$Survival=factor(df$Survival,levels=c("Healthy","Recovered","Deceased"))
  ggboxplot(df, "tax", "rela", fill = "Survival", palette = Death_color[c(3,1,2)]) +
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) + stat_compare_means(aes(group = Survival), label = "p.signif") +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))

  my_comparisons=list(c("Healthy","Recovered"),c("Recovered","Deceased"),c("Healthy","Deceased"))
ggboxplot(df,"Survival", "rela", fill = "Survival", palette = Death_color[c(3,1,2)],facet.by = "tax")+
  stat_compare_means(comparisons =my_comparisons,label = "p.signif" )
  


###Fig S3B
df=species[strep_list[1:7], c(map_D %>% rownames(), filter(metadata,Death %in% c("health")) %>% rownames()] %>%
  rownames_to_column("tax") %>%
  melt(variable.name = "Sample_ID2") %>%
  mutate(rela = value + 1e-5) %>%
  right_join(metadata, .)
df$Death[df$Death=="health"]="Healthy"
df$Death[df$Death=="1"]="Deceased"
df$Death[df$Death=="0"]="Recovered"
df=mutate(df,Survival = as.factor(Death))
df$Survival=factor(df$Survival,levels=c("Healthy","Recovered","Deceased"))
ggboxplot(df, "tax", "rela", fill = "Survival", palette = Death_color[c(3,1,2)]) +
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) + stat_compare_means(aes(group = Survival), label = "p.signif") +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))

my_comparisons=list(c("Healthy","Recovered"),c("Recovered","Deceased"),c("Healthy","Deceased"))
ggboxplot(df,"Survival", "rela", fill = "Survival", palette = Death_color[c(3,1,2)],facet.by = "tax")+
  stat_compare_means(comparisons =my_comparisons,label = "p.signif" )


## Fig S6 x
# alliviual
library(ggalluvial)
dd1 <-  species[, filter(metadata, Death == 0) %>% rownames()] %>% as.matrix() 
dd2 <- dd1[order(rowSums(dd1), decreasing = T), ]
dd3 <- dd2 %>%
  as.data.frame() %>%
  rownames_to_column("taxname")
dd3$taxname[16:nrow(dd3)] <- "Others"
dd4 <- dd3 %>%
  group_by(taxname) %>%
  summarise_each(funs = sum) %>%
  as.data.frame()
dd4$taxname <- factor(dd4$taxname, levels = rev(dd3$taxname[1:16]))
dd4 <- dd4 %>% melt()
colnames(dd4)[2] <- "Sample_ID2"
dd5 <- right_join(metadata, dd4)
dd5$sam2adm_correct[is.na(dd5$sam2adm_correct)] <- "Healthy"
#dd5$group <- ifelse(grepl("KLZ", dd5$Sample_ID), "T1", "Healthy")
dd6 <- dd5 %>%
  group_by(taxname, sam2adm_correct) %>%
  summarise(rela = mean(value)) # %>% filter(sam2adm_correct!= 21)
dd6$sam2adm_correct <- as.factor(dd6$sam2adm_correct)
dd6$sam2adm_correct <- factor(dd6$sam2adm_correct, levels = c("1", "5", "10", "14", "21", "Healthy"))
ggplot(dd6 %>% filter(sam2adm_correct != 21), aes(sam2adm_correct, rela, alluvium = taxname, stratum = taxname, fill = taxname)) +
  geom_alluvium(aes(fill = taxname, colour = taxname)) +
  geom_stratum() +
  scale_fill_manual(values = inlmisc::GetColors(16, scheme = "discrete rainbow") %>% as.character()) +
  theme_bw()



## Fig S4x
df1 <- species[, map_A_share$Sample_ID2]
df2 <- species[, map_D_share$Sample_ID2]
df21 <- df2 %>%
  as.data.frame() %>%
  rownames_to_column("tax") %>%
  melt(variable.name = "Sample_ID2") %>%
  mutate(group1 = "T2") %>%
  right_join(metadata, .)
df11 <- df1 %>%
  as.data.frame() %>%
  rownames_to_column("tax") %>%
  melt(variable.name = "Sample_ID2") %>%
  mutate(group1 = "T1") %>%
  right_join(metadata, .)
df3 <- rbind(df11, df21) %>% filter(Death != 1) # %>% filter(Individual_ID %in% map_tm_share$Individual_ID)
dd <- dplyr::select(df3 %>% filter(Severity.D. == 2), tax, group1, value)
ll=intersect(metadata1$Sample_ID2,colnames(wh_ncov_genus_p))
dd1 <- melt(species[,ll] %>% as.data.frame() %>% rownames_to_column("tax")) %>%
  mutate(group1 = "Healthy") %>%
  dplyr::select(1, 4, 3) %>%
  rbind(dd, .)
ll=species[rowMeans(species[,unique(df3$Sample_ID2)])>0.01,] %>% rownames()
ll1=df3 %>% group_by(group1,tax) %>% dplyr::summarise(aa=median(value)) %>% dcast(tax~group1) %>% filter(T2>T1) %>% .[,1] %>% intersect(.,ll)

compare_means(value ~ group, data = df3 %>% filter(Death == "Recovered") %>% filter(Severity.D. == 2), paired = T, group.by = "tax", id = "Individual_ID") %>%
  as.data.frame() %>%
  filter(p.signif != "ns") %>% filter(tax %in% intersect(ll,ll1))

A_D_up=compare_means(value ~ group1, data = df3 %>% filter(Death == 0) %>% filter(Severity.D. == 2), paired = T, group.by = "tax", id = "Individual_ID", p.adjust.method = "fdr") %>%
  as.data.frame() %>%
  filter(p.signif != "ns") %>% filter(tax %in% ll1)


df1 <- species[A_D_up$tax,map_A_share$Sample_ID2]
df3 <- species[A_D_up$tax,map_D_share$Sample_ID2]
df11 <- df1 %>%
  as.data.frame() %>%
  rownames_to_column("tax") %>%
  melt(variable.name = "Sample_ID2") %>%
  mutate(group1 = "T1") %>%
  right_join(metadata, .)
df31 <- df3 %>%
  as.data.frame() %>%
  rownames_to_column("tax") %>%
  melt(variable.name = "Sample_ID2") %>%
  mutate(group1 = "T2") %>%
  right_join(metadata, .)
df <- rbind(df11, df31)
df$value[df$value < 1e-5] <- 1e-5
df$value <- log10(df$value)
kk1 <- df %>%
  filter(Death != 1) %>%
  filter(Severity.D. == 2) %>%
  dplyr::select(group1, tax, Individual_ID, value) %>%
  dcast(Individual_ID + tax ~ group1) %>%
  filter(T2 > T1)
kk2 <- df %>%
  filter(Death != 1) %>%
  filter(Severity.D. == 2) %>%
  dplyr::select(group1, tax, Individual_ID, value) %>%
  dcast(Individual_ID + tax ~ group1) %>%
  filter(T2 < T1)

df1 <- df %>%
  filter(Death != 1) %>%
  filter(Severity.D. == 2) %>%
  dplyr::select(group1, value, tax)
dd <- species[A_D_up$tax,intersect(metadata1$Sample_ID2,colnames(wh_ncov_genus_p))] %>%as.data.frame() %>% 
  rownames_to_column("tax") %>%
  melt() %>%
  mutate(group1 = "Healthy") %>%
  dplyr::select(4, 3, 1)
dd$value <- log10(dd$value)
df2 <- rbind(df1, dd)
df2$group <- factor(df2$group, levels = c("T1", "T2", "Healthy"))
ggplot(df2 %>% filter(tax!="Lautropia mirabilis"), aes(group, value, color = group)) +
  geom_boxplot(size = 1.5) +
  facet_wrap(. ~ tax) +
  theme_bw() +
  geom_segment(data = kk2%>% filter(tax!="Lautropia mirabilis"), aes(x = "T1", y = T1, xend = "T2", yend = T2), color = "grey45", alpha = 0.25, linetype = 3, size = 0.3) +
  geom_segment(data = kk1%>% filter(tax!="Lautropia mirabilis"), aes(x = "T1", y = T1, xend = "T2", yend = T2), color = "grey45", alpha = 0.3, size = 0.3) +
  scale_color_jama() +
  stat_compare_means(comparisons = mm,label = "p.signif")+theme(text=element_text(family="Arial"))





### RNAseq analysis
df3=read.delim("rna/all_gene.tsv",row.names = 1)
ll1=colnames(df3)[colSums(df3)>100000]
library('AnnotationDbi')
library('org.Hs.eg.db')
library(tidyverse)
df3=df3[rowSums(df3>1)>5,colSums(df3)>100000]
## ll1  >100000   ll3 >300000   ll5>50000
ll=intersect(ll1,rownames(map_A)) #%>% c(.,ll1[100:102])
# ll=c(intersect(ll5,new_t2_list$T2),ll5[175:180])
# All gene count noramlization
S_raw <- df3[rowSums(df3)>10,ll]
dge <- DGEList(counts=S_raw)

comp2groups <- list(Recovered_Deceased=c("Deceased","Recovered") )
DEA_list <- list()

for (comp_var in names(comp2groups) ) {
  # comp_var <- names(comp2groups)[[1]]
  comp_i <- comp2groups[[comp_var]]
  print (comp_i)
  sample_meta=metadata[ll,]
  subsample_pheno <- dplyr::filter(sample_meta, Group %in% comp_i )
  #subsample_pheno$Group[!grepl("Healthy",subsample_pheno$Group)]="nCOV"
  subsample_pheno$Group <- factor(subsample_pheno$Group, levels=comp_i )
  Expdesign <- model.matrix(~subsample_pheno$Group)
  
  subsample_dge <- dge[,subsample_pheno$Sample_ID2]
  
  # DEA
  keep <- filterByExpr(subsample_dge, Expdesign)
  subsample_dge <- subsample_dge[keep,keep.lib.sizes=FALSE]
  subsample_dge <- calcNormFactors(subsample_dge)
  
  v <- voom(subsample_dge, Expdesign, plot=FALSE, normalize="quantile")
  
  # boxplot(v$E[,nCoV_samp])
  # boxplot(v$E)
  
  Expfit1 <- lmFit(v, Expdesign)
  Expfit2 <- eBayes(Expfit1)
  limma_res <- topTable(Expfit2, coef=tail(colnames(Expdesign), 1), number=Inf) %>% 
    tibble::rownames_to_column(var=comp_var )
  
  limma_DEA <- dplyr::filter(limma_res, adj.P.Val<=0.05, abs(logFC)>=1.5)
  #readr::write_tsv(limma_DEA, paste0("../results/", comp_var, "_limma_res.tsv") )
  
  limma_UP <- dplyr::filter(limma_DEA, logFC>0)[[comp_var]]
  limma_DN <- dplyr::filter(limma_DEA, logFC<0)[[comp_var]]
  
  limma_DEG <- limma_DEA[[comp_var]]
  print(length(limma_DEG))
  limma_DE <- v$E[limma_DEG,]
  
  
  DEA_list[["limma_res"]][[comp_var]] <- limma_res
  DEA_list[["limma_DEA"]][[comp_var]] <- limma_DEA
  DEA_list[["limma_DEG"]][[comp_var]] <- limma_DEG
  DEA_list[["limma_UPDN"]][[paste0(comp_var, "_UP")]] <- limma_UP
  DEA_list[["limma_UPDN"]][[paste0(comp_var, "_DN")]] <- limma_DN
  DEA_list[["limma_DE"]][[comp_var]] <- limma_DE
  DEA_list[["subsample_pheno"]][[comp_var]] <- subsample_pheno
  
}

library(clusterProfiler)
gene<-DEA_list$limma_DEG$Recovered_Deceased
gene.df<-bitr(gene, fromType = "SYMBOL", 
              toType = c("SYMBOL","ENTREZID","ENSEMBL"),
              OrgDb = org.Hs.eg.db)

#Go enrichment
ego_bp<-enrichGO(gene       = gene.df$ENTREZID %>% unique(),
                 OrgDb      = org.Hs.eg.db,
                 keyType    = 'ENTREZID',
                 ont        = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.01,
                 qvalueCutoff = 0.05)

barplot(ego_bp,showCategory = 15,title="The GO_BP enrichment analysis of all DEGs ")+ 
  scale_size(range=c(2, 12))+
  scale_x_discrete(labels=function(ego_bp) str_wrap(ego_bp,width = 20)).

### Fig S8x
df=v$E[DEA_list$limma_DEG$Recovered_Deceased,] %>% melt()
colnames(df)=c("symbol","Sample","value")
df1=right_join(subsample_pheno,df)
ggboxplot(df1,"Group","value",facet.by = "symbol",fill = "Group",
          palette = Death_color[c(2,1)],add = "jitter",width = 0.5)+
  ylab("Normalized gene expression")+xlab("")+theme_cleveland()