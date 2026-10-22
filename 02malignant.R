
library(dplyr)
library(ggsci)
library(ggplot2)

malignant<-readRDS('1malignant_CytoTRACE2.rds')
malignant@meta.data$sample <- as.character(malignant@meta.data$sample)
malignant@meta.data$group[malignant@meta.data$sample == 'MM3_WHO2'] <- 'WHO3'
malignant@meta.data$sample[malignant@meta.data$sample == 'MM3_WHO2'] <- 'MM3_WHO3'


#cell circle
malignant <- CellCycleScoring(malignant, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)

samples <- unique(malignant@meta.data$group)
cor_results <- data.frame(Sample = character(), Correlation = numeric(), P.value = numeric())
for (s in samples) {
  subset_cells <- malignant@meta.data$group == s
  df <- malignant@meta.data[subset_cells, c("S.Score", "G2M.Score")]
  cor_test <- cor.test(df$S.Score, df$G2M.Score, method = "pearson")
  cor_results <- rbind(cor_results, data.frame(
    Sample = s,
    Correlation = cor_test$estimate,
    P.value = cor_test$p.value
  ))
  p <- ggplot(df, aes(x = S.Score, y = G2M.Score)) +
    geom_point(alpha = 0.4) +
    geom_smooth(method = "lm", color = "red") +
    ggtitle(paste0("Sample: ", s, 
                   "\nCorrelation: ", round(cor_test$estimate, 3), 
                   ", p = ", signif(cor_test$p.value, 3))) +
    theme_minimal()
  ggsave(paste0("S_vs_G2M_", s, ".pdf"), p, width = 5, height = 5)
}
print(cor_results)


df <- malignant@meta.data[, c("niche","group","S.Score", "G2M.Score")]
df$group <- factor(df$group, levels = c("WHO3", "WHO2", "WHO1"))

library(ggplot2)
library(viridis)  
p <- ggplot(df, aes(x = S.Score, y = G2M.Score)) +
  geom_hex(bins = 50) +
  geom_smooth(aes(color = group), method = "lm", se = TRUE, linewidth = 1) +
  scale_fill_gradient(low = "gray95", high = "black", name = "Density")+
  scale_color_npg() +
  labs(title = NULL,
       x = "S Score", y = "G2M Score", color = "Group") +
  theme_minimal(base_size = 13)
ggsave(paste0("S_vs_G2M_all.pdf"), p, width = 6, height = 5)








#Hallmarker
library(GSVA)
library(tidyverse)
library(clusterProfiler)
library(msigdbr)
library(patchwork)
library(RColorBrewer)
genesets <- msigdbr(species = "Homo sapiens", category = "H") 
genesets <- subset(genesets, select = c("gs_name","gene_symbol")) %>% as.data.frame()
genesets <- split(genesets$gene_symbol, genesets$gs_name)

Idents(malignant) <- "niche" 
malignant_sub <- subset(malignant, idents=c('CN1','CN3','CN6','CN7'))
expr <- AverageExpression(malignant_sub, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,]  #选取非零基因
expr <- as.matrix(expr)
head(expr)

gsva.res <- gsva(expr, genesets, method="gsva",parallel.sz=1) #,method = "zscore"速度快'gsva'
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
rownames(gsva.res)<-substring(rownames(gsva.res), 10)#去除开头的HALLMARKER
rownames(gsva.res)<-gsub("_"," ",rownames(gsva.res))
gsva.res<-data.frame(gsva.res)
gsva.res <- gsva.res[, c('CN6', 'CN7', 'CN3', 'CN1')]
pdf("1malignant_GSVA_niche.pdf")
pheatmap::pheatmap(gsva.res, show_colnames = T, cluster_rows = T, cluster_cols = F,
                   angle_col = "45",scale='row',
                   color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(50),
                   border=F,cellwidth=30)
dev.off()






