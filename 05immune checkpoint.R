

genes <- c("PDCD1", "CD274", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "IDO1")

Tcells<-readRDS('2Tcells.rds')
Tcells@meta.data$sample <- as.character(Tcells@meta.data$sample)
Tcells@meta.data$group[Tcells@meta.data$sample == 'MM3_WHO2'] <- 'WHO3'
Tcells@meta.data$sample[Tcells@meta.data$sample == 'MM3_WHO2'] <- 'MM3_WHO3'

expr_mat <- GetAssayData(Tcells,assay = DefaultAssay(Tcells),slot = "data")

#expression percentage
result_list <- list()
for (gene in genes) {
  if (gene %in% rownames(expr_mat)) {
    expr <- expr_mat[gene, ]
    positive <- expr > 0
    prop <- tapply(positive,Tcells$group,mean,na.rm = TRUE) * 100
    tmp <- data.frame(sample = names(prop),gene = gene,positive_cells = as.numeric(prop))
    result_list[[gene]] <- tmp

  } else {
    cat("Skipping ",gene," - gene not found in expression matrix.\n",sep = "")
  }
}
gene_positive_prop <- bind_rows(result_list)
gene_positive_prop
   sample   gene positive_cells
1    WHO1  PDCD1      1.4477766
2    WHO2  PDCD1      1.3028839
3    WHO3  PDCD1      0.8550663
4    WHO1  CD274      0.4825922
5    WHO2  CD274      1.0540184
6    WHO3  CD274      0.7909363
7    WHO1  CTLA4      1.2639320
8    WHO2  CTLA4      3.8647343
9    WHO3  CTLA4      2.9927319
10   WHO1   LAG3      4.5501551
11   WHO2   LAG3      5.5775143
12   WHO3   LAG3      4.1470714
13   WHO1 HAVCR2      2.4014708
14   WHO2 HAVCR2      3.8793734
15   WHO3 HAVCR2      1.6460026
16   WHO1  TIGIT      4.6995289
17   WHO2  TIGIT     11.5356463
18   WHO3  TIGIT      3.5057717
19   WHO1   IDO1      0.2987476
20   WHO2   IDO1      0.6880398
21   WHO3   IDO1      0.4061565







library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)

df <- tribble(
  ~Method, ~Gene, ~WHO, ~Expression,
  # === PDCD1 ===
  "stRNA-seq", "PDCD1", "WHO1", 0.05814184,
  "stRNA-seq", "PDCD1", "WHO2", 0.09580126,
  "stRNA-seq", "PDCD1", "WHO3", 0.18401328,
  "snRNA-seq", "PDCD1", "WHO1", 0.35383997,
  "snRNA-seq", "PDCD1", "WHO3", 0.13812155,

  # === CD274 ===
  "stRNA-seq", "CD274", "WHO1", 0.26494179,
  "stRNA-seq", "CD274", "WHO2", 0.97783357,
  "stRNA-seq", "CD274", "WHO3", 0.46227728,
  "snRNA-seq", "CD274", "WHO1", 0.28950543,
  "snRNA-seq", "CD274", "WHO3", 0.00000000,

  # === CTLA4 ===
  "stRNA-seq", "CTLA4", "WHO1", 1.61607886,
  "stRNA-seq", "CTLA4", "WHO2", 1.15457038,
  "stRNA-seq", "CTLA4", "WHO3", 1.05471029,
  "snRNA-seq", "CTLA4", "WHO1", 1.02935263,
  "snRNA-seq", "CTLA4", "WHO3", 1.93370166,

  # === LAG3 ===
  "stRNA-seq", "LAG3", "WHO1", 0.22199612,
  "stRNA-seq", "LAG3", "WHO2", 1.17439133,
  "stRNA-seq", "LAG3", "WHO3", 1.05471029,
  "snRNA-seq", "LAG3", "WHO1", 0.28950543,
  "snRNA-seq", "LAG3", "WHO3", 0.27624309,

  # === HAVCR2 ===
  "stRNA-seq", "HAVCR2", "WHO1", 9.10250142,
  "stRNA-seq", "HAVCR2", "WHO2", 15.01106670,
  "stRNA-seq", "HAVCR2", "WHO3", 5.36780216,
  "snRNA-seq", "HAVCR2", "WHO1", 33.63892240,
  "snRNA-seq", "HAVCR2", "WHO3", 23.61878453,

  # === TIGIT ===
  "stRNA-seq", "TIGIT", "WHO1", 0.08523065,
  "stRNA-seq", "TIGIT", "WHO2", 0.16352284,
  "stRNA-seq", "TIGIT", "WHO3", 0.16157264,
  "snRNA-seq", "TIGIT", "WHO1", 0.08845999,
  "snRNA-seq", "TIGIT", "WHO3", 0.00000000,

  # === IDO1 ===
  "stRNA-seq", "IDO1", "WHO1", 0.14799741,
  "stRNA-seq", "IDO1", "WHO2", 0.23454792,
  "stRNA-seq", "IDO1", "WHO3", 0.28275212,
  "snRNA-seq", "IDO1", "WHO1", 0.78809811,
  "snRNA-seq", "IDO1", "WHO3", 0.96685083
)


df$Method <- factor(df$Method, levels = c("stRNA-seq", "snRNA-seq"))
df$WHO <- factor(df$WHO, levels = c("WHO1", "WHO2", "WHO3"))
df$group <- factor(
  interaction(df$WHO, df$Method),
  levels = c(
    "WHO1.stRNA-seq", "WHO1.snRNA-seq",
    "WHO2.stRNA-seq", "WHO2.snRNA-seq",
    "WHO3.stRNA-seq", "WHO3.snRNA-seq"
  )
)

# 绘图
ggplot(df, aes(x = Gene, y = Expression, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(
    name = "Method & WHO",
    values = c(
      "WHO1.stRNA-seq" = "#4575b4",  # 蓝
      "WHO1.snRNA-seq" = "#74add1",  # 浅蓝
      "WHO3.stRNA-seq" = "#f46d43",  # 橙
      "WHO3.snRNA-seq" = "#fdae61",  # 浅橙
      "WHO2.stRNA-seq" = "#66bd63",  # 绿
      "WHO2.snRNA-seq" = "#b2df8a"   # 浅绿
    ),
    labels = c(
      "WHO1 stRNA-seq", "WHO1 snRNA-seq",
      "WHO2 stRNA-seq", "WHO2 snRNA-seq",
      "WHO3 stRNA-seq", "WHO3 snRNA-seq"
    )
  ) +
  labs(x = NULL, y = "Expression proportion (%)", fill = "") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, face = "bold"),
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )


df_cn <- tribble(
  ~Ecotype, ~Proportion,
  "CN1", 10.85,
  "CN2", 12.37,
  "CN3", 10.45,
  "CN4", 0.79,
  "CN5", 4.99,
  "CN6", 10.36,
  "CN7", 11.36
)
ggplot(df_cn, aes(x = Ecotype, y = Proportion, group = 1)) +
  geom_line(color = "#4575b4", size = 1.2) +   
  geom_point(size = 3, color = "#4575b4") +   
  geom_text(aes(label = round(Proportion, 2)), vjust = -0.8, size = 4) +  
  labs(x = "Ecological niche", y = "Proportion (%)") +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(face = "bold")
  )


  