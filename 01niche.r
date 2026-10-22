cell_colors <- c(
  "Tumor_Multipotent"    = "#FF0049",
  "Tumor_Oligopotent"    = "#E6553A",
  "Tumor_Unipotent"      = "#FCAE91",
  "Tumor_Differentiated" = "#FDD0A2",
  "T_Cells"              = "#00F3D8",
  "Myeloid_cells"        = "#00ACEB",
  "Stromal_cells"        = "#87B400",
  "Necrotic_Cells"       = "#ABABAB"
)

samples<-c('MM1_WHO2','MM2_WHO2','MM3_WHO2','MM4_WHO3','MM5_WHO1','MM6_WHO3','MM7_WHO1','MM8_WHO1' )
pdf("0MSD_spatial_celltype.pdf", height = 6, width = 8)
for (s in samples) {
  cells_to_plot <- rownames(MSD@meta.data)[MSD@meta.data$sample == s]
  p <- DimPlot(object = MSD, reduction = "spatial", pt.size = 1,
               label = FALSE, group.by = 'cell_labels_updated', cells = cells_to_plot,raster=FALSE) +
    scale_color_manual(values = cell_colors) +
    labs(x = "spatial1", y = "spatial2", title = s) + 
    theme(panel.border = element_rect(fill = NA, color = "black", size = 1),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position = 'right',
          legend.text = element_text(color = "black", size = 16, face = "bold"))
  print(p)
}
dev.off()

#
cell_types <- levels(MSD@meta.data$cell_labels)
cell_colors <- setNames(rep("#D3D3D3", length(cell_types)), cell_types)
cell_colors["Myeloid_cells"] <- "#00ACEB"

pdf("0MSD_spatial_celltype.pdf", height = 10, width = 10)
for (s in samples) {
  cells_to_plot <- rownames(MSD@meta.data)[MSD@meta.data$sample == s]
  p <- DimPlot(object = MSD, reduction = "spatial", pt.size = 2,
               label = FALSE, group.by = 'cell_labels', cells = cells_to_plot, raster = FALSE) +NoLegend()+
    scale_color_manual(values = cell_colors) +  
    labs(x = "spatial1", y = "spatial2", title = s) + 
    theme(
      panel.border = element_blank(),axis.line=element_blank(),
      axis.text = element_blank(),axis.ticks = element_blank(),

    )
  print(p)
}
dev.off()


pdf("0MSD_RNA_snn_res.0.6.anno.spatial.pdf",height=20,width=12) 
DimPlot(object = MSD, reduction = "spatial",label = T,raster=FALSE)+scale_color_igv()+
        labs(x = "spatial1", y = "spatial2",title = "") + 
        theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
              axis.text.y = element_blank(),axis.ticks.y = element_blank(),
              axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position='right',
              legend.text=element_text(color = "black",size=16,face="bold"))
dev.off()

test<- as.data.frame(table(MSD@meta.data$cell_labels,MSD@meta.data$sample))
test$Var2 <- factor(test$Var2, levels=c("MM5_WHO1", "MM7_WHO1", "MM8_WHO1", "MM1_WHO2", "MM2_WHO2", "MM3_WHO2", "MM4_WHO3", "MM6_WHO3"), ordered=TRUE)
pdf("1MSD_celltype_per_percent_of_sample.pdf",width = 8,height=8)
ggplot(test,aes(x= Var2,y=Freq,fill= Var1))+
  geom_col( position="fill")+scale_y_continuous(expand = expansion(mult = c(0.02,0)),labels = scales::percent_format())+
  scale_x_discrete(labels= c("MM5_WHO1", "MM7_WHO1", "MM8_WHO1", "MM1_WHO2", "MM2_WHO2", "MM3_WHO2", "MM4_WHO3", "MM6_WHO3"))+ 
  labs(x='',y='Percentage')+guides(fill=guide_legend(ncol=4))+#图例分3列
  theme_classic()+
  theme(panel.background=element_rect(fill='transparent',color ="white"),
        legend.position = 'top',legend.direction = "horizontal",legend.text=element_text(color = "black",size=11,face="bold"),
        legend.title=element_blank(),legend.justification = 0.9,#position是‘top’时进行左右微调，适应字体大小
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,color = "black",size=18,face="bold"),#X轴字体
        axis.text.y = element_text(color = "black",size=15,face="bold"), #Y轴字体
        axis.title.y = element_text(color = "black",size=20,face="bold"),#Percentage的大小
        axis.ticks.x = element_blank(),axis.ticks.y = element_line(size=1.5),#XY轴刻度调整
        axis.line.x=element_line(linetype=1,color="black",size=1.5),axis.line.y=element_line(linetype=1,color="black",size=1.5))+ 
  scale_fill_npg()
dev.off()

#for Xenium Expo
samples <- c("MM1_WHO2", "MM2_WHO2", "MM3_WHO2", "MM4_WHO3", "MM5_WHO1", "MM6_WHO3", "MM7_WHO1", "MM8_WHO1")
for (s in samples) {
  cat("Processing", s, "\n")
  MM <- subset(MSD, sample == s)
  meta_df <- data.frame(cell_id = MM@meta.data$cell_id,group = MM@meta.data$cell_labels_updated)
  outname <- paste0("00", s, "_cell_labels.csv")
  write.csv(meta_df, outname, row.names = FALSE, quote = FALSE)
}

samples <- c("MM1_WHO2", "MM2_WHO2", "MM3_WHO2", "MM4_WHO3", "MM5_WHO1", "MM6_WHO3", "MM7_WHO1", "MM8_WHO1")
for (s in samples) {
  cat("Processing", s, "\n")
  MM <- subset(MSD, sample == s)
  meta_df <- data.frame(cell_id = MM@meta.data$cell_id,group = MM@meta.data$niche)
  outname <- paste0("00", s, "_niche.csv")
  write.csv(meta_df, outname, row.names = FALSE, quote = FALSE)
}


library(ggsci)
niche_colors <- c(
  CN1 = "#CE3D32FF",
  CN2 = "#5DB1DDFF",
  CN3 = "#FF8C00FF",
  CN4 = "#696969FF",
  CN5 = "#749B58FF",
  CN6 = "#40E0D0FF",
  CN7 = "#466983FF"
)

pdf("0MSD_spatial_niche.pdf", height = 15, width = 12)
for (s in samples) {
  cells_to_plot <- rownames(MSD@meta.data)[MSD@meta.data$sample == s]
  p <- DimPlot(
    object = MSD,
    reduction = "spatial",
    pt.size = 2,
    label = FALSE,
    group.by = 'niche',
    cells = cells_to_plot,
    raster = FALSE
  ) +
    scale_color_manual(values = niche_colors, drop = FALSE) +  
    labs(x = "spatial1", y = "spatial2", title = s) +
    theme(
      panel.border = element_rect(fill = NA, color = "black", size = 1),
      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      legend.position = 'right',
      legend.text = element_text(color = "black", size = 16, face = "bold")
    )
  print(p)
}
dev.off()

test <- as.data.frame(table(MSD@meta.data$niche, MSD@meta.data$group))
colnames(test) <- c("niche", "group", "Freq")
pdf("1MSD_niche_per_percent_of_group.pdf", width = 6, height = 8)
ggplot(test, aes(x = group, y = Freq, fill = niche)) +
  geom_col(position = "fill") +
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0)),
    labels = scales::percent_format()
  ) +
  scale_fill_manual(values = niche_colors) +  
  labs(x = '', y = 'Percentage') +
  guides(fill = guide_legend(ncol = 4)) +
  theme_classic() +
  theme(
    panel.background = element_rect(fill = 'transparent', color = "white"),
    legend.position = 'top',
    legend.direction = "horizontal",
    legend.text = element_text(color = "black", size = 11, face = "bold"),
    legend.title = element_blank(),
    legend.justification = 0.9,
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", size = 30, face = "bold"),
    axis.text.y = element_text(color = "black", size = 30, face = "bold"),
    axis.title.y = element_text(color = "black", size = 25, face = "bold"),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(size = 1.5),
    axis.line.x = element_line(linetype = 1, color = "black", size = 1.5),
    axis.line.y = element_line(linetype = 1, color = "black", size = 1.5)
  )
dev.off()


test<- as.data.frame(table(MSD@meta.data$niche,MSD@meta.data$sample))
test$Var2 <- factor(test$Var2, levels=c("MM5_WHO1", "MM7_WHO1", "MM8_WHO1", "MM1_WHO2", "MM2_WHO2", "MM3_WHO2", "MM4_WHO3", "MM6_WHO3"), ordered=TRUE)
colnames(test) <- c("niche", "group", "Freq")
pdf("1MSD_niche_per_percent_of_sample.pdf",width = 8,height=8)
ggplot(test, aes(x = group, y = Freq, fill = niche)) +
  geom_col(position = "fill") +
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0)),
    labels = scales::percent_format()
  ) +
  scale_fill_manual(values = niche_colors) + 
  labs(x = '', y = 'Percentage') +
  guides(fill = guide_legend(ncol = 4)) +
  theme_classic() +
  theme(
    panel.background = element_rect(fill = 'transparent', color = "white"),
    legend.position = 'top',
    legend.direction = "horizontal",
    legend.text = element_text(color = "black", size = 11, face = "bold"),
    legend.title = element_blank(),
    legend.justification = 0.9,
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", size = 30, face = "bold"),
    axis.text.y = element_text(color = "black", size = 30, face = "bold"),
    axis.title.y = element_text(color = "black", size = 25, face = "bold"),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(size = 1.5),
    axis.line.x = element_line(linetype = 1, color = "black", size = 1.5),
    axis.line.y = element_line(linetype = 1, color = "black", size = 1.5)
  )
dev.off()


#niche ----------------
library(dplyr)
library(tidyr)
library(ggplot2)

tab <- table(MSD@meta.data$cell_labels_updated, MSD@meta.data$niche)
tab <- as.data.frame(tab)
colnames(tab) <- c("cell_labels_updated", "niche", "count")

tab <- tab %>%
  group_by(niche) %>%
  mutate(percent = count / sum(count) * 100)

cell_colors <- c(
  "Tumor_Multipotent"    = "#FF0049",
  "Tumor_Oligopotent"    = "#E6553A",
  "Tumor_Unipotent"      = "#FCAE91",
  "Tumor_Differentiated" = "#FDD0A2",
  "T_Cells"              = "#00F3D8",
  "Myeloid_cells"        = "#00ACEB",
  "Stromal_cells"        = "#87B400",
  "Necrotic_Cells"       = "#ABABAB"
)

pdf("1MSD_niche_cell_percent.pdf", width = 6, height = 8)
ggplot(tab, aes(x = niche, weight = count, fill = cell_type)) +  
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Niche", y = "Percentage", fill = "Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = cell_colors)
dev.off()

tab <- table(MSD@meta.data$celltype, MSD@meta.data$niche)
tab <- as.data.frame(tab)
colnames(tab) <- c("celltype", "niche", "count")
tab <- tab %>%
  group_by(niche) %>%
  mutate(percent = count / sum(count) * 100)



test<- as.data.frame(table(malignant@meta.data$CytoTRACE2_Potency,malignant@meta.data$group))
#test$Var2 <- factor(test$Var2, levels=c("MM5_WHO1", "MM7_WHO1", "MM8_WHO1", "MM1_WHO2", "MM2_WHO2", "MM3_WHO2", "MM4_WHO3", "MM6_WHO3"), ordered=TRUE)
pdf("1malignant_CytoTRACE2_Potency_per_percent_of_sample.pdf",width = 4,height=8)
ggplot(test,aes(x= Var2,y=Freq,fill= Var1))+
  geom_col( position="fill")+scale_y_continuous(expand = expansion(mult = c(0.02,0)),labels = scales::percent_format())+
  scale_x_discrete(labels= c("MM5_WHO1", "MM7_WHO1", "MM8_WHO1", "MM1_WHO2", "MM2_WHO2", "MM3_WHO2", "MM4_WHO3", "MM6_WHO3"))+ #X轴坐标自定义
  labs(x='',y='Percentage')+guides(fill=guide_legend(ncol=2))+#图例分3列
  theme_classic()+
  theme(panel.background=element_rect(fill='transparent',color ="white"),
        legend.position = 'top',legend.direction = "horizontal",legend.text=element_text(color = "black",size=11,face="bold"),
        legend.title=element_blank(),legend.justification = 0.9,#position是‘top’时进行左右微调，适应字体大小
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,color = "black",size=18,face="bold"),#X轴字体
        axis.text.y = element_text(color = "black",size=15,face="bold"), #Y轴字体
        axis.title.y = element_text(color = "black",size=20,face="bold"),#Percentage的大小
        axis.ticks.x = element_blank(),axis.ticks.y = element_line(size=1.5),#XY轴刻度调整
        axis.line.x=element_line(linetype=1,color="black",size=1.5),axis.line.y=element_line(linetype=1,color="black",size=1.5))+ #XY轴线调整
  scale_fill_igv()
dev.off()

df <- data.frame(
  group = c("WHO1", "WHO2", "WHO3"),
  multipotent = c(1381, 3894, 24113),
  total = c(260413+1381+87+1, 332087+3894+224+9, 177602+24113+886+20)
)
df$others <- df$total - df$multipotent
table_matrix <- as.matrix(df[, c("multipotent", "others")])
rownames(table_matrix) <- df$group
chisq.test(table_matrix)


library(dplyr)
library(tidyr)
library(Seurat)
library(cowplot)
library(stringr)
library(ggplot2)


niche_colors <- c(
  CN1 = "#CE3D32FF",
  CN2 = "#5DB1DDFF",
  CN3 = "#FF8C00FF",
  CN4 = "#696969FF",
  CN5 = "#749B58FF",
  CN6 = "#40E0D0FF",
  CN7 = "#466983FF"
)

#group
df_group<-as.data.frame(table(MSD$niche,MSD$group))
colnames(df_group)<-c("niche","group","counts")
df_group_prop <- df_group %>%
  group_by(group) %>%  
  mutate(
    total_by_group = sum(counts), 
    proportion = counts / total_by_group  
  ) %>%
  ungroup() %>%
replace_na(list(proportion = 0))

p_group <- ggplot(df_group_prop, aes(x = group, y = niche, color = niche, size = proportion)) +
    geom_point(alpha = 0.7) +
    labs(
      x = "group",
      y = "Niche",
      color = "Niche",
      size = "Proportion",
      title = "group"
    ) +
    scale_size_continuous(range = c(0, 10)) +
    scale_color_manual(values = niche_colors[names(niche_colors) %in% unique(df_group$niche)]) +
    theme_classic() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 9),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      legend.position = "right"
    )

ggsave(
  filename = "1niche_WHO_group_dot.pdf",
  plot = p_group,
  width = 6,
  height = 4,
  bg = "white"
)
