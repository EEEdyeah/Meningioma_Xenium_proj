library(ggplot2)
library(ggsci)
library(Seurat)
setwd('/public/home/N202302/scRNA/MSD')
MSD<-readRDS('/public/home/N202505/project/naomoliu/result/res1/noIntegration.rds')
MSD <- NormalizeData(MSD, normalization.method = "LogNormalize", scale.factor = 10000)

MSD <- ScaleData(object = MSD,features = rownames(MSD) , vars.to.regress = c("nFeature_RNA", "nCount_RNA"), model.use = "linear")
MSD <- RunPCA(object = MSD, pc.genes = rownames(MSD), do.print = TRUE, pcs.print = 1:5, genes.print = 5)
pdf("0MSD_PCA.pdf",height=20,width=10) #UMAP聚类图，0.1分辨率
DimHeatmap(MSD, dims = 1:35, cells = 500, balanced = TRUE)
dev.off()
MSD <- FindNeighbors(MSD, dims = 1:20) #options(future.globals.maxSize = 1000 * 1024^2)
MSD <- FindClusters(MSD, resolution = c(0.5,0.6,0.8,1))
MSD <- RunUMAP(object = MSD,dims = 1:20, verbose = FALSE)

pdf("0MSD_sample.pdf",width=8,height=6) #UMAP聚类图，0.1分辨率
DimPlot(object = MSD, reduction = "umap",label = F,group.by='sample')+scale_color_igv()+#NoLegend()+
        labs(x = "UMAP1", y = "UMAP2",title = "") + 
        theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
              axis.text.y = element_blank(),axis.ticks.y = element_blank(),
              axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position='right',
              legend.text=element_text(color = "black",size=16,face="bold"))
dev.off()



pdf("0MSD_RNA_snn_res.0.6.spatial.pdf",height=20,width=12) 
DimPlot(object = MSD, reduction = "spatial",label = T,raster=FALSE,group.by='RNA_snn_res.0.6')+scale_color_igv()+#NoLegend()+
        labs(x = "spatial1", y = "spatial2",title = "") + 
        theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
              axis.text.y = element_blank(),axis.ticks.y = element_blank(),
              axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position='right',
              legend.text=element_text(color = "black",size=16,face="bold"))
dev.off()




Idents(object = MSD) <- "RNA_snn_res.0.6" 
MSD <- RenameIdents(MSD, 
    `0` = "Tumor_Cells_MM2_WHO2", 
    `1` = "Tumor_Cells_MM3_WHO2", 
    `2` = "Tumor_Cells_MM7_WHO1", 
    `3` = "Tumor_Cells_MM8_WHO1", 
    `4` = "Myeloid_cells", 
    `5` = "Tumor_Cells_MM4_WHO3", 
    `6` = "Tumor_Cells_MM6_WHO3", 
    `7` = "Myeloid_cells", 
    `8` = "Stromal_cells",
    `9` = "Tumor_Cells_MM1_WHO2", 
    `10` = "Necrotic_Cells", 
    `11` = "Stromal_cells", 
    `12` = "Tumor_Cells_MM5_WHO3", 
    `13` = "Myeloid_cells", 
    `14` = "Myeloid_cells", 
    `15` = "Stromal_cells", 
    `16` = "T_Cells", 
    `17` = "Tumor_Cells_MM3_WHO2", 
    `18` = "Tumor_Cells_MM4_WHO3", 
    `19` = "Myeloid_cells", 
    `20` = "Necrotic_Cells", 
    `21` = "Necrotic_Cells", 
    `22` = "Necrotic_Cells"
)
pdf("0MSD_RNA_snn_res.0.6.anno.pdf",height=10,width=12) 
DimPlot(object = MSD, reduction = "umap",label = F,raster=FALSE)+scale_color_npg()+
        labs(x = "UMAP1", y = "UMAP2",title = "") + 
        theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
              axis.text.y = element_blank(),axis.ticks.y = element_blank(),
              axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position='right',
              legend.text=element_text(color = "black",size=16,face="bold"))
dev.off()
MSD@meta.data$cell_labels<-MSD@active.ident



library(ggplot2)
library(ggsci)
library(Seurat)
library(patchwork)
library(dplyr)
setwd('/public/home/N202302/scRNA/MSD')

#------------------------tumor cells-----------------------------
malignant<-subset(MSD,idents='Tumor_Cells')
malignant <- NormalizeData(object = malignant, normalization.method = "LogNormalize", scale.factor = 10000)
malignant <- ScaleData(object = malignant,features = rownames(MSD) , vars.to.regress = c("nFeature_RNA", "nCount_RNA"), model.use = "linear")
malignant <- RunPCA(object = malignant, pc.genes = rownames(MSD), do.print = TRUE, pcs.print = 1:5, genes.print = 5)
pdf("1malignant_PCA.pdf",height=20,width=10) 
DimHeatmap(malignant, dims = 1:35, cells = 500, balanced = TRUE)
dev.off()
malignant <- FindNeighbors(malignant, dims = 1:20) 
malignant <- FindClusters(malignant, resolution = c(0.5))
malignant <- RunUMAP(object = malignant,dims = 1:20, verbose = FALSE)

library(CytoTRACE2)
Idents(malignant)<-'sample'
malignant <- cytotrace2(malignant,is_seurat = TRUE,slot_type = "count",
                        species = "human",seed = 970619,ncores=64)

annotation <- data.frame(phenotype=malignant@meta.data$sample) %>% set_rownames(., colnames(malignant))
plots <- plotData(cytotrace2_result=malignant, annotation=annotation, is_seurat=TRUE)
p1 <- plots$CytoTRACE2_UMAP
p2 <- plots$CytoTRACE2_Potency_UMAP
p3 <- plots$CytoTRACE2_Relative_UMAP
p4 <- plots$CytoTRACE2_Boxplot_byPheno

pdf("1malignant_CytoTRACE2.pdf",width = 12,height = 10)
(p1+p2+p3+p4) + plot_layout(ncol = 2)
dev.off()
malignant@assays$RNA@scale.data <- matrix()


#Tumor cell annotation
MSD@meta.data$cell_labels_updated <- as.character(MSD@meta.data$cell_labels)
MSD@meta.data$cell_labels_updated[ MSD@meta.data$cell_labels_updated %in% c('Tumor_Cells')]="Tumor_Differentiated"
MSD@meta.data[rownames(malignant@meta.data[malignant@meta.data$CytoTRACE2_Potency %in% c('Multipotent'),]),]$cell_labels_updated="Tumor_Multipotent"
MSD@meta.data[rownames(malignant@meta.data[malignant@meta.data$CytoTRACE2_Potency %in% c('Oligopotent'),]),]$cell_labels_updated="Tumor_Oligopotent"
MSD@meta.data[rownames(malignant@meta.data[malignant@meta.data$CytoTRACE2_Potency %in% c('Unipotent'),]),]$cell_labels_updated="Tumor_Unipotent"



#niche----------------------------------------------
set.seed(970619)
group.by<-'cell_labels_updated'
object<-MSD
niches.k = 7
neighbors.k = 25
cluster.name = "niches"
cells <- Cells(object)
group.labels <- unlist(object[[group.by]][cells, ])
groups <- sort(unique(group.labels))
cell.type.mtx <- matrix(data = 0, nrow = length(cells), ncol = length(groups))
rownames(cell.type.mtx) <- cells
colnames(cell.type.mtx) <- groups
cells.idx <- seq_along(cells)
group.idx <- match(group.labels, groups)
cell.type.mtx[cbind(cells.idx, group.idx)] <- 1

coords <-MSD@reductions$spatial@cell.embeddings
colnames(coords)<-c('x','y')

coords <- as.matrix(coords[, c("x", "y")])
neighbors <- FindNeighbors(coords, k.param = neighbors.k,compute.SNN = FALSE)
sum.mtx <- as.matrix(neighbors[["nn"]] %*% cell.type.mtx)
niche.assay <- CreateAssayObject(counts = t(sum.mtx))

object <- CreateSeuratObject(counts = niche.assay, assay = 'niches')
object <- ScaleData(object)
results <- kmeans(x = t(object[['niches']]@scale.data), centers = niches.k,nstart = 30,iter.max = 100)
MSD[['niche']]<-paste0("CN",results[["cluster"]])


saveRDS(MSD,'0MSD_new.rds')

malignant@meta.data$niche <- MSD@meta.data[rownames(malignant@meta.data), "niche"]
saveRDS(malignant,'1malignant_CytoTRACE2.rds')


library(cluster)
set.seed(123)
gap_stat <- clusGap(data.sub, FUN = kmeans, K.max = 13, B = 30)
pdf("0MSD_kmeans.pdf", height = 15, width = 15)
plot(gap_stat)
dev.off()

saveRDS(MSD,'0MSD_new.rds')

##########################################################################################################################################
MSD@meta.data$sample <- as.character(MSD@meta.data$sample)
Stromal<-readRDS('stromal_celltype3.rds')
MSD@meta.data[rownames(Stromal@meta.data[Stromal@meta.data$celltype %in% c('VECs'),]),]$cell_labels_updated="CAFs"
MSD@meta.data[rownames(Stromal@meta.data[Stromal@meta.data$celltype %in% c('Fibroblasts'),]),]$cell_labels_updated="Fibroblasts"
MSD@meta.data[rownames(Stromal@meta.data[Stromal@meta.data$celltype %in% c('EndoMTs'),]),]$cell_labels_updated="CAFs"


test<- as.data.frame(table(MSD@meta.data$cell_labels_updated,MSD@meta.data$sample))
colnames(test) <- c("celltype", "sample", "count")
test$celltype <- gsub("^Tumor.*", "Tumor", test$celltype)
test <- test %>%
  mutate(celltype = gsub("^Tumor.*", "Tumor", celltype)) %>%
  group_by(celltype, sample) %>%
  summarise(count = sum(count), .groups = "drop")
test_prop <- test %>%group_by(sample) %>%mutate(total = sum(count),prop = count / total * 100)

test<- as.data.frame(table(MSD@meta.data$cell_labels_updated,MSD@meta.data$group))
colnames(test) <- c("celltype", "group", "count")
test$celltype <- gsub("^Tumor.*", "Tumor", test$celltype)
test <- test %>%
  mutate(celltype = gsub("^Tumor.*", "Tumor", celltype)) %>%
  group_by(celltype, group) %>%
  summarise(count = sum(count), .groups = "drop")
test_prop <- test %>%group_by(group) %>%mutate(total = sum(count),prop = count / total * 100)


library(ggplot2)
library(dplyr)
data <- data.frame(
  celltype = c(
    "CAFs","CAFs","CAFs","CAFs","CAFs","CAFs","CAFs","CAFs",
    "Myeloid cells","Myeloid cells","Myeloid cells",
    "Myeloid cells","Myeloid cells","Myeloid cells","Myeloid cells","Myeloid cells",
    "T cells","T cells","T cells","T cells","T cells","T cells","T cells","T cells",
    "Tumour","Tumour","Tumour","Tumour","Tumour","Tumour","Tumour","Tumour"
  ),
  sample = c(
    "MM5_WHO1","MM7_WHO1","MM8_WHO1",
    "MM1_WHO2","MM2_WHO2",
    "MM3_WHO3","MM4_WHO3","MM6_WHO3",
    
    "MM5_WHO1","MM7_WHO1","MM8_WHO1",
    "MM1_WHO2","MM2_WHO2",
    "MM3_WHO3","MM4_WHO3","MM6_WHO3",
    
    "MM5_WHO1","MM7_WHO1","MM8_WHO1",
    "MM1_WHO2","MM2_WHO2",
    "MM3_WHO3","MM4_WHO3","MM6_WHO3",
    
    "MM5_WHO1","MM7_WHO1","MM8_WHO1",
    "MM1_WHO2","MM2_WHO2",
    "MM3_WHO3","MM4_WHO3","MM6_WHO3"
  ),
  
  proportion = c(
    2.84844502, 1.54989776, 1.38465707,
    6.56473461, 3.16034745,
    4.71842170, 15.72551920, 3.54855598,
    
    35.46743197, 44.49637099, 12.99830236,
    24.27417011, 19.49839500,
    4.18763743, 4.58576168, 6.21499315,
    
    4.06211863, 1.78378663, 0.88974170,
    5.62874251, 0.97823789,
    0.79526020, 0.74834963, 1.63012687,
    
    54.89291510, 49.85168023, 59.35943021,
    60.05115982, 73.94642796,
    81.93256780, 76.06841666, 64.90601490
  ),
  
  WHO = c(
    "WHO1","WHO1","WHO1",
    "WHO2","WHO2",
    "WHO3","WHO3","WHO3",
    
    "WHO1","WHO1","WHO1",
    "WHO2","WHO2",
    "WHO3","WHO3","WHO3",
    
    "WHO1","WHO1","WHO1",
    "WHO2","WHO2",
    "WHO3","WHO3","WHO3",
    
    "WHO1","WHO1","WHO1",
    "WHO2","WHO2",
    "WHO3","WHO3","WHO3"
  )
)

data$WHO <- factor(
  data$WHO,
  levels = c("WHO1", "WHO2", "WHO3")
)

data$celltype <- factor(
  data$celltype,
  levels = c(
    "Tumour",
    "Myeloid cells",
    "CAFs",
    "T cells"
  )
)

summary_data <- data %>%
  group_by(celltype, WHO) %>%
  summarise(
    mean = mean(proportion),
    sd = sd(proportion),
    n = n(),
    sem = sd / sqrt(n),
    .groups = "drop"
  )

WHO_colors <- c(
  "WHO1" = "#5B8FF9",
  "WHO2" = "#61DDAA",
  "WHO3" = "#F08A8A"
)


p_data <- data %>%
  group_by(celltype) %>%
  summarise(
    p = kruskal.test(proportion ~ WHO)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_label = case_when(
      p < 0.001 ~ "P < 0.001",
      TRUE ~ paste0("P = ", sprintf("%.2f", p))
    )
  )

# 查看P值
print(p_data)
#  celltype           p p_label 
#1 Tumour        0.0622 P = 0.06
#2 Myeloid cells 0.0766 P = 0.08
#3 CAFs          0.0766 P = 0.08
#4 T cells       0.249  P = 0.25

p <- ggplot(
  data,
  aes(x = WHO,y = proportion
  )
) +
  # Individual samples
  geom_jitter(
    aes(fill = WHO),
    width = 0.10,
    size = 4,
    shape = 21,
    colour = "white",
    stroke = 0.8,
    alpha = 0.9
  ) +
  # Mean ± SEM
  geom_errorbar(
    data = summary_data,
    aes(
      x = WHO,
      y = mean,
      ymin = mean - sem,
      ymax = mean + sem
    ),
    inherit.aes = FALSE,
    width = 0.10,
    linewidth = 0.9,
    colour = "#333333"
  ) +
  # Mean
  geom_point(
    data = summary_data,
    aes(
      x = WHO,
      y = mean
    ),
    inherit.aes = FALSE,
    shape = 21,
    size = 4,
    fill = "#333333",
    colour = "white",
    stroke = 0.8
  ) +
  # Colours
  scale_fill_manual(values = WHO_colors) +
  # Four panels
  facet_wrap(~ celltype,ncol = 2,scales = "free_y") 
  labs(x = NULL,y = "Cell proportion (%)") +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      size = 12,face = "bold",colour = "#222222"
    ),
    axis.text.x = element_text(
      size = 10,colour = "#222222"
    ),
    axis.text.y = element_text(
      size = 9,colour = "#222222"
    ),
    axis.title.y = element_text(
      size = 11,face = "bold"
    ),
    axis.line = element_line(
      linewidth = 0.7,colour = "#222222"
    ),
    axis.ticks = element_line(
      linewidth = 0.6,colour = "#222222"
    ),
    panel.grid.major.y = element_line(
      linewidth = 0.3,colour = "grey90"
    ),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1.2,"lines"),
    legend.position = "none",
    panel.background = element_rect(fill = "white",colour = NA),
    plot.background = element_rect(fill = "white",colour = NA
    )
  )
