library(ggplot2)
library(ggsci)
library(Seurat)

MSD<-readRDS('/public/home/N202410/naomoliu/01.snRNA_celltype_result/celltype_v5.rds')


MSD@meta.data$sample <- as.character(MSD@meta.data$sample)
MSD@meta.data$group[MSD@meta.data$sample == 'MM3'] <- 'WHO3'



pdf("0MSD_QC.pdf",height=4,width=6)
VlnPlot(MSD, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size=0, ncol = 3)
dev.off()


MSD <- NormalizeData(MSD, normalization.method = "LogNormalize", scale.factor = 10000)
MSD <- FindVariableFeatures(object = MSD,selection.method = "vst", nfeatures = 2000)
MSD <- ScaleData(object = MSD, vars.to.regress = c("nFeature_RNA", "nCount_RNA",'percent.mt'), model.use = "linear")
MSD <- RunPCA(object = MSD, pc.genes = VariableFeatures(MSD), do.print = TRUE, pcs.print = 1:5, genes.print = 5)
pdf("0MSD_PCA.pdf",height=20,width=10) 
dev.off()
MSD <- FindNeighbors(MSD, dims = 1:20) 
MSD <- FindClusters(MSD, resolution = c(0.1,0.2,0.3))
MSD <- RunUMAP(object = MSD,dims = 1:20, verbose = FALSE)


pdf("0MSD_SN_sample.pdf")
DimPlot(object = MSD, reduction = "umap",label = T,group.by='sample')+scale_color_igv()+#NoLegend()+
        labs(x = "UMAP1", y = "UMAP2",title = "") + 
        theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
              axis.text.y = element_blank(),axis.ticks.y = element_blank(),
              axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position='right',
              legend.text=element_text(color = "black",size=16,face="bold"))
dev.off()


pdf("0MSD_RNA_snn_res.0.3.pdf") 
DimPlot(object = MSD, reduction = "umap",label = T,group.by='RNA_snn_res.0.3')+scale_color_igv()+#NoLegend()+
        labs(x = "UMAP1", y = "UMAP2",title = "") + 
        theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
              axis.text.y = element_blank(),axis.ticks.y = element_blank(),
              axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position='right',
              legend.text=element_text(color = "black",size=16,face="bold"))
dev.off()
pdf("0MSD_QC.pdf",height=4,width=6)
VlnPlot(MSD, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),group.by='RNA_snn_res.0.3', pt.size=0, ncol = 3)
dev.off()


library(scCustomize)
pdf("0MSD_markers_scCustomize.pdf")
Plot_Density_Custom(seurat_object = MSD, features = "COL4A1", custom_palette = c("lightgray", "#4169E1")) + NoLegend() + NoAxes()
Plot_Density_Custom(seurat_object = MSD, features = "CD8A", custom_palette = c("lightgray", "#4169E1"))+ NoLegend() + NoAxes()
Plot_Density_Custom(seurat_object = MSD, features = "CD14", custom_palette = c("lightgray", "#4169E1")) + NoLegend() + NoAxes()
Plot_Density_Custom(seurat_object = MSD, features = "CDH1", custom_palette = c("lightgray", "#4169E1")) + NoLegend() + NoAxes()
Plot_Density_Custom(seurat_object = MSD, features = "PTPRC", custom_palette = c("lightgray", "#4169E1")) + NoLegend() + NoAxes()
dev.off()

Idents(object = MSD) <- "RNA_snn_res.0.3" 
MSD <- RenameIdents(MSD, 
    `0` = "Tumor cells", 
    `1` = "Tumor cells", 
    `2` = "Tumor cells", 
    `3` = "Tumor cells", 
    `4` = "Macrophage", 
    `5` = "Tumor cells", 
    `6` = "Macrophage", 
    `7` = "Stromal", 
    `8` = "Tumor cells", 
    `9` = "Tumor cells", 
    `10` = "T cells", 
    `11` = "Tumor cells",
    `12` = "Macrophage", 
    `13` = "Tumor cells", 
    `14` = "Tumor cells"
)
MSD@meta.data$cell_labels<-MSD@active.ident

my_colors <- c(
  "Tumor cells"    = "#E64B35FF",   
  "Macrophage"     = "#3C5488FF",  
  "Stromal"        = "#00A087FF",   
  "T cells"        = "#F39B7FFF"   
)
pdf("0MSD_RNA_snn_res.anno.pdf") 
DimPlot(
  object = MSD, 
  reduction = "umap", 
  label = FALSE
) + 
  scale_color_manual(values = my_colors) +
  labs(x = "UMAP1", y = "UMAP2", title = "") + 
  theme(
    panel.border = element_rect(fill = NA, color = "black", size = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.text = element_text(color = "black", size = 16, face = "bold")
  ) +
  NoLegend()
dev.off()

library(dplyr)

allmarkers_res.0.1 <- FindAllMarkers(MSD, logfc.threshold =1, test.use = "wilcox",assay = "sketch")
allmarkers_res.0.1_arrange <-  allmarkers_res.0.1 %>% dplyr::arrange(cluster, p_val_adj)
allmarkers_res.0.1_top5 <- allmarkers_res.0.1_arrange %>% group_by(cluster) %>% top_n(5, avg_log2FC)
 
malignant<-subset(MSD,idents=c(0:3,5))#tumor
malignant <- NormalizeData(malignant, normalization.method = "LogNormalize", scale.factor = 10000)
malignant <- FindVariableFeatures(object = malignant,selection.method = "vst", nfeatures = 2000)
malignant <- ScaleData(object = malignant, vars.to.regress = c("nFeature_RNA", "nCount_RNA",'percent.mt'), model.use = "linear")
malignant <- RunPCA(object = malignant, pc.genes = VariableFeatures(malignant), do.print = TRUE, pcs.print = 1:5, genes.print = 5)
pdf("1malignant_PCA.pdf",height=20,width=10) 
DimHeatmap(malignant, dims = 1:35, cells = 500, balanced = TRUE)
dev.off()
malignant <- FindNeighbors(malignant, dims = 20) 
malignant <- FindClusters(malignant, resolution = c(0.1))
malignant <- RunUMAP(object = malignant,dims = 1:20, verbose = FALSE)

pdf("1malignant_RNA_snn_res.0.1.pdf")
DimPlot(object = malignant, reduction = "umap",label = T,group.by='sample')+scale_color_igv()+NoLegend()+
        labs(x = "UMAP1", y = "UMAP2",title = "") + 
        theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
              axis.text.y = element_blank(),axis.ticks.y = element_blank(),
              axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position='right',
              legend.text=element_text(color = "black",size=16,face="bold"))
dev.off()

malignant@meta.data$sample <- as.character(malignant@meta.data$sample)
malignant@meta.data$group[malignant@meta.data$sample == 'MM3'] <- 'WHO3'

library(scCustomize)
pdf("1malignant_dedifferentation_genes.pdf",height=10,width=10)
Plot_Density_Custom(seurat_object = malignant, features = "SLC2A1", custom_palette = c("lightgray", "#DA0A16")) + NoLegend() + NoAxes()
Plot_Density_Custom(seurat_object = malignant, features = "COL18A1", custom_palette = c("lightgray", "#DA0A16"))+ NoLegend() + NoAxes()

dev.off()


DefaultAssay(MSD)<-'RNA'
MSD@reductions$PCA <- NULL
MSD@assays$RNA@scale.data <- matrix()







saveRDS(MSD,'0MSD.rds')

MSD<-readRDS('0MSD.rds')
MSD@meta.data$group[MSD@meta.data$sample == 'MM3'] <- 'WHO3'


pdf("0MSD_immune_markers.pdf")#marker展示
FeaturePlot(MSD,reduction = "umap", features = c('PDCD1'),cols=c("lightgrey", "#DA0A16"))+NoLegend()+NoAxes()
FeaturePlot(MSD,reduction = "umap", features = c('CD274'),cols=c("lightgrey", "#DA0A16"))+NoLegend()+NoAxes()
FeaturePlot(MSD,reduction = "umap", features = c('CTLA4'),cols=c("lightgrey", "#DA0A16"))+NoLegend()+NoAxes()
FeaturePlot(MSD,reduction = "umap", features = c('LMTK3'),cols=c("lightgrey", "red"))+NoLegend()
FeaturePlot(MSD,reduction = "umap", features = c('LAG3'),cols=c("lightgrey", "red"))+NoAxes()
FeaturePlot(MSD,reduction = "umap", features = c('HAVCR2'),cols=c("lightgrey", "red"))+NoAxes() ***
FeaturePlot(MSD,reduction = "umap", features = c('TIGIT'),cols=c("lightgrey", "red"))+NoAxes()
FeaturePlot(MSD,reduction = "umap", features = c('IDO1'),cols=c("lightgrey", "red"))+NoAxes()
dev.off()



library(ggsci)
library(monocle)
library(RColorBrewer)
library(Seurat)
library(dplyr)
library(msigdbr)
genesets <- msigdbr(species = "Homo sapiens", category = "H") 
genesets <- subset(genesets, select = c("gs_name","gene_symbol")) %>% as.data.frame()
genesets <- split(genesets$gene_symbol, genesets$gs_name)

malignant <- UCell::AddModuleScore_UCell(malignant, features = list(c('SLC2A1','COL18A1','MCTP2','CCN4','NREP')), ncores=8, name = "dedifferentation")
malignant <- AddModuleScore(malignant, features = list(c('SLC2A1','COL18A1','MCTP2','CCN4','NREP')), name = "dedifferentation2")

malignant <- UCell::AddModuleScore_UCell(malignant, features = list(c('COL18A1')), ncores=8, name = "COL18A1")
malignant <- AddModuleScore(malignant, features = list(c('COL18A1')), name = "COL18A12")


Mono_matrix <- GetAssayData(malignant, slot = "counts")
feature_ann<-data.frame(gene_id=rownames(Mono_matrix),gene_short_name=rownames(Mono_matrix))
rownames(feature_ann)<-rownames(Mono_matrix)

Mono_fd<-new("AnnotatedDataFrame", data = feature_ann)
sample_ann <- malignant@meta.data[colnames(GetAssayData(malignant, assay = "RNA")), ]
#sample_ann<-sample_ann[,-which(names(sample_ann)=='sample_name')]#不然会报错
rownames(sample_ann)<-colnames(Mono_matrix)
Mono_pd<-new("AnnotatedDataFrame", data =sample_ann)
#build new cell data set
Mono.cds<-newCellDataSet(Mono_matrix,phenoData =Mono_pd,featureData =Mono_fd,expressionFamily=negbinomial.size())


malignant<-readRDS('1malignant_CytoTRACE2.rds')
Idents(malignant)<-'niche'
ordering_genes <- FindAllMarkers(malignant, logfc.threshold = 1, test.use = "wilcox", min.pct = 0,padj=0.05, only.pos = TRUE,assay = "RNA")
pdf("./0SN_result/1malignant_dedifferentiation.pdf",height=5,width=5)
VlnPlot(malignant, features = c("signature_1dedifferentation"),group.by='group', pt.size=0)+scale_fill_npg()
dev.off()


head(pData(Mono.cds))
head(fData(Mono.cds))
#预处理
Mono.cds <- estimateSizeFactors(Mono.cds)
Mono.cds <- estimateDispersions(Mono.cds)
Mono.cds <- detectGenes(Mono.cds, min_expr = 0.1)


Mono.cds <- setOrderingFilter(Mono.cds, ordering_genes$gene)
#DDRtree 
Mono.cds <- reduceDimension(Mono.cds,method = 'DDRTree')

#计算psudotime值
Mono.cds <- orderCells(Mono.cds)
head(pData(Mono.cds))

Mono.cds=orderCells(Mono.cds,root_state = 3)


pdf("./0SN_result/1malignant_signature_1dedifferentation.pdf", width=5, height=5)
plot_cell_trajectory(Mono.cds, color_by = "signature_1dedifferentation", cell_size = 0.3) +
  scale_color_gradientn(colours = c("lightgrey", "#DA0A16"),limits = c(0, 0.7),breaks = c(0, 0.25, 0.5, 0.75)) +labs(title = 'signature_1dedifferentation') +theme_classic() +
  theme(
    panel.background = element_rect(fill = 'transparent', color = "white"),
    legend.position = c(0.95, 0.1),legend.direction = "vertical",legend.text = element_text(color = "black", size = 8, face = "bold"),legend.title = element_blank(),
    plot.title = element_text(face = 'bold', size = 12, hjust = 0.5),axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(size = 1),axis.line.x = element_line(linetype = 1, color = "black", size = 0.6),
    axis.line.y = element_line(linetype = 1, color = "black", size = 0.6)
  )
plot_cell_trajectory(Mono.cds, color_by = "HALLMARK_G2M_CHECKPOINT", cell_size = 0.3) +
  scale_color_gradientn(
    colours = c("lightgrey", "lightgrey", "lightgrey", "#DA0A16"),values = scales::rescale(c(0, 0.08, 0.10, 0.22)), limits = c(0, 0.22)) +labs(title = 'HALLMARK_G2M_CHECKPOINT') +theme_classic() +
  theme(
    panel.background = element_rect(fill = 'transparent', color = "white"),
    legend.position = c(0.95, 0.1),legend.direction = "vertical",
    legend.text = element_text(color = "black", size = 8, face = "bold"),legend.title = element_blank(),
    plot.title = element_text(face = 'bold', size = 12, hjust = 0.5),axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(size = 1),axis.line.x = element_line(linetype = 1, color = "black", size = 0.6),
    axis.line.y = element_line(linetype = 1, color = "black", size = 0.6)
  )
plot_cell_trajectory(Mono.cds, color_by = "HALLMARK_MYC_TARGETS_V1", cell_size = 0.3) +
  scale_color_gradientn(colours = c("lightgrey", "lightgrey", "lightgrey", "#DA0A16"),values = scales::rescale(c(0, 0.08, 0.12, 0.285)),limits = c(0, 0.285) ,) +labs(title = 'HALLMARK_MYC_TARGETS_V1') +theme_classic() +
  theme(
    panel.background = element_rect(fill = 'transparent', color = "white"),
    legend.position = c(0.95, 0.1),legend.direction = "vertical",
    legend.text = element_text(color = "black", size = 8, face = "bold"),legend.title = element_blank(),
    plot.title = element_text(face = 'bold', size = 12, hjust = 0.5),axis.ticks.x = element_blank(),axis.ticks.y = element_line(size = 1),
    axis.line.x = element_line(linetype = 1, color = "black", size = 0.6),axis.line.y = element_line(linetype = 1, color = "black", size = 0.6)
  )
plot_cell_trajectory(Mono.cds, color_by = "HALLMARK_MYC_TARGETS_V2", cell_size = 0.3) +
  scale_color_gradientn(
  colours = c("lightgrey", "lightgrey", "#DA0A16"),values = scales::rescale(c(0, 0.05, 0.16)),limits = c(0, 0.16) ,) +labs(title = 'HALLMARK_MYC_TARGETS_V2') +theme_classic() +
  theme(
    panel.background = element_rect(fill = 'transparent', color = "white"),
    legend.position = c(0.95, 0.1),legend.direction = "vertical",
    legend.text = element_text(color = "black", size = 8, face = "bold"),legend.title = element_blank(),
    plot.title = element_text(face = 'bold', size = 12, hjust = 0.5),axis.ticks.x = element_blank(),axis.ticks.y = element_line(size = 1),
    axis.line.x = element_line(linetype = 1, color = "black", size = 0.6),axis.line.y = element_line(linetype = 1, color = "black", size = 0.6)
  )
plot_cell_trajectory(Mono.cds, color_by = "HALLMARK_E2F_TARGETS", cell_size = 0.3) +
  scale_color_gradientn(colours = c("lightgrey", "lightgrey", "lightgrey", "#DA0A16"),values = scales::rescale(c(0, 0.08, 0.15, 0.36)),limits = c(0, 0.36) ,) +labs(title = 'HALLMARK_E2F_TARGETS') +theme_classic() +
  theme(
    panel.background = element_rect(fill = 'transparent', color = "white"),
    legend.position = c(0.95, 0.1),legend.direction = "vertical",
    legend.text = element_text(color = "black", size = 8, face = "bold"),legend.title = element_blank(),
    plot.title = element_text(face = 'bold', size = 12, hjust = 0.5),axis.ticks.x = element_blank(),axis.ticks.y = element_line(size = 1),
    axis.line.x = element_line(linetype = 1, color = "black", size = 0.6),axis.line.y = element_line(linetype = 1, color = "black", size = 0.6)
  )
plot_cell_trajectory(Mono.cds, color_by = "HALLMARK_WNT_BETA_CATENIN_SIGNALING", cell_size = 0.3) +
  scale_color_gradientn(colours = c("lightgrey", "lightgrey", "lightgrey", "#DA0A16"),values = scales::rescale(c(0, 0.08, 0.065, 0.19)),limits = c(0, 0.19) ,) +labs(title = 'HALLMARK_WNT_BETA_CATENIN_SIGNALING') +theme_classic() +
  theme(
    panel.background = element_rect(fill = 'transparent', color = "white"),
    legend.position = c(0.95, 0.1),legend.direction = "vertical",
    legend.text = element_text(color = "black", size = 8, face = "bold"),legend.title = element_blank(),
    plot.title = element_text(face = 'bold', size = 12, hjust = 0.5),axis.ticks.x = element_blank(),axis.ticks.y = element_line(size = 1),
    axis.line.x = element_line(linetype = 1, color = "black", size = 0.6),axis.line.y = element_line(linetype = 1, color = "black", size = 0.6)
  )
dev.off()


pdf("./0SN_result/1malignant_Monocle_Pseudotime.pdf",width=5,height=5)
plot_cell_trajectory(Mono.cds, color_by = "Pseudotime",cell_size = 0.75)+
                     theme(axis.ticks = element_blank(),axis.text=element_blank(),axis.title = element_text(size=12,face="bold"),
                           axis.line.y=element_line(linetype=1,color="black",size=1.5),axis.line.x=element_line(linetype=1,color="black",size=1.5) )#+NoLegend()
plot_cell_trajectory(Mono.cds, color_by = "State",cell_size = 0.75)+
                     theme(axis.ticks = element_blank(),axis.text=element_blank(),axis.title = element_text(size=12,face="bold"),
                           axis.line.y=element_line(linetype=1,color="black",size=1.5),axis.line.x=element_line(linetype=1,color="black",size=1.5) )#+NoLegend()
dev.off()



BEAM_res=BEAM(Mono.cds,branch_point = 1,cores = 8, progenitor_method = "duplicate")
BEAM_res=BEAM_res[,c("gene_short_name","pval","qval")]
save.image('./0SN_result/1malignant_monocle2.rdata')#
tmp1=plot_genes_branched_heatmap(Mono.cds[row.names(subset(BEAM_res,qval<1e-4)),],
                                 branch_point = 1,
                                 num_clusters = 5,
                                 cores = 16,
                                 branch_labels = c("trace1", "trace2"),

                                 use_gene_short_name = T,
                                 show_rownames = F,
                                 return_heatmap = T 
                                 )
pdf("./0SN_result/1malignant_BEAMheatmap.pdf",width = 5,height = 6)
tmp1$ph_res
dev.off()



#library(DOSE)
#library(GO.db)
library(org.Hs.eg.db)
library(clusterProfiler)
anno<-tmp1$annotation_row
anno$gene=rownames(anno)
#循环GO富集
for (i in 1:5){
    c1<-subset(anno,anno$Cluster== i)
    c1 = bitr(c1$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    ego_ALL <- enrichGO(gene = c1$ENTREZID,
                        OrgDb= org.Hs.eg.db,
                        ont = "BP",#MF\BP\CC\All
                        pAdjustMethod = "BH",
                        minGSSize = 1,
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05,
                        readable = TRUE)
names=paste0('./0SN_result/01malignant_BEAMheatmap_GO',i,'.pdf')
ggsave(names,barplot(ego_ALL,showCategory=20,title="Enrichment GO Top30")+theme(axis.text.y = element_text(size = 4)))
}





