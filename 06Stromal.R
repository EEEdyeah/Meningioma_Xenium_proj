
#Stromal------------------------------------------
Stromal<-subset(MSD,idents='Stromal_cells')

Stromal <- NormalizeData(object = Stromal, normalization.method = "LogNormalize", scale.factor = 10000)
Stromal <- ScaleData(object = Stromal,features = rownames(MSD) , vars.to.regress = c("nFeature_RNA", "nCount_RNA"), model.use = "linear")
Stromal <- RunPCA(object = Stromal, pc.genes = rownames(MSD), do.print = TRUE, pcs.print = 1:5, genes.print = 5)
pdf("3Stromal_PCA.pdf",height=20,width=10) #UMAP聚类图，0.1分辨率
DimHeatmap(Stromal, dims = 1:35, cells = 500, balanced = TRUE)
dev.off()
Stromal <- FindNeighbors(Stromal, dims = 1:20) #options(future.globals.maxSize = 1000 * 1024^2)
Stromal <- FindClusters(Stromal, resolution = c(0.2,0.3))
Stromal <- RunUMAP(object = Stromal,dims = 1:20, verbose = FALSE)


Idents(object = Stromal) <- "RNA_snn_res.0.2"
pdf("2Stromal_Dimplot_0.2.pdf",width=8,height=5)
DimPlot(object = Stromal, reduction = "umap",label = T)+scale_color_npg()+#NoLegend()+
       labs(x = "UMAP1", y = "UMAP2",title = "") + 
       theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
             axis.text.y = element_blank(),axis.ticks.y = element_blank(),
             axis.text.x = element_blank(),axis.ticks.x = element_blank())
dev.off()
Stromal@assays$RNA@scale.data <- matrix()

pdf("2Stromal_Dimplot_sample.pdf",width=8,height=5)
DimPlot(object = Stromal, reduction = "umap",label = T,group.by='sample')+scale_color_npg()+#NoLegend()+
       labs(x = "UMAP1", y = "UMAP2",title = "") + 
       theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
             axis.text.y = element_blank(),axis.ticks.y = element_blank(),
             axis.text.x = element_blank(),axis.ticks.x = element_blank())
dev.off()

saveRDS('stromal_celltype3.rds')

Stromal<-readRDS('stromal_celltype3.rds')
Stromal@meta.data$sample <- as.character(Stromal@meta.data$sample)
Stromal@meta.data$group[Stromal@meta.data$sample == 'MM3_WHO2'] <- 'WHO3'
Stromal@meta.data$sample[Stromal@meta.data$sample == 'MM3_WHO2'] <- 'MM3_WHO3'

VECs<-subset(Stromal,idents=c('VECs' ))
EndoMTs<-subset(Stromal,idents=c('EndoMTs'))


Idents(Stromal) <- "niche" 
Stromal<-subset(Stromal,idents=c('CN7','CN6','CN3','CN1'))
Stromal <- AddModuleScore(Stromal,list(c('ABCA1','ACVRL1','ADRA1D','AHR','APEX1','CA9','CBFB','CD1D','CD44','CDH3','CDK4','CDK6','CDKN2B','CREB1','CSNK2A1','CSNK2A2','CYP1B1','DNMT3A','EDNRB','F2R','GJA1','HMGA1','HOXB7','ICAM1','IMPDH2','ITGA5','LDHA','LOXL2','LPL','MDM2','CD99','MME','MMP11','MMP14','ABCC1','MYC','NFKB1','NFKB2','P4HA1','SERPINE1','PDGFB','PTGS2','PTHLH','RELA','CXCL5','SKP2','SNAP25','SOX9','THBD','TFPI2','HMGA2','IRS2','RIPK2','PROM1','SQSTM1','MTA1','BAG2','GDA','PSIP1','RASSF1','BACE2','MCTS1','SIGIRR','COL18A1')),nbin=25,ctrl=100,name='SP1')
Stromal <- AddModuleScore(Stromal,list(c("BIRC5","CCNB1","CDK1","CHEK1","E2F1","FOXM1","MYBL2","RRM2","TOP2A","AURKB","PEG10","UHRF1","ASF1B","H19")),nbin=25,ctrl=100,name='E2F1')

Idents(Stromal)<-'celltype'
pdf("./2Stromal_regulators.pdf",height=5,width=5)
VlnPlot(Stromal, features = c("SP11"), pt.size=0)+ggsci::scale_fill_npg()
VlnPlot(Stromal, features = c("COL18A1"), pt.size=0)+ggsci::scale_fill_npg()
VlnPlot(Stromal, features = c("COL18A1"),group.by='niche', pt.size=0)+ggsci::scale_fill_npg()
VlnPlot(Stromal, features = c("SP11"),group.by='niche', pt.size=0)+ggsci::scale_fill_npg()
dev.off()

table(Stromal@meta.data$group,Stromal@active.ident)




library(ggplot2)
library(ggalluvial)
library(reshape2)
library(scales)
library(dplyr)



data <- data.frame(
  CellType = c("EndoMTs", "Fibroblasts", "VECs"),
  WHO1 = c(3315, 44484, 4845),
  WHO2 = c(3715, 3131, 8360),
  WHO3 = c(13425, 15737, 22426)
)


df_long <- melt(data, id.vars="CellType", variable.name="Sample", value.name="Count")


df_long <- df_long %>%
  dplyr::group_by(Sample) %>%
  dplyr::mutate(Percent = Count/sum(Count))


is_alluvia_form(df_long, axes = 1:2, silent = FALSE)

# 绘图
ggplot(df_long,
       aes(x = Sample, 
           stratum = CellType, 
           alluvium = CellType,
           y = Percent,
           fill = CellType, 
           label = CellType)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", alpha = 0.5) +
  geom_stratum(width = 0.3, color = "grey") +
  scale_y_continuous(labels = percent_format()) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank()) +
  labs(y = "Cell Proportion (%)", x = NULL) +
  scale_fill_manual(values = c("#2077B4","#AEC7E8","#FF7F0D"))









