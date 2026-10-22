library(monocle)
library(RColorBrewer)
library(Seurat)
library(dplyr)
library(msigdbr)
genesets <- msigdbr(species = "Homo sapiens", category = "H") 
genesets <- subset(genesets, select = c("gs_name","gene_symbol")) %>% as.data.frame()
genesets <- split(genesets$gene_symbol, genesets$gs_name)


malignant<-readRDS('1malignant_CytoTRACE2.rds')
malignant<-subset(malignant, niche %in% c('CN1','CN3','CN6','CN7'))

malignant <- UCell::AddModuleScore_UCell(malignant, features = list(c('SLC2A1','COL18A1','MCTP2','CCN4','NREP')), ncores=8, name = "dedifferentation")
malignant <- UCell::AddModuleScore_UCell(malignant, features = genesets, ncores=8, name = "")
cell_cycle<-c("BIRC5", "LMNB1", "CCNB1", "TOP2A", "TPX2", "NUSAP1", "MKI67", "PTTG1", "UBE2T", "RRM2")
malignant <- UCell::AddModuleScore_UCell(malignant, features = list(cell_cycle), ncores=8, name = "cell_cycle")
malignant <- UCell::AddModuleScore_UCell(malignant, features = list(c('COL18A1')), ncores=8, name = "COL18A1")


Idents(malignant)<-'niche'
ordering_genes <- FindAllMarkers(malignant, logfc.threshold = 1, test.use = "wilcox", min.pct = 0,padj=0.05, only.pos = TRUE,assay = "RNA")



set.seed(970619)
malignant <- SketchData(object = malignant,ncells = 50000,method = "LeverageScore",sketched.assay = "sketch")
DefaultAssay(malignant) <- "sketch"
Mono_matrix <- GetAssayData(malignant,assay = "sketch",layer = "counts")
feature_ann<-data.frame(gene_id=rownames(Mono_matrix),gene_short_name=rownames(Mono_matrix))
rownames(feature_ann)<-rownames(Mono_matrix)

Mono_fd<-new("AnnotatedDataFrame", data = feature_ann)
sample_ann <- malignant@meta.data[colnames(GetAssayData(malignant, assay = "sketch")), ]
rownames(sample_ann)<-colnames(Mono_matrix)
Mono_pd<-new("AnnotatedDataFrame", data =sample_ann)

Mono.cds<-newCellDataSet(Mono_matrix,phenoData =Mono_pd,featureData =Mono_fd,expressionFamily=negbinomial.size())


head(pData(Mono.cds))
head(fData(Mono.cds))

Mono.cds <- estimateSizeFactors(Mono.cds)
Mono.cds <- detectGenes(Mono.cds, min_expr = 0.1)

Mono.cds <- estimateDispersions(Mono.cds)


Mono.cds <- setOrderingFilter(Mono.cds, ordering_genes$gene)

Mono.cds <- reduceDimension(Mono.cds,method = 'DDRTree')

Mono.cds <- orderCells(Mono.cds)
head(pData(Mono.cds))

Mono.cds <- orderCells(Mono.cds,root_state = 2)#修改起始位置


#------------------------------------------------
pdf("./1monocle/1malignant_Monocle_CN1367.pdf",width=5,height=5)
plot_cell_trajectory(Mono.cds, color_by = "niche",cell_size = 0.75)+
                     scale_color_manual(breaks = c("CN6", "CN7", "CN3","CN1"), values=c("#00ACEB", "#FDD0A2","#E6553A", "#FF0049")) + labs(x='DC_1',y='DC_2')+
                     theme(axis.ticks = element_blank(),axis.text=element_blank(),axis.title = element_text(size=12,face="bold"),
                           axis.line.y=element_line(linetype=1,color="black",size=1.5),axis.line.x=element_line(linetype=1,color="black",size=1.5) )#+NoLegend()
plot_cell_trajectory(Mono.cds, color_by = "Pseudotime",cell_size = 0.75)+
                     #scale_color_manual(breaks = c("CN6", "CN7", "CN3","CN1"), values=c("#00ACEB", "#FDD0A2","#E6553A", "#FF0049")) + labs(x='DC_1',y='DC_2')+
                     theme(axis.ticks = element_blank(),axis.text=element_blank(),axis.title = element_text(size=12,face="bold"),
                           axis.line.y=element_line(linetype=1,color="black",size=1.5),axis.line.x=element_line(linetype=1,color="black",size=1.5) )#+NoLegend()
plot_cell_trajectory(Mono.cds, color_by = "State",cell_size = 0.75)+
                     #scale_color_manual(breaks = c("CN6", "CN7", "CN3","CN1"), values=c("#00ACEB", "#FDD0A2","#E6553A", "#FF0049")) + labs(x='DC_1',y='DC_2')+
                     theme(axis.ticks = element_blank(),axis.text=element_blank(),axis.title = element_text(size=12,face="bold"),
                           axis.line.y=element_line(linetype=1,color="black",size=1.5),axis.line.x=element_line(linetype=1,color="black",size=1.5) )#+NoLegend()
dev.off()


pdf("./1monocle/1malignant_signature_1dedifferentation.pdf", width=5, height=5)
plot_cell_trajectory(Mono.cds, color_by = "signature_1dedifferentation", cell_size = 0.3) +
  scale_color_gradientn(colours = c("lightgrey", "#DA0A16"),limits = c(0, 1),breaks = c(0, 0.25, 0.5, 0.75, 1)) +labs(title = 'signature_1dedifferentation',x='DC_1',y='DC_2') +theme_classic() +
  theme(
    panel.background = element_rect(fill = 'transparent', color = "white"),
    legend.position = c(0.95, 0.8),legend.direction = "vertical",legend.text = element_text(color = "black", size = 8, face = "bold"),legend.title = element_blank(),
    plot.title = element_text(face = 'bold', size = 12, hjust = 0.5),axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(size = 1),axis.line.x = element_line(linetype = 1, color = "black", size = 0.6),
    axis.line.y = element_line(linetype = 1, color = "black", size = 0.6)
  )
plot_cell_trajectory(Mono.cds, color_by = "cell_cycle", cell_size = 0.3) +
  scale_color_gradientn(colours = c("lightgrey", "lightgrey", "lightgrey", "#DA0A16"),values = scales::rescale(c(-0.5, 0, 0.1, 0.9)),limits = c(-0.5, 0.9) ,) +labs(title = 'cell_cycle',x='DC_1',y='DC_2') +theme_classic() +
  theme(
    panel.background = element_rect(fill = 'transparent', color = "white"),
    legend.position = c(0.95, 0.8),legend.direction = "vertical",
    legend.text = element_text(color = "black", size = 8, face = "bold"),legend.title = element_blank(),
    plot.title = element_text(face = 'bold', size = 12, hjust = 0.5),axis.ticks.x = element_blank(),axis.ticks.y = element_line(size = 1),
    axis.line.x = element_line(linetype = 1, color = "black", size = 0.6),axis.line.y = element_line(linetype = 1, color = "black", size = 0.6)
  )
dev.off()




genes <- c("CLDN1", "PEG10", "SLC2A1", "COL18A1", "MCTP2", "CCN4", "CALCRL", "NREP","H19")
for (gene in genes) {
  expr <- exprs(Mono.cds)[gene, ]
  pData(Mono.cds)[[gene]] <- log2(expr + 0.01)
}
pdf("1malignant_degeneration_smooth.pdf", width=5, height=5)
plot_cell_trajectory(Mono.cds, color_by = "COL18A1", cell_size = 0.3) + 
  scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(100),breaks = c(0, 5)) +labs(title = 'COL18A1') +theme_classic() +
  theme(
    panel.background = element_rect(fill = 'transparent', color = "white"),
    legend.position = c(0.95, 0.8),legend.direction = "vertical",legend.text = element_text(color = "black", size = 8, face = "bold"),legend.title = element_blank(),
    plot.title = element_text(face = 'bold', size = 12, hjust = 0.5),
    axis.ticks.x = element_blank(),axis.ticks.y = element_line(size = 1),axis.line.x = element_line(linetype = 1, color = "black", size = 0.6),axis.line.y = element_line(linetype = 1, color = "black", size = 0.6)
  )
dev.off()




BEAM_res=BEAM(Mono.cds,branch_point = 1,cores = 20, progenitor_method = "duplicate")

BEAM_res=BEAM_res[,c("gene_short_name","pval","qval")]
save.image('./1monocle/1malignant_monocle2.rdata')

tmp1=plot_genes_branched_heatmap(Mono.cds[row.names(subset(BEAM_res,qval<1e-4)),],
                                 branch_point = 1,
                                 num_clusters = 5, #
                                 cores = 16,
                                 branch_labels = c("trace1", "trace2"),
                                 branch_colors = c("#7990C8", "#F05662","#979797" ),
                                 use_gene_short_name = T,
                                 show_rownames = F,
                                 return_heatmap = T 
                                 )
pdf("./1monocle/1malignant_BEAMheatmap.pdf",width = 5,height = 6)
tmp1$ph_res
dev.off()



library(org.Hs.eg.db)
library(clusterProfiler)
anno<-tmp1$annotation_row
anno$gene=rownames(anno)
#GO enrichment
for (i in 1:5){
    c1<-subset(anno,anno$Cluster== i)
    c1 = bitr(c1$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    ego_ALL <- enrichGO(gene = c1$ENTREZID,
                        OrgDb= org.Hs.eg.db,
                        ont = "BP",
                        pAdjustMethod = "BH",
                        minGSSize = 1,
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05,
                        readable = TRUE)
names=paste0('./1monocle/01malignant_BEAMheatmap_GO',i,'.pdf')
ggsave(names,barplot(ego_ALL,showCategory=20,title="Enrichment GO Top30")+theme(axis.text.y = element_text(size = 4)))
}




##################
library(circlize)
library(ComplexHeatmap)

# pseudotime 
pseudotime_raw <- pData(Mono.cds)$Pseudotime
cell_order <- order(pseudotime_raw)   

scores <- data.frame(
  row.names = colnames(Mono.cds),
  Dedifferentation = scale(pData(Mono.cds)$signature_1dedifferentation),
  E2F_TARGETS = scale(pData(Mono.cds)$HALLMARK_E2F_TARGETS),
  G2M_CHECKPOINT = scale(pData(Mono.cds)$HALLMARK_G2M_CHECKPOINT),
  Cell_cycle = scale(pData(Mono.cds)$cell_cycle),
  signature_1COL18A1 = scale(pData(Mono.cds)$signature_1COL18A1), 
  HALLMARK_MYC_TARGETS_V1 = scale(pData(Mono.cds)$HALLMARK_MYC_TARGETS_V1),
  HALLMARK_MYC_TARGETS_V2 = scale(pData(Mono.cds)$HALLMARK_MYC_TARGETS_V2),
  HALLMARK_MTORC1_SIGNALING = scale(pData(Mono.cds)$HALLMARK_MTORC1_SIGNALING),
  HALLMARK_WNT_BETA_CATENIN_SIGNALING = scale(pData(Mono.cds)$HALLMARK_WNT_BETA_CATENIN_SIGNALING),
  HALLMARK_IL6_JAK_STAT3_SIGNALING = scale(pData(Mono.cds)$HALLMARK_IL6_JAK_STAT3_SIGNALIN),
  HALLMARK_HEDGEHOG_SIGNALING = scale(pData(Mono.cds)$HALLMARK_HEDGEHOG_SIGNALING),
  HALLMARK_PI3K_AKT_MTOR_SIGNALING = scale(pData(Mono.cds)$HALLMARK_PI3K_AKT_MTOR_SIGNALING)
  )



# 
mat <- t(scores[cell_order, ])
mat <- t(scale(t(mat))) 

pseudotime <- pseudotime_raw[cell_order]
niche <- pData(Mono.cds)$niche[cell_order]
disease <- pData(Mono.cds)$CytoTRACE2_Potency[cell_order]
celltype <- pData(Mono.cds)$group[cell_order]

# annotation
top_anno <- HeatmapAnnotation(
  Pseudotime = pseudotime,
  Niche = niche,
  Disease = disease,
  Celltype = celltype,
  col = list(
    Pseudotime = colorRamp2(c(min(pseudotime), max(pseudotime)),
                            c("lightyellow", "red")),
    Niche = c("CN7"="#466983FF","CN6"="#40E0D0FF",
                "CN3"="#FF8C00FF","CN1"="#CE3D32FF"),
    Disease = c("Differentiated"="blue","Oligopotent"="pink",
                "Unipotent"="orange","Multipotent"="red"),
    Celltype = c("WHO1"="orange","WHO2"="green","WHO3"="purple")
  )
)

pdf("./1monocle/1malignant_Monocle_heatmap_pseudotime.pdf", width=20, height=5)
Heatmap(mat,
        name = "Expression",
        top_annotation = top_anno,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = TRUE,
        show_column_names = FALSE,
        column_order = 1:ncol(mat))  
dev.off()

