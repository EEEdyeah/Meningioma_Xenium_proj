

####################################################################################################################################################################################
suppressMessages({
    library(GeneNMF)
    library(Seurat)
    library(ggplot2)
    library(UCell)
    library(patchwork)
    library(Matrix)
    library(RcppML)
    library(viridis)
    library(dplyr)
    library(stringr)
    library(argparser)
    library(msigdbr)
    library(fgsea)
    library(RColorBrewer)
    })

print(sessionInfo())

mk.outdir <- function(dir) {
    if (!dir.exists(dir)){
        dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    }
    return(dir)
}

save_gg <- function(p, filename, width, height, format = c('pdf', 'png')) {
    for (d in format) {
        ggplot2::ggsave(filename = str_glue("{filename}.{d}"), plot = p, width = width, height = height, device = d)
    }
}


argv <- arg_parser('')
argv <- add_argument(argv,"--rds", help = "rds file")
argv <- add_argument(argv,"--subcelltype", help = "sub celltype",default='all')
argv <- add_argument(argv,"--group_col", type = "character",help = "", default = 'None')
argv <- add_argument(argv,"--celltype_col", type = "character", help = "")
argv <- add_argument(argv,"--outdir",help = "outdir")
argv <- add_argument(argv,"--maxGene",help = "maxGene for each Programs", default = 200)
argv <- add_argument(argv,"--knum",help = "k number for Programs", default = 10)
argv <- parse_args(argv)

rds_path <- argv$rds
maxGene <- argv$maxGene
knum <- argv$knum
outdir <- argv$outdir
mk.outdir(outdir)
subcelltype<-argv$subcelltype
celltype_col<-argv$celltype_col
group_col<-argv$group_col

# NMF
rds <- readRDS(rds_path)
if(subcelltype != 'all'){
  if (grepl(',',subcelltype)){subcelltype = unlist(strsplit(subcelltype, split = ','))}
  rds <- subset(rds, celltype_col %in% subcelltype)
}
print (head(rds))

ndim <- 15
rds <- FindVariableFeatures(rds, nfeatures = 2000)
rds <- runNMF(rds, k = ndim, assay="RNA")
rds.list <- SplitObject(rds, split.by = "sample")
geneNMF.programs <- multiNMF(rds.list, assay="RNA", slot="data", k= 4:9) # 这里不要更改，否则影响结果

geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs,
                                        nMP=knum,metric ="jaccard",
                                        max.genes=maxGene)
#ph <- plotMetaPrograms(geneNMF.metaprograms, jaccard.cutoff = c(0,1),show_rownames = TRUE, show_colnames = TRUE)
#anno_colors <- brewer.pal(n=knum, name="Paired")
#names(anno_colors) <- names(geneNMF.metaprograms$metaprograms.genes)

#ph <- plotMetaPrograms(geneNMF.metaprograms,annotation_colors = anno_colors)
ph <- plotMetaPrograms(geneNMF.metaprograms)
save_gg(ph, str_glue("{outdir}/NMF"), 11.5, 10)
saveRDS(geneNMF.metaprograms, str_glue("{outdir}/geneNMF.rds"))

# UCell
markers <- geneNMF.metaprograms$metaprograms.genes
marker_score <- AddModuleScore_UCell(rds,features=markers)
a <- colnames(marker_score@meta.data) %>% str_subset("_UCell")
pu <- FeaturePlot(marker_score,features = a,order = T, ncol = 4, cols = viridis(256),reduction = 'umap')
save_gg(pu, str_glue("{outdir}/umap.MP"), 20, 5*ceiling(knum/4))

pv = VlnPlot(marker_score,features = a,stack = TRUE,flip=TRUE)
save_gg(pv, str_glue("{outdir}/VlnPlot.MP"), 3+length(unique(rds$seurat_clusters)), 8)

p1<-VlnPlot(marker_score,features = a,stack = TRUE,flip=TRUE,group.by ='sample')
save_gg(p1, str_glue("{outdir}/VlnPlot.MP_sample"), 3+length(unique(rds$sample)), 8)
p2<-VlnPlot(marker_score,features = a,stack = TRUE,flip=TRUE,group.by =group_col)
save_gg(p2, str_glue("{outdir}/VlnPlot.MP_group"), 3+length(unique(rds@meta.data[[group_col]])), 8)


MPs <- names(markers)
Gene_Dir <- str_glue("{outdir}/MP_Gene_list")
mk.outdir(Gene_Dir)
for (i in MPs){
    tmp = data.frame("gene"= markers[[i]])
    tmp = na.omit(tmp)
    write.table(tmp,str_glue("{Gene_Dir}/{i}.genes.txt"),row.names=F,quote=F)
}


top_p <- lapply(geneNMF.metaprograms$metaprograms.genes, function(program) {
  runGSEA(program, universe=rownames(rds), category = "H")
})

for (i in 1:length(top_p)){
        for(j in 1:ncol(top_p[[i]])){
                if(is.list(top_p[[i]][[j]])){

                top_p[[i]][[j]] <- sapply(top_p[[i]][[j]], function(x) paste(x, collapse = ","))
                }
        }

        utils::write.table(as.data.frame(top_p[[i]]), file = str_glue("{Gene_Dir}/gsea_H_MP{i}.xls"),sep = "\t",quote=F,row.names=F)

}



####
library(plotly)
# 绘制散点图
fig <- plot_ly(
  data = scMarkers,
  x = ~scMarkers$pct.1, 
  y = ~scMarkers$pct_diff, 
  z = ~scMarkers$avg_log2FC, 
  marker = list(
    color = ~pct.1,                
    colorscale = c('red', 'yellow', 'green'), 
    size = ~scMarkers$avg_log2FC * 2,                 
    opacity = 1, 
    showscale = FALSE
  ),
  type = 'scatter3d', 
  mode = 'markers'
)

fig <- fig %>% add_trace(
  x = scMarkers1$pct.1, 
  y = scMarkers1$pct_diff, 
  z = scMarkers1$avg_log2FC,
  text = rownames(scMarkers1),
  type = 'scatter3d', 
  mode = 'text+markers',
  marker = list(color = 'black', size = 0.5, opacity = 0.5, symbol = 'circle'),
  textfont = list(size = 16)
)

for (i in 1:nrow(scMarkers1)) {
  fig <- fig %>% add_trace(
    x = c(0, scMarkers1$pct.1[i]), 
    y = c(0, scMarkers1$pct_diff[i]), 
    z = c(0, scMarkers1$avg_log2FC[i]),
    type = 'scatter3d',
    mode = 'lines',
    marker = list(color = 'black', size = 0.5, opacity = 0, symbol = 'circle'),
    line = list(color = 'gray', width = 2, dash = 'dash')
  )
}

fig <- fig %>% layout(
  scene = list(
    xaxis = list(showline = TRUE, zeroline = FALSE,
      title = list(
        text = 'Expression Proportion', 
        font = list(size = 15)
      )
    ),
    yaxis = list(showline = TRUE, zeroline = FALSE,
      title = list(
        text = 'Differential Expression Proportion', 
        font = list(size = 15) 
      )
    ),
    zaxis = list(showline = TRUE, zeroline = FALSE,
      title = list(
        text = 'log2FC', 
        font = list(size = 15) 
      )
    )
  )
)
fig




