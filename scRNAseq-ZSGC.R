# Full cleaned analysis script
#
# This file keeps the merged code record after removing Chinese comments and local machine paths.
# Prefer the manuscript-structured scripts in R/ for GitHub browsing.

source("R/00_packages.R")

sce@meta.data$location[is.na(sce@meta.data$location)] <- "PLF"
sce@meta.data$dataset[is.na(sce@meta.data$dataset)] <- "GSE228030"
gse$location[gse$location == "Normal"] <- "N"
tmp <- Read10X(file.path("data", "inhouse", "1492435", "N"))
tmp <- CreateSeuratObject(counts = tmp, project = "1492435_N", min.cells = 3, min.features = 200)
tmp[["percent.mt"]] <- PercentageFeatureSet(tmp, pattern = "^MT-")
tmp@meta.data$group  <- "PM"
tmp@meta.data$location  <- "N"
tmp@meta.data$dataset  <- "OWN"
tmp <- subset(tmp, subset = nFeature_RNA > 500 & nFeature_RNA < 4500 & percent.mt < 25)
s.genes=Seurat::cc.genes.updated.2019$s.genes
g2m.genes=Seurat::cc.genes.updated.2019$g2m.genes
sce <- CellCycleScoring(sce, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
DimPlot(sce,reduction = "umap")
sce <- ScaleData(sce, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(sce))
sce <- RunPCA(sce, features = VariableFeatures(sce))
#,npcs = 110)
library(harmony)
sce <- RunHarmony(sce, "orig.ident")
sce <- RunUMAP(sce,  dims = 1:20,
#min.dist = 0.8,
               reduction = "harmony")
#sce <- RunTSNE(sce,  dims = 1:20,
#               reduction = "harmony")
sce <- FindNeighbors(sce,reduction = "harmony",dims = 1:20)
sce <- FindClusters(sce, resolution = 0.2)
cluster_map <- c(
  'Epi_c1_GKN1', 'Epi_c2_TM4SF1',   'Epi_c3_PGA3',  'Epi_c4_APOA1',   'Epi_c5_GHRL',    'Epi_c6_GIF'
)
sce@meta.data$subcelltype <- cluster_map[as.character(sce@meta.data$seurat_clusters)]
Cellratio <- prop.table(table(sce$subcelltype, sce$location), margin = 2)
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
sub <- factor(Cellratio$Var2, levels = c('PM','Ascites','PT','N'))
library(ggplot2)
ggplot(Cellratio) +
  geom_bar(aes(x =sub, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+
  theme_classic() +
#scale_fill_npg()+
#scale_fill_manual(values = pal)+
  labs(x='Tissue',y = 'Ratio')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
library(Startrac)
library(ggplot2)
library(tictoc)
library(ggpubr)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(tidyverse)
library(sscVis)
in.dat <- sce@meta.data
R_oe <- calTissueDist(in.dat,
                      byPatient = F,
                      colname.cluster = "subcelltype",
                      colname.patient = "orig.ident",
                      colname.tissue = "location",
                      method = "chisq",
                      min.rowSum = 0)
R_oe
R_oe_adjusted <- pmin(R_oe, 2)
col_fun <- colorRamp2(c(0, 1,2), c("#FFE4B5", "#FFA500", "#FF4500"))
Heatmap(as.matrix(R_oe_adjusted),
        show_heatmap_legend = TRUE,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        row_names_side = 'right',
        show_column_names = TRUE,
        show_row_names = TRUE,
        col = col_fun,
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        heatmap_legend_param = list(
          title = "R_o/e",
          at =  seq(0, 2, by = 0.5),
          labels =  seq(0, 2, by = 0.5),
          legend_gp = gpar(fill = col_fun(seq(0, 2, by = 0.5)))
        ),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", R_oe_adjusted[i, j]), x, y, gp = gpar(fontsize = 8))
        }
)
pyr <- c("CAD","DHODH","UMPS","TYMS",
         "UCK1","UCK2","TK1","TK2","CDA","UCKL1","DCK",
         "SLC29A1","SLC29A2","SLC29A3","SLC29A4",
         "SLC28A1","SLC28A2","SLC28A3"
)
gene_cell_exp <- AverageExpression(sce,
                                   features = pyr,
                                   group.by = 'subcelltype',
                                   slot = 'data')
gene_cell_exp <- as.data.frame(gene_cell_exp$RNA)
library(ComplexHeatmap)
df <- data.frame(colnames(gene_cell_exp))
colnames(df) <- 'class'
library(ggsci)
mypal <- pal_npg()(6)
mypal
library(scales)
show_col(mypal)
top_anno = HeatmapAnnotation(df = df,
                             border = T,
                             show_annotation_name = F,
                             gp = gpar(col = 'black'),
                             col = list(class = c('CAPN8_Epi'= "#E64B35FF",
                                                  'TM4SF1_Epi'="#4DBBD5FF",
                                                  'TFF1_Epi'="#00A087FF",
                                                  'TFF3_Epi'="#3C5488FF",
                                                  'PGA3_Epi'="#F39B7FFF",
                                                  'UBE2C_Epi'="#8491B4FF"))
)
top_anno = HeatmapAnnotation(df = df,
                             border = T,
                             show_annotation_name = F,
                             gp = gpar(col = 'black'),
                             col = list(class = c('PT'= "#E64B35FF",
                                                  'Ascites'="#4DBBD5FF",
                                                  'M'="#00A087FF"))
)
marker_exp <- t(scale(t(gene_cell_exp),scale = T,center = T))
library(RColorBrewer)
coul <- colorRampPalette(brewer.pal(9, "OrRd"))(50)
Heatmap(marker_exp,
        cluster_rows = F,
        cluster_columns = F,
        show_column_names = T,
        show_row_names = T,
        column_title = NULL,
        heatmap_legend_param = list(
          title=' '),
        col = coul,
        border = 'black',
        rect_gp = gpar(col = "black", lwd = 1),
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10))
#top_annotation = top_anno)
cluster_map <- c(
  "23" = "Neutrophils",
  "7" = "Mast",
  "0" = "Macro_c1_CCL4",
  "4" = "Macro_c2_C1QC","6" = "Macro_c2_C1QC",
  "3" = "Macro_c3_CRIP1","12" = "Macro_c3_CRIP1",
  "11" = "Macro_c4_SPP1","18" = "Macro_c4_SPP1",
  "13" = "Macro_c5_CD1C",
  "9" = "Macro_c6_LYVE1/TIMD4","20" = "Macro_c6_LYVE1/TIMD4",
  "10" = "Macro_c7_LYVE1","15" = "Macro_c7_LYVE1", "22" = "Macro_c7_LYVE1",
  "21" = "Macro_c8_CXCL10", "19" = "Macro_c8_CXCL10",
  "1" = "Mono_c1_CCL20",
  "2" = "Mono_c2_S100A8",
  "14" = "Mono_c3_LILRA5",
  "5" = "DC_c1_CD1C", "8" = "DC_c1_CD1C",
  "16" = "DC_c2_CLEC9A",
  "17" = "DC_c3_LAMP3"
)
sce@meta.data$subcelltype <- cluster_map[as.character(sce@meta.data$seurat_clusters)]
DimPlot(sce,group.by = "subcelltype",label = T)
cluster_map <- c(
  "1" = "Meso_c1_ITLN1",
  "4" = "Meso_c2_SAA1",
  "0" = "Fib_c1_CFD",
  "2" = "Fib_c2_CXCL14",
  "5" = "Fib_c3_COL1A1",
  "6" = "Fib_c4_FTL",
  "8" = "Fib_c5_CCL11",
  "9" = "Fib_c6_CD74",
  "10" = "Fib_c7_CEMIP","11" = "Fib_c7_CEMIP",
  "3" = "Pericyte",
  "7" = "SMC"
)
sce@meta.data$subcelltype <- cluster_map[as.character(sce@meta.data$seurat_clusters)]
DimPlot(sce,group.by = "subcelltype",label = T,split.by = "dataset")
DotPlot(sce,features = c("CCL21","PROX1", #LYM
                         "HEY1","IGFBP3", #ART
                         "CD36","CA4", #CAP
                         "ACKR1", #VEIN
                         "KDR",'SPARC' #TIP
                         ))+coord_flip()
cluster_map <- c(
  "0" = "Hypoxia",
  "1" = "Tip_c1_ESM1",
  "5" = "Tip_c2_KDR",
  "4" = "Art_c1_GJA4",
  "7" = "Cap_c1_FABP4",
  "2" = "Vein_c1_IL6",
  "3" = "Vein_c2_VCAN",
  "6" = "Lym_c1_CCL21"
)
sce@meta.data$subcelltype <- cluster_map[as.character(sce@meta.data$seurat_clusters)]
DimPlot(sce,group.by = "subcelltype",label = T)
cluster_map <- c(
  "0" = "CD8_c1_IFNG",
  "1" = "CD8_c2_SLC7A5",
  "2" = "CD8_c3_GZMK",
  "12" = "CD8_c4_TLE5",
  "16" = "CD8_c5_CXCL13",
  "17" = "CD8_c6_ISG15",
  "15" = "γδT",
  "9" = "NKT",
  "3" = "CD4_c1_FOS",
  "7" = "CD4_c2_SELL","10" = "CD4_c2_SELL",
  "5" = "CD4_c3_FOXP3","19" = "CD4_c3_FOXP3",
  "14" = "CD4_c4_CXCL13",
  "8" = "CD4_c5_IL7R",
  "4" = "NK_c1_FCGR3A","11" = "NK_c1_FCGR3A",
  "6" = "NK_c2_FCER1G","13" = "NK_c2_FCER1G",
  "18" = "NK_c3_AREG"
)
sce@meta.data$subcelltype <- cluster_map[as.character(sce@meta.data$seurat_clusters)]
DimPlot(sce,group.by = "subcelltype")
celltype_color = c(
  "Mast_cells" = "#b3d1e3","Epithelial_cells" = "#7ec6da",
  "Plasma_cells"="#398caa",
  "Myeloid_cells"="#63b246",
  "Endothelial_cells"="#e13d28","Fibroblasts"="#91656d",
  "Smooth_muscle_cells"="#c0a36f",
  "T_cells_NK_cells"="#e7692f",
  "B_cells"="#cfb3cb"
)
DimPlot(sce, group.by='celltype', cols=celltype_color, pt.size=1, raster=T, label = T)+NoLegend()
colors <- c("Ascites" = "#1F77B4FF",
            "PT" = "#FF7F0EFF",
            "N" = "#2CA02CFF",
            "PM" = "#D62728FF",
            "PLF" = "#9467BDFF",
            "PBMC" = "#8C564BFF",
            "NP" = "#E377C2FF"
)
library(Startrac)
library(ggplot2)
library(tictoc)
library(ggpubr)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(tidyverse)
library(sscVis)
in.dat <- sce@meta.data
R_oe <- calTissueDist(in.dat,
                      byPatient = F,
                      colname.cluster = "subcelltype",
                      colname.patient = "orig.ident",
                      colname.tissue = "location",
                      method = "chisq",
                      min.rowSum = 0)
R_oe
R_oe_adjusted <- pmin(R_oe, 2)
col_fun <- colorRamp2(c(0, 1,2), c("#f2f4c9","#349eb7", "#1a2852"))
col_fun <- colorRamp2(c(0, 1, 2), c("#FFE4B5", "#FFA500", "#FF4500"))
R_oe_adjusted <- t(R_oe_adjusted)
Heatmap(as.matrix(R_oe_adjusted),
        show_heatmap_legend = TRUE,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        row_names_side = 'right',
        show_column_names = TRUE,
        show_row_names = TRUE,
        col = col_fun,
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        heatmap_legend_param = list(
          title = "R_o/e",
          at =  seq(0, 2, by = 0.5),
          labels =  seq(0, 2, by = 0.5),
          legend_gp = gpar(fill = col_fun(seq(0, 2, by = 0.5)))
        ),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", R_oe_adjusted[i, j]), x, y, gp = gpar(fontsize = 8))
        }
)
library(ggsci)
cl=pal_nejm("default",alpha = 0.5)(10)
Idents(sce) = "celltype"
markers <- FindAllMarkers(sce, logfc.threshold = 1, min.pct = 0.25, only.pos = T)
top5 <- markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
library(ggplot2)
p <- DotPlot(sce, features = top5$gene, #top5$gene,
             cols = cl, group.by = "celltype") #+coord_flip()
exp <- p$data
library(forcats)
exp$features.plot <- as.factor(exp$features.plot)
exp$features.plot <- fct_inorder(exp$features.plot)
exp$id <- as.factor(exp$id)
exp$id <- fct_inorder(exp$id)
ggplot(exp,aes(x=features.plot,y= id))+
  geom_point(aes(size=`pct.exp`,
                 color=`avg.exp.scaled`))+
  geom_point(aes(size=`pct.exp`,color=`avg.exp.scaled`),
             shape=21,color="black",stroke=1)+
  theme(panel.background =element_blank(),
        axis.line=element_line(colour="black", size = 1),
        panel.grid = element_blank(),
        axis.text.x=element_text(size=11,color="black",angle=90),
        axis.text.y=element_text(size=11,color="black"))+
  scale_color_gradientn(colors = colorRampPalette(c("white", "#00C1D4", "#FFED99","#FF7600"))(10))+
  labs(x=NULL,y=NULL)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyverse)
nature_col = c("#AEC7E8","#FFBB78","#9467BD","#7200DA",
               "#17BECF","#FF7F0E","#C49C94","#2CA02C","#8C564B",
               "#E377C2","#D62728","#FF9896","#98DF8A","#BCBD22","#C5B0D5",
               "#00A087FF","#3C5488FF")
getPalette = colorRampPalette(brewer.pal(9,"Set1"))(22)
nature_col = col_CD8
#scales::show_col(pal_npg(palette = c("nrc"), alpha = 1)(8))
df <- data.frame(sce@meta.data[,c('subcelltype')],
                 sce@reductions$umap@cell.embeddings[,1:2])
#df <- read.csv("umap.csv")
colnames(df)
table(df$sce.meta.data...c..subcelltype...)
cluster_order <- c('Epi_c1_GKN1', 'Epi_c2_TM4SF1',
                   'Epi_c3_PGA3',  'Epi_c4_APOA1',   'Epi_c5_GHRL',
                   'Epi_c6_GIF')
col_epi <- pal_npg()(6)
col_cluster <- setNames(col_epi,
                        c('Epi_c1_GKN1', 'Epi_c2_TM4SF1',
                          'Epi_c3_PGA3',  'Epi_c4_APOA1',   'Epi_c5_GHRL',
                          'Epi_c6_GIF'))
ggplot(df,aes(x= UMAP_1 , y = UMAP_2 ,col=factor(sce.meta.data...c..subcelltype..., levels = cluster_order))) +
  geom_point(size = 2, shape=16)+
  scale_color_manual("",values = col_cluster)+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank(),
        legend.key=element_rect(fill='white'),
        legend.text = element_text(size=15),
        legend.key.size=unit(0.6,'cm'),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  guides(color = guide_legend(override.aes = list(size=6)))+
  geom_segment(aes(x = min(UMAP_1) , y = min(UMAP_2),xend = min(UMAP_1)+3, yend = min(UMAP_2)),
               colour = "black", size=0.5,
               arrow = arrow(length = unit(0.2,"cm"),
                             type = "closed"))+
  geom_segment(aes(x = min(UMAP_1), y = min(UMAP_2),xend = min(UMAP_1),yend = min(UMAP_2)+1.5),
               colour = "black", size=0.5,arrow = arrow(length = unit(0.2,"cm"),
                                                        type = "closed")) +
  annotate("text", x = min(df$UMAP_1) +1.4, y = min(df$UMAP_2) -0.3, label = "UMAP_1",
           color="black",size = 5) +
  annotate("text", x = min(df$UMAP_1) -0.8, y = min(df$UMAP_2) + 0.8, label = "UMAP_2",
           color="black",size = 5,angle=90)
cell <- df %>%group_by(sce.meta.data...c..subcelltype...) %>%
  summarise(UMAP_1 = median(UMAP_1),
            UMAP_2 = median(UMAP_2))
rownames(cell) <- cell$sce.meta.data...c..subcelltype...
A <- cell[cluster_order,]
a <- c(1:6)
A$ID <- a
A$ID <- as.factor(A$ID)
p <- ggplot(df,aes(x= UMAP_1 , y = UMAP_2 ,col=factor(sce.meta.data...c..subcelltype..., levels = cluster_order))) +
  geom_point(size = 0.2, shape=16)+
  scale_color_manual("",values = col_cluster)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = 'none',
        axis.title = element_text(colour = 'black', size = 15, hjust = 0),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  geom_text_repel(data=A, aes(label=ID),color="black", size=6, point.padding = 0.3)
B <- A
B$x <- 1
B$lab <- c(6:1)
leg <- B %>%
  mutate(Response_round = round(5 * lab) / 5) %>%
  group_by(x, Response_round) %>%
  mutate(x = 0.1 * (seq_along(Response_round) - (0.5 * (n() + 1)))) %>%
  ungroup() %>%
  mutate(x = x + as.numeric(as.factor(x))) %>%
  ggplot(aes(x = x, y = lab)) +
  geom_point(shape = 21, size = 8, aes(x = x, y = Response_round, fill=sce.meta.data...c..subcelltype...)) +
  geom_text(aes(label = ID, x = x, y = Response_round), size = 6)+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        legend.position = 'none',
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  scale_fill_manual(values = col_cluster)+
  annotate("text", x = 1, y = 23.5, label = "Cluster",
           color="black",size = 6)+
  geom_text(aes(label = sce.meta.data...c..subcelltype..., x = x+0.001,
                y = Response_round), size = 5, hjust=0)+
  scale_x_continuous(expand=c(-0.01,0.01))
library(cowplot)
plotlist <- list(p, leg)
plot_grid(plotlist = plotlist, ncol = 2, align="hv", rel_widths = 1, hjust = 0.1)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyverse)
df <- data.frame(sce@meta.data[,c('subcelltype')],
                 sce@reductions$umap@cell.embeddings[,1:2])
#df <- read.csv("umap.csv")
colnames(df)
table(df$sce.meta.data...c..subcelltype...)
cluster_order <- c('Fib_c1_CFD', 'Fib_c2_CXCL14',
                   'Fib_c3_COL1A1', 'Fib_c4_FTL',  'Fib_c5_CCL11',
                   'Fib_c6_CD74',  'Fib_c7_CEMIP',
                   'Meso_c1_ITLN1',  'Meso_c2_SAA1',
                    'Pericyte',     'SMC',
                   'Art_c1_GJA4', 'Cap_c1_FABP4',
                    'Lym_c1_CCL21', 'Tip_c1_ESM1',
                   'Tip_c2_KDR' , 'Vein_c1_IL6',
                   'Vein_c2_VCAN', 'Hypoxia' )
col_fib <- c(pal_aaas()(9),"#c0a","#c0a36f")
col_endo <- pal_bmj()(8)
cols = c(col_fib,col_endo)
col_cluster <- setNames(cols,
                        c('Fib_c1_CFD', 'Fib_c2_CXCL14',
                          'Fib_c3_COL1A1', 'Fib_c4_FTL',  'Fib_c5_CCL11',
                          'Fib_c6_CD74',  'Fib_c7_CEMIP',
                          'Meso_c1_ITLN1',  'Meso_c2_SAA1',
                          'Pericyte',     'SMC',
                          'Art_c1_GJA4', 'Cap_c1_FABP4',
                          'Lym_c1_CCL21', 'Tip_c1_ESM1',
                          'Tip_c2_KDR' , 'Vein_c1_IL6',
                          'Vein_c2_VCAN', 'Hypoxia' ))
ggplot(df,aes(x= UMAP_1 , y = UMAP_2 ,col=factor(sce.meta.data...c..subcelltype..., levels = cluster_order))) +
  geom_point(size = 2, shape=16)+
  scale_color_manual("",values = col_cluster)+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank(),
        legend.key=element_rect(fill='white'),
        legend.text = element_text(size=15),
        legend.key.size=unit(0.6,'cm'),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  guides(color = guide_legend(override.aes = list(size=6)))+
  geom_segment(aes(x = min(UMAP_1) , y = min(UMAP_2),xend = min(UMAP_1)+3, yend = min(UMAP_2)),
               colour = "black", size=0.5,
               arrow = arrow(length = unit(0.2,"cm"),
                             type = "closed"))+
  geom_segment(aes(x = min(UMAP_1), y = min(UMAP_2),xend = min(UMAP_1),yend = min(UMAP_2)+1.5),
               colour = "black", size=0.5,arrow = arrow(length = unit(0.2,"cm"),
                                                        type = "closed")) +
  annotate("text", x = min(df$UMAP_1) +1.4, y = min(df$UMAP_2) -0.3, label = "UMAP_1",
           color="black",size = 5) +
  annotate("text", x = min(df$UMAP_1) -0.8, y = min(df$UMAP_2) + 0.8, label = "UMAP_2",
           color="black",size = 5,angle=90)
cell <- df %>%group_by(sce.meta.data...c..subcelltype...) %>%
  summarise(UMAP_1 = median(UMAP_1),
            UMAP_2 = median(UMAP_2))
rownames(cell) <- cell$sce.meta.data...c..subcelltype...
A <- cell[cluster_order,]
a <- c(1:19)
A$ID <- a
A$ID <- as.factor(A$ID)
p <- ggplot(df,aes(x= UMAP_1 , y = UMAP_2 ,col=factor(sce.meta.data...c..subcelltype..., levels = cluster_order))) +
  geom_point(size = 0.2, shape=16)+
  scale_color_manual("",values = col_cluster)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = 'none',
        axis.title = element_text(colour = 'black', size = 15, hjust = 0),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  geom_text_repel(data=A, aes(label=ID),color="black", size=6, point.padding = 0.3)
B <- A
B$x <- 1
B$lab <- c(19:1)
leg <- B %>%
  mutate(Response_round = round(5 * lab) / 5) %>%
  group_by(x, Response_round) %>%
  mutate(x = 0.1 * (seq_along(Response_round) - (0.5 * (n() + 1)))) %>%
  ungroup() %>%
  mutate(x = x + as.numeric(as.factor(x))) %>%
  ggplot(aes(x = x, y = lab)) +
  geom_point(shape = 21, size = 8, aes(x = x, y = Response_round, fill=sce.meta.data...c..subcelltype...)) +
  geom_text(aes(label = ID, x = x, y = Response_round), size = 6)+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        legend.position = 'none',
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  scale_fill_manual(values = col_cluster)+
  annotate("text", x = 1, y = 23.5, label = "Cluster",
           color="black",size = 6)+
  geom_text(aes(label = sce.meta.data...c..subcelltype..., x = x+0.001,
                y = Response_round), size = 5, hjust=0)+
  scale_x_continuous(expand=c(-0.01,0.01))
library(cowplot)
plotlist <- list(p, leg)
plot_grid(plotlist = plotlist, ncol = 2, align="hv", rel_widths = 1, hjust = 0.1)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyverse)
df <- data.frame(sce@meta.data[,c('subcelltype')],
                 sce@reductions$umap@cell.embeddings[,1:2])
#df <- read.csv("umap.csv")
colnames(df)
table(df$sce.meta.data...c..subcelltype...)
cluster_order <- c(    "Mono_c1_CCL20",
                       "Mono_c2_S100A8",
                       "Mono_c3_LILRA5",
                       "Macro_c1_CCL4",
                       "Macro_c2_C1QC",
                       "Macro_c3_CRIP1",
                       "Macro_c4_SPP1",
                       "Macro_c5_CD1C",
                       "Macro_c6_LYVE1/TIMD4",
                       "Macro_c7_LYVE1",
                       "Macro_c8_CXCL10",
                       "DC_c1_CD1C",
                       "DC_c2_CLEC9A",
                       "DC_c3_LAMP3",
                       "Neutrophils",
                       "Mast"
                )
col_mye <- pal_d3("category20c")(16)
col_cluster <- setNames(col_mye,
                        c("Mono_c1_CCL20",
                          "Mono_c2_S100A8",
                          "Mono_c3_LILRA5",
                          "Macro_c1_CCL4",
                          "Macro_c2_C1QC",
                          "Macro_c3_CRIP1",
                          "Macro_c4_SPP1",
                          "Macro_c5_CD1C",
                          "Macro_c6_LYVE1/TIMD4",
                          "Macro_c7_LYVE1",
                          "Macro_c8_CXCL10",
                          "DC_c1_CD1C",
                          "DC_c2_CLEC9A",
                          "DC_c3_LAMP3",
                          "Neutrophils",
                          "Mast"))
ggplot(df,aes(x= UMAP_1 , y = UMAP_2 ,col=factor(sce.meta.data...c..subcelltype..., levels = cluster_order))) +
  geom_point(size = 2, shape=16)+
  scale_color_manual("",values = col_cluster)+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank(),
        legend.key=element_rect(fill='white'),
        legend.text = element_text(size=15),
        legend.key.size=unit(0.6,'cm'),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  guides(color = guide_legend(override.aes = list(size=6)))+
  geom_segment(aes(x = min(UMAP_1) , y = min(UMAP_2),xend = min(UMAP_1)+3, yend = min(UMAP_2)),
               colour = "black", size=0.5,
               arrow = arrow(length = unit(0.2,"cm"),
                             type = "closed"))+
  geom_segment(aes(x = min(UMAP_1), y = min(UMAP_2),xend = min(UMAP_1),yend = min(UMAP_2)+1.5),
               colour = "black", size=0.5,arrow = arrow(length = unit(0.2,"cm"),
                                                        type = "closed")) +
  annotate("text", x = min(df$UMAP_1) +1.4, y = min(df$UMAP_2) -0.3, label = "UMAP_1",
           color="black",size = 5) +
  annotate("text", x = min(df$UMAP_1) -0.8, y = min(df$UMAP_2) + 0.8, label = "UMAP_2",
           color="black",size = 5,angle=90)
cell <- df %>%group_by(sce.meta.data...c..subcelltype...) %>%
  summarise(UMAP_1 = median(UMAP_1),
            UMAP_2 = median(UMAP_2))
rownames(cell) <- cell$sce.meta.data...c..subcelltype...
A <- cell[cluster_order,]
a <- c(1:16)
A$ID <- a
A$ID <- as.factor(A$ID)
p <- ggplot(df,aes(x= UMAP_1 , y = UMAP_2 ,col=factor(sce.meta.data...c..subcelltype..., levels = cluster_order))) +
  geom_point(size = 0.2, shape=16)+
  scale_color_manual("",values = col_cluster)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = 'none',
        axis.title = element_text(colour = 'black', size = 15, hjust = 0),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  geom_text_repel(data=A, aes(label=ID),color="black", size=6, point.padding = 0.3)
B <- A
B$x <- 1
B$lab <- c(16:1)
leg <- B %>%
  mutate(Response_round = round(5 * lab) / 5) %>%
  group_by(x, Response_round) %>%
  mutate(x = 0.1 * (seq_along(Response_round) - (0.5 * (n() + 1)))) %>%
  ungroup() %>%
  mutate(x = x + as.numeric(as.factor(x))) %>%
  ggplot(aes(x = x, y = lab)) +
  geom_point(shape = 21, size = 8, aes(x = x, y = Response_round, fill=sce.meta.data...c..subcelltype...)) +
  geom_text(aes(label = ID, x = x, y = Response_round), size = 6)+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        legend.position = 'none',
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  scale_fill_manual(values = col_cluster)+
  annotate("text", x = 1, y = 23.5, label = "Cluster",
           color="black",size = 6)+
  geom_text(aes(label = sce.meta.data...c..subcelltype..., x = x+0.001,
                y = Response_round), size = 5, hjust=0)+
  scale_x_continuous(expand=c(-0.01,0.01))
library(cowplot)
plotlist <- list(p, leg)
plot_grid(plotlist = plotlist, ncol = 2, align="hv", rel_widths = 1, hjust = 0.1)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyverse)
df <- data.frame(sce@meta.data[,c('subcelltype')],
                 sce@reductions$umap@cell.embeddings[,1:2])
#df <- read.csv("umap.csv")
colnames(df)
table(df$sce.meta.data...c..subcelltype...)
cluster_order <- c( 'CD4_c1_FOS',   'CD4_c2_SELL',  'CD4_c3_FOXP3', 'CD4_c4_CXCL13',   'CD4_c5_IL7R',
                    'CD8_c1_IFNG', 'CD8_c2_SLC7A5',   'CD8_c3_GZMK', 'CD8_c4_TLE5', 'CD8_c5_CXCL13',
                    'CD8_c6_ISG15',
                    'NK_c1_FCGR3A',  'NK_c2_FCER1G',    'NK_c3_AREG',
                    'NKT',           'γδT'
)
col_T <- pal_d3("category20b")(16)
col_cluster <- setNames(col_T,
                        c('CD4_c1_FOS',   'CD4_c2_SELL',  'CD4_c3_FOXP3', 'CD4_c4_CXCL13',   'CD4_c5_IL7R',
                          'CD8_c1_IFNG', 'CD8_c2_SLC7A5',   'CD8_c3_GZMK', 'CD8_c4_TLE5', 'CD8_c5_CXCL13',
                          'CD8_c6_ISG15',
                          'NK_c1_FCGR3A',  'NK_c2_FCER1G',    'NK_c3_AREG',
                          'NKT',           'γδT'))
ggplot(df,aes(x= UMAP_1 , y = UMAP_2 ,col=factor(sce.meta.data...c..subcelltype..., levels = cluster_order))) +
  geom_point(size = 2, shape=16)+
  scale_color_manual("",values = col_cluster)+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank(),
        legend.key=element_rect(fill='white'),
        legend.text = element_text(size=15),
        legend.key.size=unit(0.6,'cm'),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  guides(color = guide_legend(override.aes = list(size=6)))+
  geom_segment(aes(x = min(UMAP_1) , y = min(UMAP_2),xend = min(UMAP_1)+3, yend = min(UMAP_2)),
               colour = "black", size=0.5,
               arrow = arrow(length = unit(0.2,"cm"),
                             type = "closed"))+
  geom_segment(aes(x = min(UMAP_1), y = min(UMAP_2),xend = min(UMAP_1),yend = min(UMAP_2)+1.5),
               colour = "black", size=0.5,arrow = arrow(length = unit(0.2,"cm"),
                                                        type = "closed")) +
  annotate("text", x = min(df$UMAP_1) +1.4, y = min(df$UMAP_2) -0.3, label = "UMAP_1",
           color="black",size = 5) +
  annotate("text", x = min(df$UMAP_1) -0.8, y = min(df$UMAP_2) + 0.8, label = "UMAP_2",
           color="black",size = 5,angle=90)
cell <- df %>%group_by(sce.meta.data...c..subcelltype...) %>%
  summarise(UMAP_1 = median(UMAP_1),
            UMAP_2 = median(UMAP_2))
rownames(cell) <- cell$sce.meta.data...c..subcelltype...
A <- cell[cluster_order,]
a <- c(1:16)
A$ID <- a
A$ID <- as.factor(A$ID)
p <- ggplot(df,aes(x= UMAP_1 , y = UMAP_2 ,col=factor(sce.meta.data...c..subcelltype..., levels = cluster_order))) +
  geom_point(size = 0.2, shape=16)+
  scale_color_manual("",values = col_cluster)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = 'none',
        axis.title = element_text(colour = 'black', size = 15, hjust = 0),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  geom_text_repel(data=A, aes(label=ID),color="black", size=6, point.padding = 0.3)
B <- A
B$x <- 1
B$lab <- c(16:1)
leg <- B %>%
  mutate(Response_round = round(5 * lab) / 5) %>%
  group_by(x, Response_round) %>%
  mutate(x = 0.1 * (seq_along(Response_round) - (0.5 * (n() + 1)))) %>%
  ungroup() %>%
  mutate(x = x + as.numeric(as.factor(x))) %>%
  ggplot(aes(x = x, y = lab)) +
  geom_point(shape = 21, size = 8, aes(x = x, y = Response_round, fill=sce.meta.data...c..subcelltype...)) +
  geom_text(aes(label = ID, x = x, y = Response_round), size = 6)+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        legend.position = 'none',
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  scale_fill_manual(values = col_cluster)+
  annotate("text", x = 1, y = 23.5, label = "Cluster",
           color="black",size = 6)+
  geom_text(aes(label = sce.meta.data...c..subcelltype..., x = x+0.001,
                y = Response_round), size = 5, hjust=0)+
  scale_x_continuous(expand=c(-0.01,0.01))
library(cowplot)
plotlist <- list(p, leg)
plot_grid(plotlist = plotlist, ncol = 2, align="hv", rel_widths = 1, hjust = 0.1)
cluster_map <- c(
  "0" = "B_c1_TCL1A", "15" = "B_c1_TCL1A",
  "1" = "B_c2_NR4A2",
  "3" = "B_c3_MS4A1",
  "2" = "Pla_c1_SPINK4",
  "4" = "Pla_c2_IGLV2","14" = "Pla_c2_IGLV2",
  "5" = "Pla_c3_IGKV1","13" = "Pla_c3_IGKV1",
  "6" = "Pla_c4_IGLC1",
  "7" = "Pla_c5_IER3",
  "8" = "Pla_c6_IGKC",
  "9" = "Pla_c7_IGHV4",
  "10" = "Pla_c8_IGKV2",
  "11" = "Pla_c9_IGHG",
  "12" = "Pla_c10_IGKV3"
)
sce@meta.data$subcelltype <- cluster_map[as.character(sce@meta.data$seurat_clusters)]
DimPlot(sce,group.by = "subcelltype",label = T)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyverse)
df <- data.frame(sce@meta.data[,c('subcelltype')],
                 sce@reductions$umap@cell.embeddings[,1:2])
#df <- read.csv("umap.csv")
colnames(df)
table(df$sce.meta.data...c..subcelltype...)
cluster_order <- c(   "B_c1_TCL1A",
                    "B_c2_NR4A2",
                     "B_c3_MS4A1",
                    "Pla_c1_SPINK4",
                     "Pla_c2_IGLV2",
                    "Pla_c3_IGKV1",
                    "Pla_c4_IGLC1",
                    "Pla_c5_IER3",
                    "Pla_c6_IGKC",
                    "Pla_c7_IGHV4",
                     "Pla_c8_IGKV2",
                     "Pla_c9_IGHG",
                     "Pla_c10_IGKV3" )
col_B <- pal_nejm()(4)
col_pla <- pal_jco()(9)
cols = c(col_B,col_pla)
col_cluster <- setNames(cols,
                        c( "B_c1_TCL1A",
                           "B_c2_NR4A2",
                           "B_c3_MS4A1",
                           "Pla_c1_SPINK4",
                           "Pla_c2_IGLV2",
                           "Pla_c3_IGKV1",
                           "Pla_c4_IGLC1",
                           "Pla_c5_IER3",
                           "Pla_c6_IGKC",
                           "Pla_c7_IGHV4",
                           "Pla_c8_IGKV2",
                           "Pla_c9_IGHG",
                           "Pla_c10_IGKV3" ))
ggplot(df,aes(x= UMAP_1 , y = UMAP_2 ,col=factor(sce.meta.data...c..subcelltype..., levels = cluster_order))) +
  geom_point(size = 2, shape=16)+
  scale_color_manual("",values = col_cluster)+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank(),
        legend.key=element_rect(fill='white'),
        legend.text = element_text(size=15),
        legend.key.size=unit(0.6,'cm'),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  guides(color = guide_legend(override.aes = list(size=6)))+
  geom_segment(aes(x = min(UMAP_1) , y = min(UMAP_2),xend = min(UMAP_1)+3, yend = min(UMAP_2)),
               colour = "black", size=0.5,
               arrow = arrow(length = unit(0.2,"cm"),
                             type = "closed"))+
  geom_segment(aes(x = min(UMAP_1), y = min(UMAP_2),xend = min(UMAP_1),yend = min(UMAP_2)+1.5),
               colour = "black", size=0.5,arrow = arrow(length = unit(0.2,"cm"),
                                                        type = "closed")) +
  annotate("text", x = min(df$UMAP_1) +1.4, y = min(df$UMAP_2) -0.3, label = "UMAP_1",
           color="black",size = 5) +
  annotate("text", x = min(df$UMAP_1) -0.8, y = min(df$UMAP_2) + 0.8, label = "UMAP_2",
           color="black",size = 5,angle=90)
cell <- df %>%group_by(sce.meta.data...c..subcelltype...) %>%
  summarise(UMAP_1 = median(UMAP_1),
            UMAP_2 = median(UMAP_2))
rownames(cell) <- cell$sce.meta.data...c..subcelltype...
A <- cell[cluster_order,]
a <- c(1:13)
A$ID <- a
A$ID <- as.factor(A$ID)
p <- ggplot(df,aes(x= UMAP_1 , y = UMAP_2 ,col=factor(sce.meta.data...c..subcelltype..., levels = cluster_order))) +
  geom_point(size = 0.2, shape=16)+
  scale_color_manual("",values = col_cluster)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = 'none',
        axis.title = element_text(colour = 'black', size = 15, hjust = 0),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  geom_text_repel(data=A, aes(label=ID),color="black", size=6, point.padding = 0.3)
B <- A
B$x <- 1
B$lab <- c(13:1)
leg <- B %>%
  mutate(Response_round = round(5 * lab) / 5) %>%
  group_by(x, Response_round) %>%
  mutate(x = 0.1 * (seq_along(Response_round) - (0.5 * (n() + 1)))) %>%
  ungroup() %>%
  mutate(x = x + as.numeric(as.factor(x))) %>%
  ggplot(aes(x = x, y = lab)) +
  geom_point(shape = 21, size = 8, aes(x = x, y = Response_round, fill=sce.meta.data...c..subcelltype...)) +
  geom_text(aes(label = ID, x = x, y = Response_round), size = 6)+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        legend.position = 'none',
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  scale_fill_manual(values = col_cluster)+
  annotate("text", x = 1, y = 23.5, label = "Cluster",
           color="black",size = 6)+
  geom_text(aes(label = sce.meta.data...c..subcelltype..., x = x+0.001,
                y = Response_round), size = 5, hjust=0)+
  scale_x_continuous(expand=c(-0.01,0.01))
library(cowplot)
plotlist <- list(p, leg)
plot_grid(plotlist = plotlist, ncol = 2, align="hv", rel_widths = 1, hjust = 0.1)
new_order <- c(
  'Fib_c1_CFD', 'Fib_c2_CXCL14',
  'Fib_c3_COL1A1', 'Fib_c4_FTL',  'Fib_c5_CCL11',
  'Fib_c6_CD74',  'Fib_c7_CEMIP',
  'Meso_c1_ITLN1',  'Meso_c2_SAA1',
  'Pericyte',     'SMC'
)
sce$subcelltype <- factor(sce$subcelltype, levels = new_order)
col_fib <- c(pal_aaas()(9),"#c0a","#c0a36f")
DimPlot(sce,group.by = "subcelltype",raster = T,cols = col_cluster,label = TRUE)
library(Startrac)
library(ggplot2)
library(tictoc)
library(ggpubr)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(tidyverse)
library(sscVis)
in.dat <- sce@meta.data
R_oe <- calTissueDist(in.dat,
                      byPatient = F,
                      colname.cluster = "subcelltype",
                      colname.patient = "orig.ident",
                      colname.tissue = "location",
                      method = "chisq",
                      min.rowSum = 0)
R_oe
R_oe_adjusted <- pmin(R_oe, 2)
col_fun <- colorRamp2(c(0, 1, 2), c("#FFE4B5", "#FFA500", "#FF4500"))
Heatmap(as.matrix(R_oe_adjusted),
        show_heatmap_legend = TRUE,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        row_names_side = 'right',
        show_column_names = TRUE,
        show_row_names = TRUE,
        col = col_fun,
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        heatmap_legend_param = list(
          title = "R_o/e",
          at =  seq(0, 2, by = 0.5),
          labels =  seq(0, 2, by = 0.5),
          legend_gp = gpar(fill = col_fun(seq(0, 2, by = 0.5)))
        ),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", R_oe_adjusted[i, j]), x, y, gp = gpar(fontsize = 8))
        }
)
markers <-c("COL1A1", "COL10A1", "COL4A1", "MMP3", "ACTA2","TAGLN", #ECM
            "IL4","IL6", "IL11", "IL13", "TGFB1","LIF", "CSF2", #Cytokine
            "CXCL1","CXCL2", "CXCL8", "CXCL10", "CXCL12","CXCL14",
            "CCL2", "CCL8", #Chemokines
            "CFD", "C1QA", "C1QB", "C1QC", #Complement
            "HLA-DQB1","HLA-DRB1","HLA-DQA1","HLA-DRA","HLA-DMB","HLA-DOA" #MHC-II
)
markers <- as.data.frame(markers)
markers$x <- 1
markers$y <- seq(1:nrow(markers))
markers$group <- c(rep('ECM',6),rep('Cytokines',7),rep('Chemokine',8),rep('Complement',4),
                   rep("MHC-II",6))
marker_ave_exp <- AverageExpression(sce,assays = "RNA",
                                    features = markers$markers,
                                    group.by = "subcelltype",
                                    slot="data")
markers$markers=factor(markers$markers,levels = markers$markers)
data=base::apply(marker_ave_exp$RNA,1,
                 function(x) (x-mean(x))/sd(x))%>%t()%>%as.data.frame()%>%rownames_to_column('Gene')%>%reshape2::melt()
data$Gene=factor(data$Gene,levels = rev(unique(data$Gene)))
data$variable = factor(data$variable, levels = c( 'Fib_c1_CFD', 'Fib_c2_CXCL14',
                                                  'Fib_c3_COL1A1', 'Fib_c4_FTL',  'Fib_c5_CCL11',
                                                  'Fib_c6_CD74',  'Fib_c7_CEMIP',
                                                  'Meso_c1_ITLN1',  'Meso_c2_SAA1',
                                                  'Pericyte',     'SMC'))
#color_palette <- c("#7F3F00", "#B35806", "#E08214", "#FDAE61", "#FEE0B6", "#D8DAEB", "#B2ABD2", "#8073AC", "#542788", "#2D004B")
library(RColorBrewer)
color_palette <- colorRampPalette(brewer.pal(11, "RdYlBu"))(100)
p1=ggplot(data,aes(x = variable,y = Gene,fill=value))+
  geom_tile()+
  scale_y_discrete(expand = c(0,0))+
  scale_fill_gradientn(colors=rev(colorRampPalette(color_palette)(500)),
                       limits=c(-1.5,3),name="Z Score",
                       oob = scales::squish
  )+
  geom_vline(xintercept=as.numeric(cumsum(table(unique(data$variable)))+0.5),linetype=2)+
  geom_hline(yintercept=as.numeric(cumsum(c(6.5,4,8,7,6))),linetype=2)+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        text = element_text(face="bold"),
        axis.title = element_blank(),
        axis.ticks=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size = 10,colour = "black"))
p1
palette1 <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99")
p2=ggplot(markers,aes(x,reorder(y,-y),fill=group))+
  geom_tile()+
  geom_text(aes(label=markers),size=3)+
  scale_fill_manual(values = palette1)+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        text=element_text(),
        axis.text = element_blank(),
        axis.title =element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'none')+
  scale_x_continuous(expand = c(0,0))
color_labels <- c('ECM','Cytokines','Chemokine','Complement',
                  "MHC-II")
label_data <- data.frame(y = color_labels)
label_data$y=factor(label_data$y,levels = c('ECM','Cytokines','Chemokine','Complement',
                                            "MHC-II"))
p2_color_annotations <- ggplot(label_data, aes(x = 0, y = y)) +
  geom_tile(aes(fill = y), color = "white") +
  scale_fill_manual(values = palette1) +
  geom_text(aes(label = y), hjust = 0, size = 4, color = "black") +
  theme_void()
p2_color_annotations
heatmap_with_color_annotations = p2_color_annotations +p2 +p1 +
  plot_layout(ncol = 2, widths  = c(1, 3))
heatmap_with_color_annotations & theme(plot.margin = margin(0,0,0,0))
library(patchwork)
heatmap1 =p2+p1+plot_layout(ncol = 2, widths  = c(1, 3))
heatmap1 & theme(plot.margin = margin(0,0,0,0))
meso = c("WT1","DES","MSLN","KRT18","KRT19","UPK3B","CALB2","SLPI","LRRN4")
# Set color palette
# cividis：option E;
library(scCustomize)
library(viridis)
pal <- viridis(n = 10,option="C")
pal <- viridis(n = 15, option = "D", direction = -1)
FeaturePlot(object = sce,features = meso, cols = pal, order = T)
#viridis_plasma_dark_high
#viridis_plasma_light_high
#viridis_magma_light_high
#viridis_inferno_dark_high
#viridis_inferno_light_high
#viridis_dark_high
i=1
plots=list()
for (i in 1:length(meso)){
  plots[[i]]=FeaturePlot_scCustom(seurat_object = sce,
                                  colors_use = viridis_dark_high,
                                  features = meso[i])+NoLegend()+NoAxes()+
    theme(panel.border = element_rect(fill = NA,color = "black",
                                      size=1.5,linetype = "solid"))
}
library(patchwork)
p<-wrap_plots(plots, ncol = 3);p
sce = fib[, sample(1:ncol(fib),round(ncol(fib)/5)) ]
library(CytoTRACE2) #loading
results <- cytotrace2(sce,
                      species = "human",
                      is_seurat = TRUE,
                      full_model = FALSE,
                      slot_type = "counts"
)
annotation <- data.frame(phenotype = sce@meta.data$subcelltype) %>%
  set_rownames(., colnames(sce))
# plotting
plots <- plotData(cytotrace2_result = results,
                  annotation = annotation,
                  is_seurat = TRUE)
emb_1 <- sce@reductions$umap@cell.embeddings[,1]*1.3 # First dimension of your "harmony" reduction embedding
emb_2 <- sce@reductions$umap@cell.embeddings[,2]*2+2 # Second dimension of your "harmony" reduction embedding
# replace CytoTRACE2_UMAP embeddings
plots[["CytoTRACE2_UMAP"]][[1]][["data"]]["UMAP_1"] <- emb_1
plots[["CytoTRACE2_UMAP"]][[1]][["data"]]["UMAP_2"]<- emb_2
plots[["CytoTRACE2_Potency_UMAP"]][[1]][["data"]]["UMAP_1"] <- emb_1
plots[["CytoTRACE2_Potency_UMAP"]][[1]][["data"]]["UMAP_2"]<- emb_2
plots[["CytoTRACE2_Relative_UMAP"]][[1]][["data"]]["UMAP_1"] <- emb_1
plots[["CytoTRACE2_Relative_UMAP"]][[1]][["data"]]["UMAP_2"]<- emb_2
plots[["Phenotype_UMAP"]][[1]][["data"]]["UMAP_1"] <- emb_1
plots[["Phenotype_UMAP"]][[1]][["data"]]["UMAP_2"]<- emb_2
plots$CytoTRACE2_UMAP
plots$CytoTRACE2_Potency_UMAP
plots$CytoTRACE2_Relative_UMAP
plots$Phenotype_UMAP
plots$CytoTRACE2_Boxplot_byPheno
plots$CytoTRACE2_Relative_UMAP <- plots$CytoTRACE2_Relative_UMAP +
  scale_color_viridis_c(
    option = "viridis",
    name = "Relative Score"
  )
cluster_map <- c(
  "0" = "Meso_c1_ITLN1",
  "1" = "Meso_c2_SAA1"
  )
sce@meta.data$subcelltype <- cluster_map[as.character(sce@meta.data$seurat_clusters)]
DimPlot(sce,group.by = "subcelltype",label = T,cols = col_cluster)
col_data = c('#f9766e',"#b79f00",'#01ba38','#00bfc4','#619dff','#f564e3')
cluster_order <- c( "CXD","GSE130888","GSE183904","GSE206785","GSE206785","OWN"
)
col_cluster <- setNames(col_data,
                        c("CXD","GSE130888","GSE183904","GSE206785","GSE206785","OWN"))
DimPlot(sce,group.by = "dataset",label = T,cols = col_cluster)
DimPlot(sce,raster = T,group.by = "subcelltype",cols = col_cluster,label = FALSE)+NoLegend()
library(Startrac)
library(ggplot2)
library(tictoc)
library(ggpubr)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(tidyverse)
library(sscVis)
in.dat <- sce@meta.data
R_oe <- calTissueDist(in.dat,
                      byPatient = F,
                      colname.cluster = "subcelltype",
                      colname.patient = "orig.ident",
                      colname.tissue = "location",
                      method = "chisq",
                      min.rowSum = 0)
R_oe
R_oe_adjusted <- pmin(R_oe, 2)
col_fun <- colorRamp2(c(0, 1, 2), c("#FFE4B5", "#FFA500", "#FF4500"))
Heatmap(as.matrix(R_oe_adjusted),
        show_heatmap_legend = TRUE,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        row_names_side = 'right',
        show_column_names = TRUE,
        show_row_names = TRUE,
        col = col_fun,
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        heatmap_legend_param = list(
          title = "R_o/e",
          at =  seq(0, 2, by = 0.5),
          labels =  seq(0, 2, by = 0.5),
          legend_gp = gpar(fill = col_fun(seq(0, 2, by = 0.5)))
        ),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", R_oe_adjusted[i, j]), x, y, gp = gpar(fontsize = 8))
        }
)
genes <- c("FLT1","KDR","NRP1","FLT4","NRP2","TIE1","TEK","PDGFB", #Tumor vascularization
           "JAG1","EFNB2","HEY1","EPHB4","HES1", #Morphogenesis
           "NT5E","STAB1","CCL2", #Inflammation
           "HLA-A","HLA-B","HLA-C", #MHC-I
           "HLA-DMA","HLA-DMB","HLA-DOA","HLA-DPA1",
           "HLA-DPB1","HLA-DQA1","HLA-DQA2","HLA-DQB1","HLA-DRA",
           "HLA-DRB1","HLA-DRB5", #MHC-II
           "ROBO1","ROBO4","PLXND1","ROBO3","UNC5B", #Guidance receptors
           "ICAM1","VCAM1","SELE","SELP", #Endothelial Adhesion Molecules
           "IDO1","HAVCR2","CD274","PDCD1LG2",
           "FASLG","LGALS1" #Immunosuppressive Molecules
)
anergy <- c("ICAM1","VCAM1","SELE","SELP", #Endothelial Adhesion Molecules
            "IDO1","HAVCR2","CD274","PDCD1LG2",
            "FASLG","LGALS1" #Immunosuppressive Molecules
)
#trans <- c("SLC29A1","SLC29A2","SLC29A3","SLC29A4","SLC28A1","SLC28A2","SLC28A3")
library(Seurat)
library(ggplot2)
library(tidyverse)
markers <- genes
#markers <- anergy
markers <- as.data.frame(markers)
markers$x <- 1
markers$y <- seq(1:nrow(markers))
markers$group <- c(rep('Tumor vascularization',8),rep('Morphogenesis',5),rep('Inflammation',3),rep('MHC-I',3),
                   rep("MHC-II",11),rep("Guidance receptors",5),rep('Endothelial Adhesion Molecules',4),
                   rep('Immunosuppressive Molecules',6))
marker_ave_exp <- AverageExpression(sce,assays = "RNA",
                                    features = markers$markers,
                                    group.by = "subcelltype",
                                    slot="data")
markers$markers=factor(markers$markers,levels = markers$markers)
data=base::apply(marker_ave_exp$RNA,1,
                 function(x) (x-mean(x))/sd(x))%>%t()%>%as.data.frame()%>%rownames_to_column('Gene')%>%reshape2::melt()
data$Gene=factor(data$Gene,levels = rev(unique(data$Gene)))
data$variable = factor(data$variable, levels = c('Art_c1_GJA4', 'Cap_c1_FABP4',
                                                 'Lym_c1_CCL21', 'Tip_c1_ESM1',
                                                 'Tip_c2_KDR' , 'Vein_c1_IL6',
                                                 'Vein_c2_VCAN', 'Hypoxia'))
#color_palette <- c("#7F3F00", "#B35806", "#E08214", "#FDAE61", "#FEE0B6", "#D8DAEB", "#B2ABD2", "#8073AC", "#542788", "#2D004B")
library(RColorBrewer)
color_palette <- colorRampPalette(brewer.pal(11, "RdYlBu"))(100)
p1=ggplot(data,aes(x = variable,y = Gene,fill=value))+
  geom_tile()+
  scale_y_discrete(expand = c(0,0))+
  scale_fill_gradientn(colors=rev(colorRampPalette(color_palette)(500)),
                       limits=c(-2,2),name="Z Score",
                       oob = scales::squish
  )+
  geom_vline(xintercept=as.numeric(cumsum(table(unique(data$variable)))+0.5),linetype=2)+
  geom_hline(yintercept=as.numeric(cumsum(c(5.5,11,3,3,5,8))),linetype=2)+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        text = element_text(face="bold"),
        axis.title = element_blank(),
        axis.ticks=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size = 10))
p1
palette1 <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99")
p2=ggplot(markers,aes(x,reorder(y,-y),fill=group))+
  geom_tile()+
  geom_text(aes(label=markers),size=3)+
  scale_fill_manual(values = palette1)+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        text=element_text(),
        axis.text = element_blank(),
        axis.title =element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'none')+
  scale_x_continuous(expand = c(0,0))
color_labels <- c("Tumor vascularization", "Morphogenesis",'Inflammation', "MHC-I", "MHC-II", "Guidance receptors",
                  "Endothelial Adhesion Molecules", "Immunosuppressive Molecules")
#color_labels <- c("Endothelial Adhesion Molecules", "Immunosuppressive Molecules")
label_data <- data.frame(y = color_labels)
label_data$y=factor(label_data$y,levels = c("Tumor vascularization", "Morphogenesis",'Inflammation', "MHC-I", "MHC-II", "Guidance receptors",
                                            "Endothelial Adhesion Molecules", "Immunosuppressive Molecules"))
#label_data$y=factor(label_data$y,levels = c("Endothelial Adhesion Molecules", "Immunosuppressive Molecules"))
p2_color_annotations <- ggplot(label_data, aes(x = 0, y = y)) +
  geom_tile(aes(fill = y), color = "white") +
  scale_fill_manual(values = palette1) +
  geom_text(aes(label = y), hjust = 0, size = 4, color = "black") +
  theme_void()
p2_color_annotations
heatmap_with_color_annotations = p2_color_annotations +p2 +p1 +
  plot_layout(ncol = 2, widths  = c(1, 3))
heatmap_with_color_annotations & theme(plot.margin = margin(0,0,0,0))
library(patchwork)
heatmap1 =p2+p1+plot_layout(ncol = 2, widths  = c(1, 3))
heatmap1 & theme(plot.margin = margin(0,0,0,0))
DotPlot(sce,features = gene,
        group.by = "location",
)+theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradientn(values = seq(0,1,0.2),colors = c('#330066','#336699','#66CC66','#FFCC33'))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))
DimPlot(fib,group.by = "subcelltype",cols = col_cluster,raster = F,label = F,split.by = "location")
gene <- c('PDCD1','LAG3','HAVCR2','CCL5','CXCL13')
# Set color palette
# cividis：option E;
library(scCustomize)
library(viridis)
pal <- viridis(n = 10,option="C")
pal <- viridis(n = 15, option = "D", direction = -1)
FeaturePlot(object = sce,features = meso, cols = pal, order = T)
#viridis_plasma_dark_high
#viridis_plasma_light_high
#viridis_magma_light_high
#viridis_inferno_dark_high
#viridis_inferno_light_high
#viridis_dark_high
i=1
plots=list()
for (i in 1:length(gene)){
  plots[[i]]=FeaturePlot_scCustom(seurat_object = sce,
                                  colors_use = viridis_dark_high,
                                  features = gene[i])+NoLegend()+NoAxes()+
    theme(panel.border = element_rect(fill = NA,color = "black",
                                      size=1.5,linetype = "solid"))
}
library(patchwork)
p<-wrap_plots(plots, ncol = 3);p
library(Seurat)
library(ggplot2)
library(tidyverse)
markers <-c("CD8A","CD8B","CD4", #CD8/CD4
            "FCGR3A","KLRD1", #NK
            "FOXP3","IL2RA","IKZF2", #Treg
            "TCF7","SELL","LEF1","CCR7", #Naive
            "LAG3","TIGIT","PDCD1","HAVCR2","CTLA4","LAYN","ENTPD1", #EXhausted
            "GZMA","GZMB","GZMK","GNLY","IFNG","PRF1","NKG7", #Cytotoxic
            "CD28","ICOS","CD40LG","TNFRSF4","TNFRSF9","TNFRSF18", #Co-stimulatory
            "CD69","RUNX3","NR4A1", #TRM
            "CCL3","CCL5","CXCL13","IL21", #Chemokines
            "EOMES","MAF","TOX","ID2","TBX21","HOPX","ZNF683", #TFs
            "MKI67","STMN1" #Proliferating
)
markers <- as.data.frame(markers)
markers$x <- 1
markers$y <- seq(1:nrow(markers))
markers$group <- c(rep('CD8/CD4',3),rep('NK',2),rep('Treg',3),rep('Naive',4),
                   rep("Exhausted",7),rep("Cytotoxic",7),rep('Co-stimulatory',6),
                   rep("TRM",3),rep("Chemokines",4),rep("TFs",7),rep("Proliferating",2))
marker_ave_exp <- AverageExpression(sce,assays = "RNA",
                                    features = markers$markers,
                                    group.by = "subcelltype",
                                    slot="data")
markers$markers=factor(markers$markers,levels = markers$markers)
data=base::apply(marker_ave_exp$RNA,1,
                 function(x) (x-mean(x))/sd(x))%>%t()%>%as.data.frame()%>%rownames_to_column('Gene')%>%reshape2::melt()
data$Gene=factor(data$Gene,levels = rev(unique(data$Gene)))
data$variable = factor(data$variable, levels = c('CD4_c1_FOS',   'CD4_c2_SELL',  'CD4_c3_FOXP3', 'CD4_c4_CXCL13',   'CD4_c5_IL7R',
                                                 'CD8_c1_IFNG', 'CD8_c2_SLC7A5',   'CD8_c3_GZMK', 'CD8_c4_TLE5', 'CD8_c5_CXCL13',
                                                 'CD8_c6_ISG15',
                                                 'NK_c1_FCGR3A',  'NK_c2_FCER1G',    'NK_c3_AREG',
                                                 'NKT',           'γδT'))
#color_palette <- c("#7F3F00", "#B35806", "#E08214", "#FDAE61", "#FEE0B6", "#D8DAEB", "#B2ABD2", "#8073AC", "#542788", "#2D004B")
library(RColorBrewer)
color_palette <- colorRampPalette(brewer.pal(11, "RdYlBu"))(100)
p1=ggplot(data,aes(x = variable,y = Gene,fill=value))+
  geom_tile()+
  scale_y_discrete(expand = c(0,0))+
  scale_fill_gradientn(colors=rev(colorRampPalette(color_palette)(500)),
                       limits=c(-2,2),name="Z Score",
                       oob = scales::squish
  )+
  geom_vline(xintercept=as.numeric(cumsum(table(unique(data$variable)))+0.5),linetype=2)+
  geom_hline(yintercept=as.numeric(cumsum(c(2.5,7,4,3,6,7,7,4,3,2,3))),linetype=2)+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        text = element_text(face="bold"),
        axis.title = element_blank(),
        axis.ticks=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size = 10))
p1
palette1 <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99")
p2=ggplot(markers,aes(x,reorder(y,-y),fill=group))+
  geom_tile()+
  geom_text(aes(label=markers),size=3)+
  scale_fill_manual(values = palette1)+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        text=element_text(),
        axis.text = element_blank(),
        axis.title =element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'none')+
  scale_x_continuous(expand = c(0,0))
color_labels <- c("CD8/CD4", "NK", "Treg", "Naïve", "EXhausted",
                  "Cytotoxic", "Co-stimulatory", "TRM", "Chemokines",
                  "TFs", "Proliferating")
label_data <- data.frame(y = color_labels)
label_data$y=factor(label_data$y,levels = c("CD8/CD4", "NK", "Treg", "Naïve", "EXhausted",
                                            "Cytotoxic", "Co-stimulatory", "TRM", "Chemokines",
                                            "TFs", "Proliferating"))
p2_color_annotations <- ggplot(label_data, aes(x = 0, y = y)) +
  geom_tile(aes(fill = y), color = "white") +
  scale_fill_manual(values = palette1) +
  geom_text(aes(label = y), hjust = 0, size = 4, color = "black") +
  theme_void()
p2_color_annotations
heatmap_with_color_annotations = p2_color_annotations +p2 +p1 +
  plot_layout(ncol = 2, widths  = c(1, 3))
heatmap_with_color_annotations & theme(plot.margin = margin(0,0,0,0))
library(patchwork)
heatmap1 =p2+p1+plot_layout(ncol = 2, widths  = c(1, 3))
heatmap1 & theme(plot.margin = margin(0,0,0,0))
t = sce
sce1 = t[, sample(1:ncol(t),round(ncol(t)/100)) ]
library(CytoTRACE2) #loading
results <- cytotrace2(sce,
                      species = "human",
                      is_seurat = TRUE,
                      full_model = FALSE,
                      slot_type = "counts"
)
annotation <- data.frame(phenotype = sce@meta.data$subcelltype) %>%
  set_rownames(., colnames(sce))
# plotting
plots <- plotData(cytotrace2_result = results,
                  annotation = annotation,
                  is_seurat = TRUE)
emb_1 <- sce@reductions$umap@cell.embeddings[,1]+3 # First dimension of your "harmony" reduction embedding
emb_2 <- sce@reductions$umap@cell.embeddings[,2] # Second dimension of your "harmony" reduction embedding
# replace CytoTRACE2_UMAP embeddings
plots[["CytoTRACE2_UMAP"]][[1]][["data"]]["UMAP_1"] <- emb_1
plots[["CytoTRACE2_UMAP"]][[1]][["data"]]["UMAP_2"]<- emb_2
plots[["CytoTRACE2_Potency_UMAP"]][[1]][["data"]]["UMAP_1"] <- emb_1
plots[["CytoTRACE2_Potency_UMAP"]][[1]][["data"]]["UMAP_2"]<- emb_2
plots[["CytoTRACE2_Relative_UMAP"]][[1]][["data"]]["UMAP_1"] <- emb_1
plots[["CytoTRACE2_Relative_UMAP"]][[1]][["data"]]["UMAP_2"]<- emb_2
plots[["Phenotype_UMAP"]][[1]][["data"]]["UMAP_1"] <- emb_1
plots[["Phenotype_UMAP"]][[1]][["data"]]["UMAP_2"]<- emb_2
plots$CytoTRACE2_UMAP
plots$CytoTRACE2_Potency_UMAP
plots$CytoTRACE2_Relative_UMAP
plots$Phenotype_UMAP
plots$CytoTRACE2_Boxplot_byPheno
plots$CytoTRACE2_Relative_UMAP <- plots$CytoTRACE2_Relative_UMAP +
  scale_color_viridis_c(
    option = "viridis",
    name = "Relative Score"
  )
DimPlot(sce,group.by = "subcelltype",cols = col_cluster,raster = TRUE,label = TRUE)
chemo_r = c("CCR1","CCR2","CCR3","CCR4","CCR5","CCR6","CCR7","CCR8","CCR9",
            "CCR10","CXCR1","CXCR2","CXCR3","CXCR3B","CXCR4","CXCR5","CXCR6",
            "CXCR7","XCR1","CX3CR1")
DotPlot(all,features = "LYVE1",
        group.by = "celltype",
)+theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradientn(values = seq(0,1,0.2),colors = c('#330066','#336699','#66CC66','#FFCC33'))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))
FeaturePlot(sce,features = c("CCR1","CCR5","CXCR4"))
markers <-c("HLA-A","HLA-B","HLA-C", #MHC-I
            "HLA-DQB1","HLA-DRB1",
            "HLA-DQA1","HLA-DRA","HLA-DMB",
            "HLA-E","HLA-DOA", #MHC-II
            "CD274","PDCD1","HAVCR2","TIGIT","VSIR","CTLA4",
            "SIGLEC1","LAG3", #Immunocheckpoint
            "CCL2","CCL13","CCL18","CXCL1","CCL5","CXCL9",
            "CXCL10","CXCL11",
            "CCL8","CXCL2","CCL3","CCL4","CCL4L2", #Chemokine
            "IL1A","IL6","CSF1","MMP9","MMP14","TNF",
            "TGFB1","TGFB2" #Cytokine
)
markers <- as.data.frame(markers)
markers$x <- 1
markers$y <- seq(1:nrow(markers))
markers$group <- c(rep('MHC-I',3),rep('MHC-II',7),rep('Immunocheckpoint',8),rep('Chemokine',13),
                   rep("Cytokine",8))
marker_ave_exp <- AverageExpression(sce,assays = "RNA",
                                    features = markers$markers,
                                    group.by = "subcelltype",
                                    slot="data")
markers$markers=factor(markers$markers,levels = markers$markers)
data=base::apply(marker_ave_exp$RNA,1,
                 function(x) (x-mean(x))/sd(x))%>%t()%>%as.data.frame()%>%rownames_to_column('Gene')%>%reshape2::melt()
data$Gene=factor(data$Gene,levels = rev(unique(data$Gene)))
data$variable = factor(data$variable, levels = c( "Mono_c1_CCL20",
                                                  "Mono_c2_S100A8",
                                                  "Mono_c3_LILRA5",
                                                  "Macro_c1_CCL4",
                                                  "Macro_c2_C1QC",
                                                  "Macro_c3_CRIP1",
                                                  "Macro_c4_SPP1",
                                                  "Macro_c5_CD1C",
                                                  "Macro_c6_LYVE1/TIMD4",
                                                  "Macro_c7_LYVE1",
                                                  "Macro_c8_CXCL10",
                                                  "DC_c1_CD1C",
                                                  "DC_c2_CLEC9A",
                                                  "DC_c3_LAMP3",
                                                  "Neutrophils",
                                                  "Mast"))
#color_palette <- c("#7F3F00", "#B35806", "#E08214", "#FDAE61", "#FEE0B6", "#D8DAEB", "#B2ABD2", "#8073AC", "#542788", "#2D004B")
library(RColorBrewer)
color_palette <- colorRampPalette(brewer.pal(11, "RdYlBu"))(100)
p1=ggplot(data,aes(x = variable,y = Gene,fill=value))+
  geom_tile()+
  scale_y_discrete(expand = c(0,0))+
  scale_fill_gradientn(colors=rev(colorRampPalette(color_palette)(500)),
                       limits=c(-2,2),name="Z Score",
                       oob = scales::squish
  )+
  geom_vline(xintercept=as.numeric(cumsum(table(unique(data$variable)))+0.5),linetype=2)+
  geom_hline(yintercept=as.numeric(cumsum(c(8.5,13,8,7,3))),linetype=2)+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        text = element_text(face="bold"),
        axis.title = element_blank(),
        axis.ticks=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size = 10))
p1
palette1 <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99")
p2=ggplot(markers,aes(x,reorder(y,-y),fill=group))+
  geom_tile()+
  geom_text(aes(label=markers),size=3)+
  scale_fill_manual(values = palette1)+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        text=element_text(),
        axis.text = element_blank(),
        axis.title =element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'none')+
  scale_x_continuous(expand = c(0,0))
color_labels <- c("MHC-I", "MHC-II", "Immunocheckpoint", "Chemokine", "Cytokine")
label_data <- data.frame(y = color_labels)
label_data$y=factor(label_data$y,levels = c("MHC-I", "MHC-II", "Immunocheckpoint", "Chemokine", "Cytokine"))
p2_color_annotations <- ggplot(label_data, aes(x = 0, y = y)) +
  geom_tile(aes(fill = y), color = "white") +
  scale_fill_manual(values = palette1) +
  geom_text(aes(label = y), hjust = 0, size = 4, color = "black") +
  theme_void()
p2_color_annotations
heatmap_with_color_annotations = p2_color_annotations +p2 +p1 +
  plot_layout(ncol = 2, widths  = c(1, 3))
heatmap_with_color_annotations & theme(plot.margin = margin(0,0,0,0))
library(patchwork)
heatmap1 =p2+p1+plot_layout(ncol = 2, widths  = c(1, 3))
heatmap1 & theme(plot.margin = margin(0,0,0,0))
different_express_gene<- FindMarkers(sce,
                                     group.by = "subcelltype",
                                     logfc.threshold = 0.25,
                                     test.use = "wilcox",
                                     ident.1 = "Macro_c6_LYVE1/TIMD4",
                                     ident.2 = "Macro_c7_LYVE1")
different_express_gene$gene <- rownames(different_express_gene)
library(stringr)
library(EnhancedVolcano)
geneID = c("CRIP1","LYPD2", "S100A8", "FN1")
EnhancedVolcano(different_express_gene,
                lab =different_express_gene$gene,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 4.0,
                selectLab = c("CRIP1","LYPD2", "S100A8", "FN1","TIMP1","S100A10",
                              "S100A6","FCN1","S100A4",
                              "CCL4L2","CCL3","CXCL8","CCL4","EGR1","MRC1"),
                drawConnectors = TRUE,
                title = 'Macro_c6_LYVE1/TIMD4 vs Macro_c7_LYVE1')
DimPlot(sce,group)
sce = sce[, sample(1:ncol(sce),round(ncol(sce)/20)) ]
count_matrix <- sce@assays$RNA@counts
write.csv(count_matrix, file='./count_matrix.csv', row.names = T)
pyr <- c("CAD","DHODH","UMPS","TYMS",
         "UCK1","UCK2","TK1","TK2","CDA","UCKL1","DCK",
         "SLC29A1","SLC29A2","SLC29A3","SLC29A4",
         "SLC28A1","SLC28A2","SLC28A3"
)
DotPlot(epi, features = pyr,cols = c( "#00C1D4","#1E3163"))+ coord_flip()
gene_cell_exp <- AverageExpression(sce,
                                   features = anergy,
                                   group.by = 'subcelltype',
                                   slot = 'data')
gene_cell_exp <- as.data.frame(gene_cell_exp$RNA)
library(ComplexHeatmap)
df <- data.frame(colnames(gene_cell_exp))
colnames(df) <- 'class'
library(ggsci)
mypal <- pal_npg()(6)
mypal
library(scales)
show_col(mypal)
top_anno = HeatmapAnnotation(df = df,
                             border = T,
                             show_annotation_name = F,
                             gp = gpar(col = 'black'),
                             col = list(class = c('Epi_c1_GKN1'= "#E64B35FF",
                                                  'Epi_c2_TM4SF1'="#4DBBD5FF",
                                                  'Epi_c3_PGA3'="#00A087FF",
                                                  'Epi_c4_APOA1' = "#3C5488FF",
                                                  'Epi_c5_GHRL'= "#F39B7FFF",
                                                  'Epi_c6_GIF' = "#8491B4FF")
                                        ))
top_anno = HeatmapAnnotation(df = df,
                             border = T,
                             show_annotation_name = F,
                             gp = gpar(col = 'black'),
                             col = list(class = c("N" = "#2CA02CFF",
                                                  "PT" = "#FF7F0EFF",
                                                  "Ascites" = "#1F77B4FF",
                                                  "PM" = "#D62728FF"
                                                  )
                             ))
marker_exp <- t(scale(t(gene_cell_exp),scale = T,center = T))
library(RColorBrewer)
coul <- colorRampPalette(brewer.pal(9, "OrRd"))(50)
Heatmap(marker_exp,
        cluster_rows = F,
        cluster_columns = F,
        show_column_names = T,
        show_row_names = T,
        column_title = NULL,
        heatmap_legend_param = list(
          title=' '),
        col = coul,
        border = 'black',
        rect_gp = gpar(col = "black", lwd = 1),
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10))
#,
#       top_annotation = top_anno)
pm_samples <- c(
  '1422329_PT',  '1422573_PT',  '1425236_PT',  '1460016_PT',  '1483231_PT',  '1492435_PT',
  'HRR1220993/', 'HRR1220997/', 'HRR1221000/', 'HRR1221003/', 'HRR1221006/', 'HRR1221010/',
  'HRR1221013/', 'HRR1221015/', 'HRR1221022/', 'HRR1221025/', 'sample36')
sce$group <- ifelse(sce$orig.ident %in% pm_samples, "PM", "nonPM")
different_express_gene<- FindMarkers(sce,
                                     group.by = "subcelltype",
                                     logfc.threshold = 0.25,
                                     test.use = "wilcox",
                                     ident.1 = "Macro_c6_LYVE1/TIMD4",
                                     ident.2 = "Macro_c7_LYVE1")
different_express_gene$gene <- rownames(different_express_gene)
different_express_gene3<- FindMarkers(lyve1,
                                     group.by = "location",
                                     logfc.threshold = 0.25,
                                     test.use = "wilcox",
                                     ident.1 = "Ascites",
                                     ident.2 = "PM")
different_express_gene1$gene <- rownames(different_express_gene1)
myeloid_subsets <- c(
  "DC_c1_CD1C", "DC_c2_CLEC9A", "DC_c3_LAMP3",
  "Macro_c1_CCL4", "Macro_c2_C1QC", "Macro_c3_CRIP1",
  "Macro_c4_SPP1", "Macro_c5_CD1C", "Macro_c6_LYVE1/TIMD4",
  "Macro_c7_LYVE1", "Macro_c8_CXCL10",
  "Mast",
  "Mono_c1_CCL20", "Mono_c2_S100A8", "Mono_c3_LILRA5",
  "Neutrophils"
)
sce$celltype[sce$celltype %in% myeloid_subsets] <- "Myeloid_cells"
locations <- unique(sce$location)
cellchat_list <- list()
for(loc in locations) {
  cat("\nProcessing location:", loc, "\n")
  cells_loc <- which(sce$location == loc)
  data.input <- GetAssayData(sce[, cells_loc], slot = "data", assay = "RNA")
  meta <- data.frame(labels = sce$celltype[cells_loc],
                     row.names = colnames(data.input))
  meta$labels <- droplevels(factor(meta$labels))
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
  CellChatDB <- CellChatDB.human
  cellchat@DB <- CellChatDB
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat_list[[loc]] <- cellchat
  cat("Finished location:", loc, "\n")
}
save(cellchat_list, file = "mesocellchat_object.list.RData")
cell_types_by_location <- list()
for(loc in locations) {
  cellchat_obj <- cellchat_list[[loc]]
  cell_types_by_location[[loc]] <- levels(cellchat_obj@idents)
  cat(loc, " cell types:", cell_types_by_location[[loc]], "\n")
}
all_cell_types <- unique(unlist(cell_types_by_location))
cat("All cell types across locations:\n", all_cell_types, "\n")
common_cell_types <- all_cell_types
for(loc in locations) {
  cellchat_list[[loc]] <- liftCellChat(cellchat_list[[loc]],
                                       group.new = common_cell_types)
  cat("Unified cell-type identifiers for", loc, "\n")
}
cellchat_compare <- mergeCellChat(
  list(cellchat_list[["PM"]],cellchat_list[["Ascites"]]),
  add.names = c("PM","Ascites")
)
gg1 <- compareInteractions(cellchat_compare,
                           show.legend = FALSE,
                           group = c(1, 2))
gg2 <- compareInteractions(cellchat_compare,
                           show.legend = FALSE,
                           group = c(1, 2),
                           measure = "weight")
gg1 + gg2
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat_compare, weight.scale = T)
netVisual_diffInteraction(cellchat_compare, weight.scale = T, measure = "weight")
gg1 <- netVisual_heatmap(cellchat_compare)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat_compare, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2
netVisual_bubble(cellchat_compare,  targets.use = "Myeloid_cells",  comparison = c(1, 2), angle.x = 45)
pro = read.csv("proteinome.csv")
head(pro)
summary(pro)
library(VIM)
aggr(pro)
library(missForest)
pro$PG.Genes <- as.factor(pro$PG.Genes)
dat2 = missForest(pro,ntree = 100)
library(limma)
library(tidyverse)
library(openxlsx)
library(ggpubr)
library(ggthemes)
exp <- dat2$ximp
base_names <- exp[, 1]
rownames(exp) <- make.unique(as.character(base_names))
exp <- exp[, -1, drop = FALSE]
head(rownames(exp))
sample_names <- colnames(exp)
meta <- data.frame(
  sample_id = sample_names,
  Type = ifelse(grepl("neg", sample_names), "neg", "pos"),
  stringsAsFactors = FALSE
)
rownames(meta) <- sample_names
identical(rownames(meta), colnames(exp)) # check names
exp_matrix <- as.matrix(exp)
if(max(exp_matrix, na.rm = TRUE) > 50) {
  exp_matrix <- log2(exp_matrix + 1)
  print("Data were log2-transformed")
}
meta <- meta %>% mutate(contrast = as.factor(Type))
design <- model.matrix(~ 0 + contrast, data = meta)
colnames(design) <- c("neg", "pos")
fit <- lmFit(exp_matrix, design)
contrast_matrix <- makeContrasts(
  pos_vs_neg = pos - neg,
  levels = design
)
fits <- contrasts.fit(fit, contrast_matrix)
ebFit <- eBayes(fits)
limma.res <- topTable(ebFit, coef = "pos_vs_neg",
                      adjust.method = 'fdr', number = Inf)
limma.res$ID <- rownames(limma.res)
limma.res <- limma.res %>%
  filter(!is.na(P.Value)) %>%
  mutate(logP = -log10(P.Value)) %>%
  mutate(tag = "pos vs neg") %>%
  mutate(Gene = ID)
limma.res <- limma.res %>%
  mutate(group = case_when(
    (P.Value < 0.05 & logFC > 0.58) ~ "up",
    (P.Value  < 0.05 & logFC < -0.58) ~ "down",
    .default = "not sig"
  ))
result_summary <- table(limma.res$group)
print("Differential expression summary:")
print(result_summary)
write.xlsx(limma.res, "PRO_Limma_pos_vs_neg.xlsx",
           overwrite = TRUE, rowNames = FALSE)
limma.res <- limma.res %>%
  mutate(group = factor(group, levels = c("up", "down", "not sig")))
my_label <- paste0("FC:1.5 ; AdjP:0.05 ; ",
                   "Up:", ifelse("up" %in% names(result_summary), result_summary["up"], 0), " ; ",
                   "Down:", ifelse("down" %in% names(result_summary), result_summary["down"], 0))
p <- ggscatter(limma.res,
               x = "logFC", y = "logP",
               color = "group", size = 2,
               main = "Treatment vs Control",
               xlab = "log2FoldChange",
               ylab = "-log10(P.value)",
               palette = c("#D01910", "#00599F", "#CCCCCC"),
               ylim = c(-1, 10),xlim=c(-3,3)) +
  theme_base() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#222222") +
  geom_vline(xintercept = 0.58, linetype = "dashed", color = "#222222") +
  geom_vline(xintercept = -0.58, linetype = "dashed", color = "#222222") +
  labs(subtitle = my_label) +
  theme(plot.background = element_blank())
target_genes <- c("CTSG", "PADI4", "KRT6B", "MPO",
                  "CTSB","LYVE1","CD163","FOLR2","IL4I1", "GPNMB","MRC1","VSIG4")
sig_genes <- limma.res[limma.res$Gene %in% target_genes, ]
library(ggrepel)
p_annotated <- p +
  ggrepel::geom_text_repel(
    data = sig_genes,
    aes(x = logFC, y = logP, label = Gene),
    size = 4,
    color = "black",
    box.padding = 0.5,
    segment.color = "grey50",
    segment.size = 0.3,
    max.overlaps = 20
  )
p_annotated
ggsave("diff_protein.tiff", p, width = 10, height = 10, dpi = 300)
marker11$gene = rownames(marker11)
EnhancedVolcano(marker11,
                lab =marker11$gene,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 4.0,
                selectLab = c("PLAT","SERPINE1","CFB","CCN1","SAA1","IGFBP3","WFDC2","NME2","CEMIP",
                              "ITLN1","ADIRF","OGN","GNB2L1"),
                drawConnectors = TRUE,
                title = 'Meso_c2_SAA1 vs Meso_c1_ITLN1')
f = dir(file.path("data", "public", "GSE163558"))
s = c("nonPM","PM","nonPM")
t = c("N","PM","PT")
scelist = list()
for(i in 1:length(f)){
  tmp <- Read10X(paste0(f[[i]]))
  tmp <- CreateSeuratObject(counts = tmp, project = paste0(f[[i]]), min.cells = 3,
                            min.features = 200)
  tmp[["percent.mt"]] <- PercentageFeatureSet(tmp, pattern = "^MT-")
  tmp@meta.data$location  <- paste0(t[[i]])
  tmp@meta.data$dataset  <- c("GSE163558")
  tmp@meta.data$group <- paste0(s[[i]])
  scelist[i]=tmp
}
sce.all <- merge(scelist[[1]],
                 y= scelist[ -1 ])
sce.all <- subset(sce.all, subset = nFeature_RNA > 500 & nFeature_RNA < 4500 & percent.mt < 20)
VlnPlot(sce.all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "orig.ident",ncol = 3,pt.size = 0)
sce.all <- NormalizeData(sce.all)
sce.all <- FindVariableFeatures(sce.all)
sce.all <- ScaleData(sce.all)
sce.all <- RunPCA(sce.all, features = VariableFeatures(sce.all))
library(harmony)
sce.all <- RunHarmony(sce.all, "orig.ident")
names(sce.all@reductions)
sce.all <- RunUMAP(sce.all,  dims = 1:20,
                   reduction = "harmony")
sce.all <- FindNeighbors(sce.all,reduction = "harmony",dims = 1:20)
sce.all <- FindClusters(sce.all, resolution = 0.1)
DimPlot(sce.all)
saveRDS(sce.all, "GSE163558.rds")
# setwd(file.path("data", "public", "GSE228598"))
sce = qread("sce.qs")
table(sce$orig.ident)
library(stringr)
extract_gsm <- function(orig_ident_path) {
  gsm_id <- str_extract(orig_ident_path, "GSM\\d+")
  return(gsm_id)
}
sce$GSM_ID <- sapply(sce$orig.ident, extract_gsm)
head(sce$GSM_ID)
table(sce$GSM_ID)
meta_mapping <- data.frame(
  GSM_ID = c(
    'GSM7133741', 'GSM7133742', 'GSM7133743', 'GSM7133744', 'GSM7133745',
    'GSM7133746', 'GSM7133747', 'GSM7133748', 'GSM7133749', 'GSM7133750',
    'GSM7133751', 'GSM7133752', 'GSM7133753', 'GSM7133754', 'GSM7133755',
    'GSM7133756', 'GSM7133757', 'GSM7133758', 'GSM7133759', 'GSM7133760',
    'GSM7133761', 'GSM7133762', 'GSM7133763', 'GSM7133764', 'GSM7133765',
    'GSM7133766', 'GSM7133767', 'GSM7133768'
  ),
  Sample_Type = c(
    "Ascites", "PLF", "PLF", "PLF", "Ascites",
    "Ascites", "PLF", "Ascites", "PLF", "Ascites",
    "PLF","PLF","Ascites","PLF","PLF",
    "Ascites","PLF","PLF","PLF","PLF",
    "Ascites","Ascites","Ascites","Ascites","Ascites",
    "Ascites","Ascites","Ascites"
  ),
  stringsAsFactors = FALSE
)
head(meta_mapping)
sce$location <- meta_mapping$Sample_Type[match(sce$GSM_ID, meta_mapping$GSM_ID)]
head(sce$location)
table(sce$location)
sce = sce[,sce$location %in% c('Ad', 'Pt', 'As', 'PBMC')]
table(sce$location)
sce = sce[,sce$patient %in% c('MDA_Pt1', 'MDA_Pt2', 'MDA_Pt3',
                              'MDA_Pt4')]
table(sce$patient)
meta_mapping <- data.frame(
  location = c('Ad', 'Pt', 'As', 'PBMC'
  ),
  Sample_Type = c(
    "N","PT","Ascites","PBMC"
  ),
  stringsAsFactors = FALSE
)
sce$location <- meta_mapping$Sample_Type[match(sce$location, meta_mapping$location)]
sce$dataset = "GSE234129"
GSE234129 = sce
library(Seurat)
data_dir <- file.path("data", "public", "GSE239676")
readLines(gzfile("features.tsv.gz"), n = 5)
library(Seurat)
library(Matrix)
# ----------------------
# ----------------------
barcodes <- read.table(
  file = gzfile(file.path(data_dir, "barcodes.tsv.gz")),
  stringsAsFactors = FALSE,
  header = FALSE,
  col.names = "barcode"
)$barcode
# ----------------------
# ----------------------
gene_names <- readLines(gzfile(file.path(data_dir, "features.tsv.gz")))
gene_names <- gene_names[gene_names != ""]
# ----------------------
# ----------------------
count_matrix <- readMM(
  file = gzfile(file.path(data_dir, "matrix.mtx.gz"))
)
# ----------------------
# ----------------------
rownames(count_matrix) <- gene_names
colnames(count_matrix) <- barcodes
# ----------------------
# ----------------------
count_matrix_sparse <- as(count_matrix, "dgCMatrix")
# ----------------------
# ----------------------
seurat_obj <- CreateSeuratObject(
  counts = count_matrix_sparse,
  project = "GSE239676",
  min.cells = 3,
  min.features = 200
)
meta_data <- read.delim("GSE239676_meta.tsv.gz", sep = "\t", header = TRUE)
seurat_obj@meta.data <- cbind(
  seurat_obj@meta.data,
  meta_data[match(seurat_obj$orig.ident, meta_data$Patient), ]
)
table(seurat_obj$Tissue)
meta_mapping <- data.frame(
  location = c('Ad',  'As', 'Ascites',  'GC_PLF',   'N',      'NP',     'PBMC',     'PLF',
               'PM',      'Pt',      'PT'
  ),
  Sample_Type = c(
    'N',  'Ascites', 'Ascites',  'GC_PLF',   'N',      'NP',     'PBMC',     'PLF',
    'PM',      'PT',      'PT'
  ),
  stringsAsFactors = FALSE
)
sce$location <- meta_mapping$Sample_Type[match(sce$location, meta_mapping$location)]
emb <- Embeddings(sce, "harmony")
dim(emb)
range(emb, na.rm = TRUE)
sce <- RunHarmony(sce, group.by.vars = "sample_group_column")
