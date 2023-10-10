####mouse_brain_annotation_version1.h5Seurat是最初用来做分析的版本
###mouse_brain_annotation_version2.h5Seurat重新聚类cluster14后的版本
#####umap图和小提琴图
##

setwd("/jdfssz1/ST_SUPERCELLS/P20Z10200N0059/tianyingming/project/hanlei_zju/scRNA/00.data")
library(SeuratDisk)
library(Seurat)
library(dplyr)
library(ggplot2)
sc <- LoadH5Seurat("/jdfssz1/ST_SUPERCELLS/P20Z10200N0059/tianyingming/project/hanlei_zju/scRNA/00.data/mouse_brain_annotation_use.h5Seurat")
cell_order = c("Excitatory neuronal","Cholinergic inhibitory neuron","Cnr1+ inhibitory neuron","Drd1+ inhibitory neuron","Drd2+ inhibitory neuron","Npy inhibitory neuron","PV+ inhibitory neuron","SST+ inhibitory neuron","VIP+ inhibitory neuron","Zic1+ inhibitory neuron","Astrocyte","Endothelial cell","Fibroblast-like cell","Microglia","Oligodendrocyte","Oligodendrocyte precursor")
sc$celltype = ""
sc@meta.data[sc@meta.data$seurat_clusters %in% c(8,17),]$celltype <- "Astrocyte"
sc@meta.data[sc@meta.data$seurat_clusters %in% c(22),]$celltype <- "Endothelial cell"
sc@meta.data[sc@meta.data$seurat_clusters %in% c(24),]$celltype <- "Fibroblast-like cell"
sc@meta.data[sc@meta.data$seurat_clusters %in% c(14),]$celltype <- "Microglia"
sc@meta.data[sc@meta.data$seurat_clusters %in% c(16),]$celltype <- "Oligodendrocyte"
sc@meta.data[sc@meta.data$seurat_clusters %in% c(19),]$celltype <- "Oligodendrocyte precursor"
sc@meta.data[sc@meta.data$seurat_clusters %in% c(10,23),]$celltype <- "SST+ inhibitory neuron"
sc@meta.data[sc@meta.data$seurat_clusters %in% c(13),]$celltype <- "VIP+ inhibitory neuron"
sc@meta.data[sc@meta.data$seurat_clusters %in% c(6),]$celltype <- "PV+ inhibitory neuron"
sc@meta.data[sc@meta.data$seurat_clusters %in% c(11,18),]$celltype <- "Cnr1+ inhibitory neuron"
sc@meta.data[sc@meta.data$seurat_clusters %in% c(3),]$celltype <- "Drd1+ inhibitory neuron"
sc@meta.data[sc@meta.data$seurat_clusters %in% c(9),]$celltype <- "Drd2+ inhibitory neuron"
sc@meta.data[sc@meta.data$seurat_clusters %in% c(5),]$celltype <- "Zic1+ inhibitory neuron"
sc@meta.data[sc@meta.data$seurat_clusters %in% c(21),]$celltype <- "Cholinergic inhibitory neuron"
sc@meta.data[sc@meta.data$seurat_clusters %in% c(15),]$celltype <- "Npy inhibitory neuron"
sc@meta.data[sc@meta.data$seurat_clusters %in% c(0,1,2,7,12,4,11,20),]$celltype <- "Excitatory neuronal"
SaveH5Seurat(sc, filename = "/jdfssz1/ST_SUPERCELLS/P20Z10200N0059/tianyingming/project/hanlei_zju/scRNA/00.data/mouse_brain_annotation_version1.h5Seurat",overwrite = T)

sc_choose <- subset(sc,subset= seurat_clusters %in% seq(0,13,1)|seurat_clusters %in% seq(15,24,1)|celltype=="Microglia")
#sc_choose <- subset(sc,subset= seurat_clusters %in% c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,15,16,17,18,19,20,21,22,23,24)|celltype=="Microglia")
SaveH5Seurat(sc_choose, filename = "/jdfssz1/ST_SUPERCELLS/P20Z10200N0059/tianyingming/project/hanlei_zju/scRNA/00.data/mouse_brain_annotation_version2.h5Seurat",overwrite = T)


##直接调color调颜色顺序，细胞顺序
sc <- LoadH5Seurat("/jdfssz1/ST_SUPERCELLS/P20Z10200N0059/tianyingming/project/hanlei_zju/scRNA/00.data/mouse_brain_annotation_version2.h5Seurat")
colors <- c(`Excitatory neuronal`="#FEAF16",`Cholinergic inhibitory neuron`="#1CFFCE",`Cnr1+ inhibitory neuron`="#90AD1C",`Drd1+ inhibitory neuron`="#2ED9FF",`Drd2+ inhibitory neuron`="#DEA0FD",`Npy inhibitory neuron`="#AA0DFE",`PV+ inhibitory neuron`="#F8A19F",`SST+ inhibitory neuron`="#325A9B",`VIP+ inhibitory neuron`="#1C8356",`Zic1+ inhibitory neuron`="#C4451C",Astrocyte="#5A5156",`Endothelial cell`="#E4E1E3",`Fibroblast-like cell`="#F6222E",Microglia="#FE00FA",Oligodendrocyte="#16FF32",`Oligodendrocyte precursor`="#3283FE")
DimPlot(sc,reduction = "umap",label=T,group.by='celltype',label.size = 2.5,repel =T)+ 
	scale_colour_manual(values=colors)
ggsave("mouse_umap_annotation.pdf",width = 9, height = 6)


Microglia_sc= sc[,sc@meta.data$seurat_clusters %in% c(14)]
Microglia_sc <- NormalizeData(Microglia_sc, normalization.method = "LogNormalize", scale.factor = 1e4)
Microglia_sc <- FindVariableFeatures(Microglia_sc, selection.method = 'vst', nfeatures = 2000)
Microglia_sc <- ScaleData(Microglia_sc, vars.to.regress = "percent.mt")
Microglia_sc <- RunPCA(Microglia_sc, features = VariableFeatures(object = Microglia_sc))
Microglia_sc <- FindNeighbors(Microglia_sc, dims = 1:10)
Microglia_sc <- FindClusters(Microglia_sc, resolution = 1 )
DimPlot(Microglia_sc,reduction = "umap",label=T, group.by ='seurat_clusters')
ggsave("DimPlot_Microglia_umap.pdf",width = 7, height = 6)

sc.markers <- FindAllMarkers(Microglia_sc)   ##专为单细胞设计的MAST(报错)
all.markers <- sc.markers %>% select(gene, everything()) %>% subset(p_val<0.05 & abs(sc.markers$avg_log2FC) > 0.5)
top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.table(top10,"Cluster14_top10_markers.xls",sep="\t", quote=FALSE)

##再聚类
Microglia_sc@meta.data[Microglia_sc@meta.data$seurat_clusters %in% c(0,1,2,4),]$celltype <- "Excitatory neuronal"
Microglia_sc@meta.data[Microglia_sc@meta.data$seurat_clusters %in% c(3),]$celltype <- "Microglia"

cell <- rownames(Microglia_sc@meta.data)
sc@meta.data[cell,]$celltype <- Microglia_sc$celltype
sc@meta.data[cell,]$celltype_large <- Microglia_sc$celltype

SaveH5Seurat(sc, filename = "/jdfssz1/ST_SUPERCELLS/P20Z10200N0059/tianyingming/project/hanlei_zju/scRNA/00.data/mouse_brain_annotation_use.h5Seurat",overwrite = T)

colors <- c(`Excitatory neuronal`="#FEAF16",`Cholinergic inhibitory neuron`="#1CFFCE",`Cnr1+ inhibitory neuron`="#90AD1C",`Drd1+ inhibitory neuron`="#2ED9FF",`Drd2+ inhibitory neuron`="#DEA0FD",`Npy inhibitory neuron`="#AA0DFE",`PV+ inhibitory neuron`="#F8A19F",`SST+ inhibitory neuron`="#325A9B",`VIP+ inhibitory neuron`="#1C8356",`Zic1+ inhibitory neuron`="#C4451C",Astrocyte="#5A5156",`Endothelial cell`="#E4E1E3",`Fibroblast-like cell`="#F6222E",Microglia="#FE00FA",Oligodendrocyte="#16FF32",`Oligodendrocyte precursor`="#3283FE")
DimPlot(sc,reduction = "umap",label=T,group.by='celltype',label.size = 2.5,repel =T)+ 
	scale_colour_manual(values=colors)
ggsave("mouse_umap_annotation.pdf",width = 9, height = 6)


setwd("/jdfssz1/ST_SUPERCELLS/P20Z10200N0059/tianyingming/project/hanlei_zju/scRNA/tmp")
library(Seurat)
library(SeuratDisk) 
library(ggplot2)
cell_order = c("Excitatory neuronal","Cholinergic inhibitory neuron","Cnr1+ inhibitory neuron","Drd1+ inhibitory neuron","Drd2+ inhibitory neuron","Npy inhibitory neuron","PV+ inhibitory neuron","SST+ inhibitory neuron","VIP+ inhibitory neuron","Zic1+ inhibitory neuron","Astrocyte","Endothelial cell","Fibroblast-like cell","Microglia","Oligodendrocyte","Oligodendrocyte precursor")
sc <- LoadH5Seurat("/jdfssz1/ST_SUPERCELLS/P20Z10200N0059/tianyingming/project/hanlei_zju/scRNA/00.data/mouse_brain_annotation_use.h5Seurat") 
colors <- c(`Excitatory neuronal`="#FEAF16",`Cholinergic inhibitory neuron`="#B00068",`Cnr1+ inhibitory neuron`="#90AD1C",`Drd1+ inhibitory neuron`="#2ED9FF",`Drd2+ inhibitory neuron`="#DEA0FD",`Npy inhibitory neuron`="#AA0DFE",`PV+ inhibitory neuron`="#F8A19F",`SST+ inhibitory neuron`="#325A9B",`VIP+ inhibitory neuron`="#1C8356",`Zic1+ inhibitory neuron`="#C4451C",Astrocyte="#5A5156",`Endothelial cell`="#E4E1E3",`Fibroblast-like cell`="#F6222E",Microglia="#FE00FA",Oligodendrocyte="#16FF32",`Oligodendrocyte precursor`="#3283FE")
gene <- read.table("cell_marker_order",header = T,sep="\t",check.names = F)
sc$celltype <- factor(x = sc$celltype, levels = cell_order)
VlnPlot(sc, features = gene$gene, slot = "data",group.by="celltype",stack=T,fill.by = "ident",flip = T,cols = colors)
	labs(x="",y="")
ggsave("VlnPlot_celltype_test.pdf",width=12,height=14)

##检测gdf11
sc$before = ""
sc@meta.data[sc@meta.data$seurat_clusters %in% c(8,17),]$before <- "Astrocyte"
sc@meta.data[sc@meta.data$seurat_clusters %in% c(22),]$before <- "Endothelial cell"
sc@meta.data[sc@meta.data$seurat_clusters %in% c(24),]$before <- "Fibroblast-like cell"
sc@meta.data[sc@meta.data$seurat_clusters %in% c(14),]$before <- "Microglia"
sc@meta.data[sc@meta.data$seurat_clusters %in% c(16),]$before <- "Oligodendrocyte"
sc@meta.data[sc@meta.data$seurat_clusters %in% c(19),]$before <- "Oligodendrocyte precursor"
sc@meta.data[sc@meta.data$seurat_clusters %in% c(10,23),]$before <- "SST+ inhibitory neuron"
sc@meta.data[sc@meta.data$seurat_clusters %in% c(13),]$before <- "VIP+ inhibitory neuron"
sc@meta.data[sc@meta.data$seurat_clusters %in% c(6),]$before <- "PV+ inhibitory neuron"
sc@meta.data[sc@meta.data$seurat_clusters %in% c(11,18),]$before <- "Cnr1+ inhibitory neuron"
sc@meta.data[sc@meta.data$seurat_clusters %in% c(3),]$before <- "Drd1+ inhibitory neuron"
sc@meta.data[sc@meta.data$seurat_clusters %in% c(9),]$before <- "Drd2+ inhibitory neuron"
sc@meta.data[sc@meta.data$seurat_clusters %in% c(5),]$before <- "Zic1+ inhibitory neuron"
sc@meta.data[sc@meta.data$seurat_clusters %in% c(21),]$before <- "Cholinergic inhibitory neuron"
sc@meta.data[sc@meta.data$seurat_clusters %in% c(15),]$before <- "Npy inhibitory neuron"
sc@meta.data[sc@meta.data$seurat_clusters %in% c(0,1,2,7,12,4,11,20),]$before <- "Excitatory neuronal"
sc_choose <- subset(sc,subset = (before == "Excitatory neuronal"))
sc_ave<-AverageExpression(
	sc_choose,
	features=c("Gdf11"),
	group.by = "group3",
	slot = "data"
)	 ###0.0

sc_choose <- subset(sc,subset = (celltype == "Excitatory neuronal"))
sc_ave<-AverageExpression(
	sc_choose,
	features=c("Gdf11"),
	group.by = "group3",
	slot = "data"
)	###0.006992562



###检查细胞表达情况
sc_choose <- subset(sc,subset = (celltype == "Excitatory neuronal"))
sc_ex = as.matrix(GetAssayData(sc_choose,slot = "count")["Gdf11",])
colnames(sc_ex) <- "Gdf11"
sc_ex <- data.frame(sc_ex)
sc_ex$cellname <- rownames(sc_ex)
tmp <- sc_ex[sc_ex$Gdf11 >0,]

sc_choose$gdf11 <- ""
sc_choose@meta.data[rownames(tmp),]$gdf11 <- tmp$Gdf11
choose2 <- subset(sc_choose,subset = (gdf11 > 0))



####2.GDF11表达
##fos,cdkn1a及GDF11的表达量柱状图
setwd("/jdfssz1/ST_SUPERCELLS/P22Z10200N0639/tianyingming/project/hanlei_zju/scRNA/tmp")
library(SeuratDisk)
library(dplyr)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(ggpattern)
library(ggprism)
library(Rmisc)
require(reshape2)
lwd_pt <- .pt*72.27/96	##线宽
theme_set(theme_test(base_size = .pt*8, base_line_size = lwd_pt))
sc <- LoadH5Seurat("/jdfssz1/ST_SUPERCELLS/P22Z10200N0639/tianyingming/project/hanlei_zju/scRNA/00.data/mouse_brain_annotation_use.h5Seurat")
sc_choose <- subset(sc, subset = (celltype == "Oligodendrocyte"))
target_gene = c("Acvr1b","Tgfbr1", "Acvr1c","Acvr2a","Acvr2b")
sc_ex = data.frame(t(as.matrix(GetAssayData(sc_choose,slot = "data")[target_gene,])))
sc_ex2 = expm1(sc_ex)	##数据平滑处理——(log1p())和exmp1():用于计算 指数 减 1，即 exp()-1
sc_ex2$group = sc_choose@meta.data$group3
data = melt(sc_ex2)
data_sum = summarySE(data, measurevar="value", groupvars = c("group","variable"),na.rm = FALSE, conf.interval = 0.95, .drop = TRUE)
##summarySE：描述给出计数、平均值、标准差、平均值的标准误和置信区间(默认95%)。
##a data frame with count, mean, standard deviation, standard error of the mean, and confidence interval
exp = data_sum[,c("group","variable","value")]
##value :平均值
##reshape2包中的dcast函数和acast函数，两个函数都可以将长格式数据转换成宽格式数据。dcast与acast几乎没有区别，唯一的差别在于acast函数的输出结果没有行标签，dcast函数的输出结果有行标签。
exp2 = dcast(exp, variable ~ group,value.var = "value")
rownames(exp2)= exp2$variable
exp2 = exp2[,-1]
rate = 1/exp2[,1]
exp_bar = apply(exp2,2, function(x){x=x*rate})
exp_bar2 = melt(exp_bar)
data_sum$exp_adj=exp_bar2$value

genes = unique(data_sum$variable)
for (i in 1:length(genes)){
	ggplot(data = data_sum[data_sum$variable == genes[i],],mapping = aes(x = group, y =exp_adj))+ 
		geom_bar(stat = 'identity',width=0.85, position = 'dodge')+ theme_test()+
		geom_text(mapping = aes(label = round(exp_adj, 3)),position = position_nudge(y = 0.2))+
		labs(x = genes[i], y = "Relative expression",size=8)+
		geom_col_pattern(aes(pattern_type = group),pattern = c('none','stripe','none','none'),fill=c("white","white","pink","red"),pattern_fill=c("white","blue","pink","red"),colour=c("blue","blue","pink","red"))+
		scale_pattern_type_discrete(choices = c('none','stripe','none','none'))+scale_pattern_manual(values = c(`Ctrl_GFP+` = "stripe", `Ctrl_GFP-` = "none",`KO_GFP+` = "none", `KO_GFP-` = "none"))+
		geom_errorbar(aes(ymin=(exp_adj-se),ymax=(exp_adj+se)),width=0.2,position = position_dodge(width = 0.85))+
		theme(axis.text.x = element_text(size=8,angle = 45,hjust = 0.5, vjust = 0.5),axis.text.y = element_text(size=8),legend.position="none")+
  geom_signif(comparisons = list(c("Ctrl_GFP-", "KO_GFP+"),
                                 c("Ctrl_GFP+","KO_GFP+"),
                                 c("KO_GFP-","KO_GFP+")),
              map_signif_level=T,
              textsize=6,test=t.test,step_increase=0.0)
         ##如果把theme_test放在其后，theme将无效
	ggsave(paste("barplot_",genes[i],"_new2.pdf",sep=""),width = 4, height=4)
}
#差异显著性检验
##T检验
s <- sc_ex2 %>%
  gather(key = "celltype", value = "score", c(1:length(colnames(sc_ex2))-1)) %>%
  convert_as_factor(celltype)
set.seed(123)
a <- s %>% sample_n_by(group, celltype, size = 8000,replace=T)
stat.test <- a %>%
  group_by(celltype) %>%
  t_test(
    score ~ group, paired = F, var.equal= TRUE,
    alternative = 'two.sided'
  ) %>%
  select(-df, -statistic, -p)
write.csv(stat.test,"t_test.csv")
##Wilcox检验
df1 <- FindMarkers(object = sc_choose, ident.1 = "KO_GFP+",ident.2 = "KO_GFP-", min.pct = 0,logfc.threshold = 0)
df <- df1[c("Acvr1b","Tgfbr1", "Acvr1c","Acvr2a","Acvr2b"),]
write.csv(df,"diff_gene_KO_GFP+_KO_GFP-_5genes.csv")
write.csv(df1,"diff_gene_KO_GFP+_KO_GFP-.csv")
df1 <- FindMarkers(object = sc_choose, ident.1 = "KO_GFP+",ident.2 = "Ctrl_GFP-", min.pct = 0,logfc.threshold = 0)
df <- df1[c("Acvr1b","Tgfbr1", "Acvr1c","Acvr2a","Acvr2b"),]
write.csv(df,"diff_gene_KO_GFP+_Ctrl_GFP-_5genes.csv")
write.csv(df1,"diff_gene_KO_GFP+_Ctrl_GFP-.csv")
df1 <- FindMarkers(object = sc_choose, ident.1 = "KO_GFP+",ident.2 = "Ctrl_GFP+", min.pct = 0,logfc.threshold = 0)
df <- df1[c("Acvr1b","Tgfbr1", "Acvr1c","Acvr2a","Acvr2b"),]
write.csv(df,"diff_gene_KO_GFP+_Ctrl_GFP+_5genes.csv")
write.csv(df1,"diff_gene_KO_GFP+_Ctrl_GFP+.csv")

sem = data_sum[,c("group","variable","se")]
sem2 = dcast(sem, variable ~ group,value.var = "se")
rownames(sem2)= sem2$variable
sem2 = sem2[,-1]
sem3 = apply(sem2,2, function(x){x=x*rate})
sem4 = melt(sem3)
data_sum$se_adj=sem4$value

plots <- list()
genes = unique(data_sum$variable)
for (i in 1:length(genes)){
plots[[genes[i]]] = ggplot(data = data_sum[data_sum$variable == genes[i],],mapping = aes(x = group, y =exp_adj))+ 
		geom_bar(stat = 'identity',width=0.85, position = 'dodge')+ theme_test()+
		geom_text(mapping = aes(label = round(exp_adj, 3)),position = position_nudge(y = 0.2))+
		labs(x = genes[i], y = "Relative expression",size=8)+
		geom_col_pattern(aes(pattern_type = group),pattern=c('none','stripe','none','none'),fill=c("white","white","pink","red"),pattern_fill=c("white","blue","pink","red"),colour=c("blue","blue","pink","red"))+
		scale_pattern_type_discrete(choices = c('none','stripe','none','none'))+ scale_pattern_manual(values = c(`Ctrl_GFP+` = "stripe", `Ctrl_GFP-` = "none",`KO_GFP+` = "none", `KO_GFP-` = "none"))+
		geom_errorbar(aes(ymin=(exp_adj-se),ymax=(exp_adj+se)),width=0.2,position = position_dodge(width = 0.85))+
		theme(axis.text.x = element_text(size=8,angle = 45,hjust = 0.5, vjust = 0.5),axis.text.y = element_text(size=8),legend.position="none") ##如果把theme_test放在其后，theme将无效
}
p <- ggarrange(plotlist=plots, ncol=3,nrow = 1)
ggsave("barplot_3gene_split.pdf",plot=p,width = 9, height = 4)




sem = data_sum[,c("group","variable","value")]
sem2 = dcast(sem, variable ~ group,value.var = "value")
rownames(sem2)= sem2$variable
sem2 = sem2[,-1]
sem3 = apply(sem2,2, function(x){x=x*rate})
sem4 = melt(sem3)
data_sum$se_adj=sem4$value

plots <- list()
genes = unique(data_sum$variable)
for (i in 1:length(genes)){
plots[[genes[i]]] = ggplot(data = data_sum[data_sum$variable == genes[i],],mapping = aes(x = group, y =exp_adj))+ 
		geom_bar(stat = 'identity',width=0.85, position = 'dodge')+ theme_test()+
		geom_text(mapping = aes(label = round(exp_adj, 3)),position = position_nudge(y = 0.2))+
		labs(x = genes[i], y = "Relative expression",size=8)+
		geom_col_pattern(aes(pattern_type = group),pattern=c('none','stripe','none','none'),fill=c("white","white","pink","red"),pattern_fill=c("white","blue","pink","red"),colour=c("blue","blue","pink","red"))+
		scale_pattern_type_discrete(choices = c('none','stripe','none','none'))+ scale_pattern_manual(values = c(`Ctrl_GFP+` = "stripe", `Ctrl_GFP-` = "none",`KO_GFP+` = "none", `KO_GFP-` = "none"))+
		geom_errorbar(aes(ymin=(exp_adj-se),ymax=(exp_adj+se)),width=0.2,position = position_dodge(width = 0.85))+
		theme(axis.text.x = element_text(size=8,angle = 45,hjust = 0.5, vjust = 0.5),axis.text.y = element_text(size=8),legend.position="none") ##如果把theme_test放在其后，theme将无效
}
p <- ggarrange(plotlist=plots, ncol=3,nrow = 1)
ggsave("barplot_3gene_split.pdf",plot=p,width = 9, height = 4)






####3.aging基因表达热图
setwd("/jdfssz1/ST_SUPERCELLS/P20Z10200N0059/tianyingming/project/hanlei_zju/scRNA/02.analysis/4group")
library(SeuratDisk)
library(dplyr)
library(Seurat)
library(stringr)
library(ggplot2)
require(RColorBrewer)
require(circlize)
require(pheatmap)
#require(ComplexHeatmap)
sc <- LoadH5Seurat("/jdfssz1/ST_SUPERCELLS/P20Z10200N0059/tianyingming/project/hanlei_zju/scRNA/00.data/mouse_brain_annotation_use.h5Seurat")
aging_gene <- read.table("../aging_gene_use",sep="\t",header=T)
sasp_gene <- read.table("../sasp_gene_use",sep="\t",header=T)
ex_celltypes <- grep("excitatory",sc@meta.data$celltype_large,ignore.case = T,value = T)
ex_celltypes2 <- unique(ex_celltypes)
	for (i in 1:length(ex_celltypes2)){
		sc_choose <- subset(sc, subset = (celltype_large == ex_celltypes2[i]))
		sc_ave<-AverageExpression(
			sc_choose,
			group.by = "group3",
			slot = "data"
			)
		gene_ave <-as.matrix(sc_ave[[1]])
		gene_ave2 = t(scale(t(gene_ave),scale = F,center=T))
		gene_ave2[gene_ave2>=0.6]=0.6
		gene_ave2[gene_ave2<=-0.6]=-0.6
		heatmap_gene_tmp = data.frame(gene_ave2[aging_gene$Symbol,])	
		heatmap_gene_tmp2 = heatmap_gene_tmp %>% dplyr::arrange(heatmap_gene_tmp[,4])
		heatmap_gene = gene_ave2[rownames(heatmap_gene_tmp2),]

		#heatmap_gene = gene_ave2[aging_gene$Symbol,]
		
		heatmap_factor = gene_ave2[sasp_gene$Symbol,]



		pdf(paste(ex_celltypes2[i],"_heatmap_4group_ave_aging_gene.pdf",sep=""),width=4,height=6)
		p<-pheatmap(heatmap_gene,
			fontsize=6,fontsize_col=10, border_color = "NA",color = colorRampPalette(c("blue", "white","red"))(256),
			cluster_cols = F,cluster_rows = F,
			show_rownames=T,show_colnames=T,
			)
		print(p)
		dev.off()

		pdf(paste(ex_celltypes2[i],"_heatmap_4group_ave_sasp_gene.pdf",sep=""),width=4,height=5)
		p<-pheatmap(heatmap_factor,
			fontsize=6,fontsize_col=10, border_color = "NA",color = colorRampPalette(c("blue", "white","red"))(256),
			cluster_cols = F,cluster_rows = F, 
			show_rownames=T,show_colnames=T
			)
		print(p)
		dev.off()
	}










####4.统计不同细胞类型的细胞数
setwd("/jdfssz1/ST_SUPERCELLS/P20Z10200N0059/tianyingming/project/hanlei_zju/scRNA/tmp")
library(ggplot2)
library(ggrepel)
library(reshape2)
library(tidyr)
library(dplyr)
library(RColorBrewer)
sc <- LoadH5Seurat("/jdfssz1/ST_SUPERCELLS/P20Z10200N0059/tianyingming/project/hanlei_zju/scRNA/00.data/mouse_brain_annotation_version2.h5Seurat")
Case_control_fraction <- table(sc$celltype,sc$group3)
data <- melt(Case_control_fraction)
colnames(data)= c("celltype","condition","cell_num")

celltypes <- unique(data$celltype)
data$cell_sum <- 0
data$cell_percent <- 0
for (i in 1:length(celltypes)){
	data_cell <- data[data$celltype==celltypes[i],]
	data[data$celltype==celltypes[i],]$cell_sum <- sum(data_cell$cell_num)
	data[data$celltype==celltypes[i],]$cell_percent <- (data[data$celltype==celltypes[i],]$cell_num)/(data[data$celltype==celltypes[i],]$cell_sum)
}

cell_order = c("Excitatory neuronal","Cholinergic inhibitory neuron","Cnr1+ inhibitory neuron","Drd1+ inhibitory neuron","Drd2+ inhibitory neuron","Npy inhibitory neuron","PV+ inhibitory neuron","SST+ inhibitory neuron","VIP+ inhibitory neuron","Zic1+ inhibitory neuron","Astrocyte","Endothelial cell","Fibroblast-like cell","Microglia","Oligodendrocyte","Oligodendrocyte precursor")
data$celltype <- factor(x = data$celltype, levels = cell_order)

ggplot(data=data,mapping=aes(x=celltype,y=cell_percent,fill = condition))+   
	geom_bar(stat="identity",position = "stack")+theme_test() +   
	labs(x="Celltype",y="Cell number",size=12)+
	theme(axis.text=element_text(size=12))+
	#theme(axis.text=element_text(size=12),axis.text.x=element_text(angle=45))+
	geom_text(aes(label = cell_num), position = position_stack(vjust = 0.5))+
	coord_flip()+
	scale_fill_manual(values=c("brown","purple", "cadetblue2","yellow"),breaks = c("Ctrl_GFP-","Ctrl_GFP+","KO_GFP-","KO_GFP+"))
ggsave("Case_control_gfp_fraction.pdf",width = 12, height = 10)



#####5-10图
setwd("/jdfssz1/ST_SUPERCELLS/P20Z10200N0059/tianyingming/project/hanlei_zju/scRNA/03.celltype_deg/excitatory_all2")
#go_up ...gfp+ ko vs ctrl
markers_use <- read.csv("Excitatory neuronal_GFP+_KOvsCtrl_sig.csv",row.names=1)
marker_up <- markers_use[markers_use$down_up=="up",]
ids <- bitr(rownames(marker_up),fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = "org.Mm.eg.db")
GO_up <- enrichGO(gene= ids$ENTREZID, OrgDb = org.Mm.eg.db, ont='BP',pAdjustMethod = 'BH', qvalueCutoff = 0.2, pvalueCutoff = 0.05,keyType = 'ENTREZID')
go_top10 <- head(GO_up[order(GO_up$pvalue, decreasing = F), ], 10)
go_top10 <- go_top10 %>% dplyr::arrange(desc(go_top10$pvalue))
go_top10$Description <- factor(go_top10$Description,levels=go_top10$Description) 
p <- ggplot(data=go_top10,aes(x=Description,y=-log10(pvalue)))+ 
	geom_bar(stat="identity")+
	labs(x="",y="-log10(pvalue)",title="")+ 
	#scale_x_discrete(labels=function(x) str_wrap(x, width=40))+
	coord_flip()+
	theme(plot.title = element_text(hjust = 0.5))+
	theme_test()
ggsave("Excitatory neuronal_barplot_GO_up_GFP+_KOvsCtrl.pdf",plot = p,height=4,width=6)

marker_down <- markers_use[markers_use$down_up=="down",]
ids <- bitr(rownames(marker_down),fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = "org.Mm.eg.db")
GO_down <- enrichGO(gene= ids$ENTREZID, OrgDb = org.Mm.eg.db, ont='BP',pAdjustMethod = 'BH', qvalueCutoff = 0.2, pvalueCutoff = 0.05,keyType = 'ENTREZID')
go_top10 <- head(GO_down[order(GO_down$pvalue, decreasing = F), ], 10)
go_top10 <- go_top10 %>% dplyr::arrange(desc(go_top10$pvalue))
go_top10$Description <- factor(go_top10$Description,levels=go_top10$Description) 
p <- ggplot(data=go_top10,aes(x=Description,y=-log10(pvalue)))+ 
	geom_bar(stat="identity")+
	labs(x="",y="-log10(pvalue)",title="")+ 
	#scale_x_discrete(labels=function(x) str_wrap(x, width=40))+
	coord_flip()+
	theme(plot.title = element_text(hjust = 0.5))+
	theme_test()
ggsave("Excitatory neuronal_barplot_GO_down_GFP+_KOvsCtrl.pdf",plot = p,height=4,width=6)




markers_use <- read.csv("Excitatory neuronal_KO_GFP+vsGFP-_sig.csv",row.names=1)
marker_up <- markers_use[markers_use$down_up=="up",]
ids <- bitr(rownames(marker_up),fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = "org.Mm.eg.db")
GO_up <- enrichGO(gene= ids$ENTREZID, OrgDb = org.Mm.eg.db, ont='BP',pAdjustMethod = 'BH', qvalueCutoff = 0.2, pvalueCutoff = 0.05,keyType = 'ENTREZID')
go_top10 <- head(GO_up[order(GO_up$pvalue, decreasing = F), ], 10)
go_top10 <- go_top10 %>% dplyr::arrange(desc(go_top10$pvalue))
go_top10$Description <- factor(go_top10$Description,levels=go_top10$Description) 
p <- ggplot(data=go_top10,aes(x=Description,y=-log10(pvalue)))+ 
	geom_bar(stat="identity")+
	labs(x="",y="-log10(pvalue)",title="")+ 
	#scale_x_discrete(labels=function(x) str_wrap(x, width=40))+
	coord_flip()+
	theme(plot.title = element_text(hjust = 0.5))+
	theme_test()
ggsave("Excitatory neuronal_barplot_GO_up_KO_GFP+vsGFP-.pdf",plot = p,height=4,width=6)

marker_down <- markers_use[markers_use$down_up=="down",]
ids <- bitr(rownames(marker_down),fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = "org.Mm.eg.db")
GO_down <- enrichGO(gene= ids$ENTREZID, OrgDb = org.Mm.eg.db, ont='BP',pAdjustMethod = 'BH', qvalueCutoff = 0.2, pvalueCutoff = 0.05,keyType = 'ENTREZID')
go_top10 <- head(GO_down[order(GO_down$pvalue, decreasing = F), ], 10)
go_top10 <- go_top10 %>% dplyr::arrange(desc(go_top10$pvalue))
go_top10$Description <- factor(go_top10$Description,levels=go_top10$Description) 
p <- ggplot(data=go_top10,aes(x=Description,y=-log10(pvalue)))+ 
	geom_bar(stat="identity")+
	labs(x="",y="-log10(pvalue)",title="")+ 
	#scale_x_discrete(labels=function(x) str_wrap(x, width=40))+
	coord_flip()+
	theme(plot.title = element_text(hjust = 0.5))+
	theme_test()
ggsave("Excitatory neuronal_barplot_GO_down_KO_GFP+vsGFP-.pdf",plot = p,height=4,width=6)


##volcano
sc_choose.markers <- read.csv("Excitatory neuronal_GFP+_KOvsCtrl_all.csv",row.names=1)
top5 <-  sc_choose.markers %>% top_n(n = 10, wt = abs(avg_log2FC))
sc_choose.markers$condition = ifelse(sc_choose.markers$p_val_adj < 0.05 & sc_choose.markers$avg_log2FC < -0.25,"down",
	ifelse(sc_choose.markers$p_val_adj < 0.05 & sc_choose.markers$avg_log2FC > 0.25,"up","stable"))
ggplot(sc_choose.markers, aes(avg_log2FC, -log10(p_val_adj))) +
	geom_point(size=1, aes(color=condition)) +
	ylab("-log10(adjusted Pvalue)")+
	scale_color_manual(values=c("blue", "grey","red"))+
	geom_text_repel(data = top5, aes(avg_log2FC,-log10(p_val_adj),label = rownames(top5)),size=6/.pt,repel = T,min.segment.length = 0,segment.color="black",force=T)+
	geom_vline(xintercept = c(-0.25, 0.25), lty=2,col="black",lwd=0.4) +
	geom_hline(yintercept = -log10(0.05),lty=2,col="black",lwd=0.4) +
	theme_test()
ggsave("Excitatory neuronal_volcano_GFP+_KOvsCtrl.pdf",height=4,width=5)


sc_choose.markers <- read.csv("Excitatory neuronal_KO_GFP+vsGFP-_all.csv",row.names=1)
top5 <-  sc_choose.markers %>% top_n(n = 10, wt = abs(avg_log2FC))
sc_choose.markers$condition = ifelse(sc_choose.markers$p_val_adj < 0.05 & sc_choose.markers$avg_log2FC < -0.25,"down",
	ifelse(sc_choose.markers$p_val_adj < 0.05 & sc_choose.markers$avg_log2FC > 0.25,"up","stable"))
ggplot(sc_choose.markers, aes(avg_log2FC, -log10(p_val_adj))) +
	geom_point(size=1, aes(color=condition)) +
	ylab("-log10(adjusted Pvalue)")+
	scale_color_manual(values=c("blue", "grey","red"))+
	geom_text_repel(data = top5, aes(avg_log2FC,-log10(p_val_adj),label = rownames(top5)),size=6/.pt,repel = T,min.segment.length = 0,segment.color="black",force=T)+
	geom_vline(xintercept = c(-0.25, 0.25), lty=2,col="black",lwd=0.4) +
	geom_hline(yintercept = -log10(0.05),lty=2,col="black",lwd=0.4) +
	theme_test()
ggsave("Excitatory neuronal_volcano_KO_GFP+vsGFP-.pdf",height=4,width=5)







setwd("91914203.celltype_deg/excitatory_all2")
library(SeuratDisk)
library(dplyr)
library(Seurat)
library(stringr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
require(RColorBrewer)
require(circlize)
library(clusterProfiler)
library(org.Mm.eg.db)
sc <- LoadH5Seurat("/jdfssz1/ST_SUPERCELLS/P20Z10200N0059/tianyingming/project/hanlei_zju/scRNA/00.data/mouse_brain_annotation_use.h5Seurat") 
ex_celltypes <- grep("excitatory",sc@meta.data$celltype,ignore.case = T,value = T)
ex_celltypes2 <- unique(ex_celltypes)
for (i in 1:length(ex_celltypes2)){
	sc_choose <- subset(sc, subset = (celltype == ex_celltypes2[i]))		
	sc_choose_GFP_plus <- subset(sc_choose, subset = (group2 %in% c("GFP+")))
	sc_choose_GFP_plus_KO <- subset(sc_choose_GFP_plus, subset = (group1 %in% c("KO")))
	sc_choose_GFP_plus_Ctrl <- subset(sc_choose_GFP_plus, subset = (group1 %in% c("Ctrl")))
	sc_choose.markers <- FindMarkers(sc_choose_GFP_plus,slot = "data",ident.1=rownames(sc_choose_GFP_plus_KO@meta.data),ident.2=rownames(sc_choose_GFP_plus_Ctrl@meta.data),group.by="group1",logfc.threshold = -Inf) 
	sc_choose.markers$down_up <- "NA"
	sc_choose.markers[sc_choose.markers$avg_log2FC >0,]$down_up <- "up"
	sc_choose.markers[sc_choose.markers$avg_log2FC <=0,]$down_up <- "down"
	
	if(dim(sc_choose.markers)[1]>2){write.csv(sc_choose.markers,paste(ex_celltypes2[i],"_GFP+_KOvsCtrl_all.csv",sep=""),quote=F)
	markers_use <- sc_choose.markers %>% subset(p_val<0.05 & abs(sc_choose.markers$avg_log2FC) > 0.25)
	if(dim(markers_use)[1]>2){write.csv(markers_use,paste(ex_celltypes2[i],"_GFP+_KOvsCtrl_sig.csv",sep=""),quote=F)
	}else{;}
	}
	#top10 <- sc_choose.markers %>% top_n(n = 20, wt = abs(avg_log2FC)) 
	top10_up <- sc_choose.markers[sc_choose.markers$down_up=="up",] %>% top_n(n = 10, wt = avg_log2FC)
	top10_down <- sc_choose.markers[sc_choose.markers$down_up=="down",] %>% top_n(n = 10, wt = -avg_log2FC)
	top10 <- rbind(top10_up,top10_down)
	#热图
	DoHeatmap(sc_choose, features = rownames(top10), group.by = "group1",size=3)
	ggsave(paste(ex_celltypes2[i],"_heatmap_GFP+_KOvsCtrl_top10.pdf",sep=""),width = 10, height = 10)
	#go
	ids <- bitr(rownames(markers_use),fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = "org.Mm.eg.db")
	GO_BP <- enrichGO(gene= ids$ENTREZID, OrgDb = org.Mm.eg.db, ont='BP',pAdjustMethod = 'BH', qvalueCutoff = 0.2, pvalueCutoff = 0.05,keyType = 'ENTREZID')
	write.csv(GO_BP,paste(ex_celltypes2[i],"_GO_BP_GFP+_KOvsCtrl.csv",sep=""),quote = FALSE)
	#barplot(GO_all,showCategory=12,font.size=7,title="GO_all")
	#ggsave(paste(ex_celltypes2[i],"_barplot_GO_all.pdf",sep=""))
	go_top10 <- head(GO_BP[order(GO_BP$pvalue, decreasing = F), ], 10)
	go_top10 <- go_top10 %>% dplyr::arrange(desc(go_top10$pvalue))
	go_top10$Description <- factor(go_top10$Description,levels=go_top10$Description) 
	p <- ggplot(data=go_top10,aes(x=Description,y=-log10(pvalue)))+ 
		geom_bar(stat="identity")+
		labs(x="",y="-log10(pvalue)",title="")+ 
		#scale_x_discrete(labels=function(x) str_wrap(x, width=40))+
		coord_flip()+
		theme(plot.title = element_text(hjust = 0.5))+
		theme_test()
	ggsave(paste(ex_celltypes2[i],"_barplot_GO_BP_GFP+_KOvsCtrl.pdf",sep=""),plot = p,height=4,width=6)
	##volcano
	sc_choose.markers$condition = ifelse(sc_choose.markers$p_val_adj < 0.01 & sc_choose.markers$avg_log2FC < -0.25,"down",
		ifelse(sc_choose.markers$p_val_adj < 0.01 & sc_choose.markers$avg_log2FC > 0.25,"up","stable"))
	ggplot(sc_choose.markers, aes(avg_log2FC, -log10(p_val_adj))) +
		geom_point(alpha=0.6, size=0.5, aes(color=condition)) +
		ylab("-log10(Pvalue)")+
		scale_color_manual(values=c("darkblue", "grey","darkred"))+
		geom_text_repel(data = top10, aes(label = rownames(top10)),size=1.5,repel = T)+
		geom_vline(xintercept = c(-0.25, 0.25), lty=2,col="black",lwd=0.4) +
		geom_hline(yintercept = -log10(0.01),lty=2,col="black",lwd=0.4) +
		theme_test()
	ggsave(paste(ex_celltypes2[i],"_volcano_GFP+_KOvsCtrl.pdf",sep=""),height=4,width=5)
	
	
	sc_choose_GFP_minus <- subset(sc_choose, subset = (group2 %in% c("GFP-")))
	sc_choose_GFP_minus_KO <- subset(sc_choose_GFP_minus, subset = (group1 %in% c("KO")))
	sc_choose_GFP_minus_Ctrl <- subset(sc_choose_GFP_minus, subset = (group1 %in% c("Ctrl")))
	sc_choose.markers <- FindMarkers(sc_choose_GFP_minus,slot = "data",ident.1=rownames(sc_choose_GFP_minus_KO@meta.data),ident.2=rownames(sc_choose_GFP_minus_Ctrl@meta.data),group.by="group1",logfc.threshold = -Inf) 
	sc_choose.markers$down_up <- "NA"
	sc_choose.markers[sc_choose.markers$avg_log2FC >0,]$down_up <- "up"
	sc_choose.markers[sc_choose.markers$avg_log2FC <=0,]$down_up <- "down"
	if(dim(sc_choose.markers)[1]>2){write.csv(sc_choose.markers,paste(ex_celltypes2[i],"_GFP-_KOvsCtrl_all.csv",sep=""),quote=F)
	markers_use <- sc_choose.markers %>% subset(p_val<0.05 & abs(sc_choose.markers$avg_log2FC) > 0.25)
	if(dim(markers_use)[1]>2){write.csv(markers_use,paste(ex_celltypes2[i],"_GFP-_KOvsCtrl_sig.csv",sep=""),quote=F)
	}else{;}
	}
	#top10 <- sc_choose.markers %>% top_n(n = 20, wt = abs(avg_log2FC)) 
	top10_up <- sc_choose.markers[sc_choose.markers$down_up=="up",] %>% top_n(n = 10, wt = avg_log2FC)
	top10_down <- sc_choose.markers[sc_choose.markers$down_up=="down",] %>% top_n(n = 10, wt = -avg_log2FC)
	top10 <- rbind(top10_up,top10_down)
	#热图
	DoHeatmap(sc_choose, features = rownames(top10), group.by = "group1",size=3)
	ggsave(paste(ex_celltypes2[i],"_heatmap_GFP-_KOvsCtrl_top10.pdf",sep=""),width = 10, height = 10)
	#go
	ids <- bitr(rownames(markers_use),fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = "org.Mm.eg.db")
	GO_BP <- enrichGO(gene= ids$ENTREZID, OrgDb = org.Mm.eg.db, ont='BP',pAdjustMethod = 'BH', qvalueCutoff = 0.2, pvalueCutoff = 0.05,keyType = 'ENTREZID')
	write.csv(GO_BP,paste(ex_celltypes2[i],"_GO_BP_GFP-_KOvsCtrl.csv",sep=""),quote = FALSE)
	#barplot(GO_all,showCategory=12,font.size=7,title="GO_all")
	#ggsave(paste(ex_celltypes2[i],"_barplot_GO_all.pdf",sep=""))
	go_top10 <- head(GO_BP[order(GO_BP$pvalue, decreasing = F), ], 10)
	go_top10 <- go_top10 %>% dplyr::arrange(desc(go_top10$pvalue))
	go_top10$Description <- factor(go_top10$Description,levels=go_top10$Description) 
	p <- ggplot(data=go_top10,aes(x=Description,y=-log10(pvalue)))+ 
		geom_bar(stat="identity")+
		labs(x="",y="-log10(pvalue)",title="")+ 
		#scale_x_discrete(labels=function(x) str_wrap(x, width=40))+
		coord_flip()+
		theme(plot.title = element_text(hjust = 0.5))+
		theme_test()
	ggsave(paste(ex_celltypes2[i],"_barplot_GO_BP_GFP-_KOvsCtrl.pdf",sep=""),plot = p,height=4,width=6)
	##volcano
	sc_choose.markers$condition = ifelse(sc_choose.markers$p_val_adj < 0.01 & sc_choose.markers$avg_log2FC < -0.25,"down",
		ifelse(sc_choose.markers$p_val_adj < 0.01 & sc_choose.markers$avg_log2FC > 0.25,"up","stable"))
	ggplot(sc_choose.markers, aes(avg_log2FC, -log10(p_val_adj))) +
		geom_point(alpha=0.6, size=0.5, aes(color=condition)) +
		ylab("-log10(Pvalue)")+
		scale_color_manual(values=c("blue", "grey","red"))+
		geom_text_repel(data = top10, aes(label = rownames(top10)),size=1.5,repel = T)+
		geom_vline(xintercept = c(-0.25, 0.25), lty=2,col="black",lwd=0.4) +
		geom_hline(yintercept = -log10(0.01),lty=2,col="black",lwd=0.4) +
		theme_test()
	ggsave(paste(ex_celltypes2[i],"_volcano_GFP-_KOvsCtrl.pdf",sep=""),height=4,width=5)
	
	
	
	
	sc_choose_Ctrl <- subset(sc_choose, subset = (group1 %in% c("Ctrl")))
	sc_choose_Ctrl_GFP_plus <- subset(sc_choose_Ctrl, subset = (group2 %in% c("GFP+")))
	sc_choose_Ctrl_GFP_minus <- subset(sc_choose_Ctrl, subset = (group2 %in% c("GFP-")))
	sc_choose.markers <- FindMarkers(sc_choose_Ctrl,slot = "data",ident.1=rownames(sc_choose_Ctrl_GFP_plus@meta.data),ident.2=rownames(sc_choose_Ctrl_GFP_minus@meta.data),group.by="group2",logfc.threshold = -Inf) 
	sc_choose.markers$down_up <- "NA"
	sc_choose.markers[sc_choose.markers$avg_log2FC >0,]$down_up <- "up"
	sc_choose.markers[sc_choose.markers$avg_log2FC <=0,]$down_up <- "down"
	if(dim(sc_choose.markers)[1]>2){write.csv(sc_choose.markers,paste(ex_celltypes2[i],"_Ctrl_GFP+vsGFP-_all.csv",sep=""),quote=F)
	markers_use <- sc_choose.markers %>% subset(p_val<0.05 & abs(sc_choose.markers$avg_log2FC) > 0.25)
	if(dim(markers_use)[1]>2){write.csv(markers_use,paste(ex_celltypes2[i],"_Ctrl_GFP+vsGFP-_sig.csv",sep=""),quote=F)
	}else{;}
	}
	#top10 <- sc_choose.markers %>% top_n(n = 20, wt = abs(avg_log2FC)) 
	top10_up <- sc_choose.markers[sc_choose.markers$down_up=="up",] %>% top_n(n = 10, wt = avg_log2FC)
	top10_down <- sc_choose.markers[sc_choose.markers$down_up=="down",] %>% top_n(n = 10, wt = -avg_log2FC)
	top10 <- rbind(top10_up,top10_down)
	#热图
	DoHeatmap(sc_choose, features = rownames(top10), group.by = "group1",size=3)
	ggsave(paste(ex_celltypes2[i],"_heatmap_Ctrl_GFP+vsGFP-_top10.pdf",sep=""),width = 10, height = 10)
	#go
	ids <- bitr(rownames(markers_use),fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = "org.Mm.eg.db")
	GO_BP <- enrichGO(gene= ids$ENTREZID, OrgDb = org.Mm.eg.db, ont='BP',pAdjustMethod = 'BH', qvalueCutoff = 0.2, pvalueCutoff = 0.05,keyType = 'ENTREZID')
	write.csv(GO_BP,paste(ex_celltypes2[i],"_GO_BP_Ctrl_GFP+vsGFP-.csv",sep=""),quote = FALSE)
	#barplot(GO_all,showCategory=12,font.size=7,title="GO_all")
	#ggsave(paste(ex_celltypes2[i],"_barplot_GO_all.pdf",sep=""))
	go_top10 <- head(GO_BP[order(GO_BP$pvalue, decreasing = F), ], 10)
	go_top10 <- go_top10 %>% dplyr::arrange(desc(go_top10$pvalue))
	go_top10$Description <- factor(go_top10$Description,levels=go_top10$Description) 
	p <- ggplot(data=go_top10,aes(x=Description,y=-log10(pvalue)))+ 
		geom_bar(stat="identity")+
		labs(x="",y="-log10(pvalue)",title="")+ 
		#scale_x_discrete(labels=function(x) str_wrap(x, width=40))+
		coord_flip()+
		theme(plot.title = element_text(hjust = 0.5))+
		theme_test()
	ggsave(paste(ex_celltypes2[i],"_barplot_GO_BP_Ctrl_GFP+vsGFP-.pdf",sep=""),plot = p,height=4,width=6)
	##volcano
	sc_choose.markers$condition = ifelse(sc_choose.markers$p_val_adj < 0.01 & sc_choose.markers$avg_log2FC < -0.25,"down",
		ifelse(sc_choose.markers$p_val_adj < 0.01 & sc_choose.markers$avg_log2FC > 0.25,"up","stable"))
	ggplot(sc_choose.markers, aes(avg_log2FC, -log10(p_val_adj))) +
		geom_point(alpha=0.6, size=0.5, aes(color=condition)) +
		ylab("-log10(Pvalue)")+
		scale_color_manual(values=c("blue", "grey","red"))+
		geom_text_repel(data = top10, aes(label = rownames(top10)),size=1.5,repel = T)+
		geom_vline(xintercept = c(-0.25, 0.25), lty=2,col="black",lwd=0.4) +
		geom_hline(yintercept = -log10(0.01),lty=2,col="black",lwd=0.4) +
		theme_test()
	ggsave(paste(ex_celltypes2[i],"_volcano_Ctrl_GFP+vsGFP-.pdf",sep=""),height=4,width=5)


	
	sc_choose_KO <- subset(sc_choose, subset = (group1 %in% c("KO")))
	sc_choose_KO_GFP_plus <- subset(sc_choose_KO, subset = (group2 %in% c("GFP+")))
	sc_choose_KO_GFP_minus <- subset(sc_choose_KO, subset = (group2 %in% c("GFP-")))
	sc_choose.markers <- FindMarkers(sc_choose_KO,slot = "data",ident.1=rownames(sc_choose_KO_GFP_plus@meta.data),ident.2=rownames(sc_choose_KO_GFP_minus@meta.data),group.by="group2",logfc.threshold = -Inf) 
	sc_choose.markers$down_up <- "NA"
	sc_choose.markers[sc_choose.markers$avg_log2FC >0,]$down_up <- "up"
	sc_choose.markers[sc_choose.markers$avg_log2FC <=0,]$down_up <- "down"
	if(dim(sc_choose.markers)[1]>2){write.csv(sc_choose.markers,paste(ex_celltypes2[i],"_KO_GFP+vsGFP-_all.csv",sep=""),quote=F)
	markers_use <- sc_choose.markers %>% subset(p_val<0.05 & abs(sc_choose.markers$avg_log2FC) > 0.25)
	if(dim(markers_use)[1]>2){write.csv(markers_use,paste(ex_celltypes2[i],"_KO_GFP+vsGFP-_sig.csv",sep=""),quote=F)
	}else{;}
	}
	#top10 <- sc_choose.markers %>% top_n(n = 20, wt = abs(avg_log2FC)) 
	top10_up <- sc_choose.markers[sc_choose.markers$down_up=="up",] %>% top_n(n = 10, wt = avg_log2FC)
	top10_down <- sc_choose.markers[sc_choose.markers$down_up=="down",] %>% top_n(n = 10, wt = -avg_log2FC)
	top10 <- rbind(top10_up,top10_down)
	#热图
	DoHeatmap(sc_choose, features = rownames(top10), group.by = "group1",size=3)
	ggsave(paste(ex_celltypes2[i],"_heatmap_KO_GFP+vsGFP-_top10.pdf",sep=""),width = 10, height = 10)
	#go
	ids <- bitr(rownames(markers_use),fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = "org.Mm.eg.db")
	GO_BP <- enrichGO(gene= ids$ENTREZID, OrgDb = org.Mm.eg.db, ont='BP',pAdjustMethod = 'BH', qvalueCutoff = 0.2, pvalueCutoff = 0.05,keyType = 'ENTREZID')
	write.csv(GO_BP,paste(ex_celltypes2[i],"_GO_BP_Ctrl_GFP+vsGFP-.csv",sep=""),quote = FALSE)
	#barplot(GO_all,showCategory=12,font.size=7,title="GO_all")
	#ggsave(paste(ex_celltypes2[i],"_barplot_GO_all.pdf",sep=""))
	go_top10 <- head(GO_BP[order(GO_BP$pvalue, decreasing = F), ], 10)
	go_top10 <- go_top10 %>% dplyr::arrange(desc(go_top10$pvalue))
	go_top10$Description <- factor(go_top10$Description,levels=go_top10$Description) 
	p <- ggplot(data=go_top10,aes(x=Description,y=-log10(pvalue)))+ 
		geom_bar(stat="identity")+
		labs(x="",y="-log10(pvalue)",title="")+ 
		#scale_x_discrete(labels=function(x) str_wrap(x, width=40))+
		coord_flip()+
		theme(plot.title = element_text(hjust = 0.5))+
		theme_test()
	ggsave(paste(ex_celltypes2[i],"_barplot_GO_BP_KO_GFP+vsGFP-.pdf",sep=""),plot = p,height=4,width=6)
	##volcano
	sc_choose.markers$condition = ifelse(sc_choose.markers$p_val_adj < 0.01 & sc_choose.markers$avg_log2FC < -0.25,"down",
		ifelse(sc_choose.markers$p_val_adj < 0.01 & sc_choose.markers$avg_log2FC > 0.25,"up","stable"))
	ggplot(sc_choose.markers, aes(avg_log2FC, -log10(p_val_adj))) +
		geom_point(alpha=0.6, size=0.5, aes(color=condition)) +
		ylab("-log10(Pvalue)")+
		scale_color_manual(values=c("blue", "grey","red"))+
		geom_text_repel(data = top10, aes(label = rownames(top10)),size=1.5,repel = T)+
		geom_vline(xintercept = c(-0.25, 0.25), lty=2,col="black",lwd=0.4) +
		geom_hline(yintercept = -log10(0.01),lty=2,col="black",lwd=0.4) +
		theme_test()
	ggsave(paste(ex_celltypes2[i],"_volcano_KO_GFP+vsGFP-.pdf",sep=""),height=4,width=5)
}




####图11
setwd("/jdfssz1/ST_SUPERCELLS/P20Z10200N0059/tianyingming/project/hanlei_zju/scRNA/tmp")
library(SeuratDisk)
library(dplyr)
library(Seurat)
library(stringr)
library(ggplot2)
sc <- LoadH5Seurat("/jdfssz1/ST_SUPERCELLS/P20Z10200N0059/tianyingming/project/hanlei_zju/scRNA/00.data/mouse_brain_annotation_version2.h5Seurat") 
###Dotplot...KO和ctrl得单独提出来算，不然不对
sc@meta.data$group4 =""
sc@meta.data$group4 <- paste(sc@meta.data$celltype_large,sc@meta.data$group3,sep="_")
sc@meta.data$group4 <- gsub("Inhibitory neuron_KO_GFP\\+","Inhibitory neuron_KO",sc@meta.data$group4,perl=T)
sc@meta.data$group4 <- gsub("Inhibitory neuron_KO_GFP\\-","Inhibitory neuron_KO",sc@meta.data$group4,perl=T)
sc@meta.data$group4 <- gsub("Inhibitory neuron_Ctrl_GFP\\+","Inhibitory neuron_Ctrl",sc@meta.data$group4,perl=T)
sc@meta.data$group4 <- gsub("Inhibitory neuron_Ctrl_GFP\\-","Inhibitory neuron_Ctrl",sc@meta.data$group4,perl=T)

sc@meta.data$group4 <- gsub("Astrocyte_KO_GFP\\+","Astrocyte_KO",sc@meta.data$group4,perl=T)
sc@meta.data$group4 <- gsub("Astrocyte_KO_GFP\\-","Astrocyte_KO",sc@meta.data$group4,perl=T)
sc@meta.data$group4 <- gsub("Astrocyte_Ctrl_GFP\\+","Astrocyte_Ctrl",sc@meta.data$group4,perl=T)
sc@meta.data$group4 <- gsub("Astrocyte_Ctrl_GFP\\-","Astrocyte_Ctrl",sc@meta.data$group4,perl=T)

sc@meta.data$group4 <- gsub("Endothelial cell_KO_GFP\\+","Endothelial cell_KO",sc@meta.data$group4,perl=T)
sc@meta.data$group4 <- gsub("Endothelial cell_KO_GFP\\-","Endothelial cell_KO",sc@meta.data$group4,perl=T)
sc@meta.data$group4 <- gsub("Endothelial cell_Ctrl_GFP\\+","Endothelial cell_Ctrl",sc@meta.data$group4,perl=T)
sc@meta.data$group4 <- gsub("Endothelial cell_Ctrl_GFP\\-","Endothelial cell_Ctrl",sc@meta.data$group4,perl=T)

sc@meta.data$group4 <- gsub("Fibroblast-like cell_KO_GFP\\+","Fibroblast-like cell_KO",sc@meta.data$group4,perl=T)
sc@meta.data$group4 <- gsub("Fibroblast-like cell_KO_GFP\\-","Fibroblast-like cell_KO",sc@meta.data$group4,perl=T)
sc@meta.data$group4 <- gsub("Fibroblast-like cell_Ctrl_GFP\\+","Fibroblast-like cell_Ctrl",sc@meta.data$group4,perl=T)
sc@meta.data$group4 <- gsub("Fibroblast-like cell_Ctrl_GFP\\-","Fibroblast-like cell_Ctrl",sc@meta.data$group4,perl=T)

sc@meta.data$group4 <- gsub("Microglia_KO_GFP\\+","Microglia_KO",sc@meta.data$group4,perl=T)
sc@meta.data$group4 <- gsub("Microglia_KO_GFP\\-","Microglia_KO",sc@meta.data$group4,perl=T)
sc@meta.data$group4 <- gsub("Microglia_Ctrl_GFP\\+","Microglia_Ctrl",sc@meta.data$group4,perl=T)
sc@meta.data$group4 <- gsub("Microglia_Ctrl_GFP\\-","Microglia_Ctrl",sc@meta.data$group4,perl=T)

sc@meta.data$group4 <- gsub("Oligodendrocyte_KO_GFP\\+","Oligodendrocyte_KO",sc@meta.data$group4,perl=T)
sc@meta.data$group4 <- gsub("Oligodendrocyte_KO_GFP\\-","Oligodendrocyte_KO",sc@meta.data$group4,perl=T)
sc@meta.data$group4 <- gsub("Oligodendrocyte_Ctrl_GFP\\+","Oligodendrocyte_Ctrl",sc@meta.data$group4,perl=T)
sc@meta.data$group4 <- gsub("Oligodendrocyte_Ctrl_GFP\\-","Oligodendrocyte_Ctrl",sc@meta.data$group4,perl=T)

sc@meta.data$group4 <- gsub("Oligodendrocyte precursor_KO_GFP\\+","Oligodendrocyte precursor_KO",sc@meta.data$group4,perl=T)
sc@meta.data$group4 <- gsub("Oligodendrocyte precursor_KO_GFP\\-","Oligodendrocyte precursor_KO",sc@meta.data$group4,perl=T)
sc@meta.data$group4 <- gsub("Oligodendrocyte precursor_Ctrl_GFP\\+","Oligodendrocyte precursor_Ctrl",sc@meta.data$group4,perl=T)
sc@meta.data$group4 <- gsub("Oligodendrocyte precursor_Ctrl_GFP\\-","Oligodendrocyte precursor_Ctrl",sc@meta.data$group4,perl=T)

DotPlot(sc, features = "Cdkn1a",group.by="group4") + RotatedAxis()
ggsave("DotPlot_Cdkn1a.pdf")

setwd("/jdfssz1/ST_SUPERCELLS/P20Z10200N0059/tianyingming/project/hanlei_zju/scRNA/00.data")
library(SeuratDisk)
library(dplyr)
library(Seurat)
library(stringr)
library(ggplot2)
sc <- LoadH5Seurat("/jdfssz1/ST_SUPERCELLS/P20Z10200N0059/tianyingming/project/hanlei_zju/scRNA/00.data/mouse_brain_annotation_version2.h5Seurat") 
FeaturePlot(sc, features = "Cdkn1a",order=T)
ggsave("Featureplot_Cdkn1a.pdf")







####图13-14
#step1.得sc_matrix.csv
setwd("/jdfssz1/ST_SUPERCELLS/P20Z10200N0059/tianyingming/project/hanlei_zju/scRNA/06.pyscenic")
sc <- LoadH5Seurat("/jdfssz1/ST_SUPERCELLS/P20Z10200N0059/tianyingming/project/hanlei_zju/scRNA/00.data/mouse_brain_annotation_version1.h5Seurat") 
exprMat <- sc@assays$RNA@counts
exprMat <- as.matrix(exprMat)
exprMat_t <- t(exprMat)
write.csv(exprMat_t,"sc_matrix.csv")

###step2.
#!/bin/bash -l
echo "Start: `date`"
. /etc/profile.d/modules.sh

module load anaconda3/2021.05
module load R/4.0.5
python /jdfssz1/ST_SUPERCELLS/P20Z10200N0059/tianyingming/project/hanlei_zju/scRNA/06.pyscenic/scenic.wsp.py  /jdfssz1/ST_SUPERCELLS/P20Z10200N0059/tianyingming/project/hanlei_zju/scRNA/06.pyscenic/sc_matrix.csv  mouse_brain
#生成了AUC.txt,regulon.p,motifs.txt,modules.p,adjacencies.txt

###step3.将auc活性矩阵添加到meta.data中
setwd("/jdfssz1/ST_SUPERCELLS/P20Z10200N0059/tianyingming/project/hanlei_zju/scRNA/06.pyscenic")
library(SeuratDisk)
library(dplyr)
library(Seurat)
library(stringr)
library(ggplot2)
sc <- LoadH5Seurat("/jdfssz1/ST_SUPERCELLS/P20Z10200N0059/tianyingming/project/hanlei_zju/scRNA/00.data/mouse_brain_annotation_version2.h5Seurat") 
auc = read.table("mouse_brain_AUC.txt",sep="\t",header=T,check.names=F)
rownames(auc) = auc$Cell
auc.t = as.data.frame(t(auc))
auc.t=auc.t[-1,]
auc.t.f = auc.t[,colnames(sc)]

sc[["AUC"]] <- CreateAssayObject(counts = auc.t.f)
DefaultAssay(sc) <- "AUC"
SaveH5Seurat(sc, filename = "mouse_brain_auc.h5Seurat",overwrite = T)

Idents(sc)= "celltype" 
auc.marker = FindAllMarkers(sc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,test.use="wilcox",slot="counts")
write.table(auc.marker,"auc_markers.xls",sep="\t")

####step4.两两之间做差异分析
setwd("/jdfssz1/ST_SUPERCELLS/P20Z10200N0059/tianyingming/project/hanlei_zju/scRNA/06.pyscenic")
library(SeuratDisk)
library(dplyr)
library(Seurat)
library(stringr)
library(ggplot2)
sc <- LoadH5Seurat("mouse_brain_auc.h5Seurat") 
celltypes <- unique(sc@meta.data$celltype)
for (i in 1:length(celltypes)){
	sc_choose <- subset(sc, subset = (celltype == celltypes[i]))
	sc_choose_KO <- subset(sc_choose, subset = (group1 %in% c("KO")))
	sc_choose_Ctrl <- subset(sc_choose, subset = (group1 %in% c("Ctrl")))
	sc_choose.markers <- FindMarkers(sc_choose,slot = "count",ident.1=rownames(sc_choose_KO@meta.data),ident.2=rownames(sc_choose_Ctrl@meta.data),group.by="group1",logfc.threshold = -Inf,test.use="wilcox")
	write.csv(sc_choose.markers,paste(celltypes[i],"_KOvsCtrl_all.csv",sep=""),quote=F)
	markers_use <- sc_choose.markers %>% subset(p_val<0.05 & abs(sc_choose.markers$avg_log2FC) > 0.25)
	write.csv(markers_use,paste(celltypes[i],"_KOvsCtrl_sig.csv",sep=""),quote=F)
}

####step5.auc_plot.R

####step6.画热图,用step3中得到的marker画热图效果不好，用rss top5
rss.out=RunRSS(sc,
  assay = "AUC",
  slot = "counts",
  group.by = "celltype"
)
##画图辅助调顺序
rssranking = RSSRanking(
   rss.out,
   group.by="celltype",
   ggrepel_force = 1,
   ggrepel_point_padding = 0.2,
   top_genes = 6,
   plot_extended = FALSE
)
ggsave(plot=rssranking,"rssranking_top5.pdf",height=20,limitsize=F,width=10)

rss.out.top = rss.out %>% group_by(celltype)%>% top_n(n=5,wt=RSS)
write.table(rss.out.top,"rss_top5.xls",quote=F,sep="\t",row.names=F)

rss_ave<-AverageExpression(
		sc,
		assay = "AUC",
		group.by = "celltype",
		slot = "count"
		)
rss_ave1 <-as.matrix(rss_ave[[1]])


tf_gene <- read.table("../tf_sort.txt",sep="\t",header=T)
rss_heatmap = rss_ave1[tf_gene$TF,]
write.csv(rss_heatmap,"rss_ave_heatmap.csv",quote=F,sep="\t",row.names=T)
require(pheatmap)
rss_heatmap <- read.csv("rss_ave_heatmap.csv",row.names=1,check.names=F)
pdf("rsstop5_ave_heatmap_sort2.pdf",width=7,height=10)
pheatmap(rss_heatmap,
	fontsize=8,fontsize_col=10, border_color = "NA",
	cluster_cols = F,cluster_rows = F,treeheight_row=0,legend=T,
	show_rownames=T,show_colnames=T
	)
dev.off()



#####图15
setwd("/jdfssz1/ST_SUPERCELLS/P20Z10200N0059/tianyingming/project/hanlei_zju/scRNA/tmp")
library(Seurat)
library(SeuratDisk) 
library(ggplot2)
require(pheatmap)
cell_order = c("Excitatory neuronal","Cholinergic inhibitory neuron","Cnr1+ inhibitory neuron","Drd1+ inhibitory neuron","Drd2+ inhibitory neuron","Npy inhibitory neuron","PV+ inhibitory neuron","SST+ inhibitory neuron","VIP+ inhibitory neuron","Zic1+ inhibitory neuron","Astrocyte","Endothelial cell","Fibroblast-like cell","Microglia","Oligodendrocyte","Oligodendrocyte precursor")
sc <- LoadH5Seurat("/jdfssz1/ST_SUPERCELLS/P20Z10200N0059/tianyingming/project/hanlei_zju/scRNA/00.data/mouse_brain_annotation_use.h5Seurat") 
colors <- c(`Excitatory neuronal`="#FEAF16",`Cholinergic inhibitory neuron`="#B00068",`Cnr1+ inhibitory neuron`="#90AD1C",`Drd1+ inhibitory neuron`="#2ED9FF",`Drd2+ inhibitory neuron`="#DEA0FD",`Npy inhibitory neuron`="#AA0DFE",`PV+ inhibitory neuron`="#F8A19F",`SST+ inhibitory neuron`="#325A9B",`VIP+ inhibitory neuron`="#1C8356",`Zic1+ inhibitory neuron`="#C4451C",Astrocyte="#5A5156",`Endothelial cell`="#E4E1E3",`Fibroblast-like cell`="#F6222E",Microglia="#FE00FA",Oligodendrocyte="#16FF32",`Oligodendrocyte precursor`="#3283FE")
gene <- read.table("cell_marker_order",header = T,sep="\t",check.names = F)
sc$celltype <- factor(x = sc$celltype, levels = cell_order)
VlnPlot(sc, features = gene$gene, slot = "data",group.by="celltype",stack=T,fill.by = "ident",flip = T,cols = colors)
	labs(x="",y="")
ggsave("VlnPlot_celltype_test.pdf",width=12,height=14)


setwd("/jdfssz1/ST_SUPERCELLS/P20Z10200N0059/tianyingming/project/hanlei_zju/scRNA/tmp")
focus_gene <- read.table("cell_marker_order",sep="\t",header=T)
sc <- LoadH5Seurat("/jdfssz1/ST_SUPERCELLS/P20Z10200N0059/tianyingming/project/hanlei_zju/scRNA/00.data/mouse_brain_annotation_version2.h5Seurat") 
sc_ave<-AverageExpression(
	sc,
	group.by = "celltype",
	slot = "data"
)
gene_ave <-as.matrix(sc_ave[[1]])
gene_ave2 = t(scale(t(gene_ave),scale = T,center=T))
#gene_ave2=scale(gene_ave2, center = FALSE, scale = apply(gene_ave2, 2, sd, na.rm = TRUE))
#gene_ave2[gene_ave2>=0.6]=0.6
#gene_ave2[gene_ave2<=-0.6]=-0.6
#heatmap_gene_tmp = data.frame(gene_ave2[aging_gene$Symbol,])	
#heatmap_gene_tmp2 = heatmap_gene_tmp %>% dplyr::arrange(heatmap_gene_tmp[,4])
#heatmap_gene = gene_ave2[rownames(heatmap_gene_tmp2),]
#heatmap_gene = gene_ave2[aging_gene$Symbol,]

heatmap_gene = gene_ave2[focus_gene$gene,]
write.csv(heatmap_gene,"heatmap_order.csv")
heatmap_gene <- read.csv("heatmap_order.csv",row.names=1,check.names=F)	
cell_order = c("Excitatory neuronal","Cholinergic inhibitory neuron","Cnr1+ inhibitory neuron","Drd1+ inhibitory neuron","Drd2+ inhibitory neuron","Npy inhibitory neuron","PV+ inhibitory neuron","SST+ inhibitory neuron","VIP+ inhibitory neuron","Zic1+ inhibitory neuron","Astrocyte","Endothelial cell","Fibroblast-like cell","Microglia","Oligodendrocyte","Oligodendrocyte precursor")
pdf(paste("vlnplot","_heatmap2.pdf",sep=""),width=4,height=6)
p<-pheatmap(heatmap_gene,
	fontsize=6,fontsize_col=10, border_color = "NA",color = colorRampPalette(c("blue", "white","red"))(256),
	cluster_cols = F,cluster_rows = F,
	show_rownames=T,show_colnames=T,
	)
print(p)
dev.off()



####传数据
F:\scRNA\rawdata\cDNA\cDNA-3174-2-210818>certutil -hashfile DP8400020222TR_L01_71_1.fq.gz MD5

dir /b
