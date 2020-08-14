library(Seurat)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(GSVA)
library(velocyto.R)
library(pagoda2)
library(monocle)
library(data.table)
library(infercnv)
library(dplyr)
library(factoextra)
library(ROGUE)
### Global gene expression profiles and unsupervised clustering of B16 melanoma cells 

setwd("/home/kysbbubbu/B16")
merge.data <- Read10X(data.dir = "/data/MRC1_data4/kysbbubbu/B16/AGG_B16/outs/filtered_feature_bc_matrix")
pbmc <- CreateSeuratObject(merge.data, project = "merge", names.field = 2, names.delim = "-", min.cells = 5, min.features = 200)
current.cluster.ids <- c(1,2)
new.cluster.ids <- c("B16F0","B16F10")
pbmc@meta.data$orig.ident  <- plyr::mapvalues(x = pbmc@meta.data$orig.ident , from = current.cluster.ids, to = new.cluster.ids)

#Quality control 
mito.features <- grep(pattern = "^mt-", x = rownames(x = pbmc), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = pbmc, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = pbmc, slot = 'counts'))
pbmc[['percent.mito']] <- percent.mito
VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
#We decided filtering criteria based on below scatter plots
FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mito")
FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

#Figuer S1A 
ggplot(pbmc@meta.data, aes(x=nFeature_RNA, fill=orig.ident, color=orig.ident)) +geom_density(alpha=0.3)+  theme_classic()+theme(text = element_text(size=12))+theme(axis.title.x=element_text(color="black", size=12, face="bold"),axis.title.y=element_text(color="black", size=12, face="bold"),axis.text.x=element_text(color="black", size=12, face="bold"),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=12, face="bold"), strip.text = element_text(size = 12, face="bold"))+ylab("Denstiy")+xlab("Number of genes per cell")+scale_color_manual(values=c("#0000FF","#FF0000"))+scale_fill_manual(values=c("#0000FF","#FF0000"))+ geom_vline(xintercept=c(3000, 7000),linetype=3)+NoLegend()
ggplot(pbmc@meta.data, aes(x=percent.mito, fill=orig.ident, color=orig.ident)) +geom_density(alpha=0.3)+  theme_classic()+theme(text = element_text(size=12))+theme(axis.title.x=element_text(color="black", size=12, face="bold"),axis.title.y=element_text(color="black", size=12, face="bold"),axis.text.x=element_text(color="black", size=12, face="bold"),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=12, face="bold"), strip.text = element_text(size = 12, face="bold"))+ylab("Denstiy")+xlab("Percentage of mitochondrial gene expression")+scale_color_manual(values=c("#0000FF","#FF0000"))+scale_fill_manual(values=c("#0000FF","#FF0000"))+ geom_vline(xintercept=c(0.2),linetype=3)+NoLegend()
pbmc <- subset(x = pbmc, subset = nFeature_RNA > 3000 & nFeature_RNA < 7000 & percent.mito < 0.2)

#Feature selection, normalize
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 1e4)
pbmc = FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(object = pbmc, features = rownames(x = pbmc), vars.to.regress = c("nCount_RNA", "percent.mito"))

#Dimensional reduction and UMAP 
pbmc <- RunPCA(object = pbmc, features = VariableFeatures(object = pbmc), verbose = FALSE)
print(x = pbmc[['pca']], dims = 1:5, nfeatures = 10, projected = FALSE)
pbmc <- ProjectDim(object = pbmc)
pbmc <- JackStraw(object = pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(object = pbmc, dims = 1:20)
ElbowPlot(object = pbmc)
pbmc <- FindNeighbors(object = pbmc, dims = 1:12)
pbmc <- RunUMAP(pbmc, reduction.use = "pca", dims = 1:12)
pbmc <- FindClusters(object = pbmc, resolution = 0.05)
new.cluster.ids <- c("F0", "F10")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)

#FigureS1B
#Technical check 
FeaturePlot(object = pbmc, features = c("nFeature_RNA")) + ggtitle("Number of genes") + theme(plot.title = element_text(family = "arial", face = "bold", hjust = 0.5, size = 24))+ scale_colour_gradientn(colours=brewer.pal(9,"OrRd")[1:9])+ theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))

FeaturePlot(object = pbmc, features = c("nCount_RNA")) +ggtitle("Number of UMIs") + theme(plot.title = element_text(family = "arial", face = "bold", hjust = 0.5, size = 24))+ scale_colour_gradientn(colours=brewer.pal(9,"OrRd")[1:9]) +theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))

FeaturePlot(object = pbmc, features = c("percent.mito"))+ ggtitle("Mitochondria gene %") + theme(plot.title = element_text(family = "arial", face = "bold", hjust = 0.5, size = 24))+ scale_colour_gradientn(colours=brewer.pal(9,"OrRd")[1:9])+theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))

#Figure S1C
#Technical check 
VlnPlot(object = pbmc, features = c("nFeature_RNA"), pt.size = 0)+theme(text = element_text(size=12))+theme(axis.title.x=element_blank(),axis.title.y=element_text(color="black", size=12, face="bold"),axis.text.x=element_text(color="black", size=12, face="bold",angle=0,hjust=0.5),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=12, face="bold"), strip.text = element_text(size = 12, face="bold"))+NoLegend()+ggtitle("")

VlnPlot(object = pbmc, features = c("percent.mito"), pt.size = 0)+theme(text = element_text(size=12))+theme(axis.title.x=element_blank(),axis.title.y=element_text(color="black", size=12, face="bold"),axis.text.x=element_text(color="black", size=12, face="bold",angle=0,hjust=0.5),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=12, face="bold"), strip.text = element_text(size = 12, face="bold"))+NoLegend()+ggtitle("")

VlnPlot(object = pbmc, features = c("nCount_RNA") , pt.size = 0)+theme(text = element_text(size=12))+theme(axis.title.x=element_blank(),axis.title.y=element_text(color="black", size=12, face="bold"),axis.text.x=element_text(color="black", size=12, face="bold",angle=0,hjust=0.5),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=12, face="bold"), strip.text = element_text(size = 12, face="bold"))+NoLegend()+ggtitle("")

#Figure 1A
#UMAP plot by cell type
DimPlot(pbmc, group.by="orig.ident")+theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))+NoLegend()+scale_color_manual(values=c("#0000FF","#FF0000"))+scale_fill_manual(values=c("#0000FF","#FF0000"))

#Figure 1B
#UMAP plot by unsupervised clustering 
DimPlot(pbmc)+theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))+NoLegend()

#Figure 1C
#PC plot by unsupervised clustering 
DimPlot(pbmc,reduction="pca")+theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))+NoLegend()+xlab("PC1 (20.3% of variance explained)")+ylab("PC2 (11.1% of variance explained)")

#Figure 1D
#cell-to-cell Correlation
my_sample_col=as.data.frame(cbind(as.data.frame(pbmc@meta.data$orig.ident),as.data.frame(pbmc@active.ident)))
colnames(my_sample_col)=c("Type","Cluster")
Var1 = c("#0000FF","#FF0000")
names(Var1) = c("B16F0", "B16F10")
Var2 = c("#F8766D","#00BFC4")
names(Var2) = c("F0", "F10")
ann_colors = list(Type = Var1, Cluster=Var2)
my_sample_col2=my_sample_col[order(my_sample_col$Cluster),]
write_out= as.data.frame(as.matrix(pbmc@assays$RNA@scale.data))
write_out2=write_out[pbmc@assays$RNA@var.features,]
write_out3=write_out2[,order(my_sample_col$Cluster)]
dat.n=round(cor(write_out3),digits=2)
rownames(my_sample_col2)=rownames(dat.n)
pheatmap(dat.n,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(-0.3, 0.3, by = 0.006),annotation_col =my_sample_col2,annotation_row =my_sample_col2,annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_colors=ann_colors,colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100))
#F0 average correltaion calculation 
dat.n3=cor(write_out3)[1:2504,1:2504]
dat.n3[dat.n3==1] <- NA
mean(dat.n3,na.rm=TRUE)
#F10 average correltaion calculation 
dat.n3=cor(write_out3)[2505:4156,2505:4156]
dat.n3[dat.n3==1] <- NA
mean(dat.n3,na.rm=TRUE)


#FigureS2
#TPM transformation of average of scRNA-seq
avg_exp=as.data.frame(cbind(as.data.frame(rownames(pbmc@assays$RNA@counts)),rowMeans(pbmc@assays$RNA@counts[,1:2582]),rowMeans(pbmc@assays$RNA@counts[,2583:4156])))
colnames(avg_exp)=c("Gene","F0_avg","F10_avg")
avg_exp$F0_avg=log2(avg_exp$F0_avg/sum(avg_exp$F0_avg)*1000000+1)
avg_exp$F10_avg=log2(avg_exp$F10_avg/sum(avg_exp$F10_avg)*1000000+1)

#TPM from bulk RNA-seq was calculatd by Stringtie
F0_exp=read.table("/data/MRC1_data4/kysbbubbu/B16/bulk_RNA_B16/04_STRINGTIE/B16F0/B16F0.txt",sep="\t",header=TRUE)
F10_exp=read.table("/data/MRC1_data4/kysbbubbu/B16/bulk_RNA_B16/04_STRINGTIE/B16F10/B16F10.txt",sep="\t",header=TRUE)
colnames(F10_exp)
bulk_exp=merge(F0_exp, F10_exp, by="Gene.Name", all = FALSE)
?merge
bulk_exp=as.data.frame(cbind(as.data.frame(bulk_exp[,1]),bulk_exp[,9],bulk_exp[,17]))
colnames(bulk_exp)=c("Gene","F0","F10")
bulk_exp$F0=log2(bulk_exp$F0+1)
bulk_exp$F10=log2(bulk_exp$F10+1)
exp=merge(avg_exp, bulk_exp, by="Gene", all = FALSE)

#F0 (FigureS2A)
ggplot(exp, aes(F0_avg, F0)) + geom_point(color="#0000FF")  +  theme_classic() +theme(text = element_text(size=12))+theme(axis.text.x = element_text(size=12,face="bold"),axis.text.y = element_text(size=12, face="bold"), axis.title.x = element_text(size=12,face="bold"),axis.title.y = element_text(size=12, face="bold"))+ theme(legend.position = "none")+ geom_smooth(method=lm,linetype="dotted", color="black")+xlab("Average of single-cell expressions [log2(TPM+1)]")+ylab("Bulk cells [log2(TPM+1)]")+ylim(0,15)+xlim(0,15)
cor.test(exp$F0_avg,exp$F0)

#F10 (FigureS2B)
ggplot(exp, aes(F10_avg, F10)) + geom_point(color="#FF0000")  +  theme_classic() +theme(text = element_text(size=12))+theme(axis.text.x = element_text(size=12,face="bold"),axis.text.y = element_text(size=12, face="bold"), axis.title.x = element_text(size=12,face="bold"),axis.title.y = element_text(size=12, face="bold"))+ theme(legend.position = "none")+ geom_smooth(method=lm,linetype="dotted", color="black")+xlab("Average of single-cell expressions [log2(TPM+1)]")+ylab("Bulk cells [log2(TPM+1)]")+ylim(0,15)+xlim(0,15)
cor.test(exp$F10_avg,exp$F10)

### Differentially expressed genes

F0F10.markers <- FindAllMarkers(object = pbmc, return.thresh = 0.01, logfc.threshold = 0.25, only.pos = TRUE, min.pct = 0.1)
F0F10.markers=as.data.frame(F0F10.markers)
#Table S3
write(x = as.data.frame(F0F10.markers), row.names = TRUE, file = "F10F10.markers.csv")
dat.n=pbmc@assays$RNA@scale.data
F0F10.markers2=F0F10.markers[F0F10.markers$cluster=="F0",]
F0F10.markers3=F0F10.markers[F0F10.markers$cluster=="F10",]
F0F10.markers2=F0F10.markers2[order(F0F10.markers2[,6], F0F10.markers2[,2],decreasing=TRUE),]
F0F10.markers3=F0F10.markers3[order(F0F10.markers3[,6], F0F10.markers3[,2],decreasing=TRUE),]
marker=c(as.character(F0F10.markers2$gene),as.character(F0F10.markers3$gene))
dat.n2=dat.n[marker,]
my_sample_col=as.data.frame(cbind(as.data.frame(pbmc@meta.data$orig.ident),as.data.frame(pbmc@active.ident)))
colnames(my_sample_col)=c("Type","Cluster")
rownames(my_sample_col)=rownames(pbmc@meta.data)
Var1 = c("#0000FF","#FF0000")
names(Var1) = c("B16F0", "B16F10")
Var2 = c("#F8766D","#00BFC4")
names(Var2) = c("F0", "F10")
ann_colors = list(Type = Var1, Cluster=Var2)
avg.order=as.data.frame(colMeans(dat.n2[1:143,]))-as.data.frame(colMeans(dat.n2[144:274,]))
colnames(avg.order)="avg"
dat.n2=dat.n2[,order(avg.order,decreasing=TRUE)]
#Figure2A
pheatmap(dat.n2,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(-1.5, 1.5, by = 0.03),annotation_col =my_sample_col,annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_colors=ann_colors,colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100))

#Figure 2B
markers.to.plot=c("Ptgds", "Cyb5a", "Cd63", "B2m", "Lgals3", "Sparc","Met", "Tmsb4x", "Grn",  "Apoe", "Mlana")
DotPlot(pbmc, features = markers.to.plot, cols = c("grey", "red"),dot.min=0)+theme(axis.title=element_blank(), axis.text.x.bottom = element_text(size=16,face="bold"), axis.text.y.left=element_text(size = 16,face="bold")) +coord_flip()

#Figure 2C
s.gene = c("Mcm5","Pcna-ps2","Tyms","Fen1","Mcm2","Mcm4","Rrm1","Ung","Gins2","Mcm6","Cdca7","Dtl","Prim1","Uhrf1","Cenpu","Hells","Rfc2","Rpa2","Nasp","Rad51ap1","Gmnn","Wdr76","Slbp","Ccne2","Ubr7","Pold3","Msh2","Atad2","Rad51","Rrm2","Cdc45","Cdc6","Exo1","Tipin","Dscc1","Blm","Casp8ap2","Usp1","Clspn","Pola1","Chaf1b","Brip1","E2f8")
g2m.gene = c("Hmgb2","Cdk1","Nusap1","Ube2c","Birc5","Tpx2","Top2a","Ndc80","Cks2","Nuf2","Cks1b","Mki67","Tmpo","Cenpf","Tacc3","Smc4","Ccnb2","Ckap2l","Ckap2","Aurkb","Bub1","Kif11","Anp32e","Tubb4b","Gtse1","Kif20b","Hjurp","Cdca3","Cdc20","Ttk","Cdc25c","Kif2c","Rangap1","Ncapd2","Dlgap5","Cdca2","Cdca8","Ect2","Kif23","Hmmr","Aurka","Psrc1","Anln","Lbr","Ckap5","Cenpe","Ctcf","Nek2","G2e3","Gas2l3","Cbx5","Cenpa")
pbmc <- CellCycleScoring(object = pbmc, s.features = s.gene, g2m.features = g2m.gene, set.ident = FALSE)
rect2 <- data.frame (xmin=0, xmax=Inf, ymin=-Inf, ymax=Inf)
rect3 <- data.frame (xmin=-Inf, xmax=0, ymin=0, ymax=Inf)

FeatureScatter(object = pbmc, feature1 = "G2M.Score", feature2 = "S.Score",pt.size=0.7)+theme(plot.title = element_blank(), axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))+ geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="blue", alpha=0.1, inherit.aes = FALSE)+geom_rect(data=rect2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="red", alpha=0.1, inherit.aes = FALSE)+geom_rect(data=rect3, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="red", alpha=0.1, inherit.aes = FALSE)+xlab("Average expression of G2/M phase genes")+ylab("Average expression of S phase genes")+NoLegend()

FeatureScatter(object = pbmc, feature1 = "G2M.Score", feature2 = "S.Score",group.by="kdm5b",pt.size=0.7,cols=c("red","grey"))+theme(plot.title = element_blank(), axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))+ geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="blue", alpha=0.1, inherit.aes = FALSE)+geom_rect(data=rect2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="red", alpha=0.1, inherit.aes = FALSE)+geom_rect(data=rect3, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="red", alpha=0.1, inherit.aes = FALSE)+xlab("Average expression of G2/M phase genes")+ylab("Average expression of S phase genes")+NoLegend()

cluster1=c("Non-cycling","Non-cycling","Cycling", "Cycling")
cluster2=c("F0","F10","F0","F10")
value =c(1176,339,1328,1313)
data <- data.frame(cluster1,cluster2,value)
ggplot(data, aes(fill=cluster1, y=value, x=cluster2)) + geom_bar(position="fill", stat="identity")+ geom_bar(position="fill", stat="identity")+theme_classic()+theme(plot.title = element_blank(), axis.title.y =element_blank(), axis.title.x = element_text(size=14,face="bold"),axis.text.x = element_text(size=12,face="bold"), axis.text.y = element_text(size=12,angle=90,hjust=0.5,face="bold"),legend.title=element_blank(),legend.text=element_text(size=12))+ylab("Proportion")+ coord_flip()+scale_x_discrete(limits = c("F10","F0"))+scale_fill_manual(values=c("#FF9999","#99CCFF"))


cluster1=c("Non-cycling","Non-cycling","Cycling", "Cycling")
cluster2=c("Kdm5b high","Kdm5b low", "Kdm5b high", "Kdm5b low")
value =c(168,1347,64,2577)
data <- data.frame(cluster1,cluster2,value)
ggplot(data, aes(fill=cluster2, y=value, x=cluster1)) + geom_bar(position="fill", stat="identity")+ geom_bar(position="fill", stat="identity")+theme_classic()+theme(plot.title = element_blank(), axis.title.y =element_blank(), axis.title.x = element_text(size=14,face="bold"),axis.text.x = element_text(size=12,face="bold"), axis.text.y = element_text(size=12,angle=90,hjust=0.5,face="bold"),legend.title=element_blank(),legend.text=element_text(size=12))+ylab("Proportion")+ coord_flip()+scale_fill_manual(values=c("red","grey"))

cluster1=c("F0","F0","F10","F10")
cluster2=c("Kdm5b high","Kdm5b low", "Kdm5b high", "Kdm5b low")
value =c(114,2390,118,1534)
data <- data.frame(cluster1,cluster2,value)
cluster1=factor(cluster1,levels=c("F0","F10"))
ggplot(data, aes(fill=cluster2, y=value, x=cluster1)) + geom_bar(position="fill", stat="identity")+theme_classic()+theme(plot.title = element_blank(), axis.title.y =element_blank(), axis.title.x = element_text(size=12,face="bold"),axis.text.x = element_text(size=12,face="bold"), axis.text.y = element_text(size=12,angle=90,hjust=0.5,face="bold"),legend.title=element_blank(),legend.text=element_text(size=12))+ylab("Proportion")+ coord_flip()+scale_x_discrete(limits = c("F10","F0"))+scale_fill_manual(values=c("red","grey"))


###Subpopulation with intermediate characteristics of B16F0 and B16F10 cells 

#Figure S4A
#Unsupervised clustering using Various resolution test 
pbmc <- FindClusters(object = pbmc, resolution = 0.1)
DimPlot(pbmc, label=TRUE)+NoLegend()+theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))+ggtitle("Resolution=0.10")+ theme(plot.title = element_text(hjust = 0.5))
pbmc <- FindClusters(object = pbmc, resolution = 0.15)
DimPlot(pbmc, label=TRUE)+NoLegend()+theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))+ggtitle("Resolution=0.15")+ theme(plot.title = element_text(hjust = 0.5))
pbmc <- FindClusters(object = pbmc, resolution = 0.2)
DimPlot(pbmc, label=TRUE)+NoLegend()+theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))+ggtitle("Resolution=0.20")+ theme(plot.title = element_text(hjust = 0.5))
pbmc <- FindClusters(object = pbmc, resolution = 0.25)
DimPlot(pbmc, label=TRUE)+NoLegend()+theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))+ggtitle("Resolution=0.25")+ theme(plot.title = element_text(hjust = 0.5))
pbmc <- FindClusters(object = pbmc, resolution = 0.3)
DimPlot(pbmc, label=TRUE)+NoLegend()+theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))+ggtitle("Resolution=0.30(equal to C1-5 clusters)")+ theme(plot.title = element_text(hjust = 0.5))

#Kmeans silouette width calculation 
#Figure S4B
write_out= as.data.frame(as.matrix(pbmc@assays$RNA@scale.data))[pbmc@assays$RNA@var.features,]
a=fviz_nbclust(write_out, kmeans, method = "silhouette")+  labs(subtitle = "Silhouette method")
sil=as.data.frame(cbind(a$data$clusters, a$data$y)[2:8,])
sil$V1=factor(sil$V1)
sil$V2=round(sil$V2,3)
ggplot(sil, aes(y=V2, x=V1)) +  geom_bar(stat="identity", width=0.7, fill="steelblue")+theme_classic()+ theme(plot.title = element_blank(), axis.title.y =element_text(size=14,face="bold"), axis.title.x = element_text(size=14,face="bold"),axis.text.x = element_text(size=12,face="bold"), axis.text.y = element_text(size=12,hjust=0.5,face="bold"),legend.title=element_blank(),legend.text=element_text(size=12))+ylab("Average silhouette width")+ylim(0,0.06)+xlab("Number of cluster")+geom_text(aes(label=V2), hjust=0.5,vjust=2,position = position_dodge(1), size=4,color="white")

#ROGUE calculation 
#Figure S4C

current.cluster.ids <- c(0,1)
new.cluster.ids <- c("F0","F10")
pbmc@meta.data$RNA_snn_res.0.05  <- plyr::mapvalues(x = pbmc@meta.data$RNA_snn_res.0.05, from = current.cluster.ids, to = new.cluster.ids)

expr <- matr.filter(as.data.frame(pbmc@assays$RNA@counts), min.cells = 1, min.genes = 1)
ent.res <- SE_fun(expr)
rogue.value <- CalculateRogue(ent.res, platform = "UMI")
rogue.res <- rogue(expr, labels = pbmc@meta.data$RNA_snn_res.0.3, samples = pbmc@meta.data$sample, platform = "UMI", span = 0.6)
rogue.res2 <- rogue(expr, labels = pbmc@meta.data$RNA_snn_res.0.05, samples = pbmc@meta.data$sample, platform = "UMI", span = 0.6)
pbmc@meta.data$sample=rep("B16",4156)
rogue.res3 <- rogue(expr, labels = pbmc@meta.data$sample, samples = pbmc@meta.data$sample, platform = "UMI", span = 0.6)
my_data <- data.frame( 
  clusterp = c("Total", "F0","F10","C1","C2","C3","C4","C5"),
  ROGUE = c(0.545,0.562,0.572,0.629,0.669,0.602,0.562,0.666)
)
my_data$clusterp=factor(my_data$cluster,level=c("Total", "F0","F10","C1","C2","C3","C4","C5"))
ggplot(my_data, aes(y=ROGUE, x=clusterp)) +  geom_bar(stat="identity", width=0.7, fill="steelblue")+theme_classic()+ theme(plot.title = element_blank(), axis.title.y =element_text(size=14,face="bold"), axis.title.x = element_text(size=14,face="bold"),axis.text.x = element_text(size=12,face="bold"), axis.text.y = element_text(size=12,hjust=0.5,face="bold"),legend.title=element_blank(),legend.text=element_text(size=12))+ylab("ROGUE (Metric of homogeneity in each cluster)")+ ylim(0,0.8)+xlab("Cluster")+geom_text(aes(label=ROGUE), hjust=0.5,vjust=2,position = position_dodge(1), size=4,color="white")


#Figure 3A
pbmc <- FindClusters(object = pbmc, resolution = 0.3)
new.cluster.ids <- c("C5", "C3", "C1", "C4", "C2")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
my_levels <- c("C1","C2","C3","C4","C5")
UMAPPlot(pbmc) + NoLegend() +theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))

#Figure 3B
table(pbmc@active.ident,pbmc@meta.data$orig.ident)
cluster1= rep(c("B16F0","B16F10"),5)
cluster2=rep(c("C1","C2","C3","C4","C5"),each=2)
value =c(816,7,644,4,640,306,469,212,13,1045)
data <- data.frame(cluster1,cluster2,value)
ggplot(data, aes(fill=cluster1, y=value, x=cluster2)) + geom_bar(position="fill", stat="identity")+ geom_bar(position="fill", stat="identity")+theme_classic()+theme(plot.title = element_blank(), axis.title.y =element_blank(), axis.title.x = element_text(size=14,face="bold"),axis.text.x = element_text(size=12,face="bold"), axis.text.y = element_text(size=12,angle=90,hjust=0.5,face="bold"),legend.title=element_blank(),legend.text=element_text(size=12))+ylab("Proportion")+scale_fill_manual(values=c("#0000FF","#FF0000"))+ NoLegend()


#Figure 3C
current.cluster.ids <- c(0,1,2,3,4)
new.cluster.ids <- c("C5","C3","C1","C4","C2")
pbmc@meta.data$RNA_snn_res.0.3  <- plyr::mapvalues(x = pbmc@meta.data$RNA_snn_res.0.3, from = current.cluster.ids, to = new.cluster.ids)
pbmc@meta.data$RNA_snn_res.0.3=factor(pbmc@meta.data$RNA_snn_res.0.3,levels=c("C1","C2","C3","C4","C5"))
my_sample_col=as.data.frame(cbind(as.data.frame(pbmc@meta.data$orig.ident),as.data.frame(pbmc@meta.data$RNA_snn_res.0.3)))
colnames(my_sample_col)=c("Type","Cluster")
my_sample_col$Cluster=factor(my_sample_col$Cluster,level=c("C1","C2","C3","C4","C5"))
Var1 = c("#0000FF","#FF0000")
names(Var1) = c("B16F0", "B16F10")
Var2 = c("#F8766D","#A3A500","#00BF7D","#00B0F6","#E76BF3")
names(Var2) = c("C1", "C2","C3","C4","C5")
ann_colors = list(Type = Var1, Cluster=Var2)
my_sample_col2=my_sample_col[order(my_sample_col$Cluster),]
write_out= as.data.frame(as.matrix(pbmc@assays$RNA@scale.data))
write_out2=write_out[pbmc@assays$RNA@var.features,]
write_out4=write_out2[,order(my_sample_col$Cluster)]
dat.n2=round(cor(write_out4),digits=2)
rownames(my_sample_col2)=rownames(dat.n2)
pheatmap(dat.n2,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(-0.3, 0.3, by = 0.006),annotation_col =my_sample_col2,annotation_row =my_sample_col2,annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_colors=ann_colors,color=colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100))

#mean per cluster in Figure 3C 
dat.n3=cor(write_out4)[1:823,1:823]
dat.n3[dat.n3==1] <- NA
mean(dat.n3,na.rm=TRUE)
dat.n3=cor(write_out4)[824:1471,824:1471]
dat.n3[dat.n3==1] <- NA
mean(dat.n3,na.rm=TRUE)
dat.n3=cor(write_out4)[1472:2417,1472:2417]
dat.n3[dat.n3==1] <- NA
mean(dat.n3,na.rm=TRUE)
dat.n3=cor(write_out4)[2418:3098,2418:3098]
dat.n3[dat.n3==1] <- NA
mean(dat.n3,na.rm=TRUE)
dat.n3=cor(write_out4)[3099:4156,3099:4156]
dat.n3[dat.n3==1] <- NA
mean(dat.n3,na.rm=TRUE)


#TableS5
int.marker <- FindMarkers(object = pbmc, ident.1 = c("C3","C4"), ident.2 = c("C1","C2","C5"), return.thresh = 0.01, logfc.threshold = 0.25, only.pos = TRUE, min.pct = 0.1)
fwrite(x = as.data.frame(int.marker), row.names = TRUE, file = "C34_marker.csv")
int.marker2 <- FindMarkers(object = pbmc, ident.1 = c("C5"), ident.2 = c("C1","C2","C3","C4"), return.thresh = 0.01, logfc.threshold = 0.25, only.pos =TRUE, min.pct = 0.1)
fwrite(x = as.data.frame(int.marker2), row.names = TRUE, file = "C5_marker.csv")


#Figure 3D
markers.to.plot=c("Tyr", "Apoe","Pmel", "Sparc", "Serpinh1","Id2", "Tmsb4x", "Met",  "Hras", "Ccnd1", "B2m","Cdkn1a")
DotPlot(pbmc, features = markers.to.plot, cols = c("grey", "red"),dot.min=0)+theme(axis.title=element_blank(), axis.text.x.bottom = element_text(size=16,face="bold"), axis.text.y.left=element_text(size = 16,face="bold")) +coord_flip()+NoLegend()

#Figure 3E 
scale.data = as.data.frame(pbmc@assays$RNA@scale.data)
geneSets <- getGmt("gsva.gmt")
esrnaseq <- gsva(pbmc@assays$RNA@scale.data, geneSets, method="gsva", kcdf="Gaussian", mx.diff=TRUE, verbose=TRUE, parallel.sz=24)
GSVA=as.data.frame(esrnaseq)
rownames(GSVA)
pbmc@meta.data$GSVA_MITF=t(GSVA[1,])
VlnPlot(object = pbmc, features = c("GSVA_MITF"), pt.size = 0)+ theme(plot.title = element_text(family = "arial", face = "bold", hjust = 0.5, size = 24)) + NoLegend() + xlab(" ") + theme(axis.text.x=element_text(angle=0, hjust=0.5)) + ggtitle("") +ylab("GSVA z score ")+theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))

#Figure 3F
cluster1=c(rep("Non-cycling",5),rep("Cycling",5))
cluster2=c(rep(c("C1","C2","C3","C4","C5"),2))
value =c(602,3,849,2,59,221,645,97,679,999)
data <- data.frame(cluster1,cluster2,value)
ggplot(data, aes(fill=cluster1, y=value, x=cluster2)) + geom_bar(position="fill", stat="identity")+ geom_bar(position="fill", stat="identity")+theme_classic()+theme(plot.title = element_blank(), axis.title.y =element_text(size=14,face="bold"), axis.title.x = element_text(size=14,face="bold"),axis.text.x = element_text(size=12,face="bold"), axis.text.y = element_text(size=12,hjust=0.5,face="bold"),legend.title=element_blank(),legend.text=element_text(size=12))+xlab("Cluster")+ylab("Proportion")+NoLegend()+scale_fill_manual(values=c("#FF9999","#99CCFF"))

#Figure S5
dat.n=pbmc@assays$RNA@scale.data
F0F10.markers2=F0F10.markers[F0F10.markers$cluster=="F0",]
F0F10.markers3=F0F10.markers[F0F10.markers$cluster=="F10",]
F0F10.markers2=F0F10.markers2[order(F0F10.markers2[,6], F0F10.markers2[,2],decreasing=TRUE),]
F0F10.markers3=F0F10.markers3[order(F0F10.markers3[,6], F0F10.markers3[,2],decreasing=TRUE),]
marker=combine(as.character(F0F10.markers2$gene),as.character(F0F10.markers3$gene))
dat.n2=dat.n[marker,]
my_sample_col=as.data.frame(cbind(as.data.frame(pbmc@meta.data$RNA_snn_res.0.05),as.data.frame(pbmc@meta.data$RNA_snn_res.0.3)))
colnames(my_sample_col)=c("Cluster1","Cluster2")
rownames(my_sample_col)=rownames(pbmc@meta.data)
Var2 = c("#F8766D","#00CCCC")
names(Var2) = c("F0", "F10")
Var3 = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")
names(Var3) = c("C1", "C2","C3","C4","C5")
ann_colors = list(Cluster1=Var2,Cluster2=Var3)
avg.order=as.data.frame(colMeans(dat.n2[1:143,]))-as.data.frame(colMeans(dat.n2[144:274,]))
colnames(avg.order)="avg"
dat.n2=dat.n2[,order(avg.order,decreasing=TRUE)]
pheatmap(dat.n2,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(-1.5, 1.5, by = 0.03),annotation_col =my_sample_col,annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_colors=ann_colors,colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100))

#Figure S6A
FeaturePlot(object = pbmc, features = c("S.Score"))+scale_colour_gradientn(colours=brewer.pal(9,"OrRd")[1:9])+labs(title="Average expression of S phase genes") + theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))
FeaturePlot(object = pbmc, features = c("G2M.Score"))+ scale_colour_gradientn(colours=brewer.pal(9,"OrRd")[1:9])+labs(title="Average expression of G2/M phase genes")+theme(plot.title = element_text(family = "arial", face = "bold", hjust = 0.5, size = 15)) + theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))

#Figure S6B
#kdm5b status was assigned by scaled expression of kdm5b (higher or less than z score 2)
my_sample_col=as.data.frame(cbind(as.data.frame(pbmc@meta.data$RNA_snn_res.0.05),as.data.frame(pbmc@meta.data$kdm5b),as.data.frame(pbmc@meta.data$RNA_snn_res.0.3)))
colnames(my_sample_col)=c("Cluster","Kdm5b","Cluster2")
rownames(my_sample_col)=rownames(pbmc@meta.data)
colnames(my_sample_col)
Var1 = c("#F8766D","#00CCCC")
names(Var1) = c("F0", "F10")
Var2 = c("red","grey")
names(Var2) = c("Kdm5b high", "Kdm5b low")
Var3 = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")
names(Var3) = c("C1", "C2","C3","C4","C5")
ann_colors = list(Cluster = Var1, Kdm5b=Var2, Cluster2=Var3)
my_sample_col2=my_sample_col[order(my_sample_col$Cluster),]
s.g2m.gene=c(s.gene,g2m.gene)
write_out= as.data.frame(as.matrix(pbmc@assays$RNA@scale.data))
write_out2=write_out[s.g2m.gene,]
avg.order2=as.data.frame(colMeans(write_out2))
dat.n=write_out2[,order(avg.order2)]
pheatmap(dat.n,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(-1.5, 1.5, by = 0.03),annotation_col =my_sample_col,annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_colors=ann_colors,colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100))


###Inferring single-cell trajectory from B16F0 to B16F10 

data <- as(as.matrix(pbmc@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = as.data.frame(pbmc@meta.data))
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
HSMM <- newCellDataSet(data,phenoData = pd,featureData = fd,lowerDetectionLimit = 0.5,  expressionFamily = VGAM::negbinomial.size())
HSMM <- detectGenes(HSMM, min_expr = 0.1)
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
HSMM_expressed_genes <-  row.names(subset(fData(HSMM),num_cells_expressed >= 5))
clustering_DEG_genes <- differentialGeneTest(HSMM[HSMM_expressed_genes,], fullModelFormulaStr = '~RNA_snn_res.0.3',cores = 24)
HSMM_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:500]
HSMM <-  setOrderingFilter(HSMM,ordering_genes = HSMM_ordering_genes)
HSMM <- reduceDimension(HSMM, method = 'DDRTree')
HSMM <- orderCells(HSMM)
pData(HSMM)$RNA_snn_res.0.3=factor(pData(HSMM)$RNA_snn_res.0.3,levels=c("C1","C2","C3","C4","C5"))


#Figure4A
plot_cell_trajectory(HSMM, cell_size = 1, color_by = "RNA_snn_res.0.3")+ theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_blank(),legend.text=element_text(size=14))

#Figure S7A
plot_cell_trajectory(HSMM, cell_size = 1, color_by = "orig.ident")+ theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=14),legend.title=element_blank(),legend.text=element_text(size=12))+scale_colour_manual(name="orig.ident", values= c("blue", "red"))

#Figure S7B
plot_cell_trajectory(HSMM, color_by = "State") + theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=13), axis.text.y = element_text(size=13),legend.title=element_text(size=14),legend.text=element_text(size=12))

#Figure S7C
pData(HSMM)$Phase2=pbmc@meta.data$Phase
current.cluster.ids <- c("G1","G2M","S")
new.cluster.ids <- c("Non-cycling","Cycling","Cycling")
pData(HSMM)$Phase2  <- plyr::mapvalues(x = pData(HSMM)$Phase2, from = current.cluster.ids, to = new.cluster.ids)
pData(HSMM)$Phase2=factor(pData(HSMM)$Phase2,levels=c("Cycling","Non-cycling"))
plot_cell_trajectory(HSMM, cell_size = 1, color_by = "Phase2")+ theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),legend.title=element_blank(),legend.text=element_text(size=14))+scale_colour_manual(values=c("#FF9999","#99CCFF"))

HSMM <- orderCells(HSMM, root_state = 4)
#Figure S7D
plot_cell_trajectory(HSMM, cell_size = 1, color_by = "Pseudotime")+ scale_colour_gradientn(colours=brewer.pal(9,"YlOrRd")[1:9]) + theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))


#RNA velocity analysis 
ldat <- read.loom.matrices("B16merge.loom")
emat <- ldat$spliced
nmat <- ldat$unspliced
gp=plot_cell_trajectory(HSMM, cell_size = 1, color_by = "seurat_clusters")
emb <- as.data.frame(cbind(gp$data$data_dim_1, gp$data$data_dim_2))
rownames(emb)=rownames(pbmc@meta.data)
cell.dist <- as.dist(1-armaCor(t(emb)))
emat.modify=as.colnames(emat)
fwrite(x = as.data.frame(as.matrix(emat.modify)), row.names = TRUE, file = "emat.csv")
nmat.modify=colnames(nmat)
fwrite(x = as.data.frame(as.matrix(nmat.modify)), row.names = TRUE, file = "nmat.csv")
# Manually modify cell names in emat.csv / nmat.csv same as Seurat  
# Modified file saved as  ematcell.csv / nmatcell.csv
MyData <- read.csv(file="nmatcell.csv", header=FALSE, sep=",", stringsAsFactors = FALSE)
MyData2 <- read.csv(file="ematcell.csv", header=FALSE, sep=",", stringsAsFactors = FALSE)
nmat.modify2=MyData$V1
colnames(emat) = nmat.modify2
colnames(nmat) = nmat.modify2
embeddings=as.matrix(pbmc@reductions$umap@cell.embeddings)
emat.new = emat[,colnames(emat) %in% rownames(embeddings)]
nmat.new = nmat[,colnames(nmat) %in% rownames(embeddings)]
emat.new2 = emat.new[rownames(emat.new) %in% rownames(pbmc@assays$RNA@data),]
nmat.new2 = nmat.new[rownames(emat.new) %in% rownames(pbmc@assays$RNA@data),]
rvel.cd_monocle <- gene.relative.velocity.estimates(emat.new2,nmat.new2,deltaT=1,kCells=20,cell.dist=cell.dist,fit.quantile=0.05, n.cores=24,verbose=TRUE,min.nmat.emat.correlation = 0.2, min.nmat.emat.slope = 0.2)
gp = plot_cell_trajectory(HSMM, cell_size = 1, color_by = "seurat_clusters")
colors <- as.list(ggplot_build(gp)$data[[2]]$colour)
names(colors) <- rownames(emb)
par(oma=c(2,2,2,2))
par(mar=c(2,2,2,2))
par(mfrow=c(1,1))
emb2 = cbind(ggplot_build(gp)$data[[2]][,"x"],ggplot_build(gp)$data[[2]][,"y"])
rownames(emb2)=names(colors)
#Figure4B
show.velocity.on.embedding.cor(emb2,rvel.cd_monocle,n=300,scale='sqrt',cell.colors=ac(colors,alpha=0.5),cex=0.8,arrow.scale=3,show.grid.flow=T,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=F,cell.border.alpha = 0.1,n.cores=24)

#Figure 4C 
HSMM_state1 <- row.names(subset(pData(HSMM),State %in% c("1","2","3","4")))
data=pData(HSMM[,HSMM_state1])
data=data[order(data$Pseudotime),]
data$pseudo_bin=rep(c("B1","B2","B3","B4"), c(549,315,120,1170))

ggplot(data, aes(pseudo_bin, fill = seurat_clusters)) +  geom_bar(position = "fill") +  scale_y_continuous(labels = scales::percent)+theme_classic()+theme(plot.title = element_blank(), axis.title.x =element_blank(), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12,face="bold"), axis.text.y = element_text(size=12,hjust=0.5,face="bold"),legend.title=element_blank(),legend.text=element_text(size=12))+ylab("Proportion")

#Figure 4D
heatmap=plot_pseudotime_heatmap(HSMM[row.names(subset(diff_test_res_pseudotime_state, qval < 1e-100)),HSMM_state1],cores = 16,norm_method ="log",use_gene_short_name = T, show_rownames = FALSE, scale_max = 3, scale_min = -3, num_clusters = 2,return_heatmap = TRUE,hmcols=colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100))
?plot_pseudotime_heatmap
heatmap_gene=as.data.frame(rownames(HSMM[row.names(subset(diff_test_res_pseudotime_state, qval < 1e-100)),HSMM_state1][heatmap$tree_row[["order"]],]))
colnames(heatmap_gene)="gene"
heatmap_gene$cluster=rep(c("down_regulated", "up_regulated"), c(79,63))
#TableS6
fwrite(heatmap_gene , file = "monocle_gene_cluster.csv", append = FALSE, row.names = T, col.names = TRUE)
#Figure 4E

pbmc@meta.data$Pseudotime = pData(HSMM)$Pseudotime
pbmc@meta.data$State = pData(HSMM)$State
FeatureScatter(object = pbmc.state, feature1 = "Pseudotime", feature2 = "Apoe",pt.size=0.7)+theme(plot.title = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),legend.title=element_text(size=14),legend.text=element_text(size=12))+NoLegend()+ geom_smooth(method = 'loess', aes(colour = NA), colour = 'black')
FeatureScatter(object = pbmc.state, feature1 = "Pseudotime", feature2 = "Lgals3",pt.size=0.7)+theme(plot.title = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),legend.title=element_text(size=14),legend.text=element_text(size=12))+NoLegend()+ geom_smooth(method = 'loess', aes(colour = NA), colour = 'black')
FeatureScatter(object = pbmc.state, feature1 = "Pseudotime", feature2 = "Met",pt.size=0.7)+theme(plot.title = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),legend.title=element_text(size=14),legend.text=element_text(size=12))+NoLegend()+ geom_smooth(method = 'loess', aes(colour = NA), colour = 'black')
FeatureScatter(object = pbmc.state, feature1 = "Pseudotime", feature2 = "Tmsb4x",pt.size=0.7)+theme(plot.title = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),legend.title=element_text(size=14),legend.text=element_text(size=12))+NoLegend()+ geom_smooth(method = 'loess', aes(colour = NA), colour = 'black')

### Genomic alterations from B16F0 to B16F10

#Figure 5A
#prepare input files required by infercnv packages 
infercnv_obj5 = CreateInfercnvObject(raw_counts_matrix="countmatrix.txt", annotations_file="ref.txt", delim="\t",gene_order_file="gtf.txt",ref_group_names=c("C1","C2"))
infercnv_obj5 = infercnv::run(infercnv_obj5, cutoff=0.1,  out_dir="NEWC1C2_denoise",cluster_by_groups=FALSE, HMM=TRUE,denoise=TRUE,k_obs_groups=2, hclust_method="ward.D2")

#Figure S8, upper panel
infercnv_obj4 = CreateInfercnvObject(raw_counts_matrix="countmatrix.txt",annotations_file="ref2.txt",delim="\t",gene_order_file="gtf.txt",ref_group_names=c("B16F0"))
infercnv_obj4 = infercnv::run(infercnv_obj4,cutoff=0.1,out_dir="NEWF0F10_group",cluster_by_groups=FALSE)

#Figure 5B
Type <- c(rep("B16F0" ,2) , rep("B16F10" , 2) )
condition <- rep(c("Shared" , "Private" ) , 2)
condition=factor(condition, levels=c("Shared" , "Private"))
value <- c(2087,636,2087,1105)
data <- data.frame(Type,condition,value)
ggplot(data, aes(fill=condition, y=value, x=Type)) +   geom_bar(position="stack", stat="identity")+theme_classic()+ theme(plot.title = element_blank(), axis.title.y =element_text(size=13,face="bold"), axis.title.x = element_blank(),axis.text.x = element_text(size=12,face="bold"), axis.text.y = element_text(size=12,hjust=0.5,face="bold"),legend.title=element_blank(),legend.text=element_text(size=12))+ylab("Number of mutation")+xlab("Cluster")+geom_text(aes(x = Type, y = value, label = value, group = condition),  position = position_stack(vjust = .5), size=4,color="white")
Type <- c(rep("B16F0" ,2) , rep("B16F10" , 2) )
condition <- rep(c("Shared" , "Private" ) , 2)
condition=factor(condition, levels=c("Shared", "Private" ))
value <- c(21.25,12.82,21.25,12.61)
data <- data.frame(Type,condition,value)
ggplot(data, aes(fill=condition, y=value, x=Type)) +   geom_bar(position="stack", stat="identity")+theme_classic()+ theme(plot.title = element_blank(), axis.title.y =element_text(size=13,face="bold"), axis.title.x = element_blank(),axis.text.x = element_text(size=12,face="bold"), axis.text.y = element_text(size=12,hjust=0.5,face="bold"),legend.title=element_blank(),legend.text=element_text(size=12))+ylab("Proportion of genome with copy alternation (%)")+xlab("Cluster")+geom_text(aes(x = Type, y = value, label = value, group = condition),  position = position_stack(vjust = .5), size=4,color="white")
