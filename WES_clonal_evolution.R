library(clonevol)
library(sciClone)
#Prepare input for sciclone
v1 = read.table("F0_sciclone.txt", header = T)
v2 = read.table("F10_sciclone.txt", header = T)
regions = read.table("cnv.txt")
regions=regions[,1:3]
names = c("B16F0","B16F10")
c1= read.table("F0_cna.txt", header = F)
c2= read.table("F10_cna.txt", header = F)
sc2 = sciClone(vafs=list(v1,v2), 
c1=c1[(abs(c1$V4)>0.25),]
c2=c2[(abs(c2$V4)>0.25),]
sc4 = sciClone(vafs=list(v1,v2), copyNumberCalls=list(c1,c2),sampleNames=names[1:2], regionsToExclude=NULL,minimumDepth=20, maximumClusters=10, useSexChrs=FALSE,cnCallsAreLog2=TRUE)
?sciClone
writeClusterTable(sc4, "Cluster_F0F104.txt")
sc.plot2d(sc4,"Cluster_F0F104.pdf")
#From sciclone clustering result, after removal of cluster with VAFs nearly 100%
x = read.table("Cluster_F0F105.txt", header = T)
vaf.col.names <- grep('.vaf', colnames(x), value=T)
sample.names <- gsub('.vaf', '', vaf.col.names)
x[, sample.names] <- x[, vaf.col.names]
vaf.col.names <- sample.names
sample.groups <- c("B16F0","B16F10")
names(sample.groups) <- vaf.col.names
x <- x[order(x$cluster),]
clone.colors <- clone.colors <- c('#999793', '#8d4891', '#f8e356',
                                  '#fe9536', '#d7352e', '#1c96b5')
plot.pairwise(x, col.names = vaf.col.names, out.prefix = 'F0F105.pdf',colors = clone.colors, shapes = c(1,2,3,4,5,6), xMinSmall = 0, xMaxSmall= 100, yMinSmall = 0, yMaxSmall = 100)
y = infer.clonal.models(variants = x, cluster.col.name = 'cluster',vaf.col.names = vaf.col.names, sample.groups = sample.groups,subclonal.test = 'bootstrap', subclonal.test.model = 'non-parametric', num.boots = 1000, founding.cluster = 1, cluster.center = 'mean',ignore.clusters = NULL, clone.colors = clone.colors, min.cluster.vaf =  0, sum.p = 0.05, alpha = 0.05)
y <- convert.consensus.tree.clone.to.branch(y, branch.scale = 'sqrt')
plot.clonal.models(y, box.plot = TRUE, fancy.boxplot = TRUE,fancy.variant.boxplot.jitter.alpha = 1, fancy.variant.boxplot.jitter.center.color = 'grey50', fancy.variant.boxplot.base_size = 12, fancy.variant.boxplot.plot.margin = 1, fancy.variant.boxplot.vaf.suffix = '.VAF', clone.shape = 'bell',bell.event = FALSE, clone.time.step.scale = 1, bell.curve.step = 2,merged.tree.plot = FALSE, merged.tree.clone.as.branch = TRUE,mtcab.event.sep.char = ',', mtcab.show.event = FALSE,mtcab.tree.rotation.angle = 0, mtcab.branch.text.size = 0.001,mtcab.branch.width = 0.8, mtcab.node.size = 3, mtcab.node.label.size = 0.75, mtcab.node.text.size = 0.001, cell.plot = TRUE, num.cells = 100, cell.border.size = 0.25, cell.border.color = 'black', clone.grouping ='horizontal', scale.monoclonal.cell.frac = TRUE, show.score = FALSE,cell.frac.ci = TRUE, disable.cell.frac = FALSE, out.dir = 'F0F10_5',out.format = 'pdf', overwrite.output = TRUE, width = 16, height = 6,panel.widths = c(4,5,2,4))
