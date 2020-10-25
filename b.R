##########experiment####
dir <- "D:/项目/base2b/20170930"
samples <- read.table(file.path(dir, "sample.list"), header = FALSE)
samples <- data.frame(samples[-c(1:8,33),])
files <- file.path(dir,samples$samples..c.1.8..33...., "quant.sf")
names(files) <- unlist(samples)
all(file.exists(files))
tx2gene <- read.csv("D:/项目/base2b/20170930/gencode.v34.annotation.tx2gene.CSV")
library(tximport)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
lipidGene=read.csv("D:/项目/脂代谢/GENE.csv", header = FALSE)
gdict <- read.csv("D:/项目/base2b/20170930/gdict.CSV", header = FALSE)
library(DESeq2)
colData <- data.frame(condition = factor(rep(c("DL","PM2.5","control"), each = 8)))
rownames(colData) <- colnames(txi$counts)
dds <- DESeqDataSetFromTximport(txi, colData, ~condition)
counts=data.frame(counts(dds))
counts$ensgene=rownames(counts)
counts=merge(counts,gdict,by.x="ensgene",by.y="V1")
base2b <- aggregate(x = counts[,2:25], by = list(counts$V2),FUN = max)
lipidB=base2b[base2b$Group.1%in%lipidGene$V1,]
rownames(lipidB)=lipidB[,1]
lipidB=lipidB[,-1]
lipidB=lipidB[rowSums(lipidB)>1,]
dds <-DESeqDataSetFromMatrix(lipidB,colData,design=~condition)
dds <- DESeq(dds)
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
dat.pca <- PCA(t(assay(vsd)), graph = FALSE) 
fviz_pca_ind(dat.pca,
             geom.ind = "point", 
             col.ind = colData$condition, 
             addEllipses = TRUE, 
             legend.title = "Groups")
lipVst <- as.data.frame(assay(vsd))
par(cex=0.7)
par(mar=c(1,3,1,1))
n.sample=ncol(lipVst)
cols=rainbow(n.sample*1.2)
par(mfrow=c(2,1))
colnames(lipVst)=factor(rep(c("DL","PM2.5","control"), each = 8))
boxplot(lipVst,col=cols,las=2)
resDL <- results(dds, contrast=c("condition",'DL','control'), tidy=TRUE)
resDL=na.omit(resDL)
DegDL=resDL[resDL$padj<0.05&abs(resDL$log2FoldChange)>log2(1),]
resPM <- results(dds, contrast=c("condition",'PM2.5','control'), tidy=TRUE)
resPM=na.omit(resPM)
DegPM=resPM[resPM$padj<0.05&abs(resPM$log2FoldChange)>log2(1),]
write.csv(DegDL,"D:/项目/base2b/DegDL.csv")
write.csv(DegPM,"D:/项目/base2b/DegPM.csv")
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(dplyr)
resDL$threshold[resDL$padj<0.05&resDL$log2FoldChange>0]="up"
resDL$threshold[resDL$padj<0.05&resDL$log2FoldChange<(-0)]="down"
resDL$threshold[(resDL$padj>0.05)|(resDL$log2FoldChange<=0)&resDL$log2FoldChange>=(-0)]="normal"
resDL$threshold=factor(resDL$threshold,levels = c("up","down","normal"),ordered = T)
p <- ggplot(data = resDL, 
            aes(x = log2FoldChange, 
                y = -log10(padj))) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=threshold)) +
  scale_color_manual(values=c("red","blue", "grey"))+
  geom_vline(xintercept=c(-0,0),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  theme_bw()
p
library(pheatmap)
heatdata <- lipVst[DegDL$row,c(1:8,17:24)]
annotation_col = data.frame(Type = factor(rep(c("DL","contral"), c(8,8))))
rownames(annotation_col) = colnames(heatdata)
pheatmap(heatdata, 
         cluster_rows = T,
         cluster_cols = T,
         annotation_col = annotation_col,
         color =colorRampPalette(c("blue", "white","red"))(100),
         show_rownames = F,
         show_colnames = F,
         scale = "row")

resPM$threshold[resPM$padj<0.05&resPM$log2FoldChange>0]="up"
resPM$threshold[resPM$padj<0.05&resPM$log2FoldChange<(-0)]="down"
resPM$threshold[(resPM$padj>0.05)|(resPM$log2FoldChange<=0)&resPM$log2FoldChange>=(-0)]="normal"
resPM$threshold=factor(resPM$threshold,levels = c("up","down","normal"),ordered = T)
p <- ggplot(data = resPM, 
            aes(x = log2FoldChange, 
                y = -log10(padj))) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=threshold)) +
  scale_color_manual(values=c("red","blue", "grey"))+
  geom_vline(xintercept=c(-0,0),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  theme_bw()
p

heatdata <- lipVst[DegPM$row,c(9:24)]
annotation_col = data.frame(Type = factor(rep(c("PM2.5","control"), c(8,8))))
rownames(annotation_col) = colnames(heatdata)
pheatmap(heatdata, 
         cluster_rows = T,
         cluster_cols = T,
         annotation_col = annotation_col,
         color =colorRampPalette(c("blue", "white","red"))(100),
         show_rownames = F,
         show_colnames = F,
         scale = "row")

#################### LUAD########
library(SummarizedExperiment)
library(TCGAbiolinks)
library(limma)
query.exp.hg38 <- GDCquery(project = "TCGA-LUAD", 
                           data.category = "Transcriptome Profiling", 
                           data.type = "Gene Expression Quantification", 
                           workflow.type = "HTSeq - Counts")
LUADRnaseqSE <- GDCdownload(query.exp.hg38)
LUADRnaseqSE <- GDCprepare(query.exp.hg38)
rownames(LUADRnaseqSE) <- values(LUADRnaseqSE)$external_gene_name
exp.hg38.values <- assay(LUADRnaseqSE)
Exp = exp.hg38.values[,substring(colnames(exp.hg38.values),14,16)=="01A"|
                        substring(colnames(exp.hg38.values),14,16)=="11A"]
lipidExp=Exp[rownames(Exp)%in%lipidGene$V1,]
ExpCan = lipidExp[,substring(colnames(lipidExp),14,16)=="01A"]
ExpNor = lipidExp[,substring(colnames(lipidExp),14,16)=="11A"]
lipidExp=as.data.frame(t(lipidExp))
lipidExp$patients=substring(rownames(lipidExp),1,16)
lipidExp=aggregate(x = lipidExp, by = list(lipidExp$patients),FUN = max)
rownames(lipidExp)=lipidExp[,1]
lipidExp=lipidExp[,-c(1,1040)]
lipidExp=as.data.frame(t(lipidExp))
library(DESeq2)
condition <- factor(ifelse(substring(colnames(lipidExp),14,16)=="01A","cancer","normal"))
colData <- data.frame(row.names =colnames(lipidExp),condition)
dds <-DESeqDataSetFromMatrix(lipidExp,colData,design=~condition)
dds <- DESeq(dds)
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
plotPCA(vsd, "condition")
dat.pca <- PCA(t(assay(vsd)), graph = FALSE) 
fviz_pca_ind(dat.pca,
             geom.ind = "point", 
             col.ind = colData$condition, 
             addEllipses = TRUE, 
             legend.title = "Groups")

mRNAVst <- as.data.frame(assay(vsd))
par(cex=0.7)
par(mar=c(1,3,1,1))
n.sample=ncol(expVst)
cols=rainbow(n.sample*1.2)
par(mfrow=c(2,1))
boxplot(expVst,col=cols,las=2)
res <- results(dds, contrast=c("condition",'cancer','normal'),tidy=TRUE) 
res=na.omit(res)
Deg=res[res$padj<0.05&abs(res$log2FoldChange)>log2(1),]
write.csv(Deg,"D:/项目/base2b/DegLUAD.csv")
res$change <- ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > log2(1),
                        ifelse(res$log2FoldChange > log2(1) ,"up","down"),"normal") 
p <- ggplot(data = res, 
            aes(x = log2FoldChange, 
                y = -log10(padj))) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=change)) +
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-log2(1),log2(1)),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  theme_bw()
p
heatdata <- mRNAVst[Deg$row,]
heatdata=heatdata[,order(substr(colnames(heatdata),14,16))]
condition=factor(ifelse(substring(colnames(heatdata),14,16)=="01A","cancer","normal"))
annotation_col <- data.frame(row.names =colnames(heatdata),condition)

pheatmap(heatdata, 
         cluster_rows = T,
         cluster_cols = F,
         annotation_col = annotation_col,
         color =colorRampPalette(c("blue", "white","red"))(50),
         show_rownames = F,
         show_colnames = F,scale = "row")

DegUp=intersect(DegPM[DegPM$log2FoldChange>0,]$row,Deg[Deg$log2FoldChange>0,]$row)
DegDown=intersect(DegPM[DegPM$log2FoldChange<0,]$row,Deg[Deg$log2FoldChange<0,]$row)
DegUpDown=intersect(DegUp,DegDL[DegDL$log2FoldChange<0,]$row)
DegDownUp=intersect(DegDown,DegDL[DegDL$log2FoldChange>0,]$row)
keyGene=union(DegUpDown,DegDownUp)
logFC=data.frame(LUAD=Deg[Deg$row %in% keyGene,][order(Deg[Deg$row %in% keyGene,]$row),]$log2FoldChange,
                 PM=DegPM[DegPM$row %in% keyGene,][order(DegPM[DegPM$row %in% keyGene,]$row),]$log2FoldChange,
                 DL=DegDL[DegDL$row %in% keyGene,][order(DegDL[DegDL$row %in% keyGene,]$row),]$log2FoldChange) 
rownames(logFC)=Deg[Deg$row %in% keyGene,][order(Deg[Deg$row %in% keyGene,]$row),]$row
write.csv(logFC,"D:/项目/base2b/DEGs/logFC.csv")
###########临床#######
pd.all <- read.delim("TCGA-LUAD.GDC_phenotype.tsv/TCGA-LUAD.GDC_phenotype.tsv", header = T, stringsAsFactors = F)
samples=cli[,c("bcr_patient_barcode","gender","vital_status",
               "days_to_death","days_to_last_followup","race_list",
               "age_at_initial_pathologic_diagnosis",
               "tobacco_smoking_history","stage_event_pathologic_stage",
               "stage_event_tnm_categories")]
samples$survival_time<-ifelse(samples$vital_status=='Alive',
                              samples$days_to_last_followup,samples$days_to_death)
samples$status=ifelse(samples$vital_status=="Dead",1,0)
samples=samples[!duplicated(samples),]
samples=samples[!(samples$race_list==""),]
samples=samples[!(is.na(samples$tobacco_smoking_history)),]
samples=samples[!(samples$tobacco_smoking_history==5),]

samples1<-samples[match(substr(colnames(ExpCan),1,12),samples$bcr_patient_barcode),]
all(substr(rownames(samples1),1,12)==substr(colnames(ExpCan),1,12))
rownames(samples1)<-colnames(ExpCan)
samples1=samples1[na.omit(samples1$bcr_patient_barcode),]
samples2<-samples[match(substr(colnames(ExpNor),1,12),samples$bcr_patient_barcode),]
all(substr(rownames(samples2),1,12)==substr(colnames(ExpNor),1,12))
rownames(samples2)<-colnames(ExpNor)
samples2=samples2[na.omit(samples2$bcr_patient_barcode),]
ExpCan<-ExpCan[,rownames(samples1)]
ExpNor<-ExpNor[,rownames(samples2)]
mExp=cbind(ExpCan,ExpNor)
#########################clusterProfiler######
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(GOplot)
library(enrichplot)
gene=Deg[Deg$row%in%keyGene,c(1,3)]
names(gene)=c("ID","logFC")
gene$genelist=mapIds(org.Hs.eg.db,keyGene,keytype="SYMBOL",column="ENTREZID")
go=enrichGO(gene=gene$genelist,OrgDb = "org.Hs.eg.db",keyType = "ENTREZID",pvalueCutoff = 0.05,readable = T,ont="ALL")
dotplot(go, showCategory=10)
write.csv(go,"D:/项目/base2b/DEGs/GoKeyGene.csv")
kegg=enrichKEGG(genelist,organism = "hsa",pvalueCutoff = 0.05)
kegg=setReadable(kegg,org.Hs.eg.db,keyType = "ENTREZID")
dotplot(kegg, showCategory=10)
write.csv(kegg,"D:/项目/base2b/DEGs/KeggKeyGene.csv")
plotGOgraph(go)
cnetplot(go, showCategory = 10)
go1=data.frame(go)
library(stringr)
go2=go1[,c(1:3,8:9)]
go2$geneID=str_replace_all(go2$geneID,"/",",")
names(go2)=c("Category","ID","Term","adj_pval","Genes")
gene=Deg[Deg$row%in%keyGene,c(1,3)]
names(gene)=c("ID","logFC")
plot_data=list(DA=go2,GE=gene)
circ2=data.frame()
circ2=circle_dat(plot_data$DA,plot_data$GE)
reduced_circ <- reduce_overlap(circ2, overlap = 0.75)
GOBubble(reduced_circ, labels = 1)
kegg1=data.frame(kegg)
kegg1=kegg1[,c(1,2,7,8)]
kegg1$geneID=str_replace_all(kegg1$geneID,"/",",")
names(kegg1)=c("ID","Term","adj_pval","Genes")
kegg1$Category=c("KEGG pathway")
plot_data=list(DA=kegg1,GE=gene)
circ2=data.frame()
circ2=circle_dat(plot_data$DA,plot_data$GE)
GOBubble(circ2,labels = 0)

heatdata <- logFC
pheatmap(heatdata, 
         cluster_rows = T,
         cluster_cols = T,
         annotation_legend  = T, 
         show_rownames = T,
         show_colnames = T,
         display_numbers=T,
         number_color="black",
         fontsize=12,
         color = colorRampPalette(c("blue", "white","red"))(50))


###############GSE93329######
library(GEOquery)
gset <- getGEO("GSE93329", GSEMatrix =TRUE, AnnotGPL=TRUE )
exprSet <- exprs(gset[[1]])
pData <- pData(gset[[1]])
fdata<-fData(gset[[1]])
library(limma)
sign =c()
gene = fdata$gene_assignment
for (i in 1:length(fdata$gene_assignment)){
  tmp = (strsplit2(as.character(gene[i]),split = "//"))
  fdata$sign[i] = tmp[2]
}
exprSet=data.frame(exprSet)
exprSet$ID=rownames(exprSet)
exprSet =merge(exprSet,fdata[,c(1,16)],by.x="ID",by.y="ID")
exprSet<-exprSet[!(exprSet$sign==" --- "), ]
exprSet <- aggregate(x = exprSet[,2:7], by = list(exprSet$sign),FUN = max)
rownames(exprSet) =exprSet$Group.1
exprSet=exprSet[,-1]
condition =factor(rep(c("control","PM2.5"), each = 3))
colnames(exprSet) = paste(condition,1:6,sep = '_')
sampleDists <- dist(t(data.frame(exprSet)))
library(factoextra)
res <- hcut(sampleDists, k = 2, stand = TRUE)
fviz_dend(res,rect = TRUE,rect_border="cluster",rect_lty=2,lwd=1.2,
          rect_fill = T,cex = 2,color_labels_by_k=T,horiz=T)
exprSet1=exprSet[,-3]
sampleDists <- dist(t(data.frame(exprSet1)))
res <- hcut(sampleDists, k = 2, stand = TRUE)
fviz_dend(res,rect = TRUE,rect_border="cluster",rect_lty=2,lwd=1.2,
          rect_fill = T,cex = 2,color_labels_by_k=T,horiz=T)

par(cex=0.7)
par(mar=c(1,3,1,1))
n.sample=ncol(exprSet1)
cols=rainbow(n.sample*1.2)
par(mfrow=c(2,1))
boxplot(exprSet1,col=cols,las=2)
library(ggfortify)
df = as.data.frame(t(exprSet1)) 
group = factor(rep(c("control","PM2.5"),c(2,3)))
df$group = group
p5 = autoplot(prcomp(df[,1:(ncol(df)-1)]),
              data = df,colour = 'group',label=F,frame=T,frame.type='norm') 
design=model.matrix(~0+group)
colnames(design)=levels(group)
rownames(design)=colnames(exprSet1)
fit=lmFit(exprSet1,design)
cont.matrix<-makeContrasts(PM2.5-control,levels = design)
fit2=contrasts.fit(fit,cont.matrix)
fit2 <- eBayes(fit2)
tempOutput = topTable(fit2,coef=1,n=Inf,adjust="BH") 
nrDEG = na.omit(tempOutput)
deg1=nrDEG[c(" HSPH1 "," LDHA "," HSP90AA1 "," SMS "," SOAT1 "," ODC1 "," PTDSS1 ",
        " CCNC "," LBR "," ADSL "," FH "," RAN "," VDAC2 "," SREBF2 "," SLC27A3 ",
        " GLUL "," PNPLA2 "," SLC44A2 "," ESYT1 "),]
write.csv(deg1,"D:/项目/base2b/GSE93329/deg1.csv")

#########
library(GEOquery)
gset <- getGEO("GSE139294", GSEMatrix =TRUE, AnnotGPL=TRUE )
exprSet <- exprs(gset[[1]])
pData <- pData(gset[[1]])
fdata<-fData(gset[[1]])
################SNP#####
library(maftools)
LUAD_maf=GDCquery_Maf("LUAD",pipelines = "varscan2")
LUAD_maf=read.maf(maf=LUAD_maf)
write.mafSummary(maf=LUAD_maf,basename = "luad")
plotmafSummary(maf = LUAD_maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
oncoplot(maf = LUAD_maf, top = 20)
oncostrip(maf = LUAD_maf, genes = c('SLC27A3','SLC44A2', 'GLUL',"PNPLA2","RAN","LDHA","PTDSS1"))
############miRNA##########
query <- GDCquery(project = "TCGA-LUAD", 
                  data.category = "Transcriptome Profiling", 
                  data.type = "Isoform Expression Quantification")
GDCdownload(query)
mir_exp<- GDCprepare(query = query,
                     summarizedExperiment=F)
sample_id <- as.data.frame(table(mir_exp$barcode))
miRNA_id <- as.data.frame(table(mir_exp$miRNA_region))
miRNA_count <- matrix(NA,ncol = nrow(sample_id),nrow = nrow(miRNA_id))
colnames(miRNA_count) <- sample_id$Var1
rownames(miRNA_count) <- as.character(miRNA_id$Var1)
for(i in 1:nrow(sample_id)){
  temp1 <- mir_exp[which(mir_exp$barcode==as.character(sample_id[i,1])),]
  for(j in 1:nrow(miRNA_id)){
    loc <- which(temp1$miRNA_region==as.character(miRNA_id[j,1]))
    if(length(loc)>0){
      miRNA_count[j,i] <- sum(temp1[loc,3])
    }else{
      miRNA_count[j,i] <- 0
    }
  }
  print(i)
}
miRNA_count <- miRNA_count[1:2213,]
name <- substr(rownames(miRNA_count),8,nchar(rownames(miRNA_count)))
library(miRBaseVersions.db)
items <- select(miRBaseVersions.db,
                keys = name,
                keytype = "MIMAT",
                columns = c('ACCESSION',"NAME",'VERSION'))
id_name <- items[items$VERSION == 21.0, c("ACCESSION","NAME")]
miRNA_count <- cbind(id_name,miRNA_count)
rownames(miRNA_count)=miRNA_count[,2]
miRNA_count=miRNA_count[,-c(1,2)]
miExp=miRNA_count[,substring(colnames(miRNA_count),14,16)=="01A"|substring(colnames(miRNA_count),14,16)=="11A"]
colnames(miExp)=substring(colnames(miExp),1,16)
miExp=aggregate(t(miExp),by=list(colnames(miExp)),FUN=max)
rownames(miExp)=miExp[,1]
miExp=as.data.frame(t(miExp[,-1]))
miCan=miExp[,substring(colnames(miExp),14,16)=="01A"]
miNor=miExp[,substring(colnames(miExp),14,16)=="11A"]
samples1<-samples[match(substr(colnames(miCan),1,12),samples$bcr_patient_barcode),]
all(substr(rownames(samples1),1,12)==substr(colnames(miCan),1,12))
rownames(samples1)<-colnames(miCan)
samples1=samples1[na.omit(samples1$bcr_patient_barcode),]
samples2<-samples[match(substr(colnames(miNor),1,12),samples$bcr_patient_barcode),]
all(substr(rownames(samples2),1,12)==substr(colnames(miNor),1,12))
rownames(samples2)<-colnames(miNor)
samples2=samples2[na.omit(samples2$bcr_patient_barcode),]
miCan<-miCan[,rownames(samples1)]
miNor<-miNor[,rownames(samples2)]
miExp=cbind(miCan,miNor)
library(DESeq2)
condition <- factor(ifelse(substring(colnames(miExp),14,16)=="01A","cancer","normal"))
colData <- data.frame(row.names =colnames(miExp),condition)
dds <-DESeqDataSetFromMatrix(miExp,colData,design=~condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
plotPCA(vsd, "condition")
miRNAVst <- as.data.frame(assay(vsd))
res <- results(dds, contrast=c("condition",'cancer','normal'),tidy=TRUE) 
res=na.omit(res)
Deg=res[res$padj<0.05&abs(res$log2FoldChange)>log2(1),]
write.csv(Deg,"D:/项目/base2b/Deg_miRNA.csv")
hsa.miR.665=as.data.frame(t(miRNAVst["hsa-miR-665",]))
hsa.miR.665$group=ifelse(substring(colnames(miRNAVst),14,16)=="01A","cancer","normal")
colnames(hsa.miR.665)[1]="hsa-miR-665"
library(ggstatsplot)
ggbetweenstats(data=hsa.miR.665, 
               x = group, 
               y = "hsa-miR-665") 


sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)
pheatmap(assay(vsd)[select,], cluster_rows=T, show_rownames=FALSE,show_colnames=FALSE,
         cluster_cols=T, annotation_col=colData)
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)



##########
mRNAexp=as.data.frame(t(mRNAVst[c("SLC27A3","SLC44A2","GLUL","PNPLA2","RAN","LDHA","PTDSS1"),]))
mRNAexp=mRNAexp[substring(colnames(miRNAVst),1,16),]
mRNAexp=mRNAexp[substring(rownames(mRNAexp),14,16)=="01A",]
miRNAexp=as.data.frame(t(miRNAVst[Deg$row,]))
rownames(miRNAexp)=substring(rownames(miRNAexp),1,16)
miRNAexp=miRNAexp[substring(rownames(mRNAexp),1,16),]

library("psych")
cor_matrix <- corr.test(log2(mRNAexp+1),log2(miRNAexp+1))
str(cor_matrix)
library(reshape2)
library(knitr)
c.r <- melt(cor_matrix$r, value.name = "cor", as.is = T)
p_value <- melt(cor_matrix$p)[, 3]
out <- data.frame(c.r, p_value, check.names = F, stringsAsFactors = F)
out1=out[out$cor<0&out$p_value<0.05,]
write.csv(out1,"D:/项目/base2b/cor_miRNA.csv")
#############
library(ggalluvial)
network=read.csv("D:/项目/base2b/miRNA/network.csv")
p=ggplot(network,aes(axis1=miRNA,axis2=mRNA,y=Freq))+
  geom_alluvium(aes(fill=mRNA),width = 0.1,knot.pos = 0.1)+
  geom_stratum(fill="white",color="grey20",alpha=0.7,width = 1/7)+
  geom_text(stat="stratum",size=2.5,color="black",aes(label = after_stat(stratum)))+
  scale_x_discrete(limits=c("miRNA","mRNA"),expand=c(0,0))+
  xlab("")+ylab("")+
  guides(fill=F)
p
###############
setwd("D:/项目/甲基化")
library("ChAMP")
library("readr")
m6a=readr::read_tsv("TCGA-LUAD.methylation450.tsv/TCGA-LUAD.methylation450.tsv")
m6a=as.data.frame(m6a)
rownames(m6a)=m6a[,1]
m6a=m6a[,-1]
beta=as.matrix(m6a)
require(GEOquery)
require(Biobase)
GSE40279 <- getGEO("GSE40279")
beta.m <- exprs(GSE40279[[1]])
save(beta.m,file = 'Rdata/GSE40279.Rdata')
save(beta,file = 'Rdata/TCGA.Rdata')
rm(list = ls())
load('甲基化/Rdata/GSE40279.Rdata')
load('甲基化/Rdata/TCGA.Rdata')
library("impute")
beta=impute.knn(beta) 
betaData=beta$data
betaData=betaData+0.00001
sum(is.na(betaData))
betaData=betaData[,substring(colnames(betaData),14,16)=="01A"]
betaData=as.data.frame(betaData)
betaData$symbol=rownames(betaData)
beta.m=as.data.frame(beta.m)
beta.m$symbol=rownames(beta.m)
beta=merge(betaData,beta.m,by="symbol")
rownames(beta)=beta[,1]
beta=beta[,-1]
pd=data.frame(sampleID=colnames(beta),sample_type=rep(c("Tumor","Normal"),c(455,656)))
rownames(pd)=pd$sampleID

myLoad=champ.filter(beta = as.matrix(beta) ,pd = pd)
dim(myLoad$beta)
save(myLoad,file = 'Rdata/step1-output.Rdata')
rm(list = ls())
load('Rdata/step1-output.Rdata')
myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores = 100)
dim(myNorm)
QC.GUI(beta=myNorm,arraytype="450K")
num.na <- apply(myNorm,2,function(x)(sum(is.na(x))))
hist(num.na)
table(num.na)
myNorm <- myNorm[,which(num.na < 250000)]
myNorm=myNorm[,!colnames(myNorm)%in%("TCGA-64-5774-01A")]
save(myNorm,file = '甲基化/step2-champ_myNorm.Rdata')
###########
library("FactoMineR")
library("factoextra") 
dat <- t(myNorm)
group_list=ifelse(substring(colnames(myNorm),14,15)=="01","Tumor","Normal")
table(group_list)
dat.pca <- PCA(dat, graph = FALSE) 
fviz_pca_ind(dat.pca,
             geom.ind = "point", 
             col.ind = group_list, 
             addEllipses = TRUE, 
             legend.title = "Groups")
cg=names(tail(sort(apply(myNorm,1,sd)),1000))
library(pheatmap)
ac=data.frame(group=group_list)
rownames(ac)=colnames(myNorm)  
pheatmap(myNorm[cg,],show_colnames =F,show_rownames = F,
         annotation_col=ac)
pheatmap::pheatmap(cor(myNorm[cg,]),
                   annotation_col = ac,
                   show_rownames = F,
                   show_colnames = F)
#######
rm(list=ls())
load("Rdata/step2-champ_myNorm.Rdata")
group_list <- ifelse(substring(colnames(myNorm),14,15)=="01","Tumor","Normal")
myDMP <- champ.DMP(beta = myNorm,pheno=group_list)
head(myDMP[[1]])
save(myDMP,file = '甲基化/step3-output-myDMP.Rdata')
rm(list=ls())
load("Rdata/step3-output-myDMP.Rdata")
df_DMP <- myDMP$Tumor_to_Normal
lipiddf_DMP=df_DMP[df_DMP$gene%in%lipidGene$V1,]
logFC_cutoff <- 0.2
lipiddf_DMP$change[lipiddf_DMP$adj.P.Val<0.05&lipiddf_DMP$logFC>logFC_cutoff]="up"
lipiddf_DMP$change[lipiddf_DMP$adj.P.Val<0.05&lipiddf_DMP$logFC<(-logFC_cutoff)]="down"
lipiddf_DMP$change[(lipiddf_DMP$adj.P.Val>0.05)|(lipiddf_DMP$logFC<=logFC_cutoff)&lipiddf_DMP$logFC>=(-logFC_cutoff)]="normal"
lipiddf_DMP$change=factor(lipiddf_DMP$change,levels = c("up","down","normal"),ordered = T)
table(lipiddf_DMP$change) 
library(ggplot2)
this_tile <- paste0('Cutoff for logFC is ',round(0.2,3),
                    '\nThe number of hypermethylated DMPs is ',nrow(lipiddf_DMP[lipiddf_DMP$change =='up',]) ,
                    '\nThe number of hypomethylated DMPs is ',nrow(lipiddf_DMP[lipiddf_DMP$change =='down',]))

ggplot(data = lipiddf_DMP, 
       aes(x = logFC, 
           y = -log10(adj.P.Val))) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=change)) +
  scale_color_manual(values=c('red','blue','grey'))+
  geom_vline(xintercept=c(-0.2,0.2),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  ggtitle( this_tile )+
  theme_bw()
choose_gene <- data.frame(rownames(lipiddf_DMP[lipiddf_DMP$change != "normal",]))
choose_gene=lipiddf_DMP[lipiddf_DMP$change != "normal",c(14,15)]
choose_matrix <- myNorm[rownames(choose_gene),]
DMPs=lipiddf_DMP[rownames(choose_gene),]
write.csv(DMPs,"DMPs.csv")
annotation_col <- data.frame(Sample=c(ifelse(substring(colnames(choose_matrix),14,15)=="01","cancer","normal"))) 
rownames(annotation_col) <- colnames(choose_matrix)
ann_colors = list(Sample = c(cancer="green", normal="red"))
library(pheatmap)
pheatmap(choose_matrix,show_colnames = F,show_rownames = F,annotation_col = annotation_col,
         border_color=NA,
         color = colorRampPalette(colors = c("white","navy"))(50),
         annotation_colors = ann_colors)
pheatmap(choose_matrix, 
         cluster_rows = T,
         cluster_cols = T,
         annotation_col = annotation_col,
         color =colorRampPalette(c("blue", "white","red"))(100),
         show_rownames = F,
         show_colnames = F,
         scale = "row")

DMP.GUI(DMP = myDMP[[1]],beta = myNorm,pheno = group_list)
#############
library(MethylMix)
library(doParallel)
logFC=read.csv("D:/项目/base2b/DEGs/logFC.csv")
MET=choose_gene[choose_gene$gene%in%logFC$X,]
METcancer=choose_matrix[rownames(MET),substr(colnames(choose_matrix),14,15)=="01"]
METnormal=choose_matrix[rownames(MET),substr(colnames(choose_matrix),14,15)=="11"]
GEcancer=mRNAVst[logFC$X,substr(colnames(mRNAVst),14,15)=="01"]
library(ggpubr)
METcancer=as.matrix(METcancer[,intersect(colnames(GEcancer),colnames(METcancer))])
GEcancer<-as.matrix(GEcancer[,intersect(colnames(GEcancer),colnames(METcancer))])
library(MethylMix)
res <- ClusterProbes(METcancer, METnormal)
MethylMix(METcancer, GEcancer, METnormal = NULL, listOfGenes = NULL, filter = TRUE, NoNormalMode = FALSE, OutputRoot = "")
library("psych")
cor_matrix <- corr.test(t(METcancer),t(GEcancer))
str(cor_matrix)
library(reshape2)
library(knitr)
c.r <- melt(cor_matrix$r, value.name = "cor", as.is = T)
p_value <- melt(cor_matrix$p)[, 3]
out <- data.frame(c.r, p_value, check.names = F, stringsAsFactors = F)
out1=out[out$cor<(-0.2)&out$p_value<0.05,]
MET$Var1=rownames(MET)
library(dplyr)
out1=left_join(out1,MET)
out1=out1[out1$Var2==out1$gene,]
write.csv(out1,"D:/项目/base2b/cor_Methyl.csv")
################
names(pd.all)
cli=pd.all[,c("submitter_id.samples","gender.demographic","vital_status.demographic",
           "days_to_death.demographic","days_to_last_follow_up.diagnoses","race.demographic",
           "age_at_initial_pathologic_diagnosis",
           "tobacco_smoking_history","tumor_stage.diagnoses")]
cli=cli[cli$submitter_id.samples%in%colnames(choose_matrix),]
library(glmnet)
x=as.matrix(t(choose_matrix[out1$Var1,]))
y=ifelse(substring(colnames(x),14,15)=="01",1,0)
set.seed(123456)
fit3 <- cv.glmnet(x, y,family="binomial",type.measure = "deviance")
print(model_lasso)