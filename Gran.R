dir <- "E:/项目/base2b/20170930"
samples <- read.table(file.path(dir, "sample.list"), header = FALSE)
samples <- data.frame(samples[-33,])
files <- file.path(dir,samples$samples..33..., "quant.sf")
names(files) <- unlist(samples)
all(file.exists(files))
tx2gene <- read.csv("E:/项目/base2b/20170930/gencode.v34.annotation.tx2gene.CSV")
library(tximport)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
lipidGene=read.csv("E:/项目/脂代谢/GENE.csv", header = FALSE)
gdict <- read.csv("E:/项目/base2b/20170930/gdict.CSV", header = FALSE)
library(DESeq2)
colData <- data.frame(condition = factor(rep(c("DL_PM2.5","DL","PM2.5","control"), each = 8)))
rownames(colData) <- colnames(txi$counts)
dds <- DESeqDataSetFromTximport(txi, colData, ~condition)
counts=data.frame(counts(dds))
counts$ensgene=rownames(counts)
counts=merge(counts,gdict,by.x="ensgene",by.y="V1")
base2b <- aggregate(x = counts[,2:33], by = list(counts$V2),FUN = max)
lipidB=base2b[base2b$Group.1%in%lipidGene$V1,]
rownames(lipidB)=lipidB[,1]
lipidB=lipidB[,-1]
lipidB=lipidB[rowSums(lipidB)>1,]
dds <-DESeqDataSetFromMatrix(lipidB,colData,design=~condition)
dds <- DESeq(dds)
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
plotPCA(vsd, "condition")
lipVst <- as.data.frame(assay(vsd))
par(cex=0.7)
par(mar=c(1,3,1,1))
n.sample=ncol(lipVst)
cols=rainbow(n.sample*1.2)
par(mfrow=c(2,1))
boxplot(lipVst,col=cols,las=2)
resDL <- results(dds, contrast=c("condition",'DL','control'), tidy=TRUE)
resDL=na.omit(resDL)
DegDL=resDL[resDL$padj<0.05&abs(resDL$log2FoldChange)>log2(1),]
write.csv(DegDL,"E:/项目/国自然/DegDL.csv")
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
heatdata <- lipVst[DegDL$row,c(9:16,25:32)]
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
#############
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
Exp = as.data.frame(exp.hg38.values[,substring(colnames(exp.hg38.values),14,16)=="01A"|
                        substring(colnames(exp.hg38.values),14,16)=="11A"])
lipidExp=Exp[rownames(Exp)%in%lipidGene$V1,]
ExpCan = lipidExp[,substring(colnames(lipidExp),14,16)=="01A"]
ExpNor = lipidExp[,substring(colnames(lipidExp),14,16)=="11A"]
cli.query <- GDCquery(project =c("TCGA-LUAD"),file.type = "xml",
                      data.category = "Clinical") 
GDCdownload(cli.query)
cli <- GDCprepare_clinic(cli.query, clinical.info = "patient")
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
samples=rbind(samples1,samples2)
samples$group=ifelse(substring(rownames(samples),14,16)=="01A","cancer","normal")
library(DESeq2)
condition <- factor(ifelse(substring(colnames(mExp),14,16)=="01A","cancer","normal"))
colData <- data.frame(row.names =colnames(mExp),condition)
dds <-DESeqDataSetFromMatrix(mExp,colData,design=~condition)
dds <- DESeq(dds)
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
library(tableone)
stable <- CreateTableOne(vars=c("gender","race_list","age_at_initial_pathologic_diagnosis","tobacco_smoking_history"),
                         strata="group",data=samples,
                         factorVars=c("gender","race_list","tobacco_smoking_history"))
stable=print(stable,showAllLevels = TRUE)
write.csv(stable,"E:/项目/国自然/table1.csv")
plotPCA(vsd, "condition")
mRNAVst <- as.data.frame(assay(vsd))
par(cex=0.7)
par(mar=c(1,3,1,1))
n.sample=ncol(expVst)
cols=rainbow(n.sample*1.2)
par(mfrow=c(2,1))
boxplot(expVst,col=cols,las=2)
res <- results(dds, contrast=c("condition",'cancer','normal'),tidy=TRUE) 
res=na.omit(res)
Deg=res[res$padj<0.05&abs(res$log2FoldChange)>log2(2),]
write.csv(Deg,"E:/项目/国自然/DegLUAD.csv")
res$threshold[res$padj<0.05&res$log2FoldChange>1]="up"
res$threshold[res$padj<0.05&res$log2FoldChange<(-1)]="down"
res$threshold[(res$padj>0.05)|(res$log2FoldChange<=1)&res$log2FoldChange>=(-1)]="normal"
res$threshold=factor(res$threshold,levels = c("up","down","normal"),ordered = T)
p <- ggplot(data = res, 
            aes(x = log2FoldChange, 
                y = -log10(padj))) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=threshold)) +
  scale_color_manual(values=c("red","blue", "grey"))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  theme_bw()
p
heatdata <- mRNAVst[Deg$row,]
annotation_col = data.frame(Type = factor(rep(c("cancer","normal"), c(417,52))))
rownames(annotation_col) = colnames(heatdata)
pheatmap(scale(heatdata), 
         cluster_rows = T,
         cluster_cols = F,
         annotation_col = annotation_col,
         color =colorRampPalette(c("blue", "white","red"))(100),
         show_rownames = F,
         show_colnames = F,scale = "row")

DegUp=intersect(DegDL[DegDL$log2FoldChange>0,]$row,Deg[Deg$log2FoldChange<0,]$row)
DegDown=intersect(DegDL[DegDL$log2FoldChange<0,]$row,Deg[Deg$log2FoldChange>0,]$row)
keyGene=union(DegUp,DegDown)
rownames(DegDL)=DegDL$row
rownames(Deg)=Deg$row
keyGene=cbind(DegDL[keyGene,]$log2FoldChange,Deg[keyGene,]$log2FoldChange)
rownames(keyGene)=union(DegUp,DegDown)
colnames(keyGene)=c("logFC_DL","logFC_TCGA")
write.csv(keyGene,"E:/项目/国自然/keyGene.csv")
heatdata <- mRNAVst[Deg$row,]
keyGene1<-log2(keyGene+1)
annotation_col = data.frame(Type = factor(rep(c("cancer","normal"), c(426,52))))
rownames(annotation_col) = colnames(heatdata1)
pheatmap(keyGene, 
         cluster_rows = T,
         cluster_cols = T,
         color =colorRampPalette(c("blue", "white","red"))(100),
         show_rownames = T,
         show_colnames = T)

################
