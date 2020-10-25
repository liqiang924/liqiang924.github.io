###########vitro###########
dir <- "E:/项目/朱老师/vitro/quant.sf"
samples <- data.frame(samples=(paste(rep(c("AB","AL","AP","AS","HB","HL","HP","HS"),each=3),c(1:3),sep = "-")))
files <- file.path(dir,samples$samples, "quant.sf")
names(files) <- unlist(samples)
all(file.exists(files))
tx2gene <- read.csv("E:/项目/朱老师/vitro/uniqGeneName.CSV",header = F)
library(tximport)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
lipidGene=read.csv("E:/项目/脂代谢/GENE.csv", header = FALSE)
library(DESeq2)
colData <- data.frame(condition = factor(rep(c("AB","AL","AP","AS","HB","HL","HP","HS"), each = 3)))
rownames(colData) <- colnames(txi$counts)
dds <- DESeqDataSetFromTximport(txi, colData, ~condition)
dds <- DESeq(dds)
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
plotPCA(vsd, "condition")
expVst <- as.data.frame(assay(vsd))
library(factoextra)
sampleDists <- dist(t(data.frame(expVst[,c(4:6,10:12)])))
res <- hcut(sampleDists, k = 2, stand = TRUE)
fviz_dend(res,rect = TRUE,rect_border="cluster",rect_lty=2,lwd=1.0,
          rect_fill = T,cex = 1,color_labels_by_k=T,horiz=T)
sampleDists <- dist(t(data.frame(expVst[,c(7:9,10:12)])))
res <- hcut(sampleDists, k = 2, stand = TRUE)
fviz_dend(res,rect = TRUE,rect_border="cluster",rect_lty=2,lwd=1.0,
          rect_fill = T,cex = 1,color_labels_by_k=T,horiz=T)
sampleDists <- dist(t(data.frame(expVst[,c(16:18,22:24)])))
res <- hcut(sampleDists, k = 2, stand = TRUE)
fviz_dend(res,rect = TRUE,rect_border="cluster",rect_lty=2,lwd=1.0,
          rect_fill = T,cex = 1,color_labels_by_k=T,horiz=T)
sampleDists <- dist(t(data.frame(expVst[,c(19:24)])))
res <- hcut(sampleDists, k = 2, stand = TRUE)
fviz_dend(res,rect = TRUE,rect_border="cluster",rect_lty=2,lwd=1.0,
          rect_fill = T,cex = 1,color_labels_by_k=T,horiz=T)
par(cex=0.7)
par(mar=c(1,3,1,1))
n.sample=ncol(expVst)
cols=rainbow(n.sample*1.2)
par(mfrow=c(2,1))
boxplot(expVst,col=cols,las=2)
resADL <- results(dds, contrast=c("condition",'AL','AS'), tidy=TRUE)
resADL=na.omit(resADL)
lipADL=resADL[resADL$row%in%lipidGene$V1,]
DegADL=lipADL[lipADL$padj<0.05,]
resAPM <- results(dds, contrast=c("condition",'AP','AS'), tidy=TRUE)
resAPM=na.omit(resAPM)
lipAPM=resAPM[resAPM$row%in%lipidGene$V1,]
DegAPM=lipAPM[lipAPM$padj<0.05,]
resHDL <- results(dds, contrast=c("condition",'HL','HS'), tidy=TRUE)
resHDL=na.omit(resHDL)
lipHDL=resHDL[resHDL$row%in%lipidGene$V1,]
DegHDL=lipHDL[lipHDL$padj<0.05,]
resHPM <- results(dds, contrast=c("condition",'HP','HS'), tidy=TRUE)
resHPM=na.omit(resHPM)
lipHPM=resHPM[resHPM$row%in%lipidGene$V1,]
DegHPM=lipHPM[lipHPM$padj<0.05,]

###########vivo############
vivoExp=read.csv("E:/项目/朱老师/vivo/vivoExp.csv")
lipidGene=read.csv("E:/项目/脂代谢/GENE.csv", header = FALSE)
lipidVivo=vivoExp[vivoExp$GeneName%in%lipidGene$V1,]
lipidVivo=aggregate(lipidVivo[,2:38],list(lipidVivo$GeneName),max)
rownames(lipidVivo)=lipidVivo$Group.1
lipidVivo=lipidVivo[,-1]
library(DESeq2)
condition = factor(c(rep(c("CB"),5),rep(c("CD","CPD","CPM"),each=4),
                   rep(c("MB","MD","MPD","MPM"),each=5)))
colData <- data.frame(row.names =colnames(lipidVivo),condition)
dds <-DESeqDataSetFromMatrix(lipidVivo,colData,design=~condition)
dds <- DESeq(dds)
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
plotPCA(vsd, "condition")
lipVst <- as.data.frame(assay(vsd))
library(factoextra)
sampleDists <- dist(t(data.frame(lipVst[,c(1:9)])))
res <- hcut(sampleDists, k = 2, stand = TRUE)
fviz_dend(res,rect = TRUE,rect_border="cluster",rect_lty=2,lwd=1.0,
          rect_fill = T,cex = 1,color_labels_by_k=T,horiz=T)
sampleDists <- dist(t(data.frame(lipVst[,c(1:5,14:17)])))
res <- hcut(sampleDists, k = 2, stand = TRUE)
fviz_dend(res,rect = TRUE,rect_border="cluster",rect_lty=2,lwd=1.0,
          rect_fill = T,cex = 1,color_labels_by_k=T,horiz=T)
sampleDists <- dist(t(data.frame(lipVst[,c(18:27)])))
res <- hcut(sampleDists, k = 2, stand = TRUE)
fviz_dend(res,rect = TRUE,rect_border="cluster",rect_lty=2,lwd=1.0,
          rect_fill = T,cex = 1,color_labels_by_k=T,horiz=T)
sampleDists <- dist(t(data.frame(lipVst[,c(18:22,33:37)])))
res <- hcut(sampleDists, k = 2, stand = TRUE)
fviz_dend(res,rect = TRUE,rect_border="cluster",rect_lty=2,lwd=1.0,
          rect_fill = T,cex = 1,color_labels_by_k=T,horiz=T)
lipVst=lipVst[,-c(15,16)]
par(cex=0.7)
par(mar=c(1,3,1,1))
n.sample=ncol(lipVst)
cols=rainbow(n.sample*1.2)
par(mfrow=c(2,1))
boxplot(lipVst,col=cols,las=2)
resCDL <- results(dds, contrast=c("condition",'CD','CB'), tidy=TRUE)
resCDL=na.omit(resCDL)
DegCDL=resCDL[resCDL$padj<0.05,]
resCPM <- results(dds, contrast=c("condition",'CPM','CB'), tidy=TRUE)
resCPM=na.omit(resCPM)
DegCPM=resCPM[resCPM$padj<0.05,]
box=as.data.frame(t(lipVst[c("AGPAT4","ACACA","SREBF1","SREBF2","LSS","FASN"),]))
box$group=condition[-c(15,16)]
library(ggstatsplot)
ggbetweenstats(data=box[c(1:9),],x = group,y = AGPAT4) 
ggbetweenstats(data=box[c(1:9),],x = group,y = ACACA) 
ggbetweenstats(data=box[c(1:9),],x = group,y = SREBF1) 
ggbetweenstats(data=box[c(1:9),],x = group,y = SREBF2) 
ggbetweenstats(data=box[c(1:9),],x = group,y = LSS) 
ggbetweenstats(data=box[c(1:9),],x = group,y = FASN) 
ggbetweenstats(data=box[c(1:5,14:15),],x = group,y = AGPAT4) 
ggbetweenstats(data=box[c(1:5,14:15),],x = group,y = ACACA) 
ggbetweenstats(data=box[c(1:5,14:15),],x = group,y = SREBF1) 
ggbetweenstats(data=box[c(1:5,14:15),],x = group,y = SREBF2) 
ggbetweenstats(data=box[c(1:5,14:15),],x = group,y = LSS) 
ggbetweenstats(data=box[c(1:5,14:15),],x = group,y = FASN) 
ggbetweenstats(data=box[c(16:25),],x = group,y = AGPAT4) 
ggbetweenstats(data=box[c(16:25),],x = group,y = ACACA) 
ggbetweenstats(data=box[c(16:25),],x = group,y = SREBF1) 
ggbetweenstats(data=box[c(16:25),],x = group,y = SREBF2) 
ggbetweenstats(data=box[c(16:25),],x = group,y = LSS) 
ggbetweenstats(data=box[c(16:25),],x = group,y = FASN) 
ggbetweenstats(data=box[c(16:20,31:35),],x = group,y = AGPAT4) 
ggbetweenstats(data=box[c(16:20,31:35),],x = group,y = ACACA) 
ggbetweenstats(data=box[c(16:20,31:35),],x = group,y = SREBF1) 
ggbetweenstats(data=box[c(16:20,31:35),],x = group,y = SREBF2) 
ggbetweenstats(data=box[c(16:20,31:35),],x = group,y = LSS) 
ggbetweenstats(data=box[c(16:20,31:35),],x = group,y = FASN) 
write.csv(DegCDL,"E:/项目/朱老师/vivo/DegCDL.csv")
write.csv(DegCPM,"E:/项目/朱老师/vivo/DegCPM.csv")
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(dplyr)
resDL$threshold[resDL$padj<0.05&resDL$log2FoldChange>0]="up"
resDL$threshold[resDL$padj<0.05&resDL$log2FoldChange<(-0)]="down"
resDL$threshold[(resDL$padj>0.05)|(resDL$log2FoldChange<=0)&resDL$log2FoldChange>=(-0)]="normal"
resDL$threshold=factor(resDL$threshold,levels = c("up","down","normal"),ordered = T)
x_lim=max(resDL$log2FoldChange,-resDL$log2FoldChange)
theme_set(theme_bw())
p <- ggplot(resDL,aes(log2FoldChange,-1*log10(padj),
                      color =threshold))+geom_point()+
  xlim(-5,5) +  labs(x="log2(FoldChange)",y="-log10(padj)")
p <- p + scale_color_manual(values =c('up'="red","normal"="grey","down"="blue"))+
  geom_hline(yintercept=-log10(0.05),linetype=4)+
  geom_vline(xintercept=c(-log2(1),log2(1)),linetype=4)
p <- p +theme(plot.title = element_text(size = 25,face = 'bold', vjust = 0.5, hjust = 0.5),
              legend.title = element_blank(),
              legend.text = element_text(size = 18, face = 'bold'),
              legend.position = 'right',
              legend.key.size=unit(0.8,'cm'),
              axis.ticks.x=element_blank(),
              axis.text.x=element_text(size = 15,face = 'bold', vjust = 0.5, hjust = 0.5),
              axis.text.y=element_text(size = 15,face = 'bold', vjust = 0.5, hjust = 0.5),
              axis.title.x = element_text(size = 20,face = 'bold', vjust = 0.5, hjust = 0.5),
              axis.title.y = element_text(size = 20,face = 'bold', vjust = 0.5, hjust = 0.5),
              panel.background = element_rect(fill = 'transparent',colour = 'black'),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              plot.background = element_rect(fill = 'transparent',colour = 'black'))

library(pheatmap)
heatdata <- lipVst[DegDL$row,c(9:16,25:32)]
annotation_col = data.frame(Type = factor(rep(c("DL","contral"), c(8,8))))
rownames(annotation_col) = colnames(heatdata)
pheatmap(t(scale(t(heatdata))),cluster_rows = T,show_rownames = F,show_colnames = F,
         cluster_cols = T,annotation_col = annotation_col,
         color = colorRampPalette(c("navy","white","firebrick3"))(100),
         legend = T,legend_breaks = c(-2,0,2),breaks = unique(seq(-2,2,length=100)),
         cutree_rows = 2,cutree_cols = 2,
         lagend_labels = c("≤2","0","≥2"))

resPM$threshold[resPM$padj<0.05&resPM$log2FoldChange>0]="up"
resPM$threshold[resPM$padj<0.05&resPM$log2FoldChange<(-0)]="down"
resPM$threshold[(resPM$padj>0.05)|(resPM$log2FoldChange<=0)&resPM$log2FoldChange>=(-0)]="normal"
resPM$threshold=factor(resPM$threshold,levels = c("up","down","normal"),ordered = T)
x_lim=max(resPM$log2FoldChange,-resPM$log2FoldChange)
theme_set(theme_bw())
p <- ggplot(resPM,aes(log2FoldChange,-1*log10(padj),
                      color =threshold ))+geom_point()+
  xlim(-5,5) +  labs(x="log2(FoldChange)",y="-log10(padj)")
p <- p + scale_color_manual(values =c('up'="red","normal"="grey","down"="blue"))+
  geom_hline(yintercept=-log10(0.05),linetype=4)+
  geom_vline(xintercept=c(-log2(1),log2(1)),linetype=4)
p <- p +theme(plot.title = element_text(size = 25,face = 'bold', vjust = 0.5, hjust = 0.5),
              legend.title = element_blank(),
              legend.text = element_text(size = 18, face = 'bold'),
              legend.position = 'right',
              legend.key.size=unit(0.8,'cm'),
              axis.ticks.x=element_blank(),
              axis.text.x=element_text(size = 15,face = 'bold', vjust = 0.5, hjust = 0.5),
              axis.text.y=element_text(size = 15,face = 'bold', vjust = 0.5, hjust = 0.5),
              axis.title.x = element_text(size = 20,face = 'bold', vjust = 0.5, hjust = 0.5),
              axis.title.y = element_text(size = 20,face = 'bold', vjust = 0.5, hjust = 0.5),
              panel.background = element_rect(fill = 'transparent',colour = 'black'),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              plot.background = element_rect(fill = 'transparent',colour = 'black'))
p

heatdata <- lipVst[DegPM$row,c(17:24,25:32)]
annotation_col = data.frame(Type = factor(rep(c("PM2.5","control"), c(8,8))))
rownames(annotation_col) = colnames(heatdata)
pheatmap(t(scale(t(heatdata))),cluster_rows = T,show_rownames = F,show_colnames = F,
         cluster_cols = T,annotation_col = annotation_col,
         color = colorRampPalette(c("navy","white","firebrick3"))(100),
         legend = T,legend_breaks = c(-2,0,2),breaks = unique(seq(-2,2,length=100)),
         cutree_rows = 2,cutree_cols = 2,
         lagend_labels = c("≤2","0","≥2"))

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
sign =c()
patient_id = colnames(lipidExp)
for (i in 1:length(colnames(lipidExp))){
  tmp = (strsplit2(as.character(patient_id[i]),split = "-"))
  sign[i] = tmp[4]
  tmp = paste(tmp[1],tmp[2],tmp[3],sep = "-")
  patient_id[i] = tmp
}
colnames(lipidExp)=patient_id
ExpCan = lipidExp[,sign == "01A"]
ExpNor = lipidExp[,sign == "11A"]
SAMPID=colnames(ExpCan)
ExpCan=rbind(ExpCan,SAMPID)
ExpCan=data.frame(t(ExpCan))
ExpCan <- aggregate(x = ExpCan[,1:1038], by = list(ExpCan$SAMPID),FUN = max)
SAMPID=colnames(ExpNor)
ExpNor=rbind(ExpNor,SAMPID)
ExpNor=data.frame(t(ExpNor))
ExpNor <- aggregate(x = ExpNor[,1:1038], by = list(ExpNor$SAMPID),FUN = max)
cli.query <- GDCquery(project =c("TCGA-LUAD"),file.type = "xml",
                      data.category = "Clinical") 
GDCdownload(cli.query)
cli <- GDCprepare_clinic(cli.query, clinical.info = "patient")
samples=cli[,c("bcr_patient_barcode","gender","vital_status",
               "days_to_death","days_to_last_followup","race_list",
               "age_at_initial_pathologic_diagnosis",
               "tobacco_smoking_history","stage_event_pathologic_stage",
               "stage_event_tnm_categories")]
samples=samples[!(is.na(samples$days_to_last_followup)&
                    is.na(samples$days_to_death)),]
samples[samples$vital_status=="Dead","survival_time"]=
  samples[samples$vital_status=="Dead","days_to_death"]
samples[samples$vital_status=="Alive","survival_time"]=
  samples[samples$vital_status=="Alive","days_to_last_followup"]
samples$status=ifelse(samples$vital_status=="Dead",1,0)
samples=samples[!(substring(samples$stage_event_pathologic_stage,7,9)==""),]
samples$stage=ifelse(substring(samples$stage_event_pathologic_stage,7,9)=="IV"|
                       substring(samples$stage_event_pathologic_stage,7,9)=="III",2,1)
samples=samples[!duplicated(samples),]
samples=samples[!(samples$race_list==""),]
ExpCanCli=merge(ExpCan,samples,by.x="Group.1",by.y="bcr_patient_barcode")
ExpCanCli$Group.1=paste(ExpCanCli$Group.1,"01A",sep = "-")
ExpNorCli=merge(ExpNor,samples,by.x="Group.1",by.y="bcr_patient_barcode")
ExpNorCli$Group.1=paste(ExpNorCli$Group.1,"11A",sep = "-")
ExpCli=rbind(ExpCanCli,ExpNorCli)
ExpCli$group=ifelse(substring(ExpCli$Group.1,14,16)=="01A","cancer","normal")
rownames(ExpCli)=ExpCli[,1]
library(DESeq2)
ExpM=t(ExpCli[,-1])
ExpM=matrix(as.numeric(unlist(ExpM[1:1038,])),nrow=nrow(ExpM[1:1038,]), dimnames = list(rownames(ExpM[1:1038,]),colnames(ExpM[1:1038,])))
ExpM=ExpM[rowSums(ExpM)>1,]
condition <- factor(ExpCli[,1052])
colData <- data.frame(row.names =colnames(ExpM),condition)
dds <-DESeqDataSetFromMatrix(ExpM,colData,design=~condition)
dds <- DESeq(dds)
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
library(tableone)
stable <- CreateTableOne(vars=c("gender","race_list","age_at_initial_pathologic_diagnosis","tobacco_smoking_history"),
                         strata="group",data=ExpCli,
                         factorVars=c("gender","race_list","tobacco_smoking_history"))
stable=print(stable,showAllLevels = TRUE)
write.csv(stable,"E:/项目/base2b/TCGA_LUAD/table1.csv")
plotPCA(vsd, "condition")
expVst <- as.data.frame(assay(vsd))
par(cex=0.7)
par(mar=c(1,3,1,1))
n.sample=ncol(expVst)
cols=rainbow(n.sample*1.2)
par(mfrow=c(2,1))
boxplot(expVst,col=cols,las=2)
res <- results(dds, contrast=c("condition",'cancer','normal'),tidy=TRUE) 
res=na.omit(res)
Deg=res[res$padj<0.05&abs(res$log2FoldChange)>log2(1),]
write.csv(Deg,"E:/项目/base2b/TCGA_LUAD/Deg.csv")
res$threshold[res$padj<0.05&res$log2FoldChange>0]="up"
res$threshold[res$padj<0.05&res$log2FoldChange<(-0)]="down"
res$threshold[(res$padj>0.05)|(res$log2FoldChange<=0)&res$log2FoldChange>=(-0)]="normal"
res$threshold=factor(res$threshold,levels = c("up","down","normal"),ordered = T)
x_lim=max(res$log2FoldChange,-res$log2FoldChange)
theme_set(theme_bw())
p <- ggplot(res,aes(log2FoldChange,-1*log10(padj),
                    color =threshold))+geom_point()+
  xlim(-8,8) +  labs(x="log2(FoldChange)",y="-log10(padj)")
p <- p + scale_color_manual(values =c('up'="red","normal"="grey","down"="blue"))+
  geom_hline(yintercept=-log10(0.05),linetype=4)+
  geom_vline(xintercept=c(-log2(1),log2(1)),linetype=4)
p <- p +theme(plot.title = element_text(size = 25,face = 'bold', vjust = 0.5, hjust = 0.5),
              legend.title = element_blank(),
              legend.text = element_text(size = 18, face = 'bold'),
              legend.position = 'right',
              legend.key.size=unit(0.8,'cm'),
              axis.ticks.x=element_blank(),
              axis.text.x=element_text(size = 15,face = 'bold', vjust = 0.5, hjust = 0.5),
              axis.text.y=element_text(size = 15,face = 'bold', vjust = 0.5, hjust = 0.5),
              axis.title.x = element_text(size = 20,face = 'bold', vjust = 0.5, hjust = 0.5),
              axis.title.y = element_text(size = 20,face = 'bold', vjust = 0.5, hjust = 0.5),
              panel.background = element_rect(fill = 'transparent',colour = 'black'),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              plot.background = element_rect(fill = 'transparent',colour = 'black'))
p
heatdata <- expVst[Deg$row,]
annotation_col = data.frame(Type = factor(rep(c("cancer","normal"), c(436,56))))
rownames(annotation_col) = colnames(heatdata)
pheatmap(t(scale(t(heatdata))),cluster_rows = T,show_rownames = F,show_colnames = F,
         cluster_cols = F,annotation_col = annotation_col,
         color = colorRampPalette(c("navy","white","firebrick3"))(100),
         legend = T,legend_breaks = c(-2,0,2),breaks = unique(seq(-2,2,length=100)),
         cutree_rows = 2,cutree_cols = 2,
         lagend_labels = c("≤2","0","≥2"))

DegUp=intersect(DegPM[DegPM$log2FoldChange>0,]$row,Deg[Deg$log2FoldChange>0,]$row)
DegDown=intersect(DegPM[DegPM$log2FoldChange<0,]$row,Deg[Deg$log2FoldChange<0,]$row)
DegUpDown=intersect(DegUp,DegDL[DegDL$log2FoldChange<0,]$row)
DegDownUp=intersect(DegDown,DegDL[DegDL$log2FoldChange>0,]$row)
keyGene=union(DegUpDown,DegDownUp)
################LUSC###########
query.exp.hg38 <- GDCquery(project = "TCGA-LUSC", 
                           data.category = "Transcriptome Profiling", 
                           data.type = "Gene Expression Quantification", 
                           workflow.type = "HTSeq - Counts")
RnaseqSE <- GDCdownload(query.exp.hg38)
RnaseqSE <- GDCprepare(query.exp.hg38)
rownames(RnaseqSE) <- values(RnaseqSE)$external_gene_name
exp.hg38.values <- assay(RnaseqSE)
Exp = exp.hg38.values[,substring(colnames(exp.hg38.values),14,16)=="01A"|
                        substring(colnames(exp.hg38.values),14,16)=="11A"]
lipidExp=Exp[rownames(Exp)%in%lipidGene$V1,]
sign =c()
patient_id = colnames(lipidExp)
for (i in 1:length(colnames(lipidExp))){
  tmp = (strsplit2(as.character(patient_id[i]),split = "-"))
  sign[i] = tmp[4]
  tmp = paste(tmp[1],tmp[2],tmp[3],sep = "-")
  patient_id[i] = tmp
}
colnames(lipidExp)=patient_id
ExpCan = lipidExp[,sign == "01A"]
ExpNor = lipidExp[,sign == "11A"]
SAMPID=colnames(ExpCan)
ExpCan=rbind(ExpCan,SAMPID)
ExpCan=data.frame(t(ExpCan))
ExpCan <- aggregate(x = ExpCan[,1:1038], by = list(ExpCan$SAMPID),FUN = max)
SAMPID=colnames(ExpNor)
ExpNor=rbind(ExpNor,SAMPID)
ExpNor=data.frame(t(ExpNor))
ExpNor <- aggregate(x = ExpNor[,1:1038], by = list(ExpNor$SAMPID),FUN = max)
cli.query <- GDCquery(project =c("TCGA-LUSC"),file.type = "xml",
                      data.category = "Clinical") 
GDCdownload(cli.query)
cli <- GDCprepare_clinic(cli.query, clinical.info = "patient")
samples=cli[,c("bcr_patient_barcode","gender","vital_status",
               "days_to_death","days_to_last_followup","race_list",
               "age_at_initial_pathologic_diagnosis",
               "tobacco_smoking_history","stage_event_pathologic_stage",
               "stage_event_tnm_categories")]
samples=samples[!(is.na(samples$days_to_last_followup)&
                    is.na(samples$days_to_death)),]
samples[samples$vital_status=="Dead","survival_time"]=
  samples[samples$vital_status=="Dead","days_to_death"]
samples[samples$vital_status=="Alive","survival_time"]=
  samples[samples$vital_status=="Alive","days_to_last_followup"]
samples$status=ifelse(samples$vital_status=="Dead",1,0)
samples=samples[!(substring(samples$stage_event_pathologic_stage,7,9)==""),]
samples$stage=ifelse(substring(samples$stage_event_pathologic_stage,7,9)=="IV"|
                       substring(samples$stage_event_pathologic_stage,7,9)=="III",2,1)
samples=samples[!duplicated(samples),]
samples=samples[!(samples$race_list==""),]
ExpCanCli=merge(ExpCan,samples,by.x="Group.1",by.y="bcr_patient_barcode")
ExpCanCli$Group.1=paste(ExpCanCli$Group.1,"01A",sep = "-")
ExpNorCli=merge(ExpNor,samples,by.x="Group.1",by.y="bcr_patient_barcode")
ExpNorCli$Group.1=paste(ExpNorCli$Group.1,"11A",sep = "-")
ExpCli=rbind(ExpCanCli,ExpNorCli)
ExpCli$group=ifelse(substring(ExpCli$Group.1,14,16)=="01A","cancer","normal")
rownames(ExpCli)=ExpCli[,1]
library(DESeq2)
ExpM=t(ExpCli[,-1])
ExpM=matrix(as.numeric(unlist(ExpM[1:1038,])),nrow=nrow(ExpM[1:1038,]), dimnames = list(rownames(ExpM[1:1038,]),colnames(ExpM[1:1038,])))
ExpM=ExpM[rowSums(ExpM)>1,]
condition <- factor(ExpCli[,1052])
colData <- data.frame(row.names =colnames(ExpM),condition)
dds <-DESeqDataSetFromMatrix(ExpM,colData,design=~condition)
dds <- DESeq(dds)
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
library(tableone)
stable <- CreateTableOne(vars=c("gender","race_list","age_at_initial_pathologic_diagnosis","tobacco_smoking_history"),
                         strata="group",data=ExpCli,
                         factorVars=c("gender","race_list","tobacco_smoking_history"))
stable=print(stable,showAllLevels = TRUE)
write.csv(stable,"E:/项目/base2b/TCGA_LUSC/table1.csv")
plotPCA(vsd, "condition")
expVst <- as.data.frame(assay(vsd))
par(cex=0.7)
par(mar=c(1,3,1,1))
n.sample=ncol(expVst)
cols=rainbow(n.sample*1.2)
par(mfrow=c(2,1))
boxplot(expVst,col=cols,las=2)
res <- results(dds, contrast=c("condition",'cancer','normal'),tidy=TRUE) 
res=na.omit(res)
Deg=res[res$padj<0.05&abs(res$log2FoldChange)>log2(1),]
write.csv(Deg,"E:/项目/base2b/TCGA_LUSC/Deg.csv")
res$threshold[res$padj<0.05&res$log2FoldChange>0]="up"
res$threshold[res$padj<0.05&res$log2FoldChange<(-0)]="down"
res$threshold[(res$padj>0.05)|(res$log2FoldChange<=0)&res$log2FoldChange>=(-0)]="normal"
res$threshold=factor(res$threshold,levels = c("up","down","normal"),ordered = T)
x_lim=max(res$log2FoldChange,-res$log2FoldChange)
theme_set(theme_bw())
p <- ggplot(res,aes(log2FoldChange,-1*log10(padj),
                    color =threshold))+geom_point()+
  xlim(-9,9) +  labs(x="log2(FoldChange)",y="-log10(padj)")
p <- p + scale_color_manual(values =c('up'="red","normal"="grey","down"="blue"))+
  geom_hline(yintercept=-log10(0.05),linetype=4)+
  geom_vline(xintercept=c(-log2(1),log2(1)),linetype=4)
p <- p +theme(plot.title = element_text(size = 25,face = 'bold', vjust = 0.5, hjust = 0.5),
              legend.title = element_blank(),
              legend.text = element_text(size = 18, face = 'bold'),
              legend.position = 'right',
              legend.key.size=unit(0.8,'cm'),
              axis.ticks.x=element_blank(),
              axis.text.x=element_text(size = 15,face = 'bold', vjust = 0.5, hjust = 0.5),
              axis.text.y=element_text(size = 15,face = 'bold', vjust = 0.5, hjust = 0.5),
              axis.title.x = element_text(size = 20,face = 'bold', vjust = 0.5, hjust = 0.5),
              axis.title.y = element_text(size = 20,face = 'bold', vjust = 0.5, hjust = 0.5),
              panel.background = element_rect(fill = 'transparent',colour = 'black'),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              plot.background = element_rect(fill = 'transparent',colour = 'black'))
p
heatdata <- expVst[Deg$row,]
annotation_col = data.frame(Type = factor(rep(c("cancer","normal"), c(380,44))))
rownames(annotation_col) = colnames(heatdata)
pheatmap(t(scale(t(heatdata))),cluster_rows = T,show_rownames = F,show_colnames = F,
         cluster_cols = F,annotation_col = annotation_col,
         color = colorRampPalette(c("navy","white","firebrick3"))(100),
         legend = T,legend_breaks = c(-2,0,2),breaks = unique(seq(-2,2,length=100)),
         cutree_rows = 2,cutree_cols = 2,
         lagend_labels = c("≤2","0","≥2"))
#########################clusterProfiler######
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(GOplot)
genelist=mapIds(org.Hs.eg.db,keyGene,keytype="SYMBOL",column="ENTREZID")
go=enrichGO(gene=genelist,OrgDb = "org.Hs.eg.db",keyType = "ENTREZID",pvalueCutoff = 0.05,readable = T,ont="ALL")
write.csv(go,"E:/项目/base2b/DEGs/GoKeyGene.csv")
kegg=enrichKEGG(genelist,organism = "hsa",pvalueCutoff = 0.05)
kegg=setReadable(kegg,org.Hs.eg.db,keyType = "ENTREZID")
write.csv(kegg,"E:/项目/base2b/DEGs/KeggKeyGene.csv")
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
chord=chord_dat(data=circ2,genes = gene,process = unique(circ2$term))
GOBubble(circ2,labels = 2.4)
kegg1=data.frame(kegg)
kegg1=kegg1[,c(1,2,7,8)]
kegg1$geneID=str_replace_all(kegg1$geneID,"/",",")
names(kegg1)=c("ID","Term","adj_pval","Genes")
kegg1$Category=c("KEGG pathway")
plot_data=list(DA=kegg1,GE=gene)
circ2=data.frame()
circ2=circle_dat(plot_data$DA,plot_data$GE)
GOBubble(circ2,labels = 0)

logFC=read.csv("E:/项目/base2b/DEGs/logFC.csv")
heatdata <- logFC
rownames(heatdata)=heatdata[,1]
heatdata=heatdata[,-1]
pheatmap(scale(heatdata), 
         cluster_rows = T,
         cluster_cols = T,
         annotation_legend  = T, 
         show_rownames = T,
         show_colnames = T,
         scale = "row",
         color = colorRampPalette(c("navy","white","firebrick3"))(100),
         cutree_rows = 2,cutree_cols = 2,
         cellwidth = 25, cellheight = 10,
         fontsize = 10)


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
write.csv(deg1,"E:/项目/base2b/GSE93329/deg1.csv")
