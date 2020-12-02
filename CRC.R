library(SummarizedExperiment)
library(TCGAbiolinks)
library(limma)
query.exp.hg38 <- GDCquery(project = "TCGA-COAD", 
                           data.category = "Transcriptome Profiling", 
                           data.type = "Gene Expression Quantification", 
                           workflow.type = "HTSeq - Counts")
COADRnaseqSE <- GDCdownload(query.exp.hg38)
COADRnaseqSE <- GDCprepare(query.exp.hg38)
rownames(COADRnaseqSE) <- rowData(COADRnaseqSE)$external_gene_name
exp.hg38.values <- assay(COADRnaseqSE)
ExpC = exp.hg38.values[,substring(colnames(exp.hg38.values),14,16)=="01A"]
cli.query <- GDCquery(project =c("TCGA-COAD"),file.type = "xml",
                      data.category = "Clinical") 
GDCdownload(cli.query)
cliC <- GDCprepare_clinic(cli.query, clinical.info = "patient")
samplesC=cliC[,c("bcr_patient_barcode","gender","vital_status","days_to_death","weight",
               "days_to_last_followup","race_list","age_at_initial_pathologic_diagnosis",
               "stage_event_pathologic_stage","height","stage_event_tnm_categories",
               "microsatellite_instability","lymphatic_invasion","tumor_tissue_site")]
query.exp.hg38 <- GDCquery(project = "TCGA-READ", 
                           data.category = "Transcriptome Profiling", 
                           data.type = "Gene Expression Quantification", 
                           workflow.type = "HTSeq - Counts")
READRnaseqSE <- GDCdownload(query.exp.hg38)
READRnaseqSE <- GDCprepare(query.exp.hg38)
rownames(READRnaseqSE) <- values(READRnaseqSE)$external_gene_name
exp.hg38.values <- assay(READRnaseqSE)
ExpR = exp.hg38.values[,substring(colnames(exp.hg38.values),14,16)=="01A"]
cli.query <- GDCquery(project =c("TCGA-READ"),file.type = "xml",
                      data.category = "Clinical") 
GDCdownload(cli.query)
cliR <- GDCprepare_clinic(cli.query, clinical.info = "patient")
samplesR=cliR[,c("bcr_patient_barcode","gender","vital_status","days_to_death","weight",
               "days_to_last_followup","race_list","age_at_initial_pathologic_diagnosis",
               "stage_event_pathologic_stage","height","stage_event_tnm_categories",
               "microsatellite_instability","lymphatic_invasion","tumor_tissue_site")]
Exp=cbind(ExpC,ExpR)
samples=rbind(samplesC,samplesR)
samples=samples[!duplicated(samples$bcr_patient_barcode),]
MSI=read.csv("~/Downloads/项目/CRC/tcga/MSI.csv")
MSI_CLI=merge(MSI,samples,by.x="patient",by.y="bcr_patient_barcode")
MSI_CLI$survival_time<-ifelse(MSI_CLI$vital_status=='Alive',MSI_CLI$days_to_last_followup,samples$days_to_death)
MSI_CLI=MSI_CLI[!(MSI_CLI$microsatellite_instability=="NO"),]
MSI_CLI$stage=ifelse(substring(MSI_CLI$stage_event_pathologic_stage,7,9)=="IV"|
                       substring(MSI_CLI$stage_event_pathologic_stage,7,9)=="III","Advanced stage","Early stage")
library(stringr)
MSI_CLI$metastasis=ifelse(str_sub(MSI_CLI$stage_event_tnm_categories,-1)=="1"|
                          str_sub(MSI_CLI$lymphatic_invasion,1)=="YES","Metastasis","Non-metastasis")
samples1<-MSI_CLI[match(substr(colnames(Exp),1,12),MSI_CLI$patient),]
all(substr(rownames(samples1),1,12)==substr(colnames(Exp),1,12))
rownames(samples1)<-colnames(Exp)
samples1=samples1[na.omit(samples1$patient),]
ExpCli<-Exp[,rownames(samples1)]
library(tableone)
stable1 <- CreateTableOne(vars=c("gender","race_list","age_at_initial_pathologic_diagnosis",
                                 "stage","vital_status"),
                          data=samples1,strata = "metastasis",
                          factorVars=c("gender","race_list","stage","vital_status"))
stable1=print(stable1,showAllLevels = TRUE)
write.csv(stable1,"~/Downloads/项目/CRC/tcga/table1.csv")
library(DESeq2)
condition <- factor(samples1$metastasis)
colData <- data.frame(row.names =samples1$patient,condition)
dds <-DESeqDataSetFromMatrix(ExpCli,colData,design=~condition)
dds <- DESeq(dds)
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
plotPCA(vsd, "condition")
expVst <- as.data.frame(assay(vsd))
res <- results(dds, contrast=c("condition",'Metastasis','Non-metastasis'),tidy=TRUE) 
res=na.omit(res)
Deg=res[res$padj<0.05&abs(res$log2FoldChange)>0,]
write.csv(Deg,"~/Downloads/项目/CRC/tcga/Deg.csv")
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(dplyr)
res$threshold[res$padj<0.05&res$log2FoldChange>0]="up"
res$threshold[res$padj<0.05&res$log2FoldChange<(-0)]="down"
res$threshold[(res$padj>0.05)|(res$log2FoldChange<=0)&res$log2FoldChange>=(-0)]="normal"
res$threshold=factor(res$threshold,levels = c("up","down","normal"),ordered = T)
p <- ggplot(data = res, 
            aes(x = log2FoldChange, 
                y = -log10(padj))) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=threshold)) +
  scale_color_manual(values=c("red","blue", "grey"))+
  geom_vline(xintercept=c(-0,0),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  theme_bw()
for_label <- res %>% 
  filter(abs(log2FoldChange) >2& -log10(padj)> -log10(0.01))
p +
  geom_point(size = 3, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = row),
    data = for_label,
    color="black"
  )
library(pheatmap)
heatdata <- expVst[Deg$row,]
annotation_col = colData
annotation_col$patient=rownames(annotation_col)
annotation_col=annotation_col[order(annotation_col$condition),]
heatdata=heatdata[,rownames(annotation_col)]
heatdata=log2(heatdata+1)
annotation_col=as.data.frame(annotation_col[,1])
rownames(annotation_col)=colnames(heatdata)
colnames(annotation_col)=c("condition")
pheatmap(heatdata,cluster_rows = T,show_rownames = T,show_colnames = F,
         cluster_cols = F,annotation_col = annotation_col,
         scale = "row",
         color =colorRampPalette(c("blue", "white","red"))(100),
         legend = T)

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(GOplot)
genelist=mapIds(org.Hs.eg.db,Deg[Deg$log2FoldChange>0,]$row,keytype="SYMBOL",column="ENTREZID")
go=enrichGO(gene=genelist,OrgDb = "org.Hs.eg.db",keyType = "ENTREZID",pvalueCutoff = 0.05,readable = T,ont="ALL")
write.csv(go,"~/Downloads/项目/CRC/tcga/GoUp.csv")
kegg=enrichKEGG(genelist,organism = "hsa",pvalueCutoff = 0.05)
kegg=setReadable(kegg,org.Hs.eg.db,keyType = "ENTREZID")
write.csv(kegg,"~/Downloads/项目/CRC/tcga/KeggUp.csv")
library(enrichplot)
dotplot(go, showCategory=10)
dotplot(kegg, showCategory=10)
genelist=mapIds(org.Hs.eg.db,Deg[Deg$log2FoldChange<0,]$row,keytype="SYMBOL",column="ENTREZID")
go=enrichGO(gene=genelist,OrgDb = "org.Hs.eg.db",keyType = "ENTREZID",pvalueCutoff = 0.05,readable = T,ont="ALL")
write.csv(go,"~/Downloads/项目/CRC/tcga/GoDown.csv")
kegg=enrichKEGG(genelist,organism = "hsa",pvalueCutoff = 0.05)
kegg=setReadable(kegg,org.Hs.eg.db,keyType = "ENTREZID")
write.csv(kegg,"~/Downloads/项目/CRC/tcga/KeggDown.csv")
dotplot(go, showCategory=10)
dotplot(kegg, showCategory=10)
genelist=mapIds(org.Hs.eg.db,Deg$row,keytype="SYMBOL",column="ENTREZID")
go=enrichGO(gene=genelist,OrgDb = "org.Hs.eg.db",keyType = "ENTREZID",pvalueCutoff = 0.05,readable = T,ont="ALL")
write.csv(go,"~/Downloads/项目/CRC/tcga/GoAll.csv")
kegg=enrichKEGG(genelist,organism = "hsa",pvalueCutoff = 0.05)
kegg=setReadable(kegg,org.Hs.eg.db,keyType = "ENTREZID")
write.csv(kegg,"~/Downloads/项目/CRC/tcga/KeggAll.csv")
dotplot(go, showCategory=10)
dotplot(kegg, showCategory=10)

library("survival")
library("survminer")
DegExp=expVst[rownames(expVst)%in%Deg$row,]
surExp=matrix(nrow = 80,ncol = 73)
for (i in 1:length(rownames(DegExp))){
  group=ifelse(DegExp[i,]>median(as.numeric(DegExp[i,])),1,0)
  surExp[,i]=group
}
rownames(surExp)=colnames(expVst[,1:80])
colnames(surExp)=rownames(DegExp)
surExp=cbind(surExp,samples1)
surExp$status=ifelse(surExp$vital_status=="Alive",0,1)
surExp$survival_time<-ifelse(surExp$vital_status=='Alive',surExp$days_to_last_followup,surExp$days_to_death)
P=list()
HR=list()
CI=matrix(ncol = 2,nrow = 73)
colnames(CI) <- c("Lower", "Higher")
for (i in 1:73){
  fit =  coxph(Surv(as.numeric(survival_time), as.numeric(status))~surExp[,i],data=surExp)
  HR[i] <- round(exp(coef(fit)), 2)
  CI[i,] <- round(exp(confint(fit)), 2)
  P[i] <- round(coef(summary(fit))[,5], 3)
}
surTab=as.data.frame(cbind(HR, CI, P))
rownames(surTab)=colnames(surExp[,1:73])
surTab$HR=unlist(surTab$HR)
surTab$Lower=unlist(surTab$Lower)
surTab$Higher=unlist(surTab$Higher)
surTab$P=unlist(surTab$P)
surGene=surTab[surTab$P<0.05,]
write.csv(surTab,"~/Downloads/项目/CRC/tcga/sur.csv")
######WGCNA######

library(WGCNA)
library(reshape2)
library(stringr)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
exprMat <- heatdata
type = "unsigned"
corType = "pearson"
maxPOutliers = ifelse(corType=="pearson",1,0.05)
dataExpr <- as.data.frame(t(exprMat))
gsg = goodSamplesGenes(dataExpr, verbose = 3)
if (!gsg$allOK){
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
sampleTree = hclust(dist(dataExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers, 
                        networkType=type, verbose=5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=sft$powerEstimate,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
     cex=cex1, col="red")
power = sft$powerEstimate
if (is.na(power)){
  power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                 ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                        ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                               ifelse(type == "unsigned", 6, 12))       
                 )
  )
}
net = blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
                       TOMType = type, minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       corType = corType, 
                       maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                       verbose = 3)
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
##################
Deg1=res[res$padj<0.01&abs(res$log2FoldChange)>2,]
Deg1Exp=as.data.frame(t(expVst[rownames(expVst)%in%Deg1$row,]))
Deg1Exp$status=ifelse(colData$condition=="Metastasis",2,1)
library(rms)
dd=datadist(Deg1Exp)
options(datadist="dd") 
f <- lrm(status ~ GABRA3 + CCL25 + BAAT+SLC14A1+LGSN+LCN15+
           C6orf15+LINC01819+LINC00460 , data =  Deg1Exp)
nom <- nomogram(f, fun=plogis, lp=F, funlabel="Risk")
plot(nom)
cindex=rcorrcens(status~predict(f), data = Deg1Exp)
cindex

