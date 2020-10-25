library(SummarizedExperiment)
library(TCGAbiolinks)
library(limma)
query.exp.hg38 <- GDCquery(project = "TCGA-LUAD", 
                           data.category = "Transcriptome Profiling", 
                           data.type = "Gene Expression Quantification", 
                           workflow.type = "HTSeq - FPKM")
LUADRnaseqSE <- GDCdownload(query.exp.hg38)
LUADRnaseqSE <- GDCprepare(query.exp.hg38)
rownames(LUADRnaseqSE) <- values(LUADRnaseqSE)$external_gene_name
exp.hg38.values <- assay(LUADRnaseqSE)
sign =c()
patient_id = colnames(exp.hg38.values)
for (i in 1:length(colnames(exp.hg38.values))){
  tmp = (strsplit2(as.character(patient_id[i]),split = "-"))
  sign[i] = tmp[4]
  tmp = paste(tmp[1],tmp[2],tmp[3],sep = "-")
  patient_id[i] = tmp
}
colnames(exp.hg38.values) = patient_id
ExpCancer = exp.hg38.values[,sign == "01A"]
ExpCancer=data.frame(t(ExpCancer))
ExpCancer$group=ifelse(ExpCancer$PDIA3P1>median(ExpCancer$PDIA3P1),"high","low")
ExpCancer=t(ExpCancer)
colnames(ExpCancer)=ExpCancer[56494,]
ExpCancer1=matrix(as.numeric(unlist(ExpCancer[-56494,])),nrow=nrow(ExpCancer[-56494,]), dimnames = list(rownames(ExpCancer[-56494,]),colnames(ExpCancer[-56494,])))

library(edgeR)
y = DGEList(ExpCancer1,group = colnames(ExpCancer1),genes = ExpCancer1)
keep<-rowSums(cpm(y) > 0.5) >= 2
y<-y[keep, ,keep.lib.sizes=FALSE]
y=calcNormFactors(y)
design <- model.matrix(~0+factor(colnames(ExpCancer1)))  
colnames(design) <- levels(factor(colnames(ExpCancer1))) 
y <- estimateDisp(y, design, robust=TRUE)
fit <- glmQLFit(y, design, robust=TRUE)
B.LvsP <- makeContrasts(high-low, levels=design)
res <- glmQLFTest(fit, contrast=B.LvsP)
is.de <- decideTestsDGE(res)
summary(is.de)
tr <- glmTreat(fit, contrast=B.LvsP, lfc=log2(1.5))
is.de <- decideTestsDGE(tr)
summary(is.de)
deg=tr$table[abs(is.de)==1,]
write.csv(deg,"~/R/LUAD/deg.csv")
plotMD(tr, status=is.de, values=c(1,-1), col=c("red","blue"))
library(ggplot2)
res=cbind(tr$table,tr$genes)
colnames(res)[5:523]=paste(colnames(res)[5:523],1:519,sep = "_")
res$threshold[res$PValue<0.05&res$logFC>log2(1.5)]="up"
res$threshold[res$PValue<0.05&res$logFC<(-log2(1.5))]="down"
res$threshold[(res$PValue>0.05)|(res$logFC<=log2(1.5))&res$logFC>=(-log2(1.5))]="normal"
res$threshold=factor(res$threshold,levels = c("up","down","normal"),ordered = T)
x_lim=max(res$logFC,-res$logFC)
theme_set(theme_bw())
p <- ggplot(res,aes(logFC,-1*log10(PValue),
                    color =threshold))+geom_point()+
  xlim(-6,6) +  labs(x="log2(FoldChange)",y="-log10(PValue)")
p <- p + scale_color_manual(values =c('up'="red","normal"="grey","down"="blue"))+
  geom_hline(yintercept=-log10(0.05),linetype=4)+
  geom_vline(xintercept=c(-log2(1.5),log2(1.5)),linetype=4)
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

library(clusterProfiler)
library(org.Hs.eg.db)
genelist=mapIds(org.Hs.eg.db,rownames(deg),keytype="SYMBOL",column="ENTREZID")
go=enrichGO(gene=genelist,OrgDb = "org.Hs.eg.db",keyType = "ENTREZID",pvalueCutoff = 0.05,readable = T)
write.csv(go,"~/R/LUAD/Go.csv")
kegg=enrichKEGG(genelist,organism = "hsa",pvalueCutoff = 0.05)
kegg=setReadable(kegg,org.Hs.eg.db,keyType = "ENTREZID")
write.csv(kegg,"~/R/LUAD/Kegg.csv")
setwd("~/R/LUAD")
library(enrichplot)
dotplot(go, showCategory=30)+ ggtitle("dotplot for GO")
dotplot(kegg, showCategory=30)+ ggtitle("dotplot for KEGG")

kegg=read.csv("~/R/LUAD/Kegg.csv")
library(ggplot2)
p <- ggplot(kegg, aes(x=GeneRatio, y=factor(kegg$Description,levels = kegg$Description),  size=Count, color=p.adjust)) + 
  geom_point()  + 
  scale_colour_gradient(low="red",high="blue") + 
  labs(color="p.adjust",size="Count", x="GeneRatio")

p
