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
Exp = exp.hg38.values[,substring(colnames(exp.hg38.values),14,16)=="01A"|substring(colnames(exp.hg38.values),14,16)=="11A"]
clinical <- GDCquery_clinic(project = "TCGA-LUAD", type = "clinical")
samples=clinical[,c("submitter_id","tumor_stage","days_to_last_follow_up","gender","age_at_index","vital_status","days_to_death","race")]
samples=samples[!(is.na(samples$days_to_last_follow_up)&is.na(samples$days_to_death)),]
samples[samples$vital_status=="Dead","survival_time"]=samples[samples$vital_status=="Dead","days_to_death"]
samples[samples$vital_status=="Alive","survival_time"]=samples[samples$vital_status=="Alive","days_to_last_follow_up"]
samples$status=ifelse(samples$vital_status=="Dead",1,0)
sign =c()
patient_id = colnames(Exp)
for (i in 1:length(colnames(Exp))){
  tmp = (strsplit2(as.character(patient_id[i]),split = "-"))
  sign[i] = tmp[4]
  tmp = paste(tmp[1],tmp[2],tmp[3],sep = "-")
  patient_id[i] = tmp
}
colnames(Exp)=patient_id
ExpCan = Exp[,sign == "01A"]
ExpNor = Exp[,sign == "11A"]
SAMPID=colnames(ExpCan)
ExpCan=rbind(ExpCan,SAMPID)
ExpCan=data.frame(t(ExpCan))
ExpCanCli=merge(ExpCan,samples,by.x="SAMPID",by.y="submitter_id")

library(DT)
library(tableone)
Cli=ExpCanCli[,c(1,56495:56502)]
names(Cli)
vars <- c("tumor_stage","gender","age_at_index","race")
catVars <- c("tumor_stage","gender","race")
tab=CreateTableOne(vars=vars,data=Cli,factorVars = catVars)
table1=print(tab,noSpaces = TRUE)
write.csv(table1,"~/R/HCG18/table1.csv")

ExpNor =data.frame(t(ExpNor))
CanNor=rbind(ExpNor,ExpCanCli[,c(2:56494)])
rownames = rownames(CanNor)
CanNor = apply(CanNor,2,as.numeric)
CanNor = as.data.frame(CanNor)
CanNor = cbind(rownames,CanNor)
colnames(CanNor)[1] <- "SAMPID"
CanNor$group=factor(rep(c("normal","cancer"), c(58,510)))
library(ggstatsplot)
ggbetweenstats(data=CanNor[,c("HCG18","group")], 
               x = group, 
               y = HCG18) 
library(survival)
library(rms)
library(Hmisc)
library(lattice)
library(Formula)
library(ggplot2)
library(grid)
ExpCanCli[,c(2:56494)] = apply(ExpCanCli[,c(2:56494)],2,as.numeric)
ExpCanCli = as.data.frame(ExpCanCli)
ExpCanCli$ExpHCG18 = ifelse(ExpCanCli$HCG18>median(ExpCanCli$HCG18),"high","low")
fit = survival::survdiff(Surv(survival_time,status)~ExpHCG18,data=ExpCanCli)
fit1 <- coxph(Surv(survival_time,status)~ExpHCG18, data=ExpCanCli)
summary(fit1)
library(survminer)
ggsurvplot(survfit(Surv(survival_time,status)~ExpHCG18,data=ExpCanCli),data=ExpCanCli,pval = T,linetype = c('solid', 'dashed'),palette=c("red","black"))
ExpCanCli[substring(ExpCanCli$tumor_stage,7,9)=="iv","tumor_stage"]="stage IV"
ExpCanCli[substring(ExpCanCli$tumor_stage,7,9)=="iii","tumor_stage"]="stage III"
ExpCanCli[substring(ExpCanCli$tumor_stage,7,9)=="ii"|substring(ExpCanCli$tumor_stage,7,9)=="iia"|substring(ExpCanCli$tumor_stage,7,9)=="iib","tumor_stage"]="stage II"
ExpCanCli[substring(ExpCanCli$tumor_stage,7,9)=="i"|substring(ExpCanCli$tumor_stage,7,9)=="ia"|substring(ExpCanCli$tumor_stage,7,9)=="ib","tumor_stage"]="stage I"
colnames(ExpCanCli)[56495] <- "stage"
colnames(ExpCanCli)[56498] <- "age"
fit2 <- coxph(Surv(survival_time,status)~ExpHCG18+stage+gender+age, data=ExpCanCli)
summary(fit2)
ggforest(fit2, 
         data = ExpCanCli, 
         main = 'Hazard ratio of LUAD', 
         cpositions = c(0.05, 0.15, 0.35), 
         fontsize = 1,
         refLabel = 'reference', 
         noDigits = 3)

library(edgeR)
library(limma)
ExpCanCli=t(ExpCanCli)
ExpM=matrix(as.numeric(unlist(ExpCanCli[2:56494,])),nrow=nrow(ExpCanCli[2:56494,]), dimnames = list(rownames(ExpCanCli[2:56494,]),colnames(ExpCanCli[2:56494,])))
y = DGEList(ExpM,group = ExpCanCli[56504,],genes = ExpM)
keep<-rowSums(cpm(y) > 0.5) >= 2
y<-y[keep, ,keep.lib.sizes=FALSE]
y=calcNormFactors(y)
t2 = edgeR::estimateCommonDisp(y)
t3 = edgeR::exactTest(t2,pair=c("low","high"))
tableDEA <- edgeR::topTags(t3, n = nrow(t3$table))$table
LipidDeg=tableDEA[tableDEA$FDR<0.05&abs(tableDEA$logFC)>1,]
write.csv(LipidDeg[,c(578:581)],"~/Downloads/合作脂代�?/LipidDeg.csv")
library(ggplot2)
library(RColorBrewer)
tableDEA$threshold[tableDEA$FDR<0.05&tableDEA$logFC>1]="up"
tableDEA$threshold[tableDEA$FDR<0.05&tableDEA$logFC<(-1)]="down"
tableDEA$threshold[(tableDEA$FDR>0.05)|(tableDEA$logFC<=1)&tableDEA$logFC>=(-1)]="normal"
tableDEA$threshold=factor(tableDEA$threshold,levels = c("up","down","normal"),ordered = T)
x_lim=max(tableDEA$logFC,-tableDEA$logFC)
theme_set(theme_bw())
p <- ggplot(tableDEA,aes(logFC,-1*log10(FDR),
                         color =threshold ))+geom_point()+
  xlim(-8,8) +  labs(x="log2(FoldChange)",y="-log10(FDR)")
p <- p + scale_color_manual(values =c('up'="red","normal"="grey","down"="blue"))+
  geom_hline(yintercept=-log10(0.05),linetype=4)+
  geom_vline(xintercept=c(-log2(2),log2(2)),linetype=4)
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
library(dplyr)
genelist=mapIds(org.Hs.eg.db,rownames(LipidDeg),keytype="SYMBOL",column="ENTREZID")
go=enrichGO(gene=genelist,OrgDb = "org.Hs.eg.db",keyType = "ENTREZID",pvalueCutoff = 0.05,readable = T,ont="ALL")
write.csv(go,"~/Downloads/合作脂代�?/Go.csv")
kegg=enrichKEGG(genelist,organism = "hsa",pvalueCutoff = 0.05)
kegg=setReadable(kegg,org.Hs.eg.db,keyType = "ENTREZID")
write.csv(kegg,"~/Downloads/合作脂代�?/Kegg.csv")
library(enrichplot)
dotplot(go, showCategory=30)+ ggtitle("dotplot for GO")
dotplot(kegg, showCategory=30)+ ggtitle("dotplot for KEGG")
tableDEA$ENTREZID=mapIds(org.Hs.eg.db,rownames(tableDEA),keytype="SYMBOL",column="ENTREZID")
GeneGSEA=tableDEA[order(tableDEA$logFC,decreasing=T),]$ENTREZID
gmtfile <- system.file("extdata", "c5.cc.v5.0.entrez.gmt", package="clusterProfiler")
c5 <- read.gmt(gmtfile)
enrich <- enricher(GeneGSEA, TERM2GENE=c5)
glist <- tableDEA$logFC
names(glist) <- as.character(tableDEA$ENTREZID)
glist <- sort(glist,decreasing = T)
gsea <- GSEA(glist, TERM2GENE=c5, verbose=FALSE, pvalueCutoff = 0.2)
write.csv(gsea,"~/Downloads/合作脂代�?/GSEA.csv")
