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
table(substring(colnames(exp.hg38.values),14,16))
Exp = exp.hg38.values[,substring(colnames(exp.hg38.values),14,16)=="01A"|substring(colnames(exp.hg38.values),14,16)=="11A"]
group <- ifelse(substring(colnames(Exp),14,16)=="01A","cancer","normal")
LipidGene=read.csv("E:/??Ŀ/֬??л/GENE.csv",header = F)
LipidExp=Exp[rownames(Exp)%in%LipidGene$V1,]
library(edgeR)
library(limma)
y = DGEList(LipidExp,group = group,genes = LipidExp)
keep<-rowSums(cpm(y) > 0.5) >= 2
y<-y[keep, ,keep.lib.sizes=FALSE]
y=calcNormFactors(y)
t2 = edgeR::estimateCommonDisp(y)
t3 = edgeR::exactTest(t2,pair=c("normal","cancer"))
tableDEA <- edgeR::topTags(t3, n = nrow(t3$table))$table
LipidDeg=tableDEA[tableDEA$FDR<0.05&abs(tableDEA$logFC)>1,]
write.csv(LipidDeg[,c(578:581)],"~/Downloads/合作脂代???/LipidDeg.csv")
library(ggplot2)
library(RColorBrewer)
tableDEA$threshold[tableDEA$FDR<0.05&tableDEA$logFC>1]="up"
tableDEA$threshold[tableDEA$FDR<0.05&tableDEA$logFC<(-1)]="down"
tableDEA$threshold[(tableDEA$FDR>0.05)|(tableDEA$logFC<=1)&tableDEA$logFC>=(-1)]="normal"
tableDEA$threshold=factor(tableDEA$threshold,levels = c("up","down","normal"),ordered = T)
p <- ggplot(data = tableDEA, 
            aes(x = logFC, 
                y = -log10(FDR))) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=change)) +
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  theme_bw()
p

library(pheatmap)
heatdata <- LipidDeg[,-c(578:581)]
group <- ifelse(substring(colnames(heatdata),14,16)=="01A","cancer","normal")
heatdata=rbind(heatdata,group)
heatdata=heatdata[,order(heatdata[218,])]
annotation_col = data.frame(Type = factor(rep(c("cancer","normal"), c(519,58))))
rownames(annotation_col) = colnames(heatdata)
heatdata=matrix(as.numeric(unlist(heatdata[-c(218),])),nrow=nrow(heatdata[-c(218),]), dimnames = list(rownames(heatdata[-c(218),]),colnames(heatdata[-c(218),])))
pheatmap(heatdata, 
         cluster_rows = T,
         cluster_cols = F,
         annotation_col = annotation_col,
         show_rownames = F,
         show_colnames = F,
         scale = "row",
         color =colorRampPalette(c("blue", "white","red"))(100),
         cellwidth = 0.5, cellheight = 1,
         fontsize = 10)


library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(GOplot)
genelist=mapIds(org.Hs.eg.db,rownames(LipidDeg),keytype="SYMBOL",column="ENTREZID")
go=enrichGO(gene=genelist,OrgDb = "org.Hs.eg.db",keyType = "ENTREZID",pvalueCutoff = 0.05,readable = T,ont="ALL")
write.csv(go,"~/Downloads/合作脂代???/Go.csv")
kegg=enrichKEGG(genelist,organism = "hsa",pvalueCutoff = 0.05)
kegg=setReadable(kegg,org.Hs.eg.db,keyType = "ENTREZID")
go1=data.frame(go)
library(stringr)
go2=go1[,c(1:3,8:9)]
go2$geneID=str_replace_all(go2$geneID,"/",",")
names(go2)=c("Category","ID","Term","adj_pval","Genes")
gene=data.frame(genelist)
gene$SYMBOL=rownames(gene)
colnames(gene)[1]=c("ENTREZID")
gene$logFC=LipidDeg[gene$SYMBOL,]$logFC
gene=gene[,c(2,3)]
names(gene)[1]="ID"
plot_data=list(DA=go2,GE=gene)
circ2=data.frame()
circ2=circle_dat(plot_data$DA,plot_data$GE)
chord=chord_dat(data=circ2,genes = gene,process = unique(circ2$term))
GOBubble(circ2,labels = 30)
kegg1=data.frame(kegg)
kegg1=kegg1[,c(1,2,7,8)]
kegg1$geneID=str_replace_all(kegg1$geneID,"/",",")
names(kegg1)=c("ID","Term","adj_pval","Genes")
kegg1$Category=c("KEGG pathway")
plot_data=list(DA=kegg1,GE=gene)
circ2=data.frame()
circ2=circle_dat(plot_data$DA,plot_data$GE)
GOBubble(circ2,labels = 10)

write.csv(kegg,"~/Downloads/合作脂代???/Kegg.csv")
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
write.csv(gsea,"~/Downloads/合作脂代???/GSEA.csv")
library(fgsea)
library(data.table)
library(ggplot2)
pathways=read.csv("~/Downloads/合作脂代???/pathways.csv")
list_data <- list(c(as.character(pathways[1,c(2:18)])), c(as.character(pathways[2,c(2:9)])), c(as.character(pathways[3,c(2:16)])),c(as.character(pathways[4,c(2:35)])),c(as.character(pathways[5,c(2:13)])),
                  c(as.character(pathways[6,c(2:13)])), c(as.character(pathways[7,c(2:8)])), c(as.character(pathways[8,c(2:8)])),c(as.character(pathways[9,c(2:18)])),c(as.character(pathways[10,c(2:16)])))
names(list_data)=c(pathways[,1])
list_data[1]
set.seed(42)
fgseaRes <- fgsea(pathways = list_data, 
                  stats = glist,
                  minSize=15,
                  maxSize=500)
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(list_data[topPathways], glist, fgseaRes, gseaParam=0.5)

clinical <- GDCquery_clinic(project = "TCGA-LUAD", type = "clinical")
samples=clinical[,c("submitter_id","tumor_stage","days_to_last_follow_up","gender",
                    "age_at_index","vital_status","days_to_death","pack_years_smoked")]
samples=samples[!(is.na(samples$days_to_last_follow_up)&is.na(samples$days_to_death)),]
samples[samples$vital_status=="Dead","survival_time"]=samples[samples$vital_status=="Dead","days_to_death"]
samples[samples$vital_status=="Alive","survival_time"]=samples[samples$vital_status=="Alive","days_to_last_follow_up"]

tableDEA1=tableDEA[,c(1:577)]
sign =c()
patient_id = colnames(tableDEA1)
for (i in 1:length(colnames(tableDEA1))){
  tmp = (strsplit2(as.character(patient_id[i]),split = "-"))
  sign[i] = tmp[4]
  tmp = paste(tmp[1],tmp[2],tmp[3],sep = "-")
  patient_id[i] = tmp
}
colnames(tableDEA1)=patient_id
ExpCancer = tableDEA1[,sign == "01A"]
hub=ExpCancer[c("CYP2C9","DGAT1","HPGDS","INS","LPL","UGT1A6"),]
sample_id=colnames(hub)
hubgene=rbind(hub,sample_id)
hubgene=t(hubgene)
hubCli=merge(hubgene,samples,by.x="7",by.y="submitter_id")
hubCli$status=ifelse(hubCli$vital_status=="Alive",0,1)
hubCli=hubCli [which(hubCli$survival_time > 0),]
hubCli$CYP2C9=as.numeric(as.character(hubCli$CYP2C9))+0.000000001
hubCli$DGAT1=as.numeric(as.character(hubCli$DGAT1))+0.000000001
hubCli$HPGDS=as.numeric(as.character(hubCli$HPGDS))+0.000000001
hubCli$INS=as.numeric(as.character(hubCli$INS))+0.000000001
hubCli$LPL=as.numeric(as.character(hubCli$LPL))+0.000000001
hubCli$UGT1A6=as.numeric(as.character(hubCli$UGT1A6))+0.000000001

library(survival)
library(rms)
library(Hmisc)
library(lattice)
library(Formula)
library(ggplot2)
library(grid)
dd<-datadist(hubCli)
options(datadist="dd")
f=cph(Surv(survival_time,status)~CYP2C9+DGAT1+HPGDS+INS+LPL+UGT1A6,data = hubCli,
      x=T,y=T,surv=T)
survival=Survival(f)
nom=nomogram(f,fun=list(function(x) survival(1095, x),
                        function(x) survival(1825, x)),
          funlabel=c("3-year Overall Survival",
                     "5-year Overall Survival"),lp=F)
plot(nom,xfrac=.6)
library(nomogramEx)
riskFormul = nomogramEx(nomo=nom,np=2,digit=9)
riskScore=c()
for(i in 1:nrow(hubCli)){
  s1 = -2.272048277 * hubCli[i,2] + 22.720482766
  s2 = 0.118293915 * hubCli[i,3] + 0
  s3 = 0 * hubCli[i,4] ^3 + 0 * hubCli[i,4] ^2+ -7.142857143 * hubCli[i,4] + 100
  s4 = 0 * hubCli[i,5] ^3 + 0 * hubCli[i,5] ^2 + 0.311020041 * hubCli[i,5] + 0
  s5 = -0.167871493 * hubCli[i,6] + 33.574298593
  s6 = 0 * hubCli[i,7] ^2 + -0.260190391 * hubCli[i,7] + 16.912375432
  riskScore[i] = s1+s2+s3+s4+s5+s6
}
hubCli$riskScoreGroup = ifelse(riskScore>median(riskScore),"High-risk score","Low-risk score")
hubCli$status[which(hubCli$survival_time>900)] <- 0
hubCli$survival_time[which(hubCli$survival_time>900)] <- 900
hubCli$gender=ifelse(hubCli$gender=="female","Female","Male")
fit = survival::survdiff(Surv(survival_time,status)~riskScoreGroup,data=hubCli)
fit1 <- coxph(Surv(survival_time,status)~riskScoreGroup, data=hubCli)
summary(fit1)
library(survminer)
ggsurvplot(survfit(Surv(survival_time,status)~riskScoreGroup,data=hubCli),data=hubCli,pval = T,linetype = c('solid', 'dashed'),palette=c("black","red"))
hubCli=hubCli[!hubCli$tumor_stage=="not reported",]
hubCli[substring(hubCli$tumor_stage,7,9)=="iv","stage"]="Advanced stage"
hubCli[substring(hubCli$tumor_stage,7,9)=="iii","stage"]="Advanced stage"
hubCli[substring(hubCli$tumor_stage,7,9)=="ii"|substring(hubCli$tumor_stage,7,9)=="iia"|
         substring(hubCli$tumor_stage,7,9)=="iib","stage"]="Early stage"
hubCli[substring(hubCli$tumor_stage,7,9)=="i"|substring(hubCli$tumor_stage,7,9)=="ia"|
         substring(hubCli$tumor_stage,7,9)=="ib","stage"]="Early stage"
summary(is.na(hubCli$pack_years_smoked))
library(tableone)
stable1 <- CreateTableOne(vars=c("gender","age_at_index","riskScoreGroup",
                                 "stage","vital_status","pack_years_smoked"),
                          data=hubCli,factorVars=c("gender","stage","vital_status"))
stable1=print(stable1,showAllLevels = TRUE)

model <- coxph( Surv(survival_time,status) ~  riskScoreGroup+ gender + stage+
                  age_at_index , data =  hubCli )
ggforest(model,  #coxph得到的Cox回归结果
         data = hubCli,  #数据集
         cpositions = c(0.05, 0.15, 0.35),  #前三列距离
         fontsize = 1.2, #字体大小
         refLabel = 'reference', #相对变量的数值标签，也可改为1
         noDigits = 2 #保留HR值以及95%CI的小数位数
)

fp<-predict(f)
cindex=1-rcorr.cens(fp,Surv(as.numeric(hubCli$survival_time),hubCli$status))[[1]]
cindex
cal2<- calibrate(f,cmethod='KM',method='boot',u=3*365,m=166,B=1000)
plot(cal2,lwd = 2,lty = 1,errbar.col = c("#00468BFF"),xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced 3-Year OS (%)",ylab = "Observed OS (%)",
     col = c("#B2182B"),cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal2[,c('mean.predicted',"KM")],type= 'b',z = 2,col = c("#00468BFF"),pch = 16)
mtext("")
box(lwd = 1)
abline(0,1,lty= 3,lwd = 2,col =c("#224444"))

cal3<- calibrate(f,cmethod='KM',method='boot',u=5*365,m=166,B=1000)
plot(cal3,lwd = 2,lty = 1,errbar.col = c("#00468BFF"),xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced 5-Year OS (%)",ylab = "Observed OS (%)",
     col = c("#B2182B"),cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal3[,c('mean.predicted',"KM")],type= 'b',lwd = 2,col = c("#00468BFF"),pch = 16)
mtext("")
box(lwd = 1)
abline(0,1,lty= 3,lwd = 2,col =c("#224444"))

gset=getGEO("GSE13213",GSEMatrix = T,AnnotGPL = T)
expSet=exprs(gset[[1]])
pdata=pData(gset[[1]])
fdata=fData(gset[[1]])
symbol=fdata[,c(1,3)]
expSet=data.frame(expSet)
expSet$ID=rownames(expSet)
expSet=merge(expSet,symbol,by.x = "ID",by.y = "ID")
boxplot(expSet[,c(2:118)],col="blue")
expSet[,c(2:118)]=normalizeBetweenArrays(expSet[,c(2:118)])
boxplot(expSet[,c(2:118)],col="blue")
expSet<-expSet[!(expSet$`Gene symbol`==""), ]
expSet<- expSet%>%
  group_by(expSet$`Gene symbol`) %>% 
  summarise_all(max)
expSet=data.frame(expSet)
rownames(expSet)=expSet$expSet..Gene.symbol.
hubGEO=expSet[c("CYP2C9","DGAT1","HPGDS","INS","LPL","UGT1A6"),]
sample_id=colnames(hubGEO)
hubGEO=rbind(hubGEO,sample_id)
hubGEO=data.frame(t(hubGEO))
names(pdata)
samples=pdata[,c("geo_accession","Survival (days):ch1","Status:ch1","Stage (Pathological ):ch1","Smoking (BI):ch1","Sex:ch1","Age:ch1","Histology:ch1")]
samples=samples[!(is.na(samples$`Survival (days):ch1`)),]
hubCli=merge(hubGEO,samples,by.x="X7",by.y="geo_accession")
hubCli$CYP2C9=as.numeric(as.character(hubCli$CYP2C9))
hubCli$DGAT1=as.numeric(as.character(hubCli$DGAT1))
hubCli$HPGDS=as.numeric(as.character(hubCli$HPGDS))
hubCli$INS=as.numeric(as.character(hubCli$INS))
hubCli$LPL=as.numeric(as.character(hubCli$LPL))
hubCli$UGT1A6=as.numeric(as.character(hubCli$UGT1A6))
hubCli$`Survival (days):ch1`=as.numeric(as.character(hubCli$`Survival (days):ch1`))
hubCli$status=ifelse(hubCli$`Status:ch1`=="Alive",0,1)
riskScore=c()
for(i in 1:nrow(hubCli)){
  s1 = -2.272048277 * hubCli[i,2] + 22.720482766
  s2 = 0.118293915 * hubCli[i,3] + 0
  s3 = 0 * hubCli[i,4] ^3 + 0 * hubCli[i,4] ^2+ -7.142857143 * hubCli[i,4] + 100
  s4 = 0 * hubCli[i,5] ^3 + 0 * hubCli[i,5] ^2 + 0.311020041 * hubCli[i,5] + 0
  s5 = -0.167871493 * hubCli[i,6] + 33.574298593
  s6 = 0 * hubCli[i,7] ^2 + -0.260190391 * hubCli[i,7] + 16.912375432
  riskScore[i] = s1+s2+s3+s4+s5+s6
}
hubCli$riskScoreGroup = ifelse(riskScore>median(riskScore),1,0)
hubCli$status[which(hubCli$`Survival (days):ch1`>2100)] <- 0
hubCli$`Survival (days):ch1`[which(hubCli$`Survival (days):ch1`>2100)] <- 2100
fit = survival::survdiff(Surv(`Survival (days):ch1`,status)~riskScoreGroup,data=hubCli)
fit1 <- coxph(Surv(`Survival (days):ch1`,status)~riskScoreGroup, data=hubCli)
summary(fit1)
library(survminer)
ggsurvplot(survfit(Surv(`Survival (days):ch1`,status)~riskScoreGroup,data=hubCli),data=hubCli,pval = T,linetype = c('solid', 'dashed'),palette=c("black","red"))
