library(GEOquery)
gset=getGEO("GSE29265",GSEMatrix = T,AnnotGPL = T)
expSet=exprs(gset[[1]])
pdata=pData(gset[[1]])
fdata=fData(gset[[1]])
library(limma)
library(stringr)
pdata$group=ifelse(grepl(pattern ="non-tumor",pdata$title),"NT","TC")
pdata1=pdata[,c(2,60)]
exp=expSet[,colnames(expSet)%in%rownames(pdata1)]
exp=data.frame(exp)
exp$id=rownames(exp)
fdata=fdata[,c(1,3)]
exp1=merge(exp,fdata,by.x = "id",by.y = "ID")
exp1<-exp1[!(exp1$`Gene symbol`==""), ]
library(dplyr)
exp1<- exp1%>%
  group_by(exp1$`Gene symbol`) %>% 
  summarise_all(max)
table(duplicated(exp1$`Gene symbol`))
rownames(exp1)=exp1$`Gene symbol`
exp1=data.frame(t(exp1))
exp1$ID=rownames(exp1)
exp1=merge(exp1,pdata1,by.x = "ID",by.y = "geo_accession")
exp1=exp1[order(exp1[,22191]),]
Exp=t(exp1)
colnames(Exp)=Exp[22191,]
Exp=matrix(as.numeric(unlist(Exp[-c(1,22191),])),nrow=nrow(Exp[-c(1,22191),]), dimnames = list(rownames(Exp[-c(1,22191),]),colnames(Exp[-c(1,22191),])))
Exp=normalizeBetweenArrays(Exp)
design <- model.matrix(~0+factor(colnames(Exp)))  
colnames(design) <- levels(factor(colnames(Exp))) 
fit <- lmFit(Exp, design)
B.LvsP <- makeContrasts(TC-NT, levels=design)
fit2 <- contrasts.fit(fit, B.LvsP)
fit2 <- eBayes(fit2)
allDiff1=topTable(fit2,adjust='fdr',coef=1,number=Inf)
DEGs1=allDiff1[allDiff1$adj.P.Val<0.05&abs(allDiff1$logFC)>=log2(1.5),]
write.csv(DEGs1,"~/R/TC/GSE29265/DEGs1.csv")
library(ggplot2)
library(RColorBrewer)
allDiff1$threshold[allDiff1$adj.P.Val<0.05&allDiff1$logFC>=log2(1.5)]="up"
allDiff1$threshold[allDiff1$adj.P.Val<0.05&allDiff1$logFC<=-log2(1.5)]="down"
allDiff1$threshold[(allDiff1$adj.P.Val>0.05)|(allDiff1$logFC<log2(1.5))&allDiff1$logFC>-log2(1.5)]="normal"
allDiff1$threshold=factor(allDiff1$threshold,levels = c("up","down","normal"),ordered = T)
x_lim=max(allDiff1$logFC,-allDiff1$logFC)
theme_set(theme_bw())
p <- ggplot(allDiff1,aes(logFC,-1*log10(adj.P.Val),
                     color =threshold ))+geom_point()+
  xlim(-5,5) +  labs(x="log2(FoldChange)",y="-log10(FDR)")
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


library(DT)
library(tableone)

data0=pdata[,c(16,19,32,42)]
names(data0)
vars <- c("characteristics_ch1.6","characteristics_ch1.9","description","contact_country")
catVars <- c("characteristics_ch1.9","description","contact_country")
tab=CreateTableOne(vars=vars,data=data0,factorVars = catVars)
table1=print(tab,noSpaces = TRUE)
write.csv(table1,"table1.csv")


##########################
gset=getGEO("GSE33630",GSEMatrix = T,AnnotGPL = T)
expSet=exprs(gset[[1]])
pdata=pData(gset[[1]])
fdata=fData(gset[[1]])
pdata$group=ifelse(grepl(pattern ="N",pdata$title),"NT","TC")
pdata1=pdata[,c(2,35)]
exp=expSet[,colnames(expSet)%in%rownames(pdata1)]
exp=data.frame(exp)
exp$id=rownames(exp)
fdata=fdata[,c(1,3)]
exp1=merge(exp,fdata,by.x = "id",by.y = "ID")
exp1<-exp1[!(exp1$`Gene symbol`==""), ]
exp1<- exp1%>%
  group_by(exp1$`Gene symbol`) %>% 
  summarise_all(max)
table(duplicated(exp1$`Gene symbol`))
rownames(exp1)=exp1$`Gene symbol`
exp1=data.frame(t(exp1))
exp1$ID=rownames(exp1)
exp1=merge(exp1,pdata1,by.x = "ID",by.y = "geo_accession")
exp1=exp1[order(exp1[,22191]),]
Exp=t(exp1)
colnames(Exp)=Exp[22191,]
Exp=matrix(as.numeric(unlist(Exp[-c(1,22191),])),nrow=nrow(Exp[-c(1,22191),]), dimnames = list(rownames(Exp[-c(1,22191),]),colnames(Exp[-c(1,22191),])))
Exp=normalizeBetweenArrays(Exp)
design <- model.matrix(~0+factor(colnames(Exp)))  
colnames(design) <- levels(factor(colnames(Exp))) 
fit <- lmFit(Exp, design)
B.LvsP <- makeContrasts(TC-NT, levels=design)
fit2 <- contrasts.fit(fit, B.LvsP)
fit2 <- eBayes(fit2)
allDiff2=topTable(fit2,adjust='fdr',coef=1,number=Inf)
DEGs2=allDiff2[allDiff2$adj.P.Val<0.05&abs(allDiff2$logFC)>=log2(1.5),]
write.csv(DEGs2,"~/R/TC/GSE33620/DEGs2.csv")
allDiff2$threshold[allDiff2$adj.P.Val<0.05&allDiff2$logFC>=log2(1.5)]="up"
allDiff2$threshold[allDiff2$adj.P.Val<0.05&allDiff2$logFC<=-log2(1.5)]="down"
allDiff2$threshold[(allDiff2$adj.P.Val>0.05)|(allDiff2$logFC<log2(1.5))&allDiff2$logFC>-log2(1.5)]="normal"
allDiff2$threshold=factor(allDiff2$threshold,levels = c("up","down","normal"),ordered = T)
x_lim=max(allDiff2$logFC,-allDiff2$logFC)
theme_set(theme_bw())
p <- ggplot(allDiff2,aes(logFC,-1*log10(adj.P.Val),
                         color =threshold ))+geom_point()+
  xlim(-5,5) +  labs(x="log2(FoldChange)",y="-log10(FDR)")
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

data0=pdata[,c(19,30)]
names(data0)
vars <- c("description","contact_country")
catVars <- c("description","contact_country")
tab=CreateTableOne(vars=vars,data=data0,factorVars = catVars)
table1=print(tab,noSpaces = TRUE)
write.csv(table1,"~/R/TC/GSE33630/table1.csv")

######################
gset=getGEO("GSE60542",GSEMatrix = T,AnnotGPL = T)
expSet=exprs(gset[[1]])
pdata=pData(gset[[1]])
fdata=fData(gset[[1]])
sample_id=grep(pattern ="Normal|Papillary",pdata$title)
pdata1=pdata[c(sample_id),]
pdata1$group=ifelse(grepl(pattern ="Papillary",pdata1$title),"TC","NT")
pdata2=pdata1[,c(2,65)]
exp=expSet[,colnames(expSet)%in%rownames(pdata2)]
exp=data.frame(exp)
exp$id=rownames(exp)
fdata=fdata[,c(1,3)]
exp1=merge(exp,fdata,by.x = "id",by.y = "ID")
exp1<-exp1[!(exp1$`Gene symbol`==""), ]
exp1<- exp1%>%
  group_by(exp1$`Gene symbol`) %>% 
  summarise_all(max)
rownames(exp1)=exp1$`Gene symbol`
exp1=data.frame(t(exp1))
exp1$ID=rownames(exp1)
exp1=merge(exp1,pdata2,by.x = "ID",by.y = "geo_accession")
exp1=exp1[order(exp1[,22191]),]
Exp=t(exp1)
colnames(Exp)=Exp[22191,]
Exp=matrix(as.numeric(unlist(Exp[-c(1,22191),])),nrow=nrow(Exp[-c(1,22191),]), dimnames = list(rownames(Exp[-c(1,22191),]),colnames(Exp[-c(1,22191),])))
Exp=Exp[,order(colnames(Exp))]
Exp=normalizeBetweenArrays(Exp)
design <- model.matrix(~0+factor(colnames(Exp)))  
colnames(design) <- levels(factor(colnames(Exp))) 
fit <- lmFit(Exp, design)
B.LvsP <- makeContrasts(TC-NT, levels=design)
fit2 <- contrasts.fit(fit, B.LvsP)
fit2 <- eBayes(fit2)
allDiff3=topTable(fit2,adjust='fdr',coef=1,number=Inf)
DEGs3=allDiff3[allDiff3$adj.P.Val<0.05&abs(allDiff3$logFC)>=log2(1.5),]
write.csv(DEGs3,"~/R/TC/GSE60542/DEGs3.csv")
allDiff3$threshold[allDiff3$adj.P.Val<0.05&allDiff3$logFC>=log2(1.5)]="up"
allDiff3$threshold[allDiff3$adj.P.Val<0.05&allDiff3$logFC<=-log2(1.5)]="down"
allDiff3$threshold[(allDiff3$adj.P.Val>0.05)|(allDiff3$logFC<log2(1.5))&allDiff3$logFC>-log2(1.5)]="normal"
allDiff3$threshold=factor(allDiff3$threshold,levels = c("up","down","normal"),ordered = T)
x_lim=max(allDiff3$logFC,-allDiff3$logFC)
theme_set(theme_bw())
p <- ggplot(allDiff3,aes(logFC,-1*log10(adj.P.Val),
                         color =threshold ))+geom_point()+
  xlim(-5,5) +  labs(x="log2(FoldChange)",y="-log10(FDR)")
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

data0=pdata1[,c(11,14,45,65)]
names(data0)
vars <- c("characteristics_ch1.1","characteristics_ch1.4","contact_country","group")
catVars <- c("characteristics_ch1.1","characteristics_ch1.4","contact_country","group")
tab=CreateTableOne(vars=vars,data=data0,factorVars = catVars)
table1=print(tab,noSpaces = TRUE)
write.csv(table1,"~/R/TC/GSE60542/table1.csv")

#######################
gset=getGEO("GSE65144",GSEMatrix = T,AnnotGPL = T)
expSet=exprs(gset[[1]])
pdata=pData(gset[[1]])
fdata=fData(gset[[1]])
pdata$group=ifelse(grepl(pattern ="N",pdata$title),"NT","TC")
pdata1=pdata[,c(2,37)]
exp=expSet[,colnames(expSet)%in%rownames(pdata1)]
exp=data.frame(exp)
exp$id=rownames(exp)
fdata=fdata[,c(1,3)]
exp1=merge(exp,fdata,by.x = "id",by.y = "ID")
exp1<-exp1[!(exp1$`Gene symbol`==""), ]
exp1[,c(2:26)]=log2(exp1[,c(2:26)]+0.001)
exp1<- exp1%>%
  group_by(exp1$`Gene symbol`) %>% 
  summarise_all(max)
table(duplicated(exp1$`Gene symbol`))
rownames(exp1)=exp1$`Gene symbol`
exp1=data.frame(t(exp1))
exp1$ID=rownames(exp1)
exp1=merge(exp1,pdata1,by.x = "ID",by.y = "geo_accession")
exp1=exp1[order(exp1[,22191]),]
Exp=t(exp1)
colnames(Exp)=Exp[22191,]
Exp=matrix(as.numeric(unlist(Exp[-c(1,22191),])),nrow=nrow(Exp[-c(1,22191),]), dimnames = list(rownames(Exp[-c(1,22191),]),colnames(Exp[-c(1,22191),])))
Exp=normalizeBetweenArrays(Exp)
design <- model.matrix(~0+factor(colnames(Exp)))  
colnames(design) <- levels(factor(colnames(Exp))) 
fit <- lmFit(Exp, design)
B.LvsP <- makeContrasts(TC-NT, levels=design)
fit2 <- contrasts.fit(fit, B.LvsP)
fit2 <- eBayes(fit2)
allDiff4=topTable(fit2,adjust='fdr',coef=1,number=Inf)
DEGs4=allDiff4[allDiff4$adj.P.Val<0.05&abs(allDiff4$logFC)>=log2(1.5),]
write.csv(DEGs4,"~/R/TC/GSE65144/DEGs4.csv")
allDiff4$threshold[allDiff4$adj.P.Val<0.05&allDiff4$logFC>=log2(1.5)]="up"
allDiff4$threshold[allDiff4$adj.P.Val<0.05&allDiff4$logFC<=-log2(1.5)]="down"
allDiff4$threshold[(allDiff4$adj.P.Val>0.05)|(allDiff4$logFC<log2(1.5))&allDiff4$logFC>-log2(1.5)]="normal"
allDiff4$threshold=factor(allDiff4$threshold,levels = c("up","down","normal"),ordered = T)
x_lim=max(allDiff4$logFC,-allDiff4$logFC)
theme_set(theme_bw())
p <- ggplot(allDiff4,aes(logFC,-1*log10(adj.P.Val),
                         color =threshold ))+geom_point()+
  xlim(-8,8) +  labs(x="log2(FoldChange)",y="-log10(FDR)")
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

data0=pdata[,c(8,31)]
names(data0)
vars <- c("source_name_ch1","contact_country")
catVars <- c("source_name_ch1","contact_country")
tab=CreateTableOne(vars=vars,data=data0,factorVars = catVars)
table1=print(tab,noSpaces = TRUE)
write.csv(table1,"~/R/TC/GSE65144/table1.csv")


######################
up=read.csv("~/R/TC/up.csv",na="")
down=read.csv("~/R/TC/down.csv",na="")
library(RobustRankAggreg)
up<- as.list(up)
for (i in 1: length(up)){up[[i]]<- up[[i]][!is.na(up[[i]])]
up[[i]]<- as.character(up[[i]])}
count<- c(907, 1608, 1148,2598)
degUp=aggregateRanks(up,
               rmat = rankMatrix(up, N=count),
               method = "RRA")
degUp=degUp[degUp$Score<0.01,]
write.csv(degUp,"~/R/TC/degUp.csv")
down<- as.list(down)
for (i in 1: length(down)){down[[i]]<- down[[i]][!is.na(down[[i]])]
down[[i]]<- as.character(down[[i]])}
count<- c(1049, 1507, 962,2937)
degDown=aggregateRanks(down,
               rmat = rankMatrix(down, N=count), 
               method = "RRA")
degDown=degDown[degDown$Score<0.01,]
write.csv(degDown,"~/R/TC/degDown.csv")
deg=rbind(degUp,degDown)
GSE29265_diff$ID=rownames(GSE29265_diff)
GSE33630_diff$ID=rownames(GSE33630_diff)
GSE60542_diff$ID=rownames(GSE60542_diff)
GSE65144_diff$ID=rownames(GSE65144_diff)
logFC_1=merge(deg,GSE29265_diff,by.x="Name",by.y="ID")
logFC_2=merge(deg,GSE33630_diff,by.x="Name",by.y="ID")
logFC_3=merge(deg,GSE60542_diff,by.x="Name",by.y="ID")
logFC_4=merge(deg,GSE65144_diff,by.x="Name",by.y="ID")
logFC=merge(logFC_1[,c(1,3)],logFC_2[,c(1,3)],by.x="Name",by.y="Name")
logFC=merge(logFC,logFC_3[,c(1,3)],by.x="Name",by.y="Name")
colnames(logFC)=c("Name","GSE29265","GSE33630","GSE60542")
logFC=merge(logFC,logFC_4[,c(1,3)],by.x="Name",by.y="Name")
colnames(logFC)=c("Name","GSE29265","GSE33630","GSE60542","GSE65144")
logFC=logFC[order(-logFC$GSE29265),]
logFC[logFC$GSE29265<0,]=logFC[logFC$GSE29265<0,][order(logFC[logFC$GSE29265<0,]$GSE29265),]
rownames(logFC)=logFC[,1]
logFC=logFC[,-1]
write.csv(logFC,"~/R/TC/logFC.csv")
logFC=read.csv("~/R/TC/logFC.csv")
rownames(logFC)=logFC[,2]
logFC=logFC[,-c(1,2)]
deg=deg[rownames(deg)%in%rownames(logFC),]
write.csv(deg,"~/R/TC/deg.csv")
pheatmap(logFC, 
         cluster_rows = F,
         cluster_cols = F,
         annotation_legend=TRUE, 
         color =colorRampPalette(c("blue", "white","red"))(100),
         show_rownames = T,
         show_colnames = T,
         display_numbers = T,
         number_color = "black",
         cellwidth = 50, cellheight = 10,
         fontsize = 10)



###################
library(clusterProfiler)
library(org.Hs.eg.db)
list=c("FN1","TIMP1","ITGA2","MRC2","KIT")
genelist=mapIds(org.Hs.eg.db,list,keytype="SYMBOL",column="ENTREZID")
go=enrichGO(gene=genelist,OrgDb = "org.Hs.eg.db",keyType = "ENTREZID",pvalueCutoff = 0.05,readable = T)
write.csv(go,"~/R/TC/HUB GENE/Go.csv")
kegg=enrichKEGG(genelist,organism = "hsa",pvalueCutoff = 0.05)
kegg=setReadable(kegg,org.Hs.eg.db,keyType = "ENTREZID")
write.csv(kegg,"~/R/TC/HUB GENE/Kegg.csv")
library(enrichplot)
dotplot(go, showCategory=30)+ ggtitle("dotplot for GO")
dotplot(kegg, showCategory=30)+ ggtitle("dotplot for KEGG")


