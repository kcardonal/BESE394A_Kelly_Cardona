#DEA PUDP OVEREXPRESSING CELL LINES

#Load libaries
library(GenomicFeatures)
library(tidyverse)
library(openxlsx)
library(readxl)
library(HTSFilter)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(RSNNS)
library(lme4)
library(biomaRt)

###___________________0.1 DATA PREPARATION_______________________________________________________________________________
#Construct the counts matrix
file_list <- list.files("./Counts_PUDP/", pattern = "*.txt", full.names = TRUE)
#Read the data from each file into a list
counts_list <- lapply(file_list, read.table, header=TRUE)

#Create a dataframe by merging all the individual dataframes, then assign the sample names as the column names
#and the first column Geneid as row names

counts <- counts_list %>% reduce(inner_join, by='Geneid')
rownames(counts) <- counts$Geneid
#Keep only the columns containg the actual counting
counts=counts[,c(7,13,19,25,31,37,43,49,55)]
#Counts removing 2 samples: L_5816 & L_5851
#counts_7s=counts[,c(7,13,25,31,37,49,55)]
colnames(counts) <- file_list
colnames(counts) <- gsub(pattern = "./Counts_PUDP//", replacement = "", colnames(counts))
colnames(counts) <- gsub(pattern = "_gene_assigned.txt", replacement = "", colnames(counts))

#Substitute the NA with 0
counts[is.na(counts)] <- 0

#Save the counts matrix for further analysis (or just in case)
write.xlsx(counts, file = "rawcounts_PUDP_exc.xlsx", rowNames=TRUE, colNames=TRUE)

#Upload Metadata
metadataPUDP <- as.data.frame(read_excel("metadata_PUDP.xlsx"))
#metadataPUDP2 <-metadataPUDP[-c(3,7),]

###___________________0.1 DATA NORMALIZATION & PCA_______________________________________________________________________________

# combine time and karyotype factors to create a single condition for the filter on HTSFilter
metadataPUDP$Group=as.factor(metadataPUDP$Group)
#metadataPUDP2$Group=as.factor(metadataPUDP2$Group)

filtered<- HTSFilter(counts, metadataPUDP$Group, s.len=50, normalization="TMM", plot=FALSE)$filteredData
#filtered2 <-filtered[,-c(3,7)]

analysisPUDP <- DGEList(counts=filtered, group=metadataPUDP$Group)
analysisPUDP <- calcNormFactors(analysisPUDP, method = "TMM")

filter_TMM <- as.data.frame(normalizeData(filtered, "TMM"))
filter_PCA <- prcomp(t(filter_TMM))

#eigenvalues
eig <- (filter_PCA$sdev)^2
#variances in percentage
variance <- round(eig*100/sum(eig), digits=2)

filter_PCA <- as.data.frame(filter_PCA$x)
ss <- c(filter_PCA$PC1,filter_PCA$PC2)

#PCA in general
pb1=ggplot(filter_PCA, aes(x = PC1, y = PC2, color=metadataPUDP2$SampleID, clone=metadataPUDP2$Clone )) +
  geom_point(size = 3) +
  # xlim(c(min(ss), max(ss))) + ylim(c(min(ss), max(ss))) +
  xlab(paste("PCA1 (",variance[1],"%)",sep="")) +
  ylab(paste("PCA2 (",variance[2],"%)",sep=""))

pb2=ggplot(filter_PCA, aes(x = PC1, y = PC2, color=metadataPUDP2$Group, clone=metadataPUDP2$Clone )) +
  geom_point(size = 3) +
  # xlim(c(min(ss), max(ss))) + ylim(c(min(ss), max(ss))) +
  xlab(paste("PCA1 (",variance[1],"%)",sep="")) +
  ylab(paste("PCA2 (",variance[2],"%)",sep=""))

pb1 + pb2

#Save the normalized counts matrix for further analysis 

filtcounts <- as.data.frame(analysisPUDP[["counts"]])
write.xlsx(filtcounts, file = "filter_counts_PUDP.xlsx", rowNames=TRUE, colNames=TRUE)

#_______________________DE ANALYSIS_____________________________________________
# counts is the matrix of raw counts
# metadata: rownames are the same of colnames of counts.
#the 
metadataPUDP$Group=as.factor(metadataPUDP$Group)
metadataPUDP$Clone=as.factor(metadataPUDP$Clone)

#comparison factor Group 
design <- model.matrix(~ Group, metadataPUDP)

#Estimate the residual variance
analysisPUDP <- DGEList(counts=filtered, group=metadataPUDP$Group)
analysisPUDP <- estimateDisp(analysisPUDP, design, robust = T)
fit <- glmFit(analysisPUDP, design)

c <- ncol(design)
# coef should be one of  the colnames(design), otherwise contrast:

#the glmLRT() function to perform likelihood ratio tests to compare the expression 
#of genes between the different karyotypes at each time point. 
#The coef argument specifies the contrast that you want to test.

#OE_PUDP vs Empty_Control
pudp_lrt=glmLRT(fit,  coef =  "GroupOE_PUDP")
pudp <- topTags(pudp_lrt, n=Inf)[[1]]

#Adding relevant information to the DEGs, like: direction, hgcn symbol
ensembl=useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
databaseGenesEnsembl <- getBM(
  attributes=c("ensembl_gene_id", "hgnc_symbol", "chromosome_name"),
  values=rownames(filtered),
  filters="ensembl_gene_id",
  mart=ensembl
)
#resutls
PUDP_results=merge(as.data.frame(pudp), databaseGenesEnsembl,by.x="row.names", by.y="ensembl_gene_id", all.x=TRUE)
PUDP_results$DE="N.S."
PUDP_results$FDR[which(PUDP_results$FDR==0)]=min(PUDP_results$FDR[which(PUDP_results$FDR!=0)])
PUDP_results$DE[which(PUDP_results$Row.names %in% rownames(deg_pudp)[which(deg_pudp$logFC>0)])]="UP"
PUDP_results$DE[which(PUDP_results$Row.names %in% rownames(deg_pudp)[which(deg_pudp$logFC<0)])]="DOWN"
PUDP_results$DE=as.factor(PUDP_results$DE)
PUDP_results$condlabel <- ifelse(PUDP_results$FDR  < 1e-40, as.character(PUDP_results$hgnc_symbol), NA)
#results for 7 samples
#PUDP_results7s=merge(as.data.frame(pudp), databaseGenesEnsembl,by.x="row.names", by.y="ensembl_gene_id", all.x=TRUE)

######_________________________________X link genes annotation_____________________________________ 
anno <- getBM(
  attributes=c("ensembl_gene_id", "start_position", "end_position", "chromosome_name"),
  values=rownames(filtered),
  filters="ensembl_gene_id",
  mart=ensembl
)
anno$isPAR=c("no")
#chromosome:GRCh38:Y:10001 - 2781479 is shared with X: 10001 - 2781479 (PAR1)
#chromosome:GRCh38:Y:56887903 - 57217415 is shared with X: 155701383 - 156030895 (PAR2)
# PARs are annotated as X chr
anno$isPAR[which(anno$chromosome_name=="X" & anno$start_position> 10000 & anno$end_position <= 2781479)]="PAR1"
anno$isPAR[which(anno$chromosome_name=="X" & anno$start_position>= 155701383 & anno$end_position <= 156030895)]="PAR2"
par1genes=anno$ensembl_gene_id[which(anno$isPAR=="PAR1")]
par2genes=anno$ensembl_gene_id[which(anno$isPAR=="PAR2")]
pargenes=unique(c(par1genes, par2genes))

escapes=read_excel("../../MiniTasks/Saudi_revision/ImportantGenes.xlsx", 
                   sheet = "escapes_Tukianien_filtered")$GeneID
inactive=read_excel("../../MiniTasks/Saudi_revision/ImportantGenes.xlsx", 
                    sheet = "inactive_Tukianien")$GeneID
variable=read_excel("../../MiniTasks/Saudi_revision/ImportantGenes.xlsx", 
                    sheet = "variable_Tukianien")$GeneID


# DEGs with FDR<0.05 for all samples
deg_pudp=as.data.frame(PUDP_results)[which(PUDP_results$FDR < 0.05) ,]
deg_pudp$DE='UP'
deg_pudp$DE[which(deg_pudp$logFC<0)]="DOWN"
rownames(deg_pudp) = deg_pudp$Row.names
deg_pudp=deg_pudp[,-1]
deg_pudp$Category = 'N.A'
deg_pudp$Category[which(deg_pudp$chromosome_name=='X')]="No category"
deg_pudp$Category[which(rownames(deg_pudp) %in% par1genes)]="PAR1"
deg_pudp$Category[which(rownames(deg_pudp) %in% par2genes)]="PAR2"
deg_pudp$Category[which(rownames(deg_pudp) %in% variable)]="Variable"
deg_pudp$Category[which(rownames(deg_pudp) %in% escapes)]="Non PAR escape"
deg_pudp$Category[which(rownames(deg_pudp) %in% inactive)]="Inactive"

# DEGs with FDR<0.05 for only 7 samples
#deg_pudp7=as.data.frame(PUDP_results7s)[which(PUDP_results7s$FDR < 0.05) ,]
#deg_pudp7$DE='UP'
#deg_pudp7$DE[which(deg_pudp7$logFC<0)]="DOWN"
#rownames(deg_pudp7) = deg_pudp7$Row.names
#deg_pudp7=deg_pudp7[,-1]
#deg_pudp7$Category = 'N.A'
#deg_pudp7$Category[which(deg_pudp7$chromosome_name=='X')]="No category"
#deg_pudp7$Category[which(rownames(deg_pudp7) %in% par1genes)]="PAR1"
#deg_pudp7$Category[which(rownames(deg_pudp7) %in% par2genes)]="PAR2"
#deg_pudp7$Category[which(rownames(deg_pudp7) %in% variable)]="Variable"
#deg_pudp7$Category[which(rownames(deg_pudp7) %in% escapes)]="Non PAR escape"
#deg_pudp7$Category[which(rownames(deg_pudp7) %in% inactive)]="Inactive"

#Save DEGs in excel file

write.xlsx(deg_pudp, file = "DEGs_PUDP.xlsx", rowNames=TRUE)

#__________________________PLOTS________________________________________________
threshold_pval=0.05

vp <- ggplot(data=PUDP_results, aes(x=logFC, y=-log10(FDR),col=DE,label=hgnc_symbol, ensembl=Row.names)) +
  geom_point(alpha = 0.5) + 
  theme_minimal() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-LFC, LFC), linetype="dashed") +
  geom_hline(yintercept=-log10(threshold_pval), linetype="dashed") +
  geom_text_repel(aes(label = condlabel), show.legend =FALSE, max.overlaps=100)

## ---------------------------------------------------------------------------------------------------------

