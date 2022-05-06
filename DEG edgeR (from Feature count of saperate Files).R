###################### Libraries ###########################
library(edgeR) # load in your R Session 
library(limma) # plotDensities

###################### Files Reading and Merging ###########################
cancer <- read.delim("/UCEC.txt",stringsAsFactors=F,row.names=1,header=T)      #Cancer Feature Counts File
normal <- read.delim("UCEC.txt",stringsAsFactors=F,row.names=1,header=T)      #Normal Feature Counts File
#### Merge both datasets into one
merge <- merge(normal, cancer, by ="row.names", all=TRUE)
merge[is.na(merge)] <- 0 
#Checking if columns are sum of both
#head(merge[,1:25])
print("Dimentions")
dim(cancer)
dim(normal)
dim(merge)
#Write Merge File
write.csv(merge,"UCEC_merged.csv",row.names=F)
#Read Merge File
raw.data<-read.csv("UCEC_merged.csv",header=T, row.names=1)

###################### Experiment MetaData Creation ###########################
# Counting total number of columns in Normal VS Cancer
normal_col <- ncol(normal)
cancer_col <- ncol(cancer)
merge_col <- ncol(raw.data)
# Creating Reps
WT <- rep("WT",normal_col)
MU <- rep("MU",cancer_col)
# Creating Experiment
expt<-c(WT,MU)
expt<-factor(expt, levels=c("WT", "MU"))
expt

###################### Library Size Calculation ###########################
library.sizes <- colSums(raw.data[,1:merge_col])
library.sizes

###################### RPM Normalization ###########################
#rpm.data <- matrix(data=NA, nrow=nrow(raw.data), ncol=merge_col, dimnames = list(row.names(raw.data), colnames(raw.data[,1:merge_col])))
# look at what you created:
#rpm.data[1:5,1:19]
#dim(rpm.data)

#for (i in 1:merge_col) {
#    rpm.data[,i] <- raw.data[,0+i] / (library.sizes[i]/1000000)
#    }
#rpm.data[1:5,1:19]
#raw.data[1:5,1:19]
#write.csv(rpm.data, "RPKMS.csv")

###################### Differential Gene Analyses using package edgeR ###########################
deg <- DGEList(counts=as.matrix(raw.data[,1:merge_col]), lib.size=library.sizes, group=expt)
deg$ samples
deg <- calcNormFactors(deg)
deg$ samples

###################### Now do DE testing using edgeR ###########################
deg <- estimateCommonDisp(deg)
# Expression differences can be tested using an exact test with a negative binomial distribution:
nbt <- exactTest(deg, pair=c("WT","MU"))
# See what is in the deg object:
names(nbt)
#the deg we want are in the $table
nbt$ table[1:6,]
# The p.values have not been adjusted for multiple hypothesis testing; This can be done with the topTags function. It will perform FDR (False Discovery Rate correction)
nbt.corrected <- topTags(nbt,n=Inf)
nbt.corrected[1:6,] #see deg for top 5 genes
#### This result can be output by
#write.csv(nbt.corrected$table,file="NBT.csv")

# How many genes are significant at FDR p < 0.05?
sum(nbt.corrected$table$FDR<0.05)
nbttable <- nbt.corrected$table
degs <- nbttable[nbttable$FDR<0.05,]
head(degs)
dim(degs)
degs$ID<-rownames(degs)
head(degs)

# separate up-regulated and down-regulated genes
deg_up<- degs[(degs$FDR<0.05 & degs$logFC>0),]
deg_up <- deg_up[order(deg_up$logFC, decreasing=T),]
dim(deg_up)
head(deg_up)
deg_down<- degs[(degs$FDR<0.05 & degs$logFC<0),]
deg_down<- deg_down[order(deg_down$logFC),]
dim(deg_down)
head(deg_down)

#Output your degs into Saperate files
write.csv(deg_up, file="UCEC_UP.csv", row.names=F)
write.csv(deg_down, file="UCEC_DOWN.csv", row.names=F)

# combine deg with raw counts
raw.data$ID <- rownames(raw.data)
head(raw.data)
combined <- merge(raw.data, degs, by='ID', all.ID=T)
head(combined)
combined<-combined[order(combined$logFC, decreasing=T),]
head(combined)
Change Name with Cancer Name
write.csv(combined, "UCEC_DEG.csv", , row.names=F)
rm(list=ls())
gc()