## run script ginkgo_aneupl_calls.Rmd first
library(here)
library(tidyverse)
load(here("raw_data", "gnk.RData")) # file is 2KB. I saved it after generating matrix using code in ginkgo_aneupl_calls.Rmd
load(here("raw_data","gnk.c.RData"))
gnk.u_DF <- as.data.frame(gnk.u)
gnk.c_DF <- as.data.frame(gnk.c)
# do test for senescent samples first
gnk.c_sen <- gnk.c_DF[grep("sen_0_0",gnk.c_DF[,1]),]
dim(gnk.c_sen) # 120  25
# turn DF into numeric matrix to make calculations easier
rownames(gnk.c_sen) <- gnk.c_sen[,2]
gnk.c_sen <- gnk.c_sen[,c(-1,-2)]
gnk.c_sen <- as.matrix(gnk.c_sen)
class(gnk.c_sen) <- "numeric"
# make any -1 into a 1. we only want to know changes, don't care about gain or loss right now
gnk.c_sen_abs <- t(apply(gnk.c_sen, 1,abs))

# average number of gains and losses per sample
rowSums(gnk.c_sen_abs)/ncol(gnk.c_sen_abs)
# average number of gains and losses per chromosome

nullDist_freq_each_chrom <- colSums(gnk.c_sen_abs)/nrow(gnk.c_sen_abs) # 23 numbers, one for each chrom




library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
library(tidyverse)
library(ggrepel)
# calculate the length and number of genes in each chromosome
chrs <- paste0("chr", c(1:22, "X"))
chr_lens <- seqlengths(BSgenome.Hsapiens.UCSC.hg19) %>% 
  .[chrs] %>% 
  set_names(substr(names(.), 4, nchar(.)))

genes_per_chr <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene) %>% 
  data.frame() %>% 
  as_tibble() %>%
  filter(seqnames %in% chrs) %>%
  group_by(seqnames) %>%
  count() %>%
  ungroup %>%
  mutate(seqnames = as.character(seqnames),
         chr= substr(seqnames, 4, nchar(seqnames))) %>%
  dplyr::select(-seqnames) #%>% 
  #rename(nGenes = n)

chr_tbl <- data.frame(len=chr_lens, 
  aneuFreqSen=nullDist_freq_each_chrom) %>%
  rownames_to_column("chr") %>%
  as_tibble() %>%
  left_join(genes_per_chr, by="chr") %>%
  mutate(genes_per_mb = n/len * 1e6,
         trisomy = chr %in% c("13", "18", "21"))
  

ggplot(chr_tbl, aes(x = aneuFreqSen, y = genes_per_mb, color=trisomy)) + geom_point()

ggplot(chr_tbl, aes(x = aneuFreqSen, y = len, color=trisomy, size=n)) + 
  geom_point(alpha=0.5) + theme_classic() + geom_text_repel(aes(label=chr))

ggplot(chr_tbl, aes(x = aneuFreqSen, y = n, color=trisomy)) + 
  geom_point(alpha=0.5, aes(size=len)) + 
  theme_classic() + 
  geom_text_repel(aes(label=chr)) + 
  xlab("Frequency of Aneuploidy") + 
  ylab("Number of Genes") + 
  ggtitle("Aneuploidy in Senescent Cells by Chromosome")


# for chromosome 18 do 100 bootsraps to select random samples, and then take average
set.seed(1)
chrom_18_aneuploidyFreq <- NULL
for(x in 1:1000) {
    chrom_18_aneuploidyFreq <- c(chrom_18_aneuploidyFreq,mean(gnk.c_sen_abs[sample(1:120, 120, replace=T),18]))
}
# plot distributions of aneu/monoploidy frequencies where you have the single nu
plot(density(nullDist_freq_each_chrom),xlim=c(0,0.14))
lines(density(chrom_18_aneuploidyFreq))
#add the frequency of chr 18 using all 120 files
abline(v=nullDist_freq_each_chrom[18],col="red")

# cannot compute ks test because we have ties, so values are not continuous
ks.test(nullDist_freq_each_chrom,chrom_18_aneuploidyFreq,alternative="less") # p-value = 0.9989

t.test(nullDist_freq_each_chrom,chrom_18_aneuploidyFreq,alternative="less") # p-value = 6.671e-13

wilcox.test(nullDist_freq_each_chrom,chrom_18_aneuploidyFreq,alternative="less") # p-value = 9.122e-12

### now do test where the null distribution is bootstrapped 100 times

#bootstrapping on null Dist values
set.seed(1)
nullDist_boots <- NULL
for(x in 1:1000) {
    randInts <- sample(1:120, 120, replace=T)
    nullDist_boots <- c(nullDist_boots,colSums(gnk.c_sen_abs[randInts,])/length(randInts))
}

# plot bootstrapped distributions of both null and chrom 18.
plot(density(nullDist_boots))
lines(density(chrom_18_aneuploidyFreq))
#add the frequency of chr 18 using all 120 files
abline(v=nullDist_freq_each_chrom[18],col="red")

## test all chromosomes
pVals_allChroms <- NULL
for(i in 1:23) {
    set.seed(1)
    chrom_aneuploidyFreq <- NULL
    for(x in 1:1000) {
        chrom_aneuploidyFreq <- c(chrom_aneuploidyFreq,mean(gnk.c_sen_abs[sample(1:120, 120, replace=T),i]))
    }
    pVals_allChroms <- c(pVals_allChroms,t.test(nullDist_freq_each_chrom,chrom_aneuploidyFreq,alternative="less")$p.value)
}
names(pVals_allChroms) <- c(paste("chr",1:22,sep=""),"chrX")
sig_0.05_chromNums <- which(pVals_allChroms < 0.05)
pVals_allChroms %>% p.adjust(method = "BH") %>% sort
print(sig_0.05_chromNums)
#chr7 chr15 chr17 chr18 chr19 chr20  chrX
#7    15    17    18    19    20    23
sig_0.01_chromNums <- which(pVals_allChroms < 0.01)
print(sig_0.01_chromNums)
#chr17 chr18 chr19 chr20  chrX
#17    18    19    20    23
sig_0.001_chromNums <- which(pVals_allChroms < 0.001)
print(sig_0.001_chromNums)
#chr17 chr18 chr19 chr20
#17    18    19    20

### compare distributions between 2 different chromosomes.

doBootstrapping <- function(chromNum) {
    set.seed(1)
    chrom_aneuploidyFreq <- NULL
    for(x in 1:1000) {
        chrom_aneuploidyFreq <- c(chrom_aneuploidyFreq,mean(gnk.c_sen_abs[sample(1:120, 120, replace=T),chromNum]))
    }
    return(chrom_aneuploidyFreq)
}
# matrix of aneu/mono ploidy freq. Each column is a distribution of aneu/monoploidy freq for a chromosome
allChromsBoots <- sapply(seq(1:23), function(x) doBootstrapping(x))
# get matrix of all combinations to perform t-test on
t_test_comparisons <- combn(seq(1:23),2)
# name the columns so I know which comparisons are significant
colnames(t_test_comparisons) <- paste(paste("chr",t_test_comparisons[1,],sep=""), paste("chr",t_test_comparisons[2,],sep=""),sep="_")
colnames(t_test_comparisons) <- gsub("chr23","chrX", colnames(t_test_comparisons))

# do 2 sided t test
chromosome_comparison_ttest <- apply(t_test_comparisons,2, function(x) t.test(allChromsBoots[,x[1]],allChromsBoots[,x[2]])$p.value)
# adjust for multiple testing
chromosome_comparison_ttest_fdr <- p.adjust(chromosome_comparison_ttest,method="BH")
chromosome_comparison_ttest_fdr_0.05sig <- sort(chromosome_comparison_ttest_fdr[chromosome_comparison_ttest_fdr<0.05]) # 204
chromosome_comparison_ttest_fdr_0.01sig <- sort(chromosome_comparison_ttest_fdr[chromosome_comparison_ttest_fdr<0.01]) # 204

# create matrix to do heatmap
pvalue_mat <- matrix(NA,ncol=23,nrow=23)
colnames(pvalue_mat) <- paste("chr",names(nullDist_freq_each_chrom),sep="")
rownames(pvalue_mat) <- paste("chr",names(nullDist_freq_each_chrom),sep="")
# find the row and column to insert each p-value into
pval_col_row <- lapply(strsplit(names(chromosome_comparison_ttest_fdr),split="_"), function(x) {
    matchCol <- match(x[1],colnames(pvalue_mat))
    matchRow <- match(x[2],rownames(pvalue_mat))
    return(c(matchCol,matchRow))
})
# put the correct p-value in the slot in the matrix
for(x in 1:length(pval_col_row)) {
    tempIndex <- pval_col_row[[x]]
    tempRow <- tempIndex[2]
    tempCol <- tempIndex[1]
    pvalue_mat[tempRow,tempCol] <- chromosome_comparison_ttest_fdr[x]
}
#pvalue_mat[pvalue_mat>0.05] <- NA

# create matrix of p-values comparing aneuploidy frequency in different chromosomes
library(gplots)
my_palette <- colorRampPalette(c("green", "yellow", "red"))(n = 105)
pdf("ploidy_variability.pdf")
heatmap.2(pvalue_mat,dendrogram ='none',Colv=FALSE,Rowv=FALSE, col=rev(my_palette),key=TRUE, keysize=1.3, density.info='none',trace='none')
dev.off()
# linear regression on chromosomes

# get encode data for IMR90 cell type
library(AnnotationHub)
ah <- AnnotationHub()

epiFiles_RoadMap <- query(ah, c("imr90","EpigenomeRoadMap"))
# custom pick the files I want. Gapped peak, broad peak, fold change files are available, but I took pval files
epiFiles_chromImpute <- epiFiles_RoadMap[523]
epiFiles_imputed_pval <- epiFiles_RoadMap[290:319]
epiFiles_methylation <- epiFiles_RoadMap[524] # bigwig
epiFiles_chromAccesability <- epiFiles_RoadMap[285:289] # bigwig
epiFiles_UCSD_pval <- epiFiles_RoadMap[219:284] # bigwig
epiFiles_E017_pval <- epiFiles_RoadMap[119:147]# bigwig
# DNase hotspot for broad peaks fdr 0.01.
E017_DNase.hotspot.fdr0.01 <- epiFiles_RoadMap[2] # granges

# combine UCSD and E017
epiFiles_UCSD_E017 <- c(epiFiles_UCSD_pval,epiFiles_E017_pval)

#rtracklayer uses the function import bw to convert bigwig files into granges
library(rtracklayer)
# make a list of granges objects, where each element is a granges for many different chromatin marks
epiOverlaps_chrom_list <- vector(mode = "list", length = length(epiFiles_UCSD_E017$ah_id))

names(epiOverlaps_chrom_list) <- epiFiles_UCSD_E017$title
# get significant peaks for all cHipSeq data
for(x in 1:length(epiFiles_UCSD_E017$title)) {
    print(x)
    myFile <- epiFiles_UCSD_E017[[epiFiles_UCSD_E017$ah_id[x]]]
    temp_bigWig <- import.bw(myFile)
    sig_bigWig <- temp_bigWig[mcols(temp_bigWig)[,1] > -log10(0.05)]
    epiOverlaps_chrom_list[[x]] <- sig_bigWig
}

# This is a matrix where rows are different chromosomes and columns are chromatin marks. Each cell in the matrix is the number of an epigenetic mark in a chromosome.
epiOverlaps_perChrom <- sapply(epiOverlaps_chrom_list, function(x) {
    #overlap_chroms_nums <- as.numeric(table(seqnames(x)))
    #names(overlap_chroms_nums) <- names(table(seqnames(x)))
    seqnames_temp <- factor(seqnames(x),levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10", "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM"))
    overlap_chroms_nums <- table(seqnames_temp)
})
rownames(epiOverlaps_perChrom) <- unique(as.character(seqnames(epiOverlaps_chrom_list[[1]])))

# now make a matrix where I get the number of each epigenetic mark in each bin.

# first import "gnk"
source(here("scripts", "gkStats.R"))
gnk <- read_tsv(here("raw_data", "2018_01_17_ginkgo_cnv_all.txt"))#, header=T)
library(plyr)

# create a matrix where the rows are different bins and the columns are chromatin marks. Each cell in the matrix is the number of an epigenetic mark in a bin.
gnk_bins_granges <- makeGRangesFromDataFrame(gnk[,1:3])
epiOverlaps_perBin <- sapply(epiOverlaps_chrom_list, function(x) {
    binOverlaps <- findOverlaps(gnk_bins_granges,x)
    binHits <- factor(queryHits(binOverlaps),levels=seq(1,length(gnk_bins_granges)))
    binHits_freqTable <- table(binHits)
})

epiOverlaps_perBin_coords <- cbind(gnk[,1:3],epiOverlaps_perBin)

## save the 2 matrixes so I don't have to do this all over again
saveRDS(epiOverlaps_perBin_coords,file="raw_data/epiOverlaps_perBin_coords.rds")
saveRDS(epiOverlaps_perChrom,file="raw_data/epiOverlaps_perChrom.rds")

# get 2 matrixes, one for bins and one for chromosomes. But this time, instead of a matrix where columns are epigenetic signatures, the columns of our matrix are number of methylated loci, and number of chromatin accessable regions in each chromosome/bin.

roadmapFiles_meth_access_hot <- c(epiFiles_methylation,epiFiles_chromAccesability,E017_DNase.hotspot.fdr0.01)

methyl_access_hot_list <- vector(mode = "list", length = length(roadmapFiles_meth_access_hot$ah_id))

names(methyl_access_hot_list) <- roadmapFiles_meth_access_hot$title
for(x in 1:length(roadmapFiles_meth_access_hot$title)) {
    print(x)
    myFile <- roadmapFiles_meth_access_hot[[roadmapFiles_meth_access_hot$ah_id[x]]] #
    temp_bigWig <- import.bw(myFile)
    if(names(methyl_access_hot_list)[x]!="E017_WGBS_FractionalMethylation.bigwig") {
        sig_bigWig <- temp_bigWig[mcols(temp_bigWig)[,1] > -log10(0.05)]
    } else {
        sig_bigWig <- temp_bigWig
    }
    methyl_access_hot_list[[x]] <- sig_bigWig
}

methyl_access_perChrom <- sapply(methyl_access_hot_list[1:6], function(x) {
    seqnames_temp <- factor(seqnames(x),levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10", "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM"))
    overlap_chroms_nums <- table(seqnames_temp)
})

methyl_access_perBin <- sapply(methyl_access_hot_list[1:6], function(x) {
    binOverlaps <- findOverlaps(gnk_bins_granges,x)
    binHits <- factor(queryHits(binOverlaps),levels=seq(1,length(gnk_bins_granges)))
    binHits_freqTable <- table(binHits)
})

# save these matrixes so I don't have to do all this over again. Its pretty time consuming.
saveRDS(methyl_access_perBin,file="raw_data/methyl_access_perBin.rds")
saveRDS(methyl_access_perChrom,file="raw_data/methyl_access_perChrom.rds")

# get the null distribution of aneu/monoploidy frequency in each bin. We already found the null distribution of the ploidy frequency in each chromosome
gnk_binary <- as.matrix(gnk[,-c(1,2,3)])
gnk_binary[gnk_binary==2] =0
gnk_binary[gnk_binary!=0] = 1
nullDist_freq_each_bin <- rowSums(gnk_binary)/ncol(gnk_binary) # 5363 numbers, one for each bin

## read in the 2 matrixes I saved earlier where each column is an epigenetic mark

epiOverlaps_perBin_coords <- readRDS("raw_data/epiOverlaps_perBin_coords.rds")
epiOverlaps_perChrom <- readRDS("raw_data/epiOverlaps_perChrom.rds")

# scale values so the lasso works better
epiOverlaps_perBin_coords_scaled <- apply(epiOverlaps_perBin_coords[,c(-1,-2,-3)],2, function(x) scale(x,center = TRUE, scale = TRUE))
range(epiOverlaps_perBin_coords_scaled) # -1.497328 42.539047
epiOverlaps_perChrom_coords_scaled <- apply(epiOverlaps_perChrom,2, function(x) scale(x,center = TRUE, scale = TRUE))
range(epiOverlaps_perChrom_coords_scaled) # -2.129125  2.641201
## only select a subset of predictor variables
# variable selection and coefficient shrinkage. Lasso does variable selection and coefficient shrinkage. Lasso is one end and redge regression on other. ridge regression does not perform variable selection. Ridge regression you impose a penalty, like lasso, but the type of penalty is so no coefficient will ever be set to 0. In between is the elastic net. Some variable selection (if lasso is to strict and impose to much of a penalty)

# output model or model output itself
library(glmnet)
logRegElNetLooc <- function(x, y, alpha, retGlmNt=FALSE,family,type.measure,nfolds){
    #creates a matrix of coefficients using elastic net regression. Y is variability statstic. X is matrix of predictors. Family is what you want to predict. linear reg is gausian. binomial is logistic regression (categorical), and others as well. Guasian for continuous, but try others as well. type.measure is classification. Alpha is between 0 and 1. 1 is lasso, 0 is ridge regression. nlambda is the number of folds for cross validations. this is leave one out. you want to pick model that gives highest accuracy or lowest error.
    CV <- cv.glmnet(x, y, family=family, type.measure = type.measure, alpha = alpha, nlambda = 100, nfolds = nfolds) #this is Leave one out Cross validation
    fit <- glmnet(x, y,family=family, alpha = alpha, lambda = CV$lambda.1se) # lambda.1se is the lambda that gave the smallest cross validation error and using it to compute final model. This computes final model with lambda set to ideal lambda based on cross validation.
    coeffic <- coef(fit) #get the coefficients from the fit model Extract coefficients and the higher the coefficient the more important. look at magnitude and sign (pos or negatively associated with outcome).
    if(retGlmNt){
        return(fit)
    }
    return(coeffic)
}

# we could also use biostrings for gc constant, As, Ts, ATs, GCperc. ChIPseq data, sequence composition data, and chromatin access

# do glmnet with function above for UCSD epigenetic marks

# separate chrom marks done by UCSD and other place. remove chromosomes Y and M and only include UCSD data
epiOverlaps_perChrom_UCSD_noYM <- epiOverlaps_perChrom_coords_scaled[c(1:23),30:95]

#UCSD data. For some reason chromatin marks are duplicated and I am not sure what this is, so i seperated the columns into two different matrixes
epiOverlaps_perChrom_UCSD_noYM_evens <- epiOverlaps_perChrom_UCSD_noYM[,c(TRUE,FALSE)]
epiOverlaps_perChrom_UCSD_noYM_odds <- epiOverlaps_perChrom_UCSD_noYM[,c(FALSE,TRUE)]

# get chrom marks done with E017
epiOverlaps_perChrom_E017_noYM <- epiOverlaps_perChrom_coords_scaled[c(1:23),1:29]

### examine if there are any highly correlated features.
epiOverlaps_perChrom_UCSD_noYM_cor <- cor(epiOverlaps_perChrom_UCSD_noYM) # H3K9me3 is only one thats not super highly correlated with everything else
epiOverlaps_perChrom_E017_noYM_cor <- cor(epiOverlaps_perChrom_E017_noYM)

epiOverlaps_perBin_coords_E017 <- epiOverlaps_perBin_coords_scaled
epiOverlaps_perBin_coords_E017_cor <- cor(epiOverlaps_perBin_coords_E017)

epiOverlaps_perBin_coords_IMR90 <- epiOverlaps_perBin_coords_scaled[,30:95]
epiOverlaps_perBin_coords_IR90_cor <- cor(epiOverlaps_perBin_coords_IMR90)

which(apply(epiOverlaps_perBin_coords_IR90_cor,1,mean) < abs(0.5)) # H3K27me3 and H3K9me3 are only chromatin marks not highly correlated


which(apply(epiOverlaps_perBin_coords_E017_cor,1,mean) < abs(0.5)) # H3K27me3 and H3K9me3 are only chromatin marks not highly correlated

# glm net per chromosome for alll UCSD data
epiOverlaps_perChrom_UCSD_coefs <- logRegElNetLooc(epiOverlaps_perChrom_UCSD_noYM,nullDist_freq_each_chrom,alpha=0.5,family="gaussian",type.measure="deviance",nfolds=7) # grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold

# glmnet per chromosome for evens UCSD data
epiOverlaps_perChrom_UCSD_evens_coefs <- logRegElNetLooc(epiOverlaps_perChrom_UCSD_noYM_evens,nullDist_freq_each_chrom,alpha=0.5,family="gaussian",type.measure="deviance",nfolds=7) # grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold

# glmnet per chromosome for odd UCSD data
epiOverlaps_perChrom_UCSD_odds_coefs <- logRegElNetLooc(epiOverlaps_perChrom_UCSD_noYM_odds,nullDist_freq_each_chrom,alpha=0.5,family="gaussian",type.measure="deviance",nfolds=7) # grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold

# glmnet per chromosome for E017 data
epiOverlaps_perChrom_E017_coefs <- logRegElNetLooc(epiOverlaps_perChrom_E017_noYM,nullDist_freq_each_chrom,alpha=0.5,family="gaussian",type.measure="deviance",nfolds=7) # grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold

#### try binomial where you have chromosome gain or loss. 1 means that the chromosome has aneuploidy. 0 means its diploid usually
chrom_gain_loss <- rep(0,23)
chrom_gain_loss[which(pVals_allChroms < 0.01)] <- 1

## do UCSD (all, odds, and evens) and E017 as above but this time use logistic regression. I am getting errors due to cross validation
epiOverlaps_perChrom_UCSD_coefs_binomial <- logRegElNetLooc(epiOverlaps_perChrom_UCSD_noYM,chrom_gain_loss,alpha=0.5,family="binomial",type.measure="class",nfolds=3) #mse used by default. grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold

epiOverlaps_perChrom_UCSD_odds_coefs_binomial <- logRegElNetLooc(epiOverlaps_perChrom_UCSD_noYM_odds,chrom_gain_loss,alpha=0.5,family="binomial",type.measure="class",nfolds=3) # grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold

epiOverlaps_perChrom_E017_coefs_binomial <- logRegElNetLooc(epiOverlaps_perChrom_E017_noYM,chrom_gain_loss,alpha=0.5,family="binomial",type.measure="class",nfolds=3) # grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold

epiOverlaps_perChrom_E017_coefs_binomial <- logRegElNetLooc(epiOverlaps_perChrom_E017_noYM,chrom_gain_loss,alpha=0.5,family="binomial",type.measure="class",nfolds=3) #mse used by default. grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold

# instead of using lasso or rig regression just put all data points in model for multiple regression in R

library(MASS)
# first change colum names because R does not like "-". Instead use "_"
epiOverlaps_perChrom_E017_noYM_noHyp <- epiOverlaps_perChrom_E017_noYM
colnames(epiOverlaps_perChrom_E017_noYM_noHyp) <- gsub("-","_",colnames(epiOverlaps_perChrom_E017_noYM))
# there are NAs in this data because some of data is very correlated with other features.
fit <- lm(nullDist_freq_each_chrom~E017_DNase.pval.signal.txt+E017_H2A.Z.pval.signal.txt+E017_H2AK5ac.pval.signal.txt+E017_H2AK9ac.pval.signal.txt+E017_H2BK120ac.pval.signal.txt+E017_H2BK12ac.pval.signal.txt+E017_H2BK15ac.pval.signal.txt+E017_H2BK20ac.pval.signal.txt+E017_H2BK5ac.pval.signal.txt+E017_H3K14ac.pval.signal.txt+E017_H3K18ac.pval.signal.txt+E017_H3K23ac.pval.signal.txt+E017_H3K27ac.pval.signal.txt+E017_H3K27me3.pval.signal.txt+E017_H3K36me3.pval.signal.txt+E017_H3K4ac.pval.signal.txt+E017_H3K4me1.pval.signal.txt+E017_H3K4me2.pval.signal.txt+E017_H3K4me3.pval.signal.txt+E017_H3K56ac.pval.signal.txt+E017_H3K79me1.pval.signal.txt+E017_H3K79me2.pval.signal.txt+E017_H3K9ac.pval.signal.txt+E017_H3K9me1.pval.signal.txt+E017_H3K9me3.pval.signal.txt+E017_H4K20me1.pval.signal.txt+E017_H4K5ac.pval.signal.txt+E017_H4K8ac.pval.signal.txt+E017_H4K91ac.pval.signal.txt,data=as.data.frame(epiOverlaps_perChrom_E017_noYM_noHyp))
# remove NA features
fit_E017_noNA <- lm(nullDist_freq_each_chrom~E017_DNase.pval.signal.txt+E017_H2A.Z.pval.signal.txt+E017_H2AK5ac.pval.signal.txt+E017_H2AK9ac.pval.signal.txt+E017_H2BK120ac.pval.signal.txt+E017_H2BK12ac.pval.signal.txt+E017_H2BK15ac.pval.signal.txt+E017_H2BK20ac.pval.signal.txt+E017_H2BK5ac.pval.signal.txt+E017_H3K14ac.pval.signal.txt+E017_H3K18ac.pval.signal.txt+E017_H3K23ac.pval.signal.txt+E017_H3K27ac.pval.signal.txt+E017_H3K27me3.pval.signal.txt+E017_H3K36me3.pval.signal.txt+E017_H3K4ac.pval.signal.txt+E017_H3K4me1.pval.signal.txt+E017_H3K4me2.pval.signal.txt+E017_H3K4me3.pval.signal.txt+E017_H3K56ac.pval.signal.txt+E017_H3K79me1.pval.signal.txt+E017_H3K79me2.pval.signal.txt,data=as.data.frame(epiOverlaps_perChrom_E017_noYM_noHyp))

epiOverlaps_perChrom_UCSD_noYM_noHyp <- epiOverlaps_perChrom_UCSD_noYM
colnames(epiOverlaps_perChrom_UCSD_noYM_noHyp) <- gsub("-","_",colnames(epiOverlaps_perChrom_UCSD_noYM_noHyp))

fit_IMR90_chroms <- lm(nullDist_freq_each_chrom~UCSD.IMR90.H2A.Z.AY4.pval.signal.txt+UCSD.IMR90.H2A.Z.AY5.pval.signal.txt+UCSD.IMR90.H2AK5ac.YL20.pval.signal.txt+UCSD.IMR90.H2AK5ac.YL60.pval.signal.txt+UCSD.IMR90.H2AK9ac.AY6.pval.signal.txt+UCSD.IMR90.H2AK9ac.AY7.pval.signal.txt+UCSD.IMR90.H2BK120ac.SK15.pval.signal.txt+UCSD.IMR90.H2BK120ac.YL06.pval.signal.txt+UCSD.IMR90.H2BK12ac.YL138.pval.signal.txt+UCSD.IMR90.H2BK12ac.YL18.pval.signal.txt+UCSD.IMR90.H2BK12ac.YL68.pval.signal.txt+UCSD.IMR90.H2BK15ac.SK42.pval.signal.txt+UCSD.IMR90.H2BK15ac.YL19.pval.signal.txt+UCSD.IMR90.H2BK15ac.YL69.pval.signal.txt+UCSD.IMR90.H2BK20ac.SK14.pval.signal.txt+UCSD.IMR90.H2BK20ac.YL07.pval.signal.txt+UCSD.IMR90.H2BK5ac.LL229.pval.signal.txt+UCSD.IMR90.H2BK5ac.SK10.pval.signal.txt+UCSD.IMR90.H2BK5ac.YL281.pval.signal.txt+UCSD.IMR90.H3K14ac.SK17.pval.signal.txt+UCSD.IMR90.H3K14ac.YL15.pval.signal.txt+UCSD.IMR90.H3K18ac.LL230.pval.signal.txt+UCSD.IMR90.H3K18ac.SK09.pval.signal.txt+UCSD.IMR90.H3K23ac.SK18.pval.signal.txt+UCSD.IMR90.H3K23ac.YL09.pval.signal.txt+UCSD.IMR90.H3K27ac.LL231.pval.signal.txt+UCSD.IMR90.H3K27ac.LL235.pval.signal.txt+UCSD.IMR90.H3K27ac.YL58.pval.signal.txt+UCSD.IMR90.H3K27me3.LL223.pval.signal.txt+UCSD.IMR90.H3K27me3.SK05.pval.signal.txt+UCSD.IMR90.H3K36me3.LL225.pval.signal.txt+UCSD.IMR90.H3K36me3.YL67.pval.signal.txt+UCSD.IMR90.H3K4ac.LL239.pval.signal.txt+UCSD.IMR90.H3K4ac.SK16.pval.signal.txt+UCSD.IMR90.H3K4me1.LL222.pval.signal.txt+UCSD.IMR90.H3K4me1.SK04.pval.signal.txt+UCSD.IMR90.H3K4me1.SYL133.pval.signal.txt+UCSD.IMR90.H3K4me2.LL232.pval.signal.txt+UCSD.IMR90.H3K4me2.SK07.pval.signal.txt+UCSD.IMR90.H3K4me3.LL221.pval.signal.txt+UCSD.IMR90.H3K4me3.SK03.pval.signal.txt+UCSD.IMR90.H3K56ac.YL05.pval.signal.txt+UCSD.IMR90.H3K56ac.YL61.pval.signal.txt+UCSD.IMR90.H3K79me1.SK11.pval.signal.txt+UCSD.IMR90.H3K79me1.YL103.pval.signal.txt+UCSD.IMR90.H3K79me1.YL145.pval.signal.txt+UCSD.IMR90.H3K79me1.YL16.pval.signal.txt+UCSD.IMR90.H3K79me2.LL238.pval.signal.txt+UCSD.IMR90.H3K79me2.SK43.pval.signal.txt+UCSD.IMR90.H3K9ac.LL233.pval.signal.txt+UCSD.IMR90.H3K9ac.SK46.pval.signal.txt+UCSD.IMR90.H3K9me1.LL237.pval.signal.txt+UCSD.IMR90.H3K9me1.SK12.pval.signal.txt+UCSD.IMR90.H3K9me3.LL224.pval.signal.txt+UCSD.IMR90.H3K9me3.SK06.pval.signal.txt+UCSD.IMR90.H3K9me3.YL104.pval.signal.txt+UCSD.IMR90.H4K20me1.YL21.pval.signal.txt+UCSD.IMR90.H4K20me1.YL59.pval.signal.txt+UCSD.IMR90.H4K5ac.LL234.pval.signal.txt+UCSD.IMR90.H4K5ac.SK08.pval.signal.txt+UCSD.IMR90.H4K8ac.SK45.pval.signal.txt+UCSD.IMR90.H4K8ac.YL142.pval.signal.txt+UCSD.IMR90.H4K8ac.YL17.pval.signal.txt+UCSD.IMR90.H4K8ac.YL57.pval.signal.txt+UCSD.IMR90.H4K91ac.SK49.pval.signal.txt+UCSD.IMR90.H4K91ac.YL08.pval.signal.txt,data=as.data.frame(epiOverlaps_perChrom_UCSD_noYM_noHyp))
# break into evens and odds so less features are correlated
fit_IMR90_chroms_evens <- lm(nullDist_freq_each_chrom~UCSD.IMR90.H2A.Z.AY4.pval.signal.txt+UCSD.IMR90.H2AK5ac.YL20.pval.signal.txt+UCSD.IMR90.H2AK9ac.AY6.pval.signal.txt+UCSD.IMR90.H2BK120ac.SK15.pval.signal.txt+UCSD.IMR90.H2BK12ac.YL138.pval.signal.txt+UCSD.IMR90.H2BK12ac.YL68.pval.signal.txt+UCSD.IMR90.H2BK15ac.YL19.pval.signal.txt+UCSD.IMR90.H2BK20ac.SK14.pval.signal.txt+UCSD.IMR90.H2BK5ac.LL229.pval.signal.txt+UCSD.IMR90.H2BK5ac.YL281.pval.signal.txt+UCSD.IMR90.H3K14ac.YL15.pval.signal.txt+UCSD.IMR90.H3K18ac.SK09.pval.signal.txt+UCSD.IMR90.H3K23ac.YL09.pval.signal.txt+UCSD.IMR90.H3K27ac.LL235.pval.signal.txt+UCSD.IMR90.H3K27me3.LL223.pval.signal.txt+UCSD.IMR90.H3K36me3.LL225.pval.signal.txt+UCSD.IMR90.H3K4ac.LL239.pval.signal.txt+UCSD.IMR90.H3K4me1.LL222.pval.signal.txt+UCSD.IMR90.H3K4me1.SYL133.pval.signal.txt+UCSD.IMR90.H3K4me2.SK07.pval.signal.txt+UCSD.IMR90.H3K4me3.SK03.pval.signal.txt+UCSD.IMR90.H3K56ac.YL61.pval.signal.txt+UCSD.IMR90.H3K79me1.YL103.pval.signal.txt+UCSD.IMR90.H3K79me1.YL16.pval.signal.txt+UCSD.IMR90.H3K79me2.SK43.pval.signal.txt+UCSD.IMR90.H3K9ac.SK46.pval.signal.txt+UCSD.IMR90.H3K9me1.SK12.pval.signal.txt+UCSD.IMR90.H3K9me3.SK06.pval.signal.txt+UCSD.IMR90.H4K20me1.YL21.pval.signal.txt+UCSD.IMR90.H4K5ac.LL234.pval.signal.txt+UCSD.IMR90.H4K8ac.SK45.pval.signal.txt+UCSD.IMR90.H4K8ac.YL17.pval.signal.txt+UCSD.IMR90.H4K91ac.SK49.pval.signal.txt,data=as.data.frame(epiOverlaps_perChrom_UCSD_noYM_noHyp))

fit_IMR90_chroms_odds <- lm(nullDist_freq_each_chrom~UCSD.IMR90.H2A.Z.AY4.pval.signal.txt+UCSD.IMR90.H2AK5ac.YL20.pval.signal.txt+UCSD.IMR90.H2AK9ac.AY6.pval.signal.txt+UCSD.IMR90.H2BK120ac.SK15.pval.signal.txt+UCSD.IMR90.H2BK12ac.YL138.pval.signal.txt+UCSD.IMR90.H2BK12ac.YL68.pval.signal.txt+UCSD.IMR90.H2BK15ac.YL19.pval.signal.txt+UCSD.IMR90.H2BK20ac.SK14.pval.signal.txt+UCSD.IMR90.H2BK5ac.LL229.pval.signal.txt+UCSD.IMR90.H2BK5ac.YL281.pval.signal.txt+UCSD.IMR90.H3K14ac.YL15.pval.signal.txt+UCSD.IMR90.H3K18ac.SK09.pval.signal.txt+UCSD.IMR90.H3K23ac.YL09.pval.signal.txt+UCSD.IMR90.H3K27ac.LL235.pval.signal.txt+UCSD.IMR90.H3K27me3.LL223.pval.signal.txt+UCSD.IMR90.H3K36me3.LL225.pval.signal.txt+UCSD.IMR90.H3K4ac.LL239.pval.signal.txt+UCSD.IMR90.H3K4me1.LL222.pval.signal.txt+UCSD.IMR90.H3K4me1.SYL133.pval.signal.txt+UCSD.IMR90.H3K4me2.SK07.pval.signal.txt+UCSD.IMR90.H3K4me3.SK03.pval.signal.txt+UCSD.IMR90.H3K56ac.YL61.pval.signal.txt+UCSD.IMR90.H3K79me1.YL103.pval.signal.txt+UCSD.IMR90.H3K79me1.YL16.pval.signal.txt+UCSD.IMR90.H3K79me2.SK43.pval.signal.txt+UCSD.IMR90.H3K9ac.SK46.pval.signal.txt+UCSD.IMR90.H3K9me1.SK12.pval.signal.txt+UCSD.IMR90.H3K9me3.SK06.pval.signal.txt+UCSD.IMR90.H4K20me1.YL21.pval.signal.txt+UCSD.IMR90.H4K5ac.LL234.pval.signal.txt+UCSD.IMR90.H4K8ac.SK45.pval.signal.txt+UCSD.IMR90.H4K8ac.YL17.pval.signal.txt+UCSD.IMR90.H4K91ac.SK49.pval.signal.txt,data=as.data.frame(epiOverlaps_perChrom_UCSD_noYM_noHyp))

# make linear model just with E017_H3K27me3 and E017_H3K9me3 since these are not as highly correlated with any other marks
fi_E017_noNA_smallCorr <- lm(nullDist_freq_each_chrom~ E017_H3K27me3.pval.signal.txt+ E017_H3K9me3.pval.signal.txt,data=as.data.frame(epiOverlaps_perChrom_E017_noYM_noHyp))

# make linear models but using bins instead of chromosomes

epiOverlaps_perBin_coords_E017_nohyph <- epiOverlaps_perBin_coords_E017
colnames(epiOverlaps_perBin_coords_E017_nohyph) <- gsub("-","_",colnames(epiOverlaps_perBin_coords_E017_nohyph))
# linear regression model with E017 data and bins
fit_bins_E017 <- lm(nullDist_freq_each_bin~E017_DNase.pval.signal.txt+E017_H2A.Z.pval.signal.txt+E017_H2AK5ac.pval.signal.txt+E017_H2AK9ac.pval.signal.txt+E017_H2BK120ac.pval.signal.txt+E017_H2BK12ac.pval.signal.txt+E017_H2BK15ac.pval.signal.txt+E017_H2BK20ac.pval.signal.txt+E017_H2BK5ac.pval.signal.txt+E017_H3K14ac.pval.signal.txt+E017_H3K18ac.pval.signal.txt+E017_H3K23ac.pval.signal.txt+E017_H3K27ac.pval.signal.txt+E017_H3K27me3.pval.signal.txt+E017_H3K36me3.pval.signal.txt+E017_H3K4ac.pval.signal.txt+E017_H3K4me1.pval.signal.txt+E017_H3K4me2.pval.signal.txt+E017_H3K4me3.pval.signal.txt+E017_H3K56ac.pval.signal.txt+E017_H3K79me1.pval.signal.txt+E017_H3K79me2.pval.signal.txt+E017_H3K9ac.pval.signal.txt+E017_H3K9me1.pval.signal.txt+E017_H3K9me3.pval.signal.txt+E017_H4K20me1.pval.signal.txt+E017_H4K5ac.pval.signal.txt+E017_H4K8ac.pval.signal.txt+E017_H4K91ac.pval.signal.txt,data=as.data.frame(epiOverlaps_perBin_coords_E017_nohyph))

fit_IMR90_bins_evens <- lm(nullDist_freq_each_chrom~UCSD.IMR90.H2A.Z.AY4.pval.signal.txt+UCSD.IMR90.H2AK5ac.YL20.pval.signal.txt+UCSD.IMR90.H2AK9ac.AY6.pval.signal.txt+UCSD.IMR90.H2BK120ac.SK15.pval.signal.txt+UCSD.IMR90.H2BK12ac.YL138.pval.signal.txt+UCSD.IMR90.H2BK12ac.YL68.pval.signal.txt+UCSD.IMR90.H2BK15ac.YL19.pval.signal.txt+UCSD.IMR90.H2BK20ac.SK14.pval.signal.txt+UCSD.IMR90.H2BK5ac.LL229.pval.signal.txt+UCSD.IMR90.H2BK5ac.YL281.pval.signal.txt+UCSD.IMR90.H3K14ac.YL15.pval.signal.txt+UCSD.IMR90.H3K18ac.SK09.pval.signal.txt+UCSD.IMR90.H3K23ac.YL09.pval.signal.txt+UCSD.IMR90.H3K27ac.LL235.pval.signal.txt+UCSD.IMR90.H3K27me3.LL223.pval.signal.txt+UCSD.IMR90.H3K36me3.LL225.pval.signal.txt+UCSD.IMR90.H3K4ac.LL239.pval.signal.txt+UCSD.IMR90.H3K4me1.LL222.pval.signal.txt+UCSD.IMR90.H3K4me1.SYL133.pval.signal.txt+UCSD.IMR90.H3K4me2.SK07.pval.signal.txt+UCSD.IMR90.H3K4me3.SK03.pval.signal.txt+UCSD.IMR90.H3K56ac.YL61.pval.signal.txt+UCSD.IMR90.H3K79me1.YL103.pval.signal.txt+UCSD.IMR90.H3K79me1.YL16.pval.signal.txt+UCSD.IMR90.H3K79me2.SK43.pval.signal.txt+UCSD.IMR90.H3K9ac.SK46.pval.signal.txt+UCSD.IMR90.H3K9me1.SK12.pval.signal.txt+UCSD.IMR90.H3K9me3.SK06.pval.signal.txt+UCSD.IMR90.H4K20me1.YL21.pval.signal.txt+UCSD.IMR90.H4K5ac.LL234.pval.signal.txt+UCSD.IMR90.H4K8ac.SK45.pval.signal.txt+UCSD.IMR90.H4K8ac.YL17.pval.signal.txt+UCSD.IMR90.H4K91ac.SK49.pval.signal.txt,data=as.data.frame(epiOverlaps_perChrom_UCSD_noYM_noHyp))

fit_IMR90_bins_odds <- lm(nullDist_freq_each_chrom~UCSD.IMR90.H2A.Z.AY4.pval.signal.txt+UCSD.IMR90.H2AK5ac.YL20.pval.signal.txt+UCSD.IMR90.H2AK9ac.AY6.pval.signal.txt+UCSD.IMR90.H2BK120ac.SK15.pval.signal.txt+UCSD.IMR90.H2BK12ac.YL138.pval.signal.txt+UCSD.IMR90.H2BK12ac.YL68.pval.signal.txt+UCSD.IMR90.H2BK15ac.YL19.pval.signal.txt+UCSD.IMR90.H2BK20ac.SK14.pval.signal.txt+UCSD.IMR90.H2BK5ac.LL229.pval.signal.txt+UCSD.IMR90.H2BK5ac.YL281.pval.signal.txt+UCSD.IMR90.H3K14ac.YL15.pval.signal.txt+UCSD.IMR90.H3K18ac.SK09.pval.signal.txt+UCSD.IMR90.H3K23ac.YL09.pval.signal.txt+UCSD.IMR90.H3K27ac.LL235.pval.signal.txt+UCSD.IMR90.H3K27me3.LL223.pval.signal.txt+UCSD.IMR90.H3K36me3.LL225.pval.signal.txt+UCSD.IMR90.H3K4ac.LL239.pval.signal.txt+UCSD.IMR90.H3K4me1.LL222.pval.signal.txt+UCSD.IMR90.H3K4me1.SYL133.pval.signal.txt+UCSD.IMR90.H3K4me2.SK07.pval.signal.txt+UCSD.IMR90.H3K4me3.SK03.pval.signal.txt+UCSD.IMR90.H3K56ac.YL61.pval.signal.txt+UCSD.IMR90.H3K79me1.YL103.pval.signal.txt+UCSD.IMR90.H3K79me1.YL16.pval.signal.txt+UCSD.IMR90.H3K79me2.SK43.pval.signal.txt+UCSD.IMR90.H3K9ac.SK46.pval.signal.txt+UCSD.IMR90.H3K9me1.SK12.pval.signal.txt+UCSD.IMR90.H3K9me3.SK06.pval.signal.txt+UCSD.IMR90.H4K20me1.YL21.pval.signal.txt+UCSD.IMR90.H4K5ac.LL234.pval.signal.txt+UCSD.IMR90.H4K8ac.SK45.pval.signal.txt+UCSD.IMR90.H4K8ac.YL17.pval.signal.txt+UCSD.IMR90.H4K91ac.SK49.pval.signal.txt,data=as.data.frame(epiOverlaps_perChrom_UCSD_noYM_noHyp))

#make linear model just with E017_H3K27me3 and E017_H3K9me3 since these are not as highly correlated with any other marks
fi_E017_noNA_smallCorr_bin <- lm(nullDist_freq_each_bin~ E017_H3K27me3.pval.signal.txt+ E017_H3K9me3.pval.signal.txt,data=as.data.frame(epiOverlaps_perBin_coords_nohyph))


epiOverlaps_perBin_coords_IMR90_nohyph <- epiOverlaps_perBin_coords_IMR90
colnames(epiOverlaps_perBin_coords_IMR90_nohyph) <- gsub("-","_",colnames(epiOverlaps_perBin_coords_IMR90_nohyph))

fit_bins_IMR90 <- lm(nullDist_freq_each_bin~UCSD.IMR90.H2A.Z.AY4.pval.signal.txt+UCSD.IMR90.H2A.Z.AY5.pval.signal.txt+UCSD.IMR90.H2AK5ac.YL20.pval.signal.txt+UCSD.IMR90.H2AK5ac.YL60.pval.signal.txt+UCSD.IMR90.H2AK9ac.AY6.pval.signal.txt+UCSD.IMR90.H2AK9ac.AY7.pval.signal.txt+UCSD.IMR90.H2BK120ac.SK15.pval.signal.txt+UCSD.IMR90.H2BK120ac.YL06.pval.signal.txt+UCSD.IMR90.H2BK12ac.YL138.pval.signal.txt+UCSD.IMR90.H2BK12ac.YL18.pval.signal.txt+UCSD.IMR90.H2BK12ac.YL68.pval.signal.txt+UCSD.IMR90.H2BK15ac.SK42.pval.signal.txt+UCSD.IMR90.H2BK15ac.YL19.pval.signal.txt+UCSD.IMR90.H2BK15ac.YL69.pval.signal.txt+UCSD.IMR90.H2BK20ac.SK14.pval.signal.txt+UCSD.IMR90.H2BK20ac.YL07.pval.signal.txt+UCSD.IMR90.H2BK5ac.LL229.pval.signal.txt+UCSD.IMR90.H2BK5ac.SK10.pval.signal.txt+UCSD.IMR90.H2BK5ac.YL281.pval.signal.txt+UCSD.IMR90.H3K14ac.SK17.pval.signal.txt+UCSD.IMR90.H3K14ac.YL15.pval.signal.txt+UCSD.IMR90.H3K18ac.LL230.pval.signal.txt+UCSD.IMR90.H3K18ac.SK09.pval.signal.txt+UCSD.IMR90.H3K23ac.SK18.pval.signal.txt+UCSD.IMR90.H3K23ac.YL09.pval.signal.txt+UCSD.IMR90.H3K27ac.LL231.pval.signal.txt+UCSD.IMR90.H3K27ac.LL235.pval.signal.txt+UCSD.IMR90.H3K27ac.YL58.pval.signal.txt+UCSD.IMR90.H3K27me3.LL223.pval.signal.txt+UCSD.IMR90.H3K27me3.SK05.pval.signal.txt+UCSD.IMR90.H3K36me3.LL225.pval.signal.txt+UCSD.IMR90.H3K36me3.YL67.pval.signal.txt+UCSD.IMR90.H3K4ac.LL239.pval.signal.txt+UCSD.IMR90.H3K4ac.SK16.pval.signal.txt+UCSD.IMR90.H3K4me1.LL222.pval.signal.txt+UCSD.IMR90.H3K4me1.SK04.pval.signal.txt+UCSD.IMR90.H3K4me1.SYL133.pval.signal.txt+UCSD.IMR90.H3K4me2.LL232.pval.signal.txt+UCSD.IMR90.H3K4me2.SK07.pval.signal.txt+UCSD.IMR90.H3K4me3.LL221.pval.signal.txt+UCSD.IMR90.H3K4me3.SK03.pval.signal.txt+UCSD.IMR90.H3K56ac.YL05.pval.signal.txt+UCSD.IMR90.H3K56ac.YL61.pval.signal.txt+UCSD.IMR90.H3K79me1.SK11.pval.signal.txt+UCSD.IMR90.H3K79me1.YL103.pval.signal.txt+UCSD.IMR90.H3K79me1.YL145.pval.signal.txt+UCSD.IMR90.H3K79me1.YL16.pval.signal.txt+UCSD.IMR90.H3K79me2.LL238.pval.signal.txt+UCSD.IMR90.H3K79me2.SK43.pval.signal.txt+UCSD.IMR90.H3K9ac.LL233.pval.signal.txt+UCSD.IMR90.H3K9ac.SK46.pval.signal.txt+UCSD.IMR90.H3K9me1.LL237.pval.signal.txt+UCSD.IMR90.H3K9me1.SK12.pval.signal.txt+UCSD.IMR90.H3K9me3.LL224.pval.signal.txt+UCSD.IMR90.H3K9me3.SK06.pval.signal.txt+UCSD.IMR90.H3K9me3.YL104.pval.signal.txt+UCSD.IMR90.H4K20me1.YL21.pval.signal.txt+UCSD.IMR90.H4K20me1.YL59.pval.signal.txt+UCSD.IMR90.H4K5ac.LL234.pval.signal.txt+UCSD.IMR90.H4K5ac.SK08.pval.signal.txt+UCSD.IMR90.H4K8ac.SK45.pval.signal.txt+UCSD.IMR90.H4K8ac.YL142.pval.signal.txt+UCSD.IMR90.H4K8ac.YL17.pval.signal.txt+UCSD.IMR90.H4K8ac.YL57.pval.signal.txt+UCSD.IMR90.H4K91ac.SK49.pval.signal.txt+UCSD.IMR90.H4K91ac.YL08.pval.signal.txt,data=as.data.frame(epiOverlaps_perBin_coords_IMR90_nohyph))


## glmnet per bin with E017 data
epiOverlaps_perChrom_E017_coefs_linReg <- logRegElNetLooc(as.matrix(epiOverlaps_perBin_coords_E017),nullDist_freq_each_bin,alpha=0.5,family="gaussian",type.measure="deviance",nfolds=7)
epiOverlaps_perChrom_UCSD_coefs_linReg <- logRegElNetLooc(as.matrix(epiOverlaps_perBin_coords_IMR90),nullDist_freq_each_bin,alpha=0.5,family="gaussian",type.measure="deviance",nfolds=7)
epiOverlaps_perChrom_UCSD_coefs_linReg_odds <- logRegElNetLooc(as.matrix(epiOverlaps_perBin_coords_IMR90),nullDist_freq_each_bin,alpha=0.5,family="gaussian",type.measure="deviance",nfolds=7)

## now do glmnet and create linear models using methylation, chromatin accessability data
methyl_access_perBin <- readRDS(file="raw_data/methyl_access_perBin.rds")
methyl_access_perChrom <- readRDS(file="raw_data/methyl_access_perChrom.rds")

methyl_access_perBin_scaled <- scale(methyl_access_perBin,center=T,scale=T)
methyl_access_perChrom_scaled <- scale(methyl_access_perChrom,center=T,scale=T)


# look for correlations
methyl_access_perChrom_cor_perchrom <- cor(methyl_access_perChrom_scaled)
methyl_access_perBin_cor_perbin <- cor(methyl_access_perBin_scaled)

# make linear model wih null dist. For the model I used the bisulfite sequencing data. 1 of the chromatin accessibility regions, genes per chromosome scaled, and genes per MB scaled

fit_methyl_access_perchrom <- lm(nullDist_freq_each_chrom~E017_WGBS_FractionalMethylation.bigwig+UW.IMR90.ChromatinAccessibility.DS13229.pval.signal.bigwig+data.frame(genes_per_chr)[,1]+data.frame(chr_tbl)[,5],data=as.data.frame(methyl_access_perChrom_scaled[c(-24,-25),]))

fit_methyl_access <- lm(nullDist_freq_each_bin~E017_WGBS_FractionalMethylation.bigwig+UW.IMR90.ChromatinAccessibility.DS13229.pval.signal.bigwig,data=as.data.frame(methyl_access_perBin_scaled))

# make linear model wih 0 or 1

fit_methyl_access_perchrom_binom <- glm(chrom_gain_loss~E017_WGBS_FractionalMethylation.bigwig+UW.IMR90.ChromatinAccessibility.DS13229.pval.signal.bigwig+scale(data.frame(genes_per_chr)[,1])+scale(data.frame(chr_tbl)[,5]),data=as.data.frame(methyl_access_perChrom[c(-24,-25),]),family =binomial)

# lastly I will create a linear model with all features and use leaps function to get the 10 best models

allFeatures_DF_chroms <- cbind(epiOverlaps_perChrom_E017_noYM_noHyp, epiOverlaps_perChrom_UCSD_noYM_noHyp,methyl_access_perChrom_scaled[c(-24,-25),], chromPerMB=data.frame(chr_tbl)[,5],genes_per_chr=data.frame(genes_per_chr)[,1])

allFeatures_DF_bins <- cbind(epiOverlaps_perBin_coords_E017_nohyph, epiOverlaps_perBin_coords_IMR90_nohyph,methyl_access_perBin_scaled)

library(leaps)

leaps<-regsubsets(nullDist_freq_each_chrom~E017_DNase.pval.signal.txt+E017_H2A.Z.pval.signal.txt+E017_H2AK5ac.pval.signal.txt+E017_H2AK9ac.pval.signal.txt+E017_H2BK120ac.pval.signal.txt+E017_H2BK12ac.pval.signal.txt+E017_H2BK15ac.pval.signal.txt+E017_H2BK20ac.pval.signal.txt+E017_H2BK5ac.pval.signal.txt+E017_H3K14ac.pval.signal.txt+E017_H3K18ac.pval.signal.txt+E017_H3K23ac.pval.signal.txt+E017_H3K27ac.pval.signal.txt+E017_H3K27me3.pval.signal.txt+E017_H3K36me3.pval.signal.txt+E017_H3K4ac.pval.signal.txt+E017_H3K4me1.pval.signal.txt+E017_H3K4me2.pval.signal.txt+E017_H3K4me3.pval.signal.txt+E017_H3K56ac.pval.signal.txt+E017_H3K79me1.pval.signal.txt+E017_H3K79me2.pval.signal.txt+E017_H3K9ac.pval.signal.txt+E017_H3K9me1.pval.signal.txt+E017_H3K9me3.pval.signal.txt+E017_H4K20me1.pval.signal.txt+E017_H4K5ac.pval.signal.txt+E017_H4K8ac.pval.signal.txt+E017_H4K91ac.pval.signal.txt+UCSD.IMR90.H2A.Z.AY4.pval.signal.txt+UCSD.IMR90.H2A.Z.AY5.pval.signal.txt+UCSD.IMR90.H2AK5ac.YL20.pval.signal.txt+UCSD.IMR90.H2AK5ac.YL60.pval.signal.txt+UCSD.IMR90.H2AK9ac.AY6.pval.signal.txt+UCSD.IMR90.H2AK9ac.AY7.pval.signal.txt+UCSD.IMR90.H2BK120ac.SK15.pval.signal.txt+UCSD.IMR90.H2BK120ac.YL06.pval.signal.txt+UCSD.IMR90.H2BK12ac.YL138.pval.signal.txt+UCSD.IMR90.H2BK12ac.YL18.pval.signal.txt+UCSD.IMR90.H2BK12ac.YL68.pval.signal.txt+UCSD.IMR90.H2BK15ac.SK42.pval.signal.txt+UCSD.IMR90.H2BK15ac.YL19.pval.signal.txt+UCSD.IMR90.H2BK15ac.YL69.pval.signal.txt+UCSD.IMR90.H2BK20ac.SK14.pval.signal.txt+UCSD.IMR90.H2BK20ac.YL07.pval.signal.txt+UCSD.IMR90.H2BK5ac.LL229.pval.signal.txt+UCSD.IMR90.H2BK5ac.SK10.pval.signal.txt+UCSD.IMR90.H2BK5ac.YL281.pval.signal.txt+UCSD.IMR90.H3K14ac.SK17.pval.signal.txt+UCSD.IMR90.H3K14ac.YL15.pval.signal.txt+UCSD.IMR90.H3K18ac.LL230.pval.signal.txt+UCSD.IMR90.H3K18ac.SK09.pval.signal.txt+UCSD.IMR90.H3K23ac.SK18.pval.signal.txt+UCSD.IMR90.H3K23ac.YL09.pval.signal.txt+UCSD.IMR90.H3K27ac.LL231.pval.signal.txt+UCSD.IMR90.H3K27ac.LL235.pval.signal.txt+UCSD.IMR90.H3K27ac.YL58.pval.signal.txt+UCSD.IMR90.H3K27me3.LL223.pval.signal.txt+UCSD.IMR90.H3K27me3.SK05.pval.signal.txt+UCSD.IMR90.H3K36me3.LL225.pval.signal.txt+UCSD.IMR90.H3K36me3.YL67.pval.signal.txt+UCSD.IMR90.H3K4ac.LL239.pval.signal.txt+UCSD.IMR90.H3K4ac.SK16.pval.signal.txt+UCSD.IMR90.H3K4me1.LL222.pval.signal.txt+UCSD.IMR90.H3K4me1.SK04.pval.signal.txt+UCSD.IMR90.H3K4me1.SYL133.pval.signal.txt+UCSD.IMR90.H3K4me2.LL232.pval.signal.txt+UCSD.IMR90.H3K4me2.SK07.pval.signal.txt+UCSD.IMR90.H3K4me3.LL221.pval.signal.txt+UCSD.IMR90.H3K4me3.SK03.pval.signal.txt+UCSD.IMR90.H3K56ac.YL05.pval.signal.txt+UCSD.IMR90.H3K56ac.YL61.pval.signal.txt+UCSD.IMR90.H3K79me1.SK11.pval.signal.txt+UCSD.IMR90.H3K79me1.YL103.pval.signal.txt+UCSD.IMR90.H3K79me1.YL145.pval.signal.txt+UCSD.IMR90.H3K79me1.YL16.pval.signal.txt+UCSD.IMR90.H3K79me2.LL238.pval.signal.txt+UCSD.IMR90.H3K79me2.SK43.pval.signal.txt+UCSD.IMR90.H3K9ac.LL233.pval.signal.txt+UCSD.IMR90.H3K9ac.SK46.pval.signal.txt+UCSD.IMR90.H3K9me1.LL237.pval.signal.txt+UCSD.IMR90.H3K9me1.SK12.pval.signal.txt+UCSD.IMR90.H3K9me3.LL224.pval.signal.txt+UCSD.IMR90.H3K9me3.SK06.pval.signal.txt+UCSD.IMR90.H3K9me3.YL104.pval.signal.txt+UCSD.IMR90.H4K20me1.YL21.pval.signal.txt+UCSD.IMR90.H4K20me1.YL59.pval.signal.txt+UCSD.IMR90.H4K5ac.LL234.pval.signal.txt+UCSD.IMR90.H4K5ac.SK08.pval.signal.txt+UCSD.IMR90.H4K8ac.SK45.pval.signal.txt+UCSD.IMR90.H4K8ac.YL142.pval.signal.txt+UCSD.IMR90.H4K8ac.YL17.pval.signal.txt+UCSD.IMR90.H4K8ac.YL57.pval.signal.txt+UCSD.IMR90.H4K91ac.SK49.pval.signal.txt+UCSD.IMR90.H4K91ac.YL08.pval.signal.txt+E017_WGBS_FractionalMethylation.bigwig+UW.IMR90.ChromatinAccessibility.DS11759.pval.signal.bigwig+UW.IMR90.ChromatinAccessibility.DS11764.pval.signal.bigwig+UW.IMR90.ChromatinAccessibility.DS13219.pval.signal.bigwig+UW.IMR90.ChromatinAccessibility.DS13229.pval.signal.bigwig+UW.IMR90.Digital_Genomic_Footprinting.DGF.DS13219.pval.signal.bigwig+chromPerMB+genes_per_chr,data=as.data.frame(allFeatures_DF_chroms),nbest=10,really.big=T)
