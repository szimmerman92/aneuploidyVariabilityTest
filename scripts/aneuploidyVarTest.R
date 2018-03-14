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
for(x in 1:100) {
    chrom_18_aneuploidyFreq <- c(chrom_18_aneuploidyFreq,mean(gnk.c_sen_abs[sample(1:120, 100, replace=F),18]))
}
# plot distributions of aneu/monoploidy frequencies where you have the single nu
plot(density(nullDist_freq_each_chrom),xlim=c(0,0.14))
lines(density(chrom_18_aneuploidyFreq))
#add the frequency of chr 18 using all 120 files
abline(v=nullDist_freq_each_chrom[18],col="red")

# cannot compute ks test because we have ties, so values are not continuous
ks.test(nullDist_freq_each_chrom,chrom_18_aneuploidyFreq,alternative="less")

t.test(nullDist_freq_each_chrom,chrom_18_aneuploidyFreq,alternative="less") # p-value = 1.452e-13

wilcox.test(nullDist_freq_each_chrom,chrom_18_aneuploidyFreq,alternative="less") # p-value = 9.492e-13

### now do test where the null distribution is bootstrapped 100 times

#bootstrapping on null Dist values
set.seed(1)
nullDist_boots <- NULL
for(x in 1:100) {
    randInts <- sample(1:120, 100, replace=F)
    nullDist_boots <- c(nullDist_boots,colSums(gnk.c_sen_abs[randInts,])/length(randInts))
}
plot(density(nullDist_boots))
lines(density(chrom_18_aneuploidyFreq))
#add the frequency of chr 18 using all 120 files
abline(v=nullDist_freq_each_chrom[18],col="red")

## test all chromosomes
pVals_allChroms <- NULL
for(i in 1:23) {
    set.seed(1)
    chrom_aneuploidyFreq <- NULL
    for(x in 1:100) {
        chrom_aneuploidyFreq <- c(chrom_aneuploidyFreq,mean(gnk.c_sen_abs[sample(1:120, 100, replace=F),i]))
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
