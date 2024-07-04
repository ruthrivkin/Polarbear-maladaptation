setwd("/Users/ruthrivkin/Dropbox/Postdoc_2021-2024/Polar_Bears/Species_Dist_Modelling/Gradient_Forest/Subset/")
library(ggplot2)
ng <- theme(aspect.ratio=0.7,panel.background = element_blank(), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.border=element_rect(fill = NA, color = "black"),
            axis.line = element_line(linewidth=1), 
            axis.ticks=element_line(color="black"),
            axis.text=element_text(color="black",size=12, margin = 0.5), 
            axis.title=element_text(color="black",size=1), 
            axis.title.y=element_text(vjust=1,face="plain",size=15),
            axis.title.x=element_text(vjust=0.1,face="plain",size=15),
            axis.text.x=element_text(size=15),
            axis.text.y=element_text(size=15))

my.colors = colorRampPalette(c("#5E85B8","#EDF0C0","#C13127"))

#Calculate basic stats
library(hierfstat)
library(dartR)
library(StAMPP)


#load data
infile <- "../Subset/Subset.ld.hwe.raw"
myData <- read.PLINK(infile)
myData.nona <- gl.impute(myData, method = "HW")
myData.nona@pop


#Fst
fst <- stamppFst(myData.nona, nboots = 100, percent = 95, nclusters = 1)
summary(fst)
fst.matrix <- as.data.frame(fst$Fsts) #pvalues = 0 except for kb/bb = 0.12
pvalues <- fst$Pvalues

write.csv(fst.matrix, "fst.ld.matrix.csv")

#Calculate private alleles
myData <- gl.compliance.check(myData)
pa <- gl.report.pa(myData, method = "one2rest")
write.csv(pa, "privatealleles.ld.csv")

#Heterozygosity

#Plink
hetdata <- read.delim("../Subset/Subset.ld.hwe.het", sep = "")
str(hetdata)

#Prep het file
het <- hetdata %>% rename(subpop = FID)
head(het)

het <- het %>%
  mutate(
    prob.obs.hom = O.HOM./N.NM.,
    prob.exp.hom = E.HOM./N.NM.,
    prob.obs.het = 1 - prob.obs.hom,
    prob.exp.het = 1 - prob.exp.hom
  )
head(het)

write.csv(het, "Geneticdiversity_ld.csv")

het <- read.csv("Geneticdiversity_ld.csv")


meanhet <- Rmisc::summarySE(nopar, measurevar = "prob.obs.het", groupvars = "subpop", na.rm = TRUE)
meanhet

meanF <- Rmisc::summarySE(het, measurevar = "F", groupvars = "subpop", na.rm = TRUE)
meanF

#Heterozygosity
het <- read.csv("Geneticdiversity_ld.csv")
cols = my.colors(length(unique(het$subpop)))

(hetplot.subpop <- ggplot(het, aes(x=reorder(subpop, prob.obs.het, na.rm = T), y=prob.obs.het, )) + 
    geom_boxplot(outlier.size = 1, color = cols) + 
    geom_dotplot(binaxis='y', stackdir ='center', dotsize=0.01) +
    geom_jitter(shape=1, position=position_jitter(0.2)) + 
    xlab("Subpopulation") +
    ylab("Observed Heterozygosity") +
    scale_y_continuous(limits = c(0.15,0.3)) +
    ng
)
ggsave("Heterozygosity_by_subpop_hwe.ld.pdf", width = 6.78, height = 5.3, dpi = 300)

(f.subpop <- ggplot(het, aes(x=reorder(subpop, F, na.rm = T), y=F, )) + 
    geom_boxplot(outlier.size = 1, color = cols) + 
    geom_dotplot(binaxis='y', stackdir ='center', dotsize=0.01) +
    geom_jitter(shape=1, position=position_jitter(0.2)) + 
    xlab("Subpopulation") +
    ylab("Inbreeding") +
    scale_y_continuous(limits = c(-0.05,0.4)) +
    
    ng
)
ggsave("F_by_subpop_hwe.ld.pdf", width = 6.78, height = 5.3, dpi = 300)


#hierfstat

mydatahr <- gl2gi(myData, probar = FALSE, verbose = NULL)

popgen <- basic.stats(mydatahr)
popHo <- popgen$Ho 

meanHo <- colMeans(popHo) #These are the same as the plink results
meanHo <- reshape2::melt(meanHo, id.vars=c("BB","DS",        "FB",        "GB",        "KB",        "LS",        "MC",        "NB",        "NW",        "SB",        "SH",        "VM",        "WH"),
                         variable.name = "subpop", 
                         value.name = "Ho")
meanHo$subpop <- row.names(meanHo)
barplot(meanHo$Ho ~ meanHo$subpop) 


#Take a look at loci out of hwe
hwe <- hw.test(mydatahr)
hwe <- as.data.frame(hwe)
head(hwe)
#summarize loci that are significant, rename rows and move loci names to column
names(hwe) <- c("chi", "n", "pchi", "pmc")
hwe$snp <- row.names(hwe)

#summarize by chi-based p value
hwe.loc  <- filter(hwe,  hwe$pmc < 0.001)   
length(hwe.loc$snp) ## 22 snps are out of hwe
hwe.loc$snp

# summarize across pops
HWE.test <- data.frame(sapply(seppop(mydatahr), 
                              function(ls) pegas::hw.test(ls, B=0)[,3]))
HWE.test.chisq <- t(data.matrix(HWE.test))
#summarize the proportion of populations where loci are out of HWE.
alpha=0.001
Prop.loci.out.of.HWE <- data.frame(Chisq=apply(HWE.test.chisq<alpha, 1, mean))
Prop.loci.out.of.HWE  #each pop has about 5% of loci out of HWE

#Check again using FDR for multiple tests
Chisq.fdr <- matrix(p.adjust(HWE.test.chisq,method="fdr"), 
                    nrow=nrow(HWE.test.chisq))

Prop.pops.out.of.HWE <- data.frame(Chisq=apply(HWE.test.chisq<alpha, 1, mean), 
                                   Chisq.fdr=apply(Chisq.fdr<alpha, 1, mean))
Prop.pops.out.of.HWE #<1% for all pops


#Calcuclate nucleotide diversity

library(radiator)
library(tidyr)
library(readr)

#Make strata file with ind id and subset info (POP_ID and INDIVIDUALS)
vcf <- read_vcf(data = "Subset.ld.hwe.vcf", strata = "../GFAllEnvironments/strata.filtered.tsv")
pi <- pi(vcf)

pi$pi.populations

pidf <- as.data.frame(pi$pi.individuals)
pidf2 <- separate_wider_delim(pidf, cols = INDIVIDUALS, delim = "-", names = c("subpop", "IID"))

het.div.data <- merge(het, pidf2, by = c("IID", "subpop"))
het.div.data <- select(het.div.data,-c(3,12))
write.csv(het.div.data, "Geneticdiversity_ldpruned.csv")

div <- read.csv("Geneticdiversity_hwefiltered.csv")
str(div)

meanpi <- Rmisc::summarySE(div, measurevar = "PI", groupvars = "subpop", na.rm = TRUE)
meanpi


cols = my.colors(length(unique(div$subpop)))

(pi.subpop <- ggplot(div, aes(x=reorder(subpop, PI, na.rm = T), y=PI, )) + 
    geom_boxplot(outlier.size = 1, color = cols) + 
    geom_dotplot(binaxis='y', stackdir ='center', dotsize=0.01) +
    geom_jitter(shape=1, position=position_jitter(0.2)) + 
    xlab("Subpopulation") +
    ylab("Nucleotide Diversity (pi)") +
    ng
)
ggsave("Heterozygosity_by_subpop.pdf", width = 6.78, height = 5.3, dpi = 300)


#Compare transcriptomic vs rad loci
#Plink
hetdata.t <- read.delim("../Subset/Subset.trans.het", sep = "")
str(hetdata.t)

hetdata.rad <- read.delim("../Subset/Subset.RAD.het", sep = "")
str(hetdata.rad)

library(dplyr)

tran <- hetdata.t %>% rename(subpop = FID)
head(tran)

tran <- tran %>%
  mutate(Ho.tran = 1 - (O.HOM./N.NM.))  %>%
  rename(F.tran = F)

head(tran)

#Prep het file
rad <- hetdata.rad %>% rename(subpop = FID)

rad <- rad %>%
  mutate(Ho.rad = 1 - (O.HOM./N.NM.))  %>%
  rename(F.rad = F)

head(rad)

compare <- inner_join(tran, rad, by = "IID") 

compare <- dplyr::select(compare, c("subpop.x", "IID", "Ho.tran", "Ho.rad" ))



meanhet.rad <- Rmisc::summarySE(compare, measurevar = "Ho.rad", groupvars = "subpop.x", na.rm = TRUE)
meanhet.rad

meanhet.tran <- Rmisc::summarySE(compare, measurevar = "Ho.tran", groupvars = "subpop.x", na.rm = TRUE)
meanhet.tran


#Calcuclate nucleotide diversity

library(radiator)
library(tidyr)
library(readr)

#Make strata file with ind id and subset info (POP_ID and INDIVIDUALS)
radvcf <- read_vcf(data = "../Subset/Subset.RAD.vcf", strata = "../GFAllEnvironments/strata.filtered.tsv")
transvcf <- read_vcf(data = "../Subset/Subset.trans.vcf", strata = "../GFAllEnvironments/strata.filtered.tsv")

pirad <- pi(radvcf)

pitrans <- pi(transvcf)

pirad$pi.populations

pitrans$pi.populations


