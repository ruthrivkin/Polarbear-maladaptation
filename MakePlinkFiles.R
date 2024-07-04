setwd("/Users/ruthrivkin/Dropbox/Postdoc_2021-2024/Polar_Bears/Species_Dist_Modelling/SourceFiles/")

#Prune snps in LD, this is a two part process
#1) Identify snps, this makes a file with snps that are out of ld to exclude (.prune.out)
system("~/plink/plink --ped chip12.ped --map chip12.map --out forgf --geno 0.1 --maf 0.01 --indep-pairwise 10 1 0.1 --allow-extra-chr --recode --keep Subset_IDs.txt")

#2) Remove snps out of ld for gradient forest analysis
system("~/plink/plink --ped forgf.ped --map forgf.map --out forgf --missing --het --fst --recode  --allow-extra-chr --exclude forgf.prune.out --keep Subset_IDs.txt")
#recode A for .raw file and recode vcf for vcf


#Remove snps that are also out of HWE (used for gen diversity and differentiation analysis)
system("~/plink/plink --ped forgf.ped --map forgf.map --out Subset.ld.hwe --missing --het --fst --recode --hwe 0.001 midp --allow-extra-chr --exclude forgf.prune.out --keep Subset_IDs.txt")
#recode A for .raw file and recode vcf for vcf


#Calc genetic diversity and differentiation for RAD and transcriptomic SNPs
#Tanscriptomic SNP IDs determined from orignal map/ped files (eg: transcriptomic snps all start with scaffold_XX)
system("~/plink/plink --ped Subset.ld.hwe.ped --map Subset.ld.hwe.map --out Subset.trans --allow-extra-chr --missing --het --recode vcf --fst --extract TranscriptomeSNPS.txt --keep Subset_IDs.txt --within Within_Subset_IDs.txt")

system("~/plink/plink --ped Subset.ld.hwe.ped --map Subset.ld.hwe.map --out Subset.rad --allow-extra-chr --missing --het --recode vcf --fst --exclude TranscriptomeSNPS.txt --keep Subset_IDs.txt --within Within_Subset_IDs.txt")


