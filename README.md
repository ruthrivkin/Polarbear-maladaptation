---
title: "Assessing the risk of climate maladaptation for Canadian polar bears"
output: html_document
---

Assessing the risk of climate maladaptation for Canadian polar bears (doi: 10.6084/m9.figshare.25206887). This is a repository for R scripts used for performing conducting the analyses. The initial input plink genotype files can be accessed through figshare (https://figshare.com/s/543c6be97a1d1e6957e6). The description of the files  included in the repository are included below.	


1. Chip12.ped: Genotype ped file for 411 polar bears included in analysis (accessible on figshare)
2. Chip12.map: Corresponding map file from the ped. Both are required for PLINK analyses (accessible on figshare)
3. MakePlinkFiles.R:Scripts for generating plink files used in the analyses
4. Environmental_Data_Aquisition.R: R code for extracting current and future environmental data
5. PopGen.R: R code for population genetic diversity and differentiation analyses
6. snmf.R: R code for snmf ancestry distribution analysis
7. sPCA.R: R code for sPCA population structure analsysis
8. ConStruct.R: R code for conStruct analysys
9. conStruct.sh: Shell executable file for runnning conStruct
10. conStruct-xval.sh: Shell executable file for runnning conStruct cross-validation
11. AllRCP-GFcode.R: R code for Gradient Forest analysis
12. PopulationAnalysis.R: R code for genetic offset population analyses	
