# phenoscanner
phenoscanner allows users to query the PhenoScanner database of genotype-phenotype associations from inside R.

# Functions
* phenoscanner - this function allows users to query the PhenoScanner database from inside R 

# Installation
1. install.packages("devtools")
2. library(devtools) 
3. install_github("phenoscanner/phenoscanner")
4. library(phenoscanner)

# Examples 
\# SNP  
res <- phenoscanner(snpquery="rs10840293")  
head(res$results)  
res$snps  

\# Gene  
res <- phenoscanner(genequery="SWAP70")  
head(res$results)  
res$snps  

\# Region  
res <- phenoscanner(regionquery="chr11:9685624-9774538")  
head(res$results)  
res$regions  
