#' phenoscanner
#'
#' phenoscanner queries the PhenoScanner database of genotype-phenotype associations from inside R.
#' @param snpquery a vector of SNPs.
#' @param regionquery a vector of genomic regions.
#' @param genequery a vector of gene names.
#' @param catalogue the catalogue to be searched (options: None, GWAS, eQTL, pQTl, mQTL, methQTL).
#' @param pvalue the p-value threshold.
#' @param proxies the proxies database to be searched (options: None, AFR, AMR, EAS, EUR, SAS).
#' @param r2 the r2 threshold. 
#' @param build the genome build (options: 37, 38).
#' @import rjson
#' @return a list containing a data.frame of association results and a data.frame of SNP/Region/Gene information from PhenoScanner.
#' @examples
#' # SNP
#' res <- phenoscanner(snpquery="rs10840293")
#' head(res$results)
#' res$snps
#' 
#' # Gene
#' res <- phenoscanner(genequery="SWAP70")
#' head(res$results)
#' res$snps
#' 
#' # Region
#' res <- phenoscanner(regionquery="chr11:9685624-9774538")
#' head(res$results)
#' res$regions
#' @author PhenoScanner <phenoscanner@gmail.com>
#' @export

phenoscanner <- function(snpquery=NULL, genequery=NULL, regionquery=NULL, catalogue="GWAS", pvalue=1E-5, proxies="None", r2=0.8, build=37){
  cat("PhenoScanner V2\n")
  cat("Cardiovascular Epidemiology Unit\n")
  cat("University of Cambridge\n")
  cat("Email: phenoscanner@gmail.com\n\n")
  cat("Information: Each user can query a maximum of 10,000 SNPs (in batches of 100), 1,000 genes (in batches of 10) or 1,000 regions (in batches of 10) per hour. For large batch queries, please ask to download the data from www.phenoscanner.medschl.cam.ac.uk/data.\n")
  cat("Terms of use: Please refer to the terms of use when using PhenoScanner V2 (www.phenoscanner.medschl.cam.ac.uk/about). If you use the results from PhenoScanner in a publication or presentation, please cite all of the relevant references of the data used and the PhenoScanner publications: www.phenoscanner.medschl.cam.ac.uk/about/#citation.\n\n")
  if(is.null(snpquery) & is.null(regionquery) & is.null(genequery)) stop("no query has been requested")
  if((length(snpquery[1])+length(regionquery[1])+length(genequery[1]))>1) stop("only one query type allowed")
  if(!(catalogue=="None" | catalogue=="GWAS" | catalogue=="eQTL" | catalogue=="pQTL" | catalogue=="mQTL" | catalogue=="methQTL")) stop("catalogue has to be one of None, GWAS, eQTL, pQTL, mQTL or methQTL")
  if(!(proxies=="None" | proxies=="AFR" | proxies=="AMR" | proxies=="EAS" | proxies=="EUR" | proxies=="SAS")) stop("proxies has to be one of None, AFR, AMR, EAS, EUR or SAS")
  if(length(snpquery)>100) stop("a maximum of 100 SNP queries can be requested at one time")
  if(length(genequery)>10) stop("a maximum of 10 gene queries can be requested at one time")
  if(length(regionquery)>10) stop("a maximum of 10 region queries can be requested at one time")
  if(!(pvalue>0 & pvalue<=1)) stop("the p-value threshold has to be greater than 0 and less than or equal to 1")
  if(!(r2>=0.5 & r2<=1)) stop("the r2 threshold has to be greater than or equal to 0.5 and less than or equal to 1")
  if(!(build==37 | build==38)) stop("the build has to be equal to 37 or 38")
  if(!is.null(regionquery)){
    ub <- as.numeric(sub(".*-", "", sub(".*:", "",regionquery)))
    lb <- as.numeric(sub("-.*", "", sub(".*:", "",regionquery)))
    dist <- ub - lb
    if(any(dist>1000000)) stop("region query can be maximum of 1MB in size")
  }
  if(length(snpquery)>0){
    results <- data.frame()
    snps <- data.frame()
    n_queries <- length(snpquery) %/% 10
    if((length(snpquery) %% 10)>0){n_queries <- n_queries + 1}
    for(i in 1:n_queries){
      if(i < n_queries){qsnps <- paste0(snpquery[(1+10*(i-1)):(10*i)], collapse="+")}else{qsnps <- paste0(snpquery[(1+10*(i-1)):length(snpquery)], collapse="+")}
      json_file <- paste0("http://www.phenoscanner.medschl.cam.ac.uk/api/?snpquery=",qsnps,"&catalogue=",catalogue,"&p=",pvalue,"&proxies=",proxies,"&r2=",r2,"&build=",build)
      json_data <- rjson::fromJSON(file=json_file)
      if(length(json_data$results)==0 & length(json_data$snps)==0){
        print(paste0("Error: ",json_data$error))
        next
      }
      if(length(json_data$results)>0){
        fields <- json_data$results[[1]]; json_data$results[[1]] <- NULL
        if(length(json_data$results)>0){
          tables <- as.data.frame(matrix(unlist(json_data$results), ncol=length(fields), byrow=T), stringsAsFactors=F)
          names(tables) <- fields
          results <- rbind(results,tables)
          if(length(snpquery)==1){print(paste0(snpquery," -- queried"))}else{print(paste0(i," -- chunk of 10 SNPs queried"))}
        }else{if(length(snpquery)==1){print(paste0("Warning: no results found for ",snpquery))}else{print(paste0("Warning: no results found for chunk ",i))}}
      }
      if(length(json_data$snps)>0){
        fields_snps <- json_data$snps[[1]]; json_data$snps[[1]] <- NULL
        if(length(json_data$snps)>0){
          tables_snps <- as.data.frame(matrix(unlist(json_data$snps), ncol=length(fields_snps), byrow=T), stringsAsFactors=F)
          names(tables_snps) <- fields_snps
          snps <- rbind(snps,tables_snps)
          if(length(json_data$results)==0){if(length(snpquery)==1){print(paste0(snpquery," -- queried"))}else{print(paste0(i," -- chunk of 10 SNPs queried"))}}
        }else{if(length(json_data$results)==0){if(length(snpquery)==1){print(paste0("Warning: no results found for ",snpquery))}else{print(paste0("Warning: no results found for chunk ",i))}}}
      }
    }
    output <- list(snps=snps, results=results)
  }
  if(length(genequery)>0){
    results <- data.frame()
    genes <- data.frame()
    n_queries <- length(genequery)
    for(i in 1:n_queries){
      json_file <- paste0("http://www.phenoscanner.medschl.cam.ac.uk/api/?genequery=",genequery[i],"&catalogue=",catalogue,"&p=",pvalue,"&proxies=None&r2=1&build=",build)
      json_data <- rjson::fromJSON(file=json_file)
      if(length(json_data$results)==0 & length(json_data$genes)==0){
        print(paste0("Error: ", json_data$error))
        next
      }
      if(length(json_data$results)>0){
        fields <- json_data$results[[1]]; json_data$results[[1]] <- NULL
        if(length(json_data$results)>0){
          tables <- as.data.frame(matrix(unlist(json_data$results), ncol=length(fields), byrow=T), stringsAsFactors=F)
          names(tables) <- fields
          results <- rbind(results,tables)
          print(paste0(genequery," -- queried"))
        }else{print(paste0("Warning: no results found for ",genequery))}
      }
      if(length(json_data$genes)>0){
        fields_genes <- json_data$genes[[1]]; json_data$genes[[1]] <- NULL
        if(length(json_data$genes)>0){
          tables_genes <- as.data.frame(matrix(unlist(json_data$genes), ncol=length(fields_genes), byrow=T), stringsAsFactors=F)
          names(tables_genes) <- fields_genes
          genes <- rbind(genes,tables_genes)
        }
      }
    }
    output <- list(genes=genes, results=results)
  }
  if(length(regionquery)>0){
    results <- data.frame()
    regions <- data.frame()
    n_queries <- length(regionquery)
    for(i in 1:n_queries){
      json_file <- paste0("http://www.phenoscanner.medschl.cam.ac.uk/api/?regionquery=",regionquery[i],"&catalogue=",catalogue,"&p=",pvalue,"&proxies=None&r2=1&build=",build)
      json_data <- rjson::fromJSON(file=json_file)
      if(length(json_data$results)==0 & length(json_data$locations)==0){
        print(paste0("Error: ",json_data$error))
        next
      }
      if(length(json_data$results)>0){
        fields <- json_data$results[[1]]; json_data$results[[1]] <- NULL
        if(length(json_data$results)>0){
          tables <- as.data.frame(matrix(unlist(json_data$results), ncol=length(fields), byrow=T), stringsAsFactors=F)
          names(tables) <- fields
          results <- rbind(results,tables)
          print(paste0(regionquery," -- queried"))
        }else{print(paste0("Warning: no results found for ",regionquery))}
      }
      if(length(json_data$locations)>0){
        fields_regions <- json_data$locations[[1]]; json_data$locations[[1]] <- NULL
        if(length(json_data$locations)>0){
          tables_regions <- as.data.frame(matrix(unlist(json_data$locations), ncol=length(fields_regions), byrow=T), stringsAsFactors=F)
          names(tables_regions) <- fields_regions
          regions <- rbind(regions,tables_regions)
        }
      }
    }
    output <- list(regions=regions, results=results)
  }
  if(is.null(output)) stop("there is no output")
  return(output)
}
