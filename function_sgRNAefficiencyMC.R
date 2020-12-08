#### Function to detect and score binding sites for given sgRNAs in any NCBI genome ####
# Author: Vincent de Bakker
# Veening Lab, DMF, FBM, University of Lausanne, Switzerland
# vincent.debakker@unil.ch
####

sgRNAefficiencyMC <- function(sgRNAs, genes, genome, 
                              reprAct = TRUE, dist2SC = TRUE, 
                              name_by = "locus_tag", 
                              penalties = c("qi.mean.per.region", "qi", "custom"), 
                              custom.penalties = NULL, 
                              outfile = getwd(), 
                              no_cores = 1, 
                              ...){
  require(CRISPRseek)
  require(parallel)
  
  # Retrieve all gene names from GFF input
  tags <- unlist(lapply(genes$attribute, function(x){
    # add ";" at end in case locus_tag is last attribute (regular expression requires ending character)
    sub(paste0(".*?", name_by, "=(.*?);.*"), "\\1", paste0(x, ";"))
  }))
  
  # Checks
  if(any(duplicated(tags))){
    stop(paste("Duplicates of", name_by, "in genes."))
  }
  if(reprAct){ # set penalties if repression activity should be computed
    qi <- c(0.78, 0.41, 0.44, 0.63, 0.65, 0.68, 0.65, 0.63, # III (inverse counting Qi vs CRISPRseek)
            0.30, 0.25, 0.24, 0.22, 0.24,                   # II
            0.01, 0.12, 0.10, 0.07, 0.06, 0.09, 0.05)       # I
    # tapply sorts by names, so sort by 1:3 (CRISPRseek), while I:III is Qi et al regions
    qi.mean.per.region <- rep(tapply(qi, rep(c("1III", "2II", "3I"), c(8, 5, 7)), mean), c(8, 5, 7))
    penalties <- switch(penalties[1], 
                        qi = qi, 
                        qi.mean.per.region = qi.mean.per.region, 
                        custom = custom.penalties)
    if(length(penalties) != 20 | !is.numeric(penalties)){
      stop("For custom penalties, please set custom.penalties to a numeric vector of length 20")
    }
  }
  
  # Check cores specification for parallelization
  if((!is.integer(no_cores) & !is.numeric(no_cores)) | length(no_cores) != 1 | sign(no_cores) != 1){
    stop("Please set no_cores to an integer of length 1.")
  }
  
  # Binding site identification
  if(no_cores != 1){ # Parallelize or not as appropriate
    message("Making cluster...")
    # Start local socket cluster 
    cl <- makeCluster(no_cores)
    # export required objects to every core
    clusterExport(cl, c("sgRNAs", 
                        "genes", # for NT_gene, later
                        "tags",  # for NT_gene, later
                        unlist(lapply(match.call(expand.dots = FALSE)$..., all.vars)), 
                        "genome", 
                        "outfile"), 
                  envir = environment()) # from within-function environment
    # load required packages on every core
    clusterEvalQ(cl, library(CRISPRseek))
    # Identify all binding sites per sgRNA
    message("Starting parallel binding site identification...")
    ## per sgRNA, per genome
    lib_hits <- do.call(rbind, parLapply(cl, names(sgRNAs), function(y){
      do.call(rbind, lapply(seq.int(genome), function(x){
        searchHits(gRNAs = sgRNAs[names(sgRNAs) == y],
                   seqs = genome[[x]], seqname = names(genome)[x], 
                   outfile = paste0(outfile, 
                                    "_sgRNA-", y, "_genome-", names(genome)[x]), 
                   ...)
      }))
    }))
  } else{ # Parallelize or not as appropriate
    # Identify all binding sites
    message("Starting binding site identification...")
    lib_hits <- do.call(rbind, lapply(seq.int(genome), function(x){
      searchHits(gRNAs = sgRNAs,
                 seqs = genome[[x]], seqname = names(genome)[x], 
                 outfile = paste0(outfile, 
                                  "_genome-", names(genome)[x]), 
                 ...)
    }))
  }
  
  # Get NT CDS per hit
  # identify per site in which gene it lies on NT strand, if any
  #    if gene on + strand, sgRNA seq should be on - strand
  #    then complementary to coding = NT strand
  message("Identifying gene hits per sgRNA...")
  # NB force list outcome with lapply, do not use apply
  if(no_cores != 1){ # if parallel
    NT_gene <- parLapply(cl, seq.int(nrow(lib_hits)), function(x){
      # find overlapping genes per site, if any
      tmp_i <- as.numeric(lib_hits[x, "chromStart"]) - genes$end <= 0 & 
        as.numeric(lib_hits[x, "chromEnd"]) - genes$start >= 0 & 
        as.character(lib_hits[x, "strand"]) != genes$strand & 
        # also look on the right chromosome
        unlist(lapply(genes$seqid, grepl, x = as.character(lib_hits[x, "chrom"])))
      # add details of gene hits, if any
      if(any(tmp_i)){
        # could be multiple genes that are overlapped by a site
        do.call(rbind, lapply(which(tmp_i), function(y){
          if(as.numeric(lib_hits[x, "chromStart"]) < genes$start[y]){
            coverPart <- switch(genes$strand[y], 
                                `+` = "5-tail", 
                                `-` = "3-PAM")
            coverSize <- as.numeric(lib_hits[x, "chromEnd"]) - genes$start[y] + 1
          } else{
            if(as.numeric(lib_hits[x, "chromEnd"]) > genes$end[y]){
              coverPart <- switch(genes$strand[y], 
                                  `+` = "3-PAM", 
                                  `-` = "5-tail")
              coverSize <- genes$end[y] - as.numeric(lib_hits[x, "chromStart"]) + 1
            } else{
              coverPart <- "complete"
              coverSize <- 23
            }
          }
          # return details if gene hits for site
          matrix(c(tags[y], coverPart, coverSize), ncol = 3)
        }))
      } else{
        # return NAs if no gene hits for site
        matrix(NA, ncol = 3, nrow = 1)
      }
    })
  } else{ # if not parallel
    NT_gene <- lapply(seq.int(nrow(lib_hits)), function(x){
      # find overlapping genes per site, if any
      tmp_i <- as.numeric(lib_hits[x, "chromStart"]) - genes$end <= 0 & 
        as.numeric(lib_hits[x, "chromEnd"]) - genes$start >= 0 & 
        as.character(lib_hits[x, "strand"]) != genes$strand & 
        # also look on the right chromosome
        unlist(lapply(genes$seqid, grepl, x = as.character(lib_hits[x, "chrom"])))
      # add details of gene hits, if any
      if(any(tmp_i)){
        # could be multiple genes that are overlapped by a site
        do.call(rbind, lapply(which(tmp_i), function(y){
          if(as.numeric(lib_hits[x, "chromStart"]) < genes$start[y]){
            coverPart <- switch(genes$strand[y], 
                                `+` = "5-tail", 
                                `-` = "3-PAM")
            coverSize <- as.numeric(lib_hits[x, "chromEnd"]) - genes$start[y] + 1
          } else{
            if(as.numeric(lib_hits[x, "chromEnd"]) > genes$end[y]){
              coverPart <- switch(genes$strand[y], 
                                  `+` = "3-PAM", 
                                  `-` = "5-tail")
              coverSize <- genes$end[y] - as.numeric(lib_hits[x, "chromStart"]) + 1
            } else{
              coverPart <- "complete"
              coverSize <- 23
            }
          }
          # return details if gene hits for site
          matrix(c(tags[y], coverPart, coverSize), ncol = 3)
        }))
      } else{
        # return NAs if no gene hits for site
        matrix(NA, ncol = 3, nrow = 1)
      }
    })
  }
  # repeat row indices with x number of gene hits x times
  NT_gene_i <- unlist(mapply(rep, 
                             seq.int(nrow(lib_hits)), 
                             unlist(lapply(NT_gene, nrow))))
  # add details to data frame
  lib_hits <- cbind(lib_hits[NT_gene_i, ], 
                    "site" = NT_gene_i, 
                    do.call(rbind, NT_gene))
  colnames(lib_hits)[32:34] <- c("NTgene", "coverPart", "coverSize")
  # coverSize should be integer, not factor
  lib_hits$coverSize <- as.integer(as.character(lib_hits$coverSize))
  
  # Get repression activity
  if(reprAct){
    message("Calculating repression activity per site...")
    if(no_cores != 1){
      clusterExport(cl, c("penalties"), envir = environment())
      # Get retained repression per hit
      repr_act <- parApply(cl, lib_hits[, 1:20], 1, function(x){
        prod(penalties[x == 1]) # product of empty set is 1 by definition (OK for 0 mm)
      })
    } else{
      # Get retained repression per hit
      repr_act <- apply(lib_hits[, 1:20], 1, function(x){
        prod(penalties[x == 1]) # product of empty set is 1 by definition (OK for 0 mm)
      })
    }
    # add to data frame
    lib_hits$reprAct <- repr_act
  }
  
  # Get distance to start codon if sgRNA can (partially) bind within gene
  if(dist2SC){
    message("Calculating relative distances of sgRNA binding within gene hits...")
    # Get distance to start codon per hit
    #    normalize within-gene distances with feature scaling
    feat_scale <- function(x, min, max){(x - min) / (max - min)}
    if(no_cores != 1){ # if parallel
      dist_SC <- parApply(cl, lib_hits, 1, function(x){
        if(!is.na(x["NTgene"])){
          tmp_range <- genes[match(x["NTgene"], tags), c("start", "end")]
          # max is end CDS - 22 nt for binding first nt of PAM+spacer (chromStart, strand does not matter)
          #    also partial overlap (coverPart != "complete") dist_SC should always be set to [0,1]
          ifelse(genes$strand[match(x["NTgene"], tags)] == "+", 
                 # strand matters: on "-" strand, start of gene is on 3-prime end so invert distance
                 pmax(0, pmin(1, 
                              feat_scale(as.numeric(x["chromStart"]), 
                                         tmp_range[1], 
                                         tmp_range[2] - 22))), 
                 1 - pmax(0, pmin(1, 
                                  feat_scale(as.numeric(x["chromStart"]), 
                                             tmp_range[1], 
                                             tmp_range[2] - 22))))
        } else{
          NA
        }
      }) 
    } else{ # if not parallel
      dist_SC <- apply(lib_hits, 1, function(x){
        if(!is.na(x["NTgene"])){
          tmp_range <- genes[match(x["NTgene"], tags), c("start", "end")]
          # max is end CDS - 22 nt for binding first nt of PAM+spacer (chromStart, strand does not matter)
          #    also partial overlap (coverPart != "complete") dist_SC should always be set to [0,1]
          ifelse(genes$strand[match(x["NTgene"], tags)] == "+", 
                 # strand matters: on "-" strand, start of gene is on 3-prime end so invert distance
                 pmax(0, pmin(1, 
                              feat_scale(as.numeric(x["chromStart"]), 
                                         tmp_range[1], 
                                         tmp_range[2] - 22))), 
                 1 - pmax(0, pmin(1, 
                                  feat_scale(as.numeric(x["chromStart"]), 
                                             tmp_range[1], 
                                             tmp_range[2] - 22))))
        } else{
          NA
        }
      })
    }
    # add to data frame
    lib_hits$dist2SC <- unlist(dist_SC)
  }
  
  # stop cluster if required
  if(no_cores != 1){
    message("Closing cluster")
    stopCluster(cl)
  }
  
  # Return full data frame
  return(lib_hits)
}