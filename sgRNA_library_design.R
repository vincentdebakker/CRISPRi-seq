#### (Optimal) sgRNA design per annotated feature in any genome ####
# Author: Vincent de Bakker
# Veening Lab, DMF, FBM, University of Lausanne, Switzerland
# vincent.debakker@unil.ch
####

#### 1. Settings ####
input_genome <- # .gb[ff] file or NCBI assembly accession nr.
fundir <- # path to directory with file function_sgRNAefficiencyMC.R
outdir <- # path to output directory
n_cores <- 3 # number of cores to be used
#
path_ncbi_downloads <- "~" # if input_genome is accession nr.: path to directory to save genome & annotation
feature_type <- "locus_tag"
max_mismatch <- 6
reprAct_penalties <- "qi.mean.per.region"
reprAct_custom_penalties <- NULL
errorRange_maxOffreprAct <- 0.01
oligoForwardOverhang <- "TATA"
oligoReverseOverhang <- "AAAC"
spacer_length <- 20
PAM <- "NGG"
output_target_fasta <- TRUE
output_full_list <- TRUE
output_optimized_list <- TRUE


#### 2. Preliminaries ####
# check inputs
if(endsWith(input_genome, ".gb")){
  input_type <- "gbfile"
} else{
  if(any(startsWith(input_genome, c("GCA", "GCF")))){
    input_type <- "accessionnr"
  } else{
    stop("Input file should either be a .gb file or an NCBI genome assembly accession number (GCA_ or GCF_).")
  }
}
if(!file.exists(paste0(fundir, "/function_sgRNAefficiencyMC.R"))){
  stop("Please save the function file function_sgRNAefficiencyMC.R in the directory 'fundir'.")
}
if(sum(output_full_list, output_optimized_list, output_target_fasta) == 0){
  stop("You have currently selected no output. Please set at least one of the output options to TRUE.")
}
if(n_cores < 1 | !is.numeric(n_cores)){stop("n_cores should be an integer of 1 or higher.")}
if(input_type == "accessionnr"){
  db <- switch(substr(input_genome, 1, 3), 
               GCA = "genbank", 
               GCF = "refseq")
}
# detect OS
platform <- .Platform$OS.type
# install required packages if needed
required_packages_CRAN <- c("BiocManager")
required_packages_BioC <- switch(input_type, 
                                 gbfile = c("Biostrings", "genbankr", "CRISPRseek"), 
                                 accessionnr = c("Biostrings", "biomartr", "CRISPRseek"))
for(i in seq.int(required_packages_CRAN)){
  if(!requireNamespace(required_packages_CRAN[i], quietly = TRUE)){install.packages(required_packages_CRAN[i])}
}
for(i in seq.int(required_packages_BioC)){
  if(!requireNamespace(required_packages_BioC[i], quietly = TRUE)){BiocManager::install(required_packages_BioC[i])}
}
# load packages
lapply(required_packages_BioC, library, character.only = TRUE)
library(parallel)
# load user-defined functions
source(paste0(fundir, "/function_sgRNAefficiencyMC.R"))


#### 3. Get genome, features and extract feature sequences ####
if(input_type == "accessionnr"){
  if(file.exists(paste0(path_ncbi_downloads, "/_ncbi_downloads/genomes/", input_genome,"_genomic_", db, ".fna.gz"))){
    genome_path <- paste0(path_ncbi_downloads, "/_ncbi_downloads/genomes/", input_genome,"_genomic_", db, ".fna.gz")
  } else{
    genome_path <- getGenome(db = db, 
                             organism = input_genome, 
                             reference = FALSE, 
                             path = paste0(path_ncbi_downloads, "/_ncbi_downloads/genomes/"))
  }
  genome <- read_genome(genome_path)
  if(file.exists(paste0(path_ncbi_downloads, "/_ncbi_downloads/annotation/", input_genome, "_genomic_", db, ".gff.gz"))){
    gffPath <- paste0(path_ncbi_downloads, "/_ncbi_downloads/annotation/", input_genome, "_genomic_", db, ".gff.gz")
  } else{
    gffPath <- getGFF(db = db, 
                      organism = input_genome, 
                      reference = FALSE, 
                      path = paste0(path_ncbi_downloads, "/_ncbi_downloads/annotation/"))
  }
  GFF <- read_gff(gffPath)
  # extract all features of type feature_type
  genes <- GFF[unlist(lapply(GFF$attribute, grepl, pattern = feature_type)), ]
  # replace split features with same attributes by first with total range for start and end
  #   e.g. annotions of "joined feature span" 
  #       (https://www.ncbi.nlm.nih.gov/genbank/genomesubmit_annotation/)
  if(any(duplicated(genes$attribute))){
    # if any, then per unique attribute...
    genes <- do.call(rbind, lapply(split(genes, genes$attribute), function(x){
      # ...find total range of start and end
      tmp_range <- range(x[, c("start", "end")])
      # return only first row of attribute...
      res <- x[1, ]
      # ...with total chromosomal location span
      res$start <- tmp_range[1]
      res$end <- tmp_range[2]
      res
    }))
  }
  # find duplicates of same feature
  # genes_tags <- unlist(lapply(genes$attribute, function(x){
  #   # add ";" at end in case feature_type is last attribute (regular expression requires ending character)
  #   sub(paste0(".*?", feature_type, "=(.*?);.*"), "\\1", paste0(x, ";"))
  # }))
  ## add ";" at end in case feature_type is last attribute (regular expression requires ending character)
  genes$locus_tag <- sub(paste0(".*?", feature_type, "=(.*?);.*"), "\\1", paste0(genes$attribute, ";"))
  # remove duplicates
  genes <- genes[!duplicated(genes_tags), ]
  gene_tags <- genes$locus_tag
} else{
  genome_gb <- readGenBank(input_genome)
  genome <- getSeq(genome_gb)
  genes <- as.data.frame(genes(genome_gb))
  genes$seqid <- genes$seqnames
  genes_tags <- genes$locus_tag
}
genomeID <- switch(input_type, 
                   accessionnr = input_genome, 
                   gbfile = genome_gb@version["accession.version"])
# chromID <- switch(input_type, 
#                   accessionnr = "seqid", 
#                   gbfile = "seqnames")
# create DNAstringset with sequences
genes_seq <- DNAStringSet(unlist(apply(genes, 1, function(x){
  chrom <- grep(x["seqid"], names(genome))
  genome[[chrom]][x["start"]:x["end"]]
  })))
names(genes_seq) <- unique(genes_tags)
# write targets to fasta file if desired
if(output_target_fasta){writeXStringSet(genes_seq, paste0(outdir, "/", genomeID, "_maxmismatch", max_mismatch, "_targets.fasta"))}


#### 4. Find candidate sgRNAs ####
## identify all possible sgRNAs within annotated genetic elements
# multi-core socket in windows, forking otherwise
if(platform == "windows"){
  # socket built-in for function
  candidate_sgRNAs <- findgRNAs(genes_seq, annotatePaired = FALSE, 
                                n.cores.max = n_cores, 
                                enable.multicore = ifelse(n_cores > 1, TRUE, FALSE), 
                                PAM = PAM, 
                                PAM.size = nchar(PAM), 
                                gRNA.size = spacer_length)
} else{
  # looping through index instead of sequence retains feature names
  candidate_sgRNAs <- do.call(c, mclapply(seq.int(genes_seq), function(gene){
    findgRNAs(genes_seq[gene], 
              annotatePaired = FALSE, 
              n.cores.max = 1, 
              enable.multicore = FALSE, 
              PAM = PAM, 
              PAM.size = nchar(PAM), 
              gRNA.size = spacer_length)
  }, mc.cores = n_cores))
}
## retain only unique sgRNAs targeting NT strand (antisense)
# get direction (r/f) of every sgRNA
sgRNAdir <- substr(names(candidate_sgRNAs), nchar(names(candidate_sgRNAs)), nchar(names(candidate_sgRNAs)))
# get target of every sgRNA
targetnames <- sub("\\_gR[0-9]+[rf]$", "", names(candidate_sgRNAs))
# get name of every gene in GFF genes list
#genes_tags_unique <- sub(paste0(".*?", feature_type, "=(.*?);.*"), "\\1", paste0(genes$attribute, ";"))
# get strand of target gene
antisensedir <- ifelse(genes$strand[match(targetnames, genes_tags)] == "+", "r", "f")
#antisensedir <- ifelse(genes$strand[match(targetnames, genes_tags_unique)] == "+", "r", "f")
# keep only uniques of sgRNAs targeting antisense (NT) strand
candidate_sgRNAs_uNT <- unique(candidate_sgRNAs[sgRNAdir == antisensedir])


#### 5. Score sgRNA candidates ####
candidate_hits <- sgRNAefficiencyMC(sgRNAs = candidate_sgRNAs_uNT, 
                                    genes = genes, genome = genome, 
                                    reprAct = TRUE, dist2SC = TRUE, 
                                    name_by = feature_type, 
                                    penalties = reprAct_penalties, 
                                    custom.penalties = reprAct_custom_penalties, 
                                    max.mismatch = max_mismatch, 
                                    PAM = PAM, 
                                    allowed.mismatch.PAM = 1, 
                                    PAM.pattern = paste0(PAM, "$"), 
                                    no_cores = n_cores, 
                                    outfile = paste0(outdir, "/", genomeID))
# if e.g. semicolons in chrom names, will give trouble writing files
candidate_hits$chrom <- gsub(";", "", candidate_hits$chrom)
# write full, scored candidate sgRNA list to file if desired
if(output_full_list){write.csv(candidate_hits, 
                               paste0(outdir, "/", genomeID, "_maxmismatch", max_mismatch, "_candidate_sgRNAs_full.csv"), 
                               row.names = FALSE)}


#### 6. Pick optimal sgRNAs per feature ####
if(output_optimized_list){
  # optimize only for features with >0 perfect on-target sgRNAs
  ON_candidate_hits <- candidate_hits[candidate_hits$n.mismatch == 0, ]
  # if multi core
  if(n_cores > 1){
    if(platform == "windows"){
      # start cluster
      cl <- makeCluster(n_cores)
      clusterExport(cl, c("candidate_hits", "errorRange_maxOffreprAct"))
      # per feature
      optimal_sgRNAs_ls <- parLapply(cl, split(ON_candidate_hits, ON_candidate_hits$NTgene), function(target){
        # get maximum detected off-target effect per candidate sgRNA
        maxOffreprAct <- unlist(lapply(target$name, function(candidate){
          # get all binding sites
          allhits_candidate <- candidate_hits[candidate_hits$name %in% candidate, ]
          # get off-target row IDs
          offID <- !allhits_candidate$NTgene %in% target$NTgene[1]
          # if off-targets, give max reprAct, 0 otherwise (i.e. worst off-target effect)
          ifelse(sum(offID) > 0, 
                 max(allhits_candidate[offID, "reprAct"]), 
                 0)
        }))
        # add max expected off-target effect to data frame
        target$maxOffreprAct <- maxOffreprAct
        # select sgRNA with smallest maxOffreprAct + errorRange_maxOffreprAct
        opt_tmp <- target[target$maxOffreprAct <= (min(maxOffreprAct) + errorRange_maxOffreprAct), ]
        # if multiple, select on smallest dist2SC
        opt_tmp[which.min(opt_tmp$dist2SC), ]
      })
      # shut down workers
      stopCluster(cl)
    } else{
      # per feature
      optimal_sgRNAs_ls <- mclapply(split(ON_candidate_hits, ON_candidate_hits$NTgene), function(target){
        # get maximum detected off-target effect per candidate sgRNA
        maxOffreprAct <- unlist(lapply(target$name, function(candidate){
          # get all binding sites
          allhits_candidate <- candidate_hits[candidate_hits$name %in% candidate, ]
          # get off-target row IDs
          offID <- !allhits_candidate$NTgene %in% target$NTgene[1]
          # if off-targets, give max reprAct, 0 otherwise (i.e. worst off-target effect)
          ifelse(sum(offID) > 0, 
                 max(allhits_candidate[offID, "reprAct"]), 
                 0)
        }))
        # add max expected off-target effect to data frame
        target$maxOffreprAct <- maxOffreprAct
        # select sgRNA with smallest maxOffreprAct + errorRange_maxOffreprAct
        opt_tmp <- target[target$maxOffreprAct <= (min(maxOffreprAct) + errorRange_maxOffreprAct), ]
        # if multiple, select on smallest dist2SC
        opt_tmp[which.min(opt_tmp$dist2SC), ]
      }, mc.cores = n_cores)
    }
  } else{
    # if single-core
    optimal_sgRNAs_ls <- lapply(split(ON_candidate_hits, ON_candidate_hits$NTgene), function(target){
      # per candidate sgRNA
      ## get maximum detected off-target effect (reprAct [0,1])
      maxOffreprAct <- unlist(lapply(target$name, function(candidate){
        # get all binding sites
        allhits_candidate <- candidate_hits[candidate_hits$name %in% candidate, ]
        # get off-target row IDs
        offID <- !allhits_candidate$NTgene %in% target$NTgene[1]
        # if off-targets, give max reprAct, 0 otherwise (i.e. worst off-target effect)
        ifelse(sum(offID) > 0, 
               max(allhits_candidate[offID, "reprAct"]), 
               0)
      }))
      # add max expected off-target effect to data frame
      target$maxOffreprAct <- maxOffreprAct
      # select sgRNA with smallest maxOffreprAct + errorRange_maxOffreprAct
      opt_tmp <- target[target$maxOffreprAct <= (min(maxOffreprAct) + errorRange_maxOffreprAct), ]
      # if multiple, select on smallest dist2SC
      opt_tmp[which.min(opt_tmp$dist2SC), ]
    })
  }
  # collect optimal sgRNAs
  optimal_sgRNAs_df <- do.call(rbind, optimal_sgRNAs_ls)
  # get forward and reverse spacer strands
  optimal_sgRNAs_df$forward <- substr(optimal_sgRNAs_df$targetSequence, 1, spacer_length)
  optimal_sgRNAs_df$reverse <- as.character(reverseComplement(DNAStringSet(optimal_sgRNAs_df$forward)))
  # get forward and reverse oligo's
  optimal_sgRNAs_df$oligoForward <- paste0(oligoForwardOverhang, optimal_sgRNAs_df$forward)
  optimal_sgRNAs_df$oligoReverse <- paste0(oligoReverseOverhang, optimal_sgRNAs_df$reverse)
  # write to file
  write.csv(optimal_sgRNAs_df, file = paste0(outdir, "/", genomeID, "_maxmismatch", max_mismatch, "_sgRNAs_optimal.csv"))
}