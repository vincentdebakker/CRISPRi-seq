#### (Optimal) sgRNA design per annotated feature in any NCBI genome ####
# Author: Vincent de Bakker
# Veening Lab, DMF, FBM, University of Lausanne, Switzerland
# vincent.debakker@unil.ch
####

#### 1. Settings ####
wd <- "~" 
fd <- "~" 
path_ncbi_downloads <- "~"
accession_nr <- "GCA_003003495.1" 
db <- "genbank"
feature_type <- "locus_tag"
n_cores <- 1
max_mismatch <- 6
reprAct_penalties <- "qi.mean.per.region"
reprAct_custom_penalties <- NULL
errorRange_maxOffreprAct <- 0.01
spacer_length <- 20
PAM <- "NGG"
output_target_fasta <- TRUE
output_full_list <- TRUE
output_optimized_list <- TRUE


#### 2. Preliminaries ####
# install required packages if needed
required_packages_CRAN <- c("BiocManager")
required_packages_BioC <- c("Biostrings", "biomaRt", "CRISPRseek")
for(i in seq.int(required_packages_CRAN)){
  if(!requireNamespace(required_packages_CRAN[i], quietly = TRUE)){install.packages(required_packages_CRAN[i])}
}
for(i in seq.int(required_packages_BioC)){
  if(!requireNamespace(required_packages_BioC[i], quietly = TRUE)){BiocManager::install(required_packages_BioC[i])}
}
# load packages
library(Biostrings)
library(biomartr)
library(CRISPRseek)
library(parallel)
# load user-defined functions
source(paste0(fd, "/function_sgRNAefficiencyMC.R"))


#### 3. Get genome, features and extract feature sequences ####
if(file.exists(paste0(path_ncbi_downloads, "/_ncbi_downloads/genomes/", accession_nr,"_genomic_", db, ".fna.gz"))){
  genome_path <- paste0(path_ncbi_downloads, "/_ncbi_downloads/genomes/", accession_nr,"_genomic_", db, ".fna.gz")
} else{
  genome_path <- getGenome(db = db, 
                           organism = accession_nr, 
                           reference = FALSE, 
                           path = paste0(path_ncbi_downloads, "/_ncbi_downloads/genomes/"))
}
genome <- read_genome(genome_path)
if(file.exists(paste0(path_ncbi_downloads, "/_ncbi_downloads/annotation/", accession_nr, "_genomic_", db, ".gff.gz"))){
  gffPath <- paste0(path_ncbi_downloads, "/_ncbi_downloads/annotation/", accession_nr, "_genomic_", db, ".gff.gz")
} else{
  gffPath <- getGFF(db = db, 
                    organism = accession_nr, 
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
genes_tags <- unlist(lapply(genes$attribute, function(x){
  # add ";" at end in case feature_type is last attribute (regular expression requires ending character)
  sub(paste0(".*?", feature_type, "=(.*?);.*"), "\\1", paste0(x, ";"))
}))
# remove duplicates
genes <- genes[!duplicated(genes_tags), ]
# create DNAstringset with sequences
genes_seq <- DNAStringSet(unlist(apply(genes, 1, function(x){
  chrom <- grep(x["seqid"], names(genome))
  genome[[chrom]][x["start"]:x["end"]]
  })))
names(genes_seq) <- unique(genes_tags)
# write targets to fasta file if desired
if(output_target_fasta){writeXStringSet(genes_seq, paste0(wd, "/", accession_nr, "_maxmismatch", max_mismatch, "_targets.fasta"))}


#### 4. Find candidate sgRNAs ####
if(n_cores < 1 | !is.numeric(n_cores)){stop("n_cores should be an integer of 1 or higher.")}
candidate_sgRNAs <- findgRNAs(genes_seq, annotatePaired = FALSE, 
                              n.cores.max = n_cores, 
                              enable.multicore = ifelse(n_cores > 1, TRUE, FALSE), 
                              PAM = PAM, 
                              PAM.size = nchar(PAM), 
                              gRNA.size = spacer_length)
# retain only unique sgRNAs targeting NT strand (antisense)
# get direction (r/f) of every sgRNA
sgRNAdir <- substr(names(candidate_sgRNAs), nchar(names(candidate_sgRNAs)), nchar(names(candidate_sgRNAs)))
# get target of every sgRNA
targetnames <- sub("\\_gR[0-9]+[rf]$", "", names(candidate_sgRNAs))
# get name of every gene in GFF genes list
genes_tags_unique <- sub(paste0(".*?", feature_type, "=(.*?);.*"), "\\1", paste0(genes$attribute, ";"))
# get strand of target gene
antisensedir <- ifelse(genes$strand[match(targetnames, genes_tags_unique)] == "+", "r", "f")
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
                                    outfile = paste0(wd, "/", accession_nr))
# remove non-informative columns 
candidate_hits <- candidate_hits[, -match(c("forViewInUCSC", "score"), colnames(candidate_hits))]
# write full, scored candidate sgRNA list to file if desired
if(output_full_list){write.csv(candidate_hits, 
                               paste0(wd, "/", accession_nr, "_maxmismatch", max_mismatch, "_candidate_sgRNAs_full.csv"), 
                               row.names = FALSE)}


#### 6. Pick optimal sgRNAs per feature ####
if(output_optimized_list){
  # optimize only for features with >0 perfect on-target sgRNAs
  ON_candidate_hits <- candidate_hits[candidate_hits$n.mismatch == 0, ]
  # if multi core
  if(n_cores > 1){
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
  colnames(optimal_sgRNAs_df)[27] <- "TargetSequence"
  # write to file
  write.csv(optimal_sgRNAs_df, file = paste0(wd, "/", accession_nr, "_maxmismatch", max_mismatch, "_sgRNAs_optimal.csv"))
}