#### Target-site detection and reporting for given sgRNAs in any NCBI genome ####
# Author: Vincent de Bakker
# Veening Lab, DMF, FBM, University of Lausanne, Switzerland
# vincent.debakker@unil.ch
####

#### 1. Settings ####
wd <- "~" 
fd <- wd 
path_ncbi_downloads <- "~" 
sgRNA_file <- paste0(wd, "/sgRNA-library_spn-D39V.csv")
accession_nr <- "GCA_003003495.1"
out_path <- wd
max_mismatch <- 8
n_cores <- 1
feature_type <- "locus_tag"
reprAct_penalties <- "qi.mean.per.region"
reprAct_custom_penalties <- NULL
spacer_length <- 20
PAM <- "NGG"
db <- "genbank"


#### 2. Preliminaries ####
required_packages_CRAN <- c("BiocManager")
required_packages_BioC <- c("Biostrings", "CRISPRseek", "biomartr")
for(i in seq.int(required_packages_CRAN)){
  if(!requireNamespace(required_packages_CRAN[i], quietly = TRUE)){install.packages(required_packages_CRAN[i])}
}
for(i in seq.int(required_packages_BioC)){
  if(!requireNamespace(required_packages_BioC[i], quietly = TRUE)){BiocManager::install(required_packages_BioC[i])}
}
lapply(required_packages_CRAN[-1], library, character.only = TRUE)
lapply(required_packages_BioC, library, character.only = TRUE)

# load sgRNA efficiency function
source(paste0(fd, "/function_sgRNAefficiencyMC.R"))


#### 3. Read sgRNAs ####
sgRNAs <- read.csv(sgRNA_file, header = FALSE)
colnames(sgRNAs) <- c("tag", "sgRNA")
## keep only sgRNAs of length spacer_length, also removing comments as e.g. "No PAM"
sgRNAs <- sgRNAs[unlist(lapply(sgRNAs$sgRNA, nchar)) == spacer_length, ]
## check if no duplicates in either tag or sgRNA
if(any(duplicated(sgRNAs$sgRNA))){
  warning(paste("Duplicate sgRNA sequences:", paste(sgRNAs$tag[duplicated(sgRNAs$sgRNA) | duplicated(sgRNAs$sgRNA, fromLast = TRUE)], collapse = ", ")))
}
if(any(duplicated(sgRNAs$tag))){stop("Duplicate sgRNA names, please ensure all names are unique")}
## attach PAM and make sure all bases are upper case
sgRNAs$sgRNA <- apply(sgRNAs, 1, function(x){paste0(toupper(substr(x[2], 1, spacer_length)), PAM)})
# format for CRISPRseek::searchHits
gRNAs <- DNAStringSet(unlist(sgRNAs[, 2]), use.names = FALSE)
names(gRNAs) <- unlist(sgRNAs[, 1])


#### 4. Get genome, features and extract feature sequences ####
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
# extract all features with feature_type tag
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


#### 5. Find and score binding sites per sgRNA ####
effic <- sgRNAefficiencyMC(sgRNAs = gRNAs, 
                           genes = genes, genome = genome, 
                           reprAct = TRUE, dist2SC = TRUE, 
                           name_by = feature_type, 
                           penalties = reprAct_penalties, 
                           custom.penalties = reprAct_custom_penalties, 
                           outfile = paste0(out_path, "/", accession_nr), 
                           max.mismatch = max_mismatch, 
                           PAM = PAM, allowed.mismatch.PAM = 1, 
                           PAM.pattern = paste0(PAM, "$"), 
                           no_cores = n_cores)
# save results
write.csv(effic, 
          file = paste0(out_path, "/sgRNA-library_binding-sites_", accession_nr, ".csv"), 
          row.names = FALSE)