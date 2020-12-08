#### Target-site detection and reporting for given sgRNAs in any NCBI genome ####
# Author: Vincent de Bakker
# Veening Lab, DMF, FBM, University of Lausanne, Switzerland
# vincent.debakker@unil.ch
####

#### 1. Settings ####
wd <- "~" # working directory
fd <- "~" # function directory
path_ncbi_downloads <- "~" 
sgRNA_file <- "~/sgRNA_seq_all-lib-1499.xlsx"
accession_nr <- "GCA_003003495.1" # S. pneumoniae D39V: GCA_003003495.1
out_path <- wd
n_cores <- 1
PAM <- "NGG"
db <- "genbank"   # options: genbank, refseq, ensembl
max_mismatch <- 8 # I use 8; more takes much longer and is not worth it


#### 2. Preliminaries ####
required_packages_CRAN <- c("BiocManager", "readxl", "writexl")
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
sgRNAs_raw <- read_excel(sgRNA_file)
sgRNAs <- sgRNAs_raw[, c("Locustag", "sgRNA sequence")]
## would be better have sgRNA tag, not locus (target) tag, but will work
colnames(sgRNAs) <- c("tag", "sgRNA")
## keep only sgRNAs of length 20, also removing comments as e.g. "No PAM"
sgRNAs <- sgRNAs[unlist(lapply(sgRNAs$sgRNA, nchar)) == 20, ]
## filter out any leftover non-sequence (comment), check manually:
#sgRNAs[!grepl("G", sgRNAs$sgRNA) | !grepl("C", sgRNAs$sgRNA) | !grepl("A", sgRNAs$sgRNA) | !grepl("T", sgRNAs$sgRNA), ]
#sgRNAs[grep("tar", sgRNAs$sgRNA), ]
sgRNAs <- sgRNAs[-(grep("tar", sgRNAs$sgRNA)), ]
## check if no duplicates in either tag or sgRNA
if(any(duplicated(sgRNAs$sgRNA))){stop("Duplicate sgRNA sequences")}
if(any(any(duplicated(sgRNAs$tag)))){stop("Duplicate sgRNA tags")}
## attach PAM and make sure all bases are upper case
sgRNAs$sgRNA <- apply(sgRNAs, 1, function(x){paste0(toupper(substr(x[2], 1, 20)), PAM)})

# format for CRISPRseek::searchHits
gRNAs <- DNAStringSet(unlist(sgRNAs[, 2]), use.names = FALSE)
names(gRNAs) <- unlist(sgRNAs[, 1])
## check manually
#gRNAs


#### 4. Get genome, features and extract feature sequences ####
if(file.exists(paste0(path_ncbi_downloads, "_ncbi_downloads/genomes/", accession_nr,"_genomic_", db, ".fna.gz"))){
  genome_path <- paste0(path_ncbi_downloads, "_ncbi_downloads/genomes/", accession_nr,"_genomic_", db, ".fna.gz")
} else{
  genome_path <- getGenome(db = db, 
                           organism = accession_nr, 
                           reference = FALSE, 
                           path = paste0(path_ncbi_downloads, "_ncbi_downloads/genomes/"))
}
genome <- read_genome(genome_path)
if(file.exists(paste0(path_ncbi_downloads, "_ncbi_downloads/annotation/", accession_nr, "_genomic_", db, ".gff.gz"))){
  gffPath <- paste0(path_ncbi_downloads, "_ncbi_downloads/annotation/", accession_nr, "_genomic_", db, ".gff.gz")
} else{
  gffPath <- getGFF(db = db, 
                    organism = accession_nr, 
                    reference = FALSE, 
                    path = paste0(path_ncbi_downloads, "_ncbi_downloads/annotation/"))
}
GFF <- read_gff(gffPath)
# extract all features with a locus tag
genes <- GFF[unlist(lapply(GFF$attribute, grepl, pattern = "locus_tag")), ]
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
  # add ";" at end in case locus_tag is last attribute (regular expression requires ending character)
  sub(paste0(".*?", "locus_tag", "=(.*?);.*"), "\\1", paste0(x, ";"))
}))
# remove duplicates
genes <- genes[!duplicated(genes_tags), ]


#### 5. Find and score binding sites per sgRNA ####
effic <- sgRNAefficiencyMC(sgRNAs = gRNAs, 
                           genes = genes, genome = genome, 
                           reprAct = TRUE, dist2SC = TRUE, 
                           name_by = "locus_tag", 
                           penalties = "qi.mean.per.region", 
                           outfile = paste0(out_path, "/", accession_nr), 
                           max.mismatch = max_mismatch, 
                           PAM = PAM, allowed.mismatch.PAM = 1, 
                           PAM.pattern = paste0(PAM, "$"), 
                           no_cores = n_cores)
# save results
#write.csv(effic, paste0(out_path, "lib_", accession_nr, ".csv"))
write_xlsx(effic, paste0(out_path, "/lib_", accession_nr, ".xlsx"))
#write.table(effic, file = paste0(out_path, "lib_", accession_nr, ".tsv"), 
#            sep = "\t", row.names = FALSE, quote = FALSE)