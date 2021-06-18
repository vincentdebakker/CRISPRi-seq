#### (Optimal) sgRNA design per annotated feature ####
# Author: Vincent de Bakker
# Veening Lab, DMF, FBM, University of Lausanne, Switzerland
# vincent.debakker@unil.ch
####


#### 1. Settings ####
## REQUIRED ##
input_genome <- "INPUT GENOME HERE"
outdir <- "OUTPUT DIRECTORY HERE"
TINDRidir <- "TINDRi.py DIRECTORY HERE"

## OPTIONAL ##
n_sgRNA <- 1
path_ncbi_downloads <- NA
regions <- c(7, 2, 11) 
max_mismatch_cum <- c(1, 2, 11)
reprAct_penalties <- "HawkinsBsuMedian"
pen_func <- prod
errorRange_maxOffreprAct <- 0.4
avoidNonBaseNT <- TRUE
bad_seeds <- c("ACCCA", "TGGAA")
bad_seed_rule <- "ignore"
cut_sites <- c(BsmBI = "CGTCTC", tandemT = "TTTTT")
cut_site_rule <- "avoid"
oligoForwardOverhang <- "TATA"
oligoReverseOverhang <- "AAAC"
PAM <- "NGG"
filter_out_duplicates <- TRUE
output_optimized_list <- TRUE
output_all_candidates <- FALSE
output_target_fasta <- FALSE
output_sgRNAs_fasta <- FALSE
output_sites_fasta <- FALSE
output_full_list <- FALSE 
detect_offtarget_genes_full <- FALSE
keep_TINDRi_input_sgRNAs <- FALSE
keep_TINDRi_input_sites <- FALSE
keep_TINDRi_matches <- TRUE
path_python <- NULL
feature_type <- "locus_tag"


#### 2. Preliminaries ####
# starting time
message(paste0(Sys.time(), ": design pipeline started"))
start_time <- Sys.time()
# check inputs
if(any(endsWith(input_genome, c(".gb", ".gbf", ".gbff", ".gbk")))){
  input_type <- "gbfile"
} else{
  if(any(startsWith(input_genome, c("GCA", "GCF")))){
    input_type <- "accessionnr"
  } else{
    stop("Input should either be one of .gb, .gbf, .gbff or .gbk files or an NCBI genome assembly accession number (GCA_ or GCF_).")
  }
}
if(input_type == "accessionnr" & is.na(path_ncbi_downloads)){
  stop("Please provide the path to the directory to save downloaded genomes and annotations to input parameter path_ncbi_downloads.")
}
if(input_type == "gbfile" & !file.exists(input_genome)){
  stop("Input parameter input_genome: no such file. Please specify correct and complete path to a .gb, .gbf .gbff or .gbk file, or give an NCBI assembly accession number (e.g. GCA_003003495.1).")
}
if(!file.exists(paste0(TINDRidir, "/TINDRi.py"))){
  stop("Please save the Python script TINDRi.py in the directory specified for input parameter 'TINDRidir'.")
}
if(input_type == "accessionnr"){
  db <- switch(substr(input_genome, 1, 3), 
               GCA = "genbank", 
               GCF = "refseq")
  message(paste("Using", db, "data base. (Please use GCA accession for genbank and GCF for refseq.)"))
}
if(length(regions) != length(max_mismatch_cum)){
  stop("Input parameters regions and max_mismatch_cum need to be of equal length.")
}
if(!output_full_list & !output_optimized_list){
  warning("Both output_optimized_list and output_full_list are set to FALSE, which is allowed, but not very useful. Are you sure you want to continue? If not, abord the script and adjust the input parameters.")
}
if(detect_offtarget_genes_full){
  warning("Parameter detect_offtarget_genes_full is set to TRUE. This will take considerably more time, and is not required for optimal sgRNA selection. Please ensure this is desired and otherwise abord and change the input parameter.")
}
# set penalties for reprAct
pen_ls <- list("HawkinsEcoMean" = c(0.97, 0.93, 0.89, 0.91, 0.75, 0.82, 0.83, 0.87, 0.78, 0.74, 0.72, 0.66, 0.59, 0.45, 0.38, 0.47, 0.46, 0.46, 0.32, 0.19), 
               "HawkinsBsuMean" = c(0.99, 0.91, 0.84, 0.85, 0.77, 0.88, 0.85, 0.88, 0.85, 0.81, 0.79, 0.76, 0.69, 0.6, 0.56, 0.58, 0.58, 0.53, 0.4, 0.3), 
               "HawkinsEcoMedian" = c(0.98, 0.97, 0.93, 0.96, 0.78, 0.94, 0.87, 0.88, 0.88, 0.85, 0.79, 0.72, 0.69, 0.41, 0.36, 0.46, 0.39, 0.35, 0.29, 0.18), 
               "HawkinsBsuMedian" = c(0.99, 0.96, 0.94, 0.94, 0.84, 0.95, 0.92, 0.95, 0.94, 0.91, 0.9, 0.83, 0.67, 0.57, 0.53, 0.57, 0.53, 0.5, 0.46, 0.35), 
               "Qi" = c(0.78, 0.41, 0.44, 0.63, 0.65, 0.68, 0.65, 0.63,
                        0.30, 0.25, 0.24, 0.22, 0.24,
                        0.01, 0.12, 0.10, 0.07, 0.06, 0.09, 0.05))
pen_ls$QiMean <- rep(tapply(pen_ls$Qi, rep(c("1III", "2II", "3I"), c(8, 5, 7)), mean), c(8, 5, 7))
pen_ls$QiMedian <- rep(tapply(pen_ls$Qi, rep(c("1III", "2II", "3I"), c(8, 5, 7)), median), c(8, 5, 7))
if(!(reprAct_penalties %in% names(pen_ls)) & (!(is.numeric(reprAct_penalties)) | length(reprAct_penalties) != sum(regions))){
  stop(paste("For custom penalties, please set reprAct_penalties to a numeric vector of the length of the sgRNA spacer (i.e. sum(regions)). Otherwise set it to one of these standard options:", 
             paste(names(pen_ls), collapse = ", ")))
}
if(reprAct_penalties %in% names(pen_ls)){
  penalties <- pen_ls[[match(reprAct_penalties, names(pen_ls))]]
} else{
  penalties <- rev(reprAct_penalties) # custom input given from PAM-proximal to -distal
}
# install required packages if needed
message(paste0(Sys.time(), ": Loading packages (and installing if needed)..."))
required_packages_CRAN <- c("BiocManager", "reticulate")
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
invisible(suppressWarnings(suppressPackageStartupMessages(lapply(c(required_packages_BioC, 
                                                                   required_packages_CRAN[-1]), 
                                                                 library, character.only = TRUE)
)))
library(parallel)
#library(reticulate)
# load user-defined functions
#source(paste0(TINDRidir, "/function_sgRNAefficiencyMC.R"))
# check python installations
# if(is.null(path_python) & !py_available()){
#   stop("No installation of Python detected. Please ensure you have an installed instance of Python3 and if this error persists, supply the full path to the python.exe file to the path_python input parameter.")
# }
if(!is.null(path_python)){
  if(file.exists(path_python)){
    use_python(path_python)
  } else{
    stop("Path to Python3 installation invalid: python.exe was not found. Please ensure full and complete path specification if path_python is provided.")
  }
}
# check if required Python modules available
py_required_modules <- c("numpy", "numba", "time", "multiprocessing", "datetime", "os", "sys")
py_required_modules_available <- unlist(lapply(py_required_modules, py_module_available))
if(!py_available()){
  stop("Path to Python3 installation invalid: python was not found. Please ensure full and complete path specification if path_python is provided. Note that if the parameter is set to NULL, the script will look for installed python instances itself.")
}
if(any(!py_required_modules_available)){
  message("Trying to install required Python modules...")
  try(py_install(py_required_modules[!py_required_modules_available]), silent = TRUE)
}
if(any(!py_required_modules_available)){
  stop(paste("Not all Python modules are available or could be installed automatically. Please install manually:", 
             paste(py_required_modules[!py_required_modules_available], collapse = ",")))
}
# set n_cores
n_cores <- ifelse(detectCores() > 1, detectCores() - 1, 1)
# detect OS
platform <- .Platform$OS.type


#### 3. Get genome, features and extract feature sequences ####
message(paste0(Sys.time(), ": Loading genome, features and sequences..."))
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
  GFF <- suppressMessages(read_gff(gffPath))
  # extract all features of type feature_type
  genes <- GFF[unlist(lapply(GFF$attribute, grepl, pattern = feature_type)), ]
  rm(GFF)
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
  genes_tags <- genes$locus_tag
  # remove duplicates
  genes <- genes[!duplicated(genes_tags), ]
  # add gene names where possible for final output
  genes$gene <- ifelse(grepl("gene=", genes$attribute), 
                       sub(paste0(".*?", "gene", "=(.*?);.*"), "\\1", paste0(genes$attribute, ";")), 
                       genes$locus_tag)
} else{
  # made generic for multi-chromosome cases (plasmids)
  genome_txt <- readLines(input_genome)
  # chromosomes are split in .gbff files by //
  #chromsep <- grep("//", genome_txt)
  # grep also grabs non-exact matches like URLs
  chromsep <- which(genome_txt == "//")
  if(length(chromsep) > 100){message("Reading in many (", length(chromsep), ") chromosomes...")}
  # read in genbank object for each chromosome separately
  genome_gb_ls <- lapply(seq.int(length(chromsep)), function(chrom){
    chrom_txt <- if(chrom == 1){
      genome_txt[seq.int(chromsep[chrom])]
    } else{
      genome_txt[seq.int(chromsep[chrom - 1] + 1, 
                         chromsep[chrom], 
                         by = 1)]
    }
    # prokka sometimes doesn't put space between ID and length in LOCUS
    chrom_txt[1] <- sub("([0123456789]+\\sbp)", " \\1", chrom_txt[1])
    # readGenBank CANNOT HANDLE a 't' in the LOCUS line
    chrom_txt[1] <- gsub("t", "", chrom_txt[1])
    chrom_gb <- suppressWarnings(readGenBank(text = chrom_txt))
    return(chrom_gb)
  })
  # merge all chromosome information
  genome <- do.call(c, lapply(seq.int(length(genome_gb_ls)), function(chrom){
    chrom_seq <- getSeq(genome_gb_ls[[chrom]])
    chrom_loc <- strsplit(locus(genome_gb_ls[[chrom]]), " ")[[1]]
    # add chrom nr to ensure unique IDs
    names(chrom_seq) <- paste0(chrom_loc[chrom_loc != ""][2], "_", chrom)
    # for prokka files, no difference between chromosomes in any ID...
    #names(chrom_seq) <- definition(genome_gb_ls[[chrom]])
    return(chrom_seq)
  }))
  #genome <- do.call(c, lapply(genome_gb_ls, getSeq))
  # make sure resulting gene tables can be merged by selecting same columns
  genes_cols <- c("seqnames", "start", "end", "strand", "locus_tag")
  #genes_cols <- c("seqnames", "start", "end", "strand", "gene", "locus_tag", "gene_id")
  genes <- do.call(rbind, lapply(seq.int(length(genome_gb_ls)), function(chrom){
    # find features on each contig
    if(isEmpty(genes(genome_gb_ls[[chrom]]))){
      if(isEmpty(transcripts(genome_gb_ls[[chrom]])) & 
         (isEmpty(otherFeatures(genome_gb_ls[[chrom]])) | 
          all(otherFeatures(genome_gb_ls[[chrom]])$type == "assembly_gap"))){
        message("No features detected to design sgRNAs for on contig ", names(genome)[chrom])
        NULL
      } else{
        if(isEmpty(otherFeatures(genome_gb_ls[[chrom]]))){
          genes_tmp <- as.data.frame(transcripts(genome_gb_ls[[chrom]]))
          data.frame(genes_tmp[, genes_cols], 
                     "gene" = if("gene" %in% colnames(genes_tmp)){genes_tmp$gene}else{NA},
                     "seqid" = names(genome)[chrom])
        } else{
          # exclude other features without locus_tag
          othfeat <- as.data.frame(otherFeatures(genome_gb_ls[[chrom]]))
          # use product for column "gene" if is.na(gene); replace spaces by _
          othfeat$gene <- gsub(" ", "_", 
                               ifelse(is.na(othfeat$gene), othfeat$product, othfeat$gene))
          if(isEmpty(transcripts(genome_gb_ls[[chrom]]))){
            if("locus_tag" %in% colnames(othfeat)){
              data.frame(othfeat[!is.na(othfeat$locus_tag), genes_cols], 
                         "seqid" = names(genome)[chrom])
            } else{
              NULL
            }
            # data.frame(as.data.frame(otherFeatures(genome_gb_ls[[chrom]]))[, genes_cols],
            #            "seqid" = names(genome)[chrom])
          } else{
            genes_tmp <- as.data.frame(transcripts(genome_gb_ls[[chrom]]))
            if("locus_tag" %in% colnames(othfeat)){
              data.frame(rbind(othfeat[!is.na(othfeat$locus_tag), genes_cols], 
                               data.frame(genes_tmp[, genes_cols]), 
                               "gene" = if("gene" %in% colnames(genes_tmp)){genes_tmp$gene}else{NA}), 
                         "seqid" = names(genome)[chrom])
            } else{
              # superfluous?
              genes_tmp <- as.data.frame(transcripts(genome_gb_ls[[chrom]]))
              data.frame(genes_tmp[, genes_cols], 
                         "gene" = if("gene" %in% colnames(genes_tmp)){genes_tmp$gene}else{NA},
                         "seqid" = names(genome)[chrom])
            }
            # data.frame(rbind(as.data.frame(transcripts(genome_gb_ls[[chrom]]))[, genes_cols], 
            #                  as.data.frame(otherFeatures(genome_gb_ls[[chrom]]))[, genes_cols]), 
            #            "seqid" = names(genome)[chrom])
          }
        }
      }
    } else{
      # prokka does not always annotate gene flags, but if present, easier:
      genes_tmp <- as.data.frame(genes(genome_gb_ls[[chrom]]))
      data.frame(genes_tmp[, genes_cols], 
                 "gene" = if("gene" %in% colnames(genes_tmp)){genes_tmp$gene}else{NA},
                 "seqid" = names(genome)[chrom])
    }
  }))
  # genome_gb <- readGenBank(input_genome)
  # genome <- getSeq(genome_gb)
  # genes <- as.data.frame(genes(genome_gb))
  # genes$seqid <- genes$seqnames
  # genes_tags <- genes$locus_tag
}
genes <- genes[, c("seqid", "strand", "start", "end", "locus_tag", "gene")]
genomeID <- switch(input_type, 
                   accessionnr = input_genome, 
                   #gbfile = definition(genome_gb_ls[[1]]), 
                   gbfile = sub("\\..*", "", sub(".*/", "", input_genome)))
# chromID <- switch(input_type, 
#                   accessionnr = "seqid", 
#                   gbfile = "seqnames")
# create DNAstringset with sequences
genes_seq <- DNAStringSet(unlist(apply(genes, 1, function(x){
  chrom <- switch(input_type, 
                  accessionnr = grep(x["seqid"], names(genome)), 
                  gbfile = which(names(genome) == x["seqid"]))
  #chrom <- which(names(genome) == x["seqid"])
  #chrom <- grep(x["seqid"], names(genome))
  genome[[chrom]][x["start"]:x["end"]]
})))
names(genes_seq) <- genes$locus_tag
#names(genes_seq) <- unique(genes_tags)
# write targets to fasta file if desired
if(output_target_fasta){writeXStringSet(genes_seq, paste0(outdir, "/", genomeID, "_targets.fasta"))}


#### 4. Find all candidate sgRNAs ####
message(paste0(Sys.time(), ": Identifying all candidate sgRNAs..."))
# multi-core socket in windows, forking otherwise
if(platform == "windows"){
  # socket built-in for function
  invisible(capture.output(candidate_sgRNAs <- findgRNAs(genes_seq, annotatePaired = FALSE, 
                                                         n.cores.max = n_cores, 
                                                         enable.multicore = ifelse(n_cores > 1, TRUE, FALSE), 
                                                         PAM = PAM, 
                                                         PAM.size = nchar(PAM), 
                                                         gRNA.size = sum(regions))))
} else{
  # looping through index instead of sequence retains feature names
  invisible(capture.output(candidate_sgRNAs <- do.call(c, mclapply(seq.int(genes_seq), function(gene){
    if(width(genes_seq[gene]) >= sum(regions)){
      findgRNAs(genes_seq[gene], 
                annotatePaired = FALSE, 
                n.cores.max = 1, 
                enable.multicore = FALSE, 
                PAM = PAM, 
                PAM.size = nchar(PAM), 
                gRNA.size = sum(regions))
    } else{
      message("No sgRNA designed for feature ", names(genes_seq[gene]), ": feature length <", sum(regions))
    }
  }, mc.cores = n_cores))))
}
# remove, not needed anymore
rm(genes_seq)
# retain only unique sgRNAs targeting NT strand (antisense)
# get direction (r/f) of every sgRNA
sgRNAdir <- substr(names(candidate_sgRNAs), nchar(names(candidate_sgRNAs)), nchar(names(candidate_sgRNAs)))
# get target of every sgRNA
targetnames <- sub("\\_gR[0-9]+[rf]$", "", names(candidate_sgRNAs))
# get name of every gene in GFF genes list
#genes_tags_unique <- sub(paste0(".*?", feature_type, "=(.*?);.*"), "\\1", paste0(genes$attribute, ";"))
# get strand of target gene
antisensedir <- ifelse(genes$strand[match(targetnames, genes$locus_tag)] == "+", "r", "f")
# keep only sgRNAs targeting antisense (NT) strand
candidate_sgRNAs_uNT <- candidate_sgRNAs[sgRNAdir == antisensedir]
# write to file for which genes no sgRNA could be designed
targetnames <- sub("\\_gR[0-9]+[rf]$", "", names(candidate_sgRNAs_uNT))
#no_sgRNA_targets <- genes$locus_tag[!genes$locus_tag %in% targetnames]
#paste(no_sgRNA_targets, genes$gene[match(no_sgRNA_targets, genes$locus_tag)], sep = ",")
#writeLines(no_sgRNA_targets, paste0(outdir, "/", genomeID, "_no_sgRNA_targets.txt"))
write.csv(genes[!genes$locus_tag %in% targetnames, ], 
          paste0(outdir, "/", genomeID, "_no_sgRNA_targets.csv"), 
          row.names = FALSE)
# keep only unique sgRNAs; NO bc later ID'd by name, then disappears
#candidate_sgRNAs_uNT <- unique(candidate_sgRNAs[sgRNAdir == antisensedir])
#candidate_sgRNAs_uNT <- unique(candidate_sgRNAs_uNT)
# remove what's not needed anymore
rm(candidate_sgRNAs, antisensedir, targetnames, sgRNAdir)
# exclude bad seed sgRNAs if desired
if(bad_seed_rule == "exclude"){
  # bad seeds can have different lengths so lapply over them
  # location is: NNN-badseed-PAM ; e.g. N(15)-TGGAA-NGG (Cui et al)
  bad_seed_ls <- lapply(bad_seeds, function(seed){
    seed_seq <- narrow(candidate_sgRNAs_uNT, 
                       sum(regions) - nchar(seed) + 1, 
                       sum(regions))
    which(seed_seq == seed)
  })
  exclID <- unique(unlist(bad_seed_ls))
  candidate_sgRNAs_uNT <- candidate_sgRNAs_uNT[-exclID]
}
if(cut_site_rule == "exclude"){
  # only look in spacers
  spacers <- narrow(candidate_sgRNAs_uNT, 1, sum(regions))
  # but plus overhangs; both strands; same as rev.compl. on one strand
  spacers_ext <- xscat(oligoForwardOverhang, spacers, reverseComplement(DNAString(oligoReverseOverhang)))
  # also check reverse complements
  cut_sites_plusrev <- c(reverseComplement(DNAStringSet(cut_sites)), cut_sites)
  # don't use pdict; can only have one width for all sites
  #cut_site_pdict <- PDict(cut_sites_plusrev)
  cut_site_found <- vwhichPDict(cut_sites_plusrev, spacers_ext)
  exclID <- which(unlist(lapply(cut_site_found, length)) != 0)
  candidate_sgRNAs_uNT <- candidate_sgRNAs_uNT[-exclID]
}
if(output_sgRNAs_fasta){writeXStringSet(candidate_sgRNAs_uNT, paste0(outdir, "/", genomeID, "_all_sgRNA_candidates.fasta"))}
# write spacer in reverse order (PAM-proximal to -distal) for TINDRi
write.table(cbind(1:length(candidate_sgRNAs_uNT), 
                  as.character(reverse(narrow(candidate_sgRNAs_uNT, 1, sum(regions))))), 
            paste0(outdir, "/sg_candidates_within_genes.csv"), 
            sep = ",", 
            col.names = FALSE, 
            row.names = FALSE, quote = FALSE)


#### 5. Find all sgRNA binding sites ####
message(paste0(Sys.time(), ": Identifying all sgRNA binding sites..."))
# multi-core socket in windows, forking otherwise
if(platform == "windows"){
  # socket built-in for function
  all_sites <- findgRNAs(genome, annotatePaired = FALSE, 
                         n.cores.max = n_cores, 
                         enable.multicore = ifelse(n_cores > 1, TRUE, FALSE), 
                         PAM = PAM, 
                         PAM.size = nchar(PAM), 
                         gRNA.size = sum(regions))
} else{
  all_sites <- do.call(c, mclapply(seq.int(genome), function(chrom){
    findgRNAs(genome[chrom], 
              annotatePaired = FALSE, 
              n.cores.max = 1, 
              enable.multicore = FALSE, 
              PAM = PAM, 
              PAM.size = nchar(PAM), 
              gRNA.size = sum(regions))
  }, mc.cores = n_cores))
}
rm(genome)
if(output_sites_fasta){writeXStringSet(all_sites, paste0(outdir, "/", genomeID, "_all_binding_sites.fasta"))}
# write spacer in reverse order (PAM-proximal to -distal) for TINDRi
write.table(cbind(1:length(all_sites), 
                  as.character(reverse(narrow(all_sites, 1, sum(regions))))), 
            paste0(outdir, "/allsgcandidates.csv"), 
            sep = ",", 
            col.names = FALSE, 
            row.names = FALSE, quote = FALSE)


#### 6. Match candidate sgRNAs to binding sites ####
message(paste0(Sys.time(), ": Matching sgRNA candidates to binding sites..."))
if(is.null(path_python)){path_python <- "python"}
if(platform == "windows"){
  par_TINDRi <- gsub(" ", ",", paste(cumsum(regions), max_mismatch_cum, collapse = ";"))
} else{
  par_TINDRi <- gsub(" ", ",", paste(cumsum(regions), max_mismatch_cum, collapse = "\\;"))
}
system(paste(path_python, paste0(TINDRidir, "/TINDRi.py"), outdir, par_TINDRi), 
       wait = TRUE)
# remove temp files not needed anymore, unless requested to keep
if(!keep_TINDRi_input_sgRNAs){
  invisible(file.remove(paste0(outdir, "/sg_candidates_within_genes.csv")))
} else{
  invisible(file.rename(paste0(outdir, "/sg_candidates_within_genes.csv"), 
                        paste0(outdir, "/", genomeID, "_sg_candidates_within_genes.csv")))
}
if(!keep_TINDRi_input_sites){
  invisible(file.remove(paste0(outdir, "/allsgcandidates.csv")))
} else{
  invisible(file.rename(paste0(outdir, "/allsgcandidates.csv"), 
                        paste0(outdir, "/", genomeID, "_allsgcandidates.csv")))
}


#### 7. Compute required scores for sgRNA-site combinations ####
message(paste0(Sys.time(), ": Scoring sgRNA - binding site combinations..."))
# compute scores for mismatch vectors
candidate_mm <- unlist(read.table(paste0(outdir, "/missmatch_matrix_end.txt"), 
                                  colClasses = c("NULL", "character", "NULL"), 
                                  sep = "\t", 
                                  header = FALSE), 
                       use.names = FALSE)
# compute scores separately to save RAM, but takes longer
if(platform == "windows"){
  # needed anyway
  cl <- makeCluster(n_cores)
  clusterExport(cl, c("penalties", "pen_func"))
  reprAct <- unlist(parLapply(cl, candidate_mm, function(mm){
    #prod(penalties[as.integer(unlist(strsplit(mm, ","))) == 1])
    # leave penalty function customizable
    # invert penalties (to PAM-proximal to -distal)
    do.call(pen_func, list(rev(penalties)[as.integer(unlist(strsplit(mm, ","))) == 1]))
  }))
  stopCluster(cl)
  if(output_full_list){
    write(reprAct, file = paste0(outdir, "/tmp_reprAct.tsv"), sep = "\t", ncolumns = 1)
  }
  # not required otherwise, can deduce perfect matches from reprAct
  if(any(penalties == 1) | !identical(pen_func, prod) | output_full_list){
    cl <- makeCluster(n_cores)
    nmm <- unlist(parLapply(cl, candidate_mm, function(mm){
      sum(as.integer(unlist(strsplit(mm, ","))))
    }))
    stopCluster(cl)
    candidate_ID <- which(nmm == 0)
    if(output_full_list){
      write(nmm, file = paste0(outdir, "/tmp_nmm.tsv"), sep = "\t", ncolumns = 1)
    }
    rm(nmm)
  } else{
    candidate_ID <- which(reprAct == 1)
  }
} else{
  reprAct <- unlist(mclapply(candidate_mm, function(mm){
    #prod(penalties[as.integer(unlist(strsplit(mm, ","))) == 1])
    # leave penalty function customizable
    do.call(pen_func, list(rev(penalties)[as.integer(unlist(strsplit(mm, ","))) == 1]))
  }, mc.cores = n_cores))
  if(output_full_list){
    write(reprAct, file = paste0(outdir, "/tmp_reprAct.tsv"), sep = "\t", ncolumns = 1)
  }
  # not required otherwise, can deduce perfect matches from reprAct
  if(any(penalties == 1) | output_full_list){
    nmm <- unlist(mclapply(candidate_mm, function(mm){
      sum(as.integer(unlist(strsplit(mm, ","))))
    }, mc.cores = n_cores))
    candidate_ID <- which(reprAct == 1)
    if(output_full_list){
      write(nmm, file = paste0(outdir, "/tmp_nmm.tsv"), sep = "\t", ncolumns = 1)
    }
    rm(nmm)
  } else{
    candidate_ID <- which(reprAct == 1)
  }
}
rm(candidate_mm)
# required for either output
# CRISPRseek offset in within-target bp locus of sgRNA names
offset_emp <- 3


#### 8. Return all sgRNA candidates & scores ####
if(output_all_candidates){
  message(paste0(Sys.time(), ": Computing scores for all candidate sgRNAs..."))
  # for maxOffreprAct need all
  candidate_sgID <- unlist(read.table(paste0(outdir, "/missmatch_matrix_end.txt"), 
                                      colClasses = c("integer", "NULL", "NULL"), 
                                      sep = "\t", 
                                      header = FALSE), 
                           use.names = FALSE)
  # get maxOffreprAct (then rm reprAct)
  if(platform == "windows"){
    cl <- makeCluster(n_cores)
    # for all, since all were designed for >0 perfect sites
    maxOffreprAct <- unlist(parLapply(cl, split(reprAct, candidate_sgID), function(sg_ra){
      max(sg_ra[-which.max(sg_ra)])
    })) # names correspond to candidate_sgID, checked
    stopCluster(cl)
  } else{
    maxOffreprAct <- unlist(mclapply(split(reprAct, candidate_sgID), function(sg_ra){
      max(sg_ra[-which.max(sg_ra)])
    }, mc.cores = n_cores))
  }
  # get which site is maxOffreprAct for each candidate sgRNA
  if(platform == "windows"){
    cl <- makeCluster(n_cores)
    whichmaxOffreprAct <- unlist(parLapply(cl, split(reprAct, candidate_sgID), function(sg_ra){
      order(sg_ra, decreasing = TRUE)[2]
    })) # names correspond to candidate_sgID, checked
    stopCluster(cl)
  } else{
    whichmaxOffreprAct <- unlist(mclapply(split(reprAct, candidate_sgID), function(sg_ra){
      order(sg_ra, decreasing = TRUE)[2]
    }, mc.cores = n_cores))
  }
  rm(reprAct)
  # get mismatch strings
  candidate_mm <- unlist(read.table(paste0(outdir, "/missmatch_matrix_end.txt"), 
                                    colClasses = c("NULL", "character", "NULL"), 
                                    sep = "\t", 
                                    header = FALSE), 
                         use.names = FALSE)
  # invert directly to 5'-3' (could multi-core)
  #maxOffmm <- unlist(lapply(candidate_mm[whichmaxOffreprAct], function(x){intToUtf8(rev(utf8ToInt(x)))}))
  maxOffmm <- mapply(function(sg, ind){intToUtf8(rev(utf8ToInt(sg[ind])))}, 
                     sg = split(candidate_mm, candidate_sgID), 
                     ind = whichmaxOffreprAct)
  # rm from RAM
  rm(candidate_mm)
  # make data frame
  all_df <- data.frame("target" = sub("\\_gR[0-9]+[rf]$", "", names(candidate_sgRNAs_uNT[as.integer(names(maxOffreprAct))])), 
                       "target_name" = ifelse(is.na(genes$gene[match(sub("\\_gR[0-9]+[rf]$", "", names(candidate_sgRNAs_uNT[as.integer(names(maxOffreprAct))])), genes$locus_tag)]), 
                                              sub("\\_gR[0-9]+[rf]$", "", names(candidate_sgRNAs_uNT[as.integer(names(maxOffreprAct))])), 
                                              genes$gene[match(sub("\\_gR[0-9]+[rf]$", "", names(candidate_sgRNAs_uNT[as.integer(names(maxOffreprAct))])), genes$locus_tag)]), 
                       "sgRNA_name" = names(candidate_sgRNAs_uNT[as.integer(names(maxOffreprAct))]), 
                       "sgRNA_seq" = substr(candidate_sgRNAs_uNT[as.integer(names(maxOffreprAct))], 1, sum(regions)), 
                       "sgRNA_PAM" = substr(candidate_sgRNAs_uNT[as.integer(names(maxOffreprAct))], sum(regions) + 1, sum(regions) + nchar(PAM)), 
                       "sgRNA_GC" = as.numeric(letterFrequency(narrow(candidate_sgRNAs_uNT[as.integer(names(maxOffreprAct))], 1, sum(regions)), letters = "GC") / sum(regions)))
  all_df$targetlength_bp <- genes$end[match(all_df$target, genes$locus_tag)] - genes$start[match(all_df$target, genes$locus_tag)] + 1
  ## for dist2CS
  coord_tmp_on <- as.integer(sub(".*_gR([0-9]+).*", "\\1", all_df$sgRNA_name))
  strand_on <- substr(all_df$sgRNA_name, nchar(all_df$sgRNA_name), nchar(all_df$sgRNA_name))
  start_on <- ifelse(strand_on == "f", 
                     coord_tmp_on - (sum(regions) - offset_emp - 1), 
                     coord_tmp_on - offset_emp)
  end_on <- start_on + sum(regions) - 1
  ##
  all_df$dist2SC_bp <- ifelse(strand_on == "f", all_df$targetlength_bp - end_on, start_on - 1)
  all_df$dist2SC_rel <- all_df$dist2SC_bp / (all_df$targetlength_bp - sum(regions))
  all_df$maxOffreprAct <- maxOffreprAct
  ## add perfect matching site info
  # read in site IDs
  cand_site_ID <- unlist(read.table(paste0(outdir, "/missmatch_matrix_end.txt"), 
                                    colClasses = c("NULL", "NULL", "integer"), 
                                    sep = "\t", 
                                    header = FALSE), 
                         use.names = FALSE)
  # get sequence of maxOffreprAct site
  maxoff_seq <- all_sites[cand_site_ID[whichmaxOffreprAct]]
  # select zero-mismatch IDs
  candidate_sgID_sel <- candidate_sgID[candidate_ID]
  cand_site_ID_sel <- cand_site_ID[candidate_ID]
  # rm rest
  #rm(cand_site_ID)
  # compute binding site scores
  if(platform == "windows"){
    cl <- makeCluster(n_cores)
    clusterExport(cl, c("genes", "regions", "offset_emp"))
    NTgene_ls <- parLapply(cl, names(all_sites[cand_site_ID_sel]), function(site){
      chrom_tmp <- gsub(" ", "", sub("\\_gR[0-9]+[rf]$", "", site))
      strand_tmp <- ifelse(substr(site, nchar(site), nchar(site)) == "f", "+", "-")
      coord_tmp <- as.integer(sub(".*_gR([0-9]+).*", "\\1", site))
      start_tmp <- ifelse(strand_tmp == "+", 
                          yes = coord_tmp - (sum(regions) - offset_emp - 1), 
                          no = coord_tmp - offset_emp)
      end_tmp <- start_tmp + sum(regions) - 1
      tmp_i <- start_tmp - genes$end <= 0 & 
        end_tmp - genes$start >= 0 & 
        strand_tmp != genes$strand & 
        unlist(lapply(gsub(" ", "", genes$seqid), grepl, x = chrom_tmp))
      # find gene hits on site, if any (can be 0, or >1 if overlap)
      if(any(tmp_i)){
        c(paste(genes$locus_tag[which(tmp_i)], collapse = ","), 
          paste(ifelse(is.na(genes$gene[which(tmp_i)]), genes$locus_tag[which(tmp_i)], genes$gene[which(tmp_i)]), collapse = ","), 
          paste(chrom_tmp, collapse = ","), 
          paste(strand_tmp, collapse = ","), 
          mapply(paste, start_tmp, end_tmp, MoreArgs = list(sep = "..", collapse = ",")))
      } else{
        c(NA, 
          NA, 
          paste(chrom_tmp, collapse = ","), 
          paste(strand_tmp, collapse = ","), 
          mapply(paste, start_tmp, end_tmp, MoreArgs = list(sep = "..", collapse = ",")))
      }
    })
    # split per sgRNA
    hitinfo <- do.call(rbind, parLapply(cl, split(NTgene_ls, candidate_sgID_sel), function(hits){
      # cannot paste with sep or col = ","...
      gsub(" ", ",", do.call(paste, hits))
    }))
    stopCluster(cl)
  } else{
    NTgene_ls <- mclapply(names(all_sites[cand_site_ID_sel]), function(site){
      chrom_tmp <- gsub(" ", "", sub("\\_gR[0-9]+[rf]$", "", site))
      strand_tmp <- ifelse(substr(site, nchar(site), nchar(site)) == "f", "+", "-")
      coord_tmp <- as.integer(sub(".*_gR([0-9]+).*", "\\1", site))
      start_tmp <- ifelse(strand_tmp == "+", 
                          yes = coord_tmp - (sum(regions) - offset_emp - 1), 
                          no = coord_tmp - offset_emp)
      end_tmp <- start_tmp + sum(regions) - 1
      tmp_i <- start_tmp - genes$end <= 0 & 
        end_tmp - genes$start >= 0 & 
        strand_tmp != genes$strand & 
        unlist(lapply(gsub(" ", "", genes$seqid), grepl, x = chrom_tmp))
      # find gene hits on site, if any (can be 0, or >1 if overlap)
      if(any(tmp_i)){
        c(paste(genes$locus_tag[which(tmp_i)], collapse = ","), 
          paste(ifelse(is.na(genes$gene[which(tmp_i)]), genes$locus_tag[which(tmp_i)], genes$gene[which(tmp_i)]), collapse = ","), 
          paste(chrom_tmp, collapse = ","), 
          paste(strand_tmp, collapse = ","), 
          mapply(paste, start_tmp, end_tmp, MoreArgs = list(sep = "..", collapse = ",")))
      } else{
        c(NA, 
          NA, 
          paste(chrom_tmp, collapse = ","), 
          paste(strand_tmp, collapse = ","), 
          mapply(paste, start_tmp, end_tmp, MoreArgs = list(sep = "..", collapse = ",")))
      }
    }, mc.cores = n_cores)
    # split per sgRNA
    hitinfo <- do.call(rbind, mclapply(split(NTgene_ls, candidate_sgID_sel), function(hits){
      # cannot paste with sep or coll = ","...
      gsub(" ", ",", do.call(paste, hits))
    }, mc.cores = n_cores))
  }
  #rm(candidate_sgID_sel)
  colnames(hitinfo) <- c("all_targets", "all_target_names", "all_chroms", "all_strands", "all_ranges")
  ## flag cut sites
  spacers_ext <- xscat(oligoForwardOverhang, 
                       narrow(candidate_sgRNAs_uNT[as.integer(names(maxOffreprAct))], 
                              1, sum(regions)), 
                       reverseComplement(DNAString(oligoReverseOverhang)))
  cut_sites_plusrev <- c(reverseComplement(DNAStringSet(cut_sites)), cut_sites)
  names(cut_sites_plusrev) <- c(paste0("reverse_", names(cut_sites)), names(cut_sites))
  cut_site_found <- vwhichPDict(cut_sites_plusrev, spacers_ext)
  sgRNA_cut_site <- unlist(lapply(cut_site_found, function(sg_cut){
    ifelse(isEmpty(sg_cut), NA, paste(names(cut_sites_plusrev[sg_cut]), collapse = ","))
  }))
  ## flag non-standard nt
  nt_freq <- alphabetFrequency(narrow(candidate_sgRNAs_uNT[as.integer(names(maxOffreprAct))], 
                                      1, sum(regions)))
  nonbase_flag <- unlist(lapply(seq.int(nrow(nt_freq)), function(x){
    ifelse(any(nt_freq[x, -(1:4)] > 0), 
           paste(colnames(nt_freq)[-(1:4)][nt_freq[x, -(1:4)] > 0], collapse = ","), 
           NA)
  }))
  ## flag bad seeds 
  # cannot be multiple per sgRNA, use that
  bad_seed_found <- lapply(bad_seeds, function(seed){
    which(narrow(candidate_sgRNAs_uNT[as.integer(names(maxOffreprAct))], 
                 sum(regions) - nchar(seed) + 1, 
                 sum(regions)) %in% seed)
  })
  sgRNA_bad_seed <- rep(NA, nrow(all_df))
  for(i in seq.int(bad_seeds)){
    sgRNA_bad_seed[bad_seed_found[[i]]] <- bad_seeds[i]
  }
  # bind into df
  all_df <- data.frame(all_df, 
                       hitinfo[match(names(maxOffreprAct), as.integer(rownames(hitinfo))), ], 
                       "oligoForward" = paste0(oligoForwardOverhang, all_df$sgRNA_seq),
                       "oligoReverse" = paste0(oligoReverseOverhang, as.character(reverseComplement(DNAStringSet(all_df$sgRNA_seq)))), 
                       "cut_site" = sgRNA_cut_site, 
                       "bad_seed" = sgRNA_bad_seed, 
                       "maxOff_mismatches" = maxOffmm, 
                       "maxOff_seq" = as.character(maxoff_seq))
  # write to file
  write.csv(all_df, 
            file = paste0(outdir, "/", genomeID, "_sgRNAs_all.csv"), 
            row.names = FALSE)
  # remove relatively large objects
  rm(all_df, NTgene_ls, hitinfo, cand_site_ID, candidate_sgID_sel, cand_site_ID_sel, strand_on, start_on, end_on, coord_tmp_on)
}


#### 9. Select optimal sgRNA per feature ####
if(output_optimized_list){
  message(paste0(Sys.time(), ": Selecting optimal sgRNA per feature..."))
  # no need to compute if already done
  if(!output_all_candidates){
    # for maxOffreprAct need all
    candidate_sgID <- unlist(read.table(paste0(outdir, "/missmatch_matrix_end.txt"), 
                                        colClasses = c("integer", "NULL", "NULL"), 
                                        sep = "\t", 
                                        header = FALSE), 
                             use.names = FALSE)
    # get maxOffreprAct (then rm reprAct)
    if(platform == "windows"){
      cl <- makeCluster(n_cores)
      # for all, since all were designed for >0 perfect sites
      maxOffreprAct <- unlist(parLapply(cl, split(reprAct, candidate_sgID), function(sg_ra){
        max(sg_ra[-which.max(sg_ra)])
      })) # names correspond to candidate_sgID, checked
      stopCluster(cl)
    } else{
      maxOffreprAct <- unlist(mclapply(split(reprAct, candidate_sgID), function(sg_ra){
        max(sg_ra[-which.max(sg_ra)])
      }, mc.cores = n_cores))
    }
    if(platform == "windows"){
      cl <- makeCluster(n_cores)
      whichmaxOffreprAct <- unlist(parLapply(cl, split(reprAct, candidate_sgID), function(sg_ra){
        order(sg_ra, decreasing = TRUE)[2]
      })) # names correspond to candidate_sgID, checked
      stopCluster(cl)
    } else{
      whichmaxOffreprAct <- unlist(mclapply(split(reprAct, candidate_sgID), function(sg_ra){
        order(sg_ra, decreasing = TRUE)[2]
      }, mc.cores = n_cores))
    }
    rm(reprAct)
    # get mismatch strings
    candidate_mm <- unlist(read.table(paste0(outdir, "/missmatch_matrix_end.txt"), 
                                      colClasses = c("NULL", "character", "NULL"), 
                                      sep = "\t", 
                                      header = FALSE), 
                           use.names = FALSE)
    # invert directly to 5'-3' (could multi-core)
    #maxOffmm <- unlist(lapply(candidate_mm[whichmaxOffreprAct], function(x){intToUtf8(rev(utf8ToInt(x)))}))
    maxOffmm <- mapply(function(sg, ind){intToUtf8(rev(utf8ToInt(sg[ind])))}, 
                       sg = split(candidate_mm, candidate_sgID), 
                       ind = whichmaxOffreprAct)
    # rm from RAM
    rm(candidate_mm)
  }
  # sgRNA identifiers of perfect matches; take out duplicates
  candidate_sgID <- unique(candidate_sgID[candidate_ID])
  # get on-targets
  on_target <- sub("\\_gR[0-9]+[rf]$", "", names(candidate_sgRNAs_uNT[candidate_sgID]))
  # process cut_sites (before bad seeds, this is more important)
  if(cut_site_rule == "avoid"){
    # only look in spacers
    spacers <- narrow(candidate_sgRNAs_uNT[candidate_sgID], 1, sum(regions))
    # but plus overhangs; both strands; same as rev.compl. on one strand
    spacers_ext <- xscat(oligoForwardOverhang, spacers, reverseComplement(DNAString(oligoReverseOverhang)))
    # also check reverse complements
    cut_sites_plusrev <- c(reverseComplement(DNAStringSet(cut_sites)), cut_sites)
    # don't use pdict; can only have one width for all sites
    #cut_site_pdict <- PDict(cut_sites_plusrev)
    cut_site_found <- vwhichPDict(cut_sites_plusrev, spacers_ext)
    cutsite_ID <- unlist(lapply(cut_site_found, length)) != 0
    all_cs <- unlist(lapply(split(cutsite_ID, on_target), all))
    if(any(all_cs)){
      message(paste0("sgRNA with cut site selected for ", 
                     paste(names(all_cs[all_cs]), collapse = ", "), 
                     ". No alternative available."))
    }
    # redefine candidate sgRNA index including cut site processing
    # candidate_sgID & cutsite_ID same index, use!
    candidate_sgID <- unlist(lapply(split(seq.int(candidate_sgID), on_target), function(tar_cs_ind){
      if(!all(cutsite_ID[tar_cs_ind])){
        # exclude cut_site sgRNAs if alternatives available
        candidate_sgID[tar_cs_ind][!cutsite_ID[tar_cs_ind]]
      } else{
        # if different cut_sites in candidates, proceed to find "least bad" ones
        if(length(unique(cut_site_found[tar_cs_ind])) != 1){
          # base on order of input cut_sites; same index scores fwd/rev
          cs_index_sc <- lapply(cut_site_found[tar_cs_ind], function(x){
            # sort "worst" to "best" per sgRNA
            sort(ifelse(x > length(cut_sites), x - length(cut_sites), x))
          })
          # loop through sites per sgRNA, prioritizing less & "less bad" sites
          minimax_tmp <- optim_tmp <- seq.int(length(cs_index_sc))
          for(i in seq.int(max(unlist(lapply(cs_index_sc, length))))){
            score_tmp <- unlist(lapply(cs_index_sc[optim_tmp], function(x){x[i]}))
            if(any(is.na(score_tmp))){
              # is any: shorter, which is better (not possible in first loop since all had cut_site)
              minimax_tmp <- which(is.na(score_tmp))
            } else{
              # higher index score = "better"
              minimax_tmp <- which(score_tmp == max(score_tmp))
            }
            # relative ID for sgRNA IDs for this target
            optim_tmp <- optim_tmp[minimax_tmp]
          }
          # update IDs per target to include only minimax'ed cut_site sgRNAs
          candidate_sgID[tar_cs_ind][optim_tmp]
        } else{
          # if same cut_sites for all sgRNAs, no selection
          candidate_sgID[tar_cs_ind]
        }
      }
    }), 
    use.names = FALSE)
    # redefine on_target (easier than within same lapply loop)
    on_target <- sub("\\_gR[0-9]+[rf]$", "", names(candidate_sgRNAs_uNT[candidate_sgID]))
  }
  # process non-standard nucleotides
  if(avoidNonBaseNT){
    if(!hasOnlyBaseLetters(narrow(candidate_sgRNAs_uNT[candidate_sgID], 1, sum(regions)))){
      # # only have to change IDs of on_targets for which any candidate has non-base letter
      # nonbase_ID <- c(letterFrequency(narrow(candidate_sgRNAs_uNT[candidate_sgID], 1, sum(regions)), "ATGC") != sum(regions))
      # nonbase_tar <- on_target[nonbase_ID]
      # # only have to change IDs of targets of these IDs:
      # candidate_sgID[nonbase_ID]
      ##
      candidate_sgID <- unlist(lapply(split(seq.int(candidate_sgID), on_target), function(sg_seq_ID){
        # check only spacer sequences
        sg_seq <- narrow(candidate_sgRNAs_uNT[candidate_sgID[sg_seq_ID]], 1, sum(regions))
        # detect any non-base nt in all candidate spacers
        if(!hasOnlyBaseLetters(sg_seq)){
          # which candidate has non-base nt
          nonbaseletter_ID <- c(letterFrequency(sg_seq, "ATGC") != sum(regions))
          if(!all(nonbaseletter_ID)){
            candidate_sgID[sg_seq_ID][!nonbaseletter_ID]
          } else{
            message(paste0("sgRNA with non-standard nucleotide (ATGC) selected for feature ", on_target[sg_seq_ID][1]), 
                    ". No alternative available.")
            candidate_sgID[sg_seq_ID]
          }
        } else{
          # no non-base nt in any candidate
          candidate_sgID[sg_seq_ID]
        }
      }))
      on_target <- sub("\\_gR[0-9]+[rf]$", "", names(candidate_sgRNAs_uNT[candidate_sgID]))
    }
  }
  # process bad seeds
  if(bad_seed_rule == "avoid"){
    # no complete overlap of bad seeds allowed (also not logical)
    badseed_ID <- as.logical(rowSums(do.call(cbind, lapply(bad_seeds, function(seed){
      substr(candidate_sgRNAs_uNT[candidate_sgID], 
             sum(regions) - nchar(seed) + 1, 
             sum(regions)) %in% seed
    }))))
    # for any feature only bad seed sgRNAs?
    all_bs <- unlist(lapply(split(badseed_ID, on_target), all))
    if(any(all_bs)){
      message(paste0("Bad seed selected for ", paste(names(all_bs[all_bs]), collapse = ", "), 
                     ". No alternative available."))
    }
    # redefine candidate sgRNA index including bad seed processing
    # candidate_sgID & badseed_ID same index, use!
    candidate_sgID <- unlist(lapply(split(seq.int(candidate_sgID), on_target), function(tar_bs_ind){
      if(!all(badseed_ID[tar_bs_ind])){
        candidate_sgID[tar_bs_ind][!badseed_ID[tar_bs_ind]]
      } else{
        # if all bad seed for target, do not exclude
        candidate_sgID[tar_bs_ind]
      }
    }), use.names = FALSE)
    # redefine on_target (easier than within same lapply loop)
    on_target <- sub("\\_gR[0-9]+[rf]$", "", names(candidate_sgRNAs_uNT[candidate_sgID]))
  }
  # both selection rounds in one go; load only required IDs & names onto nodes
  names(candidate_sgID) <- names(candidate_sgRNAs_uNT[candidate_sgID])
  if(platform == "windows"){
    cl <- makeCluster(n_cores)
    clusterExport(cl, c("maxOffreprAct", "errorRange_maxOffreprAct", "n_sgRNA"))
    #minimax_reprAct
    optcand_name <- parLapply(cl, split(candidate_sgID, on_target), function(sel_sgID){
      # compute relevant scores
      target_mo_ra <- maxOffreprAct[match(sel_sgID, names(maxOffreprAct))]
      strand_tmp <- substr(names(sel_sgID[1]), nchar(names(sel_sgID[1])), nchar(names(sel_sgID[1])))
      coord_tmp <- as.integer(sub(".*_gR([0-9]+).*", "\\1", names(sel_sgID)))
      # prep loop
      cand_nr <- min(n_sgRNA, length(target_mo_ra))
      cand_sel <- integer(cand_nr)
      cand_ind_tmp <- rep(TRUE, length(target_mo_ra))
      for(cand in seq.int(cand_nr)){
        # filter 1: minimax off-target reprAct
        cand_ind_filt1 <- target_mo_ra[cand_ind_tmp] <= min(target_mo_ra[cand_ind_tmp]) + errorRange_maxOffreprAct
        # filter 2: min dist2SC
        mindist <- switch(strand_tmp, 
                          "r" = which.min(coord_tmp[cand_ind_tmp][cand_ind_filt1]), 
                          "f" = which.max(coord_tmp[cand_ind_tmp][cand_ind_filt1]))
        # store selected
        cand_sel[cand] <- names(sel_sgID[cand_ind_tmp][cand_ind_filt1][mindist])
        # update candidates (rm selected sgRNA)
        cand_ind_tmp[cand_ind_filt1][mindist] <- FALSE
      }
      return(cand_sel)
      # could here do
      #target_mo_ra[target_mo_ra <= (sort(target_mo_ra)[n_sgRNA] + errorRange_maxOffreprAct)]
      #sel_sgID[target_mo_ra <= (sort(target_mo_ra)[min(n_sgRNA, length(target_mo_ra))] + errorRange_maxOffreprAct)]
      # instead of
      #sel_sgID[target_mo_ra <= min(target_mo_ra) + errorRange_maxOffreprAct]
    })
    stopCluster(cl)
  } else{
    optcand_name <- mclapply(split(candidate_sgID, on_target), function(sel_sgID){
      # compute relevant scores
      target_mo_ra <- maxOffreprAct[match(sel_sgID, names(maxOffreprAct))]
      strand_tmp <- substr(names(sel_sgID[1]), nchar(names(sel_sgID[1])), nchar(names(sel_sgID[1])))
      coord_tmp <- as.integer(sub(".*_gR([0-9]+).*", "\\1", names(sel_sgID)))
      # prep loop
      cand_nr <- min(n_sgRNA, length(target_mo_ra))
      cand_sel <- integer(cand_nr)
      cand_ind_tmp <- rep(TRUE, length(target_mo_ra))
      for(cand in seq.int(cand_nr)){
        # filter 1: minimax off-target reprAct
        cand_ind_filt1 <- target_mo_ra[cand_ind_tmp] <= min(target_mo_ra[cand_ind_tmp]) + errorRange_maxOffreprAct
        # filter 2: min dist2SC
        mindist <- switch(strand_tmp, 
                          "r" = which.min(coord_tmp[cand_ind_tmp][cand_ind_filt1]), 
                          "f" = which.max(coord_tmp[cand_ind_tmp][cand_ind_filt1]))
        # store selected
        cand_sel[cand] <- names(sel_sgID[cand_ind_tmp][cand_ind_filt1][mindist])
        # update candidates (rm selected sgRNA)
        cand_ind_tmp[cand_ind_filt1][mindist] <- FALSE
      }
      return(cand_sel)
      # target_mo_ra <- maxOffreprAct[match(sel_sgID, names(maxOffreprAct))]
      # sel_sgID[target_mo_ra <= min(target_mo_ra) + errorRange_maxOffreprAct]
    }, mc.cores = n_cores)
  }
  # # get names (much faster than lapply over minimax_reprAct to get names)
  # # use here names(minimax_reprAct) instead of minimax_reprAct 1st line:
  # cand_names <- split(names(candidate_sgRNAs_uNT[unlist(minimax_reprAct)]), 
  #                     unlist(mapply(rep, 
  #                                   names(minimax_reprAct), 
  #                                   unlist(lapply(minimax_reprAct, length)))))
  # # dist2SC computation (don't load large objects on nodes)
  # if(platform == "windows"){
  #   cl <- makeCluster(n_cores)
  #   optcand_name <- parLapply(cl, cand_names, function(tar_cand){
  #     strand_tmp <- substr(tar_cand[1], nchar(tar_cand), nchar(tar_cand))
  #     coord_tmp <- as.integer(sub(".*_gR([0-9]+).*", "\\1", tar_cand))
  #     mindist <- switch(strand_tmp, 
  #                       "r" = which.min(coord_tmp), 
  #                       "f" = which.max(coord_tmp))
  #     tar_cand[mindist]
  #   })
  #   stopCluster(cl)
  # } else{
  #   optcand_name <- mclapply(cand_names, function(tar_cand){
  #     strand_tmp <- substr(tar_cand[1], nchar(tar_cand), nchar(tar_cand))
  #     coord_tmp <- as.integer(sub(".*_gR([0-9]+).*", "\\1", tar_cand))
  #     mindist <- switch(strand_tmp, 
  #                       "r" = which.min(coord_tmp), 
  #                       "f" = which.max(coord_tmp))
  #     tar_cand[mindist]
  #   }, mc.cores = n_cores)
  # }
  ## compute scores for selected sgRNAs
  sgRNA_seq <- substr(candidate_sgRNAs_uNT[match(unlist(optcand_name), names(candidate_sgRNAs_uNT))], 
                      1, sum(regions))
  sgRNA_PAM <- substr(candidate_sgRNAs_uNT[match(unlist(optcand_name), names(candidate_sgRNAs_uNT))], 
                      sum(regions) + 1, sum(regions) + nchar(PAM))
  sgRNA_GC <- as.numeric(letterFrequency(narrow(candidate_sgRNAs_uNT[match(unlist(optcand_name), names(candidate_sgRNAs_uNT))], 
                                                1, sum(regions)), letters = "GC") / sum(regions))
  # sgRNA_GC <- unlist(lapply(strsplit(sgRNA_seq, NULL), function(nt_site){
  #   sum(nt_site %in% c("C", "G")) / sum(regions)
  # }))
  ## add cut_site flags
  spacers_ext <- xscat(oligoForwardOverhang, 
                       narrow(candidate_sgRNAs_uNT[match(unlist(optcand_name), names(candidate_sgRNAs_uNT))], 
                              1, sum(regions)), 
                       reverseComplement(DNAString(oligoReverseOverhang)))
  cut_sites_plusrev <- c(reverseComplement(DNAStringSet(cut_sites)), cut_sites)
  names(cut_sites_plusrev) <- c(paste0("reverse_", names(cut_sites)), names(cut_sites))
  cut_site_found <- vwhichPDict(cut_sites_plusrev, spacers_ext)
  sgRNA_cut_site <- unlist(lapply(cut_site_found, function(sg_cut){
    ifelse(isEmpty(sg_cut), NA, paste(names(cut_sites_plusrev[sg_cut]), collapse = ","))
  }))
  ## flag bad seeds 
  # cannot be multiple per sgRNA, use that
  bad_seed_found <- lapply(bad_seeds, function(seed){
    which(narrow(candidate_sgRNAs_uNT[match(unlist(optcand_name), names(candidate_sgRNAs_uNT))], 
                 sum(regions) - nchar(seed) + 1, 
                 sum(regions)) %in% seed)
  })
  #sgRNA_bad_seed <- rep(NA, length(optcand_name))
  sgRNA_bad_seed <- rep(NA, length(unlist(optcand_name)))
  for(i in seq.int(bad_seeds)){
    sgRNA_bad_seed[bad_seed_found[[i]]] <- bad_seeds[i]
  }
  ## flag non-standard nt
  nt_freq <- alphabetFrequency(narrow(candidate_sgRNAs_uNT[match(unlist(optcand_name), names(candidate_sgRNAs_uNT))], 
                                      1, sum(regions)))
  nonbase_flag <- unlist(lapply(seq.int(nrow(nt_freq)), function(x){
    ifelse(any(nt_freq[x, -(1:4)] > 0), 
           paste(colnames(nt_freq)[-(1:4)][nt_freq[x, -(1:4)] > 0], collapse = ","), 
           NA)
  }))
  ## feature length
  tar_rep <- rep(names(optcand_name), unlist(lapply(optcand_name, length)))
  targetlength <- genes$end[match(tar_rep, genes$locus_tag)] - genes$start[match(tar_rep, genes$locus_tag)] + 1
  ## for dist2CS
  coord_tmp_on <- as.integer(sub(".*_gR([0-9]+).*", "\\1", unlist(optcand_name)))
  strand_on <- substr(unlist(optcand_name), 
                      nchar(unlist(optcand_name)), 
                      nchar(unlist(optcand_name)))
  start_on <- ifelse(strand_on == "f", 
                     coord_tmp_on - (sum(regions) - offset_emp - 1), 
                     coord_tmp_on - offset_emp)
  end_on <- start_on + sum(regions) - 1
  dist2SC_bp <- ifelse(strand_on == "f", 
                       targetlength - end_on, 
                       start_on - 1)
  dist2SC_rel <- dist2SC_bp / (targetlength - sum(regions))
  # dist2SC_rel <- ifelse(strand_on == "r", 
  #                       dist2SC_bp / (targetlength - sum(regions)), 
  #                       1 - (dist2SC_bp / (targetlength - sum(regions))))
  ##
  sel_sg_ID <- match(unlist(optcand_name), names(candidate_sgRNAs_uNT))
  maxOffreprAct_opt <- maxOffreprAct[match(sel_sg_ID, names(maxOffreprAct))]
  rm(maxOffreprAct)
  maxOffmm_opt <- maxOffmm[match(sel_sg_ID, names(maxOffmm))]
  ## for all_targets
  # read in sgRNA_IDs
  cand_sg_ID <- unlist(read.table(paste0(outdir, "/missmatch_matrix_end.txt"), 
                                  colClasses = c("integer", "NULL", "NULL"), 
                                  sep = "\t", 
                                  header = FALSE), 
                       use.names = FALSE)
  # find which lines have 0-mm && occur in optimized list
  cand_sg_ID_sel <- cand_sg_ID[candidate_ID] %in% sel_sg_ID
  # rm unwanted IDs
  cand_sg_ID <- cand_sg_ID[candidate_ID][cand_sg_ID_sel]
  # read in site IDs
  cand_site_ID <- unlist(read.table(paste0(outdir, "/missmatch_matrix_end.txt"), 
                                    colClasses = c("NULL", "NULL", "integer"), 
                                    sep = "\t", 
                                    header = FALSE), 
                         use.names = FALSE)
  # get sequence of maxOffreprAct site
  maxOffseq_opt <- all_sites[cand_site_ID[whichmaxOffreprAct][match(sel_sg_ID, names(maxOffmm))]]
  rm(maxOffmm)
  # select as sites as for cand_sg_ID_sel
  cand_site_ID_sel <- cand_site_ID[candidate_ID][cand_sg_ID_sel]
  # rm rest
  rm(cand_site_ID)
  # compute binding site scores
  chrom_tmp <- gsub(" ", "", 
                    sub("\\_gR[0-9]+[rf]$", "", names(all_sites[cand_site_ID_sel])))
  strand_tmp <- ifelse(substr(names(all_sites[cand_site_ID_sel]), 
                              nchar(names(all_sites[cand_site_ID_sel])), 
                              nchar(names(all_sites[cand_site_ID_sel]))) == "f", 
                       "+", "-")
  coord_tmp <- as.numeric(sub(".*_gR([0-9]+).*", "\\1", names(all_sites[cand_site_ID_sel])))
  start_tmp <- ifelse(strand_tmp == "+", 
                      yes = coord_tmp - (sum(regions) - offset_emp - 1), 
                      no = coord_tmp - offset_emp)
  end_tmp <- start_tmp + sum(regions) - 1
  if(platform == "windows"){
    cl <- makeCluster(n_cores)
    clusterExport(cl, c("chrom_tmp", "strand_tmp", "start_tmp", "end_tmp", "genes"))
    NTgene_ls <- parLapply(cl, seq.int(cand_site_ID_sel), function(site){
      tmp_i <- start_tmp[site] - genes$end <= 0 & 
        end_tmp[site] - genes$start >= 0 & 
        strand_tmp[site] != genes$strand & 
        unlist(lapply(gsub(" ", "", genes$seqid), grepl, x = chrom_tmp[site]))
      # find gene hits on site, if any (can be 0, or >1 if overlap)
      if(any(tmp_i)){
        c(paste(genes$locus_tag[which(tmp_i)], collapse = ","), 
          paste(ifelse(is.na(genes$gene[which(tmp_i)]), genes$locus_tag[which(tmp_i)], genes$gene[which(tmp_i)]), collapse = ","), 
          paste(chrom_tmp[site], collapse = ","), 
          paste(strand_tmp[site], collapse = ","), 
          mapply(paste, start_tmp[site], end_tmp[site], MoreArgs = list(sep = "..", collapse = ",")))
      } else{
        c(NA, 
          NA, 
          paste(chrom_tmp[site], collapse = ","), 
          paste(strand_tmp[site], collapse = ","), 
          mapply(paste, start_tmp[site], end_tmp[site], MoreArgs = list(sep = "..", collapse = ",")))
      }
    })
    # split per sgRNA
    hitinfo <- do.call(rbind, parLapply(cl, split(NTgene_ls, cand_sg_ID), function(hits){
      # cannot paste with sep or coll = ","...
      gsub(" ", ",", do.call(paste, hits))
    }))
    stopCluster(cl)
  } else{
    NTgene_ls <- mclapply(seq.int(cand_site_ID_sel), function(site){
      tmp_i <- start_tmp[site] - genes$end <= 0 & 
        end_tmp[site] - genes$start >= 0 & 
        strand_tmp[site] != genes$strand & 
        unlist(lapply(gsub(" ", "", genes$seqid), grepl, x = chrom_tmp[site]))
      # find gene hits on site, if any (can be 0, or >1 if overlap)
      if(any(tmp_i)){
        c(paste(genes$locus_tag[which(tmp_i)], collapse = ","), 
          paste(ifelse(is.na(genes$gene[which(tmp_i)]), genes$locus_tag[which(tmp_i)], genes$gene[which(tmp_i)]), collapse = ","), 
          paste(chrom_tmp[site], collapse = ","), 
          paste(strand_tmp[site], collapse = ","), 
          mapply(paste, start_tmp[site], end_tmp[site], MoreArgs = list(sep = "..", collapse = ",")))
      } else{
        c(NA, 
          NA, 
          paste(chrom_tmp[site], collapse = ","), 
          paste(strand_tmp[site], collapse = ","), 
          mapply(paste, start_tmp[site], end_tmp[site], MoreArgs = list(sep = "..", collapse = ",")))
      }
    }, mc.cores = n_cores)
    # split per sgRNA
    hitinfo <- do.call(rbind, mclapply(split(NTgene_ls, cand_sg_ID), function(hits){
      # cannot paste with sep or coll = ","...
      gsub(" ", ",", do.call(paste, hits))
    }, mc.cores = n_cores))
  }
  rm(cand_sg_ID)
  colnames(hitinfo) <- c("all_targets", "all_target_names", "all_chroms", "all_strands", "all_ranges")
  ##
  oligoForward <- paste0(oligoForwardOverhang, sgRNA_seq)
  oligoReverse <- paste0(oligoReverseOverhang, as.character(reverseComplement(DNAStringSet(sgRNA_seq))))
  ## make data frame
  opt_df <- data.frame("target" = tar_rep, 
                       "target_name" = ifelse(is.na(genes$gene[match(tar_rep, genes$locus_tag)]), 
                                              genes$locus_tag[match(tar_rep, genes$locus_tag)], 
                                              genes$gene[match(tar_rep, genes$locus_tag)]), 
                       "sgRNA_name" = unlist(optcand_name), 
                       "sgRNA_seq" = sgRNA_seq, 
                       "sgRNA_PAM" = sgRNA_PAM, 
                       "sgRNA_GC" = sgRNA_GC, 
                       "targetlength_bp" = targetlength, 
                       "dist2SC_bp" = dist2SC_bp, 
                       "dist2SC_rel" = dist2SC_rel, 
                       "maxOffreprAct" = maxOffreprAct_opt, 
                       hitinfo[match(unlist(optcand_name), names(candidate_sgRNAs_uNT[as.integer(rownames(hitinfo))])), ], 
                       "oligoForward" = oligoForward, 
                       "oligoReverse" = oligoReverse, 
                       "cut_site" = sgRNA_cut_site, 
                       "nonstandard_nt" = nonbase_flag, 
                       "bad_seed" = sgRNA_bad_seed, 
                       "maxOff_mismatches" = maxOffmm_opt, 
                       "maxOff_seq" = as.character(maxOffseq_opt))
  # exclude duplicate sgRNAs from file (same spacer+PAM, multiple targets)
  if(filter_out_duplicates){
    opt_df <- opt_df[!duplicated(opt_df$sgRNA_seq), ]
    #opt_df <- opt_df[!duplicated(paste0(opt_df$sgRNA_seq, opt_df$sgRNA_PAM)), ]
  }
  # write to file
  write.csv(opt_df, 
            file = paste0(outdir, "/", genomeID, "_sgRNAs_optimal.csv"), 
            row.names = FALSE)
  opt_designed <- length(unique(opt_df$sgRNA_seq))
  opt_unique <- length(unique(unlist(lapply(strsplit(opt_df$all_targets, ","), function(tars){tars[tars != "NA"]}))))
  message(paste(opt_designed, 
                "sgRNAs designed to target", 
                opt_unique, 
                "out of", 
                nrow(genes), 
                "unique features."))
  # remove relatively large objects
  rm(opt_df, on_target, NTgene_ls, hitinfo, cand_sg_ID_sel, optcand_name, candidate_sgID)
}


#### 10. Compile full table if requested ####
if(output_full_list){
  message(paste0(Sys.time(), ": compiling full output table (this may take a while and require much RAM, depending on settings)"))
  # reading whole file 
  full_df <- read.table(paste0(outdir, "/missmatch_matrix_end.txt"), 
                        header = FALSE, sep = "\t")
  colnames(full_df) <- c("sgRNA_ID", "mismatches", "site_ID")
  # add n_mismatch
  full_df$n_mismatch <- scan(paste0(outdir, "/tmp_nmm.tsv"), 
                             what = integer(), quiet = TRUE)
  # add reprAct
  full_df$reprAct <- scan(paste0(outdir, "/tmp_reprAct.tsv"), 
                          quiet = TRUE)
  full_df$chrom <- gsub(" ", "", sub("\\_gR[0-9]+[rf]$", "", names(all_sites[full_df$site_ID])))
  full_df$strand <- ifelse(substr(names(all_sites[full_df$site_ID]), 
                                  nchar(names(all_sites[full_df$site_ID])), 
                                  nchar(names(all_sites[full_df$site_ID]))) == "f", 
                           "+", "-")
  coord_tmp <- as.numeric(sub(".*_gR([0-9]+).*", "\\1", names(all_sites[full_df$site_ID])))
  full_df$start <- ifelse(full_df$strand == "+", 
                      yes = coord_tmp - (sum(regions) - offset_emp - 1), 
                      no = coord_tmp - offset_emp)
  rm(coord_tmp)
  full_df$end <- full_df$start + sum(regions) - 1
  full_df$sgRNA_seq <- substr(candidate_sgRNAs_uNT[full_df$sgRNA_ID], 1, sum(regions))
  full_df$site_seq <- substr(all_sites[full_df$site_ID], 1, sum(regions))
  full_df$site_PAM <- substr(all_sites[full_df$site_ID], sum(regions) + 1, sum(regions) + nchar(PAM))
  full_df$sgRNA_GC <- as.numeric(letterFrequency(candidate_sgRNAs_uNT[full_df$sgRNA_ID], letters = "GC") / sum(regions))
  full_df$target <- sub("\\_gR[0-9]+[rf]$", "", names(candidate_sgRNAs_uNT[full_df$sgRNA_ID]))
  # need lots of RAM for detecting non-template strand gene overlaps of sites
  if(detect_offtarget_genes_full){
    if(platform == "windows"){
      cl <- makeCluster(n_cores)
      clusterExport(cl, "genes")
      full_df$NTgene <- unlist(parApply(cl, full_df[, c("chrom", "strand", "start", "end")], 1, function(site){
        tmp_i <- as.integer(site["start"]) - genes$end <= 0 & 
          as.integer(site["end"]) - genes$start >= 0 & 
          site["strand"] != genes$strand & 
          unlist(lapply(gsub(" ", "", genes$seqid), grepl, x = site["chrom"]))
        # find gene hits on site, if any (can be 0, or >1 if overlap)
        if(any(tmp_i)){
          paste(genes$locus_tag[which(tmp_i)], collapse = ",")
        } else{
          NA
        }
      }), use.names = FALSE)
      stopCluster(cl)
    } else{
      full_df$NTgene <- unlist(mcapply(full_df[, c("chrom", "strand", "start", "end")], 1, function(site){
        tmp_i <- as.integer(site["start"]) - genes$end <= 0 & 
          as.integer(site["end"]) - genes$start >= 0 & 
          site["strand"] != genes$strand & 
          unlist(lapply(gsub(" ", "", genes$seqid), grepl, x = site["chrom"]))
        # find gene hits on site, if any (can be 0, or >1 if overlap)
        if(any(tmp_i)){
          paste(genes$locus_tag[which(tmp_i)], collapse = ",")
        } else{
          NA
        }
      }, mc.cores = n_cores), use.names = FALSE)
    }
  }
  # write
  write.csv(full_df, 
            file = paste0(outdir, "/", genomeID, "_sgRNAs-sites_full.csv"), 
            row.names = FALSE)
  # can remove tmp files here
  file.remove(paste0(outdir, "/tmp_nmm.tsv"))
  file.remove(paste0(outdir, "/tmp_reprAct.tsv"))
  # if(platform == "windows"){
  #   ### reading separately, if run out of RAM (much, much slower though):
  #   # cl <- makeCluster(n_cores)
  #   # clusterExport(cl, "outdir")
  #   # clusterEvalQ(cl, con <- file(paste0(outdir, "/missmatch_matrix_end.txt"), "r"))
  #   # zeromm_ls <- do.call(rbind, parLapply(cl, candidate_ID[1:1000], function(line){
  #   #   # read.table(paste0(outdir, "/missmatch_matrix_end.txt"), 
  #   #   #            skip = line - 1, nrows = 1, 
  #   #   #            sep = "\t", header = FALSE)
  #   #   scan(con, what = "character", sep = "\t", 
  #   #        skip = line - 1, nlines = 1, n = 3, 
  #   #        quiet = TRUE)
  #   # }))
  #   # clusterEvalQ(cl, close(con))
  #   # stopCluster(cl)
  # } else{
  #   
  # }
}


#### 11. Output settings file ####
message(paste0(Sys.time(), ": writing settings file"))
end_time <- Sys.time()
settings_out <- c(invisible(capture.output(timestamp())), 
                  paste("Total running time:", capture.output(end_time - start_time)), 
                  ifelse(output_optimized_list, 
                         paste(opt_designed, 
                               "sgRNAs designed to target", 
                               opt_unique, 
                               "out of", 
                               nrow(genes), 
                               "unique features."), 
                         "No optimal sgRNA selection was performed."), 
                  "", 
                  "------------------",
                  "---- SETTINGS ----", 
                  "------------------",
                  "", 
                  paste("Input genome:", input_genome), 
                  paste("Output directory:", outdir),
                  paste("TINDRi.py location:", TINDRidir), 
                  paste("Number of sgRNAs to design per feature (if available):", n_sgRNA),
                  paste("NCBI genome directory:", path_ncbi_downloads),
                  paste("Feature type:", feature_type),
                  paste("Number of cores used:", n_cores),
                  paste("Intra-spacer region lengths (PAM-proximal to PAM-distal):", paste(regions, collapse = ",")),
                  paste("Total spacer length (in bp):", sum(regions)), 
                  paste("Cumulative maximum number of mismatches allowed per region for off-target detection:", paste(max_mismatch_cum, collapse = ",")), 
                  paste("So", paste("max.", max_mismatch_cum, "mismatche(s) in the first", cumsum(regions), "nucleotides", collapse = ", ")),
                  paste("Overall maximum number of mismatches allowed for off-target detection:", max(max_mismatch_cum)),
                  paste("Penalty scoring system:", ifelse(is.numeric(reprAct_penalties), "custom", reprAct_penalties)),
                  paste("Repression activity estimation penalties (PAM-proximal to PAM-distal):", paste(rev(penalties), collapse = ",")),
                  paste("Repression activity estimation function (penalties as input):", capture.output(pen_func)), 
                  paste("Maximum off-target repression activity error range:", errorRange_maxOffreprAct), 
                  paste("Selected bad seeds:", paste(bad_seeds, collapse = ", ")), 
                  paste("Bad seed decision rule:", bad_seed_rule, "selected bad seeds (equals 'ignore' if not 'exclude' or 'avoid')"),
                  paste("Selected cut sequences:", paste(cut_sites_plusrev, collapse = ", "), "(includes reverse complement of input sequences)"), 
                  paste("Cut sequence decision rule:", cut_site_rule, "selected cut sequences (e.g. restriction enzyme cut sites or tandem bp repeats) (equals 'ignore' if not 'exclude' or 'avoid')"),
                  paste("Oligo overhang forward:", oligoForwardOverhang),
                  paste("Oligo overhang reverse:", oligoReverseOverhang),
                  paste("PAM:", PAM),
                  paste("Output filtered list with one optimal sgRNA per target feature:", output_optimized_list), 
                  paste("Output fasta file feature target sequences:", output_target_fasta),
                  paste("Output fasta file all possible sgRNA sequences (zero-mismatch for at least one feature):", output_sgRNAs_fasta),
                  paste("Output fasta file all possible sgRNA binding site sequences in the genome:", output_sites_fasta),
                  paste("Output full list of candidate sgRNAs and corresponding binding sites (including mismatches):", output_full_list), 
                  paste("Detect for all identified binding sites if within gene on non-template strand, for full list output:", detect_offtarget_genes_full),
                  paste("Keep the temporary candidate sgRNA ID & sequence file:", keep_TINDRi_input_sgRNAs),
                  paste("Keep the temporary binding site ID & sequence file:", keep_TINDRi_input_sites),
                  paste("Keep the temporary candidate sgRNA - binding site nucleotide (mis)match matrix file:", keep_TINDRi_matches), 
                  paste("Used python instance path:", path_python), 
                  "", 
                  "------------------", 
                  "-- SESSION INFO --",
                  "------------------", 
                  "", 
                  capture.output(sessionInfo()), 
                  "", 
                  capture.output(py_config()))
settings_file <- file(paste0(outdir, "/", genomeID, "_settings.txt")) 
writeLines(settings_out, settings_file)
close(settings_file)
# remove file as not needed anymore
if(!keep_TINDRi_matches){
  file.remove(paste0(outdir, "/missmatch_matrix_end.txt"))
} else{
  file.rename(paste0(outdir, "/missmatch_matrix_end.txt"), 
              paste0(outdir, "/", genomeID, "_missmatch_matrix_end.txt"))
}
# DONE
message(paste0(Sys.time(), ": design pipeline ended. Output saved in: ", outdir))