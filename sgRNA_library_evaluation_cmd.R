#!/usr/bin/Rscript

#### Target site detection and reporting for given sgRNAs ####
# Author: Vincent de Bakker
# Veening Lab, DMF, FBM, University of Lausanne, Switzerland
# vincent.debakker@unil.ch
####

#### 1. Settings ####
if(!requireNamespace("optparse", quietly = TRUE)){install.packages("optparse")}
suppressPackageStartupMessages(library(optparse))
option_list <- list(make_option(c("-g", "--input_genome"),
                                type = "character", 
                                help = "Input genome either as GenBank file or NCBI assembly accession number"), 
                    make_option(c("-s", "--sgRNA_file"), 
                                type = "character", 
                                help = "Input .csv file listing sgRNA names and spacer sequences (without PAM)"), 
                    make_option(c("-o", "--outdir"), 
                                type = "character", default = "~/", 
                                help = "Output directory [default %default]"), 
                    make_option(c("-t", "--TINDRidir"), 
                                type = "character", default = "~/", 
                                help = "Directory in which python script TINDRi.py can be found [default %default]"), 
                    make_option(c("-d", "--path_ncbi_downloads"), 
                                type = "character", default = NA, 
                                help = "Directory in which to store files downloaded from NCBI, if --input_genome is an accession number [default %default]"), 
                    make_option(c("-r", "--regions"),
                                type = "character", default = "7,2,11", 
                                help = "Character string of comma-separated integers indicating sub-spacer region sizes (PAM-proximal to PAM-distal) for which to set --max_mismatch_cum [default %default]"), 
                    make_option(c("-m", "--max_mismatch_cum"),
                                type = "character", default = "1,2,11",
                                help = "Character string of comma-separated integers indicating maximum cumulutive (PAM-proximal to PAM-distal) 
                                number of mismatches allowed with potential target sites over sub-spacer regions as defined by --regions [default %default]"), 
                    make_option(c("-p", "--reprAct_penalties"), 
                                default = "HawkinsBsuMedian", 
                                help = "Nucleotide-wise (PAM-proximal to PAM-distal) penalty scores to compute estimated repression activity of sgRNAs on each site. 
                                One of characters HawkinsEcoMean, HawkinsBsuMean, HawkinsEcoMedian, HawkinsBsuMedian, Qi, QiMean, QiMedian, 
                                or a string of integers separated by comma's of length sum(--regions) with custom penalties [default %default]"), 
                    make_option(c("-f", "--pen_func"), 
                                type = "character", default = "prod", 
                                help = "Penalty function to apply to --reprAct_penalties to compute sgRNA repression activity estimate for given binding site [default %default]"), 
                    make_option(c("-b", "--bad_seeds"), 
                                type = "character", default = "ACCCA,TGGAA", 
                                help = "Bad seeds (can be of varying lengths) to detect, input as comma-separated string [default %default]"), 
                    make_option(c("-c", "--cut_sites"), 
                                type = "character", default = "CGTCTC,TTTTT", 
                                help = "Comma-separated string indicating motifs to detect (e.g. restriction enzyme cut sites or tandem base pairs); 
                                reverse-complements are also automatically detected [default %default]"), 
                    make_option(c("-n", "--cut_sites_names"), 
                                type = "character", default = "BsmBI,tandemT", 
                                help = "String with names of --cut_sites, input in same order and separated by comma's [default %default]"), 
                    make_option(c("-F", "--oligoForwardOverhang"), 
                                type = "character", default = "TATA", 
                                help = "Forward overhang for oligonucleotide primer design (important for --cut_sites detection) [default %default]"), 
                    make_option(c("-R", "--oligoReverseOverhang"), 
                                type = "character", default = "AAAC", 
                                help = "Reverse overhang for oligonucleotide primer design (important for --cut_sites detection) [default %default]"), 
                    make_option(c("-P", "--PAM"), 
                                type = "character", default = "NGG", 
                                help = "Protospacer Adjacent Motif (PAM) sequence [default %default]"), 
                    make_option(c("-N", "--no_output_summary"), 
                                action = "store_true", default = FALSE, 
                                help = "Specifiy flag to turn OFF summary output .csv file with first- and second-highest repression activity hits per sgRNA"), 
                    make_option(c("-A", "--output_full"), 
                                action = "store_true", default = FALSE, 
                                help = "Specify flag to turn ON full output .csv file with all found sgRNA - binding sites and scores"), 
                    make_option("--output_sgRNAs_fasta", 
                                action = "store_true", default = FALSE, 
                                help = "Specify flag to turn ON output .fasta file with all sgRNA spacer sequences"), 
                    make_option("--output_sites_fasta", 
                                action = "store_true", default = FALSE, 
                                help = "Specify flag to turn ON output .fasta file with all binding site sequences"), 
                    make_option("--detect_offtarget_genes_full", 
                                action = "store_true", default = FALSE, 
                                help = "Specify flag to turn ON detection of and adding information on genes annotated on each binding site, if --output_full is turned on"), 
                    make_option("--keep_TINDRi_input_sgRNAs", 
                                action = "store_true", default = FALSE, 
                                help = "Specify flag to keep intermittently produced .csv TINDRi.py input file with all sgRNA indexes and spacer sequences"), 
                    make_option("--keep_TINDRi_input_sites", 
                                action = "store_true", default = FALSE, 
                                help = "Specify flag to keep intermittently produced .csv TINDRi.py input file with all binding site indexes and sequences"), 
                    make_option(c("-K", "--keep_TINDRi_matches"), 
                                action = "store_true", default = FALSE, 
                                help = "Specify flag to keep intermittently produced .csv TINDRi.py output file with full sgRNA-site mismatch matrix"), 
                    make_option("--path_python", 
                                type = "character", default = NULL, 
                                help = "If python 3 is installed but not in the PATH variable, specify the full path to the program here"), 
                    make_option("--feature_type", 
                                type = "character", default = "locus_tag", 
                                help = "Name of feature flag to use for feature detection. Experimental feature, in principle do not change [default %default]"))
opt <- parse_args(OptionParser(option_list = option_list))

## REQUIRED ##
input_genome <- opt$input_genome
sgRNA_file <- opt$sgRNA_file
outdir <- opt$outdir
TINDRidir <- opt$TINDRidir

## OPTIONAL ##
path_ncbi_downloads <- opt$path_ncbi_downloads
regions <- as.numeric(strsplit(opt$regions, ",")[[1]])
max_mismatch_cum <- as.numeric(strsplit(opt$max_mismatch_cum, ",")[[1]])
reprAct_penalties <- opt$reprAct_penalties
pen_func <- opt$pen_func
bad_seeds <- strsplit(opt$bad_seeds, ",")[[1]]
cut_sites <- strsplit(opt$cut_sites, ",")[[1]]
names(cut_sites) <- strsplit(opt$cut_sites_names, ",")[[1]]
oligoForwardOverhang <- opt$oligoForwardOverhang
oligoReverseOverhang <- opt$oligoReverseOverhang
PAM <- opt$PAM
output_exact <- !opt$no_output_summary
output_full <- opt$output_full
output_sgRNAs_fasta <- opt$output_sgRNAs_fasta
output_sites_fasta <- opt$output_sites_fasta
detect_offtarget_genes_full <- opt$detect_offtarget_genes_full
keep_TINDRi_input_sgRNAs <- opt$keep_TINDRi_input_sgRNAs
keep_TINDRi_input_sites <- opt$keep_TINDRi_input_sites
keep_TINDRi_matches <- opt$keep_TINDRi_matches
path_python <- opt$path_python
feature_type <- opt$feature_type


#### 2. Preliminaries ####
# check inputs
if(is.null(input_genome) | is.null(sgRNA_file)){
  stop("Please specify both input_genome and sgRNA_file")
}
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
if(!output_exact & !output_full){
  stop("Neither exact nor full binding site list for given sgRNAs specified as output. Please set at least one of output_exact or output_full to TRUE to retrieve output.")
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
#message(paste0(Sys.time(), ": Loading packages (and installing if needed)..."))
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
    return(chrom_seq)
  }))
  # make sure resulting gene tables can be merged by selecting same columns
  genes_cols <- c("seqnames", "start", "end", "strand", "locus_tag")
  genes <- do.call(rbind, lapply(seq.int(length(genome_gb_ls)), function(chrom){
    # find features on each contig
    if(isEmpty(genes(genome_gb_ls[[chrom]]))){
      if(isEmpty(transcripts(genome_gb_ls[[chrom]])) & isEmpty(otherFeatures(genome_gb_ls[[chrom]]))){
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
}
genes <- genes[, c("seqid", "strand", "start", "end", "locus_tag", "gene")]
genomeID <- switch(input_type, 
                   accessionnr = input_genome, 
                   gbfile = sub("\\..*", "", sub(".*/", "", input_genome)))
rm(genome_gb_ls, genome_txt)


#### 4. Read and prepare sgRNAs ####
sgRNAs <- read.csv(sgRNA_file, header = FALSE)
colnames(sgRNAs) <- c("tag", "sgRNA")
## keep only sgRNAs of correct length, also removing comments as e.g. "No PAM"
sgRNAs <- sgRNAs[unlist(lapply(sgRNAs$sgRNA, nchar)) == sum(regions), ]
## check if no duplicates in either tag or sgRNA
if(any(duplicated(sgRNAs$sgRNA))){
  warning(paste("Duplicate sgRNA sequences:", paste(sgRNAs$tag[duplicated(sgRNAs$sgRNA) | duplicated(sgRNAs$sgRNA, fromLast = TRUE)], collapse = ", ")))
}
if(any(duplicated(sgRNAs$tag))){stop("Duplicate sgRNA names, please ensure all names are unique")}
## attach PAM and make sure all bases are upper case
sgRNAs$sgRNA <- apply(sgRNAs, 1, function(x){paste0(toupper(substr(x[2], 1, sum(regions))), PAM)})
# format for CRISPRseek::searchHits
gRNAs <- DNAStringSet(unlist(sgRNAs[, 2]), use.names = FALSE)
names(gRNAs) <- unlist(sgRNAs[, 1])
rm(sgRNAs)
# output if wanted
if(output_sgRNAs_fasta){writeXStringSet(gRNAs, paste0(outdir, "/", genomeID, "_sgRNAs.fasta"))}
# write
write.table(cbind(1:length(gRNAs), 
                  as.character(reverse(narrow(gRNAs, 1, sum(regions))))), 
            paste0(outdir, "/sg_candidates_within_genes.csv"), 
            sep = ",", 
            col.names = FALSE, 
            row.names = FALSE, quote = FALSE)


#### 5. Find all possible binding sites ####
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
# write
write.table(cbind(1:length(all_sites), 
                  as.character(reverse(narrow(all_sites, 1, sum(regions))))), 
            paste0(outdir, "/allsgcandidates.csv"), 
            sep = ",", 
            col.names = FALSE, 
            row.names = FALSE, quote = FALSE)


#### 6. Match candidate sgRNAs to binding sites ####
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
                        paste0(outdir, "/", genomeID, "_sgRNAs.csv")))
}
if(!keep_TINDRi_input_sites){
  invisible(file.remove(paste0(outdir, "/allsgcandidates.csv")))
} else{
  invisible(file.rename(paste0(outdir, "/allsgcandidates.csv"), 
                        paste0(outdir, "/", genomeID, "_all_binding_sites.csv")))
}


#### 7. Compute scores needed in any case ####
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
    # leave penalty function customizable
    do.call(pen_func, list(rev(penalties)[as.integer(unlist(strsplit(mm, ","))) == 1]))
  }))
  stopCluster(cl)
  if(output_full){
    write(reprAct, file = paste0(outdir, "/tmp_reprAct.tsv"), sep = "\t", ncolumns = 1)
  }
  # nmismatch
  cl <- makeCluster(n_cores)
  nmm <- unlist(parLapply(cl, candidate_mm, function(mm){
    sum(as.integer(unlist(strsplit(mm, ","))))
  }))
  stopCluster(cl)
  #candidate_ID <- which(nmm == 0)
  if(output_full){
    write(nmm, file = paste0(outdir, "/tmp_nmm.tsv"), sep = "\t", ncolumns = 1)
  }
} else{
  reprAct <- unlist(mclapply(candidate_mm, function(mm){
    # leave penalty function customizable
    do.call(pen_func, list(rev(penalties)[as.integer(unlist(strsplit(mm, ","))) == 1]))
  }, mc.cores = n_cores))
  if(output_full){
    write(reprAct, file = paste0(outdir, "/tmp_reprAct.tsv"), sep = "\t", ncolumns = 1)
  }
  # nmismatch
  nmm <- unlist(mclapply(candidate_mm, function(mm){
    sum(as.integer(unlist(strsplit(mm, ","))))
  }, mc.cores = n_cores))
  #candidate_ID <- which(nmm == 0)
  if(output_full){
    write(nmm, file = paste0(outdir, "/tmp_nmm.tsv"), sep = "\t", ncolumns = 1)
  }
}
rm(candidate_mm)
# required for either output
# CRISPRseek offset in within-target bp locus of sgRNA names
offset_emp <- 3


#### 8. Get zero-mismatch sites per sgRNA plus scores ####
if(output_exact){
  # all information comes from all_sites ID, except reprAct
  sg_hits <- read.table(paste0(outdir, "/missmatch_matrix_end.txt"), 
                        colClasses = c("integer", "NULL", "integer"), 
                        sep = "\t", 
                        header = FALSE)
  # get 1st- and 2nd-best target IDs
  #names(reprAct) <- sg_hits[, 2]
  if(platform == "windows"){
    cl <- makeCluster(n_cores)
    maxreprAct <- parLapply(cl, split(data.frame(sg_hits[, 2], reprAct, nmm), sg_hits[, 1]), function(sg_ra){
      sg_max1 <- sg_ra[, 2] == max(sg_ra[, 2])
      sg_max2 <- sg_ra[, 2] == max(sg_ra[!sg_max1, 2])
      list("sites_1st" = sg_ra[sg_max1, 1], 
           "reprAct_1st" = sg_ra[sg_max1, 2][1], 
           "sites_2nd" = sg_ra[sg_max2, 1], 
           "reprAct_2nd" = sg_ra[sg_max2, 2][1], 
           "nmm_1st" = sg_ra[sg_max1, 3], 
           "nmm_2nd" = sg_ra[sg_max2, 3])
    })
    stopCluster(cl)
  } else{
    maxreprAct <- mclapply(split(data.frame(sg_hits[, 2], reprAct, nmm), sg_hits[, 1]), function(sg_ra){
      sg_max1 <- sg_ra[, 2] == max(sg_ra[, 2])
      sg_max2 <- sg_ra[, 2] == max(sg_ra[!sg_max1, 2])
      list("sites_1st" = sg_ra[sg_max1, 1], 
           "reprAct_1st" = sg_ra[sg_max1, 2][1], 
           "sites_2nd" = sg_ra[sg_max2, 1], 
           "reprAct_2nd" = sg_ra[sg_max2, 2][1], 
           "nmm_1st" = sg_ra[sg_max1, 3], 
           "nmm_2nd" = sg_ra[sg_max2, 3])
    }, mc.cores = n_cores)
  }
  rm(reprAct, nmm)
  # compute binding site scores 1st-best hits
  if(platform == "windows"){
    cl <- makeCluster(n_cores)
    clusterExport(cl, c("genes", "regions", "offset_emp"))
    NTgene_ls <- parLapply(cl, names(all_sites[unlist(lapply(maxreprAct, function(x){x[[1]]}))]), function(site){
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
    hitinfo1 <- do.call(rbind, parLapply(cl, 
                                        split(NTgene_ls, 
                                              rep(seq.int(maxreprAct), 
                                                  unlist(lapply(maxreprAct, function(x){length(x[[1]])})))), 
                                        function(hits){
      # cannot paste with sep or col = ","...
      gsub(" ", ",", do.call(paste, hits))
    }))
    stopCluster(cl)
  } else{
    NTgene_ls <- mclapply(names(all_sites[unlist(lapply(maxreprAct, function(x){x[[1]]}))]), function(site){
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
    hitinfo1 <- do.call(rbind, mclapply(split(NTgene_ls, 
                                             rep(seq.int(maxreprAct), 
                                                 unlist(lapply(maxreprAct, function(x){length(x[[1]])})))), 
                                       function(hits){
      # cannot paste with sep or coll = ","...
      gsub(" ", ",", do.call(paste, hits))
    }, mc.cores = n_cores))
  }
  colnames(hitinfo1) <- c("targets_1st", "target_names_1st", "chroms_1st", "strands_1st", "ranges_1st")
  # compute binding site scores 2nd-best hits
  site_index <- unlist(lapply(maxreprAct, function(x){
    if(isEmpty(x[[3]])){
      NA
    } else{
      names(all_sites[x[[3]]])
    }
  }))
  if(platform == "windows"){
    cl <- makeCluster(n_cores)
    clusterExport(cl, c("genes", "regions", "offset_emp"))
    NTgene_ls <- parLapply(cl, site_index, function(site){
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
    hitinfo2 <- do.call(rbind, parLapply(cl, 
                                         split(NTgene_ls, 
                                               rep(seq.int(maxreprAct), 
                                                   unlist(lapply(maxreprAct, function(x){ifelse(length(x[[3]]) == 0, 1, length(x[[3]]))})))), 
                                         function(hits){
                                           # cannot paste with sep or col = ","...
                                           gsub(" ", ",", do.call(paste, hits))
                                         }))
    stopCluster(cl)
  } else{
    NTgene_ls <- mclapply(site_index, function(site){
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
    hitinfo2 <- do.call(rbind, mclapply(split(NTgene_ls, 
                                              rep(seq.int(maxreprAct), 
                                                  unlist(lapply(maxreprAct, function(x){ifelse(length(x[[3]]) == 0, 1, length(x[[3]]))})))), 
                                        function(hits){
                                          # cannot paste with sep or coll = ","...
                                          gsub(" ", ",", do.call(paste, hits))
                                        }, mc.cores = n_cores))
  }
  colnames(hitinfo2) <- c("targets_2nd", "target_names_2nd", "chroms_2nd", "strands_2nd", "ranges_2nd")
  # get site sequences
  seqs1 <- all_sites[unlist(lapply(maxreprAct, function(x){x[[1]]}))]
  seqs2 <- unlist(DNAStringSetList(lapply(maxreprAct, function(x){
    if(isEmpty(x[[3]])){
      DNAStringSet(paste(rep("N", sum(regions) + nchar(PAM)), collapse = ""))
    } else{
      all_sites[x[[3]]]
    }
  })))
  seqs1 <- do.call(rbind, lapply(split(as.character(seqs1), 
                                       rep(seq.int(maxreprAct), 
                                           unlist(lapply(maxreprAct, function(x){length(x[[1]])})))), 
                                 function(sg_site){
                                   c(paste(substr(sg_site, 1, sum(regions)), collapse = ","), 
                                     paste(substr(sg_site, sum(regions) + 1, sum(regions) + 1 + nchar(PAM)), collapse = ","))
                                 }))
  colnames(seqs1) <- c("seq_1st", "PAM_1st")
  seqs2 <- do.call(rbind, lapply(split(as.character(seqs2), 
                                          rep(seq.int(maxreprAct), 
                                              unlist(lapply(maxreprAct, function(x){ifelse(length(x[[3]]) == 0, 1, length(x[[3]]))})))), 
                                    function(sg_site){
                                      c(paste(substr(sg_site, 1, sum(regions)), collapse = ","), 
                                        paste(substr(sg_site, sum(regions) + 1, sum(regions) + 1 + nchar(PAM)), collapse = ","))
                                    }))
  colnames(seqs2) <- c("seq_2nd", "PAM_2nd")
  # bad seed flag
  # cannot be multiple per sgRNA, use that
  bad_seed_found <- lapply(bad_seeds, function(seed){
    which(narrow(gRNAs[sort(unique(sg_hits[, 1]))], 
                 sum(regions) - nchar(seed) + 1, 
                 sum(regions)) %in% seed)
  })
  sgRNA_bad_seed <- rep(NA, length(unique(sg_hits[, 1])))
  for(i in seq.int(bad_seeds)){
    sgRNA_bad_seed[bad_seed_found[[i]]] <- bad_seeds[i]
  }
  # cut_site flag
  spacers_ext <- xscat(oligoForwardOverhang, 
                       narrow(gRNAs[sort(unique(sg_hits[, 1]))], 
                              1, sum(regions)), 
                       reverseComplement(DNAString(oligoReverseOverhang)))
  cut_sites_plusrev <- c(reverseComplement(DNAStringSet(cut_sites)), cut_sites)
  names(cut_sites_plusrev) <- c(paste0("reverse_", names(cut_sites)), names(cut_sites))
  cut_site_found <- vwhichPDict(cut_sites_plusrev, spacers_ext)
  sgRNA_cut_site <- unlist(lapply(cut_site_found, function(sg_cut){
    ifelse(isEmpty(sg_cut), NA, paste(names(cut_sites_plusrev[sg_cut]), collapse = ","))
  }))
  # make data frame (ordered by sorted unique(candidate_sgID))
  all_df <- data.frame("sgRNA_name" = names(gRNAs[sort(unique(sg_hits[, 1]))]), 
                       "sgRNA_seq" = substr(gRNAs[sort(unique(sg_hits[, 1]))], 1, sum(regions)), 
                       "sgRNA_GC" = as.numeric(letterFrequency(narrow(gRNAs[sort(unique(sg_hits[, 1]))], 1, sum(regions)), letters = "GC") / sum(regions)), 
                       "bad_seed" = sgRNA_bad_seed, 
                       "cut_site" = sgRNA_cut_site, 
                       "reprAct_1st" = unlist(lapply(maxreprAct, function(x){x[2]})), 
                       "nmismatch_1st" = unlist(lapply(maxreprAct, function(x){paste(x[[5]], collapse = ",")})), 
                       hitinfo1, 
                       seqs1, 
                       "reprAct_2nd" = unlist(lapply(maxreprAct, function(x){x[4]})), 
                       "nmismatch_2nd" = unlist(lapply(maxreprAct, function(x){paste(x[[6]], collapse = ",")})), 
                       hitinfo2, 
                       seqs2)
  # save
  write.csv(all_df, paste0(outdir, "/", genomeID, "_sgRNA_hits_summary.csv"), 
            row.names = FALSE)
  rm(all_df, hitinfo1, hitinfo2, NTgene_ls, seqs1, seqs2, maxreprAct, sg_hits)
}


#### 9. Get full table plus scores ####
if(output_full){
  # reading whole file 
  full_df <- read.table(paste0(outdir, "/missmatch_matrix_end.txt"), 
                        header = FALSE, sep = "\t")
  colnames(full_df) <- c("sgRNA_ID", "mismatches", "site_ID")
  # add sgRNA name
  full_df$sgRNA_name <- names(gRNAs[full_df$sgRNA_ID])
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
  full_df$sgRNA_seq <- substr(gRNAs[full_df$sgRNA_ID], 1, sum(regions))
  full_df$site_seq <- substr(all_sites[full_df$site_ID], 1, sum(regions))
  full_df$site_PAM <- substr(all_sites[full_df$site_ID], sum(regions) + 1, sum(regions) + nchar(PAM))
  full_df$sgRNA_GC <- as.numeric(letterFrequency(gRNAs[full_df$sgRNA_ID], letters = "GC") / sum(regions))
  # need lots of RAM for detecting non-template strand gene overlaps of sites
  if(detect_offtarget_genes_full){
    if(platform == "windows"){
      cl <- makeCluster(n_cores)
      clusterExport(cl, "genes")
      NTgene <- do.call(rbind, parLapply(cl, split(full_df[, c("chrom", "strand", "start", "end")], 
                                                   seq.int(nrow(full_df))), function(site){
        tmp_i <- as.integer(site["start"]) - genes$end <= 0 & 
          as.integer(site["end"]) - genes$start >= 0 & 
          as.character(site["strand"]) != genes$strand & 
          unlist(lapply(gsub(" ", "", genes$seqid), grepl, x = site["chrom"]))
        # find gene hits on site, if any (can be 0, or >1 if overlap)
        if(any(tmp_i)){
          c(paste(genes$locus_tag[which(tmp_i)], collapse = ","), 
            paste(ifelse(is.na(genes$gene[which(tmp_i)]), genes$locus_tag[which(tmp_i)], genes$gene[which(tmp_i)]), collapse = ","), 
            paste(site["chrom"], collapse = ","), 
            paste(site["strand"], collapse = ","), 
            mapply(paste, site["start"], site["end"], MoreArgs = list(sep = "..", collapse = ",")))
        } else{
          c(NA, 
            NA, 
            paste(site["chrom"], collapse = ","), 
            paste(site["strand"], collapse = ","), 
            mapply(paste, site["start"], site["end"], MoreArgs = list(sep = "..", collapse = ",")))
        }
      }))
      stopCluster(cl)
    } else{
      NTgene <- do.call(rbind, mclapply(split(full_df[, c("chrom", "strand", "start", "end")], seq.int(nrow(full_df))), function(site){
        tmp_i <- as.integer(site["start"]) - genes$end <= 0 & 
          as.integer(site["end"]) - genes$start >= 0 & 
          as.character(site["strand"]) != genes$strand & 
          unlist(lapply(gsub(" ", "", genes$seqid), grepl, x = site["chrom"]))
        # find gene hits on site, if any (can be 0, or >1 if overlap)
        if(any(tmp_i)){
          c(paste(genes$locus_tag[which(tmp_i)], collapse = ","), 
            paste(ifelse(is.na(genes$gene[which(tmp_i)]), genes$locus_tag[which(tmp_i)], genes$gene[which(tmp_i)]), collapse = ","), 
            paste(site["chrom"], collapse = ","), 
            paste(site["strand"], collapse = ","), 
            mapply(paste, site["start"], site["end"], MoreArgs = list(sep = "..", collapse = ",")))
        } else{
          c(NA, 
            NA, 
            paste(site["chrom"], collapse = ","), 
            paste(site["strand"], collapse = ","), 
            mapply(paste, site["start"], site["end"], MoreArgs = list(sep = "..", collapse = ",")))
        }
      }, mc.cores = n_cores))
    }
    colnames(NTgene) <- c("all_targets", "all_target_names", "all_chroms", "all_strands", "all_ranges")
    full_df <- data.frame(full_df, NTgene)
  }
  # write
  write.csv(full_df, 
            file = paste0(outdir, "/", genomeID, "_sgRNA_hits_full.csv"), 
            row.names = FALSE)
  # can remove tmp files here
  file.remove(paste0(outdir, "/tmp_nmm.tsv"))
  file.remove(paste0(outdir, "/tmp_reprAct.tsv"))
}

# remove file as not needed anymore
if(!keep_TINDRi_matches){
  file.remove(paste0(outdir, "/missmatch_matrix_end.txt"))
} else{
  file.rename(paste0(outdir, "/missmatch_matrix_end.txt"), 
              paste0(outdir, "/", genomeID, "_missmatch_matrix_end.txt"))
}