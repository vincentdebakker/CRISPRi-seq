# Genome-wide sgRNA library design and evaluation for CRISPRi-seq
Author: Vincent de Bakker  
Python code by Afonso Bravo    
See also: <https://www.veeninglab.com/crispri-seq>  
Published as a part of: 
* Liu, X., Kimmey, J.M., Matarazzo, L., de Bakker, V., Van Maele, L., Sirard, J.-C., Nizet, V., Veening, J.-W. (2021) Exploration of Bacterial Bottlenecks and _Streptococcus pneumoniae_ pathogenesis by CRISPRi-seq. _Cell Host Microbe_ **29**, 107-120.e6. https://doi.org/10.1016/j.chom.2020.10.001 (https://www.sciencedirect.com/science/article/pii/S193131282030559X).
* \<submitted\>

This repository contains R pipelines to (1) design a _de novo_  CRISPR interference (CRISPRi) single-guide RNA (sgRNA) library and (2) evaluate all potential binding sites of a given sgRNA library on any genome, either provided as GenBank file or NCBI assembly accession number. They were in principle developed for application to prokayotic (specifically bacterial) genomes. For both pipelines there is a version that can be run from within an R environment (we recommend RStudio), and a version that can be called from the command line using the Rscript command. Additionally, the repository contains an example of possible downstream expected efficiency analyses for the latter, as published in Liu _et al._ (2021). Lastly, it contains a subdirectory with the code to generate the naïve power analysis as published in \<submitted manuscript\>.

# General information
The pipelines automatically use multi-core processing, using all available logical processors except one. Input genomes can be either local GenBank files or NCBI assembly accession numbers and multi-chromosome/contig/plasmid genomes are fully supported. sgRNA candidate sequences and binding sites are detected using the [`CRISPRseek` package](https://bioconductor.org/packages/release/bioc/html/CRISPRseek.html) and match/mismatch vectors between them are defined by the here presented Tool for Identification of Nucleotide-wise, Dual-sequence, Region-delimited identity (**TINDRi.py**), which is automatically called by the pipelines. Repression activity of an sgRNA at a given site is estimated as the product of within-spacer position-specific penalties of mismatched nucleotides, and referred to as the `reprAct` score. By using the product, we assume that the observation of Qi _et al._ (2013) that individual nucleotide mismatch effects on repression activity were independent, extends to the case of >2 and non-adjacent nucleotide mismatches. Although penalties and the `reprAct` function can be customized by the user, built-in options from the literature are included with penalties from Qi _et al._ (2013) \[_E. coli_\] and Hawkins _et al._ (2020) \[_E.coli_ and _B. subtilis_\].  

The pipelines should work on all operating systems, although they have only been tested on Windows 10 and Ubuntu 18.04.5 LTS. 

# Dependencies
It is necessary to have the following software installed:
* [R version 4.x.x](https://www.r-project.org/); if the user does not plan to work from the command line, we recommend to install and run from within the [RStudio](https://www.rstudio.com/) integrated development environment.
* Python version 3.x; we advise to install the complete [Anaconda suite](https://www.anaconda.com), as this should ensure a proper working python environment.

Installation of packages and modules should be performed automatically the first time a script is invoked (so the script will then also take longer to run that one time). If errors arise nonetheless, please ensure all following dependencies / software packages and modules are installed, manually:
* R packages [CRISPRseek](https://bioconductor.org/packages/release/bioc/html/CRISPRseek.html), [Biostrings](https://doi.org/doi:10.18129/B9.bioc.Biostrings), if using a GenBank file as genome input [genbankr](https://bioconductor.org/packages/release/bioc/html/genbankr.html), if using an NCBI accession assembly number as genome input [biomartr](https://cran.r-project.org/package=biomartr), and if using the *_cmd.R scripts also [optparse](https://CRAN.R-project.org/package=optparse).
* Python modules _numpy_, _numba_, _time_, _multiprocessing_, _datetime_, _os_, _sys_.

# _de novo_ sgRNA design
The pipeline in the files **sgRNA_library_design.R** (invoked within R) and **sgRNA_library_design_cmd.R** (invoked from command line) can be run to automatically design an sgRNA library for all annotated features in a given genome.  

It identifies all possible candidate sgRNAs targeting the non-template strand for all features, as well as all possible sgRNA binding sites on either strand on the whole genome (all contigs/plasmids/chromosomes). Subsequently, it calls the TINDRi.py script, which efficiently matches all pair-wise combinations of sgRNA candidates and binding sites. It detects all relevant sites for each sgRNA based on a user-defined maximum cumulative allowed number of nucleotide mismatches per sub-spacer region, and generates a matrix scoring the mismatches for each relevant sgRNA-site combination. Using that data, the maximum mismatch-based off-target `reprAct` estimate is computed for each sgRNA, per feature. Finally, the pipeline first selects for each feature the candidate sgRNAs with the lowest expected maximum off-target repression activity (within some user-defined error range of the minimum) and secondly, of that subset, the one sgRNA with the smallest distance to the start codon. Provided that multiple candidate sgRNAs are available, the user may request to design any desired number of sgRNAs per feature.  

The pipeline can also exclude or avoid design of sgRNAs with (1) user-defined bad seed sequences, (2) certain motifs in the spacer + primer overhang (in given order of importance, since it might be more crucial to avoid certain restriction enzyme cut sites than e.g. tandem TTTTT motifs) or (3) non-standard bases (i.e., non-ATGC, as for instance some genome sequences contain N's). If desired, the pipeline can output the full list of candidate sgRNAs for all features plus their scores and additional details as well, so that the user may make their own selection. In addition, corresponding forward and reverse oligo's are designed for each selected sgRNA for library construction.  

## Command line version

### Usage
Run:
```
Rscript sgRNA_library_design_cmd.R --help
```
to get the user manual and an outline of all flag parameters than can be used. They are all either the same as, or slight adaptations of, the input parameters of **sgRNA_library_design.R** as described below in more detail. Only the first three input parameters are absolutely required. In its simplest form, the pipeline can thus be run as such: 
```
Rscript sgRNA_library_design_cmd.R -g genome.gb -o ~/output/ -t ~/wd/
```
where `genome.gb` is a GenBank file of the genome for which the library is to be designed, `~/output/` is the local directory to write the output files to and `~/wd/` is the local working directory in which **TINDRi.py** can be found.  

Depending on the genome size, number of annotated features, number of allowed mismatches in binding site identification and computer specifications (such as number of available logical processors), computation times will likely be in the order of 10-60 minutes. 

### Output
The pipeline can return multiple files, depending on user input parameters: 
* A .csv file with (an) optimal sgRNA(s) per annotated feature.
* A .csv file with all candidate sgRNAs for all annotated features.
* A .csv file with all annotated target features for which no sgRNA could be designed (due to a lack of PAM sequences on the non-template strand).
* A .txt file with statistics on the run and used input parameters.
* A .fasta file with all annotated target feature sequences for which sgRNAs were designed.
* A .fasta file with all candidate sgRNA spacer sequences. 
* A .fasta file with all sgRNA binding site sequences on the genome. 
* A .csv file with all candidate sgRNA - binding site combinations and scores. 
* A .csv file with all candidate sgRNA IDs and sequences, used as input for **TINDRi.py**.
* A .csv file with all binding site IDs and sequences, used as input for **TINDRi.py**.
* A .csv file with the full mismatch matrix for all candidate sgRNA - binding site combinations, as generated by **TINDRi.py**.
Users will in general be after either or both of the first two files. The third and fourth file are produced by default.

#### Meaning of variables in sgRNA design output files
Variable | Description
--- | ---
`target` | Locus_tag key of the feature the sgRNA was designed to target.
`target_name` | Name of the feature the sgRNA was designed to target.
`sgRNA_name` | Unique name of the sgRNA.
`sgRNA_seq` | Spacer sequence of the sgRNA (PAM-distal to -proximal).
`sgRNA_PAM` | PAM sequence of the on-target binding site of the sgRNA.
`sgRNA_GC` | GC content of the sgRNA spacer, as ratio.
`targetlength_bp` | Length of the target feature in base pairs.
`dist2SC_bp` | Distance of the on-target sgRNA binding site to the target feature start codon in base pairs.
`dist2SC_rel` | Relative distance of the on-target sgRNA binding site to the target feature start codon on the \[0,1\] interval.
`maxOffreprAct` | Maximum estimated repression activity detected, excluding the on-target binding site the sgRNA was designed for.
`all_targets` | Locus_tags keys of all target features on whose non-template strand the sgRNA has a perfect, zero-mismatch binding site.
`all_target_names` | Names of all target features on whose non-template strand the sgRNA has a perfect, zero-mismatch binding site.
`all_chroms` | Chromosome/plasmid/contig ID of the zero-mismatch binding sites.
`all_strands` | Strand of the zero-mismatch binding sites.
`all_ranges` | Chromosomal locations of the zero-mismatch binding sites, counting from the chromosome start.
`oligoForward` | Forward oligo (primer) to be ordered to construct the annealed sgRNA library pool.
`oligoReverse` | Reverse oligo (primer) to be ordered to construct the annealed sgRNA library pool.
`cut_site` | If any, name of input `cut_site` sequences detected in the sgRNA spacer; otherwise NA.
`nonstandard_nt` | If any, non-standard nucleotide (non-ATGC) detected in sgRNA spacer; otherwise NA.
`bad_seed` | If any, name of input `bad_seed` sequences detected in the sgRNA spacer; otherwise NA.
`maxOff_mismatches` | Comma-separated vector of same length as the sgRNA spacer (PAM-distal to -proximal), indicating for each nucleotide whether (1) or not (0) there was a mismatch between that nucleotide in the spacer and the corresponding nucleotide in the off-target binding site with the highest estimated maximum off-target repression activity (as given by `maxOffreprAct`).
`maxOff_seq` | Sequence (PAM-distal to -proximal) of the off-target binding site with the highest estimated maximum off-target repression activity (as given by `maxOffreprAct`).

## Within-R version

### Usage
Open the script **sgRNA_library_design.R** (we recommend in the RStudio integrated development environment), adjust the input parameters in the first code section `#### 1. Settings ####` and then run the whole script with those input parameters. In its simplest form, only `input_genome`, `outdir` and `TINDRidir` need to be specified (under `## REQUIRED ##`). Depending on the genome size, number of annotated features, number of allowed mismatches in binding site identification and computer specifications (such as number of available logical processors), computation times will likely be in the order of 10-60 minutes. 

### Input
Parameter | Description
--- | ---
`input_genome` | **required**; character. Either the full path to a GenBank file (extension .gb, .gbk, .gbf, or .gbff), or an NCBI assembly accession number (starting with either GCA_ \[GenBank\] or GCF_ \[RefSeq\], e.g. [GCA_003003495.1](https://www.ncbi.nlm.nih.gov/assembly/GCA_003003495.1/)).
`outdir` | **required** character. Path to desired output directory
`TINDRidir` | **required** character. Path to directory in which **TINDRi.py** is stored.
`path_ncbi_downloads` | character. Path to directory where downloaded files from NCBI should be stored. Only **required if** `input_genome` is a GenBank file. Default: NA.
`n_sgRNA` | integer. Desired number of sgRNAs to be designed per annotated feature. Default: 1.
`regions` | integer vector. Indicates the sub-spacer regions for which to define max_mismatch_cum. The sum of these integers must equal the sgRNA spacer length. Default: c(7, 2, 11). 
`max_mismatch_cum` | integer vector. Indicates the maximum cumulative number of mismatches allowed over regions for a binding site to be considered relevant. Should be of the same length as `regions`. For example, max_mismatch_cum of c(1, 2, 11) combined with regions of c(7, 2, 11) allows for 1 mismatch in the first 7 nucleotides (PAM-proximal to -distal), 2 mismatches in the first 7+2=9 nucleotides and 11 mismatches in the first 7+2+11=20 nucleotides (which is the whole sgRNA spacer). Binding sites that do not fulfill these conditions for a given sgRNA are not regarded as relevant and therefore excluded (in the example, a site with e.g. 2 mismatches in the first 7 nucleotides with the spacer would not be considered relevant for that sgRNA, as its repression activity is likely too low at that site). Default: c(1, 2, 11).
`reprAct_penalties` | character or numeric vector. Indicates within-spacer position-dependent penalties for `reprAct` estimation. For using built-in penalties, must be one of "HawkinsEcoMean", "HawkinsBsuMean", "HawkinsEcoMedian", "HawkinsBsuMedian", "Qi", "QiMean" or "QiMedian". These vectors are the mean, median or raw 20-nt spacer relative repression activity values found by Hawkins et al. (2020) (means and medians per nucleotide) and Qi et al. (2013) (means and medians per region as defined in their paper) for _E. coli_ (Eco) and _B. subtilis_ (Bsu). For custom penalties, please specify a numeric vector (PAM-proximal to -distal) of length `sum(regions)` (i.e., sgRNA spacer length); if `pen_func` is `prod` (recommended), the integers should be on the \[0,1\] interval. Default: "HawkinsBsuMedian" (PAM-proximal to -distal: 0.35, 0.46, 0.5, 0.53, 0.57, 0.53, 0.57, 0.67, 0.83, 0.9, 0.91, 0.94, 0.95, 0.92, 0.95, 0.84, 0.94, 0.94, 0.96, 0.99).
`pen_func` | character or function. Defines how to process reprAct_penalties of sgRNA-site mismatches to produce `reprAct` estimates; strong advice not to change. Default: prod (product).
`errorRange_maxOffreprAct` | numeric. The error allowed on maximum off-target repression activity estimates. Candidate sgRNAs with an estimated repression activity within this range from the minimum found off-target repression activity all pass the specificity filter in selecting optimal sgRNAs per feature. This is basically the parameter that weighs the importance of selecting sgRNAs based on estimated off-target activity versus distance of the binding site to the start codon. Setting it to 0 means maximum importance is placed on minimizing off-target activity, so no sgRNA selection based on distance to start codon. Default: 0.2 (20%, as this reproduced >80% of our manually designed sgRNA library of which we know it works well).
`avoidNonBaseNT` | logical. Whether or not to avoid designing sgRNAs including non-standard bases (so any other than ATGC). Relevant, as genome sequences sometimes include N's, which the user might want to avoid. sgRNAs with non-standard bases are also flagged in the output. Default: TRUE. 
`bad_seeds` | character vector. Specifies bad seed sequences to exclude, avoid or flag. Default: c("ACCCA", "TGGAA"), based on Cui et al. (2018). 
`bad_seed_rule` | character. Indicates how to treat sgRNAs with bad seed sequences given by `bad_seeds`. One of "exclude", "avoid" or "ignore", which respectively exclude candidate sgRNAs with these sequences from being designed altogether, avoid their selection if any other candidate sgRNA is available, or just flag sgRNAs with bad seeds in the output. Default: "ignore".
`cut_sites` | named character vector. Indicates in decreasing order of importance motifs to detect in sgRNA spacer + primer sequences (`oligoForwardOverhang` & `oligoReverseOverhang`). Could be e.g. specific restriction enzyme cut sites, or tandem T nucleotides to be avoided in design. Reverse-complement sequences are automatically detected too. Default: c(BsmBI = "CGTCTC", tandemT = "TTTTT"). 
`cut_site_rule` | character. Indicates how to process sgRNAs with specific motifs indicated by `cut_sites`. One of "exclude", "avoid" or "ignore", which respectively exclude, avoid or just flag sgRNAs with these sequences during design. If "avoid", design of sgRNAs with these motifs is avoided in the order of inputted `cut_sites` in the presence of alternatives. For example, if for a given feature only sgRNA candidates with BsmBI cut sites and tandem TTTTT motifs exist, the pipeline will prefer to select candidates with the tandem Ts, as per order of default `cut_sites` input. Default: "avoid". 
`oligoForwardOverhang` | character. Post-annealing 5'-3' overhang for the forward oligo to be compatible with the digested sgRNA backbone vector. Default: `TATA` (BsmBI).
`oligoReverseOverhang` | character. Post-annealing 5'-3' overhang for the reverse oligo to be compatible with the digested sgRNA backbone vector. Default: `AAAC` (BsmBI).
`PAM` | character. Proto-spacer Adjacent Motif sequence. Default: "NGG" (_S. pyogenes_ dCas9).
`filter_out_duplicates` | logical. Whether or not to filter out duplicate sgRNA spacer sequences in the output list of optimized sgRNAs. Equal sgRNAs can be designed for different annotated features, for instance in the case of gene copies. If set to FALSE, Each designed sgRNA will get its own row in the output, thus appearing multiple times if the same spacer sequence was designed multiple times. Default: TRUE.
`output_optimized_list` | logical. Whether or not to output the .csv file with the list of designed optimal sgRNAs per annotated feature. Note that setting both this parameter and `output_all_candidates` to FALSE is nonsensical. Default: TRUE.
`output_all_candidates` | logical. Whether or not to output the .csv file with the full list of candidate sgRNAs for all annotated features. Useful if the user wants to explore alternative sgRNAs to be used or to make their own selection from the full list. Default: TRUE. 
`output_target_fasta` | logical. Whether or not to output a .fasta file containing the sequences of all found annotated genetic features for which sgRNAs are designed. Default: FALSE. 
`output_sgRNAs_fasta` | logical. Whether or not to output a .fasta file containing the sequences of all candidate identified sgRNAs. Default: FALSE. 
`output_sites_fasta` | logical. Whether or not to output a .fasta file containing the sequences of all identified sgRNA binding sites on the genome. Default: FALSE. 
`output_full_list` | logical. Whether or not to output the complete list of candidate sgRNA - binding site combinations. Note that if set to TRUE, this may take considerable computation time and resources and generate a large .csv file. Not required, recommended to leave set to FALSE. Default: FALSE.
`detect_offtarget_genes_full` | logical. If `output_full_list` is set to TRUE, whether or not to also detect for each sgRNA-site combination if the site falls within an annotated feature on its non-template strand. Note that if set to TRUE, will add again much more computation time. Not required, recommended to leave set to FALSE. Default: FALSE.
`keep_TINDRi_input_sgRNAs` | logical. Whether or not to keep the intermittently generated .csv file containing all candidate sgRNA IDs and sequences (PAM-proximal to -distal), used as input for **TINDRi.py**. Default: FALSE.
`keep_TINDRi_input_sites` | logical. Whether or not to keep the intermittently generated .csv file containing all binding site IDs and sequences (PAM-proximal to -distal), used as input for **TINDRi.py**. Default: FALSE.
`keep_TINDRi_matches` | logical. Whether or not to keep the intermittently generated .csv file containing the full TINDRi.py mismatch output matrix. Note this file is generally large. Default: FALSE. 
`path_python` | character. String indicating full path to Python3 installation (e.g. ".../anaconda3/python.exe" on Windows systems). Specify only if the Python3 installation is not automatically detected by the pipeline (it will issue an error message indicating this). Default: NULL (automatic detection). 
`feature_type` | character. Feature key tag for each of whose unique identifiers the sgRNAs should be designed. Experimental pipeline feature, in principle do not alter. Default: `"locus_tag"`.

### Output
Same as described above for the command line version of the pipeline.


# sgRNA binding site evaluation
The pipeline in the files **sgRNA_library_evaluation.R** (invoked within R) and **sgRNA_library_evaluation_cmd.R** (invoked on command line) identifies and scores all binding sites up to a specified cumulative region-delimited number of mismatches in a given genome for each sgRNA of a given library.  

Much like in the _de novo_ design pipeline, **TINDRi.py** is called to compare relevant input sgRNA sequences and potential genome-wide binding sites. These results are processed into a repression activity estimate, which is reported alongside other scores.

## Command line version

### Usage
Run:
```
Rscript sgRNA_library_evaluation_cmd.R --help
```
to get the user manual and an outline of all flag parameters than can be used. They are all either the same as, or slight adaptations of, the input parameters of **sgRNA_library_evaluation.R** as described below in more detail. Only the first four input parameters are absolutely required. In its simplest form, the pipeline can thus be run as such:
```
Rscript sgRNA_library_evaluation_cmd.R -g genome.gb -s sgRNA-library_spn-D39V.csv -o ~/output/ -t ~/wd/
```
where `genome.gb` is a GenBank file of the genome for the library is to be evaluated, `sgRNA-library_spn-D39V.csv` is a .csv file with the names and spacer sequences (PAM-distal to -proximal, i.e., 5'-3') of the sgRNA library to be evaluated (here our readily available _S. pneumoniae_ D39V library, file available in this repository), `~/output/` is the local directory to write the output files to and `~/wd/` is the local working directory in which **TINDRi.py** can be found.  

Depending on the genome size, number of annotated features, number of allowed mismatches in binding site identification and computer specifications (such as number of available logical processors), computation times will likely be in the order of 10-60 minutes.

### Output
The pipeline can return multiple files, depending on user input parameters:
* A .csv file summarizing per sgRNA the binding sites with first- and second-highest repression activity.
* A .csv file with all sgRNA - binding site combinations and scores.
* A .fasta file with all sgRNA spacer sequences.
* A .fasta file with all sgRNA binding site sequences on the genome.
* A .csv file with all sgRNA IDs and sequences, used as input for **TINDRi.py**.
* A .csv file with all binding site IDs and sequences, used as input for **TINDRi.py**.
* A .csv file with the full mismatch matrix for all sgRNA - binding site combinations, as generated by **TINDRi.py**.
Users will in general be after either or both of the first two files.

#### Meaning of variables in sgRNA evaluation output files
Not all variables exist in all output files. 

Variable | Description
--- | ---
`sgRNA_ID` | Unique sgRNA identifier, used within the pipeline.
`mismatches` | Comma-separated string representing mismatches (1) and matches (0) per within-spacer position (PAM-proximal to -distal)
`site_ID` | Unique binding site identifier, used within the pipeline.
`sgRNA_name` | Unique sgRNA name.
`n_mismatch` | Total number of mismatches between given sgRNA spacer and binding site.
`reprAct` | Estimated mismatch-based repression activity of a given sgRNA on a given site.
`chrom` | Name of the chromosome on which the binding site is located.
`strand` | Strand of the chromosome on which the binding site is located.
`start` | Start base pair of the binding site location on the given chromosome.
`end` | End base pair of the binding site location on the given chromosome.
`sgRNA_seq` | Sequence of the sgRNA spacer (PAM-distal to -proximal).
`site_seq` | Sequence of the binding site (PAM-distal to -proximal).
`site_PAM` | Proto-spacer Adjacent Motif (PAM) sequence of the binding site.
`sgRNA_GC` | GC content of the sgRNA spacer, as a ratio.
`all_targets` | Locus_tag keys of all features within whose boundaries the binding site is located on the non-template strand, if any; otherwise NA.
`all_target_names` | If available, names of all features within whose boundaries the binding site is located on the non-template strand. Otherwise locus_tag keys, if any; otherwise NA.
`all_chroms` | Chromosome names on which binding site is located.
`all_strands` | DNA strands on which binding site is located.
`all_ranges` | Genomic locations in bp (formatted as start..end) of the binding site.
`bad_seed` | If any, name of input `--bad_seeds` sequences detected in the sgRNA spacer; otherwise NA.
`cut_site` | If any, name of input `--cut_sites` sequences detected in the sgRNA spacer; otherwise NA.
`reprAct_1st` / `reprAct_2nd` | Estimated mismatch-based repression activity of the sgRNA at the binding sites with first- and second-highest estimated repression activities.
`nmismatch_1st` / `nmismatch_2nd` | Total number of mismatches between the sgRNA spacer and the sequences of the binding sites with first- and second-highest estimated repression activities.
`targets_1st` / `targets_2nd` | Respectively for the binding site with the first- (`reprAct_1st`) and second-highest (`reprAct_2nd`) estimated repression activity, the locus_tag keys of features within whose boundaries the binding site is located on the non-template strand, if any; otherwise NA.
`target_names_1st` / `target_names_2nd` | If available, names of features corresponding to `targets_1st` and `targets_2nd`, otherwise locus_tag keys. If none, NA.
`chroms_1st` / `chroms_2nd` | Chromosome name on which the binding sites with first- and second-highest estimated repression activities are located.
`strands_1st` / `strands_2nd` | DNA strands on which the binding sites with first- and second-highest estimated repression activities are located.
`ranges_1st` / `ranges_2nd` | Genomic location in bp (formatted as start..end) of the binding sites with first- and second-highest estimated repression activities.
`seq_1st` / `seq_2nd` | Sequences of the binding sites with first- and second-highest estimated repression activities.
`PAM_1st` / `PAM_2nd` | PAM sequences of the binding sites with first- and second-highest estimated repression activities.

## Within-R version

### Usage
Open the script **sgRNA_library_evaluation.R** (we recommend in the RStudio integrated development environment), adjust the input parameters in the first code section `#### 1. Settings ####` and then run the whole script with those input parameters. In its simplest form, only `input_genome`, `sgRNA_file`, `outdir` and `TINDRidir` need to be specified (under `## REQUIRED ##`). Depending on the genome size, number of annotated features, number of allowed mismatches in binding site identification, and computer specifications (such as the number of available logical processors for multi-core processing), computation times may be in the order of 10-60 minutes. 

### Input
Parameter | Description
--- | ---
`input_genome` | **required**; character. Either the full path to a GenBank file (extension .gb, .gbk, .gbf, or .gbff), or an NCBI assembly accession number (starting with either GCA_ \[GenBank\] or GCF_ \[RefSeq\], e.g. [GCA_003003495.1](https://www.ncbi.nlm.nih.gov/assembly/GCA_003003495.1/)).
`sgRNA_file` | **required**; character. Full path to a header-less .csv file with two comma-separated columns, the first containing sgRNA-unique IDs (e.g. sgRNA0001, sgRNA0002, etc.) and the second containing the sgRNA spacer sequences (PAM-distal to -proximal, i.e., 5 '-3'). As **sgRNA-library_spn-D39V.csv** in this repository.
`outdir` | **required** character. Path to desired output directory.
`TINDRidir` | **required** character. Path to directory in which **TINDRi.py** is stored.
`path_ncbi_downloads` | character. Path to directory where downloaded files from NCBI should be stored. Only **required if** `input_genome` is a GenBank file. Default: NA.
`regions` | integer vector. Indicates the sub-spacer regions for which to define max_mismatch_cum. The sum of these integers must equal the sgRNA spacer length. Default: c(7, 2, 11). 
`max_mismatch_cum` | integer vector. Indicates the maximum cumulative number of mismatches allowed over regions for a binding site to be considered relevant. Should be of the same length as `regions`. For example, max_mismatch_cum of c(1, 2, 11) combined with regions of c(7, 2, 11) allows for 1 mismatch in the first 7 nucleotides (PAM-proximal to -distal), 2 mismatches in the first 7+2=9 nucleotides and 11 mismatches in the first 7+2+11=20 nucleotides (which is the whole sgRNA spacer). Binding sites that do not fulfill these conditions for a given sgRNA are not regarded as relevant and therefore excluded (in the example, a site with e.g. 2 mismatches in the first 7 nucleotides with the spacer would not be considered relevant for that sgRNA, as its repression activity is likely too low at that site). Default: c(1, 2, 11).
`reprAct_penalties` | character or numeric vector. Indicates within-spacer position-dependent penalties for `reprAct` estimation. For using built-in penalties, must be one of "HawkinsEcoMean", "HawkinsBsuMean", "HawkinsEcoMedian", "HawkinsBsuMedian", "Qi", "QiMean" or "QiMedian". These vectors are the mean, median or raw 20-nt spacer relative repression activity values found by Hawkins et al. (2020) (means and medians per nucleotide) and Qi et al. (2013) (means and medians per region as defined in their paper) for _E. coli_ (Eco) and _B. subtilis_ (Bsu). For custom penalties, please specify a numeric vector (PAM-proximal to -distal) of length `sum(regions)` (i.e., sgRNA spacer length); if `pen_func` is `prod` (recommended), the integers should be on the \[0,1\] interval. Default: "HawkinsBsuMedian" (PAM-proximal to -distal: 0.35, 0.46, 0.5, 0.53, 0.57, 0.53, 0.57, 0.67, 0.83, 0.9, 0.91, 0.94, 0.95, 0.92, 0.95, 0.84, 0.94, 0.94, 0.96, 0.99).
`pen_func` | character or function. Defines how to process reprAct_penalties of sgRNA-site mismatches to produce `reprAct` estimates; strong advice not to change. Default: prod (product).
`bad_seeds` | character vector. Specifies bad seed sequences to flag. Default: c("ACCCA", "TGGAA"), based on Cui et al. (2018). 
`cut_sites` | named character vector. Indicates the motifs to detect in sgRNA spacer + primer sequences (`oligoForwardOverhang` & `oligoReverseOverhang`). Could be e.g. specific restriction enzyme cut sites, or tandem T nucleotides one wishes to avoid. Reverse-complement sequences are automatically detected too. Default: c(BsmBI = "CGTCTC", tandemT = "TTTTT"). 
`oligoForwardOverhang` | character. Post-annealing 5'-3' overhang for the forward oligo to be compatible with the digested sgRNA backbone vector. Default: `TATA` (BsmBI).
`oligoReverseOverhang` | character. Post-annealing 5'-3' overhang for the reverse oligo to be compatible with the digested sgRNA backbone vector. Default: `AAAC` (BsmBI).
`PAM` | character. Proto-spacer Adjacent Motif sequence. Default: "NGG" (_S. pyogenes_ dCas9).
`output_exact` | logical. Whether or not to output the main results summary output table as .csv file, containing per sgRNA bad seed and cut site flags, GC content and details on the binding sites with the first- and second-highest estimated maximum repression activity. Note that setting both this and `output_full` to FALSE is nonsensical. Default: TRUE.
`output_full` | logical. Whether or not to output a .csv file with the complete list of candidate sgRNA - binding site combinations & detailed information per combination. Note that for very large sgRNA libraries, this may take considerable computation time and resources and generate a large file. Default: FALSE.
`output_sgRNAs_fasta` | logical. Whether or not to output a .fasta file containing the sequences of all sgRNAs. Default: FALSE. 
`output_sites_fasta` | logical. Whether or not to output a .fasta file containing the sequences of all identified sgRNA binding sites on the genome. Default: FALSE. 
`detect_offtarget_genes_full` | logical. If `output_full` is set to TRUE, whether or not to also detect for each sgRNA-site combination if the site falls within an annotated feature on its non-template strand. Note that if set to TRUE, this might add computation time, especially for large sgRNA libraries. Default: FALSE.
`keep_TINDRi_input_sgRNAs` | logical. Whether or not to keep the intermittently generated .csv file containing all candidate sgRNA IDs and sequences (PAM-proximal to -distal), used as input for **TINDRi.py**. Default: FALSE.
`keep_TINDRi_input_sites` | logical. Whether or not to keep the intermittently generated .csv file containing all binding site IDs and sequences (PAM-proximal to -distal), used as input for **TINDRi.py**. Default: FALSE.
`keep_TINDRi_matches` | logical. Whether or not to keep the intermittently generated .csv file containing the full TINDRi.py mismatch output matrix. Note this file is generally large. Default: FALSE.
`path_python` | character. String indicating full path to Python3 installation (e.g. ".../anaconda3/python.exe" on Windows systems). Specify only if the Python3 installation is not automatically detected by the pipeline (it will issue an error message indicating this). Default: NULL (automatic detection). 
`feature_type` | character. GenBank key tag of features to be detected. Experimental pipeline feature, in principle do not alter. Default: `"locus_tag"`.

### Output
Same as described above for the command line version of the pipeline.


# Downstream efficiency analyses
The file **Pneumococcal_sgRNA_library_efficiency_evaluation.pdf** contains downstream analyses of the results of a former version of **sgRNA_library_evaluation.R** for seven pneumococcal genomes, given the sgRNA library **sgRNA-library_spn-D39V.csv** of Liu _et al._ (2020). For the binding site results, see <https://www.veeninglab.com/crispri-seq>. Source code of the pdf: **Pneumococcal_sgRNA_library_efficiency_evaluation.Rmd**.


# References
* Liu, X., Kimmey, J.M., Matarazzo, L., de Bakker, V., Van Maele, L., Sirard, J.-C., Nizet, V., Veening, J.-W. (2021) Exploration of Bacterial Bottlenecks and _Streptococcus pneumoniae_ pathogenesis by CRISPRi-seq. _Cell Host Microbe_ **29**, 107-120.e6. https://doi.org/10.1016/j.chom.2020.10.001 (https://www.sciencedirect.com/science/article/pii/S193131282030559X).
* Qi, L.S., Larson, M.H., Gilbert, L.A., Doudna, J.A., Weissman, J.S., Arkin, A.P., Lim, W.A. (2013) Repurposing CRISPR as an RNA-guided platform for sequence-specific control of gene expression. _Cell_ **152**, 1173–1183. <https://doi.org/10.1016/j.cell.2013.02.022>
* Hawkins, J.S., Silvis, M.R., Koo, B-M., Peters, J.M., Osadnik, H., Jost, M., Hearne, C.C., Weissman, J.S., Todor, H., Gross, C.A. (2020) Mismatch-CRISPRi Reveals the Co-varying Expression-Fitness Relationships of Essential Genes in Escherichia coli and Bacillus subtilis. _Cell Syst._ **11**, 523-535.e9. <https://doi.org/10.1016/j.cels.2020.09.009>
* Cui, L., Vigouroux, A., Rousset, F., Varet, H., Khanna, V., Bikard, D. (2018) A CRISPRi screen in E. coli reveals sequence-specific toxicity of dCas9. _Nat. Commun._ **9**, 1912. <http://dx.doi.org/10.1038/s41467-018-04209-5>
* Zhu, L.J., Holmes, B.R., Aronin, N., Brodsky, M.H. (2014) CRISPRseek: A Bioconductor Package to Identify Target-Specific Guide RNAs for CRISPR-Cas9 Genome-Editing Systems. _PLoS ONE_ **9**(9): e108424. <https://doi.org/10.1371/journal.pone.0108424>
* Drost, H.-G., Paszkowski, J. (2017) Biomartr: genomic data retrieval with R. _Bioinformatics_ **33**(8), 1216–1217. <https://doi.org/10.1093/bioinformatics/btw821>
* Trevor L Davis (2020). optparse: Command Line Option Parser. R package version 1.6.6. https://CRAN.R-project.org/package=optparse
* Gabriel Becker and Michael Lawrence (2020). genbankr: Parsing GenBank files into semantically useful objects. R package version 1.16.0. <https://doi.org/doi:10.18129/B9.bioc.genbankr>
* H. Pagès, P. Aboyoun, R. Gentleman and S. DebRoy (2020). Biostrings: Efficient manipulation of biological strings. R package version 2.56.0. <https://doi.org/doi:10.18129/B9.bioc.Biostrings>