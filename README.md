# Genome-wide sgRNA library design and evaluation for CRISPRi-seq
Author: Vincent de Bakker  
Python code by Afonso Bravo    
See also: <https://www.veeninglab.com/crispri-seq>  
Published as a part of: 
* Liu, X., Kimmey, J.M., Matarazzo, L., de Bakker, V., Van Maele, L., Sirard, J.-C., Nizet, V., Veening, J.-W. (2020) Exploration of bacterial bottlenecks and _Streptococcus pneumoniae_ pathogenesis by CRISPRi-seq. _Cell Host & Microbe_. https://doi.org/10.1016/j.chom.2020.10.001
* \<submitted\>

This repository contains R pipelines to (1) design a _de novo_ sgRNA library and (2) evaluate all potential binding sites of a given sgRNA library on any NCBI genome. Additionally, it contains an example of possible downstream expected efficiency analyses for the latter. Lastly, it contains a subdirectory with the code to generate the naïve power analysis as published in \<submitted\>.

# Core function
For both types of analyses mentioned above, it is necessary to identify and score sgRNA binding sites. That is what the function `sgRNAefficiencyMC()` in the file **function_sgRNAefficiencyMC.R** does. It uses the [`CRISPRseek` package](https://bioconductor.org/packages/release/bioc/html/CRISPRseek.html) and supports multicore processing with the `parallel` package. In principle the user never has to open or adapt the script containing this function.  

This function is called by both the sgRNA design and the sgRNA evaluation scripts. 

## Input
Argument | Description
--- | ---
`sgRNAs` | DNAStringSet object containing a set of sgRNAs, passed to `gRNAs` in `CRISPRseek::searchHits()`.
`genes` | data.frame object containing gene information extracted from GFF file with columns seqid, source, type, start, end, score, strand, phase and attribute.
`genome` | DNAStringSet object containing the set of chromosomes of the genome of interest. Can contain plasmids too. Each DNAString in the set is passed to `seqs` argument of `CRISPRseek::searchHits()`.
`reprAct` | logical. Should the expected relative retained repression activity of the sgRNAs for each of their binding sites be computed? See `penalties` for details on score calculation. Default: `TRUE`.
`dist2SC` | logical. If a binding site falls (partially) within the range of an annotated feature (see `name_by`) on the genome, should the relative distance of that binding site to the start codon of that feature be computed? Default: `TRUE`.
`name_by` | character. Indicates the feature tag in `attributes` variable of `genes` to label the features by. Default: `"locus_tag"`.
`penalties` | character. One of `"qi.mean.per.region"`, `"qi"` or `"custom"`. Indicates the scoring system to quantify `reprAct`, if set to `TRUE`. `"qi"` uses the per-base mismatch penalties as estimated from Qi _et al._ (2013), Figure 5D ("Repression activity (relative)"). `"qi.mean.per.region"` uses the same penalties, averaged over the three regions indicated in that same Figure. Choose `"custom"` to define custom penalties. For every sgRNA - binding site combination, the `reprAct` score is defined as the product of the penalties for its mismatched bases. Here, we assume the observation of independence of mismatches (Qi _et al._, 2013) can be extended to >2 and non-adjacent mismatches. Default: `"qi.mean.per.region"`. 
`custom.penalties` | numeric vector of length 20, or `NULL`. If `penalties` is set to `"custom"`, indicates base pair-specific custom mismatch penalty scores ordered from from PAM-distal to PAM-proximal. Default: `NULL`.
`outfile` | character. Path to directory to store temporary results files. Default: current working directory (`getwd()`). 
`no_cores` | integer. Indicates the maximum number of logical processors to use for multicore (parallel) processing. When set to `1`, multicore processing is disabled. Parallel processing is highly recommended due to otherwise long computation times. The number of available logical processors can be found with `parallel::detectCores()`. It is advised to use `parallel::detectCores() - 2`. Default: `1`.
`...` | other arguments to be passed to `CRISPRseek::searchHits()`.

## Output
The function returns an R data.frame object in which every row is one unique binding site - sgRNA combination, with the following variables:  
Variable | Description
--- | ---
IsMismatch.pos1-20 | integer, 0 or 1. Indicates whether (1) or not (0) the base at this position in the spacer was mismatched with the binding site. 
strand | character, + or -. Indicates the strand on which strand of the chromosome the binding site is located. 
chrom | character. Indicates the chromosome on which the binding site is located. Relevant in case of plasmids in the genome. 
chromStart | numeric. First base pair of binding site, counting from the chromosome start.
chromEnd | numeric. Last base pair of binding site, counting from the chromosome start.
name | character. Name of sgRNA that binds at binding site; unique identifier. 
gRNAPlusPAM | character. Sequence of the sgRNA that binds at the site, including PAM.
targetSequence | character. Sequence of the binding site, including PAM.
n.mismatch | numeric. Total number of mismatched base pairs between the sgRNA spacer and the binding site.
site | numeric. Unique binding site identifier. One binding site overlapping e.g. two annotated features will appear as two rows in the table, with a shared site ID.
NTgene | character. If the binding site lies (partially) within an annotated feature, indicates that feature name. Otherwise `NA`.
coverPart | character. If the binding site overlaps (partially) with an annotated feature, indicates which part of the spacer does so: complete overlap, with 3' or 5' end of the spacer. Otherwise: `NA`. 
coverSize | numeric. Number of overlapping base pairs between spacer and `NTgene`, if any.
reprAct | numeric, on the \[0,1\] interval. Estimated relative retained repression activity, compared to an hypothetical zero-mismatch binding site for the sgRNA at the same distance from the transcription start site.
dist2SC | numeric, on the \[0,1\] interval. Relative distance of the binding site to the start codon of `NTgene`, if any. Here, 0 means the binding site is located on the start codon or has partial overlap with the 5'-end of the feature and 1 means the binding site is on the far end of the feature or has partial overlap with its 3'-end.

Except the last six, this output comes from calling `CRISPRseek::searchHits()`. 


# _de novo_ sgRNA design
The pipeline in the file **sgRNA_library_design.R** can be run to design an sgRNA library for the given genome and annotated features of choice. It identifies all possible candidate sgRNAs targeting the non-template strands of all features. Then it calls the core function `sgRNAefficiencyMC()` to identify for each of the candidate sgRNAs all possible binding sites on the genome. Using that data, the maximum off-target `reprAct` score is computed for each sgRNA, per feature. Finally, the pipeline first selects for each feature the candidate sgRNAs with the lowest expected maximum off-target repression activity (within some user-defined error range) and secondly, of that subset, the one sgRNA with the smallest `dist2CS`. In addition, corresponding forward and reverse oligo's are designed for each selected sgRNA for library construction.  

Candidate sgRNAs for every feature are found using `CRISPRseek::findgRNAs()` and genomes and annotations (GFF files) are imported from NCBI using the [biomartr package](https://cran.r-project.org/package=biomartr). The pipeline supports multicore processing.

## Usage
Open the script (we recommend in the RStudio integrated development environment), adjust the input parameters in the first code section `#### 1. Settings ####` and then run the whole script with those input parameters. Depending on the genome size, number of annotated features, number of allowed mismatches in binding site identification, number of cores used for parallel processing and computer specifications, computation times may be in the order of one to a few hours. Regarding the Python3 script TINDERi.py, as it is initialized from the command window, its highly advisable to download and install anaconda. This should ensure a proper working Python environment. All required dependencies can be seen at the top of the Python script as `import`. 

## Input
Parameter | Description
--- | ---
`wd` | character. Path to working directory. Output will appear here. Default: local R home directory.
`fd` | character. Path to directory containing the file **function_sgRNAefficiencyMC.R**. Default: `wd`.
`path_ncbi_downloads` | character. Path to directory where downloaded files from NCBI should be stored. Default: `wd`.
`accession_nr` | character. NCBI assembly accession number of genome for which to design sgRNAs. Default is _S. pneumoniae_ D39V: `"GCA_003003495.1"`. 
`db` | character. Data base from which to retrieve genome and annotation files. Either `"genbank"`, `"refseq"` or `"ensembl"`. Default: `"genbank"`.
`feature_type` | character. Feature key tag for each of whose unique identifiers the sgRNAs should be designed. Has to appear in the attributes columns of the GFF file. Also passed to the `name_by` argument of `sgRNAefficiencyMC()`. Defaults to `"locus_tag"`.
`n_cores` | integer. Maximum number of logical processors to use. If set to 1, multi-core processing is disabled. Must be a strictly positive integer. Available number of logical cores can be detected from within R with `parallel::detectCores()`. It is recommended to always use one or two cores less than available to allow other processes to run as well. Also passed to the `no_cores` argument of `sgRNAefficiencyMC()`. Default: `1`.
`max_mismatch` | integer. Maximum number of mismatches allowed between sgRNA spacer and binding site when searching for binding sites on the genome. Passed as `max.mismatch` argument to `sgRNAefficiencyMC()`. Default: `6`.
`reprAct_penalties` | character. Scoring system for repression activity estimation. One of `"qi.mean.per.region"`, `"qi"` or `"custom"`. Passed to the `penalties` argument of `sgRNAefficiencyMC()`, see function input description above for details. Default: `"qi.mean.per.region"`.
`reprAct_custom_penalties` | numeric vector of length 20, or `NULL`. Custom mismatch penalties for repression activity estimation. Passed to the `custom.penalties` argument of `sgRNAefficiencyMC()`, see function input description above for details. Default: `NULL`.
`errorRange_maxOffreprAct` | numeric. The error allowed on maximum off-target repression activity estimates. Candidate sgRNAs with an estimated repression activity within this range from the minimum found off-target repression activity all pass the specificity filter in selecting one optimal sgRNA per feature. Default is 1%: `0.01`.
`oligoForwardOverhang` | character. Post-annealing 5'-3' overhang for the forward oligo to be compatible with the digested sgRNA backbone vector. Default: `TATA`. 
`oligoReverseOverhang` | character. Post-annealing 5'-3' overhang for the reverse oligo to be compatible with the digested sgRNA backbone vector. Default: `AAAC`. 
`spacer_length` | integer. Length of the sgRNA spacers to be designed in bp. Must be a strictly positive integer. Default: `20`. 
`PAM` | character. Proto-spacer Adjacent Motif sequence. Default: `"NGG"`.
`output_target_fasta` | logical. Whether to output a .fasta file containing all sequences for which sgRNAs are designed. Default: `TRUE`.
`output_full_list` | logical. Whether to output a .csv file containing all found binding sites for all found candidate sgRNAs. Default: `TRUE`.
`output_optimized_list` | logical. Whether to compute and output a list of one optimal sgRNA per feature in .csv format. Default: `TRUE`.

## Output
The pipeline can return three files, depending on user input parameters: 
1. A .fasta file with all sequences for which sgRNAs were designed.
2. A .csv file with all identified binding sites for all identified candidate sgRNAs. Results as returned by `sgRNAefficiencyMC()`.
3. A .csv file with one optimal sgRNA per annotated feature. Filtered version of number 2. Includes an extra column `maxOffreprAct`, indicating the maximum off-target `reprAct` score found for that sgRNA - binding site combination. Also includes four extra columns containing forward (`forward`) and reverse (`reverse`) sgRNA base-pairing strands for sgRNA library construction, and forward (`oligoForward`) and reverse (`oligoReverse`) oligo's to be ordered to construct the annealed sgRNA library pool. 


# sgRNA binding site evaluation
The pipeline in the file **sgRNA_library_evaluation.R** identifies and scores all binding sites up to a specified number of mismatches in a given genome for each sgRNA of a given library. It is the result of a call to `sgRNAefficiencyMC()`.  

Genomes and annotations (GFF files) are imported from NCBI using the [biomartr package](https://cran.r-project.org/package=biomartr). The pipeline supports multicore processing.

## Usage
Open the script (we recommend in the RStudio integrated development environment), adjust the input parameters in the first code section `#### 1. Settings ####` and then run the whole script with those input parameters. Depending on the genome size, number of annotated features, number of allowed mismatches in binding site identification, number of cores used for parallel processing and computer specifications, computation times may be in the order of 20 minutes to a few hours. 

## Input
Parameter | Description
--- | ---
`fd` | character. Path to directory containing the file **function_sgRNAefficiencyMC.R**. Default: local R home directory.
`path_ncbi_downloads` | character. Path to directory where downloaded files from NCBI should be stored. Default: `fd`.
`sgRNA_file` | character. Path to the .csv file containing the sgRNA information. This file should strictly have no header and two columns, separated by a comma, containing the sgRNA names and sequences without PAM, respectively. Default is the example file **sgRNA-library_spn-D39V.csv** with the sgRNA library of Liu _et al._ (2020), saved in the `fd` folder: `paste0(fd, "/sgRNA-library_spn-D39V.csv")`.
`accession_nr` | character. NCBI assembly accession number of genome for which to evaluate sgRNAs. Default is _S. pneumoniae_ D39V: `"GCA_003003495.1"`. 
`out_path` | character. Output will be written here. Also passed to the `outfile` argument of `sgRNAefficiencyMC()`. Default: `fd`.
`max_mismatch` | integer. Maximum number of mismatches allowed between sgRNA spacer and binding site when searching for binding sites on the genome. Passed as `max.mismatch` argument to `sgRNAefficiencyMC()`. Default: `8`.
`n_cores` | integer. Maximum number of logical processors to use. Passed to the `no_cores` argument of `sgRNAefficiencyMC()`. Default: `1`.
`feature_type` | character. Feature key tag for each of whom to detect if binding sites fall in their range on the genome. Has to appear in the attributes columns of the GFF file. Passed to the `name_by` argument of `sgRNAefficiencyMC()`. Defaults to `"locus_tag"`.
`reprAct_penalties` | character. Scoring system for repression activity estimation. One of `"qi.mean.per.region"`, `"qi"` or `"custom"`. Passed to the `penalties` argument of `sgRNAefficiencyMC()`. Default: `"qi.mean.per.region"`.
`reprAct_custom_penalties` | numeric vector of length 20, or `NULL`. Custom mismatch penalties for repression activity estimation. Passed to the `custom.penalties` argument of `sgRNAefficiencyMC()`. Default: `NULL`.
`spacer_length` | integer. Length of the sgRNA spacers in bp. Must be a strictly positive integer. Default: `20`. 
`PAM` | character. Proto-spacer Adjacent Motif sequence. Default: `"NGG"`.
`db` | character. Data base from which to retrieve genome and annotation files. Either `"genbank"`, `"refseq"` or `"ensembl"`. Default: `"genbank"`.

## Output
The pipeline writes the results of the `sgRNAefficiencyMC()` call to a .csv file.


# Downstream efficiency analyses
The file **Pneumococcal_sgRNA_library_efficiency_evaluation.pdf** contains downstream analyses of the results of **sgRNA_library_evaluation.R** for seven penumococcal genomes, given the sgRNA library **sgRNA-library_spn-D39V.csv** of Liu _et al._ (2020). For the binding site results, see <https://www.veeninglab.com/crispri-seq>. Source code of the pdf: **Pneumococcal_sgRNA_library_efficiency_evaluation.Rmd**.


# References
* Liu, X., Kimmey, J.M., Matarazzo, L., de Bakker, V., Van Maele, L., Sirard, J.-C., Nizet, V., Veening, J.-W. (2020). Exploration of bacterial bottlenecks and _Streptococcus pneumoniae_ pathogenesis by CRISPRi-seq. _Cell Host & Microbe_. <https://doi.org/10.1016/j.chom.2020.10.001>
* Qi, L.S., Larson, M.H., Gilbert, L.A., Doudna, J.A., Weissman, J.S., Arkin, A.P., Lim, W.A. (2013) Repurposing CRISPR as an RNA-guided platform for sequence-specific control of gene expression. _Cell_ **152**, 1173–1183. <https://doi.org/10.1016/j.cell.2013.02.022>
* Zhu, L.J., Holmes, B.R., Aronin, N., Brodsky, M.H. (2014) CRISPRseek: A Bioconductor Package to Identify Target-Specific Guide RNAs for CRISPR-Cas9 Genome-Editing Systems. _PLoS ONE_ **9**(9): e108424. <https://doi.org/10.1371/journal.pone.0108424>
* Drost, H.-G., Paszkowski, J. (2017) Biomartr: genomic data retrieval with R. _Bioinformatics_ **33**(8), 1216–1217. <https://doi.org/10.1093/bioinformatics/btw821>
