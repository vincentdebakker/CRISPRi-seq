# Genome-wide sgRNA library design and evaluation for CRISPRi-seq
Author: Vincent de Bakker  
See also: <https://www.veeninglab.com/crispri-seq>  
Published as a part of: 
* Liu, X. _et al._ (2020) Exploration of bacterial bottlenecks and _Streptococcus pneumoniae_ pathogenesis by CRISPRi-seq. _Cell Host & Microbe_. https://doi.org/10.1016/j.chom.2020.10.001
* \<under review\>

This repository contains R pipelines to (1) design a _de novo_ sgRNA library and (2) evaluate all potential binding sites of a given sgRNA library on any NCBI genome. Additionally, it contains an example of possible downstream expected efficiency analyses for the latter. 

# Core function
For both types of analyses mentioned above, it is necessary to identify and score sgRNA binding sites. That is what the function `sgRNAefficiencyMC()` does. It uses the [`CRISPRseek` package](https://bioconductor.org/packages/release/bioc/html/CRISPRseek.html) and supports multicore processing with the `parallel` package.  

> **Explain reprAct and dist2CS scores!**

This function is called by both the sgRNA design and the sgRNA evaluation scripts. 

## Input
Argument | Description
--- | ---
`sgRNAs` | DNAStringSet object containing a set of sgRNAs, as in `CRISPRseek::searchHits()`.
`genes` | data.frame object containing gene information extracted from GFF file with columns seqid, source, type, start, end, score, strand, phase and attribute.
`genome` | DNAStringSet object containing the set of chromosomes of the genome of interest. Can contain plasmids too. Each DNAString in the set is passed to `seqs` argument of `CRISPRseek::searchHits()`.
`reprAct` | logical. Should the expected relative retained repression activity of the sgRNAs for each of their binding sites be computed? See `penalties` for details on score calculation. Default: `TRUE`.
`dist2SC` | logical. If a binding site falls in the range of an annotated feature (see `name_by`), should the relative distance of that binding site to the start codon of that feature be computed? Default: `TRUE`.
`name_by` | character. Indicates the feature tag in `attributes` variable of `genes` to label the features by. Default: `"locus_tag"`.
`penalties` | character. One of `"qi.mean.per.region"`, `"qi"` or `"custom"`. Indicates the scoring system to quantify `reprAct`, if set to `TRUE`. `"qi"` uses the per-base penalties as estimated from Qi _et al._ (2013), Figure 5D ("Repression activity (relative)"). `"qi.mean.per.region"` uses the same penalties, averaged over the three regions indicated in that same Figure. Choose `"custom"` to define custom penalties. For every sgRNA - binding site combination, the `reprAct` score is defined as the product of the penalties for its mismatched bases. Here, we assume the observation of independence of mismatches (Qi _et al._, 2013) can be extended to >2 and non-adjacent mismatches. Default: `"qi.mean.per.region"`. 
`custom.penalties` | numeric vector of length 20. If `penalties` is set to `"custom"`, indicates custom penalty scores from PAM-distal to PAM-proximal. Default: `NULL`.
`outfile` | character. Path to directory to store temporary results files. Default: current working directory (`getwd()`). 
`no_cores` | integer. Indicates the maximum number of logical processors to use for multicore (parallel) processing. When set to `1`, multicore processing is disabled. Parallel processing is highly recommended due to otherwise long computation times. The number of available logical processors can be found with `parallel::detectCores()`. It is advised to use `parallel::detectCores() - 2`. Default: `1`.
`...` | other arguments to be passed to `CRISPRseek::searchHits()`.

## Output


# _de novo_ sgRNA design


# sgRNA binding site evaluation


# Downstream efficiency analyses


# References
* Liu, X. _et al._ (2020). Exploration of bacterial bottlenecks and _Streptococcus pneumoniae_ pathogenesis by CRISPRi-seq. _Cell Host & Microbe_. https://doi.org/10.1016/j.chom.2020.10.001
* Zhu, L.J., Holmes, B.R., Aronin, N., Brodsky, M.H. (2014) CRISPRseek: A Bioconductor Package to Identify Target-Specific Guide RNAs for
  CRISPR-Cas9 Genome-Editing Systems. _PLoS ONE_ **9**(9): e108424. https://doi.org/10.1371/journal.pone.0108424
