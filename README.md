# Genome-wide sgRNA library design and evaluation for CRISPRi-seq
Author: Vincent de Bakker  
See also: <https://www.veeninglab.com/crispri-seq>  
Published as a part of: 
* Liu, X. _et al._ (2020) Exploration of bacterial bottlenecks and _Streptococcus pneumoniae_ pathogenesis by CRISPRi-seq. _Cell Host & Microbe_. https://doi.org/10.1016/j.chom.2020.10.001
* \<under review\>

This repository contains R pipelines to (1) design a _de novo_ sgRNA library and (2) evaluate all potential binding sites of a given sgRNA library on any NCBI genome. Additionally, it contains an example of possible downstream expected efficiency analyses for the latter. 

# Core function
For both types of analyses mentioned above, it is necessary to identify and score sgRNA binding sites. That is what the function `sgRNAefficiencyMC()` does. It uses the [`CRISPRseek` package](https://bioconductor.org/packages/release/bioc/html/CRISPRseek.html) and supports multicore processing with the `parallel` package.  

This function is called by both the sgRNA design and the sgRNA evaluation scripts. 

## Arguments
`sgRNAs` | DNAStringSet object containing a set of sgRNAs, as in `CRISPRseek::searchHits()`.
`genes` | data.frame object containing gene information extracted from GFF file with columns seqid, source, type, start, end, score, strand, phase and attribute.
`genome`
`reprAct`
`dist2SC`
`name_by`
`penalties`
`custom.penalties`
`outfile`
`no_cores`
`...` | other arguments to be passed to `CRISPRseek::searchHits()`.

## Output


# _de novo_ sgRNA design


# sgRNA binding site evaluation


# Downstream efficiency analyses


# References
* Liu, X. _et al._ (2020). Exploration of bacterial bottlenecks and _Streptococcus pneumoniae_ pathogenesis by CRISPRi-seq. _Cell Host & Microbe_. https://doi.org/10.1016/j.chom.2020.10.001
* Zhu, L.J., Holmes, B.R., Aronin, N., Brodsky, M.H. (2014) CRISPRseek: A Bioconductor Package to Identify Target-Specific Guide RNAs for
  CRISPR-Cas9 Genome-Editing Systems. _PLoS ONE_ **9**(9): e108424. doi:10.1371/journal.pone.0108424
