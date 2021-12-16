## Selection analysis on the Cactus alignment

This forder contains scripts to perform the selection analysis.

The main package used to detect sites under selection is [PHAST](http://compgen.cshl.edu/phast/). 

([neutral_model_estimation.txt](neutral_model_estimation.txt)): contains the steps used to estimate the 4-fold degenerate sites with [phyloFit](http://compgen.cshl.edu/phast/help-pages/phyloFit.txt).
Steps were taken from [PHASTCONS HOW TO](http://compgen.cshl.edu/phast/phastCons-HOWTO.html) (see 2.2 Obtaining Phylogenetic Models)

[phastCons_analysis.txt](phastCons_analysis.txt) contains the steps used to detect conserved elements (CEs) in the barn swallow genome using [phastCons](http://compgen.cshl.edu/phast/help-pages/phastCons.txt).
Steps were taken from [PHASTCONS HOW TO](http://compgen.cshl.edu/phast/phastCons-HOWTO.html) (see 4.1 Many Sites, Not Too Many Sequences)

[phyloP_analysis.txt](phyloP_analysis.txt) contains the steps used to detect conserved and accelerated sites in the barn swallow genome using [phyloP](http://compgen.cshl.edu/phast/help-pages/phyloP.txt).

Other tools used were: [hal2maf](https://github.com/ComparativeGenomicsToolkit/hal/blob/master/maf/impl/hal2maf.cpp), [maf_stream](https://github.com/joelarmstrong/maf_stream) and [msa_view](http://compgen.cshl.edu/phast/help-pages/msa_view.txt), 

