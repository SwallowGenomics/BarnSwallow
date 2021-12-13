## Recall rate-downsampling experiments joint vs per-sample call

The R file present in this folder contains the commands used to plot the number of variants identified after downsampling, comparing joint vs per-sample variant calling. <br />
Datasets are divided according to the type of variant identified: SVs (structural variants), indels (insertions-deletions), SNVs (single nucleotide variants)  <br />
`down` files contain variants number in the downsampled dataset (5x); `full` files contain variants number in the full coverage dataset. <br />
The legend can be drawn with the [ComplexHeatmap](https://jokergoo.github.io/ComplexHeatmap-reference/book/) package using this [script](https://github.com/SwallowGenomics/BarnSwallow/blob/main/Plots%20and%20figures/recall%20rate-downsampling%20experiments/joint%20vs%20per-sample%20call/legend_script.r). 