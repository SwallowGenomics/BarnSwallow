## Cactus alignment

This folder contains the species and tree input file for Cactus ([SeqFile.txt](SeqFile.txt)). <br />

Cactus was run with this command:

`cactus --logInfo --logError --binariesMode local --workDir=workDir jobStore SeqFile.txt 10_genomes.hal` <br />

The next steps were computed with [Hal Toolkit package](https://github.com/ComparativeGenomicsToolkit/hal).  <br />

The alignment was evaluated with: `HalValidate 10_genomes.hal`  <br />

To calculate how many bases of each genome aligned to bHirRus1, the next commands were used (example on one species is shown):<br />

`halAlignmentDepth --noAncestors --outWiggle coverage_10_genomes_Camarhynchus.wig --targetGenomes Camarhynchus_parvulus 10_genomes.hal Hirundo_rustica`<br />
`wig2bed < coverage_10_genomes_Camarhynchus.wig > coverage_10_genomes_Camarhynchus.bed`<br />
`grep "1" coverage_10_genomes_Camarhynchus.wi | grep -v "fixed" | wc -l`

Basewise coverage was calculated with the following command for each barn swallow chromosome (excluding the two species that didn't align):

`halAlignmentDepth --noAncestors --targetGenomes Camarhynchus_parvulus,Gallus_gallus,Lonchura_striata,Molothrus_ater,Motacilla_alba,Taeniopygia_guttata,Passer_domesticus --outWiggle coverage_10_genomes_chr1.wig --refSequence chr1 10_genomes.hal Hirundo_rustica`  <br /> 

Coverage was calculated with the following command for each barn swallow chromosome on 200kb windows (excluding the two species that didn't align):

`halAlignmentDepth --noAncestors --step 200000 --targetGenomes Camarhynchus_parvulus,Gallus_gallus,Lonchura_striata,Molothrus_ater,Motacilla_alba,Taeniopygia_guttata,Passer_domesticus --outWiggle coverage_10_genomes_chr1_200k.wig --refSequence SUPER_1 ../../10_genomes.hal Hirundo_rustica`
