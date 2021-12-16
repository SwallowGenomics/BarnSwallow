## Cactus alignment

This folder contains the species and tree input file for Cactus ([SeqFile.txt](SeqFile.txt)). <br />

Cactus was run with this command:

`cactus --logInfo --logError --binariesMode local --workDir=workDir jobStore SeqFile.txt 10_genomes.hal` <br />

The next steps were computed with [Hal Toolkit package](https://github.com/ComparativeGenomicsToolkit/hal).  <br />

The alignment was evaluated with: `HalValidate 10_genomes.hal`  <br />

Coverage for FIGURE 1 was calculated with:

`halAlignmentDepth --noAncestors --step 200000 --targetGenomes Camarhynchus_parvulus,Gallus_gallus,Lonchura_striata,Molothrus_ater,Motacilla_alba,Taeniopygia_guttata,Passer_domesticus --outWiggle coverage_10_genomes_200k.wig 10_genomes.hal Hirundo_rustica`  <br />

Basewise coverage was calculated with the following command for each barn swallow chromosome:

`halAlignmentDepth --noAncestors --targetGenomes Camarhynchus_parvulus,Gallus_gallus,Lonchura_striata,Molothrus_ater,Motacilla_alba,Taeniopygia_guttata,Passer_domesticus --outWiggle coverage_10_genomes_chr1.wig --refSequence chr1 10_genomes.hal Hirundo_rustica`  <br /> 
