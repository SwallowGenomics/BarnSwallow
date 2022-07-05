## Pangenome assembly and ortholog analysis

This folder contains the scripts and commands used to assembly the barn swallow pangenome.

[pre-processing.txt](https://github.com/SwallowGenomics/BarnSwallow/blob/main/Analyses/Pangenome/pre-processing.txt) contains the steps performed before the assembly of the pangenome, including adaptor trimming with [cutadapt](https://cutadapt.readthedocs.io/en/stable/), haplotig purging with [purge_dups](https://github.com/dfguan/purge_dups), [Genomescope2.0](http://qb.cshl.edu/genomescope/genomescope2.0/) and [Merqury](https://github.com/marbl/merqury), HiFi reads assembly with [HiFiasm](https://github.com/chhylp123/hifiasm), stats with [asm_stats](https://github.com/VGP/vgp-assembly/blob/master/pipeline/stats/asm_stats.sh) and masking with [WindowMasker](https://github.com/goeckslab/WindowMasker) and [RepeatMasker](https://github.com/rmhubley/RepeatMasker).
 
[Pangenome_assembly.txt](https://github.com/SwallowGenomics/BarnSwallow/blob/main/Analyses/Pangenome/Cactus_pangenome_pipeline.txt) containes the commands used to assemble the pangenome. The [Cactus pangenome pipeline documentation](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md) was followed.

The final pangenome was modified and indexed with [vg](https://github.com/vgteam/vg) following the GitHub solution [#132](https://github.com/vgteam/sequenceTubeMap/issues/132)
