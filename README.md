# Bioinformatic Scripts used by George S. Long

The scripts located here have been used for the following projects:
1. [Collapse of the mammoth-steppe in central Yukon as revealed by ancient environmental DNA](https://doi.org/10.1038/s41467-021-27439-6) <-- This repo isn't listed as it was was created after the study was published and results were only present in the supplementary materials. Only MapDamageLoop was used
2. [Ancient Escherichia coli](https://github.com/longg2/AncientEcoli)
3. [Ancient Brucella melitensis](https://github.com/longg2/Brucella)
4. ???

Please note that there are some ideosyncracies present in the scripts due to the specifics of the machine I'm using (take the [Core SNP Phylogeny](https://github.com/longg2/LongBioinformatics/blob/master/CoreSNPPhylogenetics.sh) as an example). You'll need to edit scripts as needed to ensure that they'll work for you.

Also, there's some very liberal use of [GNU Parallel](https://www.gnu.org/software/parallel/).

Things that need work:
1. Symbolic links aren't followed by *all* of the scripts. Need to figure out why that's the case

# Ancient Escherichia coli
The order of the scripts used is as follows:
1. AncientQC.sh
2. BWAalnMapping.sh <-- Used to map agaisnt the Human Genome (with Dedup Flag)
  * MapDamageLoop.sh
3. PanGenomeCreation.sh
4. BWAalnMapping.sh <-- Used to map agaisnt the _E.coli_ Genome (with Dedup Flag)
  * MapDamageLoop.sh
  * CoreSNPPhylogenetics.sh <-- Multiple times for global, A0, and ST4995 phylogenies. [Treemmer](https://github.com/fmenardo/Treemmer) is ran manually after the fact.
5. Manually merging the deduplicated digests
6. BWAalnMapping.sh <-- Used to map agaisnt the _K.aerogenes_ Genome (with Dedup Flag)
  * MapDamageLoop.sh

The deduplicated samtools depth results were used in the paper instead of the summarized version. This is because the summarized version was created after the paper was submitted.
