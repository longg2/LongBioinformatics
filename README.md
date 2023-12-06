# Bioinformatic Scripts used by George Long
The scripts located here have been used for the following projects:
1. [Collapse of the mammoth-steppe in central Yukon as revealed by ancient environmental DNA](https://doi.org/10.1038/s41467-021-27439-6) <-- Not cited as repo was created after the paper was published and results were only present in the supplementary materials. Only MapDamageLoop was used
2. [Ancient _Escherichia coli_](https://github.com/longg2/AncientEcoli)
3. [Ancient _Brucella melitensis_](https://github.com/longg2/Brucella)
4. [Permafrost microbial communities follow shifts in vegetation, soils, and megafauna extinctions in Late Pleistocene NW North America](https://doi.org/10.1002/edn3.493)

Please note that there are some idiosyncrasies present in the scripts due to the specifics of the machine I'm using (take the [Core SNP Phylogeny](https://github.com/longg2/LongBioinformatics/blob/master/CoreSNPPhylogenetics.sh) as an example). You'll need to edit scripts as needed to ensure that they'll work for you.

# Ancient Escherichia coli
The order of the scripts used is as follows:
1. AncientQC.sh
2. BWAalnMapping.sh <-- Used to map against the Human Genome (with Dedup Flag)
  * MapDamageLoop.sh
3. PanGenomeCreation.sh
4. BWAalnMapping.sh <-- Used to map against the _E.coli_ Genome (with Dedup Flag)
  * MapDamageLoop.sh
  * CoreSNPPhylogenetics.sh <-- Multiple times for global, A0, and ST4995 phylogenies. [Treemmer](https://github.com/fmenardo/Treemmer) is ran manually after the fact.
5. Manually merging the de-duplicated digests
6. BWAalnMapping.sh <-- Used to map against the _K.aerogenes_ Genome (with Dedup Flag)
  * MapDamageLoop.sh

# Note about _No Longer Used_ Folder
These scripts are no longer used for a variety of reasons (either created for testing only or fundamental issues with their underlying programs). They should still work, however, they're not currently being actively maintained.
