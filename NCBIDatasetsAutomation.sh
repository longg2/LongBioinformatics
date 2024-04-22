#! /usr/bin/env bash

# Getting the commands I neeed
script_name=$0
script_full_path=$(dirname $0)

source $script_full_path/lib/BasicCommands.sh # This loads the basic things I need.
usage() { printf 'NCBI Download Automation
	This script will download RefSeq (Chromosomal and Complete) genomes of your chosen
	taxa. It will also extract the dates and arrange the files in the proper order. Requires
	the NCBI datasets program and NCBIDateParsing.awk to run properly.

        -t\tTaxa of interest
        -h\tShow this help message and exit\n' 1>&2; exit 1; }

while getopts "t:h" arg; do
        case $arg in
                t)
			taxa="${OPTARG}"
                        #echo "The raw sequencing files are located in $raw"
                        ;;
                h | *)
                        usage
                        exit 0
                        ;;
        esac
done

if [ -z ${taxa+x} ]
then
	printf "No taxa listed\n"
	usage
	exit 1
fi

identifier="T${taxa}D$(date +"%Y%m%d")"
# This downloads the genomes, gff, and gbk files from NCBI
datasets download genome taxon $taxa --assembly-source RefSeq \
       	--filename RefseqGenomes${identifier}.zip \
       	--assembly-level scaffold,chromosome,complete --exclude-atypical \
       	--include genome,gff3,gbff

# Now, we'll start with extracting the dates. At least for me, these are fairly important for my phylogenies.
# The rationale for the following is that collection dates should take priority over everything else, but,
# they're not always present. If so, we'll take the earliest of either the assembly date or the biosample 
# submission date. Some of the older assemblies were uploaded years before their associated biosample

dataformat tsv genome --package RefseqGenomesT${taxa}D$(date +"%Y%m%d").zip \
	--fields accession,assminfo-biosample-publication-date,assminfo-release-date,assminfo-biosample-attribute-name,assminfo-biosample-attribute-value > ${identifier}DatesRaw.tab

dataformat tsv genome --package RefseqGenomesT${taxa}D$(date +"%Y%m%d").zip \
	--fields accession,assminfo-biosample-attribute-name,assminfo-biosample-attribute-value > ${identifier}OtherMetaData.tab
	#
# Loading the two different version in memory
mapfile -t collection < <(grep "collection_date" ${identifier}DatesRaw.tab)
mapfile -t nocollection < <(cut -f 1-3 ${identifier}DatesRaw.tab | tail -n +2 | sort -u)

# Now to get the dates extracted
$script_full_path/NCBIDateParsing.awk <(printf "%s\n" "${collection[@]}") <(printf "%s\n" "${nocollection[@]}") > ${identifier}DatesCompiled.tab

# Now that we have the dates extracted, let's get the data out in the open
mkdir -p ${identifier}Genomes
#mkdir ${identifier}GFF
#mkdir ${identifier}GBK
unzip -j RefseqGenomes${identifier}.zip ncbi_dataset/data/GCF*/*fna -d ${identifier}Genomes > /dev/null 2> /dev/null
#unzip -j RefseqGenomes${identifier}.zip ncbi_dataset/data/GCF*/*gff -d ${identifier}GFF
#unzip -j RefseqGenomes${identifier}.zip ncbi_dataset/data/GCF*/*gbff -d ${identifier}GBK

# Some quick niceties here for myself
#rename .gff3 .gff ${identifier}GFF/*.gff3
#rename .gbff .gbk ${identifier}GFF/*.gbff

printf "Data from $taxa has been successfully downloaded!\n"
