#! /usr/bin/env bash
# TODO: PE reads need to be properly handled

script_name=$0
script_full_path=$(dirname $0)

source $script_full_path/lib/BasicCommands.sh # This loads the basic things I need.
source $script_full_path/lib/CoreSNPFunctions.sh # This loads the basic things I need.

usage() { printf "Core SNP phylogeny creation V3.0
	Given a folder of sequences, create a core snp phylogeny using:	Snippy,
       	Gubbins, snp-sites, and IQTREE-2. The nucleotide model will be
	automatically chosen with ascertainment bias.

	V3 Changes: Now using an updated version of Gubbin which runs IQTREE
	in the background. Using it to keep everything consistent. Will
	continue to use a separate run of IQTREE. Decrease the default from
	10 to 8 CPU threads.

	V2 Changes: Snippy manifest could not handle bam files apparently? Due
	to this, I've created my own implementation of the script which uses
	parallel and file extension detection. Have also increased the default
	CPU threads from 8 to 10.

	NOTE: Doesn't handle PE samples for now. Assumes merged,
	assemblies, or full genomes

        -i\tWhere the sequences are
	-o\tOutput Prefix (Default = Phylogenetics)
	-r\tReference Sequence
	-g\tOutgroup Sequence
	-n\tNumber of CPU Threads to be used (Default = 10)
	-l\tLog File Name
	-m\tMinimum Read Coverage for SNP calling (Default = 10)
	-f\tFilter Percentage for Gubbins [0,100] (Default = 25)
	-b\tNumber of Bootstraps for Phylogeny (Default = 100)
        -h\tShow this help message and exit\n" 1>&2; exit 0; }
log() {	printf "Phylogeny settings for $(date):
	Log File:\t${log}
	Input folder:\t${folder}
	Output folder:\t${out}
	Reference:\t${reference}
	Outgroup:\t${outgroup}
	CPU Threads:\t${ncores}
	-------------------------------------
	Minimum Read Depth:\t${mincov}
	Filter Percentage:\t${filter}
	Bootstraps:\t${bootstrap}
	-------------------------------------\n"; exit 0;
}

log="$(date +'%Y%m%d').log"
declare -i ncores=10
filter=25
export out="Phylogenetics"
bootstrap=1000
mincov=10
while getopts "i:r:f:n:g:o:m:l:b:h" arg; do
        case $arg in
                i)
                        declare -r folder=${OPTARG}
                        ;;
                r)
                        declare -r reference=${OPTARG}
			export reference
                        ;;
                g)
                        declare outgroupFiles=${OPTARG}
			outgroup=$(echo $outgroupFiles | tr "," "\n" | sed -e "s/^.*\///" -e "s/\.f.*//" | tr "\n" "," | sed "s/,$//")
			#outgroup=$(basename $outgroup .fna)
			
			#export outgroup
                        ;;
                n)
                        declare -i ncores=${OPTARG}
                        ;;
		o)
			export out=${OPTARG}
			;;
		b)
			bootstrap=${OPTARG}
			;;
		l)
                        log=${OPTARG}
			;;
		f)
			filter=${OPTARG}
			;;
		m)
			mincov=${OPTARG}
			export mincov
			;;
                h | *)
                        usage
                        exit 0
                        ;;
        esac
done

if [ -z ${folder+x} ] || [ -z ${reference+x} ]; then
	printf "You are missing either the Input folder or the Reference Genome\n\n"
	usage
	exit 1
elif
	[ -z ${outgroup+x} ]; then
	printf "Assuming that the reference is also the outgroup\n"
	outgroup="Reference"
	#outgroup=$reference
	#outgroup=$(basename $outgroup .fna)
	export outgroup
fi

log | tee $log
printf "\n"

# Making the Directories
mkdir -p ${out}Snippy
mkdir -p ${out}Gubbins
mkdir -p ${out}IQTREE
export -f snippyParallel

# Activating conda environments so that we can switch when needed
#source ~/miniconda3/etc/profile.d/conda.sh # Activating 
#conda activate base
#conda activate pangenome

############# Running SNIPPY ###############
# Getting the SNPs
parJobs=$(echo "scale=0;var1=$ncores/8;var1"|bc) # Will round down!!!

parallel -j $parJobs --bar "snippyParallel {} $mincov $reference ${out}Snippy 2> /dev/null" ::: ${folder}/*

snippy-core --ref $reference ${out}Snippy/*

mv core* ${out}Snippy/

snippy-clean_full_aln ${out}Snippy/core.full.aln > ${out}Snippy/clean.full.aln

########## Time for Gubbins #############
echo "Running Gubbins"
#conda activate phylogenies # Need to activate

# Running Gubbins is relatively simple, if potentially long...
#run_gubbins.py --outgroup $outgroup --threads $ncores \
#	--filter-percentage $filter --tree-builder iqtree \
#       	--bootstrap $bootstrap --best-model --first-model GTRGAMMA \ 
#       	--prefix ${out}Gubbins/RecombMask ${out}Snippy/clean.full.aln > ${out}Gubbins.log
	#
#
run_gubbins.py --outgroup $outgroup --threads $ncores \
	--filter-percentage $filter \
       	--tree-builder iqtree \
	--bootstrap $bootstrap --best-model --first-model GTRGAMMA \
	--prefix ${out}Gubbins/RecombMask ${out}Snippy/clean.full.aln

if [ $? -eq 1 ]; then
	echo "Gubbins failed. Please look at the logs to figure out where"
	exit 1
fi
#
snp-sites -c ${out}Gubbins/RecombMask.filtered_polymorphic_sites.fasta > ${out}Gubbins/clean.core.aln
#
################ IQ TREE 2 ################
echo "Building the phylogeny"
iqtree2 -s ${out}Gubbins/clean.core.aln -o Reference -m MFP+ASC \
	-T AUTO --threads-max $ncores \
	-B $bootstrap --prefix ${out}IQTREE/Phylogeny

if [ $? -eq 1 ]; then
	echo "IQTREE failed. Please look at the logs to figure out where"
	exit 1
else
	echo "A core SNP phylogeny has been created!"
	exit 0
fi
