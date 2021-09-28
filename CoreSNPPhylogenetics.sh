#! /usr/bin/env bash
# TODO: PE reads need to be properly handled

script_name=$0
script_full_path=$(dirname $0)

source $script_full_path/lib/BasicCommands.sh # This loads the basic things I need.

usage() { printf 'Core SNP phylogeny creation.
	Given a folder of sequences, it will create a core snp phylogeny using:
	Snippy, Gubbins, snp-sites, and IQTREE-2. The nucleotide model will be
	automatically chosen with ascertainment bias. One thousand bootstraps
	will be used.
	NOTE: Doesnt handle PE samples for now. Assumes merged,
	assemblies, or full genomes

        -i\tWhere the sequences are
	-o\tOutput Prefix (Default = Phylogenetics)
	-r\tReference Sequence
	-n\tNumber of CPU Threads to be used (Default = 8)
	-l\tLog File Name
	-f\tFilter Percentage for Gubbins [0,100] (Default = 25)
        -h\tShow this help message and exit\n' 1>&2; exit 0; }
log() {	printf "Phylogeny settings for $(date):
	Log File:\t${log}
	Input folder:\t${folder}
	Output folder:\t${out}
	Reference:\t${db}
	CPU Threads:\t${ncores}
	-------------------------------------
	Filter Percentage:\t${filter}
	-------------------------------------\n"; exit 0;
}

log="$(date +'%Y%m%d').log"
declare -i ncores=8
filter=25
export out="Phylogenetics"
while getopts "i:r:f:n:o:l:h" arg; do
        case $arg in
                i)
                        declare -r folder=${OPTARG}
                        ;;
                r)
                        declare -r reference=${OPTARG}
                        ;;
                n)
                        declare -i ncores=${OPTARG}
                        ;;
		o)
			export out=${OPTARG}
			;;
		l)
                        log=${OPTARG}
			;;
		f)
			filter=${OPTARG}
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
fi

log | tee $log
printf "\n"

# Making the Directories
export reference
mkdir -p ${out}Snippy
mkdir -p ${out}Gubbins
mkdir -p ${out}IQTREE

# Activating conda environments so that we can switch when needed
source ~/miniconda3/etc/profile.d/conda.sh # Activating 
conda activate base
conda activate pangenome

############# Running SNIPPY ###############
# First, we need to make a manifest file for snippy-multi
find ~+/$folder -type f | $script_full_path/lib/snippyManifest.awk > ${out}snippyManifest.txt

# Running Snippy. Because of it functions, we'll be piping the log to nowhere
echo "Running snippy-multi"

cd ${out}Snippy
tmp=$(echo $OLDPWD) # So that I can work without issue...
snippy-multi $tmp/${out}snippyManifest.txt --ref $tmp/$reference --cpus $ncores --quiet | bash 2> $tmp/${out}SnippyLog.log

if [ $? -eq 1 ]; then
	echo "Snippy failed. Please look at the logs to figure out where"
	exit 1
fi

snippy-clean_full_aln core.full.aln > clean.full.aln

cd $tmp
unset tmp

########## Time for Gubbins #############
echo "Running Gubbins and snp-sites"
conda activate phylogenies # Need to activate

# Running Gubbins is relatively simple, if potentially long...
run_gubbins.py ${out}Snippy/clean.full.aln --outgroup Reference --threads $ncores \
	--filter_percentage $filter --prefix ${out}Gubbins/RecombMask > ${out}Gubbins.log

if [ $? -eq 1 ]; then
	echo "Gubbins failed. Please look at the logs to figure out where"
	exit 1
fi

snp-sites -c ${out}Gubbins/RecombMask.filtered_polymorphic_sites.fasta > ${out}Gubbins/clean.core.aln

############### IQ TREE 2 ################
echo "Building the phylogeny"
iqtree2 -s ${out}Gubbins/clean.core.aln -o Reference -m MFP+ASC \
	-T AUTO --threads-max $ncores \
	-b 1000 --prefix ${out}IQTREE/Phylogeny

if [ $? -eq 1 ]; then
	echo "IQTREE failed. Please look at the logs to figure out where"
	exit 1
else
	echo "A core SNP phylogeny has been create!"
	exit 0
fi
