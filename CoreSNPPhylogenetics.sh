#! /usr/bin/env bash
# TODO:

script_name=$0
script_full_path=$(dirname $0)

source $script_full_path/lib/BasicCommands.sh # This loads the basic things I need.

usage() { printf 'Core SNP phylogeny creation.
	Given a folder of sequences, it will create a core snp phylogeny using:
		• Snippy
		• Gubbins, and
		• IQTREE-2
	Phylogenetic model will be automated 1000 bootstraps used.
	NOTE: Doesnt handle PE samples for now. Assumes merged,
	assemblies, or full genomes

        -i\tWhere the sequences are
	-o\tOutput Prefix (Default = Phylogenetics)
	-r\tReference Sequence
	-n\tNumber of CPU Threads to be used (Default = 8)
	-l\tLog File Name
        -h\tShow this help message and exit\n' 1>&2; exit 0; }
log() {	printf "Phylogeny settings for $(date):
	Log File:\t${log}
	Input folder:\t${folder}
	Output folder:\t${out}
	Reference:\t${db}
	CPU Threads:\t${ncores}
	-------------------------------------\n"; exit 0;
}

log="$(date +'%Y%m%d').log"
declare -i ncores=8
export out="Phylogenetics"
while getopts "i:r:n:o:l:h" arg; do
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
                h | *)
                        usage
                        exit 0
                        ;;
        esac
done

if [ -z ${folder+x} ] || [ -z ${reference+x} ]; then
	echo "You are missing either the Input folder or the Reference Genome"
	usage
	exit 1
fi

log | tee $log
printf "\n"

# Making the Directories
mkdir -p ${out}Snippy
mkdir -p ${out}Gubbins
mkdir -p ${out}IQTREE

# Activating conda environments so that we can switch when needed
source ~/miniconda3/etc/profile.d/conda.sh # Activating 
conda activate pangenome

############# Running SNIPPY ###############
# First, we need to make a manifest file for snippy-multi
find ~+/$folder $script_full_path/lib/snippyManifest.txt > ${out}snippyManifest.txt

# Running Snippy. Because of it functions, we'll be piping the log to nowhere
echo "Running snippy-multi"

cd ${out}Snippy
tmp=$(echo $OLDPWD) # So that I can work without issue...
snippy-multi $tmp/${out}snippyManifest.txt --ref $reference --cpus $ncores \
	--mapqual 30 --basequal 20 | bash 2> /dev/null

if [ $? -eq 1 ]; then
	echo "Snippy failed. Please look at the logs to figure out where"
	exit 1
fi

snippy-clean_full_aln core.full.aln > clean.full.aln

cd $tmp
unset tmp

########## Time for Gubbins #############
conda activate phylogenies # Need to activate

# Running Gubbins is relatively simple, if potentially long...
run_gubbins.py ${out}Snippy/clean.full.aln --outgroup Reference --threads $ncores \
       	--prefix ${out}Gubbins/RecombMask
snp-sites -c ${out}Gubbins/RecombMask.filtered_polymorphic_sites.fasta > ${out}Gubbins/clean.core.aln

############### IQ TREE 2 ################
iqtree2 -s ${out}Gubbins/clean.core.aln -o Reference -m MFP+ASC \
	-T AUTO --threads-max $ncores \
	-b 1000 --prefix ${out}IQTREE/Phylogeny
