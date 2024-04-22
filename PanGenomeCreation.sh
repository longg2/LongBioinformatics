#! /usr/bin/env bash
######################################
### Functions that I'll be calling ###
######################################
script_name=$0
script_full_path=$(dirname $0)

# Making a small one off function to run some awk code in parallel
#modGFF3(){
#	local gff3=$1
#	awk '{if(/^##f|^# /) next; print $0} END{print "##FASTA"}' $gff3
#}

#export -f modGFF3
source $script_full_path/lib/BasicCommands.sh # This loads the basic things I need.
# These are the files and variables that will be needed.
usage() { printf 'Pan-Genome Creation V1
	Uses Prokka and Roary
        -i\tThe folder containing the genomes
	-p\tThe Prokka Output folder
	-r\tThe Roary Output folder
	-l\tThe Log file (Default: current date)
	-t\tSet of Trusted Proteins
        -n\tNumber of CPU Threads to be used
        -h\tShow this help message and exit\n' 1>&2; exit 1; }
log() {	printf "Pan-Genome Creation Settings for $(date):
	Log File:\t${log}
	Input folder:\t${folder}
	Prokka folder:\t${AnnotOut}
	Roary folder:\t${RoaryOut}
	Trusted Proteins:\t${trustedProteins}
	CPU Threads:\t${ncores}
	-------------------------------------\n"; exit 0;
}

#Bakta Database:\t${baktaDB}
#-d\tDatabase for Bakta

# The default variables
log="$(date +'%Y%m%d').log"
while getopts "i:p:r:t:l:d:n:h" arg; do
        case $arg in
                i)
                        folder=${OPTARG}
                        ;;
                l)
                        log=${OPTARG}
                        ;;
		t)
			trustedProteins=${OPTARG}
			;;
                p)
                        AnnotOut=${OPTARG}
                        ;;
                r)
                        RoaryOut=${OPTARG}
                        ;;
                n)
                        ncores=${OPTARG}
                        ;;
#		d)
#			baktaDB=${OPTARG}
#			;;
                h | *)
                        usage
                        exit 0
                        ;;
        esac
done

log | tee $log 
mkdir -p $AnnotOut

# Need to get the correct parallelization cores setup
#parcore=$((ncores/8)) # This will round down
#if [[ $parcore == 0 ]]; then
#	parcore=1
#fi

echo "Annotating the genomes"
parallel -j 1 --bar "prokka {} --outdir $AnnotOut --prefix {/.} --force --proteins $trustedProteins --cpus $ncores --kingdom bacteria --gcode 11 --quiet --compliant" ::: $folder/{*.fasta,*.fa,*.fna}
#parallel -j 1 --bar "bakta {} --output $AnnotOut --prefix {/.} --keep-contig-headers --db $baktaDB --proteins $trustedProteins --threads $ncores" ::: $folder/{*.fasta,*.fa,*.fna}
#
## Need to do some edits to ensure that the Bakta annotation is compattible with Panaroo
#echo "Modifying the GFF3s for Roary"
#mkdir -p ModGFF3
#parallel -j $ncores --bar "cat <(modGFF3 $AnnotOut/{/.}.gff3) $folder/{/.}.fna > ModGFF3/{/}" ::: $AnnotOut/*gff3
#
#echo "Now to run Roary"
roary -p $ncores -e -n -s -cd 95 -i 90 -f $RoaryOut -r $AnnotOut/*.gff

echo "Creating the pan-genome"
#panaroo -t $ncores --merge_paralogs --clean-mode strict -a core -o $RoaryOut -i ModGFF3/*.gff3 --remove-invalid-genes # Removing invalid genes ensures that we're not dealing with pseudogenes in the pan-genome
