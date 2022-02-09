#! /usr/bin/env bash
######################################
### Functions that I'll be calling ###
######################################
script_name=$0
script_full_path=$(dirname $0)

source $script_full_path/lib/BasicCommands.sh # This loads the basic things I need.
# These are the files and variables that will be needed.
usage() { printf 'Pan-Genome Creation V0.75
	Now with Prokka parallelized and not hard coded for E. coli
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
	Prokka folder:\t${ProkkaOut}
	Roary folder:\t${RoaryOut}
	Trusted Proteins:\t${trustedProteins}
	CPU Threads:\t${ncores}
	-------------------------------------\n"; exit 0;
}

# The default variables
log="$(date +'%Y%m%d').log"
while getopts "i:p:r:t:l:n:h" arg; do
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
                        ProkkaOut=${OPTARG}
                        ;;
                r)
                        RoaryOut=${OPTARG}
                        ;;
                n)
                        ncores=${OPTARG}
                        ;;
                h | *)
                        usage
                        exit 0
                        ;;
        esac
done

log | tee $log 
mkdir -p $ProkkaOut

# Need to get the correct parallelization cores setup
parcore=$((ncores/8)) # This will round down
if [[ $parcore == 0 ]]; then
	parcore=1
fi

echo "Annotating the genomes"
parallel -j 1 --bar "prokka {} --outdir $ProkkaOut --prefix {/.} --force --proteins $trustedProteins --cpus $ncores --kingdom bacteria --gcode 11 --quiet --compliant" ::: $folder/{*.fasta,*.fa,*.fna}

echo "Now to run Roary"
roary -p $ncores -e -n -s -cd 95 -i 90 -f $RoaryOut -r $ProkkaOut/*.gff
