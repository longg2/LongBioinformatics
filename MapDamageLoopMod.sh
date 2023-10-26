#! /usr/bin/env sh
# These are the files and variables that will be needed.
usage() { printf 'MapDamage Loop V1 -- Now includes a Mismatch Search, FLD Extraction, and Parallelization!
	NOTE: MapDamage does not like PE reads, R2 is removed from the MapDamage analysis.
	NOTE: For the true count of cores used, multiply n * 4.
        -i\tThe folder containing the Bam Files
	-o\tThe Output Folder Suffix
	-r\tThe Reference Sequence
	-n\tNumber of cores used (Default = 2)
        -h\tShow this help message and exit\n' 1>&2; exit 1; }
log() {	printf "MapDamage Settings for $(date):
	Log File:\t${log}
	Input folder:\t${in}
	Output Prefix:\t${out}
	Reference:\t${ref}
	CPU Threads:\t${ncores}
	-------------------------------------\n"; exit 0;
}

ncores=2
log="$(date +'%Y%m%d').log"
while getopts "i:l:o:r:n:h" arg; do
        case $arg in
                i)
                        in=${OPTARG}
                        ;;
                l)
                        log=${OPTARG}
                        ;;
                o)
                        Outfolder=${OPTARG}
                        ;;
                r)
                        ref=${OPTARG}
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
mkdir -p MapDamage$out
mkdir -p FLD$out
mkdir -p Mismatches$out

parallel -j $ncores --bar "samtools view -h -F 1 {} > tmp{%}.sam;
			mapDamage --seq-length 25 --merge-reference-sequences -t {/.} -i tmp{%}.sam -r $ref -d MapDamage$out/{/.};
			samtools stats {} | grep '^RL' | cut -f 2- > FLD${out}/{/.}.tab;
			samtools view {} | cut -f 13 | tr -d 'NM:i:' | sort | uniq -c > Mismatches${out}/{/.}.tab" ::: $in/*
rm -f tmp*.sam
