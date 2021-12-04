#! /usr/bin/env sh
# These are the files and variables that will be needed.
usage() { printf 'MapDamage Loop V0.8 -- Now includes a Mismatch Search, FLD Extraction, and Parallelization!
	NOTE: MapDamage does not like PE reads, R2 is removed from the MapDamage analysis.
	NOTE: For the true count of cores used, multiply n * 4.
        -i\tThe folder containing the Bam Files
	-o\tThe Output Folder Suffix
	-r\tThe Reference Sequence
	-n\tNumber of cores used 
        -h\tShow this help message and exit\n' 1>&2; exit 1; }

while getopts "i:o:r:n:h" arg; do
        case $arg in
                i)
                        in=${OPTARG}
                        echo "The raw sequencing files are located in $in"
                        ;;
                o)
                        Outfolder=${OPTARG}
                        echo "Output folder suffix is $Outfolder"
                        ;;
                r)
                        ref=${OPTARG}
                        echo "The reference sequence is $ref"
                        ;;
                n)
                        ncores=${OPTARG}
                        echo "Using $ncores threads for the analysis"
                        ;;
                h | *)
                        usage
                        exit 0
                        ;;
        esac
done
mkdir -p MapDamage$Outfolder
mkdir -p FLD$Outfolder
mkdir -p Mismatches$Outfolder

parallel -j $ncores --bar "samtools view -h -F 1 {} > tmp{%}.sam;
			mapDamage --merge-reference-sequences -t {/.} -i tmp{%}.sam -r $ref -d MapDamage$Outfolder/{/.};
			samtools stats {} | grep '^RL' | cut -f 2- > FLD${Outfolder}/{/.}.tab;
			samtools view {} | cut -f 13 | tr -d 'NM:i:' | sort | uniq -c > Mismatches${Outfolder}/{/.}.tab" ::: $in/*

rm -f tmp*.sam


#for file in $in/*; do
#	sample=$(basename $file .bam)
#	samtools view -h -F 1 $file > tmp.bam # The program can't deal with paired end reads....
#	mapDamage --merge-reference-sequences -i tmp.bam -r $ref -d MapDamage$Outfolder/$sample # This is the Main one.  Where I plot the
#
#	# This is the FLD extraction
#	samtools stats $file | grep "^RL" | cut -f 2- > FLD${Outfolder}/$sample.tab
#
#	# Getting Mismatches
#	samtools view $file | cut -f 13 | tr -d "NM:i:" | sort | uniq -c > Mismatches${Outfolder}/$sample.tab
#done

rm -f tmp.bam
