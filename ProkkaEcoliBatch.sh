#!/usr/bin/env bash

###########################
GenomeFolder=$1
ProkkaFolder=$2
RoaryFolder=$3
ncores=$4
##########################
mkdir -p $ProkkaFolder
#mkdir -p $RoaryFolder
mkdir -p looptmp
mkdir -p Logs

ls -1 $GenomeFolder/*.f* > list.tmp
echo "Note that there will be no GBK files as tbl2asn is too old...."

#total=$(cat list.tmp | wc -l)
#count=1
#while IFS= read -r genome; do
#	name=$(basename $genome .fna) # In case it's fna
#	name=$(basename $name .fasta)
#
#	# In case it's coming from the Masked Consensuses (they all have the same file name...)
#	echo $name
#	if echo "$name" | grep -Eq "MaskedConsensus.*"; then
#		name=$(echo $genome |cut -f 6 -d "/")
#	fi
#	prokka --outdir looptmp --force --proteins ~/EColiProject/EcoliGeneBank/K12.gbk --kingdom bacteria --gcode 11 --cpus $ncores --compliant $genome 2> Logs/$name.log
#	mv looptmp/*.gff $ProkkaFolder/$name.gff
#	echo "Strain $count out of $total Complete"
#	count=$(echo "$count + 1"|bc)
#done < list.tmp

#rm -rf looptmp* list.tmp

# Finally, we want to run roary on the results
#roary -p $ncores -f $RoaryFolder -r $ProkkaFolder/*.gff
roary -p $ncores -s -e -n -cd 95 -i 90 -f $RoaryFolder -r $ProkkaFolder/*.gff
