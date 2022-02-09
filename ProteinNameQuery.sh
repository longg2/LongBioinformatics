#! /usr/bin/env bash
file=$1
out=$2

IFS=$'\n'

for LINE in $(cat $file); do
#while IFS=$'\n' read LINE; do
	esearch -db protein -query "$LINE" | esummary |\
		xtract -pattern DocumentSummary -element Id,Caption,Title |\
		tee -a $out
	sleep 3
done
#done < $file
