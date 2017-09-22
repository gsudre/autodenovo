#!/bin/bash
# Spits out the chip for each person in the trio, to make sure they're the same
ped=$1

cut -f 2 $ped > ids.txt
while read s; do
	spath=`ls fs_c*/* | grep ${s}`;
	IFS='/' read -r -a array <<< "$spath";
	box=${array[0]};
	content=`grep Content ${box}_FinalReport.txt`;
	echo $s in $box = $content; 
done < ids.txt

rm ids.txt 
