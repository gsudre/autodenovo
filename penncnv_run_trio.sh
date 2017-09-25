#~/bin/bash
# Runs all samples in a pedigree file in the correct order
root_dir='/data/NCR_SBRB/big_fake_simplex'
out_dir="${root_dir}/penncnv"
trio=$1
ped_fname="${root_dir}/${trio}.ped"
nlines=`cat ${ped_fname} | wc -l`

if [ $nlines != 3 ]; then
	echo Pedigree file for trio needs to have exactly 3 lines! Exitting...
	return -1;
fi

cd $out_dir

while read line; do
        IFS="	" read -r -a array <<< "$line";
	id=${array[1]}
	fa=${array[2]}
	mo=${array[3]}
	sex=${array[4]}

	# figure out which box has the sample
        spath=`ls fs_c*/* | grep ${id}`;
        IFS='/' read -r -a array <<< "$spath";
        box=${array[0]};

	# figure out relationship
        if [ $fa != 0 ]; then
		child=${box}/${id};
	elif [ $sex == 1 ]; then
		father=${box}/${id};
	else
		mother=${box}/${id};
	fi;
done < $ped_fname

cat $ped_fname
echo ''
echo Found father=$father mother=$mother child=$child

module load penncnv

cd $out_dir
detect_cnv.pl -test -hmm example.hmm -pfb filtered_MAF.pfb -log ${trio}.log $father $mother $child -out ${trio}.rawcnv;
detect_cnv.pl -test -hmm example.hmm -pfb filtered_MAF.pfb -log ${trio}.log $father $mother $child -out ${trio}.adjusted.rawcnv -gcmodel filtered.h19.gcmodel;
detect_cnv.pl -trio -hmm example.hmm -pfb filtered_MAF.pfb -cnv ${trio}.rawcnv $father $mother $child -out ${trio}.triocnv;
detect_cnv.pl -trio -hmm example.hmm -pfb filtered_MAF.pfb -cnv ${trio}.adjusted.rawcnv $father $mother $child -out ${trio}.adjusted.triocnv;

