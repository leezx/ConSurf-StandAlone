#!/bin/bash
source /etc/profile
#$ -S /bin/bash
#$ -pe pvm 2
#$ -cwd
#$ -N consurf

# conda activate consurf
source /home/zz950/softwares/miniconda3/bin/activate /home/zz950/softwares/miniconda3/envs/consurf

work_dir=/home/zz950/tmpData/consurf

for i in `ls db/human/individual`
do
	j="$work_dir/db/human/individual/$i"
	out_file="$work_dir/results/$i"
	in_file="$work_dir/$i"
	#####################################
	if [ -f "$out_file" ]; then
		# echo "$out_file exists."
		continue
	else
		echo "Processing protein: $i, using consurf..."
		cp $j $in_file
		touch "$out_file"
		consurf_dir="$work_dir/results/tmp_$i"
		mkdir $consurf_dir
		python /home/zz950/softwares/ConSurf/stand_alone_consurf.py --algorithm HMMER --Maximum_Likelihood --seq $in_file --dir $consurf_dir &&\
		mv $consurf_dir/consurf_grades.txt $out_file &&\
		rm -r $consurf_dir 
		rm $in_file
		echo "$i done!"
	fi
	# break
done

# find "/home/zz950/tmpData/consurf/results/" -size 0 -type f -delete

