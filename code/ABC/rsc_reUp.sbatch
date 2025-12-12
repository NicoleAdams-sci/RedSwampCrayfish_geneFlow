#!/bin/sh
#SBATCH -C [intel18|amd20]
#SBATCH -N 1 -c 1
#SBATCH -t 32:00:00
#SBATCH --mem 2G
#SBATCH -J reup_rsc

root=/mnt/research/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/
maxjob=500

cd $root/SHELL

sq | grep rscABC | cut -f 2 -d _ | cut -f 1 -d ' ' > Qlist_$USER.txt

unset arr_str
for ((i=1; i<=$maxjob; i++)) do
	grep -w $i Qlist_$USER.txt > /dev/null
	if [[ $? -ne 0 ]]
	then
		if [[ ${#arr_str} -eq 0 ]]
		then 
			arr_str="$i"
		else
			arr_str="$arr_str,$i"
		fi	
	fi
done

if [[ ${#arr_str} -gt 0 ]]
then
	sbatch -a $arr_str ABCsims.sbatch
fi

rm Qlist_$USER.txt

sleep 30h

sbatch rsc_reUp.sbatch
