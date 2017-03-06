cd data/broadsims/00265 ##Chosen at random
leaves=($(cat original_g_trees.tree | sed -e "s/:[^),]*//g" -e "s/)[0-9.]*//g" -e "s/[(,);]/ /g" -e "s/ \+/ /g" | head -n 1 | sed -e 's/ /\n/g' | sort | paste -sd " "))
#echo ${leaves[@]}
for file in *.tree
do
	name=$(echo $file | sed "s/.tree//g")
	echo "x,y,press" > ${name}_matrix.csv
	n_line=0
	while read line
	do
		this_leaves=($(echo $line | sed -e "s/:[^),]*//g" -e "s/)[0-9.]*//g" -e "s/[(,);]/ /g" -e "s/ \+/ /g" -e 's/ /\n/g' | sort | paste -sd " "))
#		echo ${this_leaves[@]}
		nleaf=0
		nthisleaf=0
		for rleaf in "${leaves[@]}"
		do
#			echo "Comparing $rleaf with ${this_leaves[$nthisleaf]}"
			if [[ "$rleaf" == "${this_leaves[$nthisleaf]}" ]]
			then
				echo "$n_line,$nleaf,1" >> ${name}_matrix.csv
				nthisleaf=$(( $nthisleaf + 1 ))
#				echo "equal"
			else
				echo "$n_line,$nleaf,0" >> ${name}_matrix.csv
#				echo "different"
			fi
			nleaf=$(( $nleaf + 1 ))
		done
		n_line=$(( $n_line + 1 ))
	done < $file

done
