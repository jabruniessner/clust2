#!/bin/bash

if [[ ${#@} -lt 2 ]]; then
	echo "Not enough arguments"
	echo "Usage this script cluster_no cluster_out_file"
fi

a="awk -v cluster_no=$1 -f ../aux/tree_post_cutting.awk $2"


clusters=($($a))


for(( i=0 ; i<${#clusters[@]}; i++ )); do
	cluster=${clusters[$i]}	
	members=($(awk -v cluster_no=$cluster -f ../aux/print_members_of_clusters.awk cluster_average_linkage.out))
	
	for member in ${members[@]}; do
		member_string="Complexe $member : cluster $i"
		lines+=("${member_string}")
	done
done

printf "%s \n" "${lines[@]}" | sort -nk2

