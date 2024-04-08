#!/bin/awk -f
#
#Usage: this_script complexes_file clusters_file cluster_no
#

function get_members(num,b,c,	num1, num2, line_fields)
{	
	if(num>=0)
	{
		c[num]=num;
		return;
	}


	string = b[num];
	split(string, line_fields);

	num1 = line_fields[3];
	num2 = line_fields[4];
	
	get_members(num1, b, c);
	get_members(num2, b, c);
	return;

}


BEGIN{
	i=0;
	j=0;

	if(cluster_no >=0)
	{
		print cluster_no;
		exit
	}
      }

#NR==FNR && FNR>4{a[i++]=$0; next;}

#FNR<=4 && FNR==NR{print}

NF==0{next;}

FNR>2&&(-$2-1>=cluster_no){
	b[--j]=$0;
}

END{	
	split(b[j], c);
	start_num1=c[3];
	start_num2=c[4];
	get_members(start_num1, b, members)
	get_members(start_num2, b, members)
	for(member in members)
	{
		print member;
	}
}

