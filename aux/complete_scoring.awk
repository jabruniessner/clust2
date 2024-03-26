#!/usr/bin/awk -f

BEGIN{
	#Specifiying the starting max value
	max[2]=1;

	#Initializing initial values for finding max in each step
	for(i=0; i<=2; i++)
	{
		ClSize[i]=1;
	}
}

NR==1{print}
NR==2{printf("%s  increments   avg. ClSize   max. ClSize   min. ClSize\n", $0)}

NF<5{
	next;
}

NR==3{
	complex_number=$1+1;
}


NR>2{
	lines[NR]=$0;
	Distances[NR]=$NF;
	avg[NR]=complex_number/$1;
	Max_Distance=$NF;
	MAX_LINES=NR
	cluster_num = $1;
	
	ClSize[NR]=1;


}

END{
	min[MAX_LINES+1]=complex_number;
	current_min = complex_number;
	Distance_change= Max_Distance - Distances[3];
	
	for(i=3 ; i<=MAX_LINES; i++)
	{
		split(lines[i], a);
		ClSize[-(a[2]+1)]=ClSize[a[3]]+ClSize[a[4]];
		if(ClSize[-(a[2]+1)]>max[i-1])
			max[i]=ClSize[-(a[2]+1)];
		else
			max[i]=max[i-1];


	}



	#Creating the minimums array
	for(i=MAX_LINES; i>2; i--)
	{	
		min[i] = current_min;
		split(lines[i], a)
		SplitCluster1 = a[3];
		SplitCluster2 = a[4];
		min_cluster = (ClSize[a[3]]<ClSize[a[4]]) ? ClSize[a[3]] : ClSize[a[4]];
		
		if (min_cluster<current_min)
			current_min=min_cluster;
	}
	Distances[2]=Distances[3]
	for(i=3; i<=MAX_LINES; i++)
	{	
		increment = (Distances[i]-Distances[i-1])/Distance_change*100;
		printf("%s   %10.5f   %10.5f   %10i   %10i\n", lines[i], increment, avg[i], max[i], min[i]);
	}
}
