

NR>2{	
	cluster_num = -(NR-2)
	cluster_array[cluster_num]=$0;
}

END{
	
	clusters[cluster_num]=cluster_num
	
	for(i=cluster_num; length(clusters)<cluster_no; i++)
	{
		delete clusters[i]
		split(cluster_array[i], values)
		cluster1 = values[3]
		cluster2 = values[4]
		clusters[cluster1]=cluster1
		clusters[cluster2]=cluster2

	}

	
	for(cluster in clusters)
	{
		print cluster
	}

}
