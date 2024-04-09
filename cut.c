/* The clust2 program
 * Copyright (C) 2024 Jakob Niessner.
 * Contact: jabruniessner@gmail.com
 *
 *
 * This program was written at the Heidelberg Institute for theoretical studies,
 * Schloß-Wolfsbrunnenweg 35, 69118 Heidelberg, Germany
 * Under the supervision of Prof. Dr. Rebecca C. Wade
 * Contact: rebecca.wade@h-its.org
 *
 * Permission to use, copy, modify, and distribute this software and its
 * documentation with or without modifications and for any purpose and
 * without fee is hereby granted, provided that any copyright notices
 * appear in all copies and that both those copyright notices and this
 * permission notice appear in supporting documentation, and that the
 * names of the contributors or copyright holders not be used in
 * advertising or publicity pertaining to distribution of the software
 * without specific prior permission.
 *
 * THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
 * WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
 * CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
 * OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
 * OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
 * OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
 * OR PERFORMANCE OF THIS SOFTWARE.
 *
 */

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<sys/queue.h>
#include<math.h>
#include"cluster.h"
#include"pdb_structure.h"
#include"trajectories.h"
#include"math_func.h"
#include"shared_functions.h"


typedef struct tree_file_entry{
	int Node;
	int CycleNo;
	int Item1;
	int Item2;
	double Distance; 
}tree_file_entry;

typedef struct tree_list_entry{
	tree_file_entry tree_entry;
	SLIST_ENTRY(tree_list_entry) entry;
}tree_list_entry;

//Defining struct for list head
SLIST_HEAD(Tree_list, tree_list_entry);


Node* read_tree_file(char* filename, int* node_number)
{
	FILE* tree_file = fopen(filename, "r");
	const char* delimiter=" ";
	size_t buffer_size = 0;
	char* buffer = NULL;
	ssize_t read;
	
	if((read = getline(&buffer, &buffer_size, tree_file)) == -1 )
	{
		printf("There is something wrong with your tree file.");
		exit(EXIT_FAILURE);
	}

	
	if((read = getline(&buffer, &buffer_size, tree_file)) == -1 )
	{
		printf("There is something wrong with your trajectories file.");
		exit(EXIT_FAILURE);
	}

	struct Tree_list tree_list;
	SLIST_INIT(&tree_list);

	
	int node_num=0;
	while((read=getline(&buffer, &buffer_size, tree_file))>=41)
	{	

		
		tree_list_entry* list_entr = malloc(sizeof(tree_list_entry));
	        tree_file_entry* file_entry = &(list_entr->tree_entry);	
		sscanf(buffer, "%4d %7d %8d %8d   %7lf",
				&(file_entry->Node),
				&(file_entry->CycleNo),
				&(file_entry->Item1),
				&(file_entry->Item2),
				&(file_entry->Distance)
		      );
		SLIST_INSERT_HEAD(&tree_list, list_entr, entry);
		node_num++;

	}

	*node_number = node_num+1;


	tree_file_entry* nodes=malloc(node_num*sizeof(tree_file_entry));


	tree_list_entry* list_entr=NULL;
	int j = 1;
	SLIST_FOREACH(list_entr, &tree_list, entry)
	{
		nodes[node_num-j]=list_entr->tree_entry;
		j++;
	}


	while(!SLIST_EMPTY(&tree_list))
	{
		tree_list_entry *list_entry = SLIST_FIRST(&tree_list);
		SLIST_REMOVE_HEAD(&tree_list, entry);
		free(list_entry);
	}

	Node* node_array = (Node*) malloc(node_num*sizeof(Node));
	for(int i=0; i<node_num; i++)
	{
		node_array[i].left=nodes[i].Item1;
		node_array[i].right=nodes[i].Item2;
		node_array[i].distance=nodes[i].Distance;
	}
	

	fclose(tree_file);

	return node_array;


}


int main(int argc, char* argv[])
{	

	if (argc < 5 || !strcmp(argv[1], "-h") || !strcmp(argv[1], "--help"))
	{
		printf("Usage:\n");
		printf("\n");
		printf("%s tree_file complexes_file pdb_file num_of_cluster\n", argv[0]);
		printf("\n");
		printf("tree_file: This is the tree file as was output by the clust program\n");
		printf("\n");
		printf("complexes_file: This is the complexes file you would like to analyze.\n");
		printf("It needs to be in the ascii format from SDA (see https://mcm.h-its.org/sda/doc/doc_sda7/complexes_file.html)\n");
		printf("\n");
		printf("pdb_file: This is the pdb file of the ligand. i.e. the molecule that was moving during the BD simulation.\n");
		printf("The RMSD will be measured via the backbone, which are all atoms that are labelled C, N or O in the PDB file\n");
		printf("\n");
		printf("num_of_cluster: The cluster representatives to be printed into the summarized complexes file\n");
		printf("\n");
		exit(EXIT_FAILURE);
	}
      	printf("The number of input arguments is %i\n", argc);
	int cluster_num;
	sscanf(argv[4], "%i", &cluster_num);
	int num_energy_fields;
	char* tree_filename=argv[1];
	char* filename=argv[2];
	double pdb1_com[3];
	double pdb2_com[3];
	Record* record;
	char**head_lines;
	char** Energy_fields;
	int complexes = read_trajectory_file(filename, &head_lines, &record, &num_energy_fields, &Energy_fields, pdb1_com, pdb2_com);
	int complexes_input;

	printf("The number of energy fields is: %i\n", num_energy_fields);


	if(complexes<cluster_num){
		printf("Error: you chose a larger number of representatives than complexes\n");
		exit(EXIT_FAILURE);
	}

	printf("The center of mass of the first protein is: %lf\t%lf\t%lf\n", pdb1_com[0], pdb1_com[1], pdb1_com[2]);
	printf("The center of mass of the second protein is: %lf\t%lf\t%lf\n", pdb2_com[0], pdb2_com[1], pdb2_com[2]);
	printf("The number of complexes %i\n", complexes);

	pdb_struct PDB;
	char* filename2  = argv[3];
	int line_num = read_pdb_file(&PDB, filename2);
	double** structures=prepare_data(PDB.bb_atom_number, complexes, record, PDB.backbone_coordinates, pdb2_com);
	

	int** mask=prepare_mask(PDB.bb_atom_number, complexes);

	Node* Tree = read_tree_file(tree_filename, &complexes);

	
	if(Tree==NULL)
	{
		printf("Could not read Tree files");
		exit(EXIT_FAILURE);
	}

	int* clusterid = cutting_hierarchical_tree(Tree, complexes, cluster_num);
	
	
	double **cdata=malloc(sizeof(double*)*cluster_num);
	int **cmask=malloc(sizeof(double*)*cluster_num);
	
	for(int i = 0; i<cluster_num; i++)
	{
		cmask[i]=malloc(sizeof(int)*3*PDB.bb_atom_number);
		cdata[i]=malloc(sizeof(double)*3*PDB.bb_atom_number);
	}
	
	

	int num = getclustercentroids(cluster_num, complexes, 3*PDB.bb_atom_number, structures, mask, clusterid,
		cdata, cmask, 0, 'a');


	
	int* representatives;
	long* ClSize = malloc(sizeof(long)*cluster_num);
	long* ClFSize = malloc(sizeof(long)*cluster_num);
	double* AvgEnergies = malloc(sizeof(double)*cluster_num);
	double* StdEnergies = malloc(sizeof(double)*cluster_num);
	double* RepRMSD = malloc(sizeof(double)*cluster_num);
	double* CLFRMSD = malloc(sizeof(double)*cluster_num);
	double* spread = malloc(sizeof(double)*cluster_num);
	double* stddev = malloc(sizeof(double)*cluster_num);
	double* max = malloc(sizeof(double)*cluster_num);


	representatives = do_scoring(cluster_num, complexes, 3*PDB.bb_atom_number,
			structures, PDB.backbone_coordinates, pdb1_com ,clusterid, cdata, record, ClSize, ClFSize, AvgEnergies, StdEnergies, RepRMSD, CLFRMSD, spread, stddev, max);



	double**  energy_data = get_energy_data_set(complexes, num_energy_fields, record);
	double** energy_cdata = malloc(sizeof(double*)*cluster_num);
	double** energy_std_cdata = malloc(sizeof(double*)*cluster_num);
	int **   energy_cmask = malloc(sizeof(double*)*cluster_num);

	for(int i = 0; i<cluster_num; i++)
	{
		energy_cmask[i]=malloc(sizeof(int)*num_energy_fields);
		energy_cdata[i]=malloc(sizeof(double)*num_energy_fields);
		energy_std_cdata[i]=malloc(sizeof(double)*num_energy_fields);
	}

	int** energy_mask = prepare_mask(num_energy_fields, complexes);

	int num_energies = getclustercentroids(cluster_num, complexes, num_energy_fields, energy_data, mask, clusterid, energy_cdata, energy_cmask, 0, 'a' );

	get_stds_of_clusters(cluster_num, complexes, num_energy_fields, clusterid, energy_data, energy_cdata, energy_std_cdata);


	
	printf("%i\t%i\t%i\t%i\t%i\n", representatives[0]+1, representatives[1]+1, representatives[2]+1, representatives[3]+1, representatives[4]+1);
	printf("\n");
	printf("\n");
	printf("\n");


	//Computing the distance matrix for the representatives and initializing the cluster structs
	Cluster* cluster= malloc(sizeof(Cluster)*cluster_num);	
	double** representative_data=malloc(sizeof(double*)*cluster_num);
	for(int i = 0; i<cluster_num; i++)
	{	
		int rep = representatives[i];
		Record* rep_record=&record[rep];
		initializing_Cluster(&cluster[i], rep_record, i+1, rep, ClSize[i], ClFSize[i], AvgEnergies[i], StdEnergies[i],
				RepRMSD[i], CLFRMSD[i], spread[i], stddev[i], max[i], num_energy_fields, Energy_fields);
		//representative_data[i] = structures[rep];
	}

	qsort(cluster, cluster_num, sizeof(Cluster), cmp_cluster_sizes);

	print_table_header();

	
	for(int i=0; i<cluster_num; i++)
	{	
		cluster[i].no=i+1;
		print_Cluster(&cluster[i]);
	}

	printf("\n");
	printf("\n");
	printf("\n");

	
	
	//getting the correctly ordered representatives
	for(int i=0; i<cluster_num; i++)
	{	
		int rep = cluster[i].Repr-1;
		representative_data[i] = structures[rep];
	}



	double** matrix = example_distance_bb(cluster_num, 3*PDB.bb_atom_number, representative_data,  mask);

#ifdef DEBUG
	//Making a test that we understand the function correctly
	int** mask=prepare_mask(2, 2);
	double** test_data = malloc(sizeof(double*)*2);
	test_data[0] = malloc(sizeof(double)*2);
	test_data[1] = malloc(sizeof(double)*2);
	test_data[0][0] = 1.;
	test_data[0][1] = 1.;
	test_data[1][0] = 2.;
	test_data[1][1] = 2.;

	

	printf("This is a test that I understand the euclidean distance correctly!\n");
	example_distance_bb(2, 2, test_data,  mask);
#endif	
	



	print_results_description();

	

	FILE* outfile = fopen("cluster_out.ascii", "w");

	for(int i=0; i<4; i++)
		fprintf(outfile, "%s", head_lines[i]);

	//printing the head lines
	for(int i=0; i<cluster_num; i++)
	{	
		int rep = cluster[i].Repr-1;
		fprintf(outfile,"%s", record[rep].line);
	}

	fclose(outfile);


	return 0;

	
}
