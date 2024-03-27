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

typedef struct{
	int no;
	long ClSize;
	long ClFSize;
	int Repr;
	double ReprE; 
	double ClAE;
	double ClAED;
	double RepRMSD;
	double ClFRMSD; 
	double Ele;
	double ElDesE;
	double HyDesE;
	double LjE;
	double spread;
	double stddev; 
	double max;
} Cluster;


int cmp_cluster_sizes(const void* clust1, const void* clust2)
{
	Cluster* Cluster1 = (Cluster*) clust1;
	Cluster* Cluster2 = (Cluster*) clust2;
	long ClFSize1 = Cluster1->ClFSize;
	long ClFSize2 = Cluster2->ClFSize;

	return (ClFSize1 < ClFSize2) - (ClFSize1 > ClFSize2);
}



double** example_distance_bb(int nrows, int ncols, double** data, int** mask)
/* Calculate the distance matrix between genes using the Euclidean distance. */
{ int i, j;
  double* weight = malloc(ncols*sizeof(double));
  /* Set up the ragged array */
  double** matrix = malloc(nrows*sizeof(double*));
  if (matrix) {
    matrix[0] = NULL;
    for (i = 1; i < nrows; i++)
    { matrix[i] = malloc(i*sizeof(double));
      if (matrix[i]==NULL) /* Not enough memory available */
      { while (--i >= 0) free(matrix[i]);
        free(matrix);
        matrix = NULL;
        break;
      }
    }
  }
  if (!matrix) {
    free(weight);
    printf("Insufficient memory to store the distance matrix\n");
    return NULL;
  }

  printf("============ Backbone distance matrix between representatives ============\n");
  for (i = 0; i < ncols; i++) weight[i] = 1.0;
  distancematrix(nrows, ncols, data, mask, weight, 'e', 0, matrix);
  printf("   Representative:");
  for(i=0; i<nrows-1; i++) printf("%6d", i+1);
  printf("\n");
  for(i=0; i<nrows; i++)
  { printf("Representative %2d:",i+1);
    for(j=0; j<i; j++) printf(" %5.2f", sqrt(3*matrix[i][j]));
    printf("\n");
  }

  printf("\n");
  free(weight);

  
  return matrix;
  
}


Node* hierarchical_single_linkage(int nrows, int ncols, double** data, int** mask)
{
	const int nnodes = nrows-1;
	double* weight = malloc(ncols*sizeof(double));
	Node* tree;
	FILE* outfile = fopen("cluster_single_linkage.out", "w");
	for(int i=0; i<ncols; i++) weight[i]=1.0;
	fprintf(outfile, "\n");
  	fprintf(outfile,"#================ Pairwise single linkage clustering ============\n");
	/*Since we have the distance matrix herem we may as well use it. */
	tree = treecluster(nrows, ncols, data, mask, weight ,0, 'e', 's', 0);
	/*The distance matrix was modified by treecluster, so we cannot use it any
	 * more. But we still need to deallocate it here.
	 * The first row of distmatrix is a single null pointer; no need to free it.
	 */
	if(!tree)
	{ /*Indication that the treecluster routine failed*/
	  	printf("treecluster routine failed due to insufficient memory\n");
		free(weight);
		exit(EXIT_FAILURE);
	}
	fprintf(outfile,"#Node CycleNo    Item1    Item2   Distance\n");
	
	for(int i=0; i<nnodes; i++)
	{
		fprintf(outfile,"%4d %7d %8d %8d   %7.2lf\n",
			nnodes-i, i, tree[i].left, tree[i].right, sqrt(3*tree[i].distance));
	}
	
	fprintf(outfile, "\n");
	free(weight);
	fclose(outfile);
	return tree;
 
}

Node* hierarchical_maximum_linkage(int nrows, int ncols, double** data, int** mask)
{
	const int nnodes = nrows-1;
	double* weight = malloc(ncols*sizeof(double));
	Node* tree;
	FILE* outfile = fopen("cluster_maximum_linkage.out", "w");
	for(int i=0; i< ncols; i++) weight[i] = 1.0;
  	fprintf(outfile, "#================ Pairwise maximum linkage clustering ============\n");
	tree = treecluster(nrows, ncols, data, mask, weight, 0, 'e', 'm', 0);
  	/* Here, we let treecluster calculate the distance matrix for us. In that
   	* case, the treecluster routine may fail due to insufficient memory to store
   	* the distance matrix. For the small data sets in this example, that is
   	* unlikely to occur though. Let's check for it anyway:
   	*/
	if(!tree)
	{ /*Indication that the treecluster routine failed*/
    	  printf ("treecluster routine failed due to insufficient memory\n");
	  free(weight);
	  exit(EXIT_FAILURE);
	}

	fprintf(outfile,"#Node CycleNo    Item1    Item2   Distance\n");
        for(int i=0; i<nnodes; i++)
        {
                fprintf(outfile,"%4d %7d %8d %8d   %7.2lf\n",
                        nnodes-i, i, tree[i].left, tree[i].right, sqrt(3*tree[i].distance));
        }

	fprintf(outfile,"\n");
	free(weight);
	fclose(outfile);
	return tree;
}


Node* hierarchical_average_linkage(int nrows, int ncols, double** data, int** mask)
{
	const int nnodes = nrows-1;
	double* weight = malloc(ncols*sizeof(double));
	Node* tree;
	FILE* outfile = fopen("cluster_average_linkage.out", "w");
	for(int i=0; i< ncols; i++) weight[i]=1.0;
  	fprintf(outfile, "#================ Pairwise average linkage clustering ============\n");
	tree = treecluster(nrows, ncols, data, mask, weight, 0, 'e', 'a', 0);
  	/* Here, we let treecluster calculate the distance matrix for us. In that
   	* case, the treecluster routine may fail due to insufficient memory to store
   	* the distance matrix. For the small data sets in this example, that is
   	* unlikely to occur though. Let's check for it anyway:
   	*/
	if(!tree)
	{ /* Indication that the treecluster routine failed */
	  printf("treecluster routine failed due to insufficient memory\n");
	  exit(EXIT_FAILURE);
	}
	
	fprintf(outfile,"#Node CycleNo    Item1    Item2   Distance\n");
        for(int i=0; i<nnodes; i++)
        {
                fprintf(outfile,"%4d %7d %8d %8d   %7.2lf\n",
                        nnodes-i, i, tree[i].left, tree[i].right, sqrt(3*tree[i].distance));
        }
	
	fprintf(outfile,"\n");
	free(weight);
	fclose(outfile);
	return tree;
}

Node* hierarchical_centroid_linkage(int nrows, int ncols, double** data, int** mask)
{
	const int nnodes = nrows-1;
	double* weight = malloc(ncols*sizeof(double));
	FILE* outfile = fopen("cluster_centroid_linkage.out", "w");
	Node* tree;
	for(int i = 0; i<ncols; i++) weight[i]=1.0;
  	fprintf(outfile, "#================ Pairwise centroid linkage clustering ===========\n");
  	tree = treecluster(nrows, ncols, data, mask, weight, 0, 'e', 'c', 0); 
	if(!tree)
	{ /* Indication that the treecluster routine failed */
	  printf( "treecluster routine failed due to insufficient memory\n");
	  free(weight);
	  exit(EXIT_FAILURE);
	}

	fprintf(outfile,"#Node CycleNo    Item1    Item2   Distance\n");
        for(int i=0; i<nnodes; i++)
        {
                fprintf(outfile,"%4d %7d %8d %8d   %7.2lf\n",
                        nnodes-i, i, tree[i].left, tree[i].right, sqrt(3*tree[i].distance));
        }

	fprintf(outfile, "\n");
	free(weight);
	fclose(outfile);
	return tree;	
}

int* cutting_hierarchical_tree(Node* tree, int nrows, int level)
{	
	FILE* outfile = fopen("Treecuting.out", "w");
	int* clusterid = malloc(nrows*sizeof(int));
	int ok = cuttree(nrows, tree, level, clusterid);
	if(!ok)
	{
		printf ("cuttree routine failed due to insufficient memory\n");
		exit(EXIT_FAILURE);
	}
	for(int i=0; i<nrows; i++) fprintf(outfile,"Complexe %2d: cluster %2d\n", i, clusterid[i]);
	fprintf(outfile, "\n");
	fclose(outfile);
	return clusterid;
}



double**prepare_data(int number_bb_atom, int number_samples, Record* complexes, double** bb_atoms, double* pdb_com)
{
	double** structures = (double**) malloc(sizeof(double*)*number_samples);
	for(int j = 0; j< number_samples; j++)
	{
		structures[j]=(double*) malloc(sizeof(double)*3*number_bb_atom);
		for(int k = 0; k<number_bb_atom; k++)
		{	
			for(int n = 0; n<3; n++)
				structures[j][3*k+n]=bb_atoms[k][n];
			

			double* position = &(structures[j][3*k]);
			double multiplied[3];
			subtract_vectors(position, pdb_com, position);
			mat_vec_mul(complexes[j].Rot_mat, 
					position, 
					multiplied);

			add_vectors(complexes[j].Trans, multiplied, position);
		}

	}

	return structures;
}

double compute_distance_pdb_structure(int number_bb_atoms, double**bb_atoms, double* structure, double* pdb_com)
{
	double dist=0;
	for(int i=0; i< number_bb_atoms; i++)
		for(int j=0; j<3; j++)
			dist+=pow(structure[3*i+j]-(bb_atoms[i][j]-pdb_com[j]), 2);

	dist = dist/(double) number_bb_atoms;

	return sqrt(dist);
}



int**prepare_mask(int number_bb_atom, int number_samples)
{	
	int**mask=malloc(sizeof(int*)*number_samples);
	for(int j = 0; j< number_samples; j++)
	{	
        	mask[j]=(int*)malloc(sizeof(int)*3*number_bb_atom);
	for(int k = 0; k<number_bb_atom; k++)
			for(int n= 0; n<3; n++)
				mask[j][3*k+n]=1;
	}
	return mask;
}

int* do_scoring(int nclusters, int nrows, int ncolumns, 
		double** data,  double**orig_atoms, double* pdb_com ,int clusterid[], double** cdata, Record* records, long* ClSize, 
		long* ClFSize, double* AvgEnergies, double* StdEnergies, double* RepRMSD, double* CLFRMSD,
		double* spread, double* stddev, double* max)
{
	int* representatives = malloc(sizeof(int)*nclusters);
	double* currentdist = malloc(sizeof(double)*nclusters);

	
	//Setting all outputvariables initially to zero
	for(int i = 0; i<nclusters; i++)
	{
		AvgEnergies[i]=0;
		StdEnergies[i]=0;
		ClSize[i]=0;
		ClFSize[i]=0;
		CLFRMSD[i]=0;
		RepRMSD[i]=0;
		max[i]=-1;
		spread[i]=0.;
		stddev[i]=0.;
	}
	

	
	//Finding the the representatives the ClSizes and the ClFSizes
	//and summing up energies for AvgEnergies and StdEnergies
	for(int i = 0; i<nclusters; i++)
	{
		representatives[i]=-1;
		currentdist[i]=-1.;
	}

	for(int i = 0; i<nrows; i++)
	{
		int cluster = clusterid[i];
		double* centroid_coordinates = cdata[cluster];
		double* coordinates = data[i];
		double dist = compute_RMSD(coordinates, centroid_coordinates, ncolumns);
		ClSize[cluster]+=1;
		ClFSize[cluster]+=records[i].occurences;
		CLFRMSD[cluster] += compute_distance_pdb_structure((ncolumns/3), orig_atoms, coordinates, pdb_com)*records[i].occurences;
		AvgEnergies[cluster] += records[i].energy_fields[0]*records[i].occurences;


		if((currentdist[cluster]==-1.) || (currentdist[cluster] > dist))
		{
			currentdist[cluster]= dist;
			representatives[cluster] = i;
			
		}
#ifdef DEBUG
		printf("The current iteration is %i\n", i);
		printf("The clusterid is: %i\n", clusterid[i]);
		printf("the current_dist is: %lf\n", currentdist[cluster]);
		printf("The distance is %lf\n", dist);
		printf("\n");
#endif
	}




	//Dividing by the number of occurences to get the Averages
	for(int i = 0; i<nclusters; i++)
	{	
		AvgEnergies[i]/=(double) ClFSize[i];
		CLFRMSD[i]/=(double) ClFSize[i];
	}


	//The rest is about computing the standard deviations
	for(int i = 0; i<nrows; i++)
	{
		int cluster = clusterid[i];
		StdEnergies[cluster] += pow(AvgEnergies[cluster]-records[i].energy_fields[0], 2)*records[i].occurences;
	}


	//Dividing the variances and dividing the Cluster Sizes to get the Stds
	for(int i = 0; i<nclusters; i++)
	{
		StdEnergies[i]/=(double) (ClFSize[i]);
		StdEnergies[i]=sqrt(StdEnergies[i]);
	}

	//Computing the RepRMSDs
	for(int i = 0; i<nclusters; i++)
	{
		int rep = representatives[i];
		double* coordinates = data[rep]; 
		RepRMSD[i] = compute_distance_pdb_structure((ncolumns/3), orig_atoms, coordinates, pdb_com);
	}


	
	//Finding the maxs and spreads
	for(int i=0; i<nrows; i++)
	{
		int cluster = clusterid[i];
		int rep = representatives[cluster];
		double* rep_structure = data[rep];
		double* structure = data[i];
		double dist = sqrt(3)*compute_RMSD(structure, rep_structure, ncolumns);
		
		spread[cluster] += dist*records[i].occurences;

		if(max[cluster]==-1 || max[cluster]<dist)
			max[cluster]=dist;
		
		
	}



	for(int i=0; i<nclusters; i++)
	{
		spread[i]/=ClFSize[i];
	}

	//Finding the stddevs
	for(int i=0; i<nrows; i++)
	{
		int cluster = clusterid[i];
		int rep = representatives[cluster];
		double* rep_structure = data[rep];
		double* structure = data[i];
		double dist = compute_RMSD(structure, rep_structure, ncolumns);
		stddev[cluster] += pow(dist-spread[cluster], 2)*records[i].occurences;
	}



	for(int i=0; i<nclusters; i++)
	{	
		stddev[i] = stddev[i] / ClFSize[i];
		stddev[i] = sqrt(stddev[i]);
	}
	



	return representatives;
}



double** get_energy_data_set(int nrows, int energy_num, Record* records)
{
	double**energy_dat = malloc(sizeof(double*)*nrows);
	for(int k = 0; k< nrows; k++)
	{
		energy_dat[k]=malloc(sizeof(double)*energy_num);
		for(int n = 0; n< energy_num; n++)
			energy_dat[k][n] = records[k].energy_fields[n];
	}
	return energy_dat;
}

void get_stds_of_clusters(int nclusters, int nrows, int ncols, int clusterid[], double** data, double** cdata, double **std_cdata)
{	
	
	int members[nclusters];
	for(int j=0; j<nclusters; j++)
		members[j]=0;

	for(int j=0;j<nclusters; j++)
		for(int n = 0; n < ncols; n++){
			std_cdata[j][n]=0.;
		}
	for(int k=0; k<nrows; k++)
		for(int n=0; n<ncols; n++)
		{
			int cluster = clusterid[k];
			members[cluster]++;
			std_cdata[cluster][n] += pow(cdata[cluster][n]-data[k][n], 2);
		}

	for(int j=0; j<nclusters; j++)
		for(int n=0; n<ncols; n++)
			std_cdata[j][n] = sqrt(std_cdata[j][n]/members[j]); 
}

void print_table_header()
{
	printf("%7s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s\n",
			"No",
			"ClSize",
			"ClFSize",
			"Repr",
			"ReprE",
			"ClAE",
			"ClAED",
			"RepRMSD",
			"ClFRMSD",
			"Ele",
			"ElDesE",
			"HyDesE",
			"LjE",
			"spread",
			"stddev",
			"max"
			);
}

void print_Cluster(Cluster* cluster)
{
	
	printf("%7i%9li%9li%9i%9.3lf%9.3lf%9.3lf%9.3lf%9.3lf%9.3lf%9.3lf%9.3lf%9.3lf%9.3lf%9.3lf%9.3lf\n",
			cluster->no,
			cluster->ClSize,
			cluster->ClFSize,
			cluster->Repr,
			cluster->ReprE,
			cluster->ClAE,
			cluster->ClAED,
			cluster->RepRMSD,
			cluster->ClFRMSD,
			cluster->Ele,
			cluster->ElDesE,
			cluster->HyDesE,
			cluster->LjE,
			cluster->spread,
			cluster->stddev,
			cluster->max
			);
}


void initializing_Cluster(Cluster* cluster, Record* record, int cluster_num, int Repr, long ClSize, long ClFSize, double AvgEnergy, double StdEnergy, double RepRMSD, double CLFRMSD, double spread, double stddev, double max ,int num_energy_fields, char** Energy_signature)
{
			cluster->no = cluster_num;
			cluster->ClSize = ClSize;
			cluster->ClFSize = ClFSize;
			cluster->Repr = Repr+1;
			cluster->ReprE = 0;
			cluster->ClAE = AvgEnergy;
			cluster->ClAED = StdEnergy;
			cluster->RepRMSD = RepRMSD;
			cluster->ClFRMSD = CLFRMSD;
			cluster->Ele = 0;
			cluster->ElDesE = 0;
			cluster->HyDesE = 0;
			cluster->LjE = 0;
			cluster->spread = spread;
			cluster->stddev = stddev;
			cluster->max = max;

			double* Energies = &(record->energy_fields[0]);
			
			for(int i=0; i<num_energy_fields; i++)
			{	
				if(!memcmp(Energy_signature[i], "TotEn", 5))
					cluster->ReprE = Energies[i];
				else if (!memcmp(Energy_signature[i], "El", 2))
					cluster->Ele += Energies[i];
				else if (!memcmp(Energy_signature[i], "ED", 2))
					cluster->ElDesE+=Energies[i];
				else if (!memcmp(Energy_signature[i], "HD", 2))
					cluster->HyDesE+=Energies[i];
				else if (!memcmp(Energy_signature[i], "rLJ", 3))
					cluster->LjE   +=Energies[i];
				
			}

}







void print_results_description()//(int nclusters, int nows, double** cluster_energies, int* cluster_sizes, int* cluster_occurences, int* cluster_centoid_deviations)
{
	
	printf("Description of the cloumns\n");
	printf("%-10s Cluster Number\n", "No:");
	printf("%-10s Number of entries in the complexes (f55) file used in the cluster\n", "ClSize:");
	printf("%-10s Number of representative entries for the cluster\n", "ClFSize");
	printf("%-10s Representative chosen (corresponding to the line number in the complexes-f55 file)\n", "Repr:");
	printf("%-10s Total interaction energy of the chosen Representative\n", "ReprE:");
	printf("%-10s Average total energy of all cluster members weighted with number of representatvies\n", "ClAE:");
	printf("%-10s Weighted standard deviation of total energy of cluster of cluster entries in the complexes (f55) file\n", "ClAED:");
	printf("%-10s RMSD of the representative to solute 2\n", "RepRMSD:");
	printf("%-10s Average RMSD of the cluster to solute 2\n", "ClFRMSD:");
	printf("%-10s Electrostatic energy\n", "Ele:");
	printf("%-10s Electrstatic desolvation energy of the representative complex\n", "ElDesE:");
	printf("%-10s Hydrophobic desolvation energy of the representative complex\n", "HyDesE:");
	printf("%-10s Lennard-Jones energy of the representative complex\n","LjE:");
	printf("%-10s arithmetic average of rmsds for a given cluster from the representative weighted by occupancy\n", "spread:");
	printf("%-10s Stddev of the rmsds fo a given cluster, weighted by occupancy\n", "stddev:");
	printf("%-10s Maximum rmsd within one cluster from the representative\n", "max:");
}



int main(int argc, char* argv[])
{	
	printf("The number of input arguments is %i\n", argc);
	if (argc < 5 || !strcmp(argv[0], "-h") || !strcmp(argv[0], "--help"))
	{
		printf("Usage:\n");
		printf("\n");
		printf("%s complexes_file pdb_file num_of_cluster method [ncomplexes]\n", argv[0]);
		printf("\n");
		printf("complexes_file: This is the complexes file you would like to analyze.\n");
		printf("It needs to be in the ascii format from SDA (see https://mcm.h-its.org/sda/doc/doc_sda7/complexes_file.html)\n");
		printf("\n");
		printf("pdb_file: This is the pdb file of the ligand. i.e. the molecule that was moving during the BD simulation.\n");
		printf("The RMSD will be measured via the backbone, which are all atoms that are labelled CA in the PDB file\n");
		printf("All lines that have a 'C' in column 14 and an 'A' and column 15 will be considered backbone atoms\n");
		printf("\n");
		printf("num_of_cluster: The cluster representatives to be printed into the summarized complexes file\n");
		printf("\n");
		printf("method: Method to be used for clustering 'a' average, 's' single 'm' maximum and 'c' centroid linkage\n");
		printf("\n");
		printf("ncomplexes(optional): Number of complexes to use for clustering. Defaults to all complexes in the file.\n");
		exit(EXIT_FAILURE);
	}


	int cluster_num;
	sscanf(argv[3], "%i", &cluster_num);
	int num_energy_fields;
	char* filename=argv[1];
	double pdb1_com[3];
	double pdb2_com[3];
	Record* record;
	char**head_lines;
	char** Energy_fields;
	int complexes = read_trajectory_file(filename, &head_lines, &record, &num_energy_fields, &Energy_fields, pdb1_com, pdb2_com);
	int complexes_input;

	printf("The number of energy fields is: %i\n", num_energy_fields);


	if(argc==6)
	{
		sscanf(argv[5], "%i", &complexes_input);
		if(complexes_input > complexes)
		{
			printf("The number of complexes you have chosen is more than there are complexes in the file\n");
			printf("Using all complexes in the file, which is %i\n", complexes);
		}
		else
		{
			complexes = complexes_input;
			printf("Using %i complexes.\n", complexes);
		}

	}

	if(complexes<cluster_num){
		printf("Error: you chose a larger number of representatives than complexes\n");
		exit(EXIT_FAILURE);
	}

	
	printf("The center of mass of the first protein is: %lf\t%lf\t%lf\n", pdb1_com[0], pdb1_com[1], pdb1_com[2]);
	printf("The center of mass of the second protein is: %lf\t%lf\t%lf\n", pdb2_com[0], pdb2_com[1], pdb2_com[2]);
	printf("The number of complexes %i\n", complexes);

	
	pdb_struct PDB;
	char* filename2  = argv[2];
	int line_num = read_pdb_file(&PDB, filename2);
	double** structures=prepare_data(PDB.bb_atom_number, complexes, record, PDB.backbone_coordinates, pdb2_com);

	
	int** mask=prepare_mask(PDB.bb_atom_number, complexes);
	Node* Tree = NULL;


	
	switch(argv[4][0])
	{
		case 'a':
			Tree = hierarchical_average_linkage(complexes, 3*PDB.bb_atom_number, structures, mask);
			break;
		case 's':
			Tree = hierarchical_single_linkage(complexes, 3*PDB.bb_atom_number, structures, mask);
			break;
		case 'm':
			Tree = hierarchical_maximum_linkage(complexes, 3*PDB.bb_atom_number, structures, mask);
			break;
		case 'c':
			Tree = hierarchical_centroid_linkage(complexes, 3*PDB.bb_atom_number, structures, mask);
			break;
		default:
			printf("You did not specify a valid clustering method\n");
			printf("The valid methods are: (a) average, (m) maximum, (s) single, (c) centroid\n");
			printf("\n");
			printf("Try again\n");
			exit(EXIT_FAILURE);
			break;

	}


	if(Tree==NULL)
	{
		printf("Failed to do hierarchical clustering!");
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

	print_table_header();
	

	//Computing the distance matrix for the representatives and initializing the cluster structs
	Cluster* cluster= malloc(sizeof(Cluster)*cluster_num);	
	double** representative_data=malloc(sizeof(double*)*cluster_num);
	for(int i = 0; i<cluster_num; i++)
	{	
		int rep = representatives[i];
		Record* rep_record=&record[rep];
		initializing_Cluster(&cluster[i], rep_record, i+1, rep, ClSize[i], ClFSize[i], AvgEnergies[i], StdEnergies[i],
				RepRMSD[i], CLFRMSD[i], spread[i], stddev[i], max[i], num_energy_fields, Energy_fields);
		representative_data[i] = structures[rep];
	}


	qsort(cluster, cluster_num, sizeof(Cluster), cmp_cluster_sizes);

	for(int i=0; i<cluster_num; i++)
	{	
		cluster[i].no=i+1;
		print_Cluster(&cluster[i]);
	}

	printf("\n");
	printf("\n");
	printf("\n");


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
