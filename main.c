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


int main(int argc, char* argv[])
{	
	printf("The number of input arguments is %i\n", argc);
	if (argc < 5 || !strcmp(argv[1], "-h") || !strcmp(argv[1], "--help"))
	{
		printf("Usage:\n");
		printf("\n");
		printf("%s complexes_file pdb_file num_of_cluster method [ncomplexes]\n", argv[0]);
		printf("\n");
		printf("complexes_file: This is the complexes file you would like to analyze.\n");
		printf("It needs to be in the ascii format from SDA (see https://mcm.h-its.org/sda/doc/doc_sda7/complexes_file.html)\n");
		printf("\n");
		printf("pdb_file: This is the pdb file of the ligand. i.e. the molecule that was moving during the BD simulation.\n");
		printf("The RMSD will be measured via the backbone, which are all atoms that are labelled C, N or O in the PDB file\n");
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
	

	return 0;
}
