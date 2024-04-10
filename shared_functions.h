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

#ifndef SHARED_FUNC
#define SHARED_FUNC

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


int cmp_cluster_sizes(const void* clust1, const void* clust2);

int cmp_cluster_sizes_old(const void* clust1, const void* clust2);

double** example_distance_bb(int nrows, int ncols, double** data, int** mask);

Node* hierarchical_single_linkage(int nrows, int ncols, double** data, int** mask);

Node* hierarchical_maximum_linkage(int nrows, int ncols, double** data, int** mask);

Node* hierarchical_average_linkage(int nrows, int ncols, double** data, int** mask);

Node* hierarchical_centroid_linkage(int nrows, int ncols, double** data, int** mask);

int* cutting_hierarchical_tree(Node* tree, int nrows, int level);

double**prepare_data(int number_bb_atom, int number_samples, Record* complexes, double** bb_atoms, double* pdb_com);

double compute_distance_pdb_structure(int number_bb_atoms, double**bb_atoms, double* structure, double* pdb_com);

int**prepare_mask(int number_bb_atom, int number_samples);

int* do_scoring(int nclusters, int nrows, int ncolumns, 
		double** data,  double**orig_atoms, double* pdb_com ,int clusterid[], double** cdata, Record* records, long* ClSize, 
		long* ClFSize, double* AvgEnergies, double* StdEnergies, double* RepRMSD, double* CLFRMSD,
		double* spread, double* stddev, double* max);

double** get_energy_data_set(int nrows, int energy_num, Record* records);

void get_stds_of_clusters(int nclusters, int nrows, int ncols, int clusterid[], double** data, double** cdata, double **std_cdata);

void print_table_header();

void print_Cluster(Cluster* cluster);

void initializing_Cluster(Cluster* cluster, Record* record, int cluster_num, int Repr, long ClSize, long ClFSize, double AvgEnergy, double StdEnergy, double RepRMSD, double CLFRMSD, double spread, double stddev, double max ,int num_energy_fields, char** Energy_signature);

void print_results_description();//(int nclusters, int nows, double** cluster_energies, int* cluster_sizes, int* cluster_occurences, int* cluster_centoid_deviations)

#endif
