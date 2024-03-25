#ifndef PDB_H
#define PDB_H


//Defining struct for pdbfile_lines
typedef struct pdb_line{
	char * line;
	double x;
	double y;
	double z;

	SLIST_ENTRY(pdb_line) next_line;
}pdb_line;


//Defining struct for list head
SLIST_HEAD(pdb_line_list, pdb_line);


//Defining struct for pdb_file
typedef struct pdb_struct{
	char* name;
	struct pdb_line* lines;
	double** backbone_coordinates;
	long atom_number;
	long bb_atom_number;
}pdb_struct;

int read_pdb_file(pdb_struct* PDB, char Filename[]);


#endif


