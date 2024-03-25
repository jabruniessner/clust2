#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<sys/queue.h>
#include"pdb_structure.h"



int read_pdb_file(pdb_struct* PDB, char Filename[])
{
	PDB->name = Filename;

	FILE* input_file = fopen(Filename, "r");
	
	if(!input_file)
	{
		printf("ERROR: Unable to open file: %s", Filename);
		exit(EXIT_FAILURE);
	}
	
	struct pdb_line_list* line_list = malloc(sizeof(struct pdb_line_list));
	SLIST_INIT(line_list);
	
	unsigned int row = 0;
	size_t buffer_size = 0;
	ssize_t read=0;
	char* buffer=NULL;
	

	int line_num = 0;
	int backbone_atom_number = 0;
	
	while((read = getline(&buffer, &buffer_size, input_file)) !=-1 )
	{	

		//Continue if the line is an END or TER
		if(memcmp(buffer, "ATOM", 4))
		{
			continue;
		}
		
		//skip if the atom is not a backbone atom.
		if ((char*) memchr(&buffer[13], 'C', 2) != NULL ||
		    //!memcmp(&buffer[13], "C ", 2) || 
		    (char*) memchr(&buffer[13], 'O', 2) != NULL ||
		    (char*) memchr(&buffer[13], 'N', 2) != NULL
		   )
		{
			backbone_atom_number++;
		}
		pdb_line* line = malloc(sizeof(pdb_line));

		//Associating pointer of line with the buffer
		line->line = buffer;

		char a[8];
		
		memcpy(a, &buffer[30], 8*sizeof(char));
		sscanf(a, "%lf", &(line->x));

		memcpy(a, &buffer[38], 8*sizeof(char));
                sscanf(a, "%lf", &(line->y));
		
		memcpy(a, &buffer[46], 8*sizeof(char));
                sscanf(a, "%lf", &(line->z));

		SLIST_INSERT_HEAD(line_list, line, next_line);

		//Dissociating pointer from its adress
		buffer=NULL;

		line_num++;
	}
	
	pdb_line* lines =  (pdb_line*) malloc(line_num *sizeof(pdb_line));
	pdb_line* line = SLIST_FIRST(line_list);

	int line_num2 = line_num;

	PDB->backbone_coordinates = (double**) malloc(backbone_atom_number*sizeof(double*));
	PDB->atom_number=line_num;
	PDB->bb_atom_number=backbone_atom_number;
	
	SLIST_FOREACH(line, line_list, next_line)
	{
		lines[--line_num2] =*line;

		if ((char*) memchr(&(line->line[13]), 'C', 2)!=NULL ||
		    //!memcmp(&(line->line[13]), "C ", 2) ||
		    (char*) memchr(&(line->line[13]), 'O', 2)!=NULL || 
		    (char*) memchr(&(line->line[13]), 'N', 2)!=NULL
		    )
		{
			PDB->backbone_coordinates[--backbone_atom_number]=(double*)malloc(3*sizeof(double));
			PDB->backbone_coordinates[backbone_atom_number][0]=line->x;
			PDB->backbone_coordinates[backbone_atom_number][1]=line->y;
			PDB->backbone_coordinates[backbone_atom_number][2]=line->z;
		}
	}


	PDB->lines = lines;
	fclose(input_file);

	return line_num;

}
