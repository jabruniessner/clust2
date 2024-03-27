/* The clust2 program.
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
#include "trajectories.h"
#include "math_func.h"

int get_energy_signature(char* line, char***Energy_fields)
{
	int line_length = strlen(line);
	char* linecp = malloc(sizeof(char)*(line_length+1));
	strcpy(linecp, line);
	char* token = strtok(linecp, " ");
	int j=0;
	
	Energy_fields[0] = malloc(sizeof(char*)*100);

	while(memcmp(token, "TotEn", 5)){
		token = strtok(NULL, " ");
	}

	
	while( token != NULL)
	{

		if(!memcmp(token, "Occur", 5))
		{
			break;
		}

		Energy_fields[0][j]=malloc(sizeof(char)*(strlen(token)+1));
		strcpy(Energy_fields[0][j], token);
		token = strtok(NULL, " ");
		j++;




	}


	free(linecp);

	return j;
}



int read_record_from_line(char line[], Record* record, int num_energy_fields)
{	
	record->line = malloc(sizeof(char)*strlen(line));
	strcpy(record->line, line);
	
	const char* delimiter=" ";
	
	char* token = strtok(line, delimiter);
	sscanf(token, "%li", &(record->RunNb));
	

	token = strtok(NULL, delimiter);
	sscanf(token, "%li", &(record->StepNb));

	record->Trans = (double*)malloc(sizeof(double)*3);
	for(int i = 0; i<3; i++)
	{
		token=strtok(NULL, delimiter);
		sscanf(token, "%lf", &(record->Trans[i]));
	}

	for(int i = 0; i<3; i++)
        {
                token=strtok(NULL, delimiter);
                sscanf(token, "%lf", &(record->Rot_mat[0][i]));
        }

	for(int i = 0; i<3; i++)
        {
                token=strtok(NULL, delimiter);
                sscanf(token, "%lf", &(record->Rot_mat[1][i]));
        }

	//Creating the third row of the matrix
	cross(record->Rot_mat[0], record->Rot_mat[1], record->Rot_mat[2]);
	//transposing the rotation matrix in order 
	transpose(record->Rot_mat);

	record->energy_fields = (double*) malloc(num_energy_fields*sizeof(double));

	for(int i=0; i<num_energy_fields; i++)
	{
		token=strtok(NULL, delimiter);
		sscanf(token, "%lf", &(record->energy_fields[i]));
	}

	token=strtok(NULL, delimiter);
	sscanf(token, "%lf", &(record->occurences));

	//token=strtok(NULL, delimiter);
	//sscanf(token, "%lf", &(record->AvEnergy));

	//token=strtok(NULL, delimiter);
	//sscanf(token, "%lf", &(record->StdAvEnergy));
	
}



int read_trajectory_file(char* filename, char*** head_lines, Record** records, int* num_energy_fields, char*** Energy_fields ,double pdb1_com[], double pdb2_com[])
{	
	FILE* traj_file = fopen(filename, "r");

	const char* delimiter = " ";
	
	size_t buffer_size = 0;
        char * buffer = NULL;
        ssize_t read;

	head_lines[0] = malloc(sizeof(char*)*4);

	if((read = getline(&buffer, &buffer_size, traj_file)) == -1 )
	{
		printf("There is something wrong with your trajectories file.\nThis is not an SDA complexes file.\n");
	}

	printf("Made it until right before energy_signature function call\n");
	*num_energy_fields=get_energy_signature(buffer, Energy_fields);
	printf("The number of energy fields is %i\n", *num_energy_fields);

	

	//saving the first head_line
	head_lines[0][0] = (char*)malloc(sizeof(char)*strlen(buffer));
	memcpy(head_lines[0][0], buffer, strlen(buffer));
	
	if((read = getline(&buffer, &buffer_size, traj_file)) == -1 )
	{
		printf("There is something wrong with your trajectories file.\nThis is not an SDA complexes file.\n");
		exit(EXIT_FAILURE);
	}

	//saving the second head_line
        head_lines[0][1] = (char*)malloc(sizeof(char)*strlen(buffer));
        memcpy(head_lines[0][1], buffer, strlen(buffer));




	char buffer2[strlen(buffer)-1];



	if((read = getline(&buffer, &buffer_size, traj_file)) == -1 || buffer[0]!='#')
	{
		printf("The first center of mass is missing!");
		exit(EXIT_FAILURE);
	}

	//saving the third head_line
        head_lines[0][2] = (char*)malloc(sizeof(char)*strlen(buffer));
        memcpy(head_lines[0][2], buffer, strlen(buffer));


	memcpy(buffer2, &buffer[1], sizeof(char)*strlen(buffer)-1);

	
	char* token = strtok(buffer2, delimiter);
	int tok_number=0;
	
	while(token != NULL)
	{	

		sscanf(token, "%lf", &(pdb1_com[tok_number]));
		token = strtok(NULL, delimiter);
		tok_number++;

	}
	
	if((read = getline(&buffer, &buffer_size, traj_file)) == -1 || buffer[0]!='#')
        {
                printf("The second center of mass is missing!");
        }

	//saving the fourth head_line
        head_lines[0][3] = (char*)malloc(sizeof(char)*strlen(buffer));
        memcpy(head_lines[0][3], buffer, strlen(buffer));

	memcpy(buffer2, &buffer[1], sizeof(char)*24);
	token = strtok(buffer2, delimiter);
	tok_number=0;
	while(token != NULL)
        {
                sscanf(token, "%lf", &pdb2_com[tok_number++]);
                token = strtok(NULL, delimiter);

        }

	struct Record_list rec_list;
	SLIST_INIT(&rec_list);


	//Now reading the trajectories file
	int record_num=0;
	while((read=getline(&buffer, &buffer_size, traj_file))!=-1)
	{
		
		Record_list_entry* list_entr = malloc(sizeof(Record_list_entry));
		read_record_from_line(buffer, &(list_entr->record), *num_energy_fields);
		SLIST_INSERT_HEAD(&rec_list, list_entr, entry);
		record_num++;

	}

	Record* records2 = *records = (Record*)malloc(sizeof(Record)*record_num);
	
	Record_list_entry* list_entr;
	int j = 1;
	SLIST_FOREACH(list_entr, &rec_list, entry)
	{
		records2[record_num-j]=list_entr->record;
		j++;
	}


	while(!SLIST_EMPTY(&rec_list)) {

                Record_list_entry *list_entry = SLIST_FIRST(&rec_list);
                SLIST_REMOVE_HEAD(&rec_list, entry);
                free(list_entry);
        }

#ifdef DEBUG
	for(int i= 0; i<3; i++)
		print_record(records2[i], num_energy_fields);
#endif

	return record_num;

	
};


void print_record(Record record, int energy_field_num)
{
        printf("The values are: \n");
        printf("RunNb: %li\n", record.RunNb);
        printf("StepNb: %li\n", record.StepNb);
        printf("Trans X1: %lf\n", record.Trans[0]);
        printf("Trans X2: %lf\n", record.Trans[1]);
        printf("Trans X3: %lf\n", record.Trans[2]);
        printf("Rot1 x: %lf\t%lf\t%lf\n", record.Rot_mat[0][0], record.Rot_mat[1][0], record.Rot_mat[2][0]);
        printf("Rot1 y: %lf\t%lf\t%lf\n", record.Rot_mat[0][1], record.Rot_mat[1][1], record.Rot_mat[2][1]);
        printf("Rot1 z: %lf\t%lf\t%lf\n", record.Rot_mat[0][2], record.Rot_mat[1][2], record.Rot_mat[2][2]);

	printf("Energyvalues: ");
	for(int i=0; i<energy_field_num; i++)
		printf("\t%lf", record.energy_fields[i]);

	printf("\n");
        printf("Occurences: %lf\n", record.occurences);
        //printf("AvEnergy: %lf\n", record.AvEnergy);
        //printf("StdAvEnergy: %lf\n", record.StdAvEnergy);
	printf("Line:\n");
	printf("%s", record.line);
};










