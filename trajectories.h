#ifndef TRAJECTORIES_H
#define TRAJECTORIES_H

typedef struct Record{
	long RunNb;
	long StepNb;
	double* Trans;
	double Rot_mat[3][3];
	double* energy_fields;
	double occurences; 
	double AvEnergy;
	double StdAvEnergy;
	char* line;
}Record;

typedef struct Record_list_entry{
	Record record;
	SLIST_ENTRY(Record_list_entry) entry;
}Record_list_entry;

//Defining struct for list head
SLIST_HEAD(Record_list, Record_list_entry);

int read_record_from_line(char line[], Record* record,int num_energy_fields);

int get_energy_signature(char* line, char***Energy_fields);

int read_trajectory_file(char* filename, char***headlines ,Record** records, int* num_energy_fields, char*** Energy_fields , double pdb1_com[], double pdb2_com[]);

void print_record(Record record, int energy_field_num);



#endif
