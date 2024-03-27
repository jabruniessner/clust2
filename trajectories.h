/* The C clustering library.
 * Copyright (C) 2024 Jakob Niessner.
 * Contact: jabruniessner@gmail.com
 *
 *
 * This program was written at the Heidelberg Institute for theoretical studies,
 * Schloß-Wolfsbrunnenweg 35, 69118 Heidelberg, Germany. 
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


#ifndef TRAJECTORIES_H
#define TRAJECTORIES_H

typedef struct Record{
	long RunNb;
	long StepNb;
	double* Trans;
	double Rot_mat[3][3];
	double* energy_fields;
	double occurences; 
	//double AvEnergy;
	//double StdAvEnergy;
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
