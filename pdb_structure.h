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


