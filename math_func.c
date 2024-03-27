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

#include "math_func.h"
#include <math.h>
#include <stdio.h>



void cross(double x[3], double y[3], double z[3])
{
	z[0] = x[1]*y[2]-x[2]*y[1];
	z[1] = x[2]*y[0]-x[0]*y[2];
	z[2] = x[0]*y[1]-y[0]*x[1];
}

void transpose(double mat[3][3])
{
	double mat_tmp[3][3];

	for(int i = 0; i<3; i++)
		for(int j = 0; j<3; j++)
			mat_tmp[i][j]=mat[j][i];

	for(int i = 0; i<3; i++)
                for(int j = 0; j<3; j++)
                        mat[i][j]=mat_tmp[i][j];
}

void mat_vec_mul(double mat[3][3], double orig[3], double res[3])
{
	for(int i = 0; i<3; i++)
	{
		res[i] = mat[i][0]*orig[0]+mat[i][1]*orig[1]+mat[i][2]*orig[2];
	}
}


/*This routine subtracts vec2 from vec1*/
void subtract_vectors(double vec1[3], double vec2[3], double result[3])
{
	for(int i=0; i<3; i++)
	{
		result[i]=vec1[i]-vec2[i];
	}
}

/* This subroutine adds vecotrs vec1 and vec2 and saves the result in result*/
void add_vectors(double* vec1, double* vec2, double* result)
{
	for(int i = 0; i<3; i++)
	{
		result[i]=vec1[i]+vec2[i];
	}
}

double average_array(double* array, long length)
{	
	double avg = 0;
	for(long i = 0; i< length; i++)
		avg+=array[i];
	
	return avg/length;
}

double std_array(double* array, long length)
{
	double std = 0;
	double avg = average_array(array, length);
	for(long i=0; i<length; i++)
		std+=pow(array[i]-avg, 2);

	std = sqrt(std/(length-1));

	return std;

}

double compute_RMSD(double* vector1, double* vector2, long length)
{
	double RMSD=0.;
	for(long i=0; i<length; i++)
	{	
		RMSD+= pow(vector1[i]-vector2[i], 2);
	}
	

	RMSD = sqrt(RMSD/length);


	return RMSD;
}

