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


#ifndef MATH_FUNC_H
#define MATH_FUNC_H

void cross(double x[3], double y[3], double z[3]);

void transpose(double mat[3][3]);

void mat_vec_mul(double mat[3][3], double orig[3], double res[3]);

void subtract_vectors(double vec1[3], double vec2[3], double result[3]);

void add_vectors(double* vec1, double* vec2, double* result);

double average_array(double* array, long length);

double std_array(double* array, long length);

double compute_RMSD(double* vector1, double* vector2, long length);





#endif
