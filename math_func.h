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
