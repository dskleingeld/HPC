#ifndef SOLVE_H
#define SOLVE_H

#include "lu.h"

bool find_column(CompressedRowMatrix& m, int haystack_row, 
                 int needle_column, int& flat_index);

void matrix_vector_product(CompressedRowMatrix& matrix, double in_vector[], double out_vector[]);
void solve_system(CompressedRowMatrix& lu, PermutationMatrix& p, 
                  double b[], double x[]);


void permute_vector(double pb[], int length, PermutationMatrix& p);
void print_perm(PermutationMatrix& p, const size_t length);
void print_array(const double array[], const size_t length);
void print_array(const double array[], const size_t length, PermutationMatrix& p);


#endif