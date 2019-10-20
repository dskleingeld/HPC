#ifndef SOLVE_H
#define SOLVE_H

#include "lu.h"
#include "libs/dbg/dbg.h"

bool find_column(CompressedRowMatrix& m, int haystack_row, 
                 int needle_column, int& flat_index);

void matrix_vector_product(CompressedRowMatrix& matrix, double in_vector[], double out_vector[]);
void solve_system(CompressedRowMatrix& lu, PermutationMatrix& p, 
                  double b[], double x[]);


#endif