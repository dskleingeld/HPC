#ifndef LU_H
#define LU_H

#include <cstdio>
#include <ctime>
#include <iostream>
#include <cstring>
#include <cmath>

#include "libs/dbg/dbg.h"
#include "matrix.h"

/* Global variables holding the matrix data. To complete this assignment
 * you are requested to only use arrays and access these arrays with
 * subscripts. Do not use pointers.
 */
constexpr double MINIMAL_PIVOT_SIZE = 0.1;
constexpr int N_TO_OVERALLOCATE = 0;//10;
//constexpr int MAX_N_ELEMENTS = 131072; //Origional
constexpr int MAX_N_ELEMENTS = 4000000; //4m increased to prevent stop and copy
//constexpr int MAX_N_ELEMENTS = 400000000; //400m increased to prevent stop and copy
                            // 116079076; nopoly size

constexpr int MAX_N_ROWS = 16384;
constexpr int MAX_N_COLLUMNS = MAX_N_ROWS; //matrices must be square

struct CompressedRowMatrix{
  double values[MAX_N_ELEMENTS*2];
  int col_ind[MAX_N_ELEMENTS*2];
  int free; //points to the first free element in the memory array
  int old_space_end = MAX_N_ELEMENTS; //end of fromspace
  
  int row_ptr_begin[MAX_N_COLLUMNS]; //opt. move to array of structs to increase data locality
  int row_ptr_end[MAX_N_COLLUMNS]; //inclusive (points to last element)
  int row_ptr_reserved[MAX_N_COLLUMNS]; //inclusive (points to last element)

  int n_rows;
  void swap_rows(const int row, const int replacement_row);
  int n_elements_in_row(const int row_index);
  void allocate(int numb_elements, int& ptr_begin, int& ptr_reserved);
  void stop_and_copy();

  void init_memory_management(){
    free = row_ptr_end[n_rows-1]+1;
    for(int i=0; i<n_rows; i++){
      row_ptr_reserved[i] = row_ptr_end[i];
    }
  }
};

struct PermutationMatrix {
  int permuted_to_original_index[MAX_N_ROWS];
  void identity(size_t n_rows);
  void mark_swap(const int row, const int replacement_row);
};

struct DenseIndexedRow {
  double values[MAX_N_COLLUMNS] = {0};
};

void init_array(double array[], const int len, const double pattern[]);
void matrix_vector_product(CompressedRowMatrix& matrix, double in_vector[], double out_vector[]);
bool find_column(CompressedRowMatrix& m, int haystack_row, int needle_column, int& flat_index);
void lu_factorise(CompressedRowMatrix& lu, PermutationMatrix& p);

#endif