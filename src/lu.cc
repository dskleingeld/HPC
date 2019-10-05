#include <cstdio>
#include <ctime>
#include <iostream>
#include <cstring>

#include "matrix.h"

/* Code taken from the GLIBC manual.
 *
 * Subtract the ‘struct timespec’ values X and Y,
 * storing the result in RESULT.
 * Return 1 if the difference is negative, otherwise 0.
 */
static int
timespec_subtract (struct timespec *result,
                   struct timespec *x,
                   struct timespec *y)
{
  /* Perform the carry for the later subtraction by updating y. */
  if (x->tv_nsec < y->tv_nsec) {
    int nsec = (y->tv_nsec - x->tv_nsec) / 1000000000 + 1;
    y->tv_nsec -= 1000000000 * nsec;
    y->tv_sec += nsec;
  }
  if (x->tv_nsec - y->tv_nsec > 1000000000) {
    int nsec = (x->tv_nsec - y->tv_nsec) / 1000000000;
    y->tv_nsec += 1000000000 * nsec;
    y->tv_sec -= nsec;
  }

  /* Compute the time remaining to wait.
     tv_nsec is certainly positive. */
  result->tv_sec = x->tv_sec - y->tv_sec;
  result->tv_nsec = x->tv_nsec - y->tv_nsec;

  /* Return 1 if result is negative. */
  return x->tv_sec < y->tv_sec;
}


/* Global variables holding the matrix data. To complete this assignment
 * you are requested to only use arrays and access these arrays with
 * subscripts. Do not use pointers.
 */

const int max_n_elements = 131072;
const int max_n_rows = 16384;
const int max_n_columns = max_n_rows; //matrices must be square

struct CompressedRowMatrix{
  double values[max_n_elements];
  int col_ind[max_n_elements];
  int row_ptr_begin[max_n_rows];
  int row_ptr_end[max_n_rows]; //inclusive (points to last element)

  int n_rows;
};



int numb_elements_in_row(int row, CompressedRowMatrix& matrix){
  return matrix.row_ptr_end[row] - matrix.row_ptr_begin[row];
}

//out_vector needs to be all zeros for length n_rows;
void matrix_vector_product(CompressedRowMatrix& matrix, double in_vector[], double out_vector[]){
  for (int row=0; row<matrix.n_rows; row++){
    for (int element_in_row=0; element_in_row<numb_elements_in_row(row, matrix); element_in_row++){
      auto flat_index = matrix.row_ptr_begin[row] + element_in_row;
      //auto column = col_ind[flat_index];
      out_vector[row] += matrix.values[flat_index] * in_vector[row];
    }
  }
}

void print_array(double vector[], int length){
  for (int i=0; i<length; i++){
    std::cout<<vector[i]<<", ";
  }
  std::cout<<std::endl;
}


void init_array(double array[], int len, double pattern[]){
  for(int i=0; i<len; i++){
    array[i] = pattern[i%2];
  }
}

bool row_has_ok_pivot(CompressedRowMatrix& lu, int current_row, int target_row){
  for (int flat_index = lu.row_ptr_begin[current_row]; 
           flat_index <= lu.row_ptr_end[current_row]; 
           flat_index++){
  
    int column_index = lu.col_ind[flat_index];
    double value = lu.values[flat_index];
    
    if (column_index == target_row && value >= MINIMAL_PIVOT_SIZE){
      return true;
    }
  }
  return false;
}


template<typename T>
void array_to_temp(T array[], T temp[], size_t start, size_t stop){
  for (size_t i = start, j=0; i<stop; i++, j++){
    temp[i] = array[j];
  }
}

template<typename T>
void array_from_temp(T array[], T temp[], size_t start, size_t stop){
  for (size_t i = start, j=0; i<stop; i++, j++){
    array[i] = temp[j];
  }
}

//----row | in between | replacement-row | -----l
//        |1                             |2
//copy from (1-2) into memory
//copy row to the end.
//copy replacement-row back into the array
void swap_rows(CompressedRowMatrix& lu,
               int row, int replacement_row){

  double temp_values[max_n_columns];
  int temp_col_ind[max_n_columns];

  //backup into temp
  size_t n = lu.row_ptr_end[row]-lu.row_ptr_end[replacement_row];
  array_to_temp(lu.values, temp_values, lu.row_ptr_end[row], lu.row_ptr_end[replacement_row]);
  array_to_temp(lu.col_ind, temp_col_ind, lu.row_ptr_end[row], lu.row_ptr_end[replacement_row]);


}

void mark_no_swap(CompressedRowMatrix& p, int row){

}

void mark_swap(CompressedRowMatrix& p, int row, int replacement_row){

}

void init_p(CompressedRowMatrix& p){
  memset(p.row_ptr_begin, 0, max_n_rows);
  memset(p.row_ptr_end, 0, max_n_rows);
}

//TODO, what to do it there is no pivot? (all are zero)
constexpr double MINIMAL_PIVOT_SIZE = 0.1;
void apply_pivots(CompressedRowMatrix& p, CompressedRowMatrix& lu){
  init_p(p);
  for (int row=0; row<lu.n_rows; row++){
    if (row_has_ok_pivot(lu, row, row)) {
      mark_no_swap(p, row);
    } else {

      //there is no pivot for this row at its current position
      //or the pivot is too small.
      //iter through column, element is replacement_pivot 
      for (int replacement_row=row; replacement_row<lu.n_rows; replacement_row++){
        //TODO what if no row has a suitable replacement? can we leave it?
        if (row_has_ok_pivot(lu, replacement_row, row)){
          swap_rows(lu, row, replacement_row);
          mark_swap(p, row, replacement_row);
          break;
        }
      }
    }
}

//TODO what happens when there is no pivot?
void factorise_in_place(CompressedRowMatrix& a){

}

//using the opposites of the multiplies used in the 
//row operations to obtain U, we build L
void lu_factorise(CompressedRowMatrix& a, 
                  CompressedRowMatrix& lu,
                  CompressedRowMatrix& p){

  lu = a;                                               
  apply_pivots(p, lu);
  factorise_in_place(lu);
}


int
main(int argc, char **argv)
{
  if (argc != 2){
    fprintf(stderr, "usage: %s <filename>\n", argv[0]);
    return -1;
  }

  int nnz, n_rows, n_cols;
  bool ok(false);

  CompressedRowMatrix a;
  ok = load_matrix_market(argv[1], max_n_elements, max_n_rows,
                          nnz, a.n_rows, n_cols,
                          a.values, a.col_ind, a.row_ptr_begin, a.row_ptr_end);
  if (!ok){
    fprintf(stderr, "failed to load matrix.\n");
    return -1;
  }

  //solution vectors
  double solution_vector[max_n_rows];
  double pattern[] = {1., 1.};
  //double pattern[] = {.1, .1};
  //double pattern[] = {1., -1.};
  //double pattern[] = {5.,-5.};
  //double pattern[] = {100.,-100.};
  init_array(solution_vector, n_rows, pattern);

  double b[max_n_rows];
  matrix_vector_product(a, solution_vector, b);
  print_array(solution_vector, n_rows);
  //print_array(b, n_rows);

  /* For debugging, can be removed when implementation is finished. */
  //dump_nonzeros(n_rows, values, col_ind, row_ptr_begin, row_ptr_end);



  struct timespec start_time;
  clock_gettime(CLOCK_REALTIME, &start_time);

  /* Perform LU factorization here */

  struct timespec end_time;
  clock_gettime(CLOCK_REALTIME, &end_time);


  struct timespec elapsed_time;
  timespec_subtract(&elapsed_time, &end_time, &start_time);

  double elapsed = (double)elapsed_time.tv_sec +
      (double)elapsed_time.tv_nsec / 1000000000.0;
  fprintf(stderr, "elapsed time: %f s\n", elapsed);

  return 0;
}
