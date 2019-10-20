#include "solve.h"

//out_vector needs to be all zeros for length n_rows;
void matrix_vector_product(CompressedRowMatrix& matrix, double in_vector[], double out_vector[]){
  for (int row=0; row<matrix.n_rows; row++){
    for (int element_in_row=0; element_in_row<matrix.n_elements_in_row(row); element_in_row++){
      auto flat_index = matrix.row_ptr_begin[row] + element_in_row;
      //auto column = col_ind[flat_index];
      out_vector[row] += matrix.values[flat_index] * in_vector[row];
    }
  }
}

//give na lu factord matrix with permutation matrix and a vector b 
//solves the system Ax=b
//
// PAx = Pb
// LUx = Pb
//
// Ly = Pb -> solve for y
// Ux = y -> solve for x
void solve_system(CompressedRowMatrix& lu, PermutationMatrix& p, 
                  double org_b[], double x[]){

    double pb[MAX_N_ROWS];
    double y[MAX_N_ROWS];
    //permute vector b to match up LU vect.
    for(int row=0; row<lu.n_rows; row++){
        pb[row] = org_b[p.permuted_to_original_index[row]];
    }

    for (int i=0; i<4; i++) {dbg(pb[i]);}

    //forward subsitution, solves y from Ly=Pb
    for(int row=0; row<lu.n_rows; row++){
        y[row] = pb[row];

        //iterate over all columns of A
        //if A[] is zero nothing happens to c =>
        //we only need to iterate the nonzero elements
        //for(int column=0; column<row-1; column++){
        for(int flat_idx=lu.row_ptr_begin[row]; 
            flat_idx<lu.row_ptr_end[row]; flat_idx++){

            auto A = lu.values[flat_idx];
            y[row] -= A * pb[row];
        }
    }

    for (int i=0; i<4; i++) {dbg(y[i]);}   

    //back subsitution, solves x from Ux=y
    for(int row = lu.n_rows-1; row>-1; row--){
        x[row] = y[row];
        //iterate over all columns of A
        //if A[] is zero nothing happens to c =>
        //we only need to iterate the nonzero elements
        //for(int column=row+1; column<row; column++){
        int flat_idx;
        if (find_column(lu, row, row+1, flat_idx)){
            for(;flat_idx<lu.row_ptr_end[row]; flat_idx++){
                auto A = lu.values[flat_idx];
                auto column = lu.col_ind[flat_idx];
                x[row] -= A * x[column];
            }
        }
        double A;
        if (find_column(lu,row,row,flat_idx)){
            A = lu.values[flat_idx];
        } else {
            A = 1.; 
        }
        x[row] = x[row] / A;
    }
}