#include "solve.h"

//out_vector needs to be all zeros for length n_rows;
void matrix_vector_product(CompressedRowMatrix& matrix, double in_vector[], double out_vector[]){
  for (int row=0; row<matrix.n_rows; row++){
    for (int element_in_row=0; element_in_row< (matrix.row_ptr_end[row] - matrix.row_ptr_begin[row]+1); element_in_row++){
      auto flat_index = matrix.row_ptr_begin[row] + element_in_row;
      //auto column = col_ind[flat_index];
      out_vector[row] += matrix.values[flat_index] * in_vector[row];
    }
  }
}

void get_permuted_vector(const double array[], double pb[], 
                         int length, PermutationMatrix& p){
  //copy the array into pb
  for (int i=0; i<length; i++){pb[i] = array[i];}

  //create index from permutationMatrix
  int index[MAX_N_ROWS];
  for(int i=0; i<length; i++){index[i]=i;}

  for(int i=0; i<length; i++){
    auto v = p.permuted_to_original_index[i];
    auto temp = index[i];
    index[i] = index[v];
    index[v] = temp;
  }

  //permute vector b to match up LU vect.
  for(int i=0; i<length; i++){
    pb[index[i]] = array[i];
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
    double c[MAX_N_ROWS];
    //permute vector b to match up LU vect.
    get_permuted_vector(org_b, pb, lu.n_rows, p);

    //forward subsitution, solves y from Ly=Pb
    for(int row=0; row<lu.n_rows; row++){
        c[row] = pb[row];

        //iterate over all columns of A
        //if A[] is zero nothing happens to c =>
        //we only need to iterate the nonzero elements
        //for(int column=0; column<row-1; column++){
        for(int flat_idx=lu.row_ptr_begin[row]; 
            flat_idx<=lu.row_ptr_end[row]; flat_idx++){

            //only iterate till column<row-1
            auto column = lu.col_ind[flat_idx];
            if(column>row-1){break;}

            auto A = lu.values[flat_idx];
            c[row] -= A * c[column];
        }
    }

    //back subsitution, solves x from Ux=y
    for(int row = lu.n_rows-1; row>=0; row--){
        x[row] = c[row];
        //iterate over all columns of A
        //if A[] is zero nothing happens to c =>
        //we only need to iterate the nonzero elements
        //for(int column=row+1; column<row; column++){
        
        //find first non zero element in row of A starting at column row+1
        for (int flat_idx=lu.row_ptr_begin[row]; flat_idx<lu.row_ptr_end[row]; flat_idx++){
          //only truely start iterating over columns once we are past column row+1
          if (lu.col_ind[flat_idx]<row+1){ continue; }

          auto A = lu.values[flat_idx];
          auto column = lu.col_ind[flat_idx];
          x[row] -= A * x[column];
        }
        double A;
        int flat_idx;
        if (find_column(lu,row,row,flat_idx)){
          A = lu.values[flat_idx];
        } else {
          dbg(row);
          A = 1.; 
        }
        x[row] = x[row] / A;
    }
}

void print_perm(PermutationMatrix& p, const size_t length){
  std::cout<<"[ ";
  for(size_t i = 0; i<length; i++){
    printf(" %d",p.permuted_to_original_index[i]);
  }
  std::cout<<"]"<<std::endl;
}

void print_array(const double a[], const size_t length){
  std::cout<<"[ ";
  size_t i = 0;
  for (; i<length-4; i+=4){
    printf(" %.8e %.8e %.8e %.8e\n", a[i], a[i+1], a[i+2], a[i+3]);
  }
  std::cout<<" ";
  for (; i<length; i++){
    printf("%.8e", a[i]);
  }
  std::cout<<"]"<<std::endl;
}

void print_array(const double array[], const size_t length, PermutationMatrix& p){

    double pb[MAX_N_ROWS];
    get_permuted_vector(array, pb, length, p);

    std::cout<<"[ ";
    size_t i = 0;
    for (; i<length-4; i+=4){
      printf(" %.8e %.8e %.8e %.8e\n", pb[i], pb[i+1], pb[i+2], pb[i+3]);
    }
    std::cout<<" ";
    for (; i<length; i++){
      printf("%.8e", pb[i]);
    }
    std::cout<<"]"<<std::endl;
}


