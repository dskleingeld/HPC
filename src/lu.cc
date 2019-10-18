#include "lu.h"
#include "libs/dbg/dbg.h"

//regex to remove debug statements: "[ \t]+dbg\([0-z.\*]+\);\n"

void PermutationMatrix::mark_swap(const int row, const int replacement_row){
  permuted_to_original_index[row] = replacement_row;
  permuted_to_original_index[replacement_row] = row;
}

void PermutationMatrix::identity(size_t n_rows){
  for(size_t row_idx=0; row_idx<n_rows; row_idx++){
    permuted_to_original_index[row_idx] = row_idx;
  }
}

//----row | in between | replacement-row | -----l
//        |1                             |2
//copy from (1-2) into memory
//copy row to the end.
//copy replacement-row back into the array
void CompressedRowMatrix::swap_rows(const int row, const int replacement_row){

  auto row_begin = row_ptr_begin[row];
  row_ptr_begin[row] = row_ptr_begin[replacement_row];
  row_ptr_begin[replacement_row] = row_begin;

  auto row_end = row_ptr_end[row];
  row_ptr_end[row] = row_ptr_end[replacement_row];
  row_ptr_end[replacement_row] = row_end;  
}

int CompressedRowMatrix::n_elements_in_row(const int row_index){
  return row_ptr_end[row_index] - row_ptr_begin[row_index];
}

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

void print_array(const double array[], const size_t length){
  std::cout<<"array: [";
  for (size_t i=0; i<length; i++){
    std::cout<<array[i]<<", ";
  }
  std::cout<<"]"<<std::endl;
}


void init_array(double array[], const int len, const double pattern[]){
  for(int i=0; i<len; i++){
    array[i] = pattern[i%2];
  }
}

bool row_has_ok_pivot(CompressedRowMatrix& lu, const int current_row, 
                      const int diagonal, int& flat_index){
  //dbg("row_has_ok_pivot");
  //dbg(lu.row_ptr_begin[current_row]);
  //dbg(lu.row_ptr_end[current_row]);
  for (flat_index = lu.row_ptr_begin[current_row]; 
       flat_index <= lu.row_ptr_end[current_row]; 
       flat_index++){
  
    int column_index = lu.col_ind[flat_index];
    double value = lu.values[flat_index];
    
    //if there is no value that is a pivot, the value must be zero which is to small
    if (column_index == diagonal){ //opt. abort if column index > target_row
      if (value >= MINIMAL_PIVOT_SIZE){
        return true;
      }
    }
  }
  return false;
}

//TODO, what to do it there is no pivot? (all are zero)
double partial_pivot(PermutationMatrix& p, CompressedRowMatrix& lu, int pivot_row){
  int pivot_flat_index;
  for (int row=pivot_row; row<lu.n_rows; row++){
    if (row_has_ok_pivot(lu, row, row, pivot_flat_index)) {
      return lu.values[pivot_flat_index];
    } else {

      //there is no pivot for this row at its current position
      //or the pivot is too small.
      //iter through column, element is replacement_pivot 
      for (int replacement_row=row; replacement_row<lu.n_rows; replacement_row++){
        //TODO what if no row has a suitable replacement? can we leave it?
        if (row_has_ok_pivot(lu, replacement_row, row, pivot_flat_index)){
          lu.swap_rows(row, replacement_row);
          p.mark_swap(row, replacement_row);
          return lu.values[pivot_flat_index];
        }
      }
    }
  }
  //no ok row could be found
  for(int flat_index=lu.row_ptr_begin[pivot_row]; flat_index<=lu.row_ptr_end[pivot_row];
      flat_index++){

    if (lu.col_ind[flat_index] == pivot_row){
      return lu.values[flat_index];
    }
  }
  return 0.0;
}

bool find_column(CompressedRowMatrix& m, int haystack_row, 
                 int needle_column, int& flat_index){
  
  //go through the haystack
  for(flat_index = m.row_ptr_begin[haystack_row]; 
      flat_index<=m.row_ptr_end[haystack_row]; 
      flat_index++){

    auto hay_column = m.col_ind[flat_index];
    if (hay_column > needle_column){
      return false;
    } else if (hay_column == needle_column) {
      return true;
    }
  }
  return false;
}

template<typename T>
void move_array(T array[], int source, int destination, int len){
  for (int i=0; i<len; i++){
    array[destination+i] = array[source+i];
  }
}

//opt last run of stop and copy should not buffer extra space
void CompressedRowMatrix::stop_and_copy(){
  //old space at the begin of the array
  if (old_space_end == MAX_N_ELEMENTS){
    free = old_space_end + 1;
    old_space_end = 2*MAX_N_ELEMENTS;
  } else { //old space at the end of the array
    free = 0;
    old_space_end = MAX_N_ELEMENTS;
  }

  for(int row=0; row<n_rows; row++){
    //move row
    auto len = row_ptr_end[row] - row_ptr_begin[row];
    move_array(values, row_ptr_begin[row], free, len);
    
    //set pointers right
    row_ptr_begin[row] = free;
    row_ptr_end[row] = free + len;
    row_ptr_reserved[row] = free + len + N_TO_OVERALLOCATE; 

    free += len + N_TO_OVERALLOCATE;
  }
}

void CompressedRowMatrix::allocate(int numb_elements, int& ptr_begin, 
                                   int& ptr_end, int& ptr_reserved){

  if(old_space_end-free<(numb_elements+N_TO_OVERALLOCATE)){
    //do stop and copy to compact
    stop_and_copy();
  }

  ptr_begin = free;
  ptr_end = free+numb_elements;
  ptr_reserved = free + numb_elements +N_TO_OVERALLOCATE;
  free+=numb_elements +N_TO_OVERALLOCATE;
}

int DenseIndexedRow::make_sorted_col_ind(int sorted_col_ind[]){
  int k=0, i=0, j=init_columns;
  while(i<init_columns && j<added_columns){
    if(used_col_ind[i] < used_col_ind[j]){
      sorted_col_ind[k] = used_col_ind[i];
      i++; k++;
    } else if (used_col_ind[i] == used_col_ind[j]) {
      sorted_col_ind[k] = used_col_ind[i];
      i++; j++; k++;
    } else {
      sorted_col_ind[k] = used_col_ind[j];
      j++; k++;
    }
  }
  //if not done with i, copy from the i area 
  //in used_col_ind to the sorted array
  for(;i<init_columns; i++, k++){
    sorted_col_ind[k] = used_col_ind[i];
  } 
  //if not done with j, copy from the j area 
  //in used_col_ind to the sorted array
  for(;j<added_columns; j++, k++){
    sorted_col_ind[k] = used_col_ind[j];
  }
  return k;
}

//geather
void overwrite_sparse_row_with_dense(DenseIndexedRow& dense_row, 
                                     int row, CompressedRowMatrix& lu,
                                     int dense_row_offset){

  static int sorted_col_ind[MAX_N_COLLUMNS];
  auto numb_elements = dense_row.make_sorted_col_ind(sorted_col_ind);

  //+1 as lu.row_ptr_reserved points to the last element still reserved for
  //this row. thus reserved-row is one smaller then the reserved number of elements
  if(lu.row_ptr_reserved[row] - lu.row_ptr_begin[row] +1 < numb_elements){
    dbg("ALLOCATING!");
    //allocate more space
    //opt, find way to just extend reserved space if there is free space
    lu.allocate(numb_elements, lu.row_ptr_begin[row], 
                lu.row_ptr_end[row], lu.row_ptr_reserved[row]);

    //TODO: copy values not in dense_row!
  }

  //geather into sparse row
  auto flat_index = lu.row_ptr_begin[row]+dense_row_offset;
  dbg(dense_row_offset);
  dbg(numb_elements);
  for(int i=0; i<numb_elements; i++){
    auto column = sorted_col_ind[i];
    auto value = dense_row.values[column];
    if (value==0){
      dbg("skipping value");
      continue;
    }
    
    lu.values[flat_index] = value;
    dbg(value);
    dbg(column);
    lu.col_ind[flat_index] = column;
    flat_index++;
  }
  //update row ptr in case row shrunk
  lu.row_ptr_end[row] = flat_index-1;
}

//adds pivot_row, scaled by pivot, to other_row
void add_rows(CompressedRowMatrix& lu, int pivot_row, int other_row, double pivot){

  //walk until we get the flat index of the pivot column in the other row
  auto pivot_column = pivot_row;
  int pivot_column_in_other_row; //set by find column
  bool found = find_column(lu, other_row, pivot_column, pivot_column_in_other_row);
  if (!found){dbg("not found"); return;} //we are done 
  
  //set the element in the column below the pivot 
  auto mult = lu.values[pivot_column_in_other_row]/pivot;

  DenseIndexedRow dense_row; //opt make static, and reset each run
  dense_row.values[pivot_column] = mult;
  dense_row.used_col_ind[dense_row.init_columns] = pivot_column;
  dense_row.init_columns++;

  //scatter the other row to a temp row in dense form
  for(int k=pivot_column_in_other_row+1; k<=lu.row_ptr_end[other_row]; k++){
    auto value = lu.values[k]; 
    auto column = lu.col_ind[k];
    dense_row.values[column] = value;
    dense_row.used_col_ind[dense_row.init_columns] = column;
    dense_row.init_columns++;
  }
  
  //add scaled pivot_row to scatterd other row in dense form
  for(int k=lu.row_ptr_begin[pivot_row]+1; k<=lu.row_ptr_end[pivot_row]; k++){
    auto column = lu.col_ind[k];
      if (column > pivot_column){
      dense_row.values[column] -= mult*lu.values[k];
      dense_row.added_columns++;
    }
  }
  std::cout<<"dense row: "
           <<dense_row.values[0]<<"  "
           <<dense_row.values[1]<<"  "
           <<dense_row.values[2]<<"  "
           <<dense_row.values[3]<<"  "<<std::endl;
  int dense_row_offset = pivot_column_in_other_row - lu.row_ptr_begin[other_row];
  overwrite_sparse_row_with_dense(dense_row, other_row, lu, dense_row_offset);
}

//using the opposites of the multiplies used in the 
//row operations to obtain U, we build L
void lu_factorise(CompressedRowMatrix& lu,
                  PermutationMatrix& p){

  p.identity(lu.n_rows);
  auto n_columns = lu.n_rows;
  
  for(int column=0; column<n_columns; column++){
      double pivot = partial_pivot(p, lu, column);
      dump_nonzeros(lu.n_rows, lu.values, lu.col_ind, lu.row_ptr_begin, lu.row_ptr_end);

    //for each row below the current pivot
    for(int row = column+1; row<lu.n_rows; row++){
          add_rows(lu, column, row, pivot);
    }
  }
}
