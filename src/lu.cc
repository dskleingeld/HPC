#include "lu.h"

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

  auto row_reserved = row_ptr_reserved[row];
  row_ptr_reserved[row] = row_ptr_reserved[replacement_row];
  row_ptr_reserved[replacement_row] = row_reserved;
}

int CompressedRowMatrix::n_elements_in_row(const int row_index){
  return row_ptr_end[row_index] - row_ptr_begin[row_index]+1;
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
  int best_row = pivot_row;
  double largest_val;
  
  int flat_index; //set best_val to pivot row
  if (find_column(lu, pivot_row, pivot_row, flat_index)){
    largest_val = lu.values[flat_index];
  } else {
    //dbg("row has zero as pivot");
    largest_val = 0;
  }

  //find the row with the best pivot value
  for (int replacement_row=pivot_row+1; replacement_row<lu.n_rows; replacement_row++){
    for (int flat_index = lu.row_ptr_begin[replacement_row]; 
        flat_index <= lu.row_ptr_end[replacement_row]; 
        flat_index++){    
      
      int column_index = lu.col_ind[flat_index];
      double value = lu.values[flat_index];
    
      //row is better if the pivot column has a abs(value) larger then the best 
      if (column_index == pivot_row){ //opt. abort if column index > target_row

        if (std::abs(value) > std::abs(largest_val)){
          largest_val = value;
          best_row = replacement_row;
          //if(value==0){dbg("wtf!");}
        }
        break;
      }
    }
  }

  //if(std::isnan(largest_val)){dbg("nan value as pivot");}

  lu.swap_rows(pivot_row, best_row);
  p.mark_swap(pivot_row, best_row);

  //printf("Row: %d, pivot row: %d\n", pivot_row, best_row);


  if (pivot_row!=best_row){
    //printf("row: %d, pivot row: %d, pivot value: %f\n", pivot_row, best_row, largest_val);
    //dump_nonzeros(lu.n_rows, lu.values, lu.col_ind, lu.row_ptr_begin, lu.row_ptr_end);
  }

  return largest_val;
}

bool find_column(CompressedRowMatrix& m, int haystack_row, 
                 int needle_column, int& flat_index){
  
  //go through the haystack
  for(flat_index = m.row_ptr_begin[haystack_row]; 
      flat_index<=m.row_ptr_end[haystack_row]; 
      flat_index++){

    auto hay_column = m.col_ind[flat_index];
    //if (hay_column > needle_column){ //opt Re-enable
    //  return false;
    //} else if (hay_column == needle_column) {
    if (hay_column == needle_column){
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
  //dbg("STOP AND COPY");

  if (old_space_end == MAX_N_ELEMENTS){   //old space is the begin of the array
    free = old_space_end + 1;//next free element one at the begin of the second half of the array
    old_space_end = 2*MAX_N_ELEMENTS;
  } else { //old space at the end of the array
    free = 0;//next free element at begin of array
    old_space_end = MAX_N_ELEMENTS;
  }

  for(int row=0; row<n_rows; row++){
    //move row
    
    int new_idx = free;//move all data
    for(int old_idx=row_ptr_begin[row];
            old_idx<=row_ptr_end[row];  
            old_idx++){
      
      values[new_idx] = values[old_idx];
      col_ind[new_idx] = col_ind[old_idx];
      new_idx++;
    }

    //set pointers right
    row_ptr_begin[row] = free;
    row_ptr_end[row] = new_idx-1;
    row_ptr_reserved[row] = row_ptr_end[row] + N_TO_OVERALLOCATE; 

    free = row_ptr_reserved[row]+1;
  }
}

void CompressedRowMatrix::allocate(int numb_elements, int& ptr_begin, 
                                   int& ptr_reserved){

  if(free>old_space_end){//in new space
    if(free+numb_elements+N_TO_OVERALLOCATE > MAX_N_ELEMENTS*2){
      //stop_and_copy();
    }
  } else {//in old space
    if(free+numb_elements+N_TO_OVERALLOCATE > MAX_N_ELEMENTS){
      //stop_and_copy();
    }
  }

  ptr_begin = free;
  ptr_reserved = free + numb_elements + N_TO_OVERALLOCATE;
  free+=numb_elements + N_TO_OVERALLOCATE+1;
}

//geather
void overwrite_sparse_row_with_dense(DenseIndexedRow& dense_row, 
                                     int row, CompressedRowMatrix& lu){

  //count number of non zeros
  int nnz = 0;
  for(int i=0; i<lu.n_rows; i++){
    if (dense_row.values[i] != 0.0){
      nnz++;
    }
  }

  //+1 as lu.row_ptr_reserved points to the last element still reserved for
  //this row. thus reserved-row is one smaller then the reserved number of elements
  if(lu.row_ptr_reserved[row] - lu.row_ptr_begin[row] +1 < nnz){
    //dbg("allocating");
    //allocate more space
    //opt, find way to just extend reserved space if there is free space
    lu.allocate(nnz, lu.row_ptr_begin[row], 
                lu.row_ptr_reserved[row]);

    //copy values not in dense_row not needed as dense row now contains everything
  }

  //geather into sparse row
  auto flat_index = lu.row_ptr_begin[row];
  for(int column=0; column<lu.n_rows; column++){
    auto value = dense_row.values[column];
    if (value==0){ //skip zero values
      continue;
    }

    lu.values[flat_index] = value;
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
  if (!found){return;} // no pivot in this column => we are done 

  DenseIndexedRow dense_row; //opt make static, and reset each run

  //scatter elements after mult of other_row to a temp row in dense form
  for(int k=lu.row_ptr_begin[other_row]; k<=lu.row_ptr_end[other_row]; k++){
    auto value = lu.values[k]; 
    auto column = lu.col_ind[k];
    dense_row.values[column] = value;
  }
  
  auto mult = lu.values[pivot_column_in_other_row]/pivot;
  //add scaled pivot_row to scatterd other row in dense form
  for(int k=lu.row_ptr_begin[pivot_row]+1; k<=lu.row_ptr_end[pivot_row]; k++){
    auto column = lu.col_ind[k];
    if (column > pivot_column){ //TODO is this check needed?
      dense_row.values[column] -= mult*lu.values[k];
    }
  }
  //set the element in the column below the pivot 
  dense_row.values[pivot_column] = mult;

  //printf("mult: %f, dense row [%i]: %f %f %f %f\n", mult, other_row, dense_row.values[0], dense_row.values[1], dense_row.values[2], dense_row.values[3]);
  overwrite_sparse_row_with_dense(dense_row, other_row, lu);
}

//using the opposites of the multiplies used in the 
//row operations to obtain U, we build L
void lu_factorise(CompressedRowMatrix& lu,
                  PermutationMatrix& p){

  p.identity(lu.n_rows);
  auto n_columns = lu.n_rows;
  
  for(int column=0; column<n_columns; column++){
    double pivot = partial_pivot(p, lu, column);
    if (pivot==0) {continue;} //skip lines with pivot value 0

    //for each row below the current pivot
    for(int row = column+1; row<lu.n_rows; row++){
      add_rows(lu, column, row, pivot);
    }
  }
}
