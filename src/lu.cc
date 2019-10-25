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

  auto row_reserved = row_ptr_reserved[row];
  row_ptr_reserved[row] = row_ptr_reserved[replacement_row];
  row_ptr_reserved[replacement_row] = row_reserved;
}

int CompressedRowMatrix::n_elements_in_row(const int row_index){
  return row_ptr_end[row_index] - row_ptr_begin[row_index]+1;
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
  int best_row = pivot_row;
  double largest_val;
  
  int flat_index; //set best_val to pivot row
  if (find_column(lu, pivot_row, pivot_row, flat_index)){
    largest_val = lu.values[flat_index];
  } else {
    //dbg("row has zero as pivot");
    largest_val = 0;
  }

  //dbg(pivot_row);
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
          if(value==0){dbg("wtf!");}
        }
        break;
      }
    }
  }

  //if (pivot_row==101){ while(1){;}}
  if(std::isnan(largest_val)){dbg("nan value as pivot");}
  /*std::cout <<"org row: "<<pivot_row
            <<" replaced with: "<<best_row
            <<" new pivot value: "<<largest_val
            << std::endl;*/
  //printf("row: %d, pivot row: %d, pivot value: %f\n", pivot_row, best_row, largest_val);

  lu.swap_rows(pivot_row, best_row);
  p.mark_swap(pivot_row, best_row);

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
  dbg("STOP AND COPY");
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
    
    int new_idx = free;
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
    row_ptr_reserved[row] = new_idx-1 + N_TO_OVERALLOCATE; 

    free = row_ptr_reserved[row]+1;
  }
}

void CompressedRowMatrix::allocate(int numb_elements, int& ptr_begin, 
                                   int& ptr_reserved){

  if(old_space_end-free<(numb_elements+N_TO_OVERALLOCATE)){
    //do stop and copy to compact
    //stop_and_copy(); //TODO re-enable for larger matrices 
  }

  ptr_begin = free;
  ptr_reserved = free + numb_elements +N_TO_OVERALLOCATE;
  free+=numb_elements +N_TO_OVERALLOCATE+1;
}

int DenseIndexedRow::make_sorted_col_ind(int sorted_col_ind[]){
  int k=0, i=0, j=init_columns;

  while(i<init_columns && j<init_columns+added_columns){
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
  for(;j<init_columns+added_columns; j++, k++){
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

    //allocate more space
    //opt, find way to just extend reserved space if there is free space
    lu.allocate(numb_elements, lu.row_ptr_begin[row], 
                lu.row_ptr_reserved[row]);

    //copy values not in dense_row
    auto old_begin = lu.row_ptr_begin[row];
    for(int old_idx=old_begin, new_idx=lu.row_ptr_begin[row]; //opt copy less (need only everything in front of pivot)
        old_idx<=lu.row_ptr_end[row]; old_idx++, new_idx++){
      
      lu.values[new_idx] = lu.values[old_idx];
      lu.col_ind[new_idx] = lu.col_ind[old_idx];
    }
  }

  //geather into sparse row
  auto flat_index = lu.row_ptr_begin[row]+dense_row_offset;
  for(int i=0; i<numb_elements; i++){
    auto column = sorted_col_ind[i];
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

  //set the element in the column below the pivot 
  auto mult = lu.values[pivot_column_in_other_row]/pivot;

  DenseIndexedRow dense_row; //opt make static, and reset each run
  dense_row.values[pivot_column] = mult;
  dense_row.used_col_ind[0] = pivot_column;
  dense_row.init_columns=1;

  //scatter elements after mult of other_row to a temp row in dense form
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
      if (column > pivot_column){ //TODO is this check needed?
        dense_row.values[column] -= mult*lu.values[k];
        dense_row.used_col_ind[dense_row.init_columns+dense_row.added_columns] = column;
        dense_row.added_columns++;
    }
  }

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
    if (pivot==0) {continue;} //skip lines with pivot value 0

    //for each row below the current pivot
    for(int row = column+1; row<lu.n_rows; row++){
      add_rows(lu, column, row, pivot);
    }
  }
}
