#include "lu.h"
#include "solve.h"
#include <cmath>

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

int
main(int argc, char **argv)
{
  if (argc != 2){
    fprintf(stderr, "usage: %s <filename>\n", argv[0]);
    return -1;
  }

  int nnz, n_cols;
  bool ok(false);
  //must be static to make sure it is not allocated on the heap
  static CompressedRowMatrix a;
  ok = load_matrix_market(argv[1], MAX_N_ELEMENTS, MAX_N_ROWS,
                          nnz, a.n_rows, n_cols,
                          a.values, a.col_ind, a.row_ptr_begin, a.row_ptr_end);
  a.init_memory_management();

  if (!ok){
    fprintf(stderr, "failed to load matrix.\n");
    return -1;
  }
/////////////////////////////////////////////////////////////////////

  double pattern[5][2] = {
    {1., 1.},
    {.1, .1},
    {1., -1.},
    {5.,-5.},
    {100.,-100.},
  };

  double solution_vector[5][MAX_N_ROWS];
  double b_vec[5][MAX_N_ROWS];

  for(int i=0; i<5; i++){
    init_array(solution_vector[i], a.n_rows, pattern[i]);
    matrix_vector_product(a, solution_vector[i], b_vec[i]);
  }

//////////////////////////////////////////////////////////////////////

    struct timespec start_time;
    clock_gettime(CLOCK_REALTIME, &start_time);

    /* Perform LU factorization here */
    static PermutationMatrix p;
    lu_factorise(a, p);

    struct timespec end_time;
    clock_gettime(CLOCK_REALTIME, &end_time);

    //store elapsed time
    struct timespec elapsed_time;
    timespec_subtract(&elapsed_time, &end_time, &start_time);
    double factor_times = (double)elapsed_time.tv_sec
     + (double)elapsed_time.tv_nsec / 1000000000.0;  
    fprintf(stderr, "factoring took %f s\n", factor_times);

//////////////////////////////////////////////////////////////////////

  double x[MAX_N_ROWS];
  double elapsed_times[5];
  bool errors[5] = {false};
  double relative_error[5] = {0.0};
  bool any_error = false;

  for(int i=0; i<5; i++){
    struct timespec start_time;
    clock_gettime(CLOCK_REALTIME, &start_time);

    solve_system(a, p, b_vec[i], x);

    struct timespec end_time;
    clock_gettime(CLOCK_REALTIME, &end_time);

    //store elapsed time
    struct timespec elapsed_time;
    timespec_subtract(&elapsed_time, &end_time, &start_time);
    elapsed_times[i] = (double)elapsed_time.tv_sec
     + (double)elapsed_time.tv_nsec / 1000000000.0;

    //check if any errors where made
    double normed_error_size = 0;
    double normed_actual_solution = 0;
    for(int j=0; j<a.n_rows; j++){
      //TODO replace with euclidian error as in assignment
      normed_error_size += std::pow(x[j]-solution_vector[i][j], 2.0);
      normed_actual_solution += std::pow(x[j], 2.0);
      if (std::sqrt(normed_error_size)>0.1 || !std::isfinite(x[j])){
        any_error = true;
        /*std::cerr<<"INVALID SOLUTION"
                <<" \t\trow: "<<j<<" \t\tcalculated sol:"
                <<x[i]<<" \t\tcorrect sol:"<<+solution_vector[i]
                <<std::endl;*/
      }
    }
    relative_error[i] = std::sqrt(normed_error_size)/std::sqrt(normed_actual_solution);
  }

  for(int i=0; i<5; i++){
    if(any_error){
      std::cout<<"ERROR ENCOUNTERD"<<std::endl;
    } 
    fprintf(stderr, "pattern %d: %fs, relative error: %f\n", 
            i, elapsed_times[i], relative_error[i]);
    }
  std::cout<<std::endl;
  return 0;
}
