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
  PermutationMatrix p;
  ok = load_matrix_market(argv[1], MAX_N_ELEMENTS, MAX_N_ROWS,
                          nnz, a.n_rows, n_cols,
                          a.values, a.col_ind, a.row_ptr_begin, a.row_ptr_end);
  a.init_memory_management();


  if (!ok){
    fprintf(stderr, "failed to load matrix.\n");
    return -1;
  }

  //solution vectors
  double solution_vector[MAX_N_ROWS];
  //double pattern[] = {1., 1.};
  //double pattern[] = {.1, .1};
  //double pattern[] = {1., -1.};
  //double pattern[] = {5.,-5.};
  double pattern[] = {100.,-100.};
  init_array(solution_vector, a.n_rows, pattern);

  double b[MAX_N_ROWS];
  double c[MAX_N_ROWS];
  matrix_vector_product(a, solution_vector, b);
  //std::cout<<"solution vector: ";
  //print_array(b, a.n_rows);

  /* For debugging, can be removed when implementation is finished. */
  //std::cout<<"nonzeros before:"<<std::endl;
  //dump_nonzeros(a.n_rows, a.values, a.col_ind, a.row_ptr_begin, a.row_ptr_end);

  struct timespec start_time;
  clock_gettime(CLOCK_REALTIME, &start_time);

  /* Perform LU factorization here */
  lu_factorise(a, p);

  struct timespec end_time;
  clock_gettime(CLOCK_REALTIME, &end_time);

  //std::cout<<"nonzeros after:"<<std::endl;
  //dump_nonzeros(a.n_rows, a.values, a.col_ind, a.row_ptr_begin, a.row_ptr_end);

  //print_array(b, a.n_rows, p);
  solve_system(a, p, b, c);
  //check if any elements are wrong
  bool errors = false;
  for(int i=0; i<a.n_rows; i++){
    //dbg(c[i]);
    if (abs(c[i]-solution_vector[i])>0.1){
      errors = true;
      std::cerr<<"INVALID SOLUTION"
               <<" \t\trow: "<<i<<" \t\tcalculated sol:"
               <<c[i]<<" \t\tcorrect sol:"<<solution_vector[i]
               <<std::endl;
    }
  }
  if(errors == false){
    std::cout<<"NO ERRORS ENCOUNTERD, ALL IS WELL"<<std::endl;
  } else {
    dbg("ERROR ENCOUNTERD");
  }

  struct timespec elapsed_time;
  timespec_subtract(&elapsed_time, &end_time, &start_time);

  double elapsed = (double)elapsed_time.tv_sec +
  (double)elapsed_time.tv_nsec / 1000000000.0;
  fprintf(stderr, "elapsed time: %f s\n", elapsed);

  return 0;
}
