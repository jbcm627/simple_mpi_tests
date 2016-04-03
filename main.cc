// c++ includes
#include <iostream>
#include <chrono>
#include <ctime>
#include <string>
// c includes
#include <mpi.h>
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

#define STEPS 100

int main(int argc, char **argv)
{
  // Initialize the MPI environment
  MPI_Init(&argc, &argv);

  int this_rank, max_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &this_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &max_rank);

  // get result from "hostname" in a variable
  FILE *fp = popen("hostname", "r");
  char buf[1024];
  fgets(buf, sizeof(buf), fp);
  pclose(fp);
  std::string hostname(buf);

  // print MPI, OpenMP thread #, host information
  int nthreads, tid, max_threads;
  #pragma omp parallel default(shared) private(nthreads, tid, max_threads)
  {
    tid = omp_get_thread_num();
    nthreads = omp_get_num_threads();
    max_threads = omp_get_max_threads();

    std::cout << " Hello from MPI process "
      + std::to_string(this_rank) + " / " + std::to_string(max_rank)
      + " running on " + hostname
      + "   > I am running in thread " + std::to_string(tid) + " of "
      + std::to_string(nthreads) + " OpenMP threads (with "
      + std::to_string(max_threads) + " threads max)"
      + ".\n" << std::flush;
  }
  
  // Perform some computation in a NX * NY * NZ array
  int NX = 512, NY = 512, NZ = 512;
  // Each MPI proecss only handles a sub-portion of the
  // total array
  NX /= max_rank;
  int *array = new int[NX*NY*NZ];

  // time loop
  auto t_start = std::chrono::high_resolution_clock::now();

  // set array values
  int i,j,k, s;
  for(s=0; s<STEPS; ++s)
  {
    #pragma omp parallel for default(shared) private(i,j,k)
    for(i=0; i<NX; ++i)
      for(j=0; j<NY; ++j)
        for(k=0; k<NZ; ++k)
        {
          array[ i*NY*NZ + j*NZ + k ] += NX;
        }
    // end for
  }

  // print timing information
  auto t_end = std::chrono::high_resolution_clock::now();
  std::cout << "Wall clock time passed on MPI process " + std::to_string(this_rank) + ": "
    + std::to_string( std::chrono::duration<double, std::milli>(t_end-t_start).count() / 1000.0 )
    + "s.\n";

  // Finalize the MPI environment
  MPI_Finalize();
  return 0;
}
