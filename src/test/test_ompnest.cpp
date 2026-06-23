#include <cstdio>
#include <iostream>
#include <omp.h>
using namespace std;

int main(int argc, char **argv)
{
  int n = atoi(argv[1]);
  // omp_set_nested(0);
  omp_set_max_active_levels(1); // max_active_level 0: no para; 1: single layer;
  std::cout << "InPara:" << omp_in_parallel() << ", Level:" << omp_get_level() << ", nThreads:" << omp_get_num_threads()
       << std::endl;
  std::cout << "-----------------------------------\n";
#pragma omp parallel num_threads(2) if (n)
  {
#pragma omp single
    std::cout << "InPara:" << omp_in_parallel() << ", Level:" << omp_get_level() << ", nThreads:" << omp_get_num_threads()
         << std::endl;
#pragma omp parallel num_threads(3)
    {
#pragma omp single
      std::cout << "InPara:" << omp_in_parallel() << ", Level:" << omp_get_level() << ", nThreads:" << omp_get_num_threads()
           << std::endl;
    }
  }

  return 0;
}
