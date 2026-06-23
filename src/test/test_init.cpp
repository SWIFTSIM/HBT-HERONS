#include <cstdio>
#include <iostream>
#include <omp.h>
using namespace std;

class myclass
{
public:
  int x, y, z;
  myclass() : x(0), y(1) // z is uninitialized. can be anything.
  {
    printf("default called from %d\n", omp_get_thread_num());
  }
  myclass(int a, int b, int c) : x(a), y(b), z(c)
  {
    printf("initialized from %d\n", omp_get_thread_num());
  }
};
int main()
{
#pragma omp parallel num_threads(2)
  {
    myclass c = myclass(), d;
    myclass a(0, 1, 2); // every thread executes the initializer
  }
  return 0;
}