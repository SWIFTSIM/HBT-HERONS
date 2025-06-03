#include <vector>
#include <random>

#include "verify.h"
#include "mymath.h"

int main(int argc, char *argv[])
{
  // Set up repeatable RNG
  std::mt19937 rng;
  rng.seed(0);

  /* How many tests we will do and number of vectors per test. */
  const int nr_reps = 1000;
  const int number_vectors = 1000;

  /* Dimensionality of each vector in a given test. */
  const int min_dimension = 2, max_dimension = 10;
  std::uniform_int_distribution<int> number_dimensions(min_dimension, max_dimension);

  /* Maximum numerical value allowed in an individual vector entry */
  const int max_value = 1000;
  std::uniform_int_distribution<int> random_values(0, max_value);

  for (int rep_nr = 0; rep_nr < nr_reps; rep_nr += 1)
  {
    /* How many dimensions we have in the current test. */
    int vector_dimensions = number_dimensions(rng);

    /* Vector of vectors which we will randomly populate and later sum 
     * elementwise. */
    std::vector<std::vector<HBTInt>> TestVector(number_vectors);
    for(int i = 0; i < number_vectors; i++)
    {
      TestVector[i].reserve(vector_dimensions);
      for(int j = 0; j < vector_dimensions; j++)
        TestVector[i][j] = random_values(rng);
    }

    /* Results based on single thread. */
    std::vector<HBTInt> ResultSingleThread(vector_dimensions,0);
    for(int i = 0; i < number_vectors; i++)
      ResultSingleThread = SumElementwise(ResultSingleThread, TestVector[i]); 

    /* Results based on OMP implementation */
    std::vector<HBTInt> ResultMultipleThread(vector_dimensions,0);
#pragma omp parallel for reduction(SumVectorElementwise : ResultMultipleThread)
    for(int i = 0; i < number_vectors; i++)
      ResultMultipleThread = SumElementwise(ResultMultipleThread, TestVector[i]); 

    for (int i = 0; i < vector_dimensions; i++)
      verify(ResultMultipleThread[i]== ResultSingleThread[i]);
  }
  return 0;
}
