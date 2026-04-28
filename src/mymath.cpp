#include <assert.h>
#include <glob.h>
#include <iostream>

#include "mymath.h"

Timer_t global_timer;

int GetGrid(HBTReal x, HBTReal step, int dim)
{
  int i = floor(x / step);
  if (i < 0)
    i = 0;
  if (i >= dim)
    i = dim - 1;
  return i;
}

int AssignCell(const HBTxyz &Pos, const HBTxyz &step, const vector<int> &dims)
{
#define GRIDtoRank(g0, g1, g2) (((g0) * dims[1] + (g1)) * dims[2] + (g2))
#define GID(i) GetGrid(Pos[i], step[i], dims[i])
  return GRIDtoRank(GID(0), GID(1), GID(2));
#undef GID
#undef GRIDtoRank
}

int count_pattern_files(char *filename_pattern)
{
  glob_t globbuf;

  globbuf.gl_offs = 0;
  glob(filename_pattern, GLOB_ERR, NULL, &globbuf);

  globfree(&globbuf);
  return globbuf.gl_pathc;
}

#define SWAP_2(x) ((((x) & 0xff) << 8) | ((unsigned short)(x) >> 8))
#define SWAP_4(x) (((x) << 24) | (((x) << 8) & 0x00ff0000) | (((x) >> 8) & 0x0000ff00) | ((unsigned)(x) >> 24))
#define FIX_SHORT(x) (*(unsigned short *)&(x) = SWAP_2(*(unsigned short *)&(x)))
#define FIX_LONG(x) (*(unsigned *)&(x) = SWAP_4(*(unsigned *)&(x)))
// bit shift operation is invalid for 8byte+ data

/* This function is used to switch endian, for data2swap[nel] with element size mbyte */
void swap_Nbyte(void *data2swap, size_t nel, size_t mbyte)
{
  size_t i, j;
  char *data, *old_data; // by definition, sizeof(char)=1, one byte

  data = (char *)data2swap;

  switch (mbyte)
  {
  case 1:
    break;
  case 2:
    for (j = 0; j < nel; j++)
      FIX_SHORT(data[j * 2]);
    break;
  case 4:
    for (j = 0; j < nel; j++)
      FIX_LONG(data[j * 4]);
    break;
  default:
    old_data = new char[mbyte];
    for (j = 0; j < nel; j++)
    {
      memcpy(&old_data[0], &data[j * mbyte], mbyte);
      for (i = 0; i < mbyte; i++)
      {
        data[j * mbyte + i] = old_data[mbyte - i - 1];
      }
    }
    delete old_data;
  }
}

/* Find an integer factor of N that is the largest subject to x<=N**(1./dim) */
int LargestRootFactor(int N, int dim)
{
  int x = floor(pow(N, 1. / dim));
  for (; x > 0; x--)
    if (N % x == 0)
      break;
  return x;
}

/* Return a factorization of `N` into `dim` factors that are as close as possible
 * to each other*/
vector<int> ClosestFactors(int N, int dim)
{
  vector<int> factors;
  for (; dim > 0; dim--)
  {
    int x = LargestRootFactor(N, dim);
    factors.push_back(x);
    N /= x;
  }
  return factors;
}

/* Distribute a given number of tasks to available workers as equally as possible,
 * with the leading workers doing one more tasks than others if an equal division
 * of work is not possible. The tasks are assigned to worker_id as [task_begin, task_end),
 * where worker_id is in the range [0, nworkers).*/
void AssignTasks(HBTInt worker_id, HBTInt nworkers, HBTInt ntasks, HBTInt &task_begin, HBTInt &task_end)

{
  /* The (approximate) number of tasks per workers we will have. */
  HBTInt ntask_this = ntasks / nworkers;

  /* Can tasks be split equally across available workers? */
  HBTInt ntask_remainder = ntasks % nworkers;

  /* The first task for each worker. */
  task_begin = ntask_this * worker_id + std::min(ntask_remainder, worker_id);

  /* Leading workers get an extra task. */
  if (worker_id < ntask_remainder)
    ntask_this++;

  /* The last task for each worker. */
  task_end = ntask_this + task_begin;

  assert(task_end <= ntasks);
}

size_t SkipFortranBlock(FILE *fp, bool NeedByteSwap)
{
  int blocksize, blocksize2;
#define ReadBlockSize(a) fread_swap(&a, sizeof(a), 1, fp, NeedByteSwap)
  ReadBlockSize(blocksize);
  fseek(fp, blocksize, SEEK_CUR);
  ReadBlockSize(blocksize2);
  assert(blocksize == blocksize2);
  return blocksize;
#undef ReadBlockSize
}