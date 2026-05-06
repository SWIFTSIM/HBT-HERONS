#include <assert.h>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <list>
#include <numeric>
#include <sstream>
#include <string>
#include <typeinfo>

#include "../../mymath.h"

void ExchangeAndMerge(MpiWorker_t &world, std::vector<Halo_t> &Halos);
