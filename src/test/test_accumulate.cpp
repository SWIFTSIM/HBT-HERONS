#include <iostream>
#include <numeric>
#include <vector>
using namespace std;
int main()
{
  std::vector<int> x(10, 1);
  int s = 0;
  s = accumulate(x.begin(), x.end(), s);
  std::cout << s << std::endl;

  s = accumulate(x.begin() + 1, x.begin() + 5, s);
  std::cout << s << std::endl;

  int init = 100;
  int numbers[] = {10, 20, 30};

  std::cout << "using default accumulate: ";
  std::cout << std::accumulate(numbers, numbers + 3, init);
  std::cout << '\n';

  return 0;
}