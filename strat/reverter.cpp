#include <cstdlib>
#include <iomanip>
#include <iostream>

int main(int argc, char *argv[])
{
  if (argc != 2) {
    std::cerr << *argv << " N" << std::endl;
    return EXIT_FAILURE;
  }
  const int N = atoi(argv[1]);
  if (!N || (N & 1)) {
    std::cerr << argv[1] << " not positive even!" << std::endl;
    return EXIT_FAILURE;
  }
  const int S = N - 1;
  const int P = N / 2;
  const int L = 4 + S + 1;
  int w = 5;
  if (N <= 10)
    w = 1;
  else if (N <= 100)
    w = 2;
  else if (N <= 1000)
    w = 3;
  else if (N <= 10000)
    w = 4;
  const int W = 4 + P * (4 + 2 * w);

  char **lines = new char*[L]();
  lines[0] = new char[33 + 1 + 1]();
  lines[1] = new char[45 + 1 + 1]();
  lines[2] = new char[39 + 1 + 1]();
  lines[3] = new char[1 + 1 + 1]();
  for (int i = 0; i < S; ++i)
    lines[4 + i] = new char[W + 1 + 1]();
  lines[4 + S] = new char[2 + 1 + 1]();

  std::cin.getline(lines[0], 35);
  std::cin.getline(lines[1], 47);
  for (int i = 15; i < 18; ++i) {
    const char tmp = lines[1][i];
    lines[1][i] = lines[1][20 - (i - 15)];
    lines[1][20 - (i - 15)] = tmp;
  }
  std::cin.getline(lines[2], 41);
  std::cin.getline(lines[3], 3);
  for (int i = 0; i < S; ++i)
    std::cin.getline(lines[4 + S - i - 1], W + 2);
  lines[4][W - 1] = ',';
  lines[4 + S - 1][W - 1] = '\0';
  std::cin.getline(lines[4 + S], 4);

  for (int i = 0; i <= 4 + S; ++i)
    std::cout << lines[i] << std::endl;

  for (int i = 4 + S; i >= 0; --i)
    delete [] lines[i];
  delete [] lines;
  return EXIT_SUCCESS;
}
