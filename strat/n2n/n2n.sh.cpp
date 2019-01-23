#ifndef M
#define M (2u * N)
#endif // !M
#ifndef Sm
#define Sm ((M) - 1u)
#endif // !Sm
#ifndef Pm
#define Pm ((M) / 2u)
#endif // !Pm

static unsigned short stratM[Sm][Pm][C];

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>

int main(int argc, char *argv[])
{
  const unsigned maxw = 5u;
  unsigned w = maxw;
  if (M <= 10u)
    w = 1u;
  else if (M <= 100u)
    w = 2u;
  else if (M <= 1000u)
    w = 3u;
  else if (M <= 10000u)
    w = 4u;

  std::cout << "#ifdef USE_STRAT_ARRAY_DECLARATOR" << std::endl
    << "unsigned " << ((M <= 256u) ? "char  " : "short ")
    << argv[1] << std::setfill('0') << std::setw(maxw) << M << std::setfill(' ')
    << '[' << std::setw(maxw) << Sm << "][" << std::setw(maxw) << Pm << "][2] =" << std::endl
    << "#endif /* USE_STRAT_ARRAY_DECLARATOR */" << std::endl << '{' << std::endl;

  // step 0
  std::cout << "  {";
  for (unsigned j = 0u; j < N; ++j) {
    const unsigned i = 2u * j;
    stratM[0u][j][0u] = i;
    stratM[0u][j][1u] = i + 1u;
    std::cout << '{' << std::setw(w) << stratM[0u][j][0u] << ',' << std::setw(w) << stratM[0u][j][1u] << '}';
    if (j < (N - 1u))
      std::cout << ',';
  }
  std::cout << "}," << std::endl;

  if (!strcmp(argv[1], "colcyc")) {
    for (unsigned p = 0u; p < Sn; ++p) {
      const unsigned q = 2u * p + 1u;
      for (unsigned a = 0u; a < Pn; ++a) {
        const unsigned b = 2u * a;
        const unsigned i = stratN[p][a][0u];
        const unsigned j = stratN[p][a][1u];
        // NW
        stratM[q][b][0u] = 2u * i;
        stratM[q][b][1u] = 2u * j;
        // SE
        stratM[q][b + 1u][0u] = 2u * i + 1u;
        stratM[q][b + 1u][1u] = 2u * j + 1u;
        // SW
        stratM[q + 1u][b][0u] = 2u * i + 1u;
        stratM[q + 1u][b][1u] = 2u * j;
        // NE
        stratM[q + 1u][b + 1u][0u] = 2u * i;
        stratM[q + 1u][b + 1u][1u] = 2u * j + 1u;
      }
      std::cout << "  {";
      for (unsigned b = 0u; b < Pm; ++b) {
        std::cout << '{' << std::setw(w) << stratM[q][b][0u] << ',' << std::setw(w) << stratM[q][b][1u] << '}';
        if (b < (Pm - 1u))
          std::cout << ',';
      }
      std::cout << "}," << std::endl;
      std::cout << "  {";
      for (unsigned b = 0u; b < Pm; ++b) {
        std::cout << '{' << std::setw(w) << stratM[q + 1u][b][0u] << ',' << std::setw(w) << stratM[q + 1u][b][1u] << '}';
        if (b < (Pm - 1u))
          std::cout << ',';
      }
      std::cout << '}';
      if (p < (Sn - 1u))
        std::cout << ',';
      std::cout << std::endl;
    }
  }
  else { // rowcyc
    for (unsigned p = 0u; p < Sn; ++p) {
      const unsigned q = 2u * p + 1u;
      for (unsigned a = 0u; a < Pn; ++a) {
        const unsigned b = 2u * a;
        const unsigned i = stratN[p][a][0u];
        const unsigned j = stratN[p][a][1u];
        // NW
        stratM[q][b][0u] = 2u * i;
        stratM[q][b][1u] = 2u * j;
        // SE
        stratM[q][b + 1u][0u] = 2u * i + 1u;
        stratM[q][b + 1u][1u] = 2u * j + 1u;
        // NE
        stratM[q + 1u][b][0u] = 2u * i;
        stratM[q + 1u][b][1u] = 2u * j + 1u;
        // SW
        stratM[q + 1u][b + 1u][0u] = 2u * i + 1u;
        stratM[q + 1u][b + 1u][1u] = 2u * j;
      }
      std::cout << "  {";
      for (unsigned b = 0u; b < Pm; ++b) {
        std::cout << '{' << std::setw(w) << stratM[q][b][0u] << ',' << std::setw(w) << stratM[q][b][1u] << '}';
        if (b < (Pm - 1u))
          std::cout << ',';
      }
      std::cout << "}," << std::endl;
      std::cout << "  {";
      for (unsigned b = 0u; b < Pm; ++b) {
        std::cout << '{' << std::setw(w) << stratM[q + 1u][b][0u] << ',' << std::setw(w) << stratM[q + 1u][b][1u] << '}';
        if (b < (Pm - 1u))
          std::cout << ',';
      }
      std::cout << '}';
      if (p < (Sn - 1u))
        std::cout << ',';
      std::cout << std::endl;
    }
  }

  std::cout << "};" << std::endl;
  return EXIT_SUCCESS;
}
