#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>

int main(int argc, char *argv[])
{
  if (argc != 2) {
    std::cerr << *argv << " BASENAME" << std::endl;
    return EXIT_FAILURE;
  }

  static const unsigned cell_pt = 30u;
  static const unsigned size_pt = N * cell_pt;

  static const double zero = +0.0;
  static const double one = +1.0;

  // draw images
  for (unsigned s = 0u; s < S; ++s) {
    char fn[19u];
    (void)_snprintf(fn, sizeof(fn), "%s_%02u.asy", argv[1u], s);
    std::ofstream fasy(fn, std::ios_base::out | std::ios_base::trunc);

    const double one_n = one / N;
    const double one_2n = one / (2u * N);
    fasy << std::fixed << std::setprecision(16);

    // header
    fasy << "usepackage(\"amssymb\");" << std::endl;
    fasy << std::endl;
    fasy << "size(" << size_pt << ");" << std::endl;
    fasy << std::endl;
    fasy << "draw(unitsquare);" << std::endl;
    fasy << std::endl;

    // y isolines
    for (unsigned j = N - 1u; j; --j) {
      const double y = one_n * j;
      fasy << "draw((" << zero << ',' << y << ")--(" << one << ',' << y << "));" << std::endl;
    }
    fasy << std::endl;

    // x isolines
    for (unsigned i = N - 1u; i; --i) {
      const double x = one_n * i;
      fasy << "draw((" << x << ',' << zero << ")--(" << x << ',' << one << "));" << std::endl;
    }
    fasy << std::endl;

    // diagonal boxes
    // D---C
    // |   |
    // A---B
    for (unsigned i = 0u; i < N; ++i) {
      const unsigned i1 = i + 1u;
      const unsigned n_i1 = N - i1;
      const double a_x = one_n * i;
      const double a_y = one_n * n_i1;
      const double b_x = a_x + one_n;
      const double b_y = a_y;
      const double c_x = b_x;
      const double c_y = b_y + one_n;
      const double d_x = a_x;
      const double d_y = c_y;
      const double l_x = a_x + one_2n;
      const double l_y = a_y + one_2n;
      fasy << "filldraw((" <<
        a_x << ',' << a_y << ")--(" <<
        b_x << ',' << b_y << ")--(" <<
        c_x << ',' << c_y << ")--(" <<
        d_x << ',' << d_y << ")--cycle,mediumgray);" << std::endl;
      fasy << "label(\"\\huge\\boldmath$" << i << "$\",(" << l_x << ',' << l_y << "),white);" << std::endl;
    }

    // step positions
    for (unsigned i = 0u; i < P; ++i) {
      const unsigned char p = strat[s][i][0u];
      const unsigned char q = strat[s][i][1u];
      double l_x = one_n * q + one_2n;
      double l_y = one - one_n * p - one_2n;
      fasy << std::endl;
      fasy << "label(\"\\huge$\\blacksquare$\",(" << l_x << ',' << l_y << "));" << std::endl;
      l_x = one_n * p + one_2n;
      l_y = one - one_n * q - one_2n;
      fasy << "label(\"\\huge$\\blacksquare$\",(" << l_x << ',' << l_y << "),lightgray);" << std::endl;
    }

    fasy.close();
    char cmd[39u];
    (void)_snprintf(cmd, sizeof(cmd), "asy.exe -noV -f eps %s", fn);
    if (system(cmd))
      return EXIT_FAILURE;
    (void)_snprintf(cmd, sizeof(cmd), "asy.exe -noV -f pdf %s", fn);
    if (system(cmd))
      return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
