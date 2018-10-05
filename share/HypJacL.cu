#include "HypJacL.hpp"

#include "my_utils.hpp"

unsigned STRAT0 = 0u, STRAT0_STEPS = 0u, STRAT0_PAIRS = 0u;
unsigned STRAT1 = 0u, STRAT1_STEPS = 0u, STRAT1_PAIRS = 0u;

unsigned STRAT0_DTYPE strat0[STRAT0_MAX_STEPS][STRAT0_MAX_PAIRS][2u];
unsigned STRAT1_DTYPE strat1[STRAT1_MAX_STEPS][STRAT1_MAX_PAIRS][2u];

void init_strats(const char *const sdy, const char *const snp0, const unsigned n0, const char *const snp1, const unsigned n1) throw()
{
  if (!sdy || !*sdy)
    DIE("SDY empty");
  if (!snp0 || !*snp0)
    DIE("SNP0 empty");
  if (!n0 || (n0 & 1u))
    DIE("n0 not positive even");
  if (!snp1 || !*snp1)
    DIE("SNP1 empty");
  if (!n1 || (n1 & 1u))
    DIE("n1 not positive even");

  if (!strcmp(snp0, "mmstep")) {
    STRAT0 = STRAT_MMSTEP;
    STRAT0_STEPS = n0;
  }
  else if (!strcmp(snp0, "BrentL")) {
    STRAT0 = STRAT_BRENTL;
    STRAT0_STEPS = n0 - 1u;
  }
  else if (!strcmp(snp0, "colcyc")) {
    STRAT0 = STRAT_COLCYC;
    STRAT0_STEPS = n0 - 1u;
  }
  else if (!strcmp(snp0, "cycloc")) {
    STRAT0 = STRAT_CYCLOC;
    STRAT0_STEPS = n0 - 1u;
  }
  else if (!strcmp(snp0, "rowcyc")) {
    STRAT0 = STRAT_ROWCYC;
    STRAT0_STEPS = n0 - 1u;
  }
  else if (!strcmp(snp0, "cycwor")) {
    STRAT0 = STRAT_CYCWOR;
    STRAT0_STEPS = n0 - 1u;
  }
  else
    DIE("SNP0 unknown");

  STRAT0_PAIRS = (n0 >> 1u);
  (void)memset(strat0, 0, sizeof(strat0));

  if (!strcmp(snp1, "mmstep")) {
    STRAT1 = STRAT_MMSTEP;
    STRAT1_STEPS = n1;
  }
  else if (!strcmp(snp1, "BrentL")) {
    STRAT1 = STRAT_BRENTL;
    STRAT1_STEPS = n1 - 1u;
  }
  else if (!strcmp(snp1, "colcyc")) {
    STRAT1 = STRAT_COLCYC;
    STRAT1_STEPS = n1 - 1u;
  }
  else if (!strcmp(snp1, "cycloc")) {
    STRAT1 = STRAT_CYCLOC;
    STRAT1_STEPS = n1 - 1u;
  }
  else if (!strcmp(snp1, "rowcyc")) {
    STRAT1 = STRAT_ROWCYC;
    STRAT1_STEPS = n1 - 1u;
  }
  else if (!strcmp(snp1, "cycwor")) {
    STRAT1 = STRAT_CYCWOR;
    STRAT1_STEPS = n1 - 1u;
  }
  else
    DIE("SNP1 unknown");

  STRAT1_PAIRS = (n1 >> 1u);
  (void)memset(strat1, 0, sizeof(strat1));

  void *const h = strat_open(sdy);
  SYSP_CALL(h);

  const unsigned char *a0 = reinterpret_cast<const unsigned char*>(strat_ptr(h, snp0, n0));
  if (!a0) {
    SYSI_CALL(strat_close(h));
    DIE("NO SUCH STRATEGY for SNP0 & n0");
  }
  for (unsigned s = 0u; s < STRAT0_STEPS; ++s)
    for (unsigned p = 0u; p < STRAT0_PAIRS; ++p)
      for (unsigned i = 0u; i < 2u; ++i)
        strat0[s][p][i] = *a0++;

  if (n1 <= 256u) {
    const unsigned char *a1 = reinterpret_cast<const unsigned char*>(strat_ptr(h, snp1, n1));
    if (!a1) {
      SYSI_CALL(strat_close(h));
      DIE("NO SUCH STRATEGY for SNP1 & n1");
    }
    for (unsigned s = 0u; s < STRAT1_STEPS; ++s)
      for (unsigned p = 0u; p < STRAT1_PAIRS; ++p)
        for (unsigned i = 0u; i < 2u; ++i)
          strat1[s][p][i] = *a1++;
  }
  else {
    const unsigned short *a1 = reinterpret_cast<const unsigned short*>(strat_ptr(h, snp1, n1));
    if (!a1) {
      SYSI_CALL(strat_close(h));
      DIE("NO SUCH STRATEGY for SNP1 & n1");
    }
    for (unsigned s = 0u; s < STRAT1_STEPS; ++s)
      for (unsigned p = 0u; p < STRAT1_PAIRS; ++p)
        for (unsigned i = 0u; i < 2u; ++i)
          strat1[s][p][i] = *a1++;
  }

  SYSI_CALL(strat_close(h));
}
