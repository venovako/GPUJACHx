#ifndef HYPJACL_HPP
#define HYPJACL_HPP

// `Standard' maximal number of sweeps per Jacobi process
#ifndef HYPJAC_NSWEEP
#define HYPJAC_NSWEEP 30u
#endif // !HYPJAC_NSWEEP

// HEPS = DBL_EPSILON / 2 = 2^(-53) (* Lapack/RN Epsilon *)

// sqrt(HEPS)
#ifndef SQRT_HEPS
#define SQRT_HEPS 1.05367121277235087E-08
#else // SQRT_HEPS
#error SQRT_HEPS not definable externally
#endif // !SQRT_HEPS

// sqrt(2 / HEPS)
#ifndef SQRT_2_HEPS
#define SQRT_2_HEPS 1.34217728000000000E+08
#else // SQRT_2_HEPS
#error SQRT_2_HEPS not definable externally
#endif // !SQRT_2_HEPS

// sqrt(32) * HEPS
#ifndef HYPJAC_MYTOL
#define HYPJAC_MYTOL 6.28036983473510067E-16
#endif // !HYPJAC_MYTOL

// 5 / 4
#ifndef HYPJAC_COTH1_FIX
#define HYPJAC_COTH1_FIX 1.25
#else // HYPJAC_COTH1_FIX
#error HYPJAC_COTH1_FIX not definable externally
#endif // !HYPJAC_COTH1_FIX

// 1 / 2
#ifndef HYPJAC_TANH1_FIX
#define HYPJAC_TANH1_FIX 0.5
#else // HYPJAC_TANH1_FIX
#error HYPJAC_TANH1_FIX not definable externally
#endif // !HYPJAC_TANH1_FIX

// 2 / sqrt(3)
#ifndef HYPJAC_COSHI_FIX
#define HYPJAC_COSHI_FIX 1.15470053837925146
#else // HYPJAC_COSHI_FIX
#error HYPJAC_COSHI_FIX not definable externally
#endif // !HYPJAC_COSHI_FIX

// vpa(DBL_MIN / (1 - HEPS))
#ifndef MU_GAMMA
#define MU_GAMMA 2.2250738585072016301230556379557E-308
#endif // !MU_GAMMA

// vpa(sqrt(MU_GAMMA))
#ifndef MU_
#define MU_ 1.4916681462400414314624091158735E-154
#endif // !MU_

// vpa(32 * (1 + HEPS)^6)
#ifndef DELTA_32
#define DELTA_32 32.000000000000021316282072803011
#endif // !DELTA_32

// vpa(DBL_MAX / DELTA_32)
#ifndef NU_DELTA_32
#define NU_DELTA_32 5.617791046444732845753401614018E+306
#endif // !NU_DELTA_32

// vpa(sqrt(NU_DELTA_32))
#ifndef NU_
#define NU_ 2.370187977027293260705088939015E+153
#endif // !NU_

// MU / DBL_EPSILON
#ifndef MU_EPS
#define MU_EPS 1.0020841800044863889980540256751E-292
#endif // !MU_EPS

// NU / 8
#ifndef NU_8
#define NU_8 2.24711641857789464E+307
#endif // !NU_8

#ifndef HYPJAC_BLK_ORI
#define HYPJAC_BLK_ORI 8u
#else // HYPJAC_BLK_ORI
#error HYPJAC_BLK_ORI not definable externally
#endif // !HYPJAC_BLK_ORI

#ifndef HYPJAC_ACCUMV_2
#define HYPJAC_ACCUMV_2 4u
#else // HYPJAC_ACCUMV_2
#error HYPJAC_ACCUMV_2 not definable externally
#endif // !HYPJAC_ACCUMV_2

#ifndef HYPJAC_FULL_SVD
#define HYPJAC_FULL_SVD 2u
#else // HYPJAC_FULL_SVD
#error HYPJAC_FULL_SVD not definable externally
#endif // !HYPJAC_FULL_SVD

#ifndef HYPJAC_CDSORT_1
#define HYPJAC_CDSORT_1 1u
#else // HYPJAC_CDSORT_1
#error HYPJAC_CDSORT_1 not definable externally
#endif // !HYPJAC_CDSORT_1

#ifndef HYPJACL1_NCOLB
#define HYPJACL1_NCOLB 16u
#else // HYPJACL1_NCOLB
#error HYPJACL1_NCOLB not definable externally
#endif // !HYPJACL1_NCOLB

#ifndef STRAT_MMSTEP
#define STRAT_MMSTEP 1u
#else // STRAT_MMSTEP
#error STRAT_MMSTEP not definable externally
#endif // !STRAT_MMSTEP

#ifndef STRAT_BRENTL
#define STRAT_BRENTL 2u
#else // STRAT_BRENTL
#error STRAT_BRENTL not definable externally
#endif // !STRAT_BRENTL

#ifndef STRAT_COLCYC
#define STRAT_COLCYC 3u
#else // STRAT_COLCYC
#error STRAT_COLCYC not definable externally
#endif // !STRAT_COLCYC

#ifndef STRAT_CYCLOC
#define STRAT_CYCLOC 4u
#else // STRAT_CYCLOC
#error STRAT_CYCLOC not definable externally
#endif // !STRAT_CYCLOC

#ifndef STRAT_ROWCYC
#define STRAT_ROWCYC 5u
#else // STRAT_ROWCYC
#error STRAT_ROWCYC not definable externally
#endif // !STRAT_ROWCYC

#ifndef STRAT_CYCWOR
#define STRAT_CYCWOR 6u
#else // STRAT_CYCWOR
#error STRAT_CYCWOR not definable externally
#endif // !STRAT_CYCWOR

#ifndef STRAT_BLKREC
#define STRAT_BLKREC 7u
#else // STRAT_BLKREC
#error STRAT_BLKREC not definable externally
#endif // !STRAT_BLKREC

// n-1 for cyclic, n for quasi-cyclic
#ifndef STRAT0_MAX_STEPS
#define STRAT0_MAX_STEPS 32u
#endif // !STRAT0_MAX_STEPS
#ifndef STRAT1_MAX_STEPS
#define STRAT1_MAX_STEPS 1024u
#endif // !STRAT1_MAX_STEPS

// n/2 for even n
#ifndef STRAT0_MAX_PAIRS
#define STRAT0_MAX_PAIRS 16u
#endif // !STRAT0_MAX_PAIRS
#ifndef STRAT1_MAX_PAIRS
#define STRAT1_MAX_PAIRS 512u
#endif // !STRAT1_MAX_PAIRS

#ifndef STRAT0_STORAGE
#define STRAT0_STORAGE __constant__
#endif // !STRAT0_STORAGE
#ifndef STRAT1_STORAGE
#define STRAT1_STORAGE __device__
#endif // !STRAT1_STORAGE

#ifndef STRAT0_DTYPE
#define STRAT0_DTYPE char
#endif // !STRAT0_DTYPE
#ifndef STRAT1_DTYPE
#define STRAT1_DTYPE short
#endif // !STRAT1_DTYPE

extern unsigned STRAT0, STRAT0_STEPS, STRAT0_PAIRS;
extern unsigned STRAT1, STRAT1_STEPS, STRAT1_PAIRS;

extern unsigned STRAT0_DTYPE strat0[STRAT0_MAX_STEPS][STRAT0_MAX_PAIRS][2u];
extern unsigned STRAT1_DTYPE strat1[STRAT1_MAX_STEPS][STRAT1_MAX_PAIRS][2u];

extern void init_strats(const char *const sdy, const char *const snp, const unsigned n0, const unsigned n1) throw();

#endif // !HYPJACL_HPP
