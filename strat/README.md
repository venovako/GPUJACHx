# Adding a new set of parallel Jacobi strategies

## Choose a name for the new set

Make the name 6 characters in length, preferably.
It may be shorter (not tested), but not longer.
For example, let the name be `foobar`.
Follow the rules for naming variables (start with a letter, and continue with alphanumerics).

## Create the strategy tables

For each matrix order, a strategy table has to be created, preferably in a separate file.

### Create a subdirectory and the header files

In `strat` directory, create a subdirectory `foobar`.
That is where the strategy tables should reside.

For each matrix order `N`, for which a table (i.e., a strategy) exists, there should be a file named according to `%05u.h` C-style formatted output of `N`; e.g., for `N == 42`, the corresponding file should be named `00042.h`, and the strategy should be named `foobar00042`.

Each strategy is a 3D array, where the last dimension is `2` (i.e., it is a matrix of pivot pairs).  Each row of the matrix encodes a parallel step of the strategy.  See, e.g., `cycwor` subdirectory for examples.

The strategy table should be declared as
```C
#ifdef USE_STRAT_ARRAY_DECLARATOR
unsigned char  foobar00256[NPSTEP][NPAIRS][2] =
#endif /* USE_STRAT_ARRAY_DECLARATOR */
```
or as
```C
#ifdef USE_STRAT_ARRAY_DECLARATOR
unsigned short foobar01024[NPSTEP][NPAIRS][2] =
#endif /* USE_STRAT_ARRAY_DECLARATOR */
```
where `NPSTEP` is a number of parallel steps (either `N-1` or `N`), and `NPAIRS` is `N/2` (note, only *even* `N` are supported).

The element datatype is either `unsigned char` for `N <= 256`, or `unsigned short` otherwise (under assumption that `N <= 65536`).

## Include the new strategy in the dynamic library build

Depending on the operating system used, in directory `strat` edit `proc_all-lnx.sh` (Linux), `proc_all-mac.sh` (macOS), or `proc_all-win.bat` (Windows) script, and add line `./proc_dir-lnx.sh foobar` (Linux) or `./proc_dir-mac.sh foobar` (macOS) after the other similar invocations, or add `foobar` in the `for` directory list in the batch script (Windows).

### Build the dynamic library of all strategies

Call the appropriate `proc_all-lnx.sh`, `proc_all-mac.sh`, or `proc_all-win.bat`.
If all goes well, there should be a dynamic library `strat.so`, `strat.dylib`, or `strat.dll`, respectively, created in the same directory.

## Register the new strategy in the main code

In directory `shared`, edit `HypJacL.hpp`.
Assign a new ID (`STRAT_FOOBAR`, with an appropriate number, e.g., `8u`) to the strategy, and verify that the limits of `STRAT1_MAX_STEPS` and the corresponding `STRAT1_MAX_PAIRS` are still valid - won't be if the new strategy is defined for `N > 1024`.

Also, edit `HypJacL.cu`, `init_strats` function.
Add another `else if` block that handles `foobar`, i.e., `STRAT_FOOBAR`, and set `STRAT0_STEPS` and `STRAT1_STEPS` to `n0 - 1u` and `n1 - 1u`, respectively, if the strategy has `N-1` parallel steps, or (following `STRAT_MMSTEP`), to `n0` and `n1`, respectively, if the strategy has `N` parallel steps.

### Recompile and test the main code

Note that, as long as the strategy name or maximal supported matrix order does not change, any additions or changes to the new strategy set requires rebuilding of the shared library only.
