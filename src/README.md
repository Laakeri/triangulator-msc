## Installing:
1. You need to have CPLEX installed (though the code can be modified to exclude the parts for which CPLEX is needed).
2. Install MaxHS, provided in the MaxHS-3.0 directory. (Probably can be substituted with newer versions of MaxHS available in https://github.com/fbacchus/MaxHS.)
3. Edit the CPLEXDIR and MAXHSDIR in triangulator/Makefile to point to the directories where CPLEX and MaxHS are installed.
4. Use *make* in the triangulator directory.


## Usage:
The first argument denotes the problem. The instance is read from the standard input. See the instances directory for examples of input formats.

Example use:

`./triangulator tw < ../../instances/pace2017_tw_ex/ex002.gr`

Treewidth is denoted with `tw`, minimum fill-in with `minfill`, generalized hypertreewidth with `ghtw`, fractional hypertreewidth with `fhtw`, total table size with `tts`, maximum compatibility (multi-state) with `phyl_maxcomp`, binary maximum compatibility with `phyl_bin_maxcomp`, perfect phylogeny with `phyl_perf` and treelength with `tl`.
The implementation also includes utilities to count the numbers of minimal separators and potential maximal cliques, used with `count_ms` and `count_pmc` arguments (for these the input graph should be connected).


## Installing without CPLEX:
Triangulator can be modified to exclude generalized and fractional hypertreewidth and multi-state maximum compatibility to allow the usage without CPLEX.
To do this, one needs to edit only `Makefile`, `main.cpp`, `phylogeny.cpp` and `phylogeny.hpp`.
All references to `ghtw.*`, `ftw.*`, `hypercost.*` and `maxhs_iface.*` should be removed from these files.
The `CharRemoveMaxSat` functions should be removed from `phylogeny.cpp` and `phylogeny.hpp`.