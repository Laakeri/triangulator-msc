## Installing:
1. You need to have CPLEX installed (though the code can be modified to exclude the parts for which CPLEX is needed).
2. Install MaxHS, provided in the MaxHS-3.0 directory. (Probably can be substituted with newer versions of MaxHS available in https://github.com/fbacchus/MaxHS.)
3. Edit the CPLEXDIR and MAXHSDIR in triangulator/Makefile to point to the directories where CPLEX and MaxHS are installed.
4. Use *make* in the triangulator directory.


## Usage:
./triangulator tw < ../../instances/pace2017_tw_ex/ex002.gr
