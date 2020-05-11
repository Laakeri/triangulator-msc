CPLEXDIR = ../cplex
MAXHSDIR = ../MaxHS-3.0

CC = g++
CFLAGS = -DIL_STD -std=c++11 -O2 -flto -Wall -Wextra -Wshadow -g -Wfatal-errors -fPIC -march=native
INC = -I $(MAXHSDIR) -I $(CPLEXDIR)/cplex/include -I $(CPLEXDIR)/concert/include
LIB = -lm -lz -lpthread -L $(MAXHSDIR)/build/release/lib -lmaxhs -L $(CPLEXDIR)/cplex/lib/x86-64_linux/static_pic -L $(CPLEXDIR)/concert/lib/x86-64_linux/static_pic -lconcert -lilocplex -lcplex
OBJFILES = treelength.o phylogeny.o maxhs_iface.o ftw.o ghtw.o hypercost.o min-fill.o tts.o bt_algorithm.o treewidth.o enumerator.o hypergraph.o phyl_mat.o graph.o utils.o io.o staticset.o matrix.o mcs.o bitset.o

all: triangulator

triangulator: main.o $(OBJFILES)
	$(CC) $(INC) $(CFLAGS) main.o $(OBJFILES) -o triangulator $(LIB)

# Special cases to suppress warnings in external libraries
maxhs_iface.o: maxhs_iface.cpp
	$(CC) $(INC) $(CFLAGS) -Wno-shadow -MMD -c maxhs_iface.cpp -o maxhs_iface.o

ghtw.o: ghtw.cpp
	$(CC) $(INC) $(CFLAGS) -Wno-ignored-attributes -MMD -c ghtw.cpp -o ghtw.o

ftw.o: ftw.cpp
	$(CC) $(INC) $(CFLAGS) -Wno-ignored-attributes -MMD -c ftw.cpp -o ftw.o

hypercost.o: hypercost.cpp
	$(CC) $(INC) $(CFLAGS) -Wno-ignored-attributes -MMD -c hypercost.cpp -o hypercost.o
# Special cases end

%.o: %.cpp
	$(CC) $(INC) $(CFLAGS) -MMD -c $< -o $@

-include *.d

clean:
	rm -f *.o *.d

.PHONY: all clean
