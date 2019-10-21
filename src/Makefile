CC = gcc
CXX = g++
CFLAGS = -Wall -Wextra -g #-O3 #Optimisations disabled so debugger gets more info

all: lu

lu: main.o lu.o mmio.o matrix.o solve.o
	$(CXX) $(CFLAGS) -o $@ $^ -lrt

mmio.o:	mmio.c #note this is C not C++ 
	$(CC) $(CFLAGS) -c mmio.c

%.o:	%.cc #all other object files are made with C++ compiler
	$(CXX) $(CFLAGS) -c $<

clean:
	rm -f *o
	rm -f lu
