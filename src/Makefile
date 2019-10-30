CC = gcc
CXX = g++

# mcmodel=small
# Generate code for the small code model. The program and its statically defined symbols must be within 4GB of each other. Programs can be statically or dynamically linked. This is the default code model.

# mcmodel=large
# Generate code for the large code model. This makes no assumptions about addresses and sizes of sections. Programs can be statically linked only.

#CFLAGS = -Wall -Wextra -g #-O3 #Optimisations disabled so debugger gets more info
CFLAGS = -Wall -Wextra -g -mcmodel=large #-O3 #Optimisations disabled so debugger gets more info
#CFLAGS = -Wall -Wextra -g -DDBG_MACRO_DISABLE

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
