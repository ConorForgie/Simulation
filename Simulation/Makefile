# Makefile for simulation code

#
# C compiler and options for Intel
#
CC=     icpc -O3 -qopenmp -std=c++11
LIB=    -lm

#
# C compiler and options for GNU 
#
#CC=     g++ -O3 -fopenmp -std=c++11
#LIB=	-lm

#
# Object files
#
OBJ=    simulation.o

#
# Compile
#
area:   $(OBJ)
	$(CC) -o $@ $(OBJ) $(LIB)

.c.o:
	$(CC) -c $<

#
# Clean out object files and the executable.
#
clean:
	rm *.o simulation
