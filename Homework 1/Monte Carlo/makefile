PROG = prog
CC = g++
CPPFLAGS = -g -Wall -fopenmp
LINKFLAGS = -fopenmp
OBJS = MonteCarlo.o

$(PROG) : $(OBJS)
	$(CC) $(LINKFLAGS) -o $(PROG) $(OBJS)

MonteCarlo.o : MonteCarlo.cpp
	$(CC) $(CPPFLAGS) -c MonteCarlo.cpp

clean:
	rm -f core $(PROG) $(OBJS)

rebuild:
	make clean
	make
