PROG = prog
CC = g++
CPPFLAGS = -g -Wall -fopenmp
LINKFLAGS = -fopenmp
OBJS = integersort.o

$(PROG) : $(OBJS)
	$(CC) $(LINKFLAGS) -o $(PROG) $(OBJS)

integersort.o : integersort.cpp
	$(CC) $(CPPFLAGS) -c integersort.cpp

clean:
	rm -f core $(PROG) $(OBJS)

rebuild:
	make clean
	make
