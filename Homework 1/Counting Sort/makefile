PROG = prog
CC = g++
CPPFLAGS = -g -Wall -fopenmp
LINKFLAGS = -fopenmp
OBJS = countingsort.o

$(PROG) : $(OBJS)
	$(CC) $(LINKFLAGS) -o $(PROG) $(OBJS)

countingsort.o : countingsort.cpp
	$(CC) $(CPPFLAGS) -c countingsort.cpp

clean:
	rm -f core $(PROG) $(OBJS)

rebuild:
	make clean
	make
