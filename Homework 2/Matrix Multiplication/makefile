PROG = prog
CC = g++
CPPFLAGS = -g -Wall -fopenmp
LINKFLAGS = -fopenmp
OBJS = main.o

$(PROG) : $(OBJS)
	$(CC) $(LINKFLAGS) -o $(PROG) $(OBJS)

main.o : main.cpp
	$(CC) $(CPPFLAGS) -c main.cpp

clean:
	rm -f core $(PROG) $(OBJS)

rebuild:
	make clean
	make
