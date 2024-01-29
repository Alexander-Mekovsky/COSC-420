PROG = prog
CC = mpic++
CPPFLAGS = -g -Wall
OBJS = addition.o

$(PROG) : $(OBJS)
	$(CC) -o $(PROG) $(OBJS)

addition.o : addition.cpp
	$(CC) $(CPPFLAGS) -c addition.cpp

clean:
	rm -f core $(PROG) $(OBJS)

rebuild:
	make clean
	make
