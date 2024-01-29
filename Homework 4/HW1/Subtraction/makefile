PROG = prog
CC = mpic++
CPPFLAGS = -g -Wall
OBJS = subtraction.o

$(PROG) : $(OBJS)
	$(CC) -o $(PROG) $(OBJS)

subtraction.o : subtraction.cpp
	$(CC) $(CPPFLAGS) -c subtraction.cpp

clean:
	rm -f core $(PROG) $(OBJS)

rebuild:
	make clean
	make
