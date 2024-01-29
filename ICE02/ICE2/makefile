PROG = prog
CC = mpic++
CPPFLAGS = -g -Wall
OBJS = ice02.o

$(PROG) : $(OBJS)
	$(CC) -o $(PROG) $(OBJS)

ice02.o : ice02.cpp
	$(CC) $(CPPFLAGS) -c ice02.cpp

clean:
	rm -f core $(PROG) $(OBJS)

rebuild:
	make clean
	make
