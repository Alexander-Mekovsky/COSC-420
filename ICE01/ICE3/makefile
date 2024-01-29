PROG = prog
CC = mpic++
CPPFLAGS = -g -Wall
OBJS = ice3.o

$(PROG) : $(OBJS)
	$(CC) -o $(PROG) $(OBJS)

ice3.o : ice3.cpp
	$(CC) $(CPPFLAGS) -c ice3.cpp

clean:
	rm -f core $(PROG) $(OBJS)

rebuild:
	make clean
	make
