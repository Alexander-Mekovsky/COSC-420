PROG = prog
CC = mpic++
CPPFLAGS = -g -Wall
OBJS = ice2.o

$(PROG) : $(OBJS)
	$(CC) -o $(PROG) $(OBJS)

ice2.o : ice2.cpp
	$(CC) $(CPPFLAGS) -c ice2.cpp

clean:
	rm -f core $(PROG) $(OBJS)

rebuild:
	make clean
	make
