PROG = prog
CC = mpic++
CPPFLAGS = -g -Wall
OBJS = ice03.o

$(PROG) : $(OBJS)
	$(CC) -o $(PROG) $(OBJS)

ice03.o : ice03.cpp
	$(CC) $(CPPFLAGS) -c ice03.cpp

clean:
	rm -f core $(PROG) $(OBJS)

rebuild:
	make clean
	make
