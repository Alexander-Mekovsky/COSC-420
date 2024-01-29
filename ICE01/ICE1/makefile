PROG = prog
CC = mpic++
CPPFLAGS = -g -Wall
OBJS = ice1.o

$(PROG) : $(OBJS)
	$(CC) -o $(PROG) $(OBJS)

ice1.o : ice1.cpp
	$(CC) $(CPPFLAGS) -c ice1.cpp

clean:
	rm -f core $(PROG) $(OBJS)

rebuild:
	make clean
	make
