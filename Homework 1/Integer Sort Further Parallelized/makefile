PROG = prog
CC = g++
CPPFLAGS = -g -Wall -fopenmp
LINKFLAGS = -fopenmp
OBJS = integersortfurther.o

$(PROG) : $(OBJS)
	$(CC) $(LINKFLAGS) -o $(PROG) $(OBJS)

integersortfurther.o : integersortfurther.cpp
	$(CC) $(CPPFLAGS) -c integersortfurther.cpp

clean:
	rm -f core $(PROG) $(OBJS)

rebuild:
	make clean
	make
