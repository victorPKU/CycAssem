
SHELL           = /bin/sh

CC              = gcc

CC_LINKERS      = -lm

.SUFFIXES:	.c .o

.c.o:
		$(CC) -c $<


PROGRAMS = fraglink

all:		$(PROGRAMS)
 

fraglink:		fraglink.o superpose.o geometry.o mem.o mem.h superpose.h geometry.h
		$(CC) -o $@ fraglink.o superpose.o geometry.o mem.o $(CC_LINKERS)


clean:
		rm -f *.o core $(PROGRAMS)


fraglink.o:		superpose.h geometry.h
superpose.o:		superpose.h geometry.h mem.h
geometry.o:			geometry.h
mem.o:			mem.h

