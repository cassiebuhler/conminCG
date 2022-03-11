# Makefile for conmin 
# Configuration  Section

# AMPL Solvers directory
S=
# /usr/include directory for libraries
F=
# /usr/include/sys directory for system libraries
G=
# Directory with libf2c
T=

# Change the C compiler as needed
CC = gcc -I$(S) -I$(G) -I$(F) # Must be an ansi compiler

# May need to change the first flag for operating system
OPTFLAGS = -DMacOSX -USTUD -O

AR = ar

ARFLAGS = vq

#RANLIB = ranlib  # if you have ranlib, Uncomment  this line and
#                                       $(RANLIB)       below (for BSD systems)

# Default is to comple conmin only.
default: conmin

# No need to change anything below

.SUFFIXES: .o .c .ln

OBJECTS1 = linalg.o optlp.o iolp.o hash.o cputime.o 

LINTS1 = aconmin.ln linalg.ln optlp.ln iolp.ln hash.ln 

INCL1 = conmin.h myalloc.h hash.h

conmin: aconmin.o version.o libconmin.a 
	$(CC) $(OPTFLAGS) -o conmin aconmin.o version.o $S/funcadd0.o \
	libconmin.a $S/amplsolver.a $T/libf2c.a \
	-lm $(OPTFLAGS2)

libconmin.a: $(OBJECTS1) 
	-rm libconmin.a
	$(AR) $(ARFLAGS)  libconmin.a $(OBJECTS1)

lint: $(LINTS1)
	lint $(LINTS1) -lm | grep -v "op CAST" | grep -v "may sign-extend"

$(OBJECTS1): $(INCL1)

hash.o: hash.h myalloc.h

iolp.o: conmin.h hash.h myalloc.h 

linalg.o: conmin.h myalloc.h

aconmin.o: conmin.h myalloc.h

optlp.o: conmin.h conmin1.h myalloc.h

$(LINTS1): $(INCL1)

.c.o:
	$(CC) $(OPTFLAGS) -c $<

.c.ln:
	lint -c $< | grep -v "op CAST" | grep -v "may sign-extend"

clean: FORCE
	-rm libconmin.a
	-rm conmin conmin2 conmin_p
	-rm *.o
	-rm core

FORCE:

