## use "gcc" to compile source files.
#CC = gcc
##CXX = ~/mpich2-1.0.4p1-icc/bin/mpicxx
##LD = ~/mpich2-1.0.4p1-icc/bin/mpicxx
##LDFLAGS = -L ~/mpich2-1.0.4p1-icc/lib
#CXX = mpicxx
#CXXFLAGS = -g -O0 -w -I${HOME}/tnt
## the linker is also "gcc". It might be something else with other compilers.
#LD = mpicxx
## Compiler flags go here.
#CFLAGS = -g -w
## Linker flags go here. Currently there aren't any, but if we'll switch to
## code optimization, we might add "-s" here to strip debug info and symbols.
#LDFLAGS =
## use this command to erase files.
#RM = /bin/rm -f
## list of generated object files.
#OBJS = mpi_struct.o maxspd.o out.o\
#parboundary.o \
#parupdate.o \
#parflux.o \
#maxspeed.o \
#initialise_maes.o \
#initialise_blast.o \
#initialise_uniform.o \
#vanleer.o 
## program executable file name.
#PROG = ../prog
# top-level rule, to compile everything.

all: $(PROG)

$(PROG): $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o $(PROG)

# now comes a meta-rule for compiling any "C" source file.
%.o: %.cpp *.h
	$(CXX) $(CXXFLAGS) -c $<

clean:
	$(RM) *.o
tidy:
	indent -bli3 -di16 *.cpp

tar:
	tar czvf good$(egg).tgz *.cpp *.h GNUmakefile test.sh ChangeLog
