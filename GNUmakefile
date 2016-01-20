
# use "gcc" to compile source files.
export CC = gcc
#CXX = ~/mpich2-1.0.4p1-icc/bin/mpicxx
#LD = ~/mpich2-1.0.4p1-icc/bin/mpicxx
#LDFLAGS = -L ~/mpich2-1.0.4p1-icc/lib
#export MPICH_CXX = g++
export CXX = mpicxx
#export CXX = ~/mpich2-1.0.3-gcc/bin/mpicxx
export CXXFLAGS = -g -O0 -w -I${HOME}/tnt -I${HOME}/Documents/hdf5/include
# the linker is also "gcc". It might be something else with other compilers.
export LD = mpicxx
#export LD = ~/mpich2-1.0.3-gcc/bin/mpicxx
# Compiler flags go here.
export CFLAGS = -g -w -I${HOME}//Documents/hdf5/include
# Linker flags go here. Currently there aren't any, but if we'll switch to
# code optimization, we might add "-s" here to strip debug info and symbols.
export LDFLAGS =-L${HOME}/Documents/hdf5/lib -lhdf5
# use this command to erase files.
export RM = /bin/rm -f
# list of generated object files.
export OBJS = main.o  out.o\
parboundary.o \
parupdate.o \
parflux.o \
parsecond.o \
maxspeed.o \
riemann.o \
hlld.o \
sgn.o \
ctop.o \
lf.o \
outhdf5.o \
eigenvectors.o \
emf_exchange.o \
maes.o \
blast.o \
initialise_uniform.o \
cooling.o \
molcool.o \
tabfind.o \
locate.o \
orszagtang.o \
zanni.o \
zannisimple.o \
initialise_jet.o \
rungekutta.o \
minmod.o \
vanleer.o
# program executable file name.
export PROG = ../MAES

all:
	make -C src
clean:
	make -C src clean
tar:
	rm -rf out_*
	make -C src clean
	tar czvf good.tgz  src GNUmakefile test.sh
