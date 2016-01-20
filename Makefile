#CXX=/usr/local/mpich-1.2.6-gcc/bin/mpicxx
CXX=mpicxx
CXXFLAGS=-Wall -g -ansi -pedantic
CXXFLAGS=-Wall -g 
CXXFLAGS=-g 
all:
	${CXX} ${CXXFLAGS} -I${HOME}/tnt/ mpi_struct.cpp \
	maxspd.cpp out.cpp
tidy:
	indent -bli3 -di16 *.cpp
