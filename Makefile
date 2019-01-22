# Project: test
# Makefile created by Dev-C++ 4.9.9.2

CPP  = g++
CC   = gcc
RES  = 
OBJ  = ODEtest.o Embedded_Runge_Kutta_Pair_ODE_Solver.o $(RES)
LINKOBJ  = ODEtest.o Embedded_Runge_Kutta_Pair_ODE_Solver.o $(RES)
CXXFLAGS = -g 
RM = rm -f

.PHONY: all clean

all: ODEtest


clean:
	${RM} $(OBJ)

ODEtest: $(OBJ)
	$(CPP) $(LINKOBJ) -o ODEtest

ODEtest.o: ODEtest.cpp Embedded_Runge_Kutta_Pair_ODE_Solver.h
	$(CPP) -c ODEtest.cpp -o ODEtest.o $(CXXFLAGS)

Embedded_Runge_Kutta_Pair_ODE_Solver.o: Embedded_Runge_Kutta_Pair_ODE_Solver.cpp   Embedded_Runge_Kutta_Pair_ODE_Solver.h
	$(CPP) -c Embedded_Runge_Kutta_Pair_ODE_Solver.cpp -o Embedded_Runge_Kutta_Pair_ODE_Solver.o $(CXXFLAGS)
