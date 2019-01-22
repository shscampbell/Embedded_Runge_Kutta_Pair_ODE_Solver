# Embedded_Runge_Kutta_Pair_ODE_Solver

This package defines a C++ class for numerical evaluation of ordinary differential equations. It implements the [Dormand-Prince](https://en.wikipedia.org/wiki/Dormand%E2%80%93Prince_method), a fifth-order, embedded-pair, adaptive step size Runge-Kutta integrator for systems of ordinary differential equations.

# Motivation

Any system of ordinary differential equations can be written as a system of first order equations of the form
![ODE](https://latex.codecogs.com/svg.latex?%5Clarge%20%5Cfrac%7Bdx_1%7D%7Bdt%7D%3DF_1%28t%2Cx_1%2Cx_2%2C%5Cdots%2Cx_N%29%5C%5C%5C%5C%5Cmbox%7B%5C%20%5C%20%5C%20%5C%2C%7D%5Cfrac%7Bdx_2%7D%7Bdt%7D%3DF_2%28t%2Cx_1%2Cx_2%2C%5Cdots%2Cx_N%29%5C%5C%5Cmbox%7B%5C%20%5C%20%5C%20%5C%20%5C%20%5C%20%5C%20%5C%20%5C%20%5C%20%7D%5Cvdots%5C%5C%5Cmbox%7B%5C%20%5C%20%5C%2C%7D%5Cfrac%7Bdx_N%7D%7Bdt%7D%3DF_N%28t%2Cx_1%2Cx_2%2C%5Cdots%2Cx_N%29).
Here, ![t](https://latex.codecogs.com/svg.latex?t) is called the *independent variable*, and the ![N](https://latex.codecogs.com/svg.latex?N) variables ![xi](https://latex.codecogs.com/svg.latex?x_i) are called the dependent variables.

The purpose is to determine dependent variables ![xi(t)](https://latex.codecogs.com/svg.latex?x_i%28t%29) that satisfy a given set of differential equations over a range of the independent variable ![range](https://latex.codecogs.com/svg.latex?t_1%5Cleq%20t%5Cleq%20t_2) for a particular set of initial conditions ![initConds](https://latex.codecogs.com/svg.latex?x_i%28t_1%29).

A sample code demonstrating the use of the class is provided in the file ODEtest.cpp. It evaluates the functions ![x(t)](https://latex.codecogs.com/svg.latex?x%28t%29) and ![y(t)](https://latex.codecogs.com/svg.latex?y%28t%29) satisfying the differential equations
![ODEtest](https://latex.codecogs.com/svg.latex?%5Cfrac%7Bdx%7D%7Bdt%7D%3Dx+y%5C%5C%5C%5C%5Cmbox%7B%5C%20%5C%20%5C%20%5C%20%5C%2C%7D%5Cfrac%7Bdy%7D%7Bdt%7D%3Dx-y)
over the range ![ODErange](https://latex.codecogs.com/svg.latex?0%5Cleq%20t%5Cleq%2015) where ![xinit](https://latex.codecogs.com/svg.latex?x%280%29%3D3) and ![yinit](https://latex.codecogs.com/svg.latex?y%280%29%3D6). The result is directly compared to the exact solution
![soln](https://latex.codecogs.com/svg.latex?x%28t%29%3D%5B3%5Ccosh%28%5Csqrt%7B2%7Dt-%5Cln%28%5Csqrt%7B2%7D-1%29%29+6%5Csinh%28%5Csqrt%7B2%7Dt%29%5D/%5Csqrt%7B2%7D%5C%5C%5C%5C%5Cmbox%7B%5C%20%5C%20%5C%20%5C%20%5C%2C%7Dy%28t%29%3D%5B3%5Csinh%28%5Csqrt%7B2%7Dt%29+6%5Ccosh%28%5Csqrt%7B2%7Dt-%5Cln%28%5Csqrt%7B2%7D+1%29%29%5D/%5Csqrt%7B2%7D)

# Language

C++

# Installation

To compile and run the test program, open the makefile and set the C++ compiler command and flags to the desired settings for your system. Run `make` to compile the source, then run the compiled program `ODEtest`. This will generate the output file ODEtest_output.txt, which can be compared to the provided file in this repo.

# API Reference

### `Embedded_Runge_Kutta_Pair_ODE_Solver::Embedded_Runge_Kutta_Pair_ODE_Solver(double ti, double tf,const double *initConds, double (**derivFuncs)(double *, double), unsigned int nDepVars, double tolerance = 0.01, double stepSize = 1., unsigned int initLen = 0)`

### `Embedded_Runge_Kutta_Pair_ODE_Solver::Embedded_Runge_Kutta_Pair_ODE_Solver(double ti, double tf, const double *initConds, double (**derivFuncs)(double *, double), unsigned int nDepVars, int dokeepsoln, double tolerance = 0.01, double stepSize = 1., unsigned int initLen = 0)`

Perform the integration of the given differential equation with _nDepVars_ dependent variables for the independent variable ranging from _ti_ to _tf_ with initial conditions specified at _ti_.
Arguments:
	`double ti`: Initial value of the independent variable.
	`double tf`: Final value of the independent variable to integrate to.
	`const double *initConds`: pointer to an array of _nDepVars_ elements containing the initial values of each dependent variable, when the independent variable has value _ti_.
	`double (**derivFuncs)(double *, double)`: pointer to an array (of _nDepVars_ elements) of function pointers, one for each dependent variable. Each function represents the derivative of the dependent variable with respect to the independent variable, as a function of the independent and dependent variables. Each function has two arguments: a pointer to a double array of length _nDepVars_ with the values of each of the dependent variables, and a double for the value of the independent variable.
	`int nDepVars`: the number of dependent variables.
	`int dokeepsoln`: optional argument through overloading. The default behavior (if _dokeepsoln_ is not present, or if it is set to any non-zero integer) is to internally store every determined step of the solution from _ti_ to _tf_ for each dependent variable. If 0 is passed to _dokeepsoln_, then only the initial and final values of each dependent variable are stored in the class data.
	`double tolerance`: optional argument with default value 0.01. The target relative precision of the calculation.
	`double stepSize`: optional argument with default value 1. The initial step size of the algorithm.
	`unsigned int initLen`: optional argument with default value 0. This is not used if `dokeepsoln==0`. This specifies the initial size of the data array used to store the value of the independent and dependent variables at each step of the independent variable. If there are more steps than _initLen_, the data structure will be grown by 70% (an expensive operation) each time the structure fills up. When _initLen_ is 0, a default length is determined based on _stepSize_ and the integration length `tf-ti`.
Returns:
	`Embedded_Runge_Kutta_Pair_ODE_Solver`: the generated object and its solution data.
Assumptions:
	_tolerance_ is a positive number. The pointers _initConds_ and _derivFuncs_ both point to arrays of _nDepVars_ elements. Also, the functions pointed to by each element of _derivFuncs_ each have a first argument that points to an array of _nDepVars_ elements. Note that _initLen_ is ignored if `dokeepsoln==0`. Each of the functions pointed to by the array _derivFuncs_ need to be real-valued over the given range of the independent variable.

### `Embedded_Runge_Kutta_Pair_ODE_Solver::~Embedded_Runge_Kutta_Pair_ODE_Solver()`

### `double Embedded_Runge_Kutta_Pair_ODE_Solver::evalX(unsigned int xIndex, double tVal)`

Evaluate and return the dependent variable with index _xIndex_ at the independent variable value of _tVal_. If _xIndex_ is greater or equal to _nDepVars_, the number of dependent variables, then _evalX_ sets _errno_ to _ERANGE_ signifying a range error, prints an error message to _stderr_, and returns 0. If _tVal_ is sufficiently outside the independent variable range from _ti_ to _tf_, then the returned result is not likely to be accurate to within _tolerance_ and a warning will be printed to _stderr_.
Arguments:
	`unsigned int xIndex`: The index of the dependent variable to be evaluated.
	`double tVal`: The value of the independent variable where the given dependent variable will be evaluated.
Returns:
	`double`: the value of dependent variable (specified by _xIndex_) evaluated at the independent variable _tVal_.
Assumptions:
	None

###`void Embedded_Runge_Kutta_Pair_ODE_Solver::evalXAll(double tVal, double *xVals)`

Evaluate the dependent variables at the independent variable value of _tVal_. The results are provided in the array argument _xvals_, an array with a length of the number of dependent variables.
Arguments:
	`double tVal`: The value of the independent variable where the dependent variables will be evaluated.
	`double *xvals`: An array that will be populated with the values of the dependent variables evaluated at the independent variable _tVal_.
Returns:
	None
Assumptions:
	The length of the array argument _xvals_ needs to be the number of dependent variables.
Side Effects:
	The contents of _xvals_ are overwritten with the evaluated variables.

# Credits

The derivation of the algorithm for the Dormand-Prince method, and other similar methods, is found in the book [Numerical Methods for Differential Equations: A Computational Approach](https://www.amazon.com/Numerical-Methods-Differential-Equations-Computational/dp/0849394333/ref=sr_1_3?keywords=Numerical+Methods+for+Differential+Equations%3A+A+Computational+Approach&qid=1548139872&sr=8-3) by John R. Dormand.

# License

GNU General Public License v3

GNU © [Sheldon Campbell](https://www.github.com/shscampbell/)