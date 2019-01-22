#ifndef EMBEDDED_RUNGE_KUTTA_PAIR_ODE_SOLVER_H
#define EMBEDDED_RUNGE_KUTTA_PAIR_ODE_SOLVER_H

#include <cmath>
#include <cerrno>
#include <iostream>
using namespace std;

/*Numerically solve the coupled set of 1st-order ODEs dx/dt=F(x,t), 
  where x=(x1,x2,...,xN), F(x)=(F1(x,t),F2(x,t),...,FN(x,t)),
  for the functions x1(t),x2(t),...,xN(t).*/
class Embedded_Runge_Kutta_Pair_ODE_Solver {

  public:
    //Constructor arguments:
    //  ti,tf: initial and final values of independent variable.
    //  initConds: array of length nDepVars containing initial values of
    //    dependent variables (at ti).
    //  derivFuncs: array of length nDepVars of pointers to derivative
    //    functions that define the system of 1st order ODEs.
    //  nDepVars: the number of dependent variables.
    //Optional arguments with default value:
    //  tolerance: an indicator of the relative precision of the solved
    //    solution (default 0.01).
    //  stepSize: the initial step size the solver should attempt (default 1).
    //  initLen: the initial length of the solution arrays that the computer
    //    will allocate. If 0, the program estimates the needed length using a
    //    heuristic based on the initial step size used. (default 0)
    //Optional arguments (via overloading):
    //  dokeepsoln: an integer which when set to 0 causes the solution to the
    //    ODE integration to not be saved in depVars.  Only the determined
    //    values at tf are saved.  This can save a lot of memory usage if one
    //    only needs the solution evaluated at a single position, a finite
    //    number of ordered positions (continue the integration using the last
    //    result as the new initial conditions), or only wants the solution for
    //    a limited range of t (integrate up to the lower bound, then call
    //    again with the lower bound as the new initial condition, this time
    //    saving the solution).
    Embedded_Runge_Kutta_Pair_ODE_Solver(double ti, double tf,
      const double *initConds, double (**derivFuncs)(double *, double),
      unsigned int nDepVars, double tolerance = 0.01, double stepSize = 1.,
      unsigned int initLen = 0);
    Embedded_Runge_Kutta_Pair_ODE_Solver(double ti, double tf,
      const double *initConds, double (**derivFuncs)(double *, double),
      unsigned int nDepVars, int dokeepsoln, double tolerance = 0.01,
      double stepSize = 1., unsigned int initLen = 0);

    //destructor
    ~Embedded_Runge_Kutta_Pair_ODE_Solver();

    //Evaluate the dependent variable depVars[xIndex] at the independent
    //variable value of tVal.  If xIndex>=nDepVars then evalX sets errno to
    //ERANGE signifying a range error, prints an error message to stderr, and
    //returns 0. If tVal is sufficiently outside the range ti..tf, then the
    //returned result is not likely to be accurate within suggested tolerance
    //and a warning will be printed to stderr.
    double evalX(unsigned int xIndex, double tVal);

    //As with evalX except it populates xVals, an array of size nDepVars
    //containing the values of all nDepVars dependent variables evaluated at
    //tVal.
    void evalXAll(double tVal, double *xVals);


  protected:
    //Runge-Kutta ODE solver constants and variables
    static const int porder;
    static const int qorder;
    static const int stages = 7;
    static const double RKc[stages],
                        RKb[stages],
      RKbhat[stages],
      RKa[stages][stages-1];
    static const double RKorder, reduc;

    //The number of dependent variables
    unsigned int nVars;

    //The definition of the ODE via an array of pointers to functions that have
    //arguments of the nDepVars x-values (in an array) and the t-value.
    double (**funcs)(double *, double);

    //Arrays that contain the solution to the ODE, each of length nSteps
    double *indepVar, **depVars;
    unsigned int nSteps;

  private:
    double initializeStep(double step, double *vars,
      double (**derivFuncs)(double *, double), unsigned int nDepVars,
      double tolerance, short int flag);

    double doStep(double step, double (**derivFuncs)(double *, double),
      unsigned int nDepVars, const double ti, const double *initConds,
      unsigned int initCondSep, double *xVals, unsigned int xSep);

    void doQuickStep(double stepSize, double (**derivFuncs)(double *, double),
      unsigned int nDepVars, const double ti, const double *initConds,
      double *xVals);

    double newStep(double step, double err, double tolerance);

    double newStepLim(double step, double err, double tolerance);

    void checkIndepRange(double tVal);
};

#endif
