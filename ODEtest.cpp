/******************************************************************************
 * ODEtest.cpp
 *
 * A program for testing the Embedded_Runge_Kutta_Pair_ODE_Solver class.
 * Solves the differential equation system:
 *   dx/dt = x + y
 *   dy/dt = x - y
 * Outputs the relative error of the generated solutions x(t) and y(t) for
 * different values of t, output in a Mathematica list format, to a file called
 * ODEtest_output.txt.
 *****************************************************************************/
#define OUTFILE "ODEtest_output.txt"

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include "Embedded_Runge_Kutta_Pair_ODE_Solver.h"
using namespace std;

const double sqrt2 = sqrt(2);

//The derivative functions, add and subtract, used in the differential equation
//solver.
double add(double *x, double t)
{
  return x[0]+x[1];
}

double subtract(double *x, double t)
{
  return x[0]-x[1];
}

//The exact solutions, solnx and solny, to the differential equation, used to
//determine the relative error of the numerical solution.
double solnx(double t, double t0, double x0, double y0)
{
  return (x0*cosh(sqrt2*(t-t0)-log(sqrt2-1.0))+y0*sinh(sqrt2*(t-t0)))/sqrt2;
}

double solny(double t, double t0, double x0, double y0)
{
  return (x0*sinh(sqrt2*(t-t0))+y0*cosh(sqrt2*(t-t0)-log(sqrt2+1.0)))/sqrt2;
}

int main()
{
  //The range of the independent variable to integrate over.
  double ti = 0., tf = 15.;
  //The initial conditions, the value of x(ti) and y(ti).
  double initConds[2] = {3.,6.};
  //Define the derivative functions in the differential equation.
  double (*funcs[2])(double *,double) = {&add, &subtract};
  // Numerically solve the differential equation from ti to tf.
  Embedded_Runge_Kutta_Pair_ODE_Solver 
    solver(ti, tf, initConds, funcs, 2, 1e-15, 0.0001);
  fstream output(OUTFILE, fstream::out | fstream::trunc);
  double t, sign = tf-ti>0 ? 1.0 : -1.0, step=sign*0.01;
  unsigned int nSteps = (unsigned int)((tf-ti)/step)+1, i;
  double xvals[nSteps], yvals[nSteps], diffx[nSteps], diffy[nSteps], vals[2];

  //Determine the relative error of the numerical solution.
  for (i=0; i<nSteps; i++) {
    t = ti+(double)i*step;
    solver.evalXAll(t, vals);
    xvals[i] = vals[0];
    yvals[i] = vals[1];
    diffx[i] = solnx(t, ti, initConds[0], initConds[1])/vals[0] - (double)1.0;
    diffy[i] = solny(t, ti, initConds[0], initConds[1])/vals[1] - (double)1.0;
  }

  //Output the results to the OUTFILE.
  output << "diffx:=[";
  for (i=0; i<nSteps; i++) {
    t = ti+(double)i*step;
    output << "[" << t << "," << diffx[i] << "],";
  }
  //back up output position by one to overwrite the last comma.
  output.seekp((long)output.tellp()-1);
  output <<"];"<<endl;
  output << "diffy:=[";
  for (i=0; i<nSteps; i++) {
    t = ti+(double)i*step;
    output <<"["<<t<<","<<diffy[i]<<"],";
  }
  output.seekp((long)output.tellp()-1);
  output <<"];"<<endl;
  output.close();

  return 0;
}