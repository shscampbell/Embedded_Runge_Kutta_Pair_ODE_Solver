/******************************************************************************
 * Embedded_Runge_Kutta_Pair_ODE_Solver.cpp
 *
 * Implementation of the Embedded_Runge_Kutta_Pair_ODE_Solver class.
 * This class defines the data structures and methods to numerically solve
 * ordinary differential equations. This is an implementation of the Dormand-
 * Prince method, a fifth-order, embedded-pair, adaptive step size Runge-Kutta
 * integrator for systems of ordinary differential equations. The
 * differential equation is expressed in the form of a set of N coupled first-
 * order equations, as follows. Let t be the independent variable. We find a
 * set of dependent variables x1(t), x2(t), ..., xN(t) such that the
 * differential equation can be written as the N coupled equations
 * x1’(t) = F1(t, x1, x2, ..., xN), x2’(t) = F2(t, x1, x2, ..., xN), ...,
 * xN’(t) = FN(t, x1, x2, ..., xN), where x’(t) represents the derivative of x
 * with respect to t. The functions F1, F2, ..., FN are arbitrary real
 * functions of the independent and dependent variables. The fundamental
 * problem is that when given N, the functions F1 to FN, and the value of each
 * dependent variable at some initial value ti of the dependent variable, to
 * find the value of each dependent variable at some final value tf of the
 * dependent variable, to a precision specified by a given tolerance.
 *****************************************************************************/

#define MAG 10.

#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <cerrno>
#include "Embedded_Runge_Kutta_Pair_ODE_Solver.h"
using namespace std;

/************************** CONSTANTS ****************************************/

//These constants are defined by the Dormand-Prince method.
const int Embedded_Runge_Kutta_Pair_ODE_Solver::porder = 5;
const int Embedded_Runge_Kutta_Pair_ODE_Solver::qorder = 4;
const double Embedded_Runge_Kutta_Pair_ODE_Solver::RKc[] =
    {0., 0.2, 0.3, 0.8, 8./9., 1., 1.},
  Embedded_Runge_Kutta_Pair_ODE_Solver::RKb[] = 
    {5179./57600., 0., 7571./16695., 393./640., -92097./339200., 187./2100.,
    1./40.},
  Embedded_Runge_Kutta_Pair_ODE_Solver::RKbhat[] = 
    {35./384., 0., 500./1113., 125./192., -2187./6784., 11./84., 0.},
  Embedded_Runge_Kutta_Pair_ODE_Solver::RKa[][stages-1] = {
    {},
    {1./5.},
    {3./40., 9./40.},
    {44./45., -56./15., 32./9.},
    {19372./6561., -25360./2187., 64448./6561., -212./729.},
    {9017./3168., -355./33., 46732./5247., 49./176., -5103./18656.},
    {35./384., 0., 500./1113., 125./192., -2187./6784., 11./84.}
  };
const double Embedded_Runge_Kutta_Pair_ODE_Solver::RKorder = 5.;

//A constant used for determining step sizes.
const double Embedded_Runge_Kutta_Pair_ODE_Solver::reduc = 0.9;

/************************** UTILITY FUNCTIONS ********************************/

//Increase x by 70%, rounded down.
//Called by the constructors with the growArrays function.
//Arguments:
//  int x: the current integer size.
//Returns:
//  int: the new integer size.
//Assumptions:
//  None
inline int growfac(int x)
{
  return (x*17)/10;
}

//Takes a pointer called arrays pointing to an array of size narrays with
//elements that are sub-arrays of type double, each of length oldsize. The
//arrays are all increased to length newsize with the original input contents
//copied.
//Called by the constructors (unless dokeepsoln==0).
//Arguments:
//  double *&arrays: pointer to an array of double sub-arrays.
//  unsigned int narrays: the number of sub-arrays that are elements of the
//    array pointed to by arrays.
//  unsigned int oldsize: the length of each sub-array.
//  unsigned int newsize: the new length of each sub-array after the function
//    exits.
//Returns:
//  int: newsize, casted to int.
//Assumptions:
//  Assumes newsize is larger than oldsize. The pointer arrays points to an
//  array that contains narrays*oldsize elements.
int growArrays(double *&arrays, unsigned int narrays, 
  unsigned int oldsize, unsigned int newsize)
{
  unsigned int i,j,n,ulim;
  double *newarrays;

  newarrays = new double[narrays*newsize];
  if (!newarrays) {
    cerr << "out of memory" << endl;
    exit(0);
  }
  j=0;
  for (n=0; n<narrays; n++) {
    ulim = n*newsize+oldsize;
    for (i=n*newsize; i<ulim; i++,j++) {
      newarrays[i] = arrays[j];
    }
  }
  delete arrays;
  arrays = newarrays;
  return newsize;
}

//Takes a pointer to a double array called vector of the given length and
//returns the element with the maximum absolute magnitude.
//Called by norm().
//Arguments:
//  const double *vector: pointer to a double array of the given length.
//  unsigned int length: length of the array pointed to by vector.
//Returns:
//  double: the maximum element of the array.
//Assumptions:
//  The pointer vector points to an array of the given length.
inline double maxnorm(const double *vector, unsigned int length)
{
  double norm=fabs(vector[0]), tst;
  unsigned int i=1;

  do {
    if ( (tst = fabs(vector[i++])) > norm )
      norm = tst;
  } while (i<length);
  return norm;
}

//Takes a pointer to a double array called vector of the given length and
//returns the norm of the array. Currently uses the max norm, but alternative
//norms can be implemented and used instead.
//Called by doStep(). Calls maxnorm().
//Arguments:
//  const double *vector: pointer to a double array of the given length.
//  unsigned int length: length of the array pointed to by vector.
//Returns:
//  double: the norm of the array.
//Assumptions:
//  The pointer, vector, points to an array of the given length.
inline double norm(const double *vector, unsigned int length)
{
  return maxnorm(vector, length);
}

//Searches a monotonic double array called array with size elements for the
//index of the elements so that the given value val is between array[index] and
//array[index+1].
//Called by evalX() and evalXAll().
//Arguments:
//  const double *array: pointer to the monotonic double array.
//  double val: the value to locate in the array.
//  unsigned int size: the number of elements in the array.
//Returns:
//  unsigned int: The index of the array such that either val is between
//    array[index] and array[index+1], or the index is 0 or size-1 if val is
//    outside the range of the array.
//Assumptions:
//  The address pointed to by array is expected to be to an array of length
//  size. The elements of the array are assumed to be sorted, either increasing
//  or decreasing with index.
inline unsigned int findIndex(const double *array, double val,
  unsigned int size)
{
  unsigned int li=0, //left index
               ri=size-1, //right index
               mib, //middle index from binary search
               mir; //middle index from linear interpolation search
  double sign = (array[ri]-array[li])>0 ? +1.0 : -1.0;

  if (sign*val<sign*array[li]) return li;
  if (sign*val>sign*array[ri]) return ri;
  while (ri-li>1) {
    mib = (ri+li)/2;
    mir = ( (int) ( (2.0*val-(array[li]+array[ri])) *
        (double)(ri-li)/(array[ri]-array[li]) )
      + ri + li ) / 2;
    if (mib < mir) {
      if (sign*val<sign*array[mib]) {
        ri = mib;
      } else if (sign*val>sign*array[mir]) {
        li = mir;
      } else {
        li = mib;
        ri = mir;
      }
    } else {
      if (sign*val<sign*array[mir]) {
        ri = mir;
      } else if (sign*val>sign*array[mib]) {
        li = mib;
      } else {
        li = mir;
        ri = mib;
      }
    }
  }
  //val is bounded between array[ri] and array[li]. Return the smaller index.
  return li;
}

/************************** PRIVATE METHODS **********************************/

//Carry out the first Runge-Kutta step of the ODE solver, and return the
//determined suitable step size to achieve the desired accuracy of the
//solution.
//Calls doStep, newStepLim. Called by the constructor (see constructor source
//code for the intended use of this method).
//Arguments:
//  double step: The search for a suitable step size begins with the value of
//    this input step size.
//  double *vars: pointer to an array of length 2*(nDepVars+1). vars[0]
//    contains the value of the independent variable where the initial
//    conditions are defined, vars[2*i] for i from 1 to nDepVars contains the
//    initial conditions for the nDepVars dependent variables. When the
//    function returns, vars[1] will contain the new value of the independent
//    variable after the Runge-Kutta step is taken, and vars[2*i+1] will
//    contain the determined new values of the dependent variables.
//  double (**derivFuncs)(double *, double): pointer to an array (of length
//    nDepVars) of function pointers to functions that define the differential
//    equation. See details in the constructor notes.
//  unsigned int nDepVars: the number of dependent variables.
//  double tolerance: The target upper bound on the relative error of the
//    numerical solution.
//  short int flag: A flag for communicating to recursive calls of
//    initializeStep. The initial call uses flag set to 0. A recursive call
//    that increases the step size calls with flag set to 1, and a call that
//    decreases the step size sets flag to -1.
//Returns:
//  double: the value of the step size for the first Runge-Kutta step.
//Assumptions:
//  tolerance is a positive number. Any external (outside of initializeStep)
//  call to initializeStep has the flag argument set to 0. It is assumed that
//  vars points to an array of size 2*(nDepVars+1), and derivFuncs points to an
//  array of size nDepVars. The functions pointed to by each element of
//  derivFuncs each have a first argument that points to an array of nDepVars
//  elements. Each of the functions pointed to by the array derivFuncs need to
//  be real-valued over a relevant range of the independent variable.
//Side Effects:
//  The elements of vars with odd indices will be replaced with the results of
//  the calculation.
double Embedded_Runge_Kutta_Pair_ODE_Solver::initializeStep(double step,
  double *vars, double (**derivFuncs)(double *, double),
  unsigned int nDepVars, double tolerance, short int flag)
{
  double err;

  err = doStep(step, derivFuncs, nDepVars, vars[0], &vars[2], 2, &vars[3], 2);
  if (err < tolerance) {
    //the estimated error is within tolerance.
    if (flag < 0){
      //previous step sizes were too large and we have decreased it to an
      //acceptable value.
      vars[1] = vars[0] + step;
      return newStepLim(step, err, tolerance);
    } else {
      //step size is potentially too small. Try step size an order of 
      //magnitude larger.
      return initializeStep(
        step*MAG, vars, derivFuncs, nDepVars, tolerance, 1);
    }
  } else {
    //the estimated error is not within tolerance.
    if (flag > 0) {
      //previous step size worked but this one is too large.
      step /= MAG;
      err = doStep(
        step, derivFuncs, nDepVars, vars[0], &vars[2], 2, &vars[3], 2);
      vars[1] = vars[0] + step;
      return newStepLim(step, err, tolerance);
    } else {
      //step size is too large.  Try step size an order of magnitude smaller.
      return initializeStep(
        step/MAG, vars, derivFuncs, nDepVars, tolerance, -1);
    }
  }
}

//Calculate a Runge-Kutta step and return the estimated error for the step.
//Calls norm(). Called by initializeStep() and the constructors.
//Arguments:
//  const double step: the step size for the independent variable to take.
//  double (**derivFuncs)(double *, double): pointer to an array (of length
//    nDepVars) of function pointers to functions that define the differential
//    equation. See details in the constructor notes.
//  unsigned int nDepVars: the number of dependent variables.
//  const double ti: the value of the independent variable where the initial
//    conditions are defined.
//  const double *initConds: points to the initial value of the first dependent
//    variable in an array where the initial values of the dependent variables
//    are separated by initCondSep positions in the array.
//  unsigned int initCondSep: the number of array elements separation between
//    successive initial conditions in the array pointed to by initConds. For
//    example, if initial conditions are in successive array elements, then
//    initCondSep==1; if initial conditions are in every other element, then
//    initCondSep==2; if initial conditions are in every third element, then
//    initCondSep==3; and so on.
//  double *xvals: points to the position where the calculated value of the
//    first dependent variable will be stored in an array where the stored
//    values of the dependent variables are separated by xSep positions in the
//    array.
//  unsigned int xSep: the number of array elements separation between
//    dependent variable values for a fixed value of the independent variable.
//Returns:
//  double: the estimated relative error of the Runge-Kutta step.
//Assumptions:
//  The array element pointed to by initConds has at least
//  (nDepVars-1)*initCondSep subsequent elements, and the element pointed to by
//  xvals is followed by at least (nDepVars-1)*xSep elements in the array.
//  derivFuncs points to an array of size nDepVars. The functions pointed to by
//  each element of derivFuncs each have a first argument that points to an
//  array of nDepVars elements. Each of the functions pointed to by the array
//  derivFuncs need to be real-valued at the new independent variable.
//Side Effects:
//  The array elements *(xvals+i*xSep) for i from 0 to nDepVars-1 are modified
//  to be populated with the calculation results.
inline double Embedded_Runge_Kutta_Pair_ODE_Solver::doStep(const double step,
  double (**derivFuncs)(double *, double), unsigned int nDepVars,
  const double ti, const double *initConds, unsigned int initCondSep,
  double *xVals, unsigned int xSep)
{
  int xi, si, si2, xici, xvi;
  double tmid, derivVals[nDepVars][stages], xqvals[nDepVars], delta[nDepVars];

  for (si=0; si<stages; si++) {
    for (xi=0, xici=0; xi<nDepVars; xi++, xici+=initCondSep) {
      xqvals[xi]=0;
      for (si2=0; si2<si; si2++) {
        xqvals[xi] += RKa[si][si2]*derivVals[xi][si2];
      }
      xqvals[xi] = xqvals[xi]*step + initConds[xici];
    }
    tmid = ti + step*RKc[si];
    for (xi=0; xi<nDepVars; xi++) {
      derivVals[xi][si] = (*(derivFuncs[xi]))(xqvals,tmid);
    }
  }
  for (xi=0, xvi=0, xici=0; xi<nDepVars; xi++, xici+=initCondSep, xvi+=xSep) {
    xVals[xvi] = 0;
    xqvals[xi] = 0;
    for (si=0; si<stages; si++) {
      xVals[xvi] += RKbhat[si]*derivVals[xi][si];
      xqvals[xi] += RKb[si]*derivVals[xi][si];
    }
    xVals[xvi] = xVals[xvi]*step + initConds[xici];
    xqvals[xi] = xqvals[xi]*step + initConds[xici];
    //Use the relative error to measure accuracy.
    delta[xi] = xqvals[xi] / xVals[xvi] - 1.;
  }
  return norm(delta, nDepVars);
}

//A quick version of doStep that is intended to be used for calculating a
//single step, and does not estimate the relative error of the calculation.
//This is used for interpolation between previously calculated Runge-Kutta
//steps.
//Called by evalX and evalXAll.
//Arguments:
//  double step: the step size to change the independent variable by in the
//    calculation.
//  double (**derivFuncs)(double *, double): pointer to an array (of length
//    nDepVars) of function pointers to functions that define the differential
//    equation. See details in the constructor notes.
//  unsigned int nDepVars: the number of dependent variables.
//  const double ti: the independent variable at which the given initial
//    conditions are evaluated.
//  const double *initConds: pointer to an array of nDepVars elements
//    containing the initial conditions, the values of the dependent variables
//    for the independent variable at ti.
//  double *xvals: a pointer to an array of nDepVars elements that will be
//    populated with the values of the dependent variables calculated at the
//    new independent variable.
//Returns:
//  None
//Assumptions:
//  initConds and xvals both point to arrays of length nDepVars. derivFuncs
//  points to an array of size nDepVars. The functions pointed to by each
//  element of derivFuncs each have a first argument that points to an array of
//  nDepVars elements. Each of the functions pointed to by the array derivFuncs
//  need to be real-valued at the new independent variable.
//Side Effects:
//  The first nDepVars elements of the array pointed to by xvals are populated
//  with the results of the calculation.
void Embedded_Runge_Kutta_Pair_ODE_Solver::doQuickStep(double step,
  double (**derivFuncs)(double *, double), unsigned int nDepVars,
  const double ti, const double *initConds, double *xVals)
{
  int xi1, si1, si2;
  double tmid, derivVals[nDepVars][stages];

  for (si1=0; si1<stages; si1++) {  
    for (xi1=0; xi1<nDepVars; xi1++) {
      xVals[xi1]=0;
      for (si2=0; si2<si1; si2++) {
        xVals[xi1] += RKa[si1][si2]*derivVals[xi1][si2];
      }
      xVals[xi1] = xVals[xi1]*step + initConds[xi1];
    }
    tmid = ti + step*RKc[si1];
    for (xi1=0; xi1<nDepVars; xi1++) {
      derivVals[xi1][si1] = (*(derivFuncs[xi1]))(xVals,tmid);
    }
  }
  for (xi1=0; xi1<nDepVars; xi1++) {
    xVals[xi1] = 0;
    for (si1=0; si1<stages; si1++) {
      xVals[xi1] += RKbhat[si1]*derivVals[xi1][si1];
    }
    xVals[xi1] = xVals[xi1]*step + initConds[xi1];
  }
}

//The algorithm for determining after a successful Runge-Kutta step the new
//step size for the next step.
//Called by the constructors.
//Arguments:
//  double step: the step size of the previous successful Runge-Kutta step.
//  double err: the estimated relative error of the previous successful
//    Runge-Kutta step.
//  double tolerance: the target upper bound on the relative error of the
//    Runge-Kutta steps.
//Returns:
//  double: the new step size with which to begin the next Runge-Kutta step.
//Assumptions:
//  err and tolerance are positive numbers.
inline double Embedded_Runge_Kutta_Pair_ODE_Solver::newStep(double step,
  double err, double tolerance)
{
  return reduc*step*pow(tolerance/err, 1.0/((double)porder+1.0));
}

//A limited newStep that prevents the step size from increasing by more than a
//factor of 10.
//Called by the constructors and initializeStep(). Calls newStep().
//Arguments:
//  double step: the step size of the previous successful Runge-Kutta step.
//  double err: the estimated relative error of the previous successful
//    Runge-Kutta step.
//  double tolerance: the target upper bound on the relative error of the
//    Runge-Kutta steps.
//Returns:
//  double: the new step size with which to begin the next Runge-Kutta step.
//Assumptions:
//  err and tolerance are positive numbers.
inline double Embedded_Runge_Kutta_Pair_ODE_Solver::newStepLim(double step,
  double err, double tolerance)
{
  double newstep1 = newStep(step, err, tolerance);
  double newstep2 = 10*step;
  return newstep1<newstep2 ? newstep1 : newstep2;
}

//Check if the input value of the independent variable is within the range of
//the internally stored ODE solution in the
//Embedded_Runge_Kutta_Pair_ODE_Solver object. If not, output a warning message
//to standard error.
//Called by evalX() and evalXAll().
//Arguments:
//  double tVal: the value of the independent variable to be checked.
//Returns:
//  None
//Assumptions:
//  None
inline void Embedded_Runge_Kutta_Pair_ODE_Solver::checkIndepRange(double tVal)
{
  double tdiff;
  tdiff = indepVar[1]-indepVar[0];
  if (tdiff > 0) {
    if (tVal < indepVar[0]-tdiff ||
        tVal > 2*indepVar[nSteps-1]-indepVar[nSteps-2]) {
      cerr << "warning: in evalX: t=" << tVal << " is outside the range of "
           << "solution and is not likely to evaluate to within requested "
     << "tolerance." << endl;
    }
  } else /*tdiff < 0*/ {
    if (tVal > indepVar[0]-tdiff ||
        tVal < 2*indepVar[nSteps-1]-indepVar[nSteps-2]) {
      cerr << "warning: in evalX: t=" << tVal << " is outside the range of "
           << "solution and is not likely to evaluate to within requested "
           << "tolerance." << endl;
    }
  }
}

/************************** PUBLIC METHODS ***********************************/

//Perform the integration of the given differential equation with nDepVars
//dependent variables for the independent variable ranging from ti to tf with
//initial conditions specified at ti.
//Calls growfac(), growArrays(), initializeStep(), doStep(), newStep(),
//  newStepLim().
//Arguments:
//  double ti: Initial value of the independent variable.
//  double tf: Final value of the independent variable to integrate to.
//  const double *initConds: pointer to an array of nDepVars elements
//    containing the initial values of each dependent variable, when the
//    independent variable has value ti.
//  double (**derivFuncs)(double *, double): pointer to an array (of nDepVars
//    elements) of function pointers, one for each dependent variable. Each
//    function represents the derivative of the dependent variable with respect
//    to the independent variable, as a function of the independent and
//    dependent variables. Each function has two arguments: a pointer to a
//    double array of length nDepVars with the values of each of the dependent
//    variables, and a double for the value of the independent variable.
//  unsigned int nDepVars: the number of dependent variables.
//  int dokeepsoln: optional argument through overloading. The default behavior
//    (if dokeepsoln is not present, or if it is set to any non-zero integer)
//    is to internally store every determined step of the solution from ti to
//    tf for each dependent variable. If 0 is passed to dokeepsoln, then only
//    the initial and final values of each dependent variable are stored in the
//    class data.
//  double tolerance: optional argument with default value 0.01. The target
//    relative precision of the calculation.
//  double stepSize: optional argument with default value 1. The initial step
//    size of the algorithm.
//  unsigned int initLen: optional argument with default value 0. This is not
//    used if dokeepsoln==0. This specifies the initial size of the data array
//    used to store the value of the independent and dependent variables at
//    each step of the independent variable. If there are more steps than
//    initLen, the data structure will be grown by 70% (an expensive operation)
//    each time the structure fills up. When initLen is 0, a default length is
//    determined based on stepSize and the integration length tf-ti.
//Returns:
//  Embedded_Runge_Kutta_Pair_ODE_Solver: the generated object and its solution
//    data.
//Assumptions:
//  tolerance is a positive number. The pointers initConds and derivFuncs both
//  point to arrays of nDepVars elements. Also, the functions pointed to by
//  each element of derivFuncs each have a first argument that points to an
//  array of nDepVars elements. Note that initLen is ignored if dokeepsoln==0.
//  Each of the functions pointed to by the array derivFuncs need to be real-
//  valued over the given range of the independent variable.
Embedded_Runge_Kutta_Pair_ODE_Solver::Embedded_Runge_Kutta_Pair_ODE_Solver(
  double ti, double tf, const double *initConds,
  double (**derivFuncs)(double *, double), unsigned int nDepVars,
  double tolerance, double stepSize, unsigned int initLen)
{
  //vars is implemented as a two-index array using a single index.
  //indepVar[0..nSteps-1] is being put into vars[0..vlen-1].
  //depVars[i][0..nSteps-1] is being put into vars[vlen*(i+1)..vlen*(i+2)-1]
  double *vars, err;
  int iter, iter2, xiter, vlen=2;

  vars = new double[(nDepVars+1)*vlen];
  if (!vars) {
    //memory allocation failed
    cout << "out of memory" << endl;
    exit(0);
  }
  vars[0] = ti;
  for (iter=0; iter<nDepVars; iter++)
    vars[(iter+1)*vlen] = initConds[iter];
  //carry out the first iteration of the ODE solver and determine the order of
  //magnitude of suitable step size starting at suggested stepSize.
  stepSize = initializeStep(
    stepSize, vars, derivFuncs, nDepVars, tolerance, 0);

  //If the step size is constant, the number of steps needed is
  //  (ti-tf)/stepSize.
  //This is a conservative overestimate to reduce occurances of an expensive
  //regrowth operation of the vars array.
  if (initLen == 0)
    vlen = (int) (7.*abs((tf-ti)/stepSize));
  else
    vlen = initLen;
  growArrays(vars, nDepVars+1, 2, vlen);
  nSteps = 1;
  //vars[nSteps] contains the last calculated value of t. Check if it is
  //greater than tf and then increment nSteps.
  while (vars[nSteps++] < tf) {
    if (nSteps >= vlen)
      vlen = growArrays(vars, nDepVars+1, vlen, growfac(vlen));
    err = doStep(stepSize, derivFuncs, nDepVars, vars[nSteps-1],
      &vars[vlen+nSteps-1], vlen, &vars[vlen+nSteps], vlen);
    while (err > tolerance) {
      stepSize = newStep(stepSize, err, tolerance);
      err = doStep(stepSize, derivFuncs, nDepVars, vars[nSteps-1],
        &vars[vlen+nSteps-1], vlen, &vars[vlen+nSteps], vlen);
    }
    vars[nSteps] = vars[nSteps-1] + stepSize;
    stepSize = newStepLim(stepSize, err, tolerance);
  }

  //Populate the object.
  funcs = derivFuncs;
  nVars = nDepVars;
  indepVar = new double[nSteps];
  depVars = new double *[nDepVars];
  if (!indepVar || !depVars) {
    cout << "out of memory" << endl;
    exit(0);
  }
  for (iter=0; iter<nDepVars; iter++) {
    depVars[iter] = new double[nSteps];
    if (!depVars[iter]) {
      cout << "out of memory" << endl;
      exit(0);
    }
  }
  for (iter=0; iter<nSteps; iter++)
    indepVar[iter] = vars[iter];
  for (iter=0; iter<nDepVars; iter++)
    for (iter2=0, xiter=(iter+1)*vlen; iter2<nSteps; iter2++, xiter++)
      depVars[iter][iter2] = vars[xiter];
  delete vars;
}

Embedded_Runge_Kutta_Pair_ODE_Solver::Embedded_Runge_Kutta_Pair_ODE_Solver(
  double ti, double tf, const double *initConds,
  double (**derivFuncs)(double *, double), unsigned int nDepVars,
  int dokeepsoln, double tolerance, double stepSize, unsigned int initLen)
{
  double *vars, err;
  int iter;

  if (dokeepsoln != 0) {
    Embedded_Runge_Kutta_Pair_ODE_Solver(ti, tf, initConds, derivFuncs,
      nDepVars, tolerance, stepSize, initLen);
    return;
  }
  // The first 2 values of var are of t, the next 2 are of x[0], the next
  // 2 are of x[1], etc.  For each variable pair, the first value is the
  // initial condition and the second is the value after stepSize length
  // integration.
  vars = new double[(nDepVars+1)*2];
  if (!vars) {
    //memory allocation failed
    cout << "out of memory" << endl;
    exit(0);
  }
  // Put initConds in prevvars.
  vars[0] = ti;
  for (iter=0; iter<nDepVars; iter++)
    vars[(iter+1)*2] = initConds[iter];
  stepSize = initializeStep(
    stepSize, vars, derivFuncs, nDepVars, tolerance, 0);
  while (vars[1] < tf) {
    //make the solution the new initial conditions.
    for (iter=0; iter<2*(nDepVars+1); iter+=2)
      vars[iter] = vars[iter+1];
    //do the next step.
    err = doStep(
      stepSize, derivFuncs, nDepVars, vars[0], &vars[2], 2, &vars[3], 2);
    while (err > tolerance) {
      stepSize = newStep(stepSize, err, tolerance);
      err = doStep(
        stepSize, derivFuncs, nDepVars, vars[0], &vars[2], 2, &vars[3], 2);
    }
    vars[1] = vars[0] + stepSize;
    stepSize = newStepLim(stepSize, err, tolerance);
  }
  //Populate the object, keeping only the final 2 points which bind tf.
  nSteps = 2;
  funcs = derivFuncs;
  nVars = nDepVars;
  indepVar = new double[2];
  depVars = new double *[nDepVars];
  if (!indepVar || !depVars) {
    cout << "out of memory" << endl;
    exit(0);
  }
  for (iter=0; iter<nDepVars; iter++) {
    depVars[iter] = new double[2];
    if (!depVars[iter]) {
      cout << "out of memory" << endl;
      exit(0);
    }
  }
  indepVar[0] = vars[0];
  indepVar[1] = vars[1];
  for (iter=0; iter<nDepVars; iter++){
    depVars[iter][0] = vars[2*iter+2];
    depVars[iter][1] = vars[2*iter+3];
  }
  delete vars;
}

Embedded_Runge_Kutta_Pair_ODE_Solver::~Embedded_Runge_Kutta_Pair_ODE_Solver()
{
  unsigned int iter;
  delete indepVar;
  for (iter=0; iter<nVars; iter++) {
    delete depVars[iter];
  }
  delete depVars;
}

//Evaluate and return the dependent variable with index xIndex at the
//independent variable value of tVal. If xIndex is greater or equal to
//nDepVars, the number of dependent variables, then evalX sets errno to ERANGE
//signifying a range error, prints an error message to stderr, and returns 0.
//If tVal is sufficiently outside the independent variable range from ti to tf,
//then the returned result is not likely to be accurate to within tolerance and
//a warning will be printed to stderr.
//Calls checkIndepRange(), findIndex(), doQuickStep().
//Arguments:
//  unsigned int xIndex: The index of the dependent variable to be evaluated.
//  double tVal: The value of the independent variable where the given
//    dependent variable will be evaluated.
//Returns:
//  double: the value of dependent variable (specified by xIndex) evaluated at
//    the independent variable tVal.
//Assumptions:
//  None
double Embedded_Runge_Kutta_Pair_ODE_Solver::evalX(unsigned int xIndex,
  double tVal)
{
  unsigned int iter, titer;
  double xVals[nVars], initConds[nVars];

  if (xIndex >= nVars) {
    cerr << "invalid independent variable " << xIndex << " in call to evalX."
   << endl;
    errno = ERANGE;
    return 0;
  }
  checkIndepRange(tVal);
  //find the index titer of indepVar so that
  //  indepVar[titer] <= tVal <= indepVar[titer+1].
  titer = findIndex(indepVar, tVal, nSteps);
  for (iter=0; iter<nVars; iter++)
    initConds[iter] = depVars[iter][titer];
  doQuickStep(tVal-indepVar[titer], funcs, nVars, tVal, initConds, xVals);
  return xVals[xIndex];
}

//Evaluate the dependent variables at the independent variable value of tVal.
//The results are provided in the array argument xvals, an array with a length
//of the number of dependent variables.
//Calls checkIndepRange(), findIndex(), doQuickStep().
//Arguments:
//  double tVal: The value of the independent variable where the dependent
//    variables will be evaluated.
//  double *xvals: An array that will be populated with the values of the
//    dependent variables evaluated at the independent variable tVal.
//Returns:
//  None
//Assumptions:
//  The length of the array argument xvals needs to be the number of dependent
//    variables.
//Side Effects:
//  The contents of xvals are overwritten with the evaluated variables.
void Embedded_Runge_Kutta_Pair_ODE_Solver::evalXAll(double tVal,
  double *xVals)
{
  unsigned int iter, titer;
  double initConds[nVars];

  checkIndepRange(tVal);
  //find the index titer of indepVar so that
  //  indepVar[titer] <= tVal <= indepVar[titer+1].
  titer = findIndex(indepVar, tVal, nSteps);
  for (iter=0; iter<nVars; iter++)
    initConds[iter] = depVars[iter][titer];
  doQuickStep(tVal-indepVar[titer], funcs, nVars, tVal, initConds, xVals);
}
