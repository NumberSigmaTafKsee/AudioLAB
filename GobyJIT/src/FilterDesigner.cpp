// BEGIN FINAL FILE: FilterDesign.cpp
/*
 *
 *  Filter Design Utility
 *  Source
 *
 *  Brent Lehman
 *  16 March 2007
 *
 * DESolver.hpp
 * DESolver.cpp
 */


////////////////////////////////////////////////////////////////////
//                                                                //
//  The idea is that an optimization algorithm passes a bunch of  //
//  different filter specifications to the function               //
//  "EnergyFunction" below.  That function is supposed to         //
//  compute an "error" or "cost" value for each specification     //
//  it receives, which the algorithm uses to decide on other      //
//  filter specifications to try.  Over the course of several     //
//  thousand different specifications, the algorithm will         //
//  eventually converge on a single best one.  This one has the   //
//  lowest error value of all possible specifications.  Thus,     //
//  you effectively tell the optimization algorithm what it's     //
//  looking for through code that you put into EnergyFunction.    //
//                                                                //
//  Look for a note in the code like this one to see what part    //
//  you need to change for your own uses.                         //
//                                                                //
////////////////////////////////////////////////////////////////////


#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <conio.h>
#include <math.h>
#include <time.h>
#include "DESolver.h"


#define kIntIsOdd(x) (((x) & 0x00000001) == 1)


double FilterSolver::EnergyFunction(double testSolution[], bool& bAtSolution)
{
  unsigned i;
  double   tempReal;
  double   tempImag;

  // You probably will want to keep this if statement and its contents
  if (mMinimumPhase)
  {
    // Make sure there are no zeros outside the unit circle
    unsigned lastEvenZero = (mNumZeros & 0xfffffffe) - 1;
    for (i = 0; i <= lastEvenZero; i+=2)
    {
      tempReal = testSolution[i];
      tempImag = testSolution[i+1];
      if ((tempReal*tempReal + tempImag*tempImag) > 1.0)
      {
        return 1.0e+300;
      }
    }

    if (kIntIsOdd(mNumZeros))
    {
      tempReal = testSolution[mNumZeros - 1];
      if ((tempReal * tempReal) > 1.0)
      {
        return 1.0e+300;
      }
    }
  }

  // Make sure there are no poles on or outside the unit circle
  // You probably will want to keep this too
  unsigned lastEvenPole = mNumZeros + (mNumPoles & 0xfffffffe) - 2;
  for (i = mNumZeros; i <= lastEvenPole; i+=2)
  {
    tempReal = testSolution[i];
    tempImag = testSolution[i+1];
    if ((tempReal*tempReal + tempImag*tempImag) > 0.999999999)
    {
      return 1.0e+300;
    }
  }

  // If you keep the for loop above, keep this too
  if (kIntIsOdd(mNumPoles))
  {
    tempReal = testSolution[mNumZeros + mNumPoles - 1];
    if ((tempReal * tempReal) > 1.0)
    {
      return 1.0e+300;
    }
  }

  double* evenZeros = &(testSolution[0]);
  double* evenPoles = &(testSolution[mNumZeros]);
  double* oddZero   = NULL;
  double* oddPole   = NULL;
  double  gain = testSolution[mNumZeros + mNumPoles];

  if (kIntIsOdd(mNumZeros))
  {
    oddZero = &(testSolution[mNumZeros - 1]);
  }

  if (kIntIsOdd(mNumPoles))
  {
    oddPole = &(testSolution[mNumZeros + mNumPoles - 1]);
  }

  ComputeSpectrum(evenZeros, mNumZeros & 0xfffffffe, oddZero,
                  evenPoles, mNumPoles & 0xfffffffe, oddPole,
                  gain, &mSpectrum);

  unsigned numPoints = mSpectrum.mNumValues;

/////////////////////////////////////////////////////////////////
//                                                             //
//   Use the impulse response, held in the variable            //
//   "mSpectrum", to compute a score for the solution that     //
//   has been passed into this function.  You probably don't   //
//   want to touch any of the code above this point, but       //
//   from here to the end of this function, it's all you!      //
//                                                             //
/////////////////////////////////////////////////////////////////

  #define kLnTwoToThe127 88.02969193111305
  #define kRecipLn10      0.4342944819032518

  // Compute square sum of errors for magnitude
  double magnitudeError = 0.0;
  double magnitude = 0.0;
  double logMagnitude = 0.0;
  tempReal = mSpectrum.mReals[0];
  tempImag = mSpectrum.mImags[0];
  magnitude = tempReal*tempReal + tempImag*tempImag;
  double baseMagnitude = 0.0;
  if (0.0 == magnitude)
  {
    baseMagnitude = -kLnTwoToThe127;
  }
  else
  {
    baseMagnitude = log(magnitude) * kRecipLn10;
    baseMagnitude *= 0.5;
  }

  for (i = 0; i < numPoints; i++)
  {
    tempReal = mSpectrum.mReals[i];
    tempImag = mSpectrum.mImags[i];
    magnitude = tempReal*tempReal + tempImag*tempImag;
    if (0.0 == magnitude)
    {
      logMagnitude = -kLnTwoToThe127;
    }
    else
    {
      logMagnitude = log(magnitude) * kRecipLn10;
      logMagnitude *= 0.5;  // Half the log because it's mag squared
    }

    logMagnitude -= baseMagnitude;
    magnitudeError += logMagnitude * logMagnitude;
  }

  // Compute errors for phase
  double phaseError = 0.0;
  double phase = 0.0;
  double componentError = 0.0;
  double degree = 1;//((mNumZeros + 1) & 0xfffffffe) - 1;
  double angleSpacing = -3.141592653589793 * 0.5 / numPoints * degree;
  double targetPhase = 0.0;
  double oldPhase = 0.0;
  double phaseDifference = 0;
  double totalPhaseTraversal = 0.0;
  double traversalError = 0.0;
  for (i = 0; i < (numPoints - 5); i++)
  {
    tempReal = mSpectrum.mReals[i];
    tempImag = mSpectrum.mImags[i];
    oldPhase = phase;
    phase = atan2(tempImag, tempReal);
    phaseDifference = phase - oldPhase;
    if (phaseDifference > 3.141592653589793)
    {
      phaseDifference -= 3.141592653589793;
      phaseDifference -= 3.141592653589793;
    }
    else if (phaseDifference < -3.141592653589793)
    {
      phaseDifference += 3.141592653589793;
      phaseDifference += 3.141592653589793;
    }
    totalPhaseTraversal += phaseDifference;
    componentError = cosh(200.0*(phaseDifference - angleSpacing)) - 0.5;
    phaseError += componentError * componentError;
    targetPhase += angleSpacing;
    if (targetPhase < -3.141592653589793)
    {
      targetPhase += 3.141592653589793;
      targetPhase += 3.141592653589793;
    }
  }

  traversalError = totalPhaseTraversal - angleSpacing * numPoints;
  traversalError *= traversalError;

  double baseMagnitudeError = baseMagnitude * baseMagnitude;

  // Compute weighted sum of the two subtotals
  // Take square root
  return sqrt(baseMagnitudeError*1.0 + magnitudeError*100.0 +
              phaseError*400.0 + traversalError*4000000.0);
}


///////////////////////////////
int main(int argc, char** argv)
{
  srand((unsigned)time(NULL));

  unsigned numZeros;
  unsigned numPoles;
  bool     minimumPhase;

  if (argc < 4)
  {
    printf("Usage: FilterDesign.exe <minimumPhase?> <numZeros> <numPoles>\n");
    return 0;
  }
  else
  {
    if (0 == atoi(argv[1]))
    {
      minimumPhase = false;
    }
    else
    {
      minimumPhase = true;
    }

    numZeros = (unsigned)atoi(argv[2]);
    if (0 == numZeros)
    {
      numZeros = 1;
    }

    numPoles = (unsigned)atoi(argv[3]);
  }

  unsigned vectorLength   = numZeros + numPoles + 1;
  unsigned populationSize = vectorLength * 10;
  FilterSolver theSolver(vectorLength, populationSize, 200,
                         numZeros, numPoles, minimumPhase);

  double* minimumSolution = new double[vectorLength];
  unsigned i;
  if (minimumPhase)
  {
    for (i = 0; i < numZeros; i++)
    {
      minimumSolution[i] = -sqrt(0.5);
    }
  }
  else
  {
    for (i = 0; i < numZeros; i++)
    {
      minimumSolution[i] = -10.0;
    }
  }

  for (; i < (vectorLength - 1); i++)
  {
    minimumSolution[i] = -sqrt(0.5);
  }

  minimumSolution[vectorLength - 1] = 0.0;

  double* maximumSolution = new double[vectorLength];
  if (minimumPhase)
  {
    for (i = 0; i < numZeros; i++)
    {
      maximumSolution[i] = sqrt(0.5);
    }
  }
  else
  {
    for (i = 0; i < numZeros; i++)
    {
      maximumSolution[i] = 10.0;
    }
  }

  for (i = 0; i < (vectorLength - 1); i++)
  {
    maximumSolution[i] = sqrt(0.5);
  }

  maximumSolution[vectorLength - 1] = 2.0;

  theSolver.Setup(minimumSolution, maximumSolution, 0, 0.5, 0.75);
  theSolver.Solve(1000000);

  double* bestSolution = theSolver.Solution();
  printf("\nZeros:\n");
  unsigned numEvenZeros = numZeros & 0xfffffffe;
  for (i = 0; i < numEvenZeros; i+=2)
  {
    printf("%.10f +/- %.10fi\n", bestSolution[i], bestSolution[i+1]);
  }

  if (kIntIsOdd(numZeros))
  {
    printf("%.10f\n", bestSolution[numZeros-1]);
  }

  printf("Poles:\n");
  unsigned lastEvenPole = numZeros + (numPoles & 0xfffffffe) - 2;
  for (i = numZeros; i <= lastEvenPole; i+=2)
  {
    printf("%.10f +/- %.10fi\n", bestSolution[i], bestSolution[i+1]);
  }

  unsigned numRoots = numZeros + numPoles;
  if (kIntIsOdd(numPoles))
  {
    printf("%.10f\n", bestSolution[numRoots-1]);
  }

  double gain = bestSolution[numRoots];
  printf("Gain: %.10f\n", gain);

  _getch();
  unsigned j;
  ASpectrum* spectrum = theSolver.Spectrum();
  double logMagnitude;
  printf("Magnitude Response, millibels:\n");
  for (i = 0; i < 20; i++)
  {
    for (j = 0; j < 10; j++)
    {
      logMagnitude = kRecipLn10 *
         log(spectrum->mReals[i*10 + j] * spectrum->mReals[i*10 + j] +
             spectrum->mImags[i*10 + j] * spectrum->mImags[i*10 + j]);
      if (logMagnitude < -9.999)
      {
        logMagnitude = -9.999;
      }
      printf("%+5.0f ", logMagnitude*1000);
    }
    printf("\n");
  }

  _getch();
  double phase;
  printf("Phase Response, milliradians:\n");
  for (i = 0; i < 20; i++)
  {
    for (j = 0; j < 10; j++)
    {
      phase = atan2(spectrum->mImags[i*10 + j], spectrum->mReals[i*10 + j]);
      printf("%+5.0f ", phase*1000);
    }
    printf("\n");
  }

  _getch();
  printf("Biquad Sections:\n");
  unsigned numBiquadSections =
    (numZeros > numPoles) ? ((numZeros + 1) >> 1) : ((numPoles + 1) >> 1);
  double x0, x1, x2;
  double y0, y1, y2;
  if (numZeros >=2)
  {
    x0 = (bestSolution[0]*bestSolution[0] + bestSolution[1]*bestSolution[1]) *
         gain;
    x1 = 2.0 * bestSolution[0] * gain;
    x2 = gain;
  }
  else if (1 == numZeros)
  {
    x0 = bestSolution[0] * gain;
    x1 = gain;
    x2 = 0.0;
  }
  else
  {
    x0 = gain;
    x1 = 0.0;
    x2 = 0.0;
  }

  if (numPoles >= 2)
  {
    y0 = (bestSolution[numZeros]*bestSolution[numZeros] +
          bestSolution[numZeros+1]*bestSolution[numZeros+1]);
    y1 = 2.0 * bestSolution[numZeros];
    y2 = 1.0;
  }
  else if (1 == numPoles)
  {
    y0 = bestSolution[numZeros];
    y1 = 1.0;
    y2 = 0.0;
  }
  else
  {
    y0 = 1.0;
    y1 = 0.0;
    y2 = 0.0;
  }

  x0 /= y0;
  x1 /= y0;
  x2 /= y0;
  y1 /= y0;
  y2 /= y0;

  printf("y[n] = %.10fx[n]", x0);
  if (numZeros > 0)
  {
    printf(" + %.10fx[n-1]", x1);
  }
  if (numZeros > 1)
  {
    printf(" + %.10fx[n-2]", x2);
  }
  printf("\n");

  if (numPoles > 0)
  {
    printf("                   + %.10fy[n-1]", y1);
  }
  if (numPoles > 1)
  {
    printf(" + &.10fy[n-2]", y2);
  }
  if (numPoles > 0)
  {
    printf("\n");
  }

  int numRemainingZeros = numZeros - 2;
  int numRemainingPoles = numPoles - 2;
  for (i = 1; i < numBiquadSections; i++)
  {
    if (numRemainingZeros >= 2)
    {
      x0 = (bestSolution[i*2]   * bestSolution[i*2] +
            bestSolution[i*2+1] * bestSolution[i*2+1]);
      x1 = -2.0 * bestSolution[i*2];
      x2 = 1.0;
    }
    else if (numRemainingZeros >= 1)
    {
      x0 = bestSolution[i*2];
      x1 = 1.0;
      x2 = 0.0;
    }
    else
    {
      x0 = 1.0;
      x1 = 0.0;
      x2 = 0.0;
    }

    if (numRemainingPoles >= 2)
    {
      y0 = (bestSolution[i*2+numZeros]   * bestSolution[i*2+numZeros] +
            bestSolution[i*2+numZeros+1] * bestSolution[i*2+numZeros+1]);
      y1 = -2.0 * bestSolution[i*2+numZeros];
      y2 = 1.0;
    }
    else if (numRemainingPoles >= 1)
    {
      y0 = bestSolution[i*2+numZeros];
      y1 = 1.0;
      y2 = 0.0;
    }
    else
    {
      y0 = 1.0;
      y1 = 0.0;
      y2 = 0.0;
    }

    x0 /= y0;
    x1 /= y0;
    x2 /= y0;
    y1 /= y0;
    y2 /= y0;

    printf("y[n] = %.10fx[n]", x0);
    if (numRemainingZeros > 0)
    {
      printf(" + %.10fx[n-1]", x1);
    }
    if (numRemainingZeros > 1)
    {
      printf(" + %.10fx[n-2]", x2);
    }
    printf("\n");

    if (numRemainingPoles > 0)
    {
      printf("                   + %.10fy[n-1]", -y1);
    }
    if (numRemainingPoles > 1)
    {
      printf(" + %.10fy[n-2]", -y2);
    }
    if (numRemainingPoles > 0)
    {
      printf("\n");
    }

    numRemainingZeros -= 2;
    numRemainingPoles -= 2;
  }

  _getch();
  printf("Full Expansion:\n");
  double* xpolynomial = new double[numRoots + 1];
  memset(xpolynomial, 0, sizeof(double) * (numRoots + 1));
  xpolynomial[0] = 1.0;
  if (numZeros >= 2)
  {
    xpolynomial[0] = bestSolution[0] * bestSolution[0] +
                     bestSolution[1] * bestSolution[1];
    xpolynomial[1] = -2.0 * bestSolution[0];
    xpolynomial[2] = 1.0;
  }
  else if (numZeros == 1)
  {
    xpolynomial[0] = bestSolution[0];
    xpolynomial[1] = 1.0;
  }
  else
  {
    xpolynomial[0] = 1.0;
  }

  for (i  = 2, numRemainingZeros = numZeros; numRemainingZeros >= 2;
       i += 2, numRemainingZeros-=2)
  {
    x2 = 1.0;
    x1 = -2.0 * bestSolution[i];
    x0 = bestSolution[i]   * bestSolution[i] +
         bestSolution[i+1] * bestSolution[i+1];
    for (j = numRoots; j > 1; j--)
    {
      xpolynomial[j] = xpolynomial[j-2] + xpolynomial[j-1] * x1 +
                       xpolynomial[j] * x0;
    }
    xpolynomial[1]  = xpolynomial[0] * x1 + xpolynomial[1] * x0;
    xpolynomial[0] *= x0;
  }

  if (numRemainingZeros > 0)
  {
    x1 = 1.0;
    x0 = bestSolution[numZeros-1];
    for (j = numRoots; j > 0; j--)
    {
      xpolynomial[j] = xpolynomial[j-1] + xpolynomial[j] * x0;
    }
    xpolynomial[0] *= x0;
  }

  double* ypolynomial = new double[numRoots + 1];
  memset(ypolynomial, 0, sizeof(double) * (numRoots + 1));
  ypolynomial[0] = 1.0;
  if (numPoles >= 2)
  {
    ypolynomial[0] = bestSolution[numZeros]   * bestSolution[numZeros] +
                     bestSolution[numZeros+1] * bestSolution[numZeros+1];
    ypolynomial[1] = -2.0 * bestSolution[numZeros];
    ypolynomial[2] = 1.0;
  }
  else if (numPoles == 1)
  {
    ypolynomial[0] = bestSolution[numZeros];
    ypolynomial[1] = 1.0;
  }
  else
  {
    xpolynomial[0] = 1.0;
  }

  for (i  = 2, numRemainingPoles = numPoles; numRemainingPoles >= 2;
       i += 2, numRemainingPoles-=2)
  {
    y2 = 1.0;
    y1 = -2.0 * bestSolution[numZeros+i];
    y0 = bestSolution[numZeros+i]   * bestSolution[numZeros+i] +
         bestSolution[numZeros+i+1] * bestSolution[numZeros+i+1];
    for (j = numRoots; j > 1; j--)
    {
      ypolynomial[j] = ypolynomial[j-2] + ypolynomial[j-1] * y1 +
                       ypolynomial[j] * y0;
    }
    ypolynomial[1]  = ypolynomial[0] * y1 + ypolynomial[1] * y0;
    ypolynomial[0] *= y0;
  }

  if (numRemainingPoles > 0)
  {
    y1 = 1.0;
    y0 = bestSolution[numZeros+numPoles-1];
    for (j = numRoots; j > 0; j--)
    {
      ypolynomial[j] = ypolynomial[j-1] + ypolynomial[j] * y0;
    }
    ypolynomial[0] *= y0;
  }

  y0 = ypolynomial[0];
  for (i = 0; i <= numRoots; i++)
  {
    xpolynomial[i] /= y0;
    ypolynomial[i] /= y0;
  }

  printf("y[n] = %.10fx[n]", xpolynomial[0]*gain);
  for (i = 1; i <= numZeros; i++)
  {
    printf(" + %.10fx[n-%d]", xpolynomial[i]*gain, i);
    if ((i % 3) == 2)
    {
      printf("\n");
    }
  }

  if ((i % 3) != 0)
  {
    printf("\n");
  }

  if (numPoles > 0)
  {
    printf("                 ");
  }

  for (i = 1; i <= numPoles; i++)
  {
    printf(" + %.10fy[n-%d]", -ypolynomial[i], i);
    if ((i % 3) == 2)
    {
      printf("\n");
    }
  }

  if ((i % 3) != 0)
  {
    printf("\n");
  }

  delete[] minimumSolution;
  delete[] maximumSolution;
  delete[] xpolynomial;
  delete[] ypolynomial;
}


bool ComputeSpectrum(double* evenZeros, unsigned numEvenZeros, double* oddZero,
                     double* evenPoles, unsigned numEvenPoles, double* oddPole,
                     double gain, ASpectrum* spectrum)
{
  unsigned i, j;

  // For equally spaced points on the unit circle
  unsigned numPoints = spectrum->mNumValues;
  double   spacingAngle = 3.141592653589793 / (numPoints - 1);
  double   pointArgument = 0.0;
  double   pointReal = 0.0;
  double   pointImag = 0.0;
  double   rootReal = 0.0;
  double   rootImag = 0.0;
  double   differenceReal = 0.0;
  double   differenceImag = 0.0;
  double   responseReal = 1.0;
  double   responseImag = 0.0;
  double   recipSquareMagnitude = 0.0;
  double   recipReal = 0.0;
  double   recipImag = 0.0;
  double   tempRealReal = 0.0;
  double   tempRealImag = 0.0;
  double   tempImagReal = 0.0;
  double   tempImagImag = 0.0;

  for (i = 0; i < numPoints; i++)
  {
    responseReal = 1.0;
    responseImag = 0.0;

    // The imaginary component is negated because we're using 1/z, not z
    pointReal =  cos(pointArgument);
    pointImag = -sin(pointArgument);

    // For each even zero
    for (j = 0; j < numEvenZeros; j+=2)
    {
      rootReal = evenZeros[j];
      rootImag = evenZeros[j + 1];
      // Compute distance from that zero to that point
      differenceReal = pointReal - rootReal;
      differenceImag = pointImag - rootImag;
      // Multiply that distance by the accumulating product
      tempRealReal = responseReal * differenceReal;
      tempRealImag = responseReal * differenceImag;
      tempImagReal = responseImag * differenceReal;
      tempImagImag = responseImag * differenceImag;
      responseReal = tempRealReal - tempImagImag;
      responseImag = tempRealImag + tempImagReal;
      // Do the same with the conjugate root
      differenceImag = pointImag + rootImag;
      tempRealReal = responseReal * differenceReal;
      tempRealImag = responseReal * differenceImag;
      tempImagReal = responseImag * differenceReal;
      tempImagImag = responseImag * differenceImag;
      responseReal = tempRealReal - tempImagImag;
      responseImag = tempRealImag + tempImagReal;
      // The following way is little faster, if any
      // response *= (1/z - r) * (1/z - conj(r))
      //          *= r*conj(r) - (r + conj(r))/z + 1/(z*z)
      //          *= real(r)*real(r) + imag(r)*imag(r) - 2*real(r)/z + 1/(z*z)
      //          *= ... - 2*real(r)*conj(z) + conj(z)*conj(z)
      //          *= ... - 2*real(r)*real(z) + 2i*real(r)*imag(z) +
      //             real(z)*real(z) - 2i*real(z)*imag(z) + imag(z)*imag(z)
      //          *= real(r)*real(r) + imag(r)*imag(r) - 2*real(r)*real(z) +
      //             real(z)*real(z) + imag(z)*imag(z) +
      //              2i * imag(z) * (real(r) - real(z))
      //          *= (real(r) - real(z))^2  + imag(r)^2 + imag(z)^2 +
      //              2i * imag(z) * (real(r) - real(z))
      // This ends up being 8 multiplications, 6 additions
    }

    if (NULL != oddZero)
    {
      rootReal = *oddZero;
      // Compute distance from that zero to that point
      differenceReal = pointReal - rootReal;
      differenceImag = pointImag;
      // Multiply that distance by the accumulating product
      tempRealReal = responseReal * differenceReal;
      tempRealImag = responseReal * differenceImag;
      tempImagReal = responseImag * differenceReal;
      tempImagImag = responseImag * differenceImag;
      responseReal = tempRealReal - tempImagImag;
      responseImag = tempRealImag + tempImagReal;
    }

    // For each pole
    for (j = 0; j < numEvenPoles; j+=2)
    {
      rootReal = evenPoles[j];
      rootImag = evenPoles[j + 1];
      // Compute distance from that pole to that point
      differenceReal = pointReal - rootReal;
      differenceImag = pointImag - rootImag;
      // Multiply the reciprocal of that distance by the accumulating product
      recipSquareMagnitude = 1.0 / (differenceReal * differenceReal +
                                    differenceImag * differenceImag);
      recipReal =  differenceReal * recipSquareMagnitude;
      recipImag = -differenceImag * recipSquareMagnitude;
      tempRealReal = responseReal * recipReal;
      tempRealImag = responseReal * recipImag;
      tempImagReal = responseImag * recipReal;
      tempImagImag = responseImag * recipImag;
      responseReal = tempRealReal - tempImagImag;
      responseImag = tempRealImag + tempImagReal;
      // Do the same with the conjugate root
      differenceImag = pointImag + rootImag;
      recipSquareMagnitude = 1.0 / (differenceReal * differenceReal +
                                    differenceImag * differenceImag);
      recipReal =  differenceReal * recipSquareMagnitude;
      recipImag = -differenceImag * recipSquareMagnitude;
      tempRealReal = responseReal * recipReal;
      tempRealImag = responseReal * recipImag;
      tempImagReal = responseImag * recipReal;
      tempImagImag = responseImag * recipImag;
      responseReal = tempRealReal - tempImagImag;
      responseImag = tempRealImag + tempImagReal;
    }

    if (NULL != oddPole)
    {
      rootReal = *oddPole;
      // Compute distance from that point to that zero
      differenceReal = pointReal - rootReal;
      differenceImag = pointImag;
      // Multiply the reciprocal of that distance by the accumulating product
      recipSquareMagnitude = 1.0 / (differenceReal * differenceReal +
                                    differenceImag * differenceImag);
      recipReal =  differenceReal * recipSquareMagnitude;
      recipImag = -differenceImag * recipSquareMagnitude;
      tempRealReal = responseReal * recipReal;
      tempRealImag = responseReal * recipImag;
      tempImagReal = responseImag * recipReal;
      tempImagImag = responseImag * recipImag;
      responseReal = tempRealReal - tempImagImag;
      responseImag = tempRealImag + tempImagReal;
    }

    // Multiply by the gain
    responseReal *= gain;
    responseImag *= gain;

    spectrum->mReals[i] = responseReal;
    spectrum->mImags[i] = responseImag;

    pointArgument += spacingAngle;
  }

  return true;
}

// Half-band lowpass
/*
  #define kLnTwoToThe127 88.02969193111305
  #define kRecipLn10      0.4342944819032518

  // Compute square sum of errors for bottom half band
  unsigned numLoBandPoints = numPoints >> 1;
  double loBandError = 0.0;
  double magnitude = 0.0;
  double logMagnitude = 0.0;
  for (i = 0; i < numLoBandPoints; i++)
  {
    tempReal = mSpectrum.mReals[i];
    tempImag = mSpectrum.mImags[i];
    magnitude = tempReal*tempReal + tempImag*tempImag;
    if (0.0 == magnitude)
    {
      logMagnitude = -kLnTwoToThe127;
    }
    else
    {
      logMagnitude = log(magnitude) * kRecipLn10;
      logMagnitude *= 0.5;  // Half the log because it's mag squared
    }

    loBandError += logMagnitude * logMagnitude;
  }

  // Compute errors for top half of band
  double hiBandError = 0.0;
  for ( ; i < numPoints; i++)
  {
    tempReal = mSpectrum.mReals[i];
    tempImag = mSpectrum.mImags[i];
    magnitude = tempReal*tempReal + tempImag*tempImag;
    hiBandError += magnitude; // Already a squared value
  }

  // Compute weighted sum of the two subtotals
  // Take square root
  return sqrt(loBandError + 5000.0 * hiBandError);
*/