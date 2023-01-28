#include <cstdlib>
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h> 
#include "FilterMath.h"

/*====================================================
  acosh() - calcs the inverse hyperbolic cosine of x
  Prototype:  double acosh(double x);
  Return:     inverse hyperbolic cosine value
  Arguments:  x - argument of acosh
====================================================*/
double acosh (double x)
{ 
   return log (x + sqrt (x * x - 1.0));
}

/*====================================================
  arcsc() - calcs Jacobian elliptic arctan function
              using arithmetic-geometric mean method
  Prototype:  double arcsc(double u,double k);
  Return:     Jacobian elliptic arctan function
  Arguments:  u - the value of the elliptic integral
              k - the modulus of the integral
====================================================*/
double arcsc (double u, double k)
{ 
   int     i, L;
   double  A, B, BT, Y;

  /*  Initialize the starting values. */
  A = 1.0;
  B = k;
  Y = 1.0 / u;
  L = 0;

  //  Iterate until error is small enough
  for (i = 0; i < MAX_TERMS; i++)
  { 
     BT = A * B;
      A = A + B;
      B = 2 * sqrt(BT);
      Y = Y - BT / Y;
      if (Y == 0.0)
      { 
         Y = sqrt (BT) * ERR_SMALL;
      }
      if (fabs (A - B) < (A * ERR_SMALL))
      { 
         break;
      }
      L = 2 * L;
      if(Y < 0)
      { 
         L++;
      }
  }

  if(Y < 0)
  { 
     L++;
  }
  return ((atan (A / Y) + PI * L) / A);
}

/*====================================================
  asinh() - calcs the inverse hyperbolic sine of x
  Prototype:  double asinh(double x);
  Return:     inverse hyperbolic sine value
  Arguments:  x - argument of asinh
====================================================*/
double asinh (double x)
{ 
   return log(x + sqrt (x * x + 1.0));
}

/*====================================================
  Ellip_Funcs() - calcs Jacobian elliptic functions
              using arithmetic-geometric mean method
  Prototype:  void Ellip_Funcs (double u, double k,
                    double *sn, double *cn, double *dn);
  Return:     none (values return via arg list)
  Arguments:  u - the value of the elliptic integral
              k - the modulus of the integral
              sn - Jacobian elliptic sine function
              cn - Jacobian elliptic cosine function
              dn - Jacobian elliptic differenc func
====================================================*/
void Ellip_Funcs (double u, double k, double *sn, double *cn, double *dn)
{ 
   int     i, imax;   /*  Loop counter and max value. */
   double  A[MAX_TERMS],B[MAX_TERMS],  /* Array storage*/
           C[MAX_TERMS],P[MAX_TERMS];  /* values. */

  // Square the modulus as required in this method
  k = k * k;

  // Initialize the starting values
  A[0] = 1.0;
  B[0] = sqrt (1-k);
  C[0] = sqrt (k);

  //  Iterate until error is small enough
  for (i = 1; i < MAX_TERMS; i++)
  { 
     A[i] = (A[i - 1] + B[i - 1]) / 2.0;
     B[i] = sqrt (A[i - 1] * B[i - 1]);
     C[i] = (A[i - 1] - B[i - 1]) / 2.0;
     if (C[i] < ERR_SMALL)
     { 
        break;
     }
  }

  //  Get last value of i
  if (i == MAX_TERMS)
  { 
     imax = i - 1;
  }
  else
  { 
     imax = i;
  }
  
  // Determine phase and unwrap it
  P[imax] = pow (2.0, imax) * A[imax] * u;
  for (i = imax; i > 0 ; i--)
  { 
     P[i - 1] = (asin(C[i] * sin (P[i]) / A[i]) + P[i]) / 2.0;
  }
  
  //  Place values in memory locations specified. */
  *sn = sin(P[0]);
  *cn = cos(P[0]);
  *dn = sqrt(1 - k * (*sn) * (*sn));
}

/*====================================================
  Ellip_Integral() - calcs complete elliptic integral
              using arithmetic-geometric mean method
  Prototype:  void Ellip_Integral(double k);
  Return:     complete elliptic integral value
  Arguments:  k - the modulus of the integral
====================================================*/
double Ellip_Integral (double k)
{ 
   int     i;      // Loop counter
   double  A[MAX_TERMS], B[MAX_TERMS], C[MAX_TERMS]; // Array storage values

  //  Square the modulus as required by this method
  k = k * k;
  //  Initialize the starting values
  A[0] = 1.0;
  B[0] = sqrt (1 - k);
  C[0] = sqrt (k);

  //  Iterate until error is small enough
  for (i = 1; i < MAX_TERMS ; i++)
  { 
     A[i] = (A[i - 1] + B[i - 1])/2.0;
     B[i] = sqrt (A[i - 1] * B[i - 1]);
     C[i] = (A[i - 1] - B[i - 1]) / 2.0;
     if (C[i] < ERR_SMALL)
     { 
        break;
     }
  }
  return PI / (2.0 * A[i]);
}

/*====================================================
  Io () - calcs the mod Bessel function of first kind
  Prototype:  double Io (double value);
  Return:     modified Bessel function
  Arguments:  value - the argument of Bessel function
====================================================*/
double Io (double value)
{ 
   int     i, converge;    // Loop counter and indicator
   double  Iold, Inew, J, K;  // Iteration constants

  //  Initialize values, set convergence to false
  Iold = 1.0;
  J = 1.0;
  K = value / 2.0;
  converge = 0;

  //  Use series expansion definition of Bessel
  for (i = 1; i < MAX_TERMS; i++)
  { 
     J *= K / i;
     Inew = Iold + J * J;
     if((Inew - Iold) < ERR_SMALL)
     { 
        converge = 1; 
        break;
     }
     Iold = Inew;
  }
  if (!converge) 
  { 
     return 0.0;
  }
  return Inew;
}

/*====================================================
  Calc_Ellipt_Coefs() - calcs normal elliptic coefs
  Prototype:  int Calc_Ellipt_Coefs(Filt_Params *FP);
  Return:     error value
  Arguments:  FP - ptr to struct holding filter params
====================================================*/
int Calc_Ellipt_Coefs (Filt_Params *FP)
{ int     m,a,b,            /*  loop counters and indices */
          odd;              /*  odd order indicator */
  double  kernel,ratio,     /*  internal constants  */
          epsilon,kn,rt,    /*  internal constants  */
          sigma,omega,zero, /*  pole and zero       */
          CEIrt,CEIkn,vo,fm,/*  elliptic constants  */
          SP,CP,DP,SN,CN,DN;/*  ellip sin,cos,dif   */

  /*  Check for NULL ptrs and zero order. */
  if(!FP->acoefs)
  { return ERR_NULL;}
  if(!FP->bcoefs)
  { return ERR_NULL;}
  if(FP->order == 0)
  { return ERR_VALUE;}

  /*  Make calculations of necessary constants. */
  epsilon = sqrt (pow (10.0, -0.1 * FP->apass1) - 1.0 );

  /*  Calc freq ratio based on filter selectivity. */
  switch (FP->select)
  { case 'L': ratio = FP->wstop1 / FP->wpass1;
              break;
    case 'H': ratio = FP->wpass1 / FP->wstop1;
              break;
    case 'P': ratio = (FP->wstop2 - FP->wstop1) / (FP->wpass2 - FP->wpass1);
              break;
    case 'S': ratio = (FP->wpass2-FP->wpass1) / (FP->wstop2-FP->wstop1);
              break;
    default:  return ERR_FILTER;
  }

  //  Calc values used in future calculations
  kernel = (pow (10, -0.1 * FP->astop1) - 1.0 ) / (pow (10.0, -0.1 * FP->apass1) -1.0);
  rt = 1.0 / ratio;
  kn = 1.0 / sqrt (kernel);

  //  Calc Ellip_Integrals, vo and Ellip_Funcs
  CEIrt = Ellip_Integral (rt);
  CEIkn = Ellip_Integral (kn);
  vo = (CEIrt / (CEIkn * FP->order)) * arcsc (1.0/ epsilon, kn);
  Ellip_Funcs (vo, sqrt (1.0 - rt * rt), &SP, &CP, &DP);
  
  //  Start indices at 0
  a = 0; 
  b = 0;

  //  Handle odd order if necessary
  odd = FP->order % 2;
  
  if (odd)
  { 
     FP->gain = 1.0;   /* Gain set to 1.0 for odd. */
     FP->acoefs[a++] = 0.0;
     FP->acoefs[a++] = 0.0;
     FP->acoefs[a++] = SP * CP / (1.0 - SP *SP);
     FP->bcoefs[b++] = 0.0;
     FP->bcoefs[b++] = 1.0;
     FP->bcoefs[b++] = SP * CP / (1 - SP * SP);
  } else {                 /*  Gain adjusted for even. */
     FP->gain = pow (10.0, 0.05 * FP->apass1);
  }

  //  Handle all quadratic terms
  for (m = 0; m < FP->order / 2; m++)
  { 
     // Make intermediate calculations
     fm = CEIrt * (2 * m + 1 + odd) / FP->order;
     Ellip_Funcs (fm, rt, &SN, &CN, &DN);
     // Calc real and imag coordinates of poles
     sigma = -1.0 * CN * DN * SP * CP / (1.0 - DN * DN * SP * SP);
     omega = SN * DP / (1.0 - DN * DN * SP * SP);

     //  Calculate the zero location
     zero = 1.0 / (rt * SN);

     // Set the quadratic coefs. 
     FP->acoefs[a++] = 1.0;
     FP->acoefs[a++] = 0.0;
     FP->acoefs[a++] = zero * zero;
     FP->bcoefs[b++] = 1.0;
     FP->bcoefs[b++] = -2.0 * sigma;
     FP->bcoefs[b++] = sigma * sigma + omega * omega;
    
     //  Update the gain
     FP->gain *= ((sigma * sigma + omega * omega) / (zero * zero));
  }
  return ERR_NONE;
}