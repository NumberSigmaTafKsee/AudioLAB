float chebyshev(float x, float A[], int order)
{
   // To = 1
   // T1 = x
   // Tn = 2.x.Tn-1 - Tn-2
   // out = sum(Ai*Ti(x)) , i C {1,..,order}
   float Tn_2 = 1.0f;
   float Tn_1 = x;
   float Tn;
   float out = A[0]*Tn_1;

   for(int n=2;n<=order;n++)
   {
      Tn    =   2.0f*x*Tn_1 - Tn_2;
      out    +=   A[n-1]*Tn;
      Tn_2 =   Tn_1;
      Tn_1 =  Tn;
   }
   return out;
}
