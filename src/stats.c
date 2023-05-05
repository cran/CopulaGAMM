/******************************************************************/
/*                                                                */
/*   Library of functions for goodness-of-fit tests               */
/*                                                                */
/*   Bruno Remillard, Nov 29, 2021                                */
/******************************************************************/




#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>



   double maxi(double a, double b)
   {
      if (a > b)
         return a;
      else
         return b;
   }


void stats(double *y, int *n, double *Fn, double *Fm)
{
int i,j;
double sum1,sum2,yy;


for(i=0;i<n[0];i++)
{
  sum1 = 0.0;
  sum2 = 0.0;
  yy = y[i];


  for(j=0;j<n[0];j++)
  {
      sum1 += (y[j] <= yy);
      sum2 += (y[j] <  yy);
  }
Fn[i]  = sum1/((double) n[0]);
Fm[i] = sum2/((double) n[0]);
}

}




