#include <stdio.h>
#include <math.h>
#include "local.h"

#define PI 3.14159265

double sin_degrees (double x)
{
   double ret;
   ret = sin(x * PI / 180.0);
   printf("The sine of %lf is %lf degrees.", x, ret);
   return(ret);
}