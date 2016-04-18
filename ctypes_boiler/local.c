#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "local.h"

#define PI 3.14159265359

/*
simplest example, args and ret by value
*/

double sin_degrees (double x)
{
   double ret;
   ret = sin(x * PI / 180.0);
   return(ret);
}


/*
args and ret by reference
*/

double *sin_degrees_p (double *x)
{
   // don't return local pointer (use static)
   static double *ret;
   ret = malloc(sizeof(double) * 1);
   ret[0] = sin((*x) * PI / 180.0);
   // printf("The sine of %lf is %lf, %lf.", *x, *ret, ret[0]);
   return(ret);
}


/*
numpy array args
*/

double sum(double * array, int n, int m){
    double total = 0.0;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            total = total + array[i + m * j];
        }
    }

    printf("(%d,%d) = %lf\n", 1, 2, array[2 + m * 1]);
    
    return(total);
}