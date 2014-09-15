// we have a 300x300 grid,
// 90 000x90 000 correlation matrix
// with 8 100 000 000 elements
// 4 bytes each would make 32 GB
// Harry only has 16GB RAM, so don't save it
// we have two options:
//    1. make C sparse ( make_C then use av )
//    2. compute C's elements every time, don't save it ( set_dims then use av_no_C)

#include "correlation.hpp"

#include <iostream>
#include <map>
#include <math.h>       /* exp */
#include <omp.h>

using namespace std;


C_matrix::C_matrix(){
  this->sigma2 = pow(0.2, 2);
  this->L2 = pow(3000,2);
}

C_matrix::C_matrix(double sigma, double L){
  this->sigma2 = pow(sigma, 2);
  this->L2 = pow(L,2);
}

double C_matrix::correlation ( double* x1, int l1, double* x2, int l2){
  double sigma2 = this->sigma2;
  double L2 = this->L2;
  double D = 2;  // for (int i = 0; i < N; ++i)
  double sum = pow(x1[0] - x2[0], 2) + pow(x1[1] - x2[1], 2);
  return sigma2 * exp(-D*sum/L2);
}

double C_matrix::get_elem(int ind_i, int ind_j){
  unsigned long long int i = (unsigned long long int) ind_i;
  unsigned long long int j = (unsigned long long int) ind_j;
  unsigned long long int unique_key;
  double out;

  if (j >= i){
    unique_key = (i<<32) | j;
  }
  else{
    unique_key = (j<<32) | i;
  }

  if(this->C_sparse.count(unique_key) > 0){
    out = this->C_sparse[unique_key];
  }
  else out = 0.0;

  return out;
}

void C_matrix::make_C(int nrow, int ncol, double Lx, double Ly){

  // we want a unique map key from each i,j pair
  // using i, j and unique_key as 64 bit ints (unsigned long long int)
  // the lower 32 bits are i, higher 32 bits are j
  // masking shouldn't be necessary, low8bits = i & 0xFF ??
  if (pow(2,31) < nrow*ncol )
  {
    cout << "ERROR: key masking not possible" << endl;
  }

  unsigned long long int i,j;
  unsigned long long int unique_key;

  unsigned long long int Nrow = (unsigned long long int) nrow;
  unsigned long long int Ncol = (unsigned long long int) ncol;

  double x1[2];
  double x2[2];
  double val;

  for (unsigned long long int row1 = 0; row1 < Nrow; row1++)
  {
    for (unsigned long long int col1 = 0; col1 < Ncol; col1++)
    {
      for (unsigned long long int row2 = 0; row2 < Nrow; row2++)
      {
        for (unsigned long long int col2 = 0; col2 < Ncol; col2++)
        {

          i = row1*Ncol+col1;
          j = row2*Ncol+col2;

          // if on or above the diagonal
          if (j >= i){

            x1[0] = (row1+0.5)*Lx;
            x1[1] = (col1+0.5)*Ly;
            x2[0] = (row2+0.5)*Lx;
            x2[1] = (col2+0.5)*Ly;
            val = correlation(x1, 2, x2, 2);

            // save a lot of memory by ignoring small values
            if (val > 1.0e-10){
              unique_key = (i<<32) | j;
              this->C_sparse[unique_key] = val;
            }
            // NOTE: if j<i look for c_sparse[j,i] not c_sparse[i,j], C_sparse.count(key)>0 if contains

          }

        }
      }
    }
  }

}


void C_matrix::av(int n, double *in, int nin, double *out, int nout)
{


  unsigned long long int i,j;
  unsigned long long int unique_key;

  unsigned long long int N = (unsigned long long int) n;

  // C_matrix has N=nrow*col rows and columns
  for (i = 0; i < N; i++)
  {
    out[i] = 0;
  }


  for (i = 0; i < N; i++)
  {
      for (j = 0; j < N; j++)
      {
        // C_sparse is symmetric, but we only stored values on or above diagonal
        if (j >= i){
          unique_key = (i<<32) | j;
        }
        else{
          unique_key = (j<<32) | i;
        }

        if(this->C_sparse.count(unique_key) > 0){
          out[i] += in[j]*(this->C_sparse[unique_key]); // out[i] += in[j]*C_matrix[i][j];
        }

      }
  }

  // cout << "calling av" << endl;

}

void C_matrix::set_dims(int nrow, int ncol, double Lx, double Ly){
  this->nrow = nrow;
  this->ncol = ncol;
  this->Lx = Lx;
  this->Ly = Ly;
}


 // sets out to A * in (matrix multiplication)
void C_matrix::av_no_C(int n, double *in, int nin, double *out, int nout)
{

  int nrow = this->nrow;
  int ncol = this->ncol;
  double Lx = this->Lx;
  double Ly = this->Ly;

  #pragma omp parallel
  {
    // private(x1,x2,val,i,j, k row1, row2,col1, col2)

    double x1[2];
    double x2[2];
    double val;
    int i,j;
    int N = n;
    // C_matrix has N=nrow*col rows and columns
    #pragma omp for
    for (int k = 0; k < N; ++k)
    {
      out[k] = 0;
    }

    #pragma omp barrier

    #pragma omp for
    for (int row1 = 0; row1 < nrow; row1++)
    {
      for (int col1 = 0; col1 < ncol; col1++)
      {
        // C_matrix is symmetric
        // out[row1*ncol+col1] = 0;
        for (int row2 = 0; row2 < nrow; row2++)
        {
          for (int col2 = 0; col2 < ncol; col2++)
          {
            i = row1*ncol+col1;
            j = row2*ncol+col2;

            // the if statement won't work with default openMP

            // if on or above the diagonal
            // if (j >= i){

              x1[0] = (row1+0.5)*Lx;
              x1[1] = (col1+0.5)*Ly;
              x2[0] = (row2+0.5)*Lx;
              x2[1] = (col2+0.5)*Ly;
              val = correlation(x1, 2, x2, 2);

              out[i] += in[j]*val;

              // if above diagonal set reflection
              // if (j > i){
                // out[j] += in[i]*val;
              // }

            // }


          }
        }

      }
    }

  }
  // cout << "calling av" << endl;

}
