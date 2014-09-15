// double correlation( double* x1, int l1, double* x2, int l2);

#include <map>
using namespace std;
class C_matrix {
  private:
  	double sigma2;
  	double L2;
    map<unsigned long long int,double> C_sparse;
    int nrow;
    int ncol;
    double Lx;
    double Ly;
    double N;
  public:
  	// C_matrix(int n) : N(n){ };
  	C_matrix();
  	C_matrix(double sigma, double L);
  	double get_elem(int ind_i, int ind_j);
  	double correlation( double* x1, int l1, double* x2, int l2);
    void make_C(int nrow, int ncol, double Lx, double Ly);
    void av(int n, double *in, int nin, double *out, int nout);
    void set_dims(int nrow, int ncol, double Lx, double Ly);
    void av_no_C(int n, double *in, int nin, double *out, int nout);
    // void av(int n, double *in, double *out);
    // void av_no_C(int n, double *in, double *out);

};