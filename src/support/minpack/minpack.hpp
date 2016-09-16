#include "defines.h"

void chkder ( int m, int n, real_t x[], real_t fvec[], real_t fjac[],
  int ldfjac, real_t xp[], real_t fvecp[], int mode, real_t err[] );
void dogleg ( int n, real_t r[], int lr, real_t diag[], real_t qtb[],
  real_t delta, real_t x[], real_t wa1[], real_t wa2[] );
real_t enorm ( int n, real_t x[] );
void fdjac1 ( void fcn (void* usr_data, int n, real_t x[], real_t f[], int *iflag ),
  void* usr_data, int n, real_t x[], real_t fvec[], real_t fjac[], int ldfjac, int *iflag,
  int ml, int mu, real_t epsfcn, real_t wa1[], real_t wa2[] );
void fdjac2 ( void fcn ( int m, int n, real_t x[], real_t fvec[], int *iflag ),
  int m, int n, real_t x[], real_t fvec[], real_t fjac[], int ldfjac,
  int *iflag, real_t epsfcn, real_t wa[] );
int hybrd ( void fcn (void* usr_data, int n, real_t x[], real_t fvec[], int *iflag ),
  void* usr_data, int n, real_t x[], real_t fvec[], real_t xtol, int maxfev, int ml,
  int mu, real_t epsfcn, real_t diag[], int mode, real_t factor, int nprint,
  int nfev, real_t fjac[], int ldfjac, real_t r[], int lr, real_t qtf[],
  real_t wa1[], real_t wa2[], real_t wa3[], real_t wa4[] );
int hybrd1 ( void fcn (void* usr_data, int n, real_t x[], real_t fvec[], int *iflag ),
             void* usr_data, int n, real_t x[], real_t fvec[], real_t tol, real_t wa[], int lwa );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
void qform ( int m, int n, real_t q[], int ldq, real_t wa[] );
void qrfac ( int m, int n, real_t a[], int lda, bool pivot, int ipvt[],
  int lipvt, real_t rdiag[], real_t acnorm[], real_t wa[] );
void r1mpyq ( int m, int n, real_t a[], int lda, real_t v[], real_t w[] );
bool r1updt ( int m, int n, real_t s[], int ls, real_t u[], real_t v[], real_t w[] );
real_t r8_abs ( real_t x );
real_t r8_epsilon ( );
real_t r8_huge ( );
real_t r8_max ( real_t x, real_t y );
real_t r8_min ( real_t x, real_t y );
real_t r8_tiny ( );
real_t r8_uniform_01 ( int *seed );
void r8mat_print ( int m, int n, real_t a[], char* title );
void r8mat_print_some ( int m, int n, real_t a[], int ilo, int jlo, int ihi,
  int jhi, char* title );
void r8vec_print ( int n, real_t a[], char* title );
void timestamp ( );
