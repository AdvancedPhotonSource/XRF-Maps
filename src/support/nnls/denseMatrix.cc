// File: denseMatrix.cc -- implements dense matrix functionality 

// Author: Suvrit Sra <suvrit@tuebingen.mpg.de>
// (c) Copyright 2010   Suvrit Sra
// Max-Planck-Institute for Biological Cybernetics

// denseMatrix.cc - implements dense matrix functionality
// Copyright (C) 2010 Suvrit Sra (suvrit@tuebingen.mpg.de)

// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

#include "denseMatrix.h"
#include <cstdio>

#ifdef _HAVE_CBLAS
#include <cblas.h>
#endif 

using namespace nsNNLS;

int denseMatrix::load(const char* fn, bool asbin)
{
  if (asbin)
    return load_as_bin(fn);
  else
    return load_as_txt(fn);
}

int denseMatrix::load_as_bin(const char* fn)
{
  FILE* fp = fopen (fn, "rb");
  if (fp == 0) {
    fprintf(stderr, "denseMatrix:: Error opening `%s' for reading\n", fn);
    return -1;
  }

  size_t m;
  size_t n;

  // read the dimension info
  if (fread(&m,  sizeof(size_t), 1, fp) != 1) {
    fprintf(stderr, "denseMatrix:: Error reading M (numrows) from %s\n", fn);
    return -1;
  }

  if (fread(&n,  sizeof(size_t), 1, fp) != 1) {
    fprintf(stderr, "denseMatrix:: Error reading N (numcols) from %s\n", fn);
    return -2;
  }

  fprintf(stderr, "Found a %zu x %zu matrix in %s\n", m, n, fn);

  setsize(m, n);
  size = m*n;
  data = new double[size];
  if (!data) {
    fprintf(stderr, "Could not allocate memory for matrix\n");
    return -4;
  }
    
  if (fread(data, sizeof(double), size, fp) != size) {
    fprintf(stderr, "denseMatrix:: Error reading data from %s\n", fn);
    return -3;
  }

  fclose(fp);
  return 0;
}

int denseMatrix::load_as_txt(const char* fn)
{
  FILE* fp = fopen(fn, "r");
  size_t m, n;
  int r;

  r = fscanf(fp, "%zu %zu\n", &m, &n);
  setsize(m, n);
  size = m*n;
  data = new double[size];
  for (size_t i = 0; i < size; i++) 
    r = fscanf(fp, "%lf", &data[i]);

  fclose(fp);

  return 0;
}

/// Returns 'r'-th row into pre-alloced vector
int denseMatrix::get_row (size_t i, vector*& v)
{
  return -1;
}

/// Returns 'c'-th col as a vector
int denseMatrix::get_col (size_t j, vector*& c)
{
  return -1;
}

/// Returns main or second diagonal (if p == true)
int denseMatrix::get_diag(bool p, vector*& d)
{
  return -1;
}
    
/// Sets the specified row to the given vector
int denseMatrix::set_row(size_t i, vector*& r)
{
  return -1;
}

/// Sets the specified col to the given vector
int denseMatrix::set_col(size_t j, vector*& c)
{
  return -1;
}

/// Sets the specified diagonal to the given vector
int denseMatrix::set_diag(bool p, vector*&d)
{
  return -1;
}

/// Vector l_p norms for this matrix, p > 0
double denseMatrix::norm (double p)
{
  return -1;
}

/// Vector l_p norms, p is 'l1', 'l2', 'fro', 'inf'
double denseMatrix::norm (const char*  p)
{
  return -1;
}

/// r = a*row(i) + r
int    denseMatrix::row_daxpy(size_t i, double a, vector* r)
{
  return -1;
}

/// c = a*col(j) + c
int  denseMatrix::col_daxpy(size_t j, double a, vector* c)
{
  return -1;
}

/// Let r := this * x or  this^T * x depending on tranA
int denseMatrix::dot (bool transp, vector* x, vector*r)
{
  // Replace by call to BLAS library, in case it is available; otherwise
  // the under-optimized code below will be invoked

#ifdef _HAVE_CBLAS
  if (transp)
    cblas_dgemv(CblasRowMajor, CblasTrans, nrows(), ncols(), 
                1.0, data, ncols(), x->getData(), 1, 0.0, r->getData(), 1);
  else
    cblas_dgemv(CblasRowMajor, CblasNoTrans, nrows(), ncols(), 
                1.0, data, ncols(), x->getData(), 1, 0.0, r->getData(), 1);
#else
  double* pr = r->getData();
  double* px = x->getData();
  double* rp;
  // M*x
  if (!transp) {
    r->zeroOut();
    for (size_t i = 0; i < nrows(); i++) {
      rp = &data[i*ncols()];
      for (size_t j = 0; j < ncols(); j++) {
        pr[i] += rp[j] * px[j];
      }
    }
  } else {                      // M'*x
    r->zeroOut();
    for (size_t i = 0; i < nrows(); i++) {
      rp = &data[i*ncols()];
      for (size_t j = 0; j < ncols(); j++) {
        pr[j] += rp[j] * px[i];
      }
    }
  }
#endif
  return 0;
}
