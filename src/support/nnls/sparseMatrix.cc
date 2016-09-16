// File: sparseMatrix.cc -- implements sparse matrix functionality 

// Author: Suvrit Sra <suvrit@tuebingen.mpg.de>
// (c) Copyright 2010   Suvrit Sra
// Max-Planck-Institute for Biological Cybernetics

// sparseMatrix.cc - implements sparse matrix functionality
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

#include "sparseMatrix.h"
#include <cstdio>
#include <string>

using namespace nsNNLS;

int sparseMatrix::load(const char* fn, bool asbin)
{
  if (asbin)
    return load_as_bin(fn);
  else 
    return load_as_txt(fn);
  return 0;
}

int sparseMatrix::load_as_bin(const char* fn)
{
  FILE* fp = fopen (fn, "rb");
  if (fp == 0) {
    fprintf(stderr, "sparseMatrix:: Error opening `%s' for reading\n", fn);
    return -1;
  }

  size_t m;
  size_t n;
  size_t nz;

  // read the dimension info
  if (fread(&m,  sizeof(size_t), 1, fp) != 1) {
    fprintf(stderr, "sparseMatrix:: Error reading M (numrows) from %s\n", fn);
    return -2;
  }

  if (fread(&n,  sizeof(size_t), 1, fp) != 1) {
    fprintf(stderr, "sparseMatrix:: Error reading N (numcols) from %s\n", fn);
    return -2;
  }

  if (fread(&nz, sizeof(size_t), 1, fp) != 1) {
    fprintf(stderr, "sparseMatrix:: Error reading #NNZ (nnz) from %s\n", fn);
    return -2;
  }

  setsize(m, n);
  nnz = nz; 
  size = nnz;

  // Allocate the arrays
  cols = new size_t[n+1];
  ridx = new size_t[nz];
  data = new double[nz];
   
  // Read the colptrs
  if (fread(cols, sizeof(size_t), n+1, fp) != n+1) {
    fprintf(stderr, "sparseMatrix:: Error reading colptrs from %s\n", fn);
    return -3;
  }

  // Read the row indices
  // if nz happens to be very big, we should break up the reads -- TODO
  if (fread(ridx, sizeof(size_t), nz, fp) != nz) {
    fprintf(stderr, "sparseMatrix:: Error reading rowindices from %s\n", fn);
    return -3;
  }

  // Read the nonzeros themselves
  if (fread(data, sizeof(double), nz, fp) != nz) {
    fprintf(stderr, "sparseMatrix:: Error reading nonzeros from %s\n", fn);
    return -3;
  }

  fclose(fp);
  return 0;
}

int sparseMatrix::load_as_txt(const char* fn)
{
  std::string prefix(fn);
  // Set up the file names that we are gonna open
  std::string dim(prefix + "_dim");
  std::string rows_file(prefix + "_row_ccs");
  std::string cols_file(prefix + "_col_ccs");
  std::string nz_file (prefix + "_txx_nz");
  
  fprintf(stderr, "sparseMatrix: Nonzeros will be loaded from %s\n", nz_file.c_str());
  
  FILE *fp;
  fp = fopen(dim.c_str(), "r");
  if (!fp) {
    fprintf(stderr, "sparseMatrix:: Could not open `%s'", dim.c_str());
    return -1;
  }
  
  size_t m, n;  size_t nz;
  ssize_t r;
  r = fscanf(fp, "%zu %zu %zu", &m, &n, &nz);
  fclose(fp);

  fprintf(stderr, "Found (%zu,%zu,%zu) matrix\n",m,n,nz);
  setsize(m, n);
  nnz = nz;   size = nnz;

  // Allocate the arrays
  cols = new size_t[n+1];
  ridx = new size_t[nz];
  data = new double[nz];

  cols[n] = nz; // Add the convenience value

  // Now read in all the columnpointers
  fp = fopen(cols_file.c_str(), "r");
  if (!fp) {
    fprintf(stderr, "sparseMatrix:: Could not open `%s'", cols_file.c_str());
    return -2;
  }

  for (size_t i = 0; i < n; i++) 
    r = fscanf(fp, "%zu", &cols[i]);
  fclose(fp);
  
  // Now read in all the rowindices
  fp = fopen(rows_file.c_str(), "r");
  if (!fp) {
    fprintf(stderr, "sparseMatrix:: Could not open `%s'", rows_file.c_str());
    return -3;
  }
  for (size_t i = 0; i <nz; i++)
    r = fscanf(fp, "%zu", &ridx[i]);

  fclose(fp);
  
  // Now read in all the values
  fp = fopen(nz_file.c_str(), "r");
  if (!fp) {
    fprintf(stderr, "sparseMatrix:: Could not open `%s'", nz_file.c_str());
    return -4;
  }
  for (size_t i = 0; i <nz; i++)
    r = fscanf(fp, "%lf", &data[i]);
  fclose(fp);

  return 0;
}

/// Get the (i,j) entry of the matrix
double sparseMatrix::operator()   (size_t i, size_t j)
{
  // size_t sz = cols[j+1]-cols[j];
  // int idx = binary_search(m_rowindx + m_colptrs[j], i, sz);
  // if (idx != -1) 
  //   return (m_val + m_colptrs[j])[idx];
  return 0.0;
}

/// Get the (i,j) entry of the matrix
double sparseMatrix::get (size_t i, size_t j)
{
  // size_t sz = cols[j+1]-cols[j];
  // int idx = binary_search(m_rowindx + m_colptrs[j], i, sz);
  // if (idx != -1) 
  //   return (m_val + m_colptrs[j])[idx];
  return 0.0;
}

/// Set the (i,j) entry of the matrix. If entry does not exist, function bombs.
int sparseMatrix::set (size_t i, size_t j, double val)
{
  return -1;
}
    
/// Returns 'r'-th row into pre-alloced vector
int sparseMatrix::get_row (size_t i, vector*& r)
{
  return -1;
}

/// Returns 'c'-th col as a vector
int sparseMatrix::get_col (size_t j, vector*& c)
{
  return -1;
}

/// Returns main or second diagonal (if p == true)
int sparseMatrix::get_diag(bool p, vector*& d) 
{
  return -1;
}

/// Vector l_p norms for this matrix, p > 0
double sparseMatrix::norm (double p)
{
  return -1;
}

/// Vector l_p norms, p is 'l1', 'l2', 'fro', 'inf'
double sparseMatrix::norm (const char*  p)
{
  return -1;
}

/// r = a*row(i) + r
int  sparseMatrix::row_daxpy(size_t i, double a, vector* r)
{
  return -1;
}

/// c = a*col(j) + c
int  sparseMatrix::col_daxpy(size_t j, double a, vector* c)
{
  return -1;
}

/// Let r := this * x or  this^T * x depending on tranA
int sparseMatrix::dot (bool transp, vector* x, vector*r)
{
  double* px = x->getData();
  double* pr = r->getData();

  if (!transp) {
    r->zeroOut();
    for (size_t i = 0; i < ncols(); i++) {
      for (size_t j = cols[i]; j < cols[i+1]; j++)
        pr[ridx[j]] += data[j] * px[i];
    }
  } else {
    double yi;
    for (size_t i = 0; i < ncols(); i++) {
      yi = 0;
      for (size_t j = cols[i]; j < cols[i+1]; j++) {
        yi += data[j] * px[ ridx[j] ];
      }
      pr[i] = yi;  
    }
  }
  return 0;
}
