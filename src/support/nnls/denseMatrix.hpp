// File: denseMatrix.h -*- c++ -*-
// Author: Suvrit Sra <suvrit@tuebingen.mpg.de>
// (c) Copyright 2010   Suvrit Sra
// Max-Planck-Institute for Biological Cybernetics

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

#ifndef _MY_DENSE_H
#define _MY_DENSE_H

#include "matrix.h"
#include <cassert>

namespace nsNNLS {

#define forall(x) for (size_t x = 0; x < size; x++)
template<typename _T>
  class denseMatrix : public matrix<_T> {
    size_t size;
    _T* data;
    bool external;              // is data managed externally
  public:
    denseMatrix(): matrix<_T>() { external = true; data = 0;}

    denseMatrix (size_t r, size_t c) : matrix<_T>(r, c)
    {
      assert (r > 0 && c > 0);
      data = new _T[r*c];
      size = r*c;
      external = false;
    }

    ~denseMatrix () { if (!external) delete[] data;}
    
    denseMatrix(size_t r, size_t c, _T* data) : matrix<_T>(r, c)
    { external = true; size=r*c; this->data=data;}


  public:
    /// Get the (i,j) entry of the matrix
    _T operator() (size_t i, size_t j) { return data[i*this->ncols()+j];}

    /// Get the (i,j) entry of the matrix
    _T get (size_t i, size_t j) { return data[i*this->ncols()+j];}

    /// Set the (i,j) entry of the matrix. Not all matrices can support this.
    void set (size_t i, size_t j, _T val) { data[i*this->ncols()+j] = val; }
    
    
    /// Returns 'r'-th row into pre-alloced vector
    int get_row (size_t, vector<_T>*&){return -1;}
    /// Returns 'c'-th col as a vector
    int get_col (size_t, vector<_T>*&){return -1;}
    /// Returns main or second diagonal (if p == true)
    int get_diag(bool p, vector<_T>*& d){return -1;}
    
    /// Sets the specified row to the given vector
    int set_row(size_t r, vector<_T>*&){return -1;}
    /// Sets the specified col to the given vector
    int set_col(size_t c, vector<_T>*&){return -1;}
    /// Sets the specified diagonal to the given vector
    int set_diag(bool p, vector<_T>*&){return -1;}

    /// Vector l_p norms for this matrix, p > 0
    _T norm (_T p){return (_T)-1;}
    /// Vector l_p norms, p is 'l1', 'l2', 'fro', 'inf'
    _T norm (const char*  p){return (_T)-1;}

    /// Apply an arbitrary function elementwise to this matrix. Changes the matrix.
    int apply (_T (* fn)(_T)) { forall(i) data[i] = fn(data[i]); return 0;}

    /// Scale the matrix so that x_ij := s * x_ij
    int scale (_T s) { forall(i) data[i] *= s; return 0;}

    /// Adds a const 's' so that x_ij := s + x_ij
    int add_const(_T s) { forall(i) data[i] += s; return 0;}

    /// r = a*row(i) + r
     int    row_daxpy(size_t i, _T a, vector<_T>* r){return -1;}
    /// c = a*col(j) + c
     int  col_daxpy(size_t j, _T a, vector<_T>* c) {return -1;}

    /// Let r := this * x or  this^T * x depending on tranA
    int dot (bool transp, vector<_T>* x, vector<_T>*r)
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
      _T* pr = r->getData();
      _T* px = x->getData();
      _T* rp;
      // M*x
      if (!transp) {
        r->zeroOut();
        for (size_t i = 0; i < this->nrows(); i++) {
          rp = &data[i*this->ncols()];
          for (size_t j = 0; j < this->ncols(); j++) {
            pr[i] += rp[j] * px[j];
          }
        }
      } else {                      // M'*x
        r->zeroOut();
        for (size_t i = 0; i < this->nrows(); i++) {
          rp = &data[i*this->ncols()];
          for (size_t j = 0; j < this->ncols(); j++) {
            pr[j] += rp[j] * px[i];
          }
        }
      }
    #endif
      return 0;
    }

    size_t memoryUsage() { return this->nrows()*this->ncols()*sizeof(_T);}
  };
}

#endif
