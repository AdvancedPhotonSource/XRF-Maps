// File: sparseMatrix.h -*- c++ -*-
// Author: Suvrit Sra <suvrit@tuebingen.mpg.de>
// (c) Copyright 2010   Suvrit Sra
// Max-Planck-Institute for Biological Cybernetics

// sparseMatrix.h -  CCS sparse matrix class (double precision)
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

#ifndef _MY_SPARSE_H
#define _MY_SPARSE_H

#include "matrix.h"
#include <cassert>

#define forall(x) for (size_t x = 0; x < size; x++)

namespace nsNNLS {
  template<typename _T>
  class sparseMatrix : public matrix<_T> {
  private:
    _T* data;
    size_t* ridx;
    size_t* cols;
    size_t size;
    size_t nnz;
    bool external;
  public:
    sparseMatrix() : matrix<_T>() {external = true;}
    sparseMatrix (size_t r, size_t c, size_t nnz) : matrix<_T>(r, c) {
      assert (r > 0 && c > 0);
      this->nnz = nnz; this->size = nnz;
      cols = new size_t [c+1];
      ridx = new size_t [nnz];
      data = new _T [nnz];
    }

    sparseMatrix(size_t r, size_t c, size_t nnz, size_t* ridx, size_t* cptr, _T* val) : matrix<_T>(r, c)
    { this->nnz = nnz; this->size = nnz; this->ridx = ridx; this->cols = cptr; this->data = val;}

    ~sparseMatrix () { if (!external) { delete [] ridx; delete[] cols; delete [] data;}}

  public:

    /// Get the (i,j) entry of the matrix
    _T operator()   (size_t i, size_t j)
    {
      // size_t sz = cols[j+1]-cols[j];
      // int idx = binary_search(m_rowindx + m_colptrs[j], i, sz);
      // if (idx != -1)
      //   return (m_val + m_colptrs[j])[idx];
      return 0.0;
    }

    /// Get the (i,j) entry of the matrix
    _T get (size_t i, size_t j)
    {
      // size_t sz = cols[j+1]-cols[j];
      // int idx = binary_search(m_rowindx + m_colptrs[j], i, sz);
      // if (idx != -1)
      //   return (m_val + m_colptrs[j])[idx];
      return 0.0;
    }

    size_t getNnz() const { return nnz;}

    /// Set the (i,j) entry of the matrix. If entry does not exist, function bombs.
    int set (size_t i, size_t j, _T val) {  return -1;  }
    
    /// Returns 'r'-th row into pre-alloced vector
    int get_row (size_t, vector<_T>*&){  return -1;  }

    /// Returns 'c'-th col as a vector
    int get_col (size_t, vector<_T>*&){  return -1;  }

    /// Returns main or second diagonal (if p == true)
    int get_diag(bool p, vector<_T>*& d){  return -1;  }

    
    /// Sets the specified row to the given vector
    int set_row(size_t r, vector<_T>*&) { assert("Unsupported"); return -1;}

    /// Sets the specified col to the given vector
    int set_col(size_t c, vector<_T>*&) { assert("Unsupported"); return -1;}

    /// Sets the specified diagonal to the given vector
    int set_diag(bool p, vector<_T>*&) { assert("Unsupported"); return -1;}

    /// Vector l_p norms for this matrix, p > 0
    _T norm (_T p){  return -1;  }
    /// Vector l_p norms, p is 'l1', 'l2', 'fro', 'inf'
    _T norm (const char*  p){  return -1;  }

    /// Apply an arbitrary function elementwise to this matrix. Changes the matrix.
    int apply (_T (* fn)(_T)) { forall(i) data[i] = fn(data[i]); return 0;};

    /// Scale the matrix so that x_ij := s * x_ij
    int scale (_T s) { forall(i) data[i] *= s; return 0;}

    /// Adds a const 's' so that x_ij := s + x_ij
    int add_const(_T s) { forall(i) data[i] += s; return 0;};

    /// r = a*row(i) + r
    int    row_daxpy(size_t i, _T a, vector<_T>* r){  return -1;  }
    /// c = a*col(j) + c
    int  col_daxpy(size_t j, _T a, vector<_T>* c){  return -1;  }

    /// Let r := this * x or  this^T * x depending on tranA
    int dot (bool transp, vector<_T>* x, vector<_T>*r)
    {
      _T* px = x->getData();
      _T* pr = r->getData();

      if (!transp) {
        r->zeroOut();
        for (size_t i = 0; i < this->ncols(); i++) {
          for (size_t j = cols[i]; j < cols[i+1]; j++)
            pr[ridx[j]] += data[j] * px[i];
        }
      } else {
        _T yi;
        for (size_t i = 0; i < this->ncols(); i++) {
          yi = 0;
          for (size_t j = cols[i]; j < cols[i+1]; j++) {
            yi += data[j] * px[ ridx[j] ];
          }
          pr[i] = yi;
        }
      }
      return 0;
    }


    size_t memoryUsage() { this->memory = nnz*sizeof(_T) + nnz*sizeof(size_t) + (this->ncols()+1)*sizeof(size_t); return this->memory;}
  };
}

#endif 
