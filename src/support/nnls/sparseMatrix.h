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
  
  class sparseMatrix : public matrix {
  private:
    double* data;
    size_t* ridx;
    size_t* cols;
    size_t size;
    size_t nnz;
    bool external;
  public:
    sparseMatrix() {external = true;}
    sparseMatrix (size_t r, size_t c, size_t nnz) : matrix(r, c) {
      assert (r > 0 && c > 0);
      this->nnz = nnz; this->size = nnz;
      cols = new size_t [c+1];
      ridx = new size_t [nnz];
      data = new double [nnz];
    }

    sparseMatrix(size_t r, size_t c, size_t nnz, size_t* ridx, size_t* cptr, double* val) : matrix(r, c)
    { this->nnz = nnz; this->size = nnz; this->ridx = ridx; this->cols = cptr; this->data = val;}

    int load(const char* fn, bool asbin);

    ~sparseMatrix () { if (!external) { delete [] ridx; delete[] cols; delete [] data;}}


  private:
    int load_as_bin(const char* fn);
    int load_as_txt(const char* fn);
  public:

    /// Get the (i,j) entry of the matrix
    double operator()   (size_t i, size_t j);

    /// Get the (i,j) entry of the matrix
    double get (size_t i, size_t j);

    size_t getNnz() const { return nnz;}
    /// Set the (i,j) entry of the matrix. If entry does not exist, function bombs.
    int set (size_t i, size_t j, double val);
    
    /// Returns 'r'-th row into pre-alloced vector
    int get_row (size_t, vector*&);
    /// Returns 'c'-th col as a vector
    int get_col (size_t, vector*&);
    /// Returns main or second diagonal (if p == true)
    int get_diag(bool p, vector*& d);
    
    /// Sets the specified row to the given vector
    int set_row(size_t r, vector*&) { assert("Unsupported"); return -1;}

    /// Sets the specified col to the given vector
    int set_col(size_t c, vector*&) { assert("Unsupported"); return -1;}

    /// Sets the specified diagonal to the given vector
    int set_diag(bool p, vector*&) { assert("Unsupported"); return -1;}

    /// Vector l_p norms for this matrix, p > 0
    double norm (double p);
    /// Vector l_p norms, p is 'l1', 'l2', 'fro', 'inf'
    double norm (const char*  p);

    /// Apply an arbitrary function elementwise to this matrix. Changes the matrix.
    int apply (double (* fn)(double)) { forall(i) data[i] = fn(data[i]); return 0;};

    /// Scale the matrix so that x_ij := s * x_ij
    int scale (double s) { forall(i) data[i] *= s; return 0;}

    /// Adds a const 's' so that x_ij := s + x_ij
    int add_const(double s) { forall(i) data[i] += s; return 0;};

    /// r = a*row(i) + r
    int    row_daxpy(size_t i, double a, vector* r);
    /// c = a*col(j) + c
    int  col_daxpy(size_t j, double a, vector* c);

    /// Let r := this * x or  this^T * x depending on tranA
    int dot (bool transp, vector* x, vector*r);


    size_t memoryUsage() { memory = nnz*sizeof(double) + nnz*sizeof(size_t) + (ncols()+1)*sizeof(size_t); return memory;}
  };
}

#endif 
