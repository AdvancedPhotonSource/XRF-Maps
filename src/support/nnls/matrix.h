// File: matrix.h        -*- c++ -*-
// Author: Suvrit Sra
// Time-stamp: <04 March 2011 07:35:51 PM CET --  suvrit>
// Base class for both sparse and dense matrices....

// matrix.h  - minimal base class for dense and sparse matrices
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


#ifndef _MY_MATRIX_BASE_H
#define _MY_MATRIX_BASE_H

#include <cstring>

#include "vector.h"

namespace nsNNLS {
  class matrix {
    size_t M;               /// Num of rows
    size_t N;               /// Num of colums

  protected:
    size_t memory;
  public:

    static const bool TRAN = true;
    static const bool NOTRAN = false;

    matrix() {}

    /// Create an empty (unallocated) matrix of specified size.
    matrix (size_t rows, size_t cols) { M = rows;   N = cols;   }
    
    virtual ~matrix() {;}


    virtual int load(const char* f, bool asbin) = 0;

    void setsize(size_t r, size_t c) { M = r; N = c;}
    /// Returns number of rows. 
    size_t nrows() const { return M;}

    /// Returns number of colums
    size_t ncols() const { return N;}

    /// Sets the row and column dimensionality of the matrix
    void matrix_setsize(size_t m, size_t n) { M = m; N = n; }

    /// Get the (i,j) entry of the matrix
    virtual double operator()   (size_t i, size_t j            ) = 0;
    /// Get the (i,j) entry of the matrix
    virtual double get          (size_t i, size_t j            ) = 0;
    /// Set the (i,j) entry of the matrix. Not all matrices can support this.
    virtual int set            (size_t i, size_t j, double val) = 0;

    /// Returns 'r'-th row as a vector
    virtual int get_row (size_t, vector*&) = 0;
    /// Returns 'c'-th col as a vector
    virtual int get_col (size_t, vector*&) = 0;
    /// Returns main or second diagonal (if p == true)
    virtual int get_diag(bool p, vector*& d  ) = 0;
    
    /// Sets the specified row to the given vector
    virtual int set_row(size_t r, vector*&) = 0;
    /// Sets the specified col to the given vector
    virtual int set_col(size_t c, vector*&) = 0;
    /// Sets the specified diagonal to the given vector
    virtual int set_diag(bool p, vector*&)  = 0;

    /// Vector l_p norms for this matrix, p > 0
    virtual double norm (double p) = 0;
    /// Vector l_p norms, p is 'l1', 'l2', 'fro', 'inf'
    virtual double norm (const char*  p) = 0;

    /// Apply an arbitrary function elementwise to this matrix. Changes the matrix.
    virtual int apply (double (* fn)(double)) = 0;

    /// Scale the matrix so that x_ij := s * x_ij
    virtual int scale (double s) = 0;

    /// Adds a const 's' so that x_ij := s + x_ij
    virtual int add_const(double s) = 0;

    /// r = a*row(i) + r
    virtual int    row_daxpy(size_t i, double a, vector* r) = 0;
    /// c = a*col(j) + c
    virtual int  col_daxpy(size_t j, double a, vector* c) = 0;

    /// Let r := this * x or  this^T * x depending on tranA
    virtual int dot (bool transp, vector* x, vector*r) = 0;

    /// get statistics about storage
    virtual size_t memoryUsage() = 0;
  };

}
#endif
