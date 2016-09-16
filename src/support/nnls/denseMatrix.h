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

  class denseMatrix : public matrix {
    size_t size;
    double* data;
    bool external;              // is data managed externally
  public:
    denseMatrix() { external = true; data = 0;}

    denseMatrix (size_t r, size_t c) : matrix(r, c) { 
      assert (r > 0 && c > 0);
      data = new double[r*c];
      size = r*c;
      external = false;
    }

    ~denseMatrix () { if (!external) delete[] data;}
    
    denseMatrix(size_t r, size_t c, double* data) : matrix(r, c) 
    { external = true; size=r*c; this->data=data;}

    int load(const char* fn, bool asbin);

  private:
    int load_as_bin(const char*);
    int load_as_txt(const char*);
  public:
    /// Get the (i,j) entry of the matrix
    double operator()   (size_t i, size_t j) { return data[i*ncols()+j];}

    /// Get the (i,j) entry of the matrix
    double get (size_t i, size_t j) { return data[i*ncols()+j];}

    /// Set the (i,j) entry of the matrix. Not all matrices can support this.
    int set (size_t i, size_t j, double val) { data[i*ncols()+j] = val; return 0;}
    
    
    /// Returns 'r'-th row into pre-alloced vector
    int get_row (size_t, vector*&);
    /// Returns 'c'-th col as a vector
    int get_col (size_t, vector*&);
    /// Returns main or second diagonal (if p == true)
    int get_diag(bool p, vector*& d);
    
    /// Sets the specified row to the given vector
    int set_row(size_t r, vector*&);
    /// Sets the specified col to the given vector
    int set_col(size_t c, vector*&);
    /// Sets the specified diagonal to the given vector
    int set_diag(bool p, vector*&);

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

    size_t memoryUsage() { return nrows()*ncols()*sizeof(double);}
  };
}

#endif
