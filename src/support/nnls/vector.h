// File: vector.h        -*- c++ -*-
// Author: Suvrit Sra
// Time-stamp: <31 March 2010 04:04:01 PM CEST --  suvrit>


// vector.h - a simple vector class to have a cleaner interface
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



#ifndef _MY_VECTOR_H
#define _MY_VECTOR_H

#include <cstring>
#include <cmath>

// Vectors are dense -- usu. faster this way
// Does minimal error checking

namespace nsNNLS {
  class vector {
    size_t size;
    double* data;
    bool external;
#define forall(x) for (size_t x = 0; x < size; x++)

  public:
    vector() { external = true;}

    vector(size_t sz, double* v) { external = true; size = sz; data = v;}

    // Create a new allocated vector
    vector (size_t t) { external = false; size = t; data = new double[size]; memset(data, 0, sizeof(double)*size); }

    ~vector() { if (!external) delete[] data;}

    void   setSize(size_t s) { size = s;}
    size_t length() const { return size;}
    void   setAll(double c) { for (size_t i = 0; i < size; i++) data[i] = c;}
    void   zeroOut() { memset(data, 0, sizeof(double)*size);}
    double* getData() { return data;}
    double operator ()(vector* v, size_t j) { return v->data[j];}
    double  get(size_t j) { return data[j]; }

    void set(size_t j, double v) { data[j] = v;}
    // this'*b
    double ddot(vector*b) { double *pb = b->getData(); double r = 0; forall(i) r += data[i]*pb[i]; return r;}

    // norm(this) --- not doing blas style that avoids overflow
    double norm2() { double n = 0.0; forall(i) n += data[i]*data[i]; return sqrt(n);}

    double norm1() { double n = 0.0; forall(i) n += fabs(data[i]); return n;}

    double norminf() { double n = 0.0; double t; forall(i) {t = fabs(data[i]); if (t > n) n = t;} return n;}

    // add: this = this + b
    void add(vector* b)
    { double* bd = b->getData();
      forall(i) data[i] += bd[i];
    }

    void sub(vector* b)
    { double* pb = b->getData();
      forall(i) data[i] -= pb[i];
    }

    void copy(vector* b) { double* pb = b->getData(); memcpy(data, pb, sizeof(double)*size);}

    // x = x + a*y (x == this is used)
    void scalePlusAdd(double a, vector* y) 
    { if (a == 0) add(y); 
      else {
        double* py = y->getData();
        forall(i) data[i] = data[i] + a*py[i];
      }
    }
  };
}

#endif 
