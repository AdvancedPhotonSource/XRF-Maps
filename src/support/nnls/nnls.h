// File: nnls.h -*- c++ -*-
// Author: Suvrit Sra
// Time-stamp: <08 March 2011 12:23:00 PM CET --  suvrit>
// Header file for the NNLS

// nnls.h - base class for nnls solver
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
#include "sparseMatrix.h"
#include <ctime>

namespace nsNNLS {
  class nnls {
    // The variables used during compute time
  private:                      
    vector* x;                  // The solution -- also current iterate
    vector* x0;                 // Starting value
    vector* oldx;               // Previous iterate
    vector* g;                  // Current gradient
    vector* oldg;               // Previous gradient
    vector* xdelta;             // x - oldx
    vector* gdelta;             // g - oldg
    vector* refx;               // iterate from 'M' steps ago
    vector* refg;               // reference gradient from M steps ago

    matrix* A;                  
    vector* b;
    vector* ax;                 // vector to hold A*x
    size_t* fset;               // fixed set 
    size_t fssize;              // sizeof fixed set

    // The parameters of the solver
  private:
    int maxit;               // maximum num. iters
    int   M;                    // max num. of null iterations
    double decay;               // parameter to make diminishing scalar to decay by
    double beta;                // diminishing scalar
    double pgtol;               // projected gradient tolerance
    double sigma;               // constant for descent condition

    // The solution and statistics variables
  private:
    struct out_ {
      clock_t start;
      size_t memory;
      vector* obj;
      int  iter;
      vector* time;
      vector* pgnorms;
      double npg;                // inf-norm of projected gradient
    } out;

    // The helper functions used by the solver
  private:                      
    int    initialize();            // computes / sets x0, g0, oldx, oldg, etc.
    int    checkTermination();      // embodies various termination criteria
    void   showStatus();            // 
    int    cleanUp();               // memory deallocation and friends
    void   findFixedVariables();    // compute fixed set (binding set)
    void   computeXandGradDelta();  //  
    void   computeObjGrad();        // compute both together to sav time
    double computeBBStep();         // compute step-size, take care of decay
    void  checkDescentUpdateBeta(); // check the wolfe-condition
    double normProjectedGradient(); // compute inf-norm of chopped gradient

    // The actual interface to the world!
  public:
    nnls() { x = 0; A = 0; b = 0; maxit = 0;}
  
    nnls (matrix* A, vector* b, int maxit) {
      this->x = 0; this->A = A; this->b = b;
      this->maxit = maxit; this->x0 = 0;
      out.obj = 0; out.iter = -1; out.time = 0;
      fset = 0; out.memory = 0;
      // convergence controlling parameters
      M = 100; beta = 1.0; decay = 0.9; pgtol = 1e-3;  sigma = .01;
    }

    nnls (matrix* A, vector* b, vector* x0, int maxit)  { nnls(A, b, maxit); this->x0 = x0;}

    // The various accessors and mutators (or whatever one calls 'em!)

    double  getDecay() const   { return decay; }
    int     getM()     const   { return M;     }
    double  getBeta()  const   { return beta;  }
    double  getObj()           { return out.obj->get(out.iter-1);}
    double  getPgTol() const   { return pgtol; }
    vector* getSolution()      { return x;     }
    size_t* getFset()          { return fset;  }
    size_t  getMaxit() const   { return maxit; }
    double  getSigma() const   { return sigma; }

    void  setDecay(double d) { decay = d;  }
    void  setM(int m)        { M = m;      }
    void  setBeta(double b)  { beta = b;   }
    void  setPgTol(double pg){ pgtol = pg; }
    void  setMaxit(size_t m) { maxit = m;  }
    void  setSigma(double s) { sigma = s; }

    void  setData(matrix* A, vector* b)  { this->A = A; this->b = b;}

    // The functions that actually launch the ship, and land it!
    int     optimize();
    int     saveStats(const char*fn);
    double  getOptimizationTime() { return out.time->get(out.iter);}
  };
}
