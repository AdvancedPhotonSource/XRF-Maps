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


#include "denseMatrix.hpp"
#include "sparseMatrix.hpp"
#include <ctime>

namespace nsNNLS {
template <typename _T>
  class nnls {
    // The variables used during compute time
  private:                      
    vector<_T>* x;                  // The solution -- also current iterate
    vector<_T>* x0;                 // Starting value
    vector<_T>* oldx;               // Previous iterate
    vector<_T>* g;                  // Current gradient
    vector<_T>* oldg;               // Previous gradient
    vector<_T>* xdelta;             // x - oldx
    vector<_T>* gdelta;             // g - oldg
    vector<_T>* refx;               // iterate from 'M' steps ago
    vector<_T>* refg;               // reference gradient from M steps ago

    matrix<_T>* A;
    vector<_T>* b;
    vector<_T>* ax;                 // vector to hold A*x
    size_t* fset;               // fixed set 
    size_t fssize;              // sizeof fixed set

    // The parameters of the solver
  private:
    int maxit;               // maximum num. iters
    int   M;                    // max num. of null iterations
    _T decay;               // parameter to make diminishing scalar to decay by
    _T beta;                // diminishing scalar
    _T pgtol;               // projected gradient tolerance
    _T sigma;               // constant for descent condition

    // The solution and statistics variables
  private:
    struct out_ {
      clock_t start;
      size_t memory;
      vector<_T>* obj;
      int  iter;
      vector<_T>* time;
      vector<_T>* pgnorms;
      _T npg;                // inf-norm of projected gradient
    } out;

    // The helper functions used by the solver
  private:                      
    int    initialize()
    {
      size_t n = A->ncols();
      x = new vector<_T>(n);
      g = new vector<_T>(n);
      refx = new vector<_T>(n);
      refg = new vector<_T>(n);
      oldx = new vector<_T>(n);
      oldg = new vector<_T>(n);
      xdelta = new vector<_T>(n);
      gdelta =  new vector<_T>(n);
      out.memory += 8*n*sizeof(_T);

      ax = new vector<_T>(A->nrows());
      out.memory += sizeof(_T)*ax->length();

      fset = (size_t*) malloc(sizeof(size_t)*n);
      out.memory += sizeof(size_t)*n;

      x->setAll(.5);

      out.obj = new vector<_T>(maxit+1);
      out.pgnorms = new vector<_T>(maxit+1);
      out.time = new vector<_T>(maxit+1);

      out.memory += sizeof(_T)*3*(maxit+1);

      if (x0) {
        x->copy(x0);
        A->dot(matrix<_T>::NOTRAN, x, ax); ax->sub(b);
        A->dot(matrix<_T>::TRAN, ax, g);
      } else {
        // Initial gradient = A'(0 - b), since x0 = 0
        A->dot(matrix<_T>::TRAN, b, oldg);
        _T* dat = oldg->getData();
        for (size_t i = 0; i < oldg->length(); i++)
          dat[i] = -dat[i];
      }

      // old gradient = A'*(ax - b)
      A->dot(matrix<_T>::NOTRAN, x, ax);  ax->sub(b);
      A->dot(matrix<_T>::TRAN, ax, g);

      // Set the reference iterations
      refx->copy(x);
      refg->copy(g);

      return 0;
    }

    int    checkTermination()
    {
      int term = 0;
      if (out.iter >= maxit)
        term = 1;

      out.pgnorms->set(out.iter, normProjectedGradient());
      out.npg = out.pgnorms->get(out.iter);
      if (out.npg < pgtol)
        term = 2;
      return term;
    }

    void   showStatus()
    {
      fprintf(stderr, "%05d\t %010E\t %010E\n", out.iter, out.obj->get(out.iter), out.npg);
    }

    int    cleanUp()
    {
            delete ax;
            delete oldx;
            delete oldg;
            delete xdelta;
            delete gdelta;
            free(fset);
            delete refx;
            return 0;
    }

    void   findFixedVariables()
    {
      size_t n = A->ncols();
      _T *vars = x->getData();
      _T* grads = g->getData();

      // Clean out fixed set first
      memset(fset, 0, sizeof(size_t)*n);
      fssize = 0;

      for (size_t i = 0; i < n; i++) {
        if (vars[i] == 0 && grads[i] > 0)
          fset[fssize++] = i;
      }
    }

    void   computeXandGradDelta()
    {
      size_t n = A->ncols();
      _T *px = x->getData();
      _T *poldx = oldx->getData();
      _T* pg = g->getData();
      _T *poldg = oldg->getData();
      _T* pxd = xdelta->getData();
      _T* pgd = gdelta->getData();

      for (size_t i = 0; i < n; i++) {
        pxd[i] = px[i] - poldx[i];
        pgd[i] = pg[i] - poldg[i];
      }

      // Zero out fixed variables
      for (size_t i = 0; i < fssize; i++) {
        pxd[fset[i]] = 0;
        pgd[fset[i]] = 0;
      }

      // Save oldx and oldg
      memcpy(poldg, pg, n*sizeof(_T));
      memcpy(poldx, px, n*sizeof(_T));
    }

    void   computeObjGrad()
    {
      A->dot(matrix<_T>::NOTRAN, x, ax);
      ax->sub(b);        // ax = ax - b
      _T d;
      d = ax->norm2();
      _T tm  = (_T) (clock() - out.start) / CLOCKS_PER_SEC;
      out.obj->set(out.iter, 0.5*d*d);
      out.time->set(out.iter, tm);
      A->dot(matrix<_T>::TRAN, ax, g);          // A'(ax)
    }

    _T computeBBStep()
    {
      // Use xdelta and gdelta to compute BB step
      _T step = 0.0;
      _T nr;
      _T dr;

      if (out.iter % 2) {
        nr = xdelta->ddot(xdelta);
        dr = xdelta->ddot(gdelta);
      } else {
        nr = xdelta->ddot(gdelta);
        dr = gdelta->ddot(gdelta);
      }

      //fprintf(stderr, " nr / dr = %lf / %lf\n", nr, dr);
      if (nr == 0) return 0;
      step = nr / dr;

      step *= beta;
      return step;
    }

    void  checkDescentUpdateBeta()
    {

       // if(out.iter == 0)
       //     return;

      _T d = 0;
      _T* grad = g->getData();
      _T* prx  = refx->getData();
      _T* px   = x->getData();

      // compute sigma*<grad, refx - x>
      for (size_t i = 0; i < g->length(); i++)
        d += grad[i]*(prx[i]-px[i]);
      d *= sigma;

      if (out.iter >= M)
        d = out.obj->get(out.iter - M) - out.obj->get(out.iter) - d;
      else
        d = out.obj->get(out.iter) - d;

      //fprintf(stderr, "Descent value: %g - %g - %g, beta=%f\n", out.obj->get(out.iter -M), out.obj->get(out.iter - 2), d, beta*decay);

      // check if suff. descent failed
      if (d < 0) {
        beta *= decay;
      } else {                      // update reference iteration
        refx->copy(x);
        refg->copy(g);
      }
    }


    _T normProjectedGradient()
    {
      size_t n = A->ncols();
      _T* grads = g->getData();
      _T pg = 0.0;
      size_t ctr = 0;
      // compute the norm of the gradient for all variables not in fset
      for (size_t i = 0; i < n; i++) {
        if (i==fset[ctr]) {ctr++; continue;}
        pg = std::max(pg, std::abs(grads[i]));
      }
      return pg;
    }

    // The actual interface to the world!
  public:
    nnls() { x = 0; A = 0; b = 0; maxit = 0;}
  
    nnls (matrix<_T>* A, vector<_T>* b, int maxit) {
      this->x = 0; this->A = A; this->b = b;
      this->maxit = maxit; this->x0 = 0;
      out.obj = 0; out.iter = -1; out.time = 0;
      fset = 0; out.memory = 0;
      // convergence controlling parameters
      M = 100; beta = 1.0; decay = 0.9; pgtol = 1e-3;  sigma = .01;
    }

    nnls (matrix<_T>* A, vector<_T>* b, vector<_T>* x0, int maxit)  { nnls(A, b, maxit); this->x0 = x0;}

    ~nnls()
    {
      delete x;
      delete g;
      delete refg;
      delete out.obj;
      delete out.pgnorms;
      delete out.time;
    }

    // The various accessors and mutators (or whatever one calls 'em!)

    _T  getDecay() const   { return decay; }
    int     getM()     const   { return M;     }
    _T  getBeta()  const   { return beta;  }
    _T  getObj()           { return out.obj->get(out.iter-1);}
    _T  getPgTol() const   { return pgtol; }
    vector<_T>* getSolution()      { return x;     }
    size_t* getFset()          { return fset;  }
    size_t  getMaxit() const   { return maxit; }
    _T  getSigma() const   { return sigma; }

    void  setDecay(_T d) { decay = d;  }
    void  setM(int m)        { M = m;      }
    void  setBeta(_T b)  { beta = b;   }
    void  setPgTol(_T pg){ pgtol = pg; }
    void  setMaxit(size_t m) { maxit = m;  }
    void  setSigma(_T s) { sigma = s; }

    void  setData(matrix<_T>* A, vector<_T>* b)  { this->A = A; this->b = b;}

    // The functions that actually launch the ship, and land it!
    int     optimize()
    {
      out.start = clock();
      initialize();

      int term = 0;
      _T step;
      _T *px = x->getData();

      //fprintf(stderr, "Iter           Obj           ||g||\n");
      //fprintf(stderr, "----------------------------------\n");
      while (!term) {
        out.iter++;
        term = checkTermination();

        findFixedVariables();
        computeXandGradDelta();

        step = computeBBStep();

        x->scalePlusAdd(-step, g); // x = x - step*gradient

        // project
        for (size_t i = 0; i < x->length(); i++)
          if (px[i] < 0) px[i] = 0.0;

        computeObjGrad();
        // check the descent condition
        if (out.iter % M == 0) {
          checkDescentUpdateBeta();
        }
        //if (out.iter % 10 == 0)
        //  showStatus();
      }
      //showStatus();
      cleanUp();
      fprintf(stderr, "-");
      return out.iter;
    }

    int     saveStats(const char*fn)
    {
      FILE* fp = fopen(fn, "w");
      fprintf(fp, "NNLS: Solver (v) 1.0\n(c) 2010 Suvrit Sra\n");
      fprintf(fp, "Matrix = %zu bytes, rhs = %zu bytes\n", A->memoryUsage(), sizeof(_T)*A->nrows());
      fprintf(fp, "Extra Memory used=%zu bytes\nTime\t\tObj\n",out.memory);

      for (int i = 0; i < out.iter; i++) {
        fprintf(fp, "%.3lf\t%.8g\t%.8g\n", out.time->get(i), out.obj->get(i),out.pgnorms->get(i));
      }
      fclose(fp);
      return 0;
    }
    _T  getOptimizationTime() { return out.time->get(out.iter);}
  };

}
