// File: nnls.cc 
// Author: Suvrit Sra
// Time-stamp: <08 March 2011 01:40:04 PM CET --  suvrit>
// nnls solver

// Copyright (C) 2009, 2010 Suvrit Sra (suvrit@tuebingen.mpg.de)
// Copyright Max-Planck-Institute for Biological Cybernetics

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


#include "nnls.h"
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <algorithm>

using namespace nsNNLS;

int nnls::optimize()
{
  out.start = clock();
  initialize();

  int term = 0;
  double step;
  double *px = x->getData();

  fprintf(stderr, "Iter           Obj           ||g||\n");
  fprintf(stderr, "----------------------------------\n");
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
    if (out.iter % 10 == 0)
      showStatus();
  }
  showStatus();
  cleanUp();
  return 0;
}

int nnls::saveStats(const char* fn)
{
  FILE* fp = fopen(fn, "w");
  fprintf(fp, "NNLS: Solver (v) 1.0\n(c) 2010 Suvrit Sra\n");
  fprintf(fp, "Matrix = %zu bytes, rhs = %zu bytes\n", A->memoryUsage(), sizeof(double)*A->nrows());
  fprintf(fp, "Extra Memory used=%zu bytes\nTime\t\tObj\n",out.memory);

  for (int i = 0; i < out.iter; i++) {
    fprintf(fp, "%.3lf\t%.8g\t%.8g\n", out.time->get(i), out.obj->get(i),out.pgnorms->get(i));
  }
  fclose(fp);
  return 0;
}

void nnls::computeObjGrad()
{
  A->dot(matrix::NOTRAN, x, ax);
  ax->sub(b);        // ax = ax - b
  double d;
  d = ax->norm2();
  double tm  = (double) (clock() - out.start) / CLOCKS_PER_SEC; 
  out.obj->set(out.iter, 0.5*d*d);
  out.time->set(out.iter, tm);
  A->dot(matrix::TRAN, ax, g);          // A'(ax)
}

void nnls::findFixedVariables()
{
  size_t n = A->ncols();
  double *vars = x->getData();
  double* grads = g->getData();

  // Clean out fixed set first
  memset(fset, 0, sizeof(size_t)*n);
  fssize = 0;

  for (size_t i = 0; i < n; i++) {
    if (vars[i] == 0 && grads[i] > 0)
      fset[fssize++] = i;
  }
}

void nnls::computeXandGradDelta()
{
  size_t n = A->ncols();
  double *px = x->getData();
  double *poldx = oldx->getData();
  double* pg = g->getData();
  double *poldg = oldg->getData();
  double* pxd = xdelta->getData();
  double* pgd = gdelta->getData();

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
  memcpy(poldg, pg, n*sizeof(double));
  memcpy(poldx, px, n*sizeof(double));
}

// Have to yet do error checking for too-small or too-large step-sizes
double nnls::computeBBStep()
{
  // Use xdelta and gdelta to compute BB step
  double step = 0.0;
  double nr;
  double dr;  

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

void nnls::checkDescentUpdateBeta()
{

  double d = 0;
  double* grad = g->getData();
  double* prx  = refx->getData();
  double* px   = x->getData();

  // compute sigma*<grad, refx - x>
  for (size_t i = 0; i < g->length(); i++)
    d += grad[i]*(prx[i]-px[i]);
  d *= sigma;
  
  d = out.obj->get(out.iter - M) - out.obj->get(out.iter) - d; 

  //fprintf(stderr, "Descent value: %g - %g - %g, beta=%f\n", out.obj->get(out.iter -M), out.obj->get(out.iter - 2), d, beta*decay);

  // check if suff. descent failed
  if (d < 0) {
    beta *= decay;
  } else {                      // update reference iteration
    refx->copy(x);
    refg->copy(g);
  }
}

int nnls::initialize()
{
  size_t n = A->ncols();
  x = new vector(n);
  g = new vector(n);
  refx = new vector(n);
  refg = new vector(n);
  oldx = new vector(n);
  oldg = new vector(n);
  xdelta = new vector(n);
  gdelta =  new vector(n);
  out.memory += 8*n*sizeof(double);

  ax = new vector(A->nrows());
  out.memory += sizeof(double)*ax->length();

  fset = (size_t*) malloc(sizeof(size_t)*n);
  out.memory += sizeof(size_t)*n;

  x->setAll(.5);

  out.obj = new vector(maxit+1);
  out.pgnorms = new vector(maxit+1);
  out.time = new vector(maxit+1);

  out.memory += sizeof(double)*3*(maxit+1);

  if (x0) {
    x->copy(x0);
    A->dot(matrix::NOTRAN, x, ax); ax->sub(b);
    A->dot(matrix::TRAN, ax, g);
  } else {
    // Initial gradient = A'(0 - b), since x0 = 0
    A->dot(matrix::TRAN, b, oldg);
    double* dat = oldg->getData();
    for (size_t i = 0; i < oldg->length(); i++)
      dat[i] = -dat[i];
  }

  // old gradient = A'*(ax - b)
  A->dot(matrix::NOTRAN, x, ax);  ax->sub(b);
  A->dot(matrix::TRAN, ax, g);

  // Set the reference iterations
  refx->copy(x);
  refg->copy(g);

  return 0;
}

int nnls::checkTermination()
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

double nnls::normProjectedGradient()
{
  size_t n = A->ncols();
  double* grads = g->getData();
  double pg = 0.0;
  size_t ctr = 0;
  // compute the norm of the gradient for all variables not in fset
  for (size_t i = 0; i < n; i++) {
    if (i==fset[ctr]) {ctr++; continue;}
    pg = std::max(pg, fabs(grads[i]));
  }
  return pg;
}

void nnls::showStatus()
{
  fprintf(stderr, "%05d\t %010E\t %010E\n", out.iter, out.obj->get(out.iter), out.npg);
}

int nnls::cleanUp()
{
  delete ax;
  delete oldx;
  delete oldg;
  delete xdelta;
  delete gdelta;
  delete fset;
  delete refx;
  return 0;
}
