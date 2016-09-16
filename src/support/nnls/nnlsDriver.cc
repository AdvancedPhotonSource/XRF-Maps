// File: nnlsDriver.cc
// Author: Suvrit Sra
// Time-stamp: <08 March 2011 12:19:07 PM CET --  suvrit>

// nnlsDriver - Driver routine to invoke libnnls
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

#include <ctime>
#include <cstdlib>
#include <string>
#include <cstdio>
#include <cstdarg>

#include "nnls.h"


int error(int ve, FILE* fp, const char* fmt, ...)
{
  if (ve > 0) {
    va_list ap;
    va_start(ap, fmt);
    vfprintf(fp, fmt, ap);
    va_end(ap);
  } 
  return 0;
}

void showHelp(char* msg) {
  error(1, stderr, "NNLS Driver version 0.1 -- Suvrit Sra <suvrit@gmail.com>\n");
  error(1, stderr, "Usage: %s <dD|sS|> A b maxit\n", msg);
  error(1, stderr, "Example: ./nnls S A b 10\nSolves NNLS(A, b) for a SPARSE matrix A, loaded in BIN format\n");
  error(1, stderr, "Upper case for BIN mode matrices, lower case for txt\n");
}

using namespace nsNNLS;

matrix* load (const char* fname, char type);
vector* load_vector(const char* fil, size_t sz, char dos);
int     write_vector(const char* fil, vector* v);
int     write_binvector(const char* fil, vector* v);

int main (int argc, char** argv) 
{
  matrix* A;
  vector* b;
  vector* x;
  nnls*   solver;
  int     flag;

  if (argc != 5) {
    showHelp(argv[0]);
     return -1;
  }

  char dos = argv[1][0];
  error(1, stderr, "Loading...`%s'\n", argv[2]);
  A = load(argv[2], dos);
  error(1, stderr, "Loaded: (%zu,%zu)\n", A->nrows(), A->ncols());
  b = load_vector(argv[3], A->nrows(), dos);
  solver = new nnls(A, b, atoi(argv[4]));

  error(1, stderr, "Optimizing...\n");
  flag = solver->optimize();
  error(1, stderr, "Done!\n");

  error(1, stderr, "Optimization time: %.3lf seconds\n", solver->getOptimizationTime());

  if (flag < 0) {
    error(1, stderr, "NNLS: Solver terminated with error flag: %d\n", flag);
    return flag;
  }

  x = solver->getSolution();

  std::string solfile = std::string(argv[2]) + ".solution";
  write_binvector(solfile.c_str(), x);
  std::string solstat = std::string(argv[2]) + ".stats";
  solver->saveStats(solstat.c_str());
  error(1, stderr, "Solution left in: %s\n", solfile.c_str());
  error(1, stderr, "Statistics left in: %s\n", solstat.c_str());
  return 0;
}

matrix* load (const char* fname, char type)
{
  matrix* mat = 0;
  bool asbin = false;
  char* mtype;

  switch (type) {
  case 'D': 
    mat = new denseMatrix();
    asbin = true;
    mtype = "Dense BIN matrix";
    break;
  case 'd':
    mat = new denseMatrix();
    mtype = "Dense TXT matrix";
    break;
  case 'S': 
    mat = new sparseMatrix();
    mtype = "Sparse BIN matrix";
    asbin = true;
    break;
  case 's':
    mat = new sparseMatrix();
    mtype = "Sparse TXT matrix";
    break;
  default:
    error(1, stderr, "Unrecognized matrix type requested\n");
    return 0;
  }
  error(1, stderr, "Trying to load matrix of type: `%s'\n", mtype);
  if (mat->load(fname, asbin) < 0) { delete mat; mat = 0;}
  return mat;
}

vector* load_vector(const char* fn, size_t sz, char dos)
{
  vector* v = new vector(sz);
  FILE* fp;
  bool isdense = false;
  if (dos == 's' || dos == 'd') isdense = true;
  if (isdense) {
    fp = fopen(fn, "r");
    if (!fp) {
      error(1, stderr, "could not load vector (length %zu) from: `%s'\n", sz, fn);
      delete v;
      return 0;
    }
    double *pv = v->getData();
    size_t i = 0;

    size_t m, n;
    ssize_t rv;
    // read size info, and skip
    rv = fscanf(fp, "%zu %zu", &m, &n);
    if (n != 1 || m != sz) {
      error(1, stderr, "could not load vector (length %zu) from: `%s'\n", sz, fn);
      delete v;
      return 0;
    }
    while (i < sz) {
      rv = fscanf(fp, "%lf", &pv[i]);
      i++;
    }

    fclose(fp);
  } else {
    fp = fopen(fn, "rb");
    if (!fp) {
      error(1, stderr, "could not load vector (length %zu) from: `%s'\n", sz, fn);
      delete v;
      return 0;
    }
    double* pv = v->getData();
    size_t m, n, nb;
    nb=fread(&m, sizeof(size_t), 1, fp);
    nb=fread(&n, sizeof(size_t), 1, fp);
    if (n != 1) {
      error(1, stderr, "expected vector, found object of size (%zu %zu) in %s\n", m, n, fn);
      delete v;
      return 0;
    }
    error(1, stderr, "found vector (%d,%d) in %s\n", m, n, fn);
    nb=fread(pv, sizeof(double), m, fp);
    fclose(fp);
  }
  //fprintf(stderr, "v[0] = %lf\t v[end]=%lf\n", v->get(0), v->get(v->length()-1));
  return v;
}

int write_binvector(const char* fn, vector* v)
{
  FILE* fp = fopen(fn, "wb");
  if (!fp) {
    error(1, stderr, "could not write vector to: %s\n", fn);
    return -1;
  }

  size_t m, n;
  m = v->length();
  n = 1;

  double *pv = v->getData();
  fwrite(&m, sizeof(size_t), 1, fp);
  fwrite(&n, sizeof(size_t), 1, fp);
  fwrite(pv, sizeof(double), v->length(), fp);
  fclose(fp);
  return 0;
}


int write_vector(const char* fn, vector* v)
{
  FILE* fp = fopen(fn, "w");
  if (!fp) {
    error(1, stderr, "could not write vector to: %s\n", fn);
    return -1;
  }

  double *pv = v->getData();
  for (size_t i = 0; i < v->length(); i++)
    fprintf(fp, "%lf\n", pv[i]);

  fclose(fp);
  return 0;
}
