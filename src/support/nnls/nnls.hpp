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

// --------
// Arthur Glowacki
// Argonne National Lab
// Dec 2017 : Modified to make it template class and use Eigen data structures
#include <Eigen/Core>

namespace nsNNLS 
{
	template <typename _T>
	class nnls 
	{
	public:
		nnls() 
		{
			this->x = nullptr;
			maxit = 100;
		}

		nnls(Eigen::Matrix<_T, Eigen::Dynamic, Eigen::Dynamic> *A, Eigen::Array<_T, Eigen::Dynamic, Eigen::RowMajor> *b, int maxit)
		{
			this->A = A;
			this->b = b;
			this->maxit = maxit; 
			this->x0 = nullptr;
			out.iter = -1;
			fset = 0;
			// convergence controlling parameters
			M = 100;
			beta = 1.0;
			decay = 0.9;
			pgtol = 1e-3;
			sigma = .01;
		}

		nnls(Eigen::Matrix<_T, Eigen::Dynamic, Eigen::Dynamic> *A, Eigen::Array<_T, Eigen::Dynamic, Eigen::RowMajor> *b, Eigen::Array<_T, Eigen::Dynamic, Eigen::RowMajor>* x0, int maxit)
		{
			nnls(A, b, maxit);
			this->x0 = x0; 
		}

		~nnls()
		{
			
		}

		// The various accessors and mutators (or whatever one calls 'em!)
		_T getDecay() const { return decay; }
		int getM()     const { return M; }
		_T getBeta()  const { return beta; }
		_T getObj() { return out.obj[out.iter - 1]; }
		_T getPgTol() const { return pgtol; }
		Eigen::Array<_T, Eigen::Dynamic, Eigen::RowMajor>* getSolution() { return &x; }
		size_t* getFset() { return fset; }
		size_t getMaxit() const { return maxit; }
		_T getSigma() const { return sigma; }

		void setDecay(_T d) { decay = d; }
		void setM(int m) { M = m; }
		void setBeta(_T b) { beta = b; }
		void setPgTol(_T pg) { pgtol = pg; }
		void setMaxit(size_t m) { maxit = m; }
		void setSigma(_T s) { sigma = s; }

		void setData(Eigen::Matrix<_T, Eigen::Dynamic, Eigen::Dynamic>* A, Eigen::Array<_T, Eigen::Dynamic, Eigen::RowMajor>* b) { this->A = A; this->b = b; }

		// The functions that actually launch the ship, and land it!
		int optimize()
		{
			initialize();

			int term = 0;
			_T step;

			while (!term) 
			{
				out.iter++;
				term = checkTermination();

				findFixedVariables();
				computeXandGradDelta();

				step = computeBBStep();

				x = x - (step * gradient); // x = x - step*gradient
											// project
				for (size_t i = 0; i < x.rows(); i++)
				{
					if (x[i] < 0)
					{
						x[i] = 0.0;
					}
				}

				computeObjGrad();
				// check the descent condition
				if (out.iter % M == 0) 
				{
					checkDescentUpdateBeta();
				}
			}
			fprintf(stderr, "-");
			return out.iter;
		}
		
    // The variables used during compute time
	private:                      
		Eigen::Array<_T, Eigen::Dynamic, Eigen::RowMajor> x;                  // The solution -- also current iterate
		Eigen::Array<_T, Eigen::Dynamic, Eigen::RowMajor>* x0;                 // Starting value
		Eigen::Array<_T, Eigen::Dynamic, Eigen::RowMajor> oldx;               // Previous iterate
		Eigen::Array<_T, Eigen::Dynamic, Eigen::RowMajor> gradient;                  // Current gradient
		Eigen::Array<_T, Eigen::Dynamic, Eigen::RowMajor> oldg;               // Previous gradient
		Eigen::Array<_T, Eigen::Dynamic, Eigen::RowMajor> xdelta;             // x - oldx
		Eigen::Array<_T, Eigen::Dynamic, Eigen::RowMajor> gdelta;             // g - oldg
		Eigen::Array<_T, Eigen::Dynamic, Eigen::RowMajor> refx;               // iterate from 'M' steps ago
		Eigen::Array<_T, Eigen::Dynamic, Eigen::RowMajor> refg;               // reference gradient from M steps ago

		Eigen::Matrix<_T, Eigen::Dynamic, Eigen::Dynamic> *A;
		Eigen::Array<_T, Eigen::Dynamic, Eigen::RowMajor> *b;
		Eigen::Array<_T, Eigen::Dynamic, Eigen::RowMajor> ax;                 // vector to hold A*x
		size_t* fset;               // fixed set 
		size_t fssize;              // sizeof fixed set

		// The parameters of the solver
		int maxit;               // maximum num. iters
		int   M;                    // max num. of null iterations
		_T decay;               // parameter to make diminishing scalar to decay by
		_T beta;                // diminishing scalar
		_T pgtol;               // projected gradient tolerance
		_T sigma;               // constant for descent condition

		// The solution and statistics variables
		struct out_
		{
			Eigen::Array<_T, Eigen::Dynamic, Eigen::RowMajor> obj;
			int  iter;
			Eigen::Array<_T, Eigen::Dynamic, Eigen::RowMajor> pgnorms;
			_T npg;                // inf-norm of projected gradient
		} out;

		// The helper functions used by the solver                      
		int initialize()
		{
			size_t n = A->cols();
			x.resize(n);
			gradient.resize(n);
			refx.resize(n);
			refg.resize(n);
			
			oldx.resize(n);
			oldx.setZero();
			oldg.resize(n);

			xdelta.resize(n);
			gdelta.resize(n);
			
			ax.resize(A->rows());

			fset = (size_t*) malloc(sizeof(size_t)*n);

			x.setConstant(.5);

			out.obj.resize(maxit+1);
			out.obj.setZero();
			out.pgnorms.resize(maxit+1);
			out.pgnorms.setZero();

			if (x0) 
			{
				oldg.setZero();
				x = (*x0);
				ax = (*A) * x.matrix();
				ax -= (*b);
				gradient = A->transpose() * ax.matrix();
			}
			else
			{
				// Initial gradient = A'(0 - b), since x0 = 0
				oldg = A->transpose() * b->matrix();
				oldg *= -1.0;
			}

			// old gradient = A'*(ax - b)
			ax = (*A) * x.matrix();
			ax -= (*b);
			gradient = A->transpose() * ax.matrix();

			// Set the reference iterations
			refx = x;
			refg = gradient;

			return 0;
		}

		int checkTermination()
		{
			int term = 0;
			if (out.iter >= maxit)
			{
				term = 1;
			}
			out.pgnorms[out.iter] = normProjectedGradient();
			out.npg = out.pgnorms[out.iter];
			if (out.npg < pgtol)
			{
				term = 2;
			}
			return term;
		}

		void findFixedVariables()
		{
			size_t n = A->cols();

			// Clean out fixed set first
			memset(fset, 0, sizeof(size_t)*n);
			fssize = 0;

			for (size_t i = 0; i < n; i++)
			{
				if (x[i] == 0 && gradient[i] > 0)
				{
					fset[fssize++] = i;
				}
			}
		}

		void computeXandGradDelta()
		{
			size_t n = A->cols();

			for (size_t i = 0; i < n; i++) 
			{
				xdelta[i] = x[i] - oldx[i];
				gdelta[i] = gradient[i] - oldg[i];
			}

			// Zero out fixed variables
			for (size_t i = 0; i < fssize; i++) 
			{
				xdelta[fset[i]] = 0;
				gdelta[fset[i]] = 0;
			}

			// Save oldx and oldg
			oldx = x;
			oldg = gradient;
		}

		void computeObjGrad()
		{
			ax = (*A) * x.matrix();
			ax -= (*b);        // ax = ax - b
			_T d = (ax*ax).sum();
			d = std::sqrt(d);
			out.obj[out.iter] =  (0.5 * d * d);
			gradient = A->transpose() * ax.matrix(); // A'(ax)
		}

		_T computeBBStep()
		{
			// Use xdelta and gdelta to compute BB step
			_T step = 0.0;
			_T nr = 0.0;
			_T dr = 0.0;

			if (out.iter % 2) 
			{
				for (int i = 0; i < xdelta.rows(); i++)
				{
					nr += xdelta[i] * xdelta[i];
					dr += xdelta[i] * gdelta[i];
				}
			}
			else 
			{
				for (int i = 0; i < xdelta.rows(); i++)
				{
					nr += xdelta[i] * gdelta[i];
					dr += gdelta[i] * gdelta[i];
				}
			}

			//fprintf(stderr, " nr / dr = %lf / %lf\n", nr, dr);
			if (nr == 0)
			{
				return 0;
			}
			step = nr / dr;

			step *= beta;
			return step;
		}

		void checkDescentUpdateBeta()
		{
			// if(out.iter == 0)
			//     return;

			_T d = 0;

			// compute sigma*<grad, refx - x>
			for (size_t i = 0; i < gradient.rows(); i++)
			{
				d += gradient[i] * (refx[i] - x[i]);
			}
			d *= sigma;

			if (out.iter >= M)
			{
				d = out.obj[out.iter - M] - out.obj[out.iter] - d;
			}
			else
			{
				d = out.obj[out.iter] - d;
			}

			// check if suff. descent failed
			if (d < 0) 
			{
				beta *= decay;
			} 
			else
			{                      // update reference iteration
				refx = x;
				refg = gradient;
			}
		}


		_T normProjectedGradient()
		{
			size_t n = A->cols();
			_T pg = 0.0;
			size_t ctr = 0;
			// compute the norm of the gradient for all variables not in fset
			for (size_t i = 0; i < n; i++) 
			{
				if (i==fset[ctr]) 
				{
					ctr++; continue;
				}
				pg = std::max(pg, std::abs(gradient[i]));
			}
			return pg;
		}
	};
}
