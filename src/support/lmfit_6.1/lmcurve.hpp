/*
 * Library:   lmfit (Levenberg-Marquardt least squares fitting)
 *
 * File:      lmcurve.h
 *
 * Contents:  Declares lmcurve, a simplified API for curve fitting
 *            using the generic Levenberg-Marquardt routine lmmin.
 *
 * Copyright: Joachim Wuttke, Forschungszentrum Juelich GmbH (2004-2013)
 *
 * License:   see ../COPYING (FreeBSD)
 *
 * Homepage:  apps.jcns.fz-juelich.de/lmfit
 *
 * Note to programmers: Don't patch and fork, but copy and variate!
 *   If you need to compute residues differently, then please do not patch
 * lmcurve.h, but copy it to a differently named file, and change lmcurve()
 * into a differently named function declaration, like we have done in
 * lmcurve_tyd.h.
 */

#ifndef LMCURVE_H
#define LMCURVE_H
#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
#define __BEGIN_DECLS extern "C" {
#define __END_DECLS }
#else
#define __BEGIN_DECLS /* empty */
#define __END_DECLS   /* empty */
#endif

#include <lmstruct.hpp>

__BEGIN_DECLS

template <typename _T>
typedef struct {
    const _T* t;
    const _T* y;
    _T (*f)(const _T t, const _T* par);
} lmcurve_data_struct;

template <typename _T>
void lmcurve_evaluate(
    const _T* par, const int m_dat, const void* data, _T* fvec,
    int* info)
{
    lmcurve_data_struct* D = (lmcurve_data_struct*)data;
    int i;
    for (i = 0; i < m_dat; i++)
        fvec[i] = D->y[i] - D->f(D->t[i], par);
}

template <typename _T>
void lmcurve(
    const int n_par, _T* par, const int m_dat,
    const _T* t, const _T* y,
    _T (*f)(_T t, const _T* par),
    const lm_control_struct* control, lm_status_struct* status)
{
    lmcurve_data_struct data = { t, y, f };

    lmmin(n_par, par, m_dat, (const void*)&data, lmcurve_evaluate,
          control, status);
}

__END_DECLS
#endif /* LMCURVE_H */
