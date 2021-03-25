/*
 * Library:   lmfit (Levenberg-Marquardt least squares fitting)
 *
 * File:      lmcurve_tyd.h
 *
 * Contents:  Declares lmcurve_tyd(), a variant of lmcurve() that weighs
 *            data points y(t) with the inverse of the standard deviations dy.
 *
 * Copyright: Joachim Wuttke, Forschungszentrum Juelich GmbH (2004-2013)
 *
 * License:   see ../COPYING (FreeBSD)
 *
 * Homepage:  apps.jcns.fz-juelich.de/lmfit
 */

#ifndef LMCURVETYD_H
#define LMCURVETYD_H
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
    const _T* dy;
    _T (*f)(const _T t, const _T* par);
} lmcurve_tyd_data_struct;


template <typename _T>
void lmcurve_tyd_evaluate(
    const _T* par, const int m_dat, const void* data, _T* fvec,
    int* info)
{
    lmcurve_tyd_data_struct* D = (lmcurve_tyd_data_struct*)data;
    int i;
    for (i = 0; i < m_dat; i++)
        fvec[i] = ( D->y[i] - D->f(D->t[i], par) ) / D->dy[i];
}


template <typename _T>
void lmcurve_tyd(
    const int n_par, _T* par, const int m_dat,
    const _T* t, const _T* y, const _T* dy,
    _T (*f)(_T t, const _T* par),
    const lm_control_struct* control, lm_status_struct* status)
{
    lmcurve_tyd_data_struct data = { t, y, dy, f };

    lmmin(n_par, par, m_dat, (const void*)&data, lmcurve_tyd_evaluate,
          control, status);
}

__END_DECLS
#endif /* LMCURVETYD_H */
