/*-----------------------------------------------------------------------------
 * Copyright (c) 2023, UChicago Argonne, LLC
 * See LICENSE file.
 *---------------------------------------------------------------------------*/

#ifndef Correlation_Coefficient_H
#define Correlation_Coefficient_H

#include <Eigen/Core>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template<typename _T>
_T find_coefficient(const Eigen::Array<_T, Eigen::Dynamic, Eigen::RowMajor>& x_arr, const Eigen::Array<_T, Eigen::Dynamic, Eigen::RowMajor>& y_arr)
{
    Eigen::Array<_T, Eigen::Dynamic, Eigen::RowMajor> xy_arr = x_arr * y_arr;
    Eigen::Array<_T, Eigen::Dynamic, Eigen::RowMajor> xx_arr = x_arr * x_arr;
    Eigen::Array<_T, Eigen::Dynamic, Eigen::RowMajor> yy_arr = y_arr * y_arr;
    _T sum_x = x_arr.sum();
    _T sum_y = y_arr.sum();
    _T sum_xy = xy_arr.sum();
    _T sqsum_x = xx_arr.sum();
    _T sqsum_y = yy_arr.sum();
    _T n = (_T)x_arr.size();
    _T corr = (n * sum_xy - sum_x * sum_y) / sqrt((n * sqsum_x - sum_x * sum_x) * (n * sqsum_y - sum_y * sum_y));
    return corr;
}

/*---------------------------------------------------------------------------*/

#endif

/*---------------------------------------------------------------------------*/
