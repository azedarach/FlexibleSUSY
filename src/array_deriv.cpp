// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

/**
 * @file  src/array_deriv.cpp
 *
 * @brief The following implementation is an extension of
 * gsl_deriv_central() in GNU Scientific Library 1.16 to array-valued
 * functions
 *
 */

#include <gsl/gsl_math.h>
#include "array_deriv.hpp"

namespace flexiblesusy {

using namespace std;
using namespace Eigen;

namespace {

void
central_deriv(function<ArrayXd(double x)> f,
	      double x, double h,
	      ArrayXd& result, ArrayXd& abserr_round, ArrayXd& abserr_trunc)
{
    /* Compute the derivative using the 5-point rule (x-h, x-h/2, x,
       x+h/2, x+h). Note that the central point is not used.

       Compute the error using the difference between the 5-point and
       the 3-point rule (x-h,x,x+h). Again the central point is not
       used. */

    ArrayXd fm1 = f(x - h);
    ArrayXd fp1 = f(x + h);

    ArrayXd fmh = f(x - h / 2);
    ArrayXd fph = f(x + h / 2);

    ArrayXd r3 = 0.5 * (fp1 - fm1);
    ArrayXd r5 = (4.0 / 3.0) * (fph - fmh) - (1.0 / 3.0) * r3;

    ArrayXd e3 = (fp1.abs() + fm1.abs()) * GSL_DBL_EPSILON;
    ArrayXd e5 = (fph.abs() + fmh.abs()) * (2.0 * GSL_DBL_EPSILON) + e3;

    /* The next term is due to finite precision in x+h = O (eps * x) */

    ArrayXd dy = r3.abs().cwiseMax(r5.abs()) *
		 (fabs(x) / h / fabs(h) * GSL_DBL_EPSILON);

    /* The truncation error in the r5 approximation itself is O(h^4).
       However, for safety, we estimate the error from r5-r3, which is
       O(h^2).  By scaling h we will minimise this estimated error, not
       the actual truncation error in r5. */

    result = r5 / h;
    abserr_trunc = ((r5 - r3) / h).abs(); /* Estimated truncation error O(h^2) */
    abserr_round = (e5 / h).abs() + dy;   /* Rounding error (cancellations) */
}

} // namespace

void
array_deriv_central(function<ArrayXd(double x)> f,
		    double x, double h,
		    ArrayXd& r_0, ArrayXd& error)
{
    ArrayXd round, trunc;
    central_deriv(f, x, h, r_0, round, trunc);
    error = round + trunc;
    double round_sum = round.sum(), trunc_sum = trunc.sum();

    if (round_sum < trunc_sum && round_sum > 0 && trunc_sum > 0)
    {
	ArrayXd r_opt, round_opt, trunc_opt, error_opt;

	/* Compute an optimised stepsize to minimize the total error,
	   using the scaling of the truncation error (O(h^2)) and
	   rounding error (O(1/h)). */

	double h_opt = h * pow(round_sum / (2.0 * trunc_sum), 1.0 / 3.0);
	central_deriv(f, x, h_opt, r_opt, round_opt, trunc_opt);
	error_opt = round_opt + trunc_opt;

	/* Check that the new error is smaller, and that the new derivative
	   is consistent with the error bounds of the original estimate. */

	double error_sum = error.sum();
	if (error_opt.sum() < error_sum &&
	    (r_opt - r_0).abs().sum() < 4.0 * error_sum)
        {
	    r_0 = r_opt;
	    error = error_opt;
        }
    }
}

} // namespace flexiblesusy
