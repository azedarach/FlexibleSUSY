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

#ifndef ARRAY_DERIV_H
#define ARRAY_DERIV_H

#include <functional>
#include <Eigen/Dense>

namespace flexiblesusy {

/**
 * Takes the numerical derivative of an array-valued function using an
 * adaptive central difference algorithm.  The error in each component
 * of the derivative is in general larger than what would result from
 * applying gsl_deriv_central() on the same component of f.
 *
 * @param f	 function to differentiate
 * @param x	 argument to derivative of f
 * @param h	 initial step size
 * @param result derivative of f
 * @param abserr array of error estimates
 *
 */
void array_deriv_central(std::function<Eigen::ArrayXd(double x)> f,
			 double x, double h,
			 Eigen::ArrayXd& result, Eigen::ArrayXd& abserr);

}

#endif
