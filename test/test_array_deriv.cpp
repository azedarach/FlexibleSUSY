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

#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include "array_deriv.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_array_deriv

#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/mpl/list.hpp>

using namespace std;
using namespace Eigen;
using namespace flexiblesusy;

ArrayXd foo(double x)
{
    return (ArrayXd(2) << 7 / sqrt(x), log(x)).finished();
}

ArrayXd bar(double x)
{
    return (ArrayXd(4) << sin(x), exp(x), atan(x), 11 / x).finished();
}

struct Param_gsl {
    function<ArrayXd(double x)> f;
    size_t i;
};

double f_wrap(double x, void *params)
{
    Param_gsl *p = (Param_gsl *)params;
    return p->f(x)[p->i];
}

void gsl_array_deriv_central(function<ArrayXd(double x)> f,
			     double x, double h,
			     ArrayXd& result, ArrayXd& abserr)
{
    Param_gsl p;
    p.f = f;

    gsl_function f_gsl;
    f_gsl.function = f_wrap;
    f_gsl.params = &p;

    size_t n = f(x).size();
    result.resize(n);
    abserr.resize(n);

    for (size_t i = 0; i < n; i++) {
	p.i = i;
	gsl_deriv_central(&f_gsl, x, h, &result[i], &abserr[i]);
    }
}

template<ArrayXd f_(double), int sign_>
struct Test_deriv {
    enum { sign = sign_ };
    ArrayXd (*fp())(double) { return f_; };
};

typedef boost::mpl::list<
    Test_deriv<foo,  1>,
    Test_deriv<bar,  1>,
    Test_deriv<bar, -1>
> deriv_tests;

BOOST_AUTO_TEST_CASE_TEMPLATE(test_array_deriv_central, T, deriv_tests)
{
    int sign = T::sign;

    for (size_t n = 10000; n; n--) {
	double x = sign * (1e2 * rand() / RAND_MAX + 1);

	ArrayXd derivs_arr, abserrs_arr;
	array_deriv_central
	(T().fp(), x, GSL_SQRT_DBL_EPSILON, derivs_arr, abserrs_arr);

	ArrayXd derivs_gsl, abserrs_gsl;
	gsl_array_deriv_central
	(T().fp(), x, GSL_SQRT_DBL_EPSILON, derivs_gsl, abserrs_gsl);

	BOOST_CHECK(((derivs_arr - derivs_gsl).abs() < abserrs_arr).all());
	BOOST_WARN (((derivs_arr - derivs_gsl).abs() < abserrs_gsl).all());
    }
}
