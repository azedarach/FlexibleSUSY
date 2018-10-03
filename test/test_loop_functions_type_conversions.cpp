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

#include "loop_functions/type_conversion.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_loop_functions_type_conversions

#include <boost/test/unit_test.hpp>

#include <complex>
#include <type_traits>

#include <iostream>

using namespace flexiblesusy;
using namespace loop_functions;

BOOST_AUTO_TEST_CASE( test_complexification )
{
   BOOST_TEST((std::is_same<Complexification<int>::type, std::complex<int> >::value));
   BOOST_TEST((std::is_same<Complexification<float>::type, std::complex<float> >::value));
   BOOST_TEST((std::is_same<Complexification<double>::type, std::complex<double> >::value));
   BOOST_TEST((std::is_same<Complexification<long double>::type, std::complex<long double> >::value));
   BOOST_TEST((std::is_same<Complexification<std::complex<float> >::type, std::complex<float> >::value));
   BOOST_TEST((std::is_same<Complexification<std::complex<double> >::type, std::complex<double> >::value));
   BOOST_TEST((std::is_same<Complexification<std::complex<long double> >::type, std::complex<long double> >::value));
}

BOOST_AUTO_TEST_CASE( test_int_as_double_promotion )
{
   BOOST_TEST((std::is_same<Promote_args<int>::type, double>::value));
   BOOST_TEST((std::is_same<Promote_args<int, int>::type, double>::value));
   BOOST_TEST((std::is_same<Promote_args<int, int, int>::type, double>::value));
   BOOST_TEST((std::is_same<Promote_args<int, int, int, int>::type, double>::value));
   BOOST_TEST((std::is_same<Promote_args<int, int, int,
              int, int, int, int, int, int, int>::type, double>::value));

   BOOST_TEST((std::is_same<Promote_args<float, int>::type, double>::value));
   BOOST_TEST((std::is_same<Promote_args<int, float>::type, double>::value));

   BOOST_TEST((std::is_same<Promote_args<int, float, float>::type, double>::value));
   BOOST_TEST((std::is_same<Promote_args<float, int, float>::type, double>::value));
   BOOST_TEST((std::is_same<Promote_args<float, float, int>::type, double>::value));
   BOOST_TEST((std::is_same<Promote_args<int, int, float>::type, double>::value));
   BOOST_TEST((std::is_same<Promote_args<int, float, int>::type, double>::value));
   BOOST_TEST((std::is_same<Promote_args<float, int, int>::type, double>::value));

   BOOST_TEST((std::is_same<Promote_args<int, float, double>::type, double>::value));
   BOOST_TEST((std::is_same<Promote_args<int, double, float>::type, double>::value));
   BOOST_TEST((std::is_same<Promote_args<float, int, double>::type, double>::value));
   BOOST_TEST((std::is_same<Promote_args<float, double, int>::type, double>::value));
   BOOST_TEST((std::is_same<Promote_args<double, int, float>::type, double>::value));
   BOOST_TEST((std::is_same<Promote_args<double, float, int>::type, double>::value));
}

BOOST_AUTO_TEST_CASE( test_float_promotion )
{
   BOOST_TEST((std::is_same<Promote_args<float>::type, float>::value));
   BOOST_TEST((std::is_same<Promote_args<float, float>::type, float>::value));
   BOOST_TEST((std::is_same<Promote_args<float, float, float>::type, float>::value));
}

BOOST_AUTO_TEST_CASE( test_double_promotion )
{
   BOOST_TEST((std::is_same<Promote_args<double>::type, double>::value));
   BOOST_TEST((std::is_same<Promote_args<double, double>::type, double>::value));
   BOOST_TEST((std::is_same<Promote_args<double, double, double>::type, double>::value));

   BOOST_TEST((std::is_same<Promote_args<double, int>::type, double>::value));
   BOOST_TEST((std::is_same<Promote_args<int, double>::type, double>::value));
   BOOST_TEST((std::is_same<Promote_args<float, double>::type, double>::value));
   BOOST_TEST((std::is_same<Promote_args<double, float>::type, double>::value));

   BOOST_TEST((std::is_same<Promote_args<double, float, float>::type, double>::value));
   BOOST_TEST((std::is_same<Promote_args<float, double, float>::type, double>::value));
   BOOST_TEST((std::is_same<Promote_args<float, float, double>::type, double>::value));
   BOOST_TEST((std::is_same<Promote_args<double, double, float>::type, double>::value));
   BOOST_TEST((std::is_same<Promote_args<double, float, double>::type, double>::value));
   BOOST_TEST((std::is_same<Promote_args<float, double, double>::type, double>::value));

   BOOST_TEST((std::is_same<Promote_args<int, double, double>::type, double>::value));
   BOOST_TEST((std::is_same<Promote_args<double, int, double>::type, double>::value));
   BOOST_TEST((std::is_same<Promote_args<double, double, int>::type, double>::value));
   BOOST_TEST((std::is_same<Promote_args<int, int, double>::type, double>::value));
   BOOST_TEST((std::is_same<Promote_args<int, double, int>::type, double>::value));
   BOOST_TEST((std::is_same<Promote_args<double, int, int>::type, double>::value));

}

BOOST_AUTO_TEST_CASE( test_long_double_promotion )
{
   BOOST_TEST((std::is_same<Promote_args<long double>::type, long double>::value));
   BOOST_TEST((std::is_same<Promote_args<long double, long double>::type, long double>::value));
   BOOST_TEST((std::is_same<Promote_args<long double, long double, long double>::type, long double>::value));

   BOOST_TEST((std::is_same<Promote_args<int, long double>::type, long double>::value));
   BOOST_TEST((std::is_same<Promote_args<long double, int>::type, long double>::value));
   BOOST_TEST((std::is_same<Promote_args<float, long double>::type, long double>::value));
   BOOST_TEST((std::is_same<Promote_args<long double, float>::type, long double>::value));
   BOOST_TEST((std::is_same<Promote_args<double, long double>::type, long double>::value));
   BOOST_TEST((std::is_same<Promote_args<long double, double>::type, long double>::value));

   BOOST_TEST((std::is_same<Promote_args<long double, double, double>::type, long double>::value));
   BOOST_TEST((std::is_same<Promote_args<double, long double, double>::type, long double>::value));
   BOOST_TEST((std::is_same<Promote_args<double, double, long double>::type, long double>::value));
   BOOST_TEST((std::is_same<Promote_args<long double, long double, double>::type, long double>::value));
   BOOST_TEST((std::is_same<Promote_args<long double, double, long double>::type, long double>::value));
   BOOST_TEST((std::is_same<Promote_args<double, long double, long double>::type, long double>::value));

   BOOST_TEST((std::is_same<Promote_args<long double, float, float>::type, long double>::value));
   BOOST_TEST((std::is_same<Promote_args<float, long double, float>::type, long double>::value));
   BOOST_TEST((std::is_same<Promote_args<float, float, long double>::type, long double>::value));
   BOOST_TEST((std::is_same<Promote_args<long double, long double, float>::type, long double>::value));
   BOOST_TEST((std::is_same<Promote_args<long double, float, long double>::type, long double>::value));
   BOOST_TEST((std::is_same<Promote_args<float, long double, long double>::type, long double>::value));

   BOOST_TEST((std::is_same<Promote_args<long double, float, double>::type, long double>::value));
   BOOST_TEST((std::is_same<Promote_args<long double, double, float>::type, long double>::value));
   BOOST_TEST((std::is_same<Promote_args<float, long double, double>::type, long double>::value));
   BOOST_TEST((std::is_same<Promote_args<float, double, long double>::type, long double>::value));
   BOOST_TEST((std::is_same<Promote_args<double, long double, float>::type, long double>::value));
   BOOST_TEST((std::is_same<Promote_args<double, float, long double>::type, long double>::value));

   BOOST_TEST((std::is_same<Promote_args<int, double, long double>::type, long double>::value));
   BOOST_TEST((std::is_same<Promote_args<int, long double, double>::type, long double>::value));
   BOOST_TEST((std::is_same<Promote_args<double, int, long double>::type, long double>::value));
   BOOST_TEST((std::is_same<Promote_args<double, long double, int>::type, long double>::value));
   BOOST_TEST((std::is_same<Promote_args<long double, int, double>::type, long double>::value));
   BOOST_TEST((std::is_same<Promote_args<long double, double, int>::type, long double>::value));
}

BOOST_AUTO_TEST_CASE( test_complex_float_promotion )
{
   BOOST_TEST((std::is_same<Promote_args<std::complex<float> >::type, std::complex<float> >::value));
   BOOST_TEST((std::is_same<Promote_args<std::complex<float>, std::complex<float> >::type,
               std::complex<float> >::value));
   BOOST_TEST((std::is_same<Promote_args<std::complex<float>, std::complex<float>, std::complex<float> >::type,
               std::complex<float> >::value));

   BOOST_TEST((std::is_same<Promote_args<std::complex<float>, float>::type, std::complex<float> >::value));
   BOOST_TEST((std::is_same<Promote_args<float, std::complex<float> >::type, std::complex<float> >::value));
   BOOST_TEST((std::is_same<Promote_args<std::complex<float>, float, float>::type, std::complex<float> >::value));
   BOOST_TEST((std::is_same<Promote_args<float, std::complex<float>, float>::type, std::complex<float> >::value));
   BOOST_TEST((std::is_same<Promote_args<float, float, std::complex<float> >::type, std::complex<float> >::value));
}

BOOST_AUTO_TEST_CASE( test_int_as_double_complex_promotion )
{
   BOOST_TEST((std::is_same<Promote_args<std::complex<float>, int>::type, std::complex<double> >::value));
   BOOST_TEST((std::is_same<Promote_args<int, std::complex<float> >::type, std::complex<double> >::value));

   BOOST_TEST((std::is_same<Promote_args<int, std::complex<float>, float>::type, std::complex<double> >::value));
   BOOST_TEST((std::is_same<Promote_args<float, std::complex<float>, int>::type, std::complex<double> >::value));
}

BOOST_AUTO_TEST_CASE( test_complex_double_promotion )
{
   BOOST_TEST((std::is_same<Promote_args<std::complex<double> >::type, std::complex<double> >::value));
   BOOST_TEST((std::is_same<Promote_args<std::complex<double>, std::complex<double> >::type,
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<Promote_args<std::complex<double>, std::complex<double>, std::complex<double> >::type,
               std::complex<double> >::value));

   BOOST_TEST((std::is_same<Promote_args<std::complex<double>, std::complex<float> >::type,
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<Promote_args<std::complex<float>, std::complex<double> >::type,
               std::complex<double> >::value));

   BOOST_TEST((std::is_same<Promote_args<std::complex<double>, int>::type, std::complex<double> >::value));
   BOOST_TEST((std::is_same<Promote_args<int, std::complex<double> >::type, std::complex<double> >::value));
   BOOST_TEST((std::is_same<Promote_args<std::complex<double>, float>::type, std::complex<double> >::value));
   BOOST_TEST((std::is_same<Promote_args<float, std::complex<double> >::type, std::complex<double> >::value));
   BOOST_TEST((std::is_same<Promote_args<std::complex<double>, double>::type, std::complex<double> >::value));
   BOOST_TEST((std::is_same<Promote_args<double, std::complex<double> >::type, std::complex<double> >::value));
   BOOST_TEST((std::is_same<Promote_args<std::complex<float>, double>::type, std::complex<double> >::value));
   BOOST_TEST((std::is_same<Promote_args<double, std::complex<float> >::type, std::complex<double> >::value));

   BOOST_TEST((std::is_same<Promote_args<double, std::complex<float>, std::complex<float> >::type,
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<Promote_args<std::complex<float>, double, std::complex<float> >::type,
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<Promote_args<std::complex<float>, std::complex<float>, double>::type,
               std::complex<double> >::value));

   BOOST_TEST((std::is_same<Promote_args<std::complex<double>, float, int>::type, std::complex<double> >::value));
   BOOST_TEST((std::is_same<Promote_args<std::complex<double>, int, float>::type, std::complex<double> >::value));
   BOOST_TEST((std::is_same<Promote_args<float, std::complex<double>, int>::type, std::complex<double> >::value));
   BOOST_TEST((std::is_same<Promote_args<float, int, std::complex<double> >::type, std::complex<double> >::value));
   BOOST_TEST((std::is_same<Promote_args<int, std::complex<double>, float>::type, std::complex<double> >::value));
   BOOST_TEST((std::is_same<Promote_args<int, float, std::complex<double> >::type, std::complex<double> >::value));
}

BOOST_AUTO_TEST_CASE( test_complex_long_double_promotion )
{
   BOOST_TEST((std::is_same<Promote_args<std::complex<long double> >::type, std::complex<long double> >::value));
   BOOST_TEST((std::is_same<Promote_args<std::complex<long double>, std::complex<long double> >::type,
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<Promote_args<std::complex<long double>, std::complex<long double>,
               std::complex<long double> >::type, std::complex<long double> >::value));

   BOOST_TEST((std::is_same<Promote_args<std::complex<long double>, std::complex<float> >::type,
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<Promote_args<std::complex<float>, std::complex<long double> >::type,
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<Promote_args<std::complex<long double>, std::complex<double> >::type,
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<Promote_args<std::complex<double>, std::complex<long double> >::type,
               std::complex<long double> >::value));

   BOOST_TEST((std::is_same<Promote_args<std::complex<long double>, int>::type, std::complex<long double> >::value));
   BOOST_TEST((std::is_same<Promote_args<int, std::complex<long double> >::type, std::complex<long double> >::value));
   BOOST_TEST((std::is_same<Promote_args<std::complex<long double>, float>::type, std::complex<long double> >::value));
   BOOST_TEST((std::is_same<Promote_args<float, std::complex<long double> >::type, std::complex<long double> >::value));
   BOOST_TEST((std::is_same<Promote_args<std::complex<long double>, double>::type, std::complex<long double> >::value));
   BOOST_TEST((std::is_same<Promote_args<double, std::complex<long double> >::type, std::complex<long double> >::value));
   BOOST_TEST((std::is_same<Promote_args<std::complex<float>, long double>::type, std::complex<long double> >::value));
   BOOST_TEST((std::is_same<Promote_args<long double, std::complex<float> >::type, std::complex<long double> >::value));
   BOOST_TEST((std::is_same<Promote_args<std::complex<double>, long double>::type, std::complex<long double> >::value));
   BOOST_TEST((std::is_same<Promote_args<long double, std::complex<double> >::type, std::complex<long double> >::value));

   BOOST_TEST((std::is_same<Promote_args<long double, std::complex<float>, std::complex<float> >::type,
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<Promote_args<std::complex<float>, long double, std::complex<float> >::type,
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<Promote_args<std::complex<float>, std::complex<float>, long double>::type,
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<Promote_args<long double, std::complex<double>, std::complex<double> >::type,
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<Promote_args<std::complex<double>, long double, std::complex<double> >::type,
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<Promote_args<std::complex<double>, std::complex<double>, long double>::type,
               std::complex<long double> >::value));

   BOOST_TEST((std::is_same<Promote_args<std::complex<long double>, float, int>::type, std::complex<long double> >::value));
   BOOST_TEST((std::is_same<Promote_args<std::complex<long double>, int, float>::type, std::complex<long double> >::value));
   BOOST_TEST((std::is_same<Promote_args<float, std::complex<long double>, int>::type, std::complex<long double> >::value));
   BOOST_TEST((std::is_same<Promote_args<float, int, std::complex<long double> >::type, std::complex<long double> >::value));
   BOOST_TEST((std::is_same<Promote_args<int, std::complex<long double>, float>::type, std::complex<long double> >::value));
   BOOST_TEST((std::is_same<Promote_args<int, float, std::complex<long double> >::type, std::complex<long double> >::value));

   BOOST_TEST((std::is_same<Promote_args<std::complex<long double>, float, double>::type, std::complex<long double> >::value));
   BOOST_TEST((std::is_same<Promote_args<std::complex<long double>, double, float>::type, std::complex<long double> >::value));
   BOOST_TEST((std::is_same<Promote_args<float, std::complex<long double>, double>::type, std::complex<long double> >::value));
   BOOST_TEST((std::is_same<Promote_args<float, double, std::complex<long double> >::type, std::complex<long double> >::value));
   BOOST_TEST((std::is_same<Promote_args<double, std::complex<long double>, float>::type, std::complex<long double> >::value));
   BOOST_TEST((std::is_same<Promote_args<double, float, std::complex<long double> >::type, std::complex<long double> >::value));
}
