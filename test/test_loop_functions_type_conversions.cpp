#include "loop_functions/type_conversion.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_loop_functions_type_conversions

#include <boost/test/unit_test.hpp>

#include <complex>
#include <type_traits>

#include <iostream>

using namespace flexiblesusy;
using namespace loop_functions;

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

BOOST_AUTO_TEST_CASE( test_complexification )
{
   BOOST_TEST((std::is_same<Complexification<int>::type, std::complex<int> >::value));
   BOOST_TEST((std::is_same<Complexification<float>::type, std::complex<float> >::value));
   BOOST_TEST((std::is_same<Complexification<double>::type, std::complex<double> >::value));
   BOOST_TEST((std::is_same<Complexification<long double>::type, std::complex<long double> >::value));
}
