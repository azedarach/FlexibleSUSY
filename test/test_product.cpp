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

#include <iostream>
#include <cmath>
#include "product.hpp"
#include "wrappers.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_product

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

using namespace std;
using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE(test_product)
{
    BOOST_CHECK_EQUAL((product<size_t, 1, 10>([](size_t i){ return i; })),
		      3628800);
    BOOST_CHECK_EQUAL((product<char, -5, -1>([](char i){ return i; })),
		      -120);

    double c = 1.0/6;
    BOOST_CHECK_CLOSE_FRACTION(Pi * c * PRODUCT(i, 1, 1000, 1 - Sqr(c/i)),
			       1.0/2, 3e-5);
    BOOST_CHECK_EQUAL(PRODUCT(char, i, -5, -1, i), -120);
}

BOOST_AUTO_TEST_CASE(test_uproduct)
{
    BOOST_CHECK_EQUAL((uproduct<size_t, 1, 10>([](size_t i){ return i; })),
		      3628800);
    BOOST_CHECK_EQUAL((uproduct<char, -5, -1>([](char i){ return i; })),
		      -120);

    double c = 1.0/6;
    BOOST_CHECK_CLOSE_FRACTION(Pi * c * UPRODUCT(i, 1, 1000, 1 - Sqr(c/i)),
			       1.0/2, 3e-5);
    BOOST_CHECK_EQUAL(UPRODUCT(char, i, -5, -1, i), -120);

    // monitor evaluation order
    BOOST_CHECK_EQUAL(UPRODUCT(i,1,100,(cout << "i=" << i << "\n", pow(-1,i))),
		      1);
}
