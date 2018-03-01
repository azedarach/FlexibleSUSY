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
#include <limits>

#include "numerics.h"
#include "numerics2.hpp"
#include "pv2.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_pv2

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

struct Values {
   Values(double p_, double m1_, double m2_, double q_)
      : p(p_), m1(m1_), m2(m2_), q(q_) {}
   double p{}, m1{}, m2{}, q{};
};

constexpr double sqr(double a) noexcept { return a*a; }

const double scale  = 100;
const double scale2 = sqr(scale);
const double p  = 91.0;
const double p2 = sqr(p);

BOOST_AUTO_TEST_CASE( test_ReA0 )
{
   BOOST_CHECK_EQUAL(flexiblesusy::a0(0., scale2), 0.);
   BOOST_CHECK_CLOSE_FRACTION(softsusy::a0(p, scale), flexiblesusy::a0(p2, scale2), 1e-15);
   BOOST_CHECK_CLOSE_FRACTION(softsusy::a0(0, scale), flexiblesusy::a0(0., scale2), 1e-15);
   BOOST_CHECK_CLOSE_FRACTION(softsusy::a0(1e-4, scale), flexiblesusy::a0(1e-8, scale2), 1e-15);
   BOOST_CHECK_CLOSE_FRACTION(softsusy::a0(1e-5, scale), flexiblesusy::a0(1e-10, scale2), 1e-15);
}

BOOST_AUTO_TEST_CASE( test_ReB0_values )
{
   const std::vector<Values> vals = {
      Values(0.  , 1.  , 0., 1.),
      Values(1e-1, 1.  , 0., 1.),
      Values(1e-2, 1.  , 0., 1.),
      Values(1e-3, 1.  , 0., 1.),
      Values(1e-4, 1.  , 0., 1.),
      Values(1e-5, 1.  , 0., 1.),
      Values(0.  , 1.  , 1e-15, 1.),
      Values(1e-1, 1.  , 1e-15, 1.),
      Values(1e-2, 1.  , 1e-15, 1.),
      Values(1e-3, 1.  , 1e-15, 1.),
      Values(1e-4, 1.  , 1e-15, 1.),
      Values(1e-5, 1.  , 1e-15, 1.),
      Values(0.  , 1e20, 0., 1.),
      Values(1e-1, 1e20, 0., 1.),
      Values(1e-2, 1e20, 0., 1.),
      Values(1e-3, 1e20, 0., 1.),
      Values(1e-4, 1e20, 0., 1.),
      Values(1e-5, 1e20, 0., 1.),
      Values(0.  , 0.  , 1., 1.),
      Values(1e-1, 0.  , 1., 1.),
      Values(1e-2, 0.  , 1., 1.),
      Values(1e-3, 0.  , 1., 1.),
      Values(1e-4, 0.  , 1., 1.),
      Values(1e-5, 0.  , 1., 1.),
      Values(0.  , 1e20, 1., 1.),
      Values(1e-1, 1e20, 1., 1.),
      Values(1e-2, 1e20, 1., 1.),
      Values(1e-3, 1e20, 1., 1.),
      Values(1e-4, 1e20, 1., 1.),
      Values(1e-5, 1e20, 1., 1.),
      Values(0.  , 0., 1e20, 1.),
      Values(1e-1, 0., 1e20, 1.),
      Values(1e-2, 0., 1e20, 1.),
      Values(1e-3, 0., 1e20, 1.),
      Values(1e-4, 0., 1e20, 1.),
      Values(1e-5, 0., 1e20, 1.),
      Values(1.  , 1.  , 1., 1.),
      Values(1.  , 2.  , 3., 4.)
   };

   for (const auto v: vals) {
      const auto p = v.p;
      const auto m1 = v.m1;
      const auto m2 = v.m2;
      const auto q = v.q;

      const auto v1 = softsusy::b0(p,m1,m2,q);
      const auto v2 = flexiblesusy::b0(sqr(p),sqr(m1),sqr(m2),sqr(q));

      BOOST_TEST_MESSAGE("testing p = " << p << ", m1 = " << m1
                         << ", m2 = " << m2 << ", q = " << q);
      BOOST_CHECK_CLOSE_FRACTION(v1, v2, 8e-5);
   }
}

BOOST_AUTO_TEST_CASE( test_ReB1_values )
{
   const std::vector<Values> vals = {
      Values(0.  , 1.  , 0., 1.),
      Values(1e-1, 1.  , 0., 1.),
      Values(1e-2, 1.  , 0., 1.),
      Values(1e-3, 1.  , 0., 1.),
      Values(1e-4, 1.  , 0., 1.),
      Values(1e-5, 1.  , 0., 1.),
      Values(0.  , 1.  , 1e-15, 1.),
      Values(1e-1, 1.  , 1e-15, 1.),
      Values(1e-2, 1.  , 1e-15, 1.),
      Values(1e-3, 1.  , 1e-15, 1.),
      Values(1e-4, 1.  , 1e-15, 1.),
      Values(1e-5, 1.  , 1e-15, 1.),
      Values(0.  , 1e20, 0., 1.),
      Values(1e-1, 1e20, 0., 1.),
      Values(1e-2, 1e20, 0., 1.),
      Values(1e-3, 1e20, 0., 1.),
      Values(1e-4, 1e20, 0., 1.),
      Values(1e-5, 1e20, 0., 1.),
      Values(0.  , 0.  , 1., 1.),
      Values(1e-1, 0.  , 1., 1.),
      Values(1e-2, 0.  , 1., 1.),
      Values(1e-3, 0.  , 1., 1.),
      Values(1e-4, 0.  , 1., 1.),
      Values(1e-5, 0.  , 1., 1.),
      Values(0.  , 1e20, 1., 1.),
      Values(1e-1, 1e20, 1., 1.),
      Values(1e-2, 1e20, 1., 1.),
      Values(1e-3, 1e20, 1., 1.),
      Values(1e-4, 1e20, 1., 1.),
      Values(1e-5, 1e20, 1., 1.),
      Values(0.  , 0., 1e20, 1.),
      Values(1e-1, 0., 1e20, 1.),
      Values(1e-2, 0., 1e20, 1.),
      Values(1e-3, 0., 1e20, 1.),
      Values(1e-4, 0., 1e20, 1.),
      Values(1e-5, 0., 1e20, 1.),
      Values(1.  , 1.  , 1., 1.)
   };

   for (const auto v: vals) {
      const auto p = v.p;
      const auto m1 = v.m1;
      const auto m2 = v.m2;
      const auto q = v.q;

      const auto v1 = softsusy::b1(p,m1,m2,q);
      const auto v2 = flexiblesusy::b1(sqr(p),sqr(m1),sqr(m2),sqr(q));

      BOOST_TEST_MESSAGE("testing p = " << p << ", m1 = " << m1
                         << ", m2 = " << m2 << ", q = " << q);
      BOOST_CHECK_CLOSE_FRACTION(v1, v2, 1e-12);
   }
}
