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

#include "config.h"
#include "loop_functions/one_loop_one_point_functions.hpp"
#include "numerics.h"
#include "numerics2.hpp"
#include "pv.hpp"
#include "wrappers.hpp"

#include "stopwatch.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_one_point_functions

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <array>
#include <vector>
#include <type_traits>

#ifdef ENABLE_RANDOM
#include <random>
#else
#include <cstdlib>
#endif

using namespace flexiblesusy;

double get_random_real(double lower_bound, double upper_bound)
{
#ifdef ENABLE_RANDOM
   static std::default_random_engine generator;
   static std::uniform_real_distribution<> dist(0., 1.);
   const double r = dist(generator);
#else
   const double r = 1. * std::rand() / RAND_MAX;
#endif
   return lower_bound + r * (upper_bound - lower_bound);
}

BOOST_AUTO_TEST_CASE( test_A0_return_type_real_args )
{
   int i;
   float f;
   double d;
   long double ld;

   BOOST_TEST((std::is_same<decltype(loop_functions::A0(i, i)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(f, f)),
               std::complex<float> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(d, d)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(ld, ld)),
               std::complex<long double> >::value));

   BOOST_TEST((std::is_same<decltype(loop_functions::A0(i, f)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(f, i)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(i, d)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(d, i)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(i, ld)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(ld, i)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(f, d)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(d, f)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(f, ld)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(ld, f)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(d, ld)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(ld, d)),
               std::complex<long double> >::value));

   BOOST_TEST((std::is_same<decltype(loop_functions::A0(i, i, i)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(f, f, f)),
               std::complex<float> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(d, d, d)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(ld, ld, ld)),
               std::complex<long double> >::value));

   BOOST_TEST((std::is_same<decltype(loop_functions::A0(i, i, f)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(i, f, i)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(f, i, i)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(i, i, d)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(i, d, i)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(d, i, i)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(i, i, ld)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(i, ld, i)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(ld, i, i)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(f, f, i)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(f, i, f)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(i, f, f)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(f, f, d)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(f, d, f)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(d, f, f)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(f, f, ld)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(f, ld, f)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(ld, f, f)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(d, d, i)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(d, i, d)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(i, d, d)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(d, d, f)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(d, f, d)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(f, d, d)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(d, d, ld)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(d, ld, d)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(ld, d, d)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(ld, ld, i)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(ld, i, ld)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(i, ld, ld)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(ld, ld, f)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(ld, f, ld)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(f, ld, ld)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(ld, ld, d)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(ld, d, ld)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(d, ld, ld)),
               std::complex<long double> >::value));

   BOOST_TEST((std::is_same<decltype(loop_functions::A0(i, f, d)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(i, f, ld)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(i, d, f)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(i, d, ld)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(i, ld, f)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(i, ld, d)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(f, i, d)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(f, i, ld)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(f, d, i)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(f, d, ld)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(f, ld, i)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(f, ld, d)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(d, i, f)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(d, i, ld)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(d, f, i)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(d, f, ld)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(d, ld, i)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(d, ld, f)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(ld, i, f)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(ld, i, d)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(ld, f, i)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(ld, f, d)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(ld, d, i)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(ld, d, f)),
               std::complex<long double> >::value));
}

BOOST_AUTO_TEST_CASE( test_A0_return_type_complex_args )
{
   std::complex<float> cf;
   std::complex<double> cd;
   std::complex<long double> cld;

   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cf, cf)),
               std::complex<float> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cd, cd)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cld, cld)),
               std::complex<long double> >::value));

   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cf, cd)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cf, cld)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cd, cf)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cd, cld)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cld, cf)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cld, cd)),
               std::complex<long double> >::value));

   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cf, cf, cf)),
               std::complex<float> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cd, cd, cd)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cld, cld, cld)),
               std::complex<long double> >::value));

   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cf, cf, cd)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cf, cd, cf)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cd, cf, cf)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cf, cf, cld)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cf, cld, cf)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cld, cf, cf)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cd, cd, cf)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cd, cf, cd)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cf, cd, cd)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cd, cd, cld)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cd, cld, cd)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cld, cd, cd)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cld, cld, cf)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cld, cf, cld)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cf, cld, cld)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cld, cld, cd)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cld, cd, cld)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cd, cld, cld)),
               std::complex<long double> >::value));

   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cf, cd, cld)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cf, cld, cd)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cd, cf, cld)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cd, cld, cf)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cld, cf, cd)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cld, cd, cf)),
               std::complex<long double> >::value));
}

BOOST_AUTO_TEST_CASE( test_A0_return_type_mixed_args )
{
   int i;
   float f;
   double d;
   long double ld;

   std::complex<float> cf;
   std::complex<double> cd;
   std::complex<long double> cld;

   BOOST_TEST((std::is_same<decltype(loop_functions::A0(i, cf)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(i, cd)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(i, cld)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(f, cf)),
               std::complex<float> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(f, cd)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(f, cld)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(d, cf)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(d, cd)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(d, cld)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(ld, cf)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(ld, cd)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(ld, cld)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cf, i)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cf, f)),
               std::complex<float> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cf, d)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cf, ld)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cd, i)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cd, f)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cd, d)),
               std::complex<double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cd, ld)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cld, i)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cld, f)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cld, d)),
               std::complex<long double> >::value));
   BOOST_TEST((std::is_same<decltype(loop_functions::A0(cld, ld)),
               std::complex<long double> >::value));
}

BOOST_AUTO_TEST_CASE( test_A0_zero_complex_mass_limit_real_scale )
{
   const std::size_t n_points = 100;
   const double smallest_mass_sq = 1.e-14;
   const double largest_mass_sq = 1.;

   const std::array<double, 8> angles{{-0.75 * Pi, -0.5 * Pi, -0.25 * Pi, 0.,
            0.25 * Pi, 0.5 * Pi, 0.75 * Pi, Pi }};
   std::vector<double> mass_sq(n_points);
   for (std::size_t i = 0; i < n_points - 1; ++i) {
      mass_sq[i] =
         std::exp(std::log(largest_mass_sq)
                  + i * (std::log(smallest_mass_sq) - std::log(largest_mass_sq))
                  / (n_points - 2.));
   }
   mass_sq[n_points - 1] = 0.;

   const double scale_sq = 1.;

   const auto approx_bound = [](const std::complex<double>& m2a,
                                const std::complex<double>& m2b,
                                const std::complex<double>& q2) {
      const auto delta = std::abs(m2a - m2b);
      const auto la = (std::abs(m2a) == 0. ? 0. : m2a * fast_log(m2a));
      const auto lb = (std::abs(m2b) == 0. ? 0. : m2b * fast_log(m2b));

      return 2. * (delta + std::abs(lb - la + m2a * fast_log(q2)
                                    - m2b * fast_log(q2)));
   };

   for (const auto theta : angles) {
      for (std::size_t i = 1; i < n_points; ++i) {
         const auto m_sq_im1 = std::polar(mass_sq[i - 1], theta);
         const auto m_sq_i = std::polar(mass_sq[i], theta);
         const auto val_im1 = loop_functions::A0(m_sq_im1, scale_sq);
         const auto val_i = loop_functions::A0(m_sq_i, scale_sq);

         BOOST_CHECK_LE(std::abs(val_i - val_im1),
                        approx_bound(m_sq_im1, m_sq_i, scale_sq));

         if (i == n_points - 1) {
            BOOST_CHECK_EQUAL(val_i, std::complex<double>(0., 0.));
         }
      }
   }
}

BOOST_AUTO_TEST_CASE( test_A0_negative_real_mass_value )
{
   const double tol = 1.e-15;

   BOOST_CHECK_CLOSE_FRACTION(Re(loop_functions::A0(-1, 1)), -1., tol);
   BOOST_CHECK_CLOSE_FRACTION(Im(loop_functions::A0(-1, 1)), Pi, tol);

   BOOST_CHECK_CLOSE_FRACTION(Re(loop_functions::A0(-10, 1)),
                              -10. * (1 - std::log(10.)), tol);
   BOOST_CHECK_CLOSE_FRACTION(Im(loop_functions::A0(-10, 1)),
                              10. * Pi, tol);
}

BOOST_AUTO_TEST_CASE( test_softsusy_ReA0_pos_real_args_value )
{
   const double tol = 1.e-14;
   BOOST_CHECK_CLOSE_FRACTION(softsusy::a0(std::sqrt(2), std::sqrt(1)),
                              Re(loop_functions::A0(2, 1)), tol);
}

BOOST_AUTO_TEST_CASE( test_softsusy_ReA0_pos_real_args_time )
{
   const double min_mass_sq = 1.e-10;
   const double max_mass_sq = 1.e12;
   const double min_scale_sq = 1.e-10;
   const double max_scale_sq = 1.e12;

   const std::size_t n_calls = 10000;
   std::vector<double> mass_sq(n_calls);
   std::vector<double> scale_sq(n_calls);
   for (std::size_t i = 0; i < n_calls; ++i) {
      mass_sq[i] = get_random_real(min_mass_sq, max_mass_sq);
      scale_sq[i] = get_random_real(min_scale_sq, max_scale_sq);
   }

   std::vector<double> mass(n_calls);
   std::vector<double> scale(n_calls);
   for (std::size_t i = 0; i < n_calls; ++i) {
      mass[i] = std::sqrt(mass_sq[i]);
      scale[i] = std::sqrt(scale_sq[i]);
   }

   Stopwatch stopwatch;
   stopwatch.start();
   for (std::size_t i = 0; i < n_calls; ++i) {
      volatile const auto val = softsusy::a0(mass[i], scale[i]);
   }
   stopwatch.stop();
   const double ss_time = stopwatch.get_time_in_seconds();

   stopwatch.start();
   for (std::size_t i = 0; i < n_calls; ++i) {
      volatile const auto val = Re(loop_functions::A0(mass_sq[i], scale_sq[i]));
   }
   stopwatch.stop();
   const double fs_time = stopwatch.get_time_in_seconds();

   BOOST_TEST_MESSAGE("Calculating A0 with positive arguments " << n_calls
                      << " times with\n"
                      "Softsusy     : " << ss_time << "s\n"
                      "FlexibleSUSY : " << fs_time << "s\n");

   BOOST_CHECK_LT(fs_time / ss_time, 1.2);
}

#ifdef ENABLE_LOOPTOOLS

BOOST_AUTO_TEST_CASE( test_looptools_A0_pos_real_args_value )
{
   const double tol = 1.e-14;
   BOOST_CHECK_LT(std::abs(loop_functions::A0(2, 1) -
                           passarino_veltman::A0(2, 1)), tol);
}

#endif
