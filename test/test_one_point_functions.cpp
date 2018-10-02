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

#include "loop_functions/one_loop_one_point_functions.hpp"
#include "numerics.h"
#include "pv.hpp"
#include "wrappers.hpp"

#include "stopwatch.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_one_point_functions

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

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

BOOST_AUTO_TEST_CASE( test_softsusy_ReA0_pos_real_args_value )
{
   const double tol = 1.e-14;
   BOOST_CHECK_CLOSE_FRACTION(softsusy::a0(std::sqrt(2), std::sqrt(1)), Re(loop_functions::A0(2, 1)), tol);
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

   BOOST_CHECK_GT(ss_time, fs_time);
}
