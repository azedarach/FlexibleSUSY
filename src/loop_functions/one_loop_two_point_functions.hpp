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

#ifndef LOOP_FUNCTIONS_ONE_LOOP_ONE_POINT_FUNCTIONS_H
#define LOOP_FUNCTIONS_ONE_LOOP_ONE_POINT_FUNCTIONS_H

#include "type_conversion.hpp"
#include "wrappers.hpp"

#include <cmath>
#include <complex>

namespace flexiblesusy {

namespace loop_functions {

namespace detail {

template <typename T>
T log_branch_correction(T a, T b)
{
   const auto im_a = Im(a);
   const auto im_b = Im(b);
   const auto im_ab = Im(a * b);

   if (im_a < 0 && im_b < 0 && im_ab > 0) {
      return 1;
   } else if (im_a > 0 && im_b > 0 && im_ab < 0) {
      return -1;
   }

   return 0;
}

// @note no error checking as yet for p_sq == 0, assumed to have been
//       handled in public interface
// @todo handle types of numeric constants such as Pi
template <typename T>
typename std::enable_if<Is_complex<T>::value, typename Complexification<T>::type>::type
B0_momentum_terms_massless(T p_sq) {
   using complex_type = typename Complexification<T>::type;

   return (Re(p_sq) > 0 ? -fast_log(p_sq) + complex_type(0, Pi) + 2 :
           -fast_log(-p_sq) + 2);
}

template <typename T>
typename std::enable_if<!Is_complex<T>::value, typename Complexification<T>::type>::type
B0_momentum_terms_massless(T p_sq)
{
   using complex_type = typename Complexification<T>::type;
   using std::log;

   return (p_sq < 0 ? -log(-p_sq) + 2 : complex_type(-log(p_sq) + 2, Pi));
}

template <typename T>
typename std::enable_if<Is_complex<T>::value, typename Complexification<T>::type>::type
B0_momentum_terms_one_zero_mass(T p_sq, T m_sq)
{
}

template <typename T>
typename std::enable_if<!Is_complex<T>::value, typename Complexification<T>::type>::type
B0_momentum_terms_one_zero_mass(T p_sq, T m_sq)
{
   if (is_zero(p_sq)) {
      return 1;
   }

   const auto diff = m_sq - p_sq;
   const auto s1 = p_sq / m_sq;
   const auto s2 = diff / m_sq;
}

// @note assumes at least one argument is non-zero
template <typename T>
typename std::enable_if<!Is_complex<T>::value, typename Complexification<T>::type>::type
B0_momentum_terms(T p_sq, T m1_sq, T m2_sq)
{
   if (is_zero(m1_sq)) {
      if (is_zero(m2_sq)) {
         return B0_momentum_terms_massless(p_sq);
      } else {
         return B0_momentum_terms_one_zero_mass(p_sq, m2_sq);
      }
   } else if (is_zero(m2_sq)) {
      return B0_momentum_terms_one_zero_mass(p_sq, m1_sq);
   } else if (is_equal(m1_sq, m2_sq)) {
      return B0_momentum_terms_equal_masses(p_sq, m1_sq);
   }
   return B0_momentum_terms_distinct_masses(p_sq, m1_sq, m2_sq);
}

template <typename T>
typename std::enable_if<Is_complex<T>::value, typename Complexification<T>::type>::type
B0_momentum_terms(T p_sq, T m1_sq, T m2_sq)
{
   if (is_zero(Im(p_sq)) && is_zero(Im(m1_sq)) && is_zero(Im(m2_sq))) {
      return B0_momentum_terms(Re(p_sq), Re(m1_sq), Re(m2_sq));
   }

   if (is_zero(abs(m1_sq))) {
      if (is_zero(abs(m2_sq))) {
         return B0_momentum_terms_massless(p_sq);
      } else {
         return B0_momentum_terms_one_zero_mass(p_sq, m2_sq);
      }
   } else if (is_zero(abs(m2_sq))) {
      return B0_momentum_terms_one_zero_mass(p_sq, m1_sq);
   } else if (is_equal(m1_sq, m2_sq)) {
      return B0_momentum_terms_equal_masses(p_sq, m1_sq);
   }
   return B0_momentum_terms_distinct_masses(p_sq, m1_sq, m2_sq);
}

template <typename T>
typename std::enable_if<!Is_complex<T>::value, typename Complexification<T>::type>::type
B0_impl(T p_sq, T m1_sq, T m2_sq, T scale_sq, T divergence)
{
   using complex_type = typename Complexification<T>::type;

   const auto value = B0_momentum_terms(p_sq, m1_sq, m2_sq);

   const complex_type log_arg = 0;

   if (is_zero(m1_sq)) {
      log_arg = is_zero(m2_sq) ? 1 : m2_sq * m2_sq;
   } else if (is_zero(m2_sq)) {
      log_arg = m1_sq * m1_sq;
   } else {
      log_arg = m1_sq * m2_sq;
   }

   log_arg /= (scale_sq * scale_sq);

   const auto mass_log =
      is_zero(abs(log_arg)) ? 0 : -fast_log(log_arg) / 2;

   return divergence + value + mass_log;
}

template <typename T>
typename std::enable_if<Is_complex<T>::value, typename Complexification<T>::type>::type
B0_impl(T p_sq, T m1_sq, T m2_sq, T scale_sq, T divergence)
{
   using complex_type = typename Complexification<T>::type;
   using std::abs;

   const auto value = B0_momentum_terms(p_sq, m1_sq, m2_sq);

   T eta = 0;
   T log_arg = 0;

   if (is_zero(abs(m1_sq))) {
      log_arg = is_zero(abs(m2_sq)) ? 1 : m2_sq * m2_sq;
   } else if (is_zero(abs(m2_sq))) {
      log_arg = m1_sq * m1_sq;
   } else {
      log_arg = m1_sq * m2_sq;
      eta = log_branch_correction(m1_sq, m2_sq);
   }

   log_arg /= (scale_sq * scale_sq);
   const auto mass_log =
   is_zero(abs(log_arg)) ? 0 : -fast_log(log_arg) / 2 + eta * complex_type(0, Pi);

   return divergence + value + mass_log;
}

} // namespace detail

template <typename T1, typename T2, typename T3, typename T4, typename T5>
typename Promote_and_complexify<T1, T2, T3, T4, T5>::type
B0(T1 p_sq, T2 m1_sq, T3 m2_sq, T4 scale_sq, T5 divergence)
{
   using result_type = typename Promote_and_complexify<T1, T2, T3, T4, T5>::type;
   using argument_type = typename Promote_args<T1, T2, T3, T4, T5>::type;
   using std::abs;

   if (is_zero(abs(p_sq)) && is_zero(abs(m1_sq)) && is_zero(abs(m2_sq))) {
      return result_type(divergence);
   }

   return detail::B0_impl<argument_type>(p_sq, m1_sq, m2_sq, scale_sq, divergence);
}

template <typename T1, typename T2, typename T3, typename T4>
typename Promote_and_complexify<T1, T2, T3, T4>::type
B0(T1 p_sq, T2 m1_sq, T3 m2_sq, T4 scale_sq)
{
   return B0(p_sq, m1_sq, m2_sq, scale_sq, T1(0));
}

} // namespace loop_functions

} // namespace flexiblesusy

#endif
