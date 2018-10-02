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

#include "evaluation_policies.hpp"
#include "numerics2.hpp"
#include "type_conversion.hpp"

#include <cmath>
#include <complex>

#include <iostream>

namespace flexiblesusy {

namespace loop_functions {

namespace detail {

// template <typename MassSq, typename ScaleSq>
// std::enable_if<is_complex<MassSq>, MassSq>
// A0_impl(MassSq m_sq, ScaleSq scale_sq)
// {
//    using std::abs;
//    using std::arg;
//    using std::log;

//    return m_sq * (divergence + MassSq(1) - log(m_sq / scale_sq));
// }

template <typename T1, typename T2, typename T3, class Evaluation_policy>
Complex_promotion<T1,T2,T3>::type
A0_impl(T1 m_sq, T2 scale_sq, T3 divergence, const Evaluation_policy& policy)
{
   using std::log;
   const auto value = m_sq * (divergence + result_type(1) - log(m_sq / scale_sq));
   return value;
}

} // namespace detail

template <typename T1, typename T2, typename T3, class Evaluation_policy>
typename Complex_promotion<T1, T2, T3>::type
A0(T1 m_sq, T2 scale_sq, T3 divergence, const Evaluation_policy& policy)
{
   return detail::A0_impl(m_sq, scale_sq, divergence, policy);
}

template <typename T1, typename T2, typename T3>
typename Complex_promotion<T1, T2, T3>::type
A0(T1 m_sq, T2 scale_sq, T3 divergence)
{
   return detail::A0_impl(m_sq, scale_sq, divergence,
                          Default_evaluation_policy());
}

template <typename T1, typename T2>
typename Complex_promotion<T1,T2>::type
A0(T1 m_sq, T2 scale_sq)
{
   return detail::A0_impl(m_sq, scale_sq, T1(0),
                          Default_evaluation_policy());
}

} // namespace loop_functions

} // namespace flexiblesusy

#endif
