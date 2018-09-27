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

#include "numerics2.hpp"

#include <cmath>
#include <complex>

namespace flexiblesusy {

namespace loop_functions {

namespace detail {

template <typename MassSq, typename ScaleSq>
std::enable_if<is_complex<MassSq>, MassSq>
A0_impl(MassSq m_sq, ScaleSq scale_sq)
{
   using std::abs;
   using std::arg;
   using std::log;

   return m_sq * (divergence + MassSq(1) - log(m_sq / scale_sq));
}

template <typename Real>
Real A0_impl(Real m_sq, Real scale_sq, Real divergence)
{
   using std::log;
   return m_sq * (divergence + Real(1) - log(m_sq / scale_sq));
}

} // namespace detail

template <typename T1, typename T2>
Real A0(T1 m_sq, T2 scale_sq)
{
   return detail::A0_impl(m_sq, scale_sq, T1(0));
}

} // namespace loop_functions

} // namespace flexiblesusy
