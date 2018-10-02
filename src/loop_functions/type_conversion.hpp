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

#ifndef LOOP_FUNCTIONS_TYPE_CONVERSION_H
#define LOOP_FUNCTIONS_TYPE_CONVERSION_H

#include <complex>
#include <type_traits>

namespace flexiblesusy {

namespace loop_functions {

namespace detail {

template <template <typename, typename> class Promotion_policy, typename... Args>
struct Promote_args_helper;

template <template <typename, typename> class Promotion_policy, typename T>
struct Promote_args_helper<Promotion_policy, T> {
   using type = typename Promotion_policy<T, T>::type;
};

template <template <typename, typename> class Promotion_policy,
          typename First, typename Second, typename... Rest>
struct Promote_args_helper<Promotion_policy, First, Second, Rest...> {
   using type =
      typename Promotion_policy<
      First,
      typename Promote_args_helper<Promotion_policy, Second, Rest...>::type>::type;
};

} // namespace detail

template <typename T1, typename T2>
struct Default_promotion_rule {
   using type = std::common_type<T1, T2>;
};

template <>
struct Default_promotion_rule<int, int> {
   using type = double;
};

template <typename T>
struct Default_promotion_rule<int, T> {
   using type = Default_promotion_rule<double, T>;
};

template <typename T>
struct Default_promotion_rule<T, int> {
   using type = Default_promotion_rule<T, double>;
};

template <typename... Types>
struct Promote_args {
   using type =
      typename detail::Promote_args_helper<Default_promotion_rule, Types...>::type;
};

template <typename T>
struct Complexification {
   using type = std::complex<T>;
};

// helper function
template <typename... Types>
struct Complex_promotion {
   using type =
      typename Complexification<typename Promote_args<Types...>::type>::type;
};

} // namespace loop_functions

} // namespace flexiblesusy

#endif
