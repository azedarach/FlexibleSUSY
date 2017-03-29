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

#ifndef LRS_TENSOR_H
#define LRS_TENSOR_H

#include <Eigen/Core>
#include <complex>

namespace flexiblesusy {

/// class to represent a one component of LRS_tensor
template <typename T> struct LRS_tensor_one;
/// class to represent a zero component of LRS_tensor
template <typename T> struct LRS_tensor_zero;

template <>
template <typename Scalar, int M, int N>
struct LRS_tensor_one<Eigen::Matrix<Scalar,M,N>> {
   using Matrix_t = Eigen::Matrix<Scalar,M,N>;
   static Matrix_t get() { return Matrix_t::Ones(); }
};

template <>
template <typename Scalar>
struct LRS_tensor_one<std::complex<Scalar>> {
   using T = std::complex<Scalar>;
   static T get() { return T(1,0); }
};

template <>
struct LRS_tensor_one<double> {
   static double get() { return 1.; }
};

template <>
template <typename Scalar, int M, int N>
struct LRS_tensor_zero<Eigen::Matrix<Scalar,M,N>> {
   using Matrix_t = Eigen::Matrix<Scalar,M,N>;
   static Matrix_t get() { return Matrix_t::Zero(); }
};

template <>
template <typename Scalar>
struct LRS_tensor_zero<std::complex<Scalar>> {
   using T = std::complex<Scalar>;
   static T get() { return T{}; }
};

template <>
struct LRS_tensor_zero<double> {
   static double get() { return 0.; }
};

/**
 * @class LRS_tensor
 * @brief collection of 3 matrices
 */
template <typename T>
class LRS_tensor {
public:
   using Tuple_t = std::tuple<T,T,T>;

   LRS_tensor()
      : lrs(std::make_tuple(LRS_tensor_zero<T>::get(),LRS_tensor_zero<T>::get(),LRS_tensor_zero<T>::get())) {}
   explicit LRS_tensor(const Tuple_t& t)
      : lrs(t) {}

   const T& L() const { return std::get<0>(lrs); }
   const T& R() const { return std::get<1>(lrs); }
   const T& S() const { return std::get<2>(lrs); }
   T& L() { return std::get<0>(lrs); }
   T& R() { return std::get<1>(lrs); }
   T& S() { return std::get<2>(lrs); }
   Tuple_t tuple() const { return lrs; }
   static LRS_tensor PL() { return LRS_tensor(std::make_tuple(LRS_tensor_one<T>::get(),LRS_tensor_zero<T>::get(),LRS_tensor_zero<T>::get())); }
   static LRS_tensor PR() { return LRS_tensor(std::make_tuple(LRS_tensor_zero<T>::get(),LRS_tensor_one<T>::get(),LRS_tensor_zero<T>::get())); }
   static LRS_tensor PS() { return LRS_tensor(std::make_tuple(LRS_tensor_zero<T>::get(),LRS_tensor_zero<T>::get(),LRS_tensor_one<T>::get())); }

   LRS_tensor operator+(const LRS_tensor& rhs) const
   {
      return LRS_tensor(std::make_tuple(L()+rhs.L(), R()+rhs.R(), S()+rhs.S()));
   }

   LRS_tensor operator-(const LRS_tensor& rhs) const
   {
      return LRS_tensor(std::make_tuple(L()-rhs.L(), R()-rhs.R(), S()-rhs.S()));
   }

   LRS_tensor operator-() const
   {
      return LRS_tensor(std::make_tuple(-L(), -R(), -S()));
   }

   LRS_tensor& operator+=(const LRS_tensor& rhs)
   {
      L() += rhs.L();
      R() += rhs.R();
      S() += rhs.S();
      return *this;
   }

   template <typename U>
   friend LRS_tensor operator*(const LRS_tensor& lhs, U rhs)
   {
      return LRS_tensor(std::make_tuple(lhs.L()*rhs, lhs.R()*rhs, lhs.S()*rhs));
   }

   template <typename U>
   friend LRS_tensor operator*(U lhs, const LRS_tensor& rhs)
   {
      return LRS_tensor(std::make_tuple(lhs*rhs.L(), lhs*rhs.R(), lhs*rhs.S()));
   }

private:
   Tuple_t lrs; ///< L, R, S components
};

} // namespace flexiblesusy

#endif
