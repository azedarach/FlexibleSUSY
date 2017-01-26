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

#ifndef TWO_LOOP_FUNCTIONS
#define TWO_LOOP_FUNCTIONS

#include "cextensions.hpp"

#ifndef PVATTR
#define PVATTR noexcept ATTR(const)
#endif

namespace flexiblesusy {
namespace two_loop_functions {

// tadpole functions

double TfSS(double, double, double, double) PVATTR;
double TfSSS(double, double, double, double) PVATTR;
double TfSSSS(double, double, double, double, double) PVATTR;
double TfSSFF(double, double, double, double, double) PVATTR;
double TfSSFbFb(double, double, double, double, double) PVATTR;
double TfFFFbS(double, double, double, double, double) PVATTR;
double TfFFbFS(double, double, double, double, double) PVATTR;
double TfFbFbFbS(double, double, double, double, double) PVATTR;
double TfSV(double, double) PVATTR;
double TfFV(double, double) PVATTR;

} // namespace two_loop_functions
} // namespace flexiblesusy

#endif // TWO_LOOP_FUNCTIONS
