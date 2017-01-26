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

#include "two_loop_functions.hpp"

extern "C" {
// tadpole functions
extern double __pole2lfunctions_MOD_tfss(const double&, const double&, const double&, const double&);
extern double __pole2lfunctions_MOD_tfsss(const double&, const double&, const double&, const double&);
extern double __pole2lfunctions_MOD_tfssss(const double&, const double&, const double&, const double&, const double&);
extern double __pole2lfunctions_MOD_tfssff(const double&, const double&, const double&, const double&, const double&);
extern double __pole2lfunctions_MOD_tfssfbfb(const double&, const double&, const double&, const double&, const double&);
extern double __pole2lfunctions_MOD_tffffbs(const double&, const double&, const double&, const double&, const double&);
extern double __pole2lfunctions_MOD_tfffbfs(const double&, const double&, const double&, const double&, const double&);
extern double __pole2lfunctions_MOD_tffbfbfbs(const double&, const double&, const double&, const double&, const double&);
extern double __pole2lfunctions_MOD_tfsv(const double&, const double&);
extern double __pole2lfunctions_MOD_tffv(const double&, const double&);
}

namespace flexiblesusy {
namespace two_loop_functions {

double TfSS(double x, double y, double z, double scl2) noexcept
{
   return __pole2lfunctions_MOD_tfss(x, y, z, scl2);
}

double TfSSS(double x, double y, double z, double scl2) noexcept
{
   return __pole2lfunctions_MOD_tfsss(x, y, z, scl2);
}

double TfSSSS(double x, double y, double z, double u, double scl2) noexcept
{
   return __pole2lfunctions_MOD_tfssss(x, y, z, u, scl2);
}

double TfSSFF(double x, double y, double z, double u, double scl2) noexcept
{
   return __pole2lfunctions_MOD_tfssff(x, y, z, u, scl2);
}

double TfSSFbFb(double x, double y, double z, double u, double scl2) noexcept
{
   return __pole2lfunctions_MOD_tfssfbfb(x, y, z, u, scl2);
}

double TfFFFbS(double x, double y, double z, double u, double scl2) noexcept
{
   return __pole2lfunctions_MOD_tffffbs(x, y, z, u, scl2);
}

double TfFFbFS(double x, double y, double z, double u, double scl2) noexcept
{
   return __pole2lfunctions_MOD_tfffbfs(x, y, z, u, scl2);
}

double TfFbFbFbS(double x, double y, double z, double u, double scl2) noexcept
{
   return __pole2lfunctions_MOD_tffbfbfbs(x, y, z, u, scl2);
}

double TfSV(double x, double scl2) noexcept
{
   return __pole2lfunctions_MOD_tfsv(x, scl2);
}

double TfFV(double x, double scl2) noexcept
{
   return __pole2lfunctions_MOD_tffv(x, scl2);
}

} // namespace two_loop_functions
} // namespace flexiblesusy
