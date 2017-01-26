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
extern void __pole2lfunctions_MOD_ctfss(double&, const double&, const double&, const double&, const double&);
extern void __pole2lfunctions_MOD_ctfsss(double&, const double&, const double&, const double&, const double&);
extern void __pole2lfunctions_MOD_ctfssss(double&, const double&, const double&, const double&, const double&, const double&);
extern void __pole2lfunctions_MOD_ctfssff(double&, const double&, const double&, const double&, const double&, const double&);
extern void __pole2lfunctions_MOD_ctfssfbfb(double&, const double&, const double&, const double&, const double&, const double&);
extern void __pole2lfunctions_MOD_ctffffbs(double&, const double&, const double&, const double&, const double&, const double&);
extern void __pole2lfunctions_MOD_ctfffbfs(double&, const double&, const double&, const double&, const double&, const double&);
extern void __pole2lfunctions_MOD_ctffbfbfbs(double&, const double&, const double&, const double&, const double&, const double&);
extern void __pole2lfunctions_MOD_ctfsv(double&, const double&, const double&);
extern void __pole2lfunctions_MOD_ctffv(double&, const double&, const double&);
}

namespace flexiblesusy {
namespace two_loop_functions {

double TfSS(double x, double y, double z, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_ctfss(out, x, y, z, scl2);
   return out;
}

double TfSSS(double x, double y, double z, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_ctfsss(out, x, y, z, scl2);
   return out;
}

double TfSSSS(double x, double y, double z, double u, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_ctfssss(out, x, y, z, u, scl2);
   return out;
}

double TfSSFF(double x, double y, double z, double u, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_ctfssff(out, x, y, z, u, scl2);
   return out;
}

double TfSSFbFb(double x, double y, double z, double u, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_ctfssfbfb(out, x, y, z, u, scl2);
   return out;
}

double TfFFFbS(double x, double y, double z, double u, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_ctffffbs(out, x, y, z, u, scl2);
   return out;
}

double TfFFbFS(double x, double y, double z, double u, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_ctfffbfs(out, x, y, z, u, scl2);
   return out;
}

double TfFbFbFbS(double x, double y, double z, double u, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_ctffbfbfbs(out, x, y, z, u, scl2);
   return out;
}

double TfSV(double x, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_ctfsv(out, x, scl2);
   return out;
}

double TfFV(double x, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_ctffv(out, x, scl2);
   return out;
}

} // namespace two_loop_functions
} // namespace flexiblesusy
