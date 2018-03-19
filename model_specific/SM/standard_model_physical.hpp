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

// File generated at Sat 13 Feb 2016 15:47:12

#ifndef STANDARD_MODEL_PHYSICAL_H
#define STANDARD_MODEL_PHYSICAL_H

#include "linalg2.hpp"
#include <Eigen/Core>

#include <iosfwd>
#include <string>

namespace flexiblesusy {

namespace standard_model {

struct Standard_model_physical {
   Standard_model_physical();
   void clear();
   void convert_to_hk();   ///< converts pole masses to HK convention
   void convert_to_slha(); ///< converts pole masses to SLHA convention
   Eigen::ArrayXd get() const; ///< returns array with all masses and mixings
   void set(const Eigen::ArrayXd&); ///< set all masses and mixings
   Eigen::ArrayXd get_masses() const; ///< returns array with all masses
   void set_masses(const Eigen::ArrayXd&); ///< set all masses
   void print(std::ostream&) const;

   double MVG;
   double MHp;
   Eigen::Array<double,3,1> MFv;
   double MAh;
   double Mhh;
   double MVP;
   double MVZ;
   Eigen::Array<double,3,1> MFd;
   Eigen::Array<double,3,1> MFu;
   Eigen::Array<double,3,1> MFe;
   double MVWp;
   Eigen::Array<double,2,1> MVPVZ;

   Eigen::Matrix<std::complex<double>,3,3> Vd;
   Eigen::Matrix<std::complex<double>,3,3> Ud;
   Eigen::Matrix<std::complex<double>,3,3> Vu;
   Eigen::Matrix<std::complex<double>,3,3> Uu;
   Eigen::Matrix<std::complex<double>,3,3> Ve;
   Eigen::Matrix<std::complex<double>,3,3> Ue;
   Eigen::Matrix<double,2,2> ZZ;

   /// @todo adapt to squared masses
   double get_MVG() const { return MVG; }
   double get_MHp() const { return MHp; }
   Eigen::Array<double,3,1> get_MFv() const { return MFv; }
   double get_MFv(int i) const { return MFv(i); }
   double get_MAh() const { return MAh; }
   double get_Mhh() const { return Mhh; }
   Eigen::Array<double,3,1> get_MFd() const { return MFd; }
   double get_MFd(int i) const { return MFd(i); }
   Eigen::Array<double,3,1> get_MFu() const { return MFu; }
   double get_MFu(int i) const { return MFu(i); }
   Eigen::Array<double,3,1> get_MFe() const { return MFe; }
   double get_MFe(int i) const { return MFe(i); }
   double get_MVWp() const { return MVWp; }
   double get_MVP() const { return MVP; }
   double get_MVZ() const { return MVZ; }
};

std::ostream& operator<<(std::ostream&, const Standard_model_physical&);

} // namespace standard_model

} // namespace flexiblesusy

#endif
