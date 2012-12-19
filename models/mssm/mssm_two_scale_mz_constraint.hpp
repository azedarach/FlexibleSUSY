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

#ifndef MSSM_MZ_CONSTRAINT_H
#define MSSM_MZ_CONSTRAINT_H

#include "two_scale_constraint.hpp"

class Two_scale;
template<class T> class Mssm;

/**
 * @class Mssm_mz_constraint
 * @brief MSSM low-energy constraint at the Z mass MZ
 *
 * This class represents the low-energy constraint of the MSSM at the
 * Z mass MZ.  The apply() function calculates the threshold
 * corrections to the gauge and Yukawa couplings.  It is assumed that
 * the MSSM model class is filled with the low-energy data set (see
 * MssmSoftsusy::setData).
 */

class Mssm_mz_constraint : public Constraint<Two_scale> {
public:
   Mssm_mz_constraint(Mssm<Two_scale>* mssm_,
                      double tanBeta_);
   virtual ~Mssm_mz_constraint();
   virtual void apply();
   virtual double get_scale() const;

private:
   Mssm<Two_scale>* mssm;
   double tanBeta;
   double scale;

   void update_scale();
};

#endif
