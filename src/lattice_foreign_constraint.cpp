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

#include "lattice_foreign_constraint.hpp"

namespace flexiblesusy {

ForeignConstraint::ForeignConstraint(size_t nrows) :
    SingleSiteConstraint(), nr(nrows)
{
}

void ForeignConstraint::init(RGFlow<Lattice> *flow, size_t theory, size_t site)
{
    SingleSiteConstraint::init(flow, theory, site);
    rows.resize(nr, f->efts[T].w->width);
    rhss.resize(nr);
}

void ForeignConstraint::copy_rows()
{
    for (size_t j = 0; j < rows.cols(); j++)
	for (size_t r = 0; r < nr; r++)
	    A(r,j) = rows(r,j)*u(j);

    for (size_t r = 0; r < nr; r++)
	z(r) = rhss(r);
}

ForeignMatching::ForeignMatching(size_t nrows) :
    InterTheoryConstraint(), nr(nrows)
{
}

void ForeignMatching::init(RGFlow<Lattice> *flow, size_t lower_theory)
{
    InterTheoryConstraint::init(flow, lower_theory);
    rows.resize(nr, f->efts[TL].w->width + f->efts[TL+1].w->width);
    rhss.resize(nr);
}

void ForeignMatching::copy_rows()
{
    for (size_t j = 0; j < f->efts[TL  ].w->width; j++)
	for (size_t r = 0; r < nr; r++)
	    A(r,0,j) = u(0,j) * rows(r,j);
    for (size_t j = 0; j < f->efts[TL+1].w->width; j++)
	for (size_t r = 0; r < nr; r++)
	    A(r,1,j) = u(1,j) * rows(r,j+f->efts[TL].w->width);

    for (size_t r = 0; r < nr; r++)
	z(r) = rhss(r);
}

} // namespace flexiblesusy
