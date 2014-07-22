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

#ifndef LATTICE_FOREIGN_CONSTRAINT_H
#define LATTICE_FOREIGN_CONSTRAINT_H


#include "lattice_constraint.hpp"

namespace flexiblesusy {

class ForeignConstraint : public SingleSiteConstraint {
public:
    ForeignConstraint(size_t nrows) :
	SingleSiteConstraint(), nr(nrows) {}
    void init(RGFlow<Lattice> *flow, size_t theory, size_t site) {
	SingleSiteConstraint::init(flow, theory, site);
	row.resize(f->efts[T].w->width);
    }
    void alloc_rows() { ralloc(nr, true); }
    using SingleSiteConstraint::init;
protected:
    const std::vector<double>& x() { return f->xbuf(T, mbegin); }
    void copy_row(size_t r) {
	for (size_t j = 0; j < row.size(); j++) A(r,j) = row[j]*u(j);
	z(r) = rhs;
    }
    RVec row;
    Real rhs;
private:
    size_t nr;
};

class ForeignMatching : public InterTheoryConstraint {
public:
    ForeignMatching(size_t nrows) :
	InterTheoryConstraint(), nr(nrows) {}
    void init(RGFlow<Lattice> *flow, size_t lower_theory) {
	InterTheoryConstraint::init(flow, lower_theory);
	row.resize(f->efts[TL].w->width + f->efts[TL+1].w->width);
    }
    void alloc_rows() { ralloc(nr, true); }
    using InterTheoryConstraint::init;
protected:
    const std::vector<double>& w() { return f->xbuf(TL  , m(0)); }
    const std::vector<double>& x() { return f->xbuf(TL+1, m(1)); }
    void copy_row(size_t r) {
	RVec::const_iterator p = row.begin();
	for (size_t j = 0; j < f->efts[TL  ].w->width; j++)
	    A(r,0,j) = u(0,j) * *p++;
	for (size_t j = 0; j < f->efts[TL+1].w->width; j++)
	    A(r,1,j) = u(1,j) * *p++;
	z(r) = rhs;
    }
    RVec row;
    Real rhs;
private:
    size_t nr;
};

}

#endif // LATTICE_FOREIGN_CONSTRAINT_H
