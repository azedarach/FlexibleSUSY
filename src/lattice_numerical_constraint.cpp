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

#include "lattice_numerical_constraint.hpp"

namespace flexiblesusy {

using namespace std;


void NumericalConstraint::init
(RGFlow<Lattice> *flow, size_t theory, size_t site)
{
    ForeignConstraint::init(flow, theory, site);
    x_local.resize(row.size());
    depends_on.resize(row.size());
    for (auto i: nonzeros) depends_on[i] = true;
    vector<size_t>().swap(nonzeros); // forces deallocation
}

void NumericalConstraint::operator()()
{
    x_local = x();
    Real BC0 = c(&x_local[0]);
    for (j = 0; j < row.size(); j++) {
	Real dc;
	Real xj = x_local[j];
	if (depends_on[j]) {
	    Real abserr;
	    gsl_deriv_central(&F_gsl, xj, deriv_epsilon*u(j), &dc, &abserr);
	    x_local[j] = xj;
	}
	else dc = 0;
	BC0 -= xj * (row[j] = dc);
    }
    rhs = -BC0;
    copy_row(0);
}

double NumericalConstraint::c_wrap(double xj, void *params)
{
    NumericalConstraint *self = static_cast<NumericalConstraint *>(params);
    self->x_local[self->j] = xj;
    return self->c(&self->x_local[0]);
}

void NumericalMatching::init(RGFlow<Lattice> *flow, size_t lower_theory)
{
    ForeignMatching::init(flow, lower_theory);
    w_local.resize(f->efts[TL  ].w->width);
    x_local.resize(f->efts[TL+1].w->width);
    depends_on.resize(row.size());
    for (auto i: nonzerosL) depends_on[i               ] = true;
    for (auto i: nonzerosH) depends_on[i+w_local.size()] = true;
    vector<size_t>().swap(nonzerosL); // forces deallocation
    vector<size_t>().swap(nonzerosH); // forces deallocation
}

void NumericalMatching::operator()()
{
    w_local = w(); x_local = x();
    Real MC0 = c(&w_local[0], &x_local[0]);
    for (j = 0; j < depends_on.size(); j++) {
	Real dc;
	Real wxj = wx_local(j);
	if (depends_on[j]) {
	    Real uj = j < w_local.size() ? u(0,j) : u(1,j-w_local.size());
	    Real abserr;
	    gsl_deriv_central(&F_gsl, wxj, deriv_epsilon*uj, &dc, &abserr);
	    wx_local(j) = wxj;
	}
	else dc = 0;
	MC0 -= wxj * (row[j] = dc);
    }
    rhs = -MC0;
    copy_row(0);
}

double NumericalMatching::c_wrap(double wxj, void *params)
{
    NumericalMatching *self = static_cast<NumericalMatching *>(params);
    self->wx_local(self->j) = wxj;
    return self->c(&self->w_local[0], &self->x_local[0]);
}

}
