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

#include "array_deriv.hpp"
#include "lattice_numerical_constraint.hpp"

namespace flexiblesusy {

using namespace std;
using namespace Eigen;

const vector<size_t> NumericalConstraintCommon::empty_vector;

NumericalConstraint::NumericalConstraint
(size_t nrows, vector<size_t> dependence, Real epsilon) :
    ForeignConstraint(nrows), nonzeros(dependence), deriv_epsilon(epsilon)
{
}

void NumericalConstraint::init
(RGFlow<Lattice> *flow, size_t theory, size_t site)
{
    ForeignConstraint::init(flow, theory, site);
    depends_on.resize(rows.cols());
    if (nonzeros.empty()) fill(depends_on.begin(), depends_on.end(), true);
    else for (auto i: nonzeros) depends_on[i] = true;
    vector<size_t>().swap(nonzeros); // forces deallocation
}

void NumericalConstraint::operator()()
{
    vector<double> x_local = x();
    VectorXd BC0 = cs(&x_local[0]);
    for (size_t j = 0; j < rows.cols(); j++) {
	ArrayXd dcs(rows.rows());
	Real xj = x_local[j];
	if (depends_on[j]) {
	    ArrayXd abserr;
	    array_deriv_central(
		[&](double xj) {
		    x_local[j] = xj;
		    return cs(&x_local[0]);
		},
		xj, deriv_epsilon*u(j), dcs, abserr);
	    x_local[j] = xj;
	}
	else dcs.setZero();
	BC0 -= xj * (rows.col(j) = dcs);
    }
    rhss = -BC0;
    copy_rows();
}

AnyNumericalConstraint::AnyNumericalConstraint
(size_t nrows,
 function<ArrayXd(const AnyNumericalConstraint *, const double *x)> fxn,
 vector<size_t> dependence, double epsilon) :
    NumericalConstraint(nrows, dependence, epsilon), fxn_(fxn)
{
}

AnyNumericalConstraint::AnyNumericalConstraint
(function<double(const AnyNumericalConstraint *, const double *x)> fxn,
 vector<size_t> dependence, double epsilon) :
    AnyNumericalConstraint(
	1,
	[=](const AnyNumericalConstraint *self, const double *x) {
	    return (ArrayXd(1) << fxn(self, x)).finished();
	},
	dependence, epsilon)
{
}

NumericalMatching::NumericalMatching
(size_t nrows, vector<size_t> depL, vector<size_t> depH,
 Real epsilon) :
    ForeignMatching(nrows), nonzerosL(depL), nonzerosH(depH),
    deriv_epsilon(epsilon)
{
}

void NumericalMatching::init(RGFlow<Lattice> *flow, size_t lower_theory)
{
    ForeignMatching::init(flow, lower_theory);
    depends_on.resize(rows.cols());
    if (nonzerosL.empty())
	fill_n(depends_on.begin(), f->efts[TL].w->width, true);
    else for (auto i: nonzerosL) depends_on[i			  ] = true;
    if (nonzerosH.empty())
	fill(depends_on.begin()+f->efts[TL].w->width, depends_on.end(), true);
    else for (auto i: nonzerosH) depends_on[i+f->efts[TL].w->width] = true;
    vector<size_t>().swap(nonzerosL); // forces deallocation
    vector<size_t>().swap(nonzerosH); // forces deallocation
}

namespace {

double& cat(vector<double>& w, vector<double>& x, size_t j)
{
    return j < w.size() ? w[j] : x[j - w.size()];
}

}

void NumericalMatching::operator()()
{
    vector<double> w_local = w(), x_local = x();
    VectorXd MC0 = cs(&w_local[0], &x_local[0]);
    for (size_t j = 0; j < depends_on.size(); j++) {
	ArrayXd dcs(rows.rows());
	Real wxj = cat(w_local, x_local, j);
	if (depends_on[j]) {
	    Real uj = j < w_local.size() ? u(0,j) : u(1,j-w_local.size());
	    ArrayXd abserr;
	    array_deriv_central(
		[&](double wxj) {
		    cat(w_local, x_local, j) = wxj;
		    return cs(&w_local[0], &x_local[0]);
		},
		wxj, deriv_epsilon*uj, dcs, abserr);
	    cat(w_local, x_local, j) = wxj;
	}
	else dcs.setZero();
	MC0 -= wxj * (rows.col(j) = dcs);
    }
    rhss = -MC0;
    copy_rows();
}

AnyNumericalMatching::AnyNumericalMatching
(size_t nrows,
 function<ArrayXd(const AnyNumericalMatching *,
		  const double *w, const double *x)> fxn,
 vector<size_t> depL, vector<size_t> depH, double epsilon) :
    NumericalMatching(nrows, depL, depH, epsilon), fxn_(fxn)
{
}

AnyNumericalMatching::AnyNumericalMatching
(function<double(const AnyNumericalMatching *,
		 const double *w, const double *x)> fxn,
 vector<size_t> depL, vector<size_t> depH, double epsilon) :
    AnyNumericalMatching(
	1,
	[=](const AnyNumericalMatching *sf, const double *w, const double *x) {
	    return (ArrayXd(1) << fxn(sf, w, x)).finished();
	},
	depL, depH, epsilon)
{
}

}
