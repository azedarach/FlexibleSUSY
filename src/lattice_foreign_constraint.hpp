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
