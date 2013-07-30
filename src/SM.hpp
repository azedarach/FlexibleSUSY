#ifndef SM_HPP
#define SM_HPP


#include "mathdefs.hpp"
#include "small_matrices.hpp"

namespace flexiblesusy {

const Real modVub_cent  = 4.31e-3; // PDG 2007

CM33 standard_unitary_matrix(Real s12, Real s13, Real s23, Real delta);
CM33 standard_VCKM(Real gamma, Real Vub);
CM33 standard_VCKM(Real gamma);

}

#endif // SM_HPP
