#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CMSSM_self_energies

#include <boost/test/unit_test.hpp>
#include "wrappers.hpp"
#include "CMSSM_mass_eigenstates.hpp"

using namespace flexiblesusy;

CMSSM_mass_eigenstates setup()
{
   CMSSM_input_parameters input;
   input.m0 = 125.;
   input.m12 = 500.;
   input.TanBeta = 10.;
   input.SignMu = 1;
   input.Azero = 100.;

   const double ALPHASMZ = 0.1176;
   const double ALPHAMZ = 1.0 / 127.918;
   const double sinthWsq = 0.23122;
   const double alpha1 = 5 * ALPHAMZ / (3 * (1 - sinthWsq));
   const double alpha2 = ALPHAMZ / sinthWsq;
   const double g1 = sqrt(4 * Pi * alpha1);
   const double g2 = sqrt(4 * Pi * alpha2);
   const double g3 = sqrt(4 * Pi * ALPHASMZ);
   const double tanBeta = input.TanBeta;
   const double sinBeta = sin(atan(tanBeta));
   const double cosBeta = cos(atan(tanBeta));
   const double M12 = input.m12;
   const double m0 = input.m0;
   const double a0 = input.Azero;
   const double root2 = sqrt(2.0);
   const double vev = 246.0;
   const double vu = vev * sinBeta;
   const double vd = vev * cosBeta;
   const double susyMu = input.SignMu * 120.0;
   const double BMu = Sqr(2.0 * susyMu);
   const double scale = 91.;
   Eigen::Matrix<double,3,3> Yu(Eigen::Matrix<double,3,3>::Zero()),
      Yd(Eigen::Matrix<double,3,3>::Zero()),
      Ye(Eigen::Matrix<double,3,3>::Zero()),
      mm0(Eigen::Matrix<double,3,3>::Zero());
   Yu(2,2) = 165.0   * root2 / (vev * sinBeta);
   Yd(2,2) = 2.9     * root2 / (vev * cosBeta);
   Ye(2,2) = 1.77699 * root2 / (vev * cosBeta);
   mm0 = Sqr(m0) * Eigen::Matrix<double,3,3>::Identity();

   CMSSM_mass_eigenstates m;

   m.set_input_parameters(input);
   m.set_scale(1000.);
   m.set_loops(2);
   m.set_thresholds(2);
   m.set_g1(g1);
   m.set_g2(g2);
   m.set_g3(g3);
   m.set_Yu(Yu);
   m.set_Yd(Yd);
   m.set_Ye(Ye);
   m.set_MassB(M12);
   m.set_MassG(M12);
   m.set_MassWB(M12);
   m.set_mq2(mm0);
   m.set_ml2(mm0);
   m.set_md2(mm0);
   m.set_mu2(mm0);
   m.set_me2(mm0);
   m.set_mHd2(Sqr(m0));
   m.set_mHu2(Sqr(m0));
   m.set_TYu(a0 * Yu);
   m.set_TYd(a0 * Yd);
   m.set_TYe(a0 * Ye);
   m.set_Mu(susyMu);
   m.set_BMu(BMu);
   m.set_vu(vu);
   m.set_vd(vd);

   m.calculate_DRbar_masses();

   return m;
}

#define CHECK_SCALAR_SE(particle,dim)                                   \
   do {                                                                 \
      const double p = 100.;                                            \
      const auto se = m.self_energy_ ## particle ## _1loop(p);          \
                                                                        \
      BOOST_TEST_MESSAGE("comparing self-energies for " # particle);    \
                                                                        \
      for (int i = 0; i < dim; i++) {                                   \
         for (int k = 0; k < dim; k++) {                                \
            BOOST_CHECK_CLOSE_FRACTION(Re(m.self_energy_ ## particle ## _1loop(p,i,k)), Re(se(i,k)), 1e-10); \
            BOOST_CHECK_CLOSE_FRACTION(Im(m.self_energy_ ## particle ## _1loop(p,i,k)), Im(se(i,k)), 1e-10); \
         }                                                              \
      }                                                                 \
   } while (false)

#define CHECK_FERMION_SE(particle,dim)                                  \
   do {                                                                 \
      const double p = 100.;                                            \
      const auto seS = m.self_energy_ ## particle ## _1loop_1(p);       \
      const auto seL = m.self_energy_ ## particle ## _1loop_PL(p);      \
      const auto seR = m.self_energy_ ## particle ## _1loop_PR(p);      \
      const auto se  = m.self_energy_ ## particle ## _1loop(p);         \
                                                                        \
      BOOST_TEST_MESSAGE("comparing self-energies for " # particle);    \
                                                                        \
      for (int i = 0; i < dim; i++) {                                   \
         for (int k = 0; k < dim; k++) {                                \
            BOOST_CHECK_CLOSE_FRACTION(Re(m.self_energy_ ## particle ## _1loop_1 (p,i,k)), Re(seS(i,k))   , 1e-10); \
            BOOST_CHECK_CLOSE_FRACTION(Im(m.self_energy_ ## particle ## _1loop_1 (p,i,k)), Im(seS(i,k))   , 1e-10); \
            BOOST_CHECK_CLOSE_FRACTION(Re(m.self_energy_ ## particle ## _1loop_1 (p,i,k)), Re(se.S()(i,k)), 1e-10); \
            BOOST_CHECK_CLOSE_FRACTION(Im(m.self_energy_ ## particle ## _1loop_1 (p,i,k)), Im(se.S()(i,k)), 1e-10); \
            BOOST_CHECK_CLOSE_FRACTION(Re(m.self_energy_ ## particle ## _1loop   (p,i,k).S()), Re(se.S()(i,k)), 1e-10); \
            BOOST_CHECK_CLOSE_FRACTION(Im(m.self_energy_ ## particle ## _1loop   (p,i,k).S()), Im(se.S()(i,k)), 1e-10); \
            BOOST_CHECK_CLOSE_FRACTION(Re(m.self_energy_ ## particle ## _1loop_PL(p,i,k)), Re(seL(i,k))   , 1e-10); \
            BOOST_CHECK_CLOSE_FRACTION(Im(m.self_energy_ ## particle ## _1loop_PL(p,i,k)), Im(seL(i,k))   , 1e-10); \
            BOOST_CHECK_CLOSE_FRACTION(Re(m.self_energy_ ## particle ## _1loop_PL(p,i,k)), Re(se.L()(i,k)), 1e-10); \
            BOOST_CHECK_CLOSE_FRACTION(Im(m.self_energy_ ## particle ## _1loop_PL(p,i,k)), Im(se.L()(i,k)), 1e-10); \
            BOOST_CHECK_CLOSE_FRACTION(Re(m.self_energy_ ## particle ## _1loop   (p,i,k).L()), Re(se.L()(i,k)), 1e-10); \
            BOOST_CHECK_CLOSE_FRACTION(Im(m.self_energy_ ## particle ## _1loop   (p,i,k).L()), Im(se.L()(i,k)), 1e-10); \
            BOOST_CHECK_CLOSE_FRACTION(Re(m.self_energy_ ## particle ## _1loop_PR(p,i,k)), Re(seR(i,k))   , 1e-10); \
            BOOST_CHECK_CLOSE_FRACTION(Im(m.self_energy_ ## particle ## _1loop_PR(p,i,k)), Im(seR(i,k))   , 1e-10); \
            BOOST_CHECK_CLOSE_FRACTION(Re(m.self_energy_ ## particle ## _1loop_PR(p,i,k)), Re(se.R()(i,k)), 1e-10); \
            BOOST_CHECK_CLOSE_FRACTION(Im(m.self_energy_ ## particle ## _1loop_PR(p,i,k)), Im(se.R()(i,k)), 1e-10); \
            BOOST_CHECK_CLOSE_FRACTION(Re(m.self_energy_ ## particle ## _1loop   (p,i,k).R()), Re(se.R()(i,k)), 1e-10); \
            BOOST_CHECK_CLOSE_FRACTION(Im(m.self_energy_ ## particle ## _1loop   (p,i,k).R()), Im(se.R()(i,k)), 1e-10); \
         }                                                              \
      }                                                                 \
   } while (false)

BOOST_AUTO_TEST_CASE(test_CMSSM_self_energies)
{
   const auto m = setup();

   CHECK_SCALAR_SE(hh,2);
   CHECK_SCALAR_SE(Ah,2);
   CHECK_SCALAR_SE(Hpm,2);
   CHECK_SCALAR_SE(Su,6);
   CHECK_SCALAR_SE(Sd,6);
   CHECK_SCALAR_SE(Se,6);
   CHECK_SCALAR_SE(Sv,3);

   CHECK_FERMION_SE(Fe,3);
   CHECK_FERMION_SE(Fu,3);
   CHECK_FERMION_SE(Fd,3);
   CHECK_FERMION_SE(Fv,3);
   CHECK_FERMION_SE(Chi,4);
   CHECK_FERMION_SE(Cha,2);
}
