#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CMSSM_database

#include <boost/test/unit_test.hpp>

#include "CMSSM_mass_eigenstates.hpp"
#include "CMSSM_input_parameters.hpp"
#include "CMSSM_utilities.hpp"
#include "ew_input.hpp"
#include "test.h"
#include "lowe.h"

using namespace flexiblesusy;

void setup_CMSSM(CMSSM_mass_eigenstates& m, const CMSSM_input_parameters& input)
{
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
   const double scale = Electroweak_constants::MZ;

   Eigen::Matrix<double,3,3> Yu(Eigen::Matrix<double,3,3>::Zero()),
      Yd(Eigen::Matrix<double,3,3>::Zero()),
      Ye(Eigen::Matrix<double,3,3>::Zero()),
      mm0(Eigen::Matrix<double,3,3>::Zero());
   Yu(2,2) = 165.0   * root2 / (vev * sinBeta);
   Yd(2,2) = 2.9     * root2 / (vev * cosBeta);
   Ye(2,2) = 1.77699 * root2 / (vev * cosBeta);
   mm0 = Sqr(m0) * Eigen::Matrix<double,3,3>::Identity();

   m.set_input_parameters(input);
   m.set_scale(scale);
   m.set_loops(2);
   m.set_thresholds(3);
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

   m.do_force_output(true);

   m.calculate_DRbar_masses();
   m.solve_ewsb();
   m.calculate_spectrum();
}

template <class T, class P>
void test_parameter_equality(const T& a, const P& b)
{
   TEST_EQUALITY(a.get_loops(), b.get_loops());
   TEST_EQUALITY(a.get_scale(), b.get_scale());
   TEST_EQUALITY(a.get_thresholds(), b.get_thresholds());
   TEST_EQUALITY(a.get_ewsb_loop_order(), b.get_ewsb_loop_order());
   TEST_EQUALITY(a.get_pole_mass_loop_order(), b.get_pole_mass_loop_order());
   TEST_EQUALITY(a.get_precision(), b.get_precision());
   TEST_EQUALITY(a.get_ewsb_iteration_precision(), b.get_ewsb_iteration_precision());

   TEST_EQUALITY(a.get_number_of_parameters(), b.get_number_of_parameters());

   TEST_EQUALITY(a.get_g1(), b.get_g1());
   TEST_EQUALITY(a.get_g2(), b.get_g2());
   TEST_EQUALITY(a.get_g3(), b.get_g3());

   TEST_EQUALITY(a.get_Yu(), b.get_Yu());
   TEST_EQUALITY(a.get_Yd(), b.get_Yd());
   TEST_EQUALITY(a.get_Ye(), b.get_Ye());

   TEST_EQUALITY(a.get_MassB(), b.get_MassB());
   TEST_EQUALITY(a.get_MassWB(), b.get_MassWB());
   TEST_EQUALITY(a.get_MassG(), b.get_MassG());

   TEST_EQUALITY(a.get_mHd2(), b.get_mHd2());
   TEST_EQUALITY(a.get_mHu2(), b.get_mHu2());
   TEST_EQUALITY(a.get_mq2(), b.get_mq2());
   TEST_EQUALITY(a.get_mu2(), b.get_mu2());
   TEST_EQUALITY(a.get_md2(), b.get_md2());
   TEST_EQUALITY(a.get_ml2(), b.get_ml2());
   TEST_EQUALITY(a.get_me2(), b.get_me2());

   TEST_EQUALITY(a.get_TYu(), b.get_TYu());
   TEST_EQUALITY(a.get_TYd(), b.get_TYd());
   TEST_EQUALITY(a.get_TYe(), b.get_TYe());

   TEST_EQUALITY(a.get_Mu(), b.get_Mu());
   TEST_CLOSE(a.get_BMu(), b.get_BMu(), 1e-10);

   TEST_EQUALITY(a.get_vu(), b.get_vu());
   TEST_EQUALITY(a.get_vd(), b.get_vd());
}

template <class T, class P>
void test_mass_equality(const T& a, const P& b)
{
   TEST_EQUALITY(a.Mhh, b.Mhh);
   TEST_EQUALITY(a.MAh, b.MAh);
   TEST_EQUALITY(a.MHpm, b.MHpm);
   TEST_EQUALITY(a.MCha, b.MCha);
   TEST_EQUALITY(a.MChi, b.MChi);
   TEST_EQUALITY(a.MGlu, b.MGlu);
   TEST_EQUALITY(a.MSu, b.MSu);
   TEST_EQUALITY(a.MSd, b.MSd);
   TEST_EQUALITY(a.MSe, b.MSe);
   TEST_EQUALITY(a.MSv, b.MSv);
   TEST_EQUALITY(a.MFu, b.MFu);
   TEST_EQUALITY(a.MFd, b.MFd);
   TEST_EQUALITY(a.MFe, b.MFe);
   TEST_EQUALITY(a.MFv, b.MFv);
   TEST_EQUALITY(a.MVP, b.MVP);
   TEST_EQUALITY(a.MVG, b.MVG);
   TEST_EQUALITY(a.MVWm, b.MVWm);
   TEST_EQUALITY(a.MVZ, b.MVZ);

   TEST_EQUALITY(a.ZH, b.ZH);
   TEST_EQUALITY(a.ZA, b.ZA);
   TEST_EQUALITY(a.ZP, b.ZP);
   TEST_EQUALITY(a.ZN, b.ZN);
   TEST_EQUALITY(a.UM, b.UM);
   TEST_EQUALITY(a.UP, b.UP);
   TEST_EQUALITY(a.ZU, b.ZU);
   TEST_EQUALITY(a.ZD, b.ZD);
   TEST_EQUALITY(a.ZE, b.ZE);
   TEST_EQUALITY(a.ZV, b.ZV);
   TEST_EQUALITY(a.ZEL, b.ZEL);
   TEST_EQUALITY(a.ZER, b.ZER);
   TEST_EQUALITY(a.ZUL, b.ZUL);
   TEST_EQUALITY(a.ZUR, b.ZUR);
   TEST_EQUALITY(a.ZDL, b.ZDL);
   TEST_EQUALITY(a.ZDR, b.ZDR);
}

template <class T, class P>
void test_input_parameter_equality(const T& a, const P& b)
{
   TEST_EQUALITY(a.m0, b.m0);
   TEST_EQUALITY(a.m12, b.m12);
   TEST_EQUALITY(a.TanBeta, b.TanBeta);
   TEST_EQUALITY(a.SignMu, b.SignMu);
   TEST_EQUALITY(a.Azero, b.Azero);
}

BOOST_AUTO_TEST_CASE( test_CMSSM_read_write )
{
   softsusy::QedQcd qedqcd1, qedqcd2;
   qedqcd1.setAlpha(softsusy::ALPHA, 1./127.);
   qedqcd1.toMz();

   CMSSM_input_parameters pp;
   pp.m0 = 125.;
   pp.m12 = 500.;
   pp.TanBeta = 10.;
   pp.SignMu = 1;
   pp.Azero = 0.1;

   CMSSM_mass_eigenstates model(pp);
   setup_CMSSM(model, pp);

   const std::string db_file("test/test_CMSSM_database.db");

   CMSSM_database::to_database(db_file, model, &qedqcd1);

   const CMSSM_mass_eigenstates tmp(CMSSM_database::from_database(db_file, -1, &qedqcd2));

   gErrors = 0;
   test_parameter_equality(model, tmp);
   test_mass_equality(model.get_physical(), tmp.get_physical());
   test_input_parameter_equality(model.get_input(), tmp.get_input());
   BOOST_REQUIRE(gErrors == 0);
   BOOST_REQUIRE((qedqcd1.display_input() - qedqcd2.display_input()).maxCoeff() < 1e-10);
}
