
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_numerics

#include <boost/test/unit_test.hpp>
#include "numerics2.hpp"
#include "stopwatch.hpp"
#include <complex>
#include <cstdlib>

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE(test_is_equal)
{
   short s = (short)0;
   int i = 0;
   long l = 0L;
   long long ll = 0LL;
   unsigned short us = 0;
   unsigned int ui = 0;
   unsigned long ul = 0UL;
   unsigned long long ull = 0ULL;

   BOOST_CHECK(is_zero(s));
   BOOST_CHECK(is_zero(i));
   BOOST_CHECK(is_zero(l));
   BOOST_CHECK(is_zero(ll));
   BOOST_CHECK(is_zero(us));
   BOOST_CHECK(is_zero(ui));
   BOOST_CHECK(is_zero(ul));
   BOOST_CHECK(is_zero(ull));

   // BOOST_CHECK(is_equal(s, (short)0));
   BOOST_CHECK(is_equal(i, 0));
   BOOST_CHECK(is_equal(l, 0L));
   BOOST_CHECK(is_equal(ll, 0LL));
   // BOOST_CHECK(is_equal(us, (unsigned short)0));
   BOOST_CHECK(is_equal(ui, 0U));
   BOOST_CHECK(is_equal(ul, 0UL));
   BOOST_CHECK(is_equal(ull, 0ULL));

   // BOOST_CHECK(is_equal_rel(s, (short)0));
   BOOST_CHECK(is_equal_rel(i, 0));
   BOOST_CHECK(is_equal_rel(l, 0L));
   BOOST_CHECK(is_equal_rel(ll, 0LL));
   // BOOST_CHECK(is_equal_rel(us, (unsigned short)0));
   BOOST_CHECK(is_equal_rel(ui, 0U));
   BOOST_CHECK(is_equal_rel(ul, 0UL));
   BOOST_CHECK(is_equal_rel(ull, 0ULL));
}

template <long N>
std::array<std::complex<double>, N> make_logs(double max = 1.)
{
   std::array<std::complex<double>, N> vals;

   for (auto& v: vals)
      v = std::complex<double>(max*std::rand()/RAND_MAX, max*std::rand()/RAND_MAX);

   return vals;
}

BOOST_AUTO_TEST_CASE(test_complex_log)
{
   const auto vals = make_logs<100000>();

   for (auto v: vals) {
      BOOST_CHECK_CLOSE_FRACTION(std::real(std::log(v)), std::real(fast_log(v)), 1e-9);
      BOOST_CHECK_CLOSE_FRACTION(std::imag(std::log(v)), std::imag(fast_log(v)), 1e-9);
   }
}

BOOST_AUTO_TEST_CASE(test_complex_log_bench)
{
   const auto vals = make_logs<100000>();

   Stopwatch sw;

   sw.start();
   for (auto v: vals)
      volatile auto res = std::log(v);
   sw.stop();
   const double t_log = sw.get_time_in_seconds();

   sw.start();
   for (auto v: vals)
      volatile auto res = fast_log(v);
   sw.stop();
   const double t_fslog = sw.get_time_in_seconds();

   BOOST_TEST_MESSAGE("time for log     : " << t_log);
   BOOST_TEST_MESSAGE("time for fast_log: " << t_fslog);

   BOOST_CHECK_LT(t_fslog, t_log);
}

BOOST_AUTO_TEST_CASE( test_complex_log1p )
{
   const double tol = 1.e-14;

   BOOST_CHECK_EQUAL(complex_log1p(std::complex<double>(0, 0)), std::complex<double>(0., 0.));

   const auto vals = make_logs<100000>(0.5);
   for (auto v: vals) {
      const auto x = std::real(v);
      const auto y = std::imag(v);
      const auto theta = std::atan2(y, 1 + x);

      BOOST_CHECK_CLOSE_FRACTION(0.5 * std::log1p(x * x + 2. * x + y * y),
                                 std::real(complex_log1p(v)), tol);
      BOOST_CHECK_CLOSE_FRACTION(theta, std::imag(complex_log1p(v)), tol);
   }

   // precalculated checks
   const std::array<double, 81> real_parts{{
         -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.15, -0.15,
            -0.15, -0.15, -0.15, -0.15, -0.15, -0.15, -0.15, -0.1, -0.1, -0.1,
            -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.05, -0.05, -0.05, -0.05,
            -0.05, -0.05, -0.05, -0.05, -0.05, 0., 0., 0., 0., 0., 0., 0., 0.,
            0., 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.1, 0.1,
            0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.15, 0.15, 0.15, 0.15, 0.15,
            0.15, 0.15, 0.15, 0.15, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2
            }};
   const std::array<double, 81> imag_parts{{
         -0.2, -0.15, -0.1, -0.05, 0., 0.05, 0.1, 0.15, 0.2, -0.2, -0.15,
            -0.1, -0.05, 0., 0.05, 0.1, 0.15, 0.2, -0.2, -0.15, -0.1, -0.05, 0.,
            0.05, 0.1, 0.15, 0.2, -0.2, -0.15, -0.1, -0.05, 0., 0.05, 0.1, 0.15,
            0.2, -0.2, -0.15, -0.1, -0.05, 0., 0.05, 0.1, 0.15, 0.2, -0.2, -0.15,
            -0.1, -0.05, 0., 0.05, 0.1, 0.15, 0.2, -0.2, -0.15, -0.1, -0.05, 0.,
            0.05, 0.1, 0.15, 0.2, -0.2, -0.15, -0.1, -0.05, 0., 0.05, 0.1, 0.15,
            0.2, -0.2, -0.15, -0.1, -0.05, 0., 0.05, 0.1, 0.15, 0.2
            }};
   const std::array<std::complex<double>, 81> expected{{
         std::complex<double>(-0.19283124040599233447599202428962, -0.24497866312686415417208248121128),
            std::complex<double>(-0.20586736056087988906260934001492, -0.18534794999569476488602596122854),
            std::complex<double>(-0.21539145804622712869086806728861, -0.12435499454676143503135484916387),
            std::complex<double>(-0.22119423110638109425932641876191, -0.06241880999595734847397911298551),
            std::complex<double>(-0.22314355131420975576629509030983, 0),
            std::complex<double>(-0.22119423110638109425932641876191, 0.06241880999595734847397911298551),
            std::complex<double>(-0.21539145804622712869086806728861, 0.12435499454676143503135484916387),
            std::complex<double>(-0.20586736056087988906260934001492, 0.18534794999569476488602596122854),
            std::complex<double>(-0.19283124040599233447599202428962, 0.24497866312686415417208248121128),
            std::complex<double>(-0.13557638525028518175914935781664, -0.23109066719589707990301793905020),
            std::complex<double>(-0.14718553030128876823039337786139, -0.17467219900823969307190196932964),
            std::complex<double>(-0.15564596904545734046530344473629, -0.11710874456686428298901824354029),
            std::complex<double>(-0.16079181206373113783459183261940, -0.05875582271572269290218330245177),
            std::complex<double>(-0.16251892949777491318568895826941, 0),
            std::complex<double>(-0.16079181206373113783459183261940, 0.05875582271572269290218330245177),
            std::complex<double>(-0.15564596904545734046530344473629, 0.11710874456686428298901824354029),
            std::complex<double>(-0.14718553030128876823039337786139, 0.17467219900823969307190196932964),
            std::complex<double>(-0.13557638525028518175914935781664, 0.23109066719589707990301793905020),
            std::complex<double>(-0.08125946474888745659284447913471, -0.21866894587394196204217375024994),
            std::complex<double>(-0.09166102856376907985593050370429, -0.16514867741462683827912828964394),
            std::complex<double>(-0.09922546936191912737599370743657, -0.11065722117389564655913987221006),
            std::complex<double>(-0.103819682389122250807720522133694, -0.055498505245716835557198148092237),
            std::complex<double>(-0.10536051565782630122750098083931, 0),
            std::complex<double>(-0.103819682389122250807720522133694, 0.055498505245716835557198148092237),
            std::complex<double>(-0.09922546936191912737599370743657, 0.11065722117389564655913987221006),
            std::complex<double>(-0.09166102856376907985593050370429, 0.16514867741462683827912828964394),
            std::complex<double>(-0.08125946474888745659284447913471, 0.21866894587394196204217375024994),
            std::complex<double>(-0.02960967982998561181684383917893, -0.20749622643520266494202316381465),
            std::complex<double>(-0.03898077073485592924218001328464, -0.15660187698201535512227632471472),
            std::complex<double>(-0.045783596762745241588789480810176, -0.104876938730233895818336167534754),
            std::complex<double>(-0.049910167641105465307063996096948, -0.052583061610941717974868773085595),
            std::complex<double>(-0.051293294387550533426196144254687, 0),
            std::complex<double>(-0.049910167641105465307063996096948, 0.052583061610941717974868773085595),
            std::complex<double>(-0.045783596762745241588789480810176, 0.104876938730233895818336167534754),
            std::complex<double>(-0.03898077073485592924218001328464, 0.15660187698201535512227632471472),
            std::complex<double>(-0.02960967982998561181684383917893, 0.20749622643520266494202316381465),
            std::complex<double>(0.01961035657664064813460044828556, -0.19739555984988075837004976519479),
            std::complex<double>(0.01112530446740987942153752842165, -0.14888994760949725058653039165587),
            std::complex<double>(0.004975165426584041424107678772130, -0.099668652491162027378446119878021),
            std::complex<double>(0.001248440099293599490115967532134, -0.049958395721942761410006287034845),
            std::complex<double>(0, 0),
            std::complex<double>(0.001248440099293599490115967532134, 0.049958395721942761410006287034845),
            std::complex<double>(0.004975165426584041424107678772130, 0.099668652491162027378446119878021),
            std::complex<double>(0.01112530446740987942153752842165, 0.14888994760949725058653039165587),
            std::complex<double>(0.01961035657664064813460044828556, 0.19739555984988075837004976519479),
            std::complex<double>(0.06660942189311135616134211443420, -0.18822150530477076117135427387118),
            std::complex<double>(0.05889151782819172726939705473526, -0.14189705460416392281285161710255),
            std::complex<double>(0.053304867529129113024063860805478, -0.094951706342756319757251984940447),
            std::complex<double>(0.049922667484858069424903514305770, -0.047583103276983395802032716016486),
            std::complex<double>(0.048790164169432003065374404223165, 0),
            std::complex<double>(0.049922667484858069424903514305770, 0.047583103276983395802032716016486),
            std::complex<double>(0.053304867529129113024063860805478, 0.094951706342756319757251984940447),
            std::complex<double>(0.05889151782819172726939705473526, 0.14189705460416392281285161710255),
            std::complex<double>(0.06660942189311135616134211443420, 0.18822150530477076117135427387118),
            std::complex<double>(0.11157177565710487788314754515492, -0.17985349979247827058855299725611),
            std::complex<double>(0.10452231346735406028117974897498, -0.13552771398550073213150543559916),
            std::complex<double>(0.099425429372582595066319157757531, -0.090659887200745113498386675308466),
            std::complex<double>(0.096342171914750604923516901322534, -0.045423279421577015023084402434214),
            std::complex<double>(0.095310179804324860043952123280765, 0),
            std::complex<double>(0.096342171914750604923516901322534, 0.045423279421577015023084402434214),
            std::complex<double>(0.099425429372582595066319157757531, 0.090659887200745113498386675308466),
            std::complex<double>(0.10452231346735406028117974897498, 0.13552771398550073213150543559916),
            std::complex<double>(0.11157177565710487788314754515492, 0.17985349979247827058855299725611),
            std::complex<double>(0.15466062377763104405381431657515, -0.17219081452293902421865161002832),
            std::complex<double>(0.14819700652690121438336907483417, -0.12970253715591210695114112663540),
            std::complex<double>(0.14352843852893127652490183115882, -0.08673833867598511181090989298473),
            std::complex<double>(0.14070622971909276564600672071417, -0.04345089539153084207317434412599),
            std::complex<double>(0.13976194237515869737152925566766, 0),
            std::complex<double>(0.14070622971909276564600672071417, 0.04345089539153084207317434412599),
            std::complex<double>(0.14352843852893127652490183115882, 0.08673833867598511181090989298473),
            std::complex<double>(0.14819700652690121438336907483417, 0.12970253715591210695114112663540),
            std::complex<double>(0.15466062377763104405381431657515, 0.17219081452293902421865161002832),
            std::complex<double>(0.19602104388801184758328850228954, -0.16514867741462683827912828964394),
            std::complex<double>(0.19007365006193725328714504817574, -0.12435499454676143503135484916387),
            std::complex<double>(0.18578177821624151687402422810969, -0.08314123188844122991066831465078),
            std::complex<double>(0.18318885970005878730306407749280, -0.04164257909858842386017806470075),
            std::complex<double>(0.18232155679395462621171802515451, 0),
            std::complex<double>(0.18318885970005878730306407749280, 0.04164257909858842386017806470075),
            std::complex<double>(0.18578177821624151687402422810969, 0.08314123188844122991066831465078),
            std::complex<double>(0.19007365006193725328714504817574, 0.12435499454676143503135484916387),
            std::complex<double>(0.19602104388801184758328850228954, 0.16514867741462683827912828964394)
            }};

   for (std::size_t i = 0; i < real_parts.size(); ++i) {
      BOOST_CHECK_CLOSE_FRACTION(std::real(complex_log1p(std::complex<double>(real_parts[i], imag_parts[i]))),
                                 std::real(expected[i]), 1.e-15);
      BOOST_CHECK_CLOSE_FRACTION(std::imag(complex_log1p(std::complex<double>(real_parts[i], imag_parts[i]))),
                                 std::imag(expected[i]), 1.e-15);
   }
}
