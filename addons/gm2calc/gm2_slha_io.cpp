// ====================================================================
// This file is part of GM2Calc.
//
// GM2Calc is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// GM2Calc is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with GM2Calc.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#include "gm2_slha_io.hpp"
#include "ffunctions.hpp"
#include "MSSMNoFV_onshell.hpp"

#include <cmath>
#include <fstream>
#include <limits>

#include <Eigen/Core>

namespace flexiblesusy {
namespace gm2calc {

/**
 * @brief reads from source
 *
 * If source is "-", then read_from_stream() is called.  Otherwise,
 * read_from_file() is called.
 *
 * @param source string that specifies the source
 */
void GM2_slha_io::read_from_source(const std::string& source)
{
   if (source == "-")
      read_from_stream(std::cin);
   else
      read_from_file(source);
}

/**
 * @brief opens SLHA input file and reads the content
 * @param file_name SLHA input file name
 */
void GM2_slha_io::read_from_file(const std::string& file_name)
{
   std::ifstream ifs(file_name);
   if (ifs.good()) {
      data.clear();
      data.read(ifs);
   } else {
      std::ostringstream msg;
      msg << "cannot read SLHA file: \"" << file_name << "\"";
      throw ReadError(msg.str());
   }
}

/**
 * @brief reads SLHA data from a stream
 * @param istr input stream
 */
void GM2_slha_io::read_from_stream(std::istream& istr)
{
   data.read(istr);
}

double GM2_slha_io::read_entry(const std::string& block_name, int key) const
{
   SLHAea::Coll::const_iterator block =
      data.find(data.cbegin(), data.cend(), block_name);

   double entry = 0.;
   const SLHAea::Block::key_type keys(1, std::to_string(key));
   SLHAea::Block::const_iterator line;

   while (block != data.cend()) {
      line = block->find(keys);

      if (line != block->end() && line->is_data_line() && line->size() > 1)
         entry = convert_to<double>(line->at(1));

      ++block;
      block = data.find(block, data.cend(), block_name);
   }

   return entry;
}

/**
 * Reads scale definition from SLHA block.
 *
 * @param block_name block name
 *
 * @return scale (or 0 if no scale is defined)
 */
double GM2_slha_io::read_scale(const std::string& block_name) const
{
   if (!block_exists(block_name))
      return 0.;

   double scale = 0.;

   for (SLHAea::Block::const_iterator line = data.at(block_name).cbegin(),
        end = data.at(block_name).cend(); line != end; ++line) {
      if (!line->is_data_line()) {
         if (line->size() > 3 &&
             to_lower((*line)[0]) == "block" && (*line)[2] == "Q=")
            scale = convert_to<double>((*line)[3]);
         break;
      }
   }

   return scale;
}

bool GM2_slha_io::block_exists(const std::string& block_name) const
{
   return data.find(block_name) != data.cend();
}

std::string GM2_slha_io::to_lower(const std::string& str)
{
   std::string lower(str.size(), ' ');
   std::transform(str.begin(), str.end(), lower.begin(), ::tolower);
   return lower;
}

/**
 * Applies processor to each (key, value) pair of a SLHA block.
 * Non-data lines are ignored.
 *
 * @param block_name block name
 * @param processor tuple processor to be applied
 *
 * @return scale (or 0 if no scale is defined)
 */
double GM2_slha_io::read_block(const std::string& block_name, const Tuple_processor& processor) const
{
   SLHAea::Coll::const_iterator block =
      data.find(data.cbegin(), data.cend(), block_name);

   double scale = 0.;

   while (block != data.cend()) {
      for (SLHAea::Block::const_iterator line = block->cbegin(),
              end = block->cend(); line != end; ++line) {
         if (!line->is_data_line()) {
            // read scale from block definition
            if (line->size() > 3 &&
                to_lower((*line)[0]) == "block" && (*line)[2] == "Q=")
               scale = convert_to<double>((*line)[3]);
            continue;
         }

         if (line->size() >= 2) {
            const int key = convert_to<int>((*line)[0]);
            const double value = convert_to<double>((*line)[1]);
            processor(key, value);
         }
      }

      ++block;
      block = data.find(block, data.cend(), block_name);
   }

   return scale;
}

void GM2_slha_io::set_block(const std::ostringstream& lines, Position position)
{
   SLHAea::Block block;
   block.str(lines.str());
   data.erase(block.name());
   if (position == front)
      data.push_front(block);
   else
      data.push_back(block);
}

double read_scale(const GM2_slha_io& slha_io)
{
   char const * const drbar_blocks[] =
      { "Yu", "Yd", "Ye", "Ae", "Ad", "Au", "HMIX", "MSOFT" };

   double scale = 0.;

   for (unsigned i = 0; i < sizeof(drbar_blocks)/sizeof(*drbar_blocks); i++) {
      const double block_scale = slha_io.read_scale(drbar_blocks[i]);
      if (!is_zero(block_scale)) {
         scale = block_scale;
         break;
      }
   }

   if (is_zero(scale)) {
      std::cerr << "could not find renormalization scale in any"
         " SLHA block.\n";
   }

   return scale;
}

void fill_alpha_s(const GM2_slha_io& slha_io, MSSMNoFV_onshell& model)
{
   const double alpha_S = slha_io.read_entry("SMINPUTS", 3);
   model.set_g3(std::sqrt(4*M_PI*alpha_S));
}

void fill_drbar_parameters(const GM2_slha_io& slha_io, MSSMNoFV_onshell& model)
{
   {
      Eigen::Matrix<double,3,3> Ae(Eigen::Matrix<double,3,3>::Zero());
      slha_io.read_block("AE", Ae);
      model.set_Ae(Ae);
   }
   {
      Eigen::Matrix<double,3,3> Au(Eigen::Matrix<double,3,3>::Zero());
      slha_io.read_block("AU", Au);
      model.set_Au(Au);
   }
   {
      Eigen::Matrix<double,3,3> Ad(Eigen::Matrix<double,3,3>::Zero());
      slha_io.read_block("AD", Ad);
      model.set_Ad(Ad);
   }

   model.set_Mu(slha_io.read_entry("HMIX", 1));
   model.set_mHd2(slha_io.read_entry("MSOFT", 21));
   model.set_mHu2(slha_io.read_entry("MSOFT", 22));
   model.set_ml2(0, 0, signed_sqr(slha_io.read_entry("MSOFT", 31)));
   model.set_ml2(1, 1, signed_sqr(slha_io.read_entry("MSOFT", 32)));
   model.set_ml2(2, 2, signed_sqr(slha_io.read_entry("MSOFT", 33)));
   model.set_me2(0, 0, signed_sqr(slha_io.read_entry("MSOFT", 34)));
   model.set_me2(1, 1, signed_sqr(slha_io.read_entry("MSOFT", 35)));
   model.set_me2(2, 2, signed_sqr(slha_io.read_entry("MSOFT", 36)));
   model.set_mq2(0, 0, signed_sqr(slha_io.read_entry("MSOFT", 41)));
   model.set_mq2(1, 1, signed_sqr(slha_io.read_entry("MSOFT", 42)));
   model.set_mq2(2, 2, signed_sqr(slha_io.read_entry("MSOFT", 43)));
   model.set_mu2(0, 0, signed_sqr(slha_io.read_entry("MSOFT", 44)));
   model.set_mu2(1, 1, signed_sqr(slha_io.read_entry("MSOFT", 45)));
   model.set_mu2(2, 2, signed_sqr(slha_io.read_entry("MSOFT", 46)));
   model.set_md2(0, 0, signed_sqr(slha_io.read_entry("MSOFT", 47)));
   model.set_md2(1, 1, signed_sqr(slha_io.read_entry("MSOFT", 48)));
   model.set_md2(2, 2, signed_sqr(slha_io.read_entry("MSOFT", 49)));
   model.set_MassB(slha_io.read_entry("MSOFT", 1));
   model.set_MassWB(slha_io.read_entry("MSOFT", 2));
   model.set_MassG(slha_io.read_entry("MSOFT", 3));

   const double tanb = slha_io.read_entry("HMIX", 2);
   const double MA2_drbar = slha_io.read_entry("HMIX", 4);
   const double MW = model.get_MW();
   const double MZ = model.get_MZ();
   const double cW = MW/MZ;
   const double sW = std::sqrt(1. - cW*cW);
   const double vev = 2. * model.get_MW() * sW / model.get_EL();
   const double sinb = tanb / std::sqrt(1 + tanb*tanb);
   const double cosb = 1.   / std::sqrt(1 + tanb*tanb);

   model.set_vd(vev * cosb);
   model.set_vu(vev * sinb);
   model.set_BMu(MA2_drbar * sinb * cosb);

   model.set_scale(read_scale(slha_io));
}

void fill_fermion_pole_masses_from_sminputs(const GM2_slha_io& slha_io, MSSMNoFV_onshell_physical& physical)
{
   physical.MFd = slha_io.read_entry("SMINPUTS", 21);
   physical.MFs = slha_io.read_entry("SMINPUTS", 23);
   physical.MFb = slha_io.read_entry("SMINPUTS", 5);
   physical.MFu = slha_io.read_entry("SMINPUTS", 22);
   physical.MFc = slha_io.read_entry("SMINPUTS", 24);
   physical.MFt = slha_io.read_entry("SMINPUTS", 6);
   physical.MFe = slha_io.read_entry("SMINPUTS", 11);
   physical.MFm = slha_io.read_entry("SMINPUTS", 13);
   physical.MFtau = slha_io.read_entry("SMINPUTS", 7);
}

void fill_physical(const GM2_slha_io& slha_io, MSSMNoFV_onshell_physical& physical)
{
   // read MW from MASS[24].  If not given there, read from
   // SMINPUTS[9]
   double MW = slha_io.read_entry("MASS", 24);
   if (is_zero(MW))
      MW = slha_io.read_entry("SMINPUTS", 9);

   fill_fermion_pole_masses_from_sminputs(slha_io, physical);

   physical.MVWm = MW;
   physical.MVZ = slha_io.read_entry("SMINPUTS", 4);
   physical.MSveL = slha_io.read_entry("MASS", 1000012);
   physical.MSvmL = slha_io.read_entry("MASS", 1000014);
   physical.MSvtL = slha_io.read_entry("MASS", 1000016);
   physical.MSd(0) = slha_io.read_entry("MASS", 1000001);
   physical.MSd(1) = slha_io.read_entry("MASS", 2000001);
   physical.MSu(0) = slha_io.read_entry("MASS", 1000002);
   physical.MSu(1) = slha_io.read_entry("MASS", 2000002);
   physical.MSe(0) = slha_io.read_entry("MASS", 1000011);
   physical.MSe(1) = slha_io.read_entry("MASS", 2000011);
   physical.MSm(0) = slha_io.read_entry("MASS", 1000013);
   physical.MSm(1) = slha_io.read_entry("MASS", 2000013);
   physical.MStau(0) = slha_io.read_entry("MASS", 1000015);
   physical.MStau(1) = slha_io.read_entry("MASS", 2000015);
   physical.MSs(0) = slha_io.read_entry("MASS", 1000003);
   physical.MSs(1) = slha_io.read_entry("MASS", 2000003);
   physical.MSc(0) = slha_io.read_entry("MASS", 1000004);
   physical.MSc(1) = slha_io.read_entry("MASS", 2000004);
   physical.MSb(0) = slha_io.read_entry("MASS", 1000005);
   physical.MSb(1) = slha_io.read_entry("MASS", 2000005);
   physical.MSt(0) = slha_io.read_entry("MASS", 1000006);
   physical.MSt(1) = slha_io.read_entry("MASS", 2000006);
   physical.Mhh(0) = slha_io.read_entry("MASS", 25);
   physical.Mhh(1) = slha_io.read_entry("MASS", 35);
   physical.MAh(1) = slha_io.read_entry("MASS", 36);
   physical.MHpm(1) = slha_io.read_entry("MASS", 37);
   physical.MChi(0) = slha_io.read_entry("MASS", 1000022);
   physical.MChi(1) = slha_io.read_entry("MASS", 1000023);
   physical.MChi(2) = slha_io.read_entry("MASS", 1000025);
   physical.MChi(3) = slha_io.read_entry("MASS", 1000035);
   physical.MCha(0) = slha_io.read_entry("MASS", 1000024);
   physical.MCha(1) = slha_io.read_entry("MASS", 1000037);
}

void fill_pole_masses_from_sminputs_and_mass(const GM2_slha_io& slha_io, MSSMNoFV_onshell& model)
{
   MSSMNoFV_onshell_physical physical_hk;
   fill_physical(slha_io, physical_hk);
   physical_hk.convert_to_hk();
   model.get_physical() = physical_hk;
}

void fill_pole_masses_from_sminputs(const GM2_slha_io& slha_io, MSSMNoFV_onshell& model)
{
   fill_fermion_pole_masses_from_sminputs(slha_io, model.get_physical());

   model.get_physical().MVWm = slha_io.read_entry("SMINPUTS", 9);
   model.get_physical().MVZ = slha_io.read_entry("SMINPUTS", 4);
}

void fill_gm2_specific_alphas(const GM2_slha_io& slha_io, MSSMNoFV_onshell& model)
{
   const double alpha_MZ = std::abs(slha_io.read_entry("FlexibleSUSYGM2Input", 1));
   const double alpha_thompson = std::abs(slha_io.read_entry("FlexibleSUSYGM2Input", 2));

   if (alpha_MZ > std::numeric_limits<double>::epsilon())
      model.set_alpha_MZ(alpha_MZ);

   if (alpha_thompson > std::numeric_limits<double>::epsilon())
      model.set_alpha_thompson(alpha_thompson);
}

void fill_gm2_specific_onshell_parameters(const GM2_slha_io& slha_io, MSSMNoFV_onshell& model)
{
   const double tanb = slha_io.read_entry("FlexibleSUSYGM2Input", 3);
   const double MW = model.get_MW();
   const double MZ = model.get_MZ();
   const double cW = MW/MZ;
   const double sW = std::sqrt(1. - cW*cW);
   const double vev = 2. * MW * sW / model.get_EL();
   const double sinb = tanb / std::sqrt(1 + tanb*tanb);
   const double cosb = 1.   / std::sqrt(1 + tanb*tanb);

   model.set_vd(vev * cosb);
   model.set_vu(vev * sinb);

   model.set_scale(         slha_io.read_entry("FlexibleSUSYGM2Input", 0));
   model.set_alpha_MZ(      slha_io.read_entry("FlexibleSUSYGM2Input", 1));
   model.set_alpha_thompson(slha_io.read_entry("FlexibleSUSYGM2Input", 2));
   model.set_Mu(            slha_io.read_entry("FlexibleSUSYGM2Input", 4));
   model.set_MassB(         slha_io.read_entry("FlexibleSUSYGM2Input", 5));
   model.set_MassWB(        slha_io.read_entry("FlexibleSUSYGM2Input", 6));
   model.set_MassG(         slha_io.read_entry("FlexibleSUSYGM2Input", 7));
   model.set_BMu(           slha_io.read_entry("FlexibleSUSYGM2Input", 8));
   model.set_ml2(0, 0,      signed_sqr(slha_io.read_entry("FlexibleSUSYGM2Input", 9 )));
   model.set_ml2(1, 1,      signed_sqr(slha_io.read_entry("FlexibleSUSYGM2Input", 10)));
   model.set_ml2(2, 2,      signed_sqr(slha_io.read_entry("FlexibleSUSYGM2Input", 11)));
   model.set_me2(0, 0,      signed_sqr(slha_io.read_entry("FlexibleSUSYGM2Input", 12)));
   model.set_me2(1, 1,      signed_sqr(slha_io.read_entry("FlexibleSUSYGM2Input", 13)));
   model.set_me2(2, 2,      signed_sqr(slha_io.read_entry("FlexibleSUSYGM2Input", 14)));
   model.set_mq2(0, 0,      signed_sqr(slha_io.read_entry("FlexibleSUSYGM2Input", 15)));
   model.set_mq2(1, 1,      signed_sqr(slha_io.read_entry("FlexibleSUSYGM2Input", 16)));
   model.set_mq2(2, 2,      signed_sqr(slha_io.read_entry("FlexibleSUSYGM2Input", 17)));
   model.set_mu2(0, 0,      signed_sqr(slha_io.read_entry("FlexibleSUSYGM2Input", 18)));
   model.set_mu2(1, 1,      signed_sqr(slha_io.read_entry("FlexibleSUSYGM2Input", 19)));
   model.set_mu2(2, 2,      signed_sqr(slha_io.read_entry("FlexibleSUSYGM2Input", 20)));
   model.set_md2(0, 0,      signed_sqr(slha_io.read_entry("FlexibleSUSYGM2Input", 21)));
   model.set_md2(1, 1,      signed_sqr(slha_io.read_entry("FlexibleSUSYGM2Input", 22)));
   model.set_md2(2, 2,      signed_sqr(slha_io.read_entry("FlexibleSUSYGM2Input", 23)));
   model.set_Ae( 0, 0,      slha_io.read_entry("FlexibleSUSYGM2Input", 24));
   model.set_Ae( 1, 1,      slha_io.read_entry("FlexibleSUSYGM2Input", 25));
   model.set_Ae( 2, 2,      slha_io.read_entry("FlexibleSUSYGM2Input", 26));
   model.set_Ad( 0, 0,      slha_io.read_entry("FlexibleSUSYGM2Input", 27));
   model.set_Ad( 1, 1,      slha_io.read_entry("FlexibleSUSYGM2Input", 28));
   model.set_Ad( 2, 2,      slha_io.read_entry("FlexibleSUSYGM2Input", 29));
   model.set_Au( 0, 0,      slha_io.read_entry("FlexibleSUSYGM2Input", 30));
   model.set_Au( 1, 1,      slha_io.read_entry("FlexibleSUSYGM2Input", 31));
   model.set_Au( 2, 2,      slha_io.read_entry("FlexibleSUSYGM2Input", 32));
}

void fill_gm2calc(const GM2_slha_io& slha_io, MSSMNoFV_onshell& model)
{
   fill_pole_masses_from_sminputs(slha_io, model);
   fill_alpha_s(slha_io, model);
   fill_gm2_specific_alphas(slha_io, model);
   fill_gm2_specific_onshell_parameters(slha_io, model);
}

void fill_slha(const GM2_slha_io& slha_io, MSSMNoFV_onshell& model)
{
   fill_pole_masses_from_sminputs_and_mass(slha_io, model);
   fill_alpha_s(slha_io, model);
   fill_drbar_parameters(slha_io, model);
   fill_gm2_specific_alphas(slha_io, model);
}

} // namespace gm2calc
} // namespace flexiblesusy