(* ::Package:: *)

(* :Copyright:

   ====================================================================
   This file is part of FlexibleSUSY.

   FlexibleSUSY is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published
   by the Free Software Foundation, either version 3 of the License,
   or (at your option) any later version.

   FlexibleSUSY is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with FlexibleSUSY.  If not, see
   <http://www.gnu.org/licenses/>.
   ====================================================================

*)

BeginPackage["FToFConversionInNucleus`",
    {"SARAH`", "TextFormatting`", "TreeMasses`", "Vertices`", "CXXDiagrams`"}
];

FToFConversionInNucleusCreateInterface::usage = "";

Begin["Private`"];
(*
(* we pass the generation index but a particle might also have other indices, like color, so we need to construct full vector of indices *)
CreateIndexArray[particle_, idxName_String, generationIndex_String] :=
    Module[{numberOfIndices = CXXDiagrams`NumberOfFieldIndices[particle]},
        "std::array<int, " <> ToString @ numberOfIndices <>
            "> " <> idxName <> " = {" <>
            If[TreeMasses`GetDimension[particle] =!= 1,
                generationIndex <>
                    If[numberOfIndices1 =!= 1,
                        StringJoin @ Table[", 0", {numberOfIndices-1}],
                        ""] <> " ",
                If[numberOfIndices =!= 0,
                    StringJoin @ Riffle[Table[" 0", {numberOfIndices}], ","] <> " ",
                    ""]
            ] <> "};\n"
    ];
*)
FToFConversionInNucleusCreateInterface[{inFermion_, outFermion_, nucleus_}] :=
    Module[{prototype, definition},

        prototype =
            "double calculate_" <> CXXNameOfField[inFermion] <> "_to_" <>
                CXXNameOfField[outFermion] <> "_in_nucleus (\n" <>
                If[TreeMasses`GetDimension[inFermion] =!= 1,
                    "int generationIndex1, ",
                    " "
                ] <>
                If[TreeMasses`GetDimension[outFermion] =!= 1,
                    " int generationIndex2, ",
                    " "
                ] <>
                "const " <> FlexibleSUSY`FSModelName <> "_f_to_f_conversion::Nucleus nucleus,\n" <>
                "const " <> FlexibleSUSY`FSModelName <> "_mass_eigenstates& model, const softsusy::QedQcd& qedqcd)";

        definition =
            prototype <> " {\n" <>
            IndentText[
                "\n" <>
                FlexibleSUSY`FSModelName <> "_mass_eigenstates model_ = model;\n" <>
                "EvaluationContext context{ model_ };\n" <>

                "// get Fermi constant from Les Houches input file\n" <>
                "const auto GF {qedqcd.displayFermiConstant()};\n" <>

                "const auto photon_exchange = calculate_" <> CXXNameOfField[inFermion] <> "_" <>
                    CXXNameOfField[outFermion] <> "_" <> CXXNameOfField[SARAH`VP] <> "_form_factors (" <>
                    If[TreeMasses`GetDimension[inFermion] =!= 1, "generationIndex1, ", " "] <>
                    If[TreeMasses`GetDimension[outFermion] =!= 1, " generationIndex2, ", " "] <>
                    "model);\n" <>

                "\n// translate from the convention of Hisano, Moroi & Tobe to Kitano, Koike & Okada\n" <>
                (* TODO: check the statement below *)
                "// Hisano defines form factors A2 through a matrix element in eq. 14\n" <>
                "// Kitano uses a lagrangian with F_munu. There is a factor of 2 from translation\n" <>
                "// because Fmunu = qeps - eps q\n" <>
                "const auto A2L = -0.5 * photon_exchange[2]/(4.*GF/sqrt(2.)) *(-1.);\n" <>
                "const auto A2R = -0.5 * photon_exchange[3]/(4.*GF/sqrt(2.)) *(-1.);\n" <>

                "\n// ------ penguins ------\n" <>
                "// 2 up and 1 down quark in proton (gp couplings)\n" <>
                "// 1 up and 2 down in neutron (gn couplings)\n" <>

                "\n// mediator: massless vectors\n" <>
                "\n// construct 4-fermion operators from A1 form factors\n" <>
                "// i q^2 A1 * (- i gmunu/q^2) * (-i Qq e) = GF/sqrt2 * gpV\n" <>
                "// VP\n" <>

                "const auto uEMVectorCurrent = vectorCurrent<typename Fu::lorentz_conjugate, Fu, typename VP::lorentz_conjugate>(model);\n" <>
                "const auto dEMVectorCurrent = vectorCurrent<typename Fd::lorentz_conjugate, Fd, typename VP::lorentz_conjugate>(model);\n\n" <>

                "auto gpLV = -sqrt(2.0)/GF * (2.*uEMVectorCurrent + dEMVectorCurrent) * photon_exchange[0] *(-1.);\n" <>
                "auto gpRV = -sqrt(2.0)/GF * (2.*uEMVectorCurrent + dEMVectorCurrent) * photon_exchange[1] *(-1.);\n" <>
                "auto gnLV = -sqrt(2.0)/GF * (uEMVectorCurrent + 2.*dEMVectorCurrent) * photon_exchange[0] *(-1.);\n" <>
                "auto gnRV = -sqrt(2.0)/GF * (uEMVectorCurrent + 2.*dEMVectorCurrent) * photon_exchange[1] *(-1.);\n" <>

                "\n// mediator: massive vectors\n" <>
                StringJoin @ Map[
                    ("\n// " <> CXXNameOfField[#] <> "\n" <>
                        "const auto " <> CXXNameOfField[#] <> "_exchange = create_massive_penguin_amp<" <> CXXNameOfField[#] <>">(" <>
                        If[TreeMasses`GetDimension[inFermion] =!= 1, "generationIndex1, ", " "] <>
                        If[TreeMasses`GetDimension[outFermion] =!= 1, " generationIndex2, ", " "] <>
                        "model, qedqcd);\n" <>
                        "gpLV += 2.*" <> CXXNameOfField[#] <> "_exchange[0] + " <> CXXNameOfField[#] <> "_exchange[2];\n" <>
                        "gpRV += 2.*" <> CXXNameOfField[#] <> "_exchange[1] + " <> CXXNameOfField[#] <> "_exchange[3];\n" <>
                        "gnLV += " <> CXXNameOfField[#] <> "_exchange[0] + 2.*" <> CXXNameOfField[#] <> "_exchange[2];\n" <>
                        "gnRV += " <> CXXNameOfField[#] <> "_exchange[1] + 2.*" <> CXXNameOfField[#] <> "_exchange[3];\n")&,

                    (* create a list of massive and electrically neutral gauge bosons *)
                    Select[GetVectorBosons[], !(IsMassless[#] || IsElectricallyCharged[#])&]
                ] <>

                (* TODO: add contributions from scalar penguins *)
                "\n// mediator: massive scalars \n\n" <>

                "gpLV += sqrt(2.0)/GF * 0.;\n" <>
                "gpRV += sqrt(2.0)/GF * 0.;\n" <>
                "gnLV += sqrt(2.0)/GF * 0.;\n" <>
                "gnRV += sqrt(2.0)/GF * 0.;\n" <>

                "\n// ------ boxes ------\n\n" <>

                "gpLV += 0.;\n" <>
                "gpRV += 0.;\n" <>
                "gnLV += 0.;\n" <>
                "gnRV += 0.;\n" <>

                "const auto nuclear_form_factors = get_overlap_integrals(flexiblesusy::" <>
                    ToString[FlexibleSUSY`FSModelName] <> "_f_to_f_conversion::Nucleus::" <> SymbolName[nucleus] <>
                    ", qedqcd" <> ");\n" <>

                "\nconst auto left {A2R*nuclear_form_factors.D + gpLV*nuclear_form_factors.Vp + gnLV*nuclear_form_factors.Vn};\n" <>
                "const auto right {A2L*nuclear_form_factors.D + gpRV*nuclear_form_factors.Vp + gnRV*nuclear_form_factors.Vn};\n" <>

                "\n// eq. 14 of Kitano, Koike and Okada\n" <>
                "return 2.*pow(GF,2)*(std::norm(left) + std::norm(right));\n"
            ] <>
            "}\n";
    
        {prototype <> ";", definition}
    ];

End[];
EndPackage[];
