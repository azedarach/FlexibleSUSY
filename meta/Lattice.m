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

BeginPackage["Lattice`", {
    "SARAH`",
    "BetaFunction`",
    "Parametrization`",
    "Traces`",
    "TreeMasses`",
    "Phases`",
    "Vertices`",
    "SelfEnergies`",
    "LatticeUtils`",
    "CConversion`",
    "TextFormatting`",
    "WriteOut`",
    "Makefile`"}]

WriteLatticeCode::usage;
ParametrizeBetaFunctions::usage;

Begin["`Private`"]

Unprotect[Real];

Real /: Format[x_Real, CForm] :=
    Format[ToString@NumberForm[x, Ceiling[MachinePrecision],
			       NumberFormat -> CDoubleFormat],
	   OutputForm];

CDoubleFormat[m_, _, ""] :=
    StringReplace[m, RegularExpression["\\.$"] -> ".0"];

CDoubleFormat[m_, _, e_] :=
    Row[{StringReplace[m, RegularExpression["\\.$"] -> ""], "e", e}];

Protect[Real];

Format[x, CForm] := Format["x", OutputForm];

Format[  x[i_], CForm] := Format[  "x["<>ToString[CForm[i]]<>"]", OutputForm];

Format[ dx[i_], CForm] := Format[ "dx["<>ToString[CForm[i]]<>"]", OutputForm];

Format[ddx[i_,j_], CForm] := Format["ddx(" <> ToString[CForm[i]] <> "," <>
				    ToString[CForm[j]] <> ")", OutputForm];

Format[a, CForm] := Format["a", OutputForm];

Format[Lattice`Private`p2, CForm] := Format["p2", OutputForm];

Format[scl0, CForm] := Format["scl0", OutputForm];

Format[scl[], CForm] := Format["scl()", OutputForm];

Format[drv[ap_, p_], CForm] :=
    Format[DrvToCFormString[drv[ap, p]], OutputForm];

DrvToCFormString[drv[ap_, p_]] :=
    "d"<>ToString[CForm[HoldForm[ap]]] <> "d"<>ToString[CForm[HoldForm[p]]];

Format[Lattice`Private`Re[z_], CForm] :=
    Format["RE" <> ToValidCSymbolString[z], OutputForm];
Format[Lattice`Private`Im[z_], CForm] :=
    Format["IM" <> ToValidCSymbolString[z], OutputForm];

Format[Lattice`Private`M[f_[{i__}]], CForm] :=
    Format[Lattice`Private`M[f][i], CForm];

Format[Lattice`Private`M2[f_[{i__}]], CForm] :=
    Format[Lattice`Private`M2[f][i], CForm];

Format[Lattice`Private`M[f_], CForm] :=
    Format["M" <> ToValidCSymbolString[f], OutputForm];

Format[Lattice`Private`M2[f_], CForm] :=
    Format["M2" <> ToValidCSymbolString[f], OutputForm];

Format[Lattice`Private`SUM, CForm] := Format["SUM", OutputForm];

Format[Lattice`Private`USUM, CForm] := Format["USUM", OutputForm];

Format[Lattice`Private`LispAnd, CForm] := Format["LispAnd", OutputForm];

Format[Lattice`Private`Complex, CForm] :=
    Format["std::complex<double>", OutputForm];

Format[InCScope[scope_, z_], CForm] :=
    Format[If[ValueQ@CContext[scope], CContext[scope], ""] <>
	   ToString[CForm[z]], OutputForm];

WriteLatticeCode[
    sarahAbbrs_List, betaFunctions_List, anomDims_List,
    fsMassMatrices_, fsNPointFunctions_,
    vertexRules_,
    phases_,
    fsEwsbEquations_,
    gaugeCouplingRules_, vevRules_, otherParameterRules_, templateRules_,
    modelName_, templateDir_, outputDir_] :=

Block[{
	drv,
	ToEnumSymbol,
	DeclaredRealQ,
	DependenceNode,
	CScope
    },
    Format[d:drv[_, _], CForm] := Format[DrvToCFormString[d], OutputForm];

Module[{
	parameterRules = Join[
	    {t -> Re[t]}, gaugeCouplingRules, vevRules, otherParameterRules],
	gaugeCouplings = RealVariables[gaugeCouplingRules],
	betaFunctionRules, betaFunctionDerivRules,
	abbrRules, abbrDerivRules,
        trivialAbbrRules, nonTrivialAbbrRules,
	parameters, enumRules, enumParameters,
	abbrDecls, abbrDefs,
	betaDecls, betaDefs,
	defChunks, nDefChunks,
	replacementFiles, cFiles,
	betaCFile, betaCFiles,
	massMatrices = fsMassMatrices /. sarahOperatorReplacementRules,
	matrixDefs, matrixStmts,
	eigenVarsDefs, eigenVarsStmts,
	dependenceNumDecls, dependenceNumDefs,
	vertexDecls, vertexDefs,
	nPointFunctions = fsNPointFunctions //. restoreMassPowerRules /.
	    FlexibleSUSY`M -> Lattice`Private`M,
	replaceGhosts =
	    SelfEnergies`Private`ReplaceGhosts[FlexibleSUSY`FSEigenstates],
	nPointDecls, nPointDefs,
	phaseDefs,
	vevs = Union @ Variables @ vevRules[[All,2]],
	treeEwsbConstraints = fsEwsbEquations /. sarahOperatorReplacementRules,
	softHiggsMasses,
	treeEwsbEquations, shiftHiggsMasses,
	ewsbConstraints, ewsbEquations, ewsbDep, ewsbList,
	fixTsusy, tsusyConstraint = (scl[])^4 - Lattice`Private`M2[Global`Su[{1}]] Lattice`Private`M2[Global`Su[{6}]] /. sarahOperatorReplacementRules,
    },
    DeclaredRealQ[a | scl0 | scl[]] := True;
    DeclaredRealQ[_] := False;
    DependenceNode[(SARAH`A0|SARAH`B0|SARAH`B1|SARAH`B00|SARAH`B22|
		    SARAH`F0|SARAH`G0|SARAH`H0)[m__]] := {Re[t], m};
    SetDependenceNode[scl[], Re[t]];
    parameters = RealVariables[parameterRules];
    enumRules = EnumRules[parameters];
    enumParameters = EnumParameters[enumRules];
    SetToEnumSymbol /@ enumRules;
    {betaFunctionRules, abbrRules} =
	ParametrizeBetaFunctions[betaFunctions, sarahAbbrs, parameterRules];
    DoneLn[
    abbrDerivRules = DifferentiateRules[abbrRules, parameters, abbrRules];
    {trivialAbbrRules, nonTrivialAbbrRules} =
	SeparateTrivialRules@Flatten[{abbrRules, abbrDerivRules}];
    {abbrDecls, abbrDefs} = CFxnsToCCode[
        AbbrRuleToC /@ nonTrivialAbbrRules, Global`$flexiblesusyCSrcChunkSize],
    "Differentiating abbreviations... "];
    betaFunctionRules = Transpose[
	ScaleByA[#, gaugeCouplings]& /@
	SortBy[betaFunctionRules,
	       (# /. {{BETA[_Integer, p_] -> _, ___}, ___} :>
		Position[parameters, p])&]];
    {betaDecls, betaDefs} = CFxnsToCCode[Flatten[
	BetaFunctionRulesToC[betaFunctionRules, enumRules, abbrRules]],
	Global`$flexiblesusyCSrcChunkSize];
    {dependenceNumDecls, dependenceNumDefs} = CFxnsToCCode[
	DepNumRuleToC /@ ParametrizeDependenceNums[
	    FindDependenceNums[massMatrices], parameterRules]];
    {matrixDefs, matrixStmts} = CMatricesToCCode[
	MatrixToC /@ Cases[parameterRules, HoldPattern[_ -> _?MatrixQ]]];
    {eigenVarDefs, eigenVarStmts} = CMatricesToCCode[
	FSMassMatrixToC /@
	ParametrizeMasses[massMatrices, parameterRules]];
    {vertexDecls, vertexDefs} = CFxnsToCCode[
	VertexRuleToC /@
	ParametrizeVertexRules[vertexRules, parameterRules]];
    {nPointDecls, nPointDefs} = CFxnsToCCode[
	NPointFunctionToC@ParametrizeNPointFunction[#, replaceGhosts]& /@
	nPointFunctions];
    phaseDefs = Phases`CreatePhasesDefinition[phases];
    replacementFiles = {
	{FileNameJoin[{templateDir, "lattice_info.hpp.in"}],
	 FileNameJoin[{outputDir, modelName <> "_lattice_info.hpp"}]},
	{FileNameJoin[{templateDir, "lattice_model.hpp.in"}],
	 FileNameJoin[{outputDir, modelName <> "_lattice_model.hpp"}]},
	{FileNameJoin[{templateDir, "lattice_model.cpp.in"}],
	 FileNameJoin[{outputDir, modelName <> "_lattice_model.cpp"}]},
	{FileNameJoin[{templateDir, "lattice_model_interactions.cpp.in"}],
	 FileNameJoin[{outputDir, modelName <> "_lattice_model_interactions.cpp"}]}};
(* if EWSB is demanded *)
    softHiggsMasses = SoftHiggsMasses[treeEwsbConstraints];
    treeEwsbEquations = ParametrizeEWSBEquations[
	treeEwsbConstraints, ConditionPositiveVevs[vevs], softHiggsMasses,
	parameterRules];
    shiftHiggsMasses = CShiftHiggsMasses[treeEwsbEquations];
    ewsbConstraints = EWSBConstraintsWithCorrections[treeEwsbConstraints];
    ewsbEquations = ParametrizeEWSBEquations[
	ewsbConstraints, ConditionPositiveVevs[vevs], softHiggsMasses,
	parameterRules];
    {ewsbDep, ewsbList} = CNConstraintsToCCode @
			  EWSBConditionsToC[ewsbEquations];
    ewsbDep = ToString[ewsbDep];
    fixTsusy = CNConstraintToCCode[NConstraintToC[tsusyConstraint /. parameterRules]];
    replacementFiles = Join[replacementFiles, {
	{FileNameJoin[{templateDir, "lattice_susy_scale_constraint.hpp.in"}],
	 FileNameJoin[{outputDir, modelName <> "_lattice_susy_scale_constraint.hpp"}]},
	{FileNameJoin[{templateDir, "lattice_susy_scale_constraint.cpp.in"}],
	 FileNameJoin[{outputDir, modelName <> "_lattice_susy_scale_constraint.cpp"}]},
	{FileNameJoin[{templateDir, "lattice_tsusy_constraint.hpp.in"}],
	 FileNameJoin[{outputDir, modelName <> "_lattice_tsusy_constraint.hpp"}]},
	{FileNameJoin[{templateDir, "lattice_tsusy_constraint.cpp.in"}],
	 FileNameJoin[{outputDir, modelName <> "_lattice_tsusy_constraint.cpp"}]},
	{FileNameJoin[{templateDir, "lattice_ewsb_constraint.hpp.in"}],
	 FileNameJoin[{outputDir, modelName <> "_lattice_ewsb_constraint.hpp"}]},
	{FileNameJoin[{templateDir, "lattice_ewsb_constraint.cpp.in"}],
	 FileNameJoin[{outputDir, modelName <> "_lattice_ewsb_constraint.cpp"}]}}];
(* end if *)
    WriteOut`ReplaceInFiles[replacementFiles,
	Join[templateRules, {
	"@abbrDecls@"	    -> IndentText[abbrDecls, 2],
	"@betaDecls@"	    -> IndentText[betaDecls, 2],
	"@dependenceNumDecls@" -> IndentText[dependenceNumDecls, 4],
	"@dependenceNumDefs@"  -> WrapText[StringJoin@dependenceNumDefs],
	"@eigenVarDefs@"    -> IndentText[eigenVarDefs, 4],
	"@eigenVarStmts@"   -> WrapText[eigenVarStmts],
	"@enumParameters@"  -> WrapText@IndentText[enumParameters, 2],
	"@nRowsEwsb@"	    -> ToString@Length[ewsbEquations],
	"@ewsbDep@"	    -> StringTrim@WrapText@IndentText[ewsbDep, 4],
	"@ewsbList@"	    -> StringTrim@WrapText@IndentText[ewsbList, 4],
	"@fixTsusy@"        -> StringTrim@WrapText@IndentText[fixTsusy, 2],
	"@matrixDefs@"	    -> IndentText[matrixDefs, 4],
	"@matrixStmts@"	    -> WrapText[matrixStmts],
	"@nPointDecls@"	    -> IndentText[nPointDecls, 4],
	"@nPointDefs@"	    -> WrapText[StringJoin@nPointDefs],
	"@phaseDefs@"	    -> IndentText[phaseDefs, 4],
	"@shiftHiggsMasses@"-> WrapText@IndentText[shiftHiggsMasses, 2],
	"@vertexDecls@"	    -> IndentText[vertexDecls, 4],
	"@vertexDefs@"	    -> WrapText[StringJoin@vertexDefs]
    }]];
    defChunks = Join[abbrDefs, betaDefs];
    nDefChunks = Length[defChunks];
    betaCFiles = MapIndexed[(
	betaCFile = FileNameJoin[
	    {outputDir, modelName <> "_lattice_model_betafunctions_" <>
	     IntegerString[First[#2], 10, StringLength@ToString[nDefChunks]] <>
	     ".cpp"}];
	WriteOut`ReplaceInFiles[{
	    {FileNameJoin[{templateDir, "lattice_model_betafunctions.cpp.in"}],
	     betaCFile}},
	    Join[templateRules, {
		"@abbrBetaDefs@"    -> WrapText[#1]
	    }]];
	betaCFile)&,
    defChunks];
    cFiles = Join[
	Select[replacementFiles[[All,2]], StringMatchQ[#, ___ ~~ ".cpp"]&],
	betaCFiles];
    WriteMakefile[templateDir, outputDir, cFiles, templateRules]
]];

ParametrizeEWSBEquations[
    constraints_List, assumptions_List, outputVars_, parameterRules_] :=
Module[{
	pConstraints = constraints /. parameterRules,
	pAssumptions = assumptions /. parameterRules,
	vars = Union@Variables[outputVars /. parameterRules],
	eqs
    },
    eqs = Refine[Reduce[Join[# == 0& /@ pConstraints, pAssumptions], vars],
		 pAssumptions];
    If[MatchQ[eqs, HoldPattern[And[(_?(MemberQ[vars, #]&) == _)...]]],
       List @@ eqs,
       Print["Lattice`ParametrizeEWSBEquations[] failed to understand the result from Reduce[]: ", eqs];
       Abort[]]
];

EWSBConditionsToC[pEquations_List] := EWSBConditionToC /@ pEquations;

EWSBConditionToC[lhs_ == rhs_] := NConstraintToC[rhs - lhs];

CShiftHiggsMasses[treeEwsbEquations_List] :=
    StringJoin[CShiftHiggsMass /@ treeEwsbEquations];

CShiftHiggsMass[lhs_ == rhs_] := Module[{
	clhs = CExpToCFormString @ ToCExp[lhs, x],
	crhs = CExpToCFormString @ ToCExp[rhs, x]
    },
    clhs <> " = " <> crhs <> ";\n"
];

NConstraintToC[constraint_] :=
    CNConstraint[
	Dependence -> DependenceList[constraint],
	Expr -> ToCExp[constraint, x]
    ];

DependenceList[expr_] :=
    SortEnums[ToEnumSymbol@PrivatizeReIm@ToCExp[#]& /@ FindDependence[expr]];

SortEnums[enums_List] :=
    SortBy[enums,
	   ToExpression @ StringReplace[
	       ToString[#], RegularExpression["^l([[:digit:]]+).*$"] :> "$1"] &
    ];

EWSBConstraintsWithCorrections[treeConstraints_List] := Module[{
	tadpoles = (# /. {field_, index_} :>
		    Symbol["Tadpole"<>ToString[field]][index-1])& /@
		   FlexibleSUSY`Private`CreateHiggsToEWSBEqAssociation[]
    },
    treeConstraints - a tadpoles
];

SoftHiggsMasses[ewsbConstraints_] := Union @
    Select[SARAH`ListSoftBreakingScalarMasses, !FreeQ[ewsbConstraints, #]&];

ConditionPositiveVevs[vevs_List] := # > 0& /@ vevs;

CNConstraintsToCCode[cncs_List] := Module[{
	dependences, dependence,
	exprs, csl
    },
    {dependences, exprs} = Transpose[({Dependence, Expr} /. List@@#)& /@ cncs];
    dependence = SortEnums[Union @@ dependences];
    csl = StringJoin @
	  Riffle[Block[{CContext},
		       CContext["CLASSNAME::Interactions::"] = "I.";
		       CExpToCFormString @ ReCExp[#]]& /@ exprs,
		 ",\n"];
    {dependence, csl}
];

CNConstraintToCCode[cnc_] := Module[{
	dependence, expr
    },
    {dependence, expr} = {Dependence, Expr} /. List@@cnc;
    StringJoin[
	"AnyNumericalConstraint(\n",
	"  [&](const AnyNumericalConstraint *self, const double *x) {\n",
	"    double a    = self->f->a;\n",
	"    double scl0 = self->f->scl0;\n",
	"    CLASSNAME::Interactions I;\n",
	"    I.set(Eigen::Map<const Eigen::VectorXd>(x, eftWidth), scl0);\n",
	"    return ",
	Block[{CContext},
	    CContext["CLASSNAME::Interactions::"] = "I.";
	    CExpToCFormString @ ReCExp[expr]], ";\n",
	"  },\n",
	"  ", ToString[dependence], "\n",
	")"]
];

restoreMassPowerRules = {
    (f:SARAH`A0|SARAH`B0|SARAH`B1|SARAH`B00|SARAH`B22|
     SARAH`F0|SARAH`G0|SARAH`H0)[a___,FlexibleSUSY`M[b_],c___] :>
    f[a,Lattice`Private`M2[b],c]
};

ParametrizeNPointFunction[h_[field_, expr_], replaceGhosts_] :=
    h[field,
      expr /. p^2 -> Lattice`Private`p2 /. C -> 1 /.
      sarahOperatorReplacementRules /.
      cp:_SARAH`Cp|_SARAH`Cp[_] :> CVertexFunction[cp] /. replaceGhosts //.
      matrixOpRules];

ParametrizeVertexRules[vertexRules_, parameterRules_] := Module[{
	scalarParameterRules =
	    DeleteCases[parameterRules, HoldPattern[_ -> _?MatrixQ]]
    },
    ParametrizeVertexRule[#, scalarParameterRules]& /@ vertexRules
];

ParametrizeVertexRule[lhs_ -> rhs_, parameterRules_] :=
    (lhs /. sarahOperatorReplacementRules) ->
    (rhs /. sarahOperatorReplacementRules /. parameterRules //. matrixOpRules);

ParametrizeDependenceNums[depNums_List, parameterRules_] :=
    ParametrizeDependenceNum[#, parameterRules]& /@ depNums;

ParametrizeDependenceNum[lhs_ -> rhs_, parameterRules_] :=
    lhs -> (rhs /. parameterRules);

ParametrizeMasses[massMatrices_, parameterRules_] := Module[{
	matrices = Cases[parameterRules, HoldPattern[_ -> _?MatrixQ]][[All,1]]
    },
    ReleaseHold[Hold[massMatrices] /.
	(m_[i__Integer] /; MemberQ[matrices, m]) :> m[[i]] /.
	parameterRules //. matrixOpRules]
];

CMatricesToCCode[cmatrices_] :=
    StringJoin /@ Transpose[
	({{EigenDef, "\n"}, {SetStmt, "\n"}} /. List@@#)& /@
	Flatten[cmatrices]];

FSMassMatrixToC[FSMassMatrix[{m_?RealQ}, f_, _]] :=
    MassToC[m, f, "double"];

FSMassMatrixToC[FSMassMatrix[{m_}, f_, _]] :=
    MassToC[m, f, "std::complex<double>"];

FSMassMatrixToC[FSMassMatrix[m_?RealMatrixQ, f_, {u_Symbol, v_Symbol}]] :=
    SVDToC[m, f, u, v, "double"];

FSMassMatrixToC[FSMassMatrix[m_?MatrixQ, f_, {u_Symbol, v_Symbol}]] :=
    SVDToC[m, f, u, v, "std::complex<double>"];

FSMassMatrixToC[FSMassMatrix[m_?RealMatrixQ, f_, z_Symbol]] :=
    HermitianToC[m, f, z, "double"];

FSMassMatrixToC[FSMassMatrix[m_?MatrixQ, f_, z_Symbol]] :=
    HermitianToC[m, f, z, "std::complex<double>"];

MatrixToC[symbol_ -> m_?RealMatrixQ] :=
    MatrixToC[m, symbol, "double"];

MatrixToC[symbol_ -> m_?MatrixQ] :=
    MatrixToC[m, symbol, "std::complex<double>"];

MatrixToC[m_, symbol_, scalarType_] := Module[{
	d1, d2, name = CExpToCFormString@ToCExp[symbol]
    },
    SetDependenceNode[symbol[_,_], m];
    {d1, d2} = Dimensions[m];
    If[scalarType === "double", DeclaredRealQ[symbol[_,_]] := True];
    CMatrix[
	EigenDef ->
	    CEigenMatrixType[scalarType, d1, d2] <> " " <> name <> ";",
	SetStmt ->
	    CSetMatrix[m, scalarType === "double", name]
    ]
];

MassToC[m_, f_, cType_] := Module[{
	ev = ToCMassName[f]
    },
    SetDependenceNode[MassN[f], m];
    (* DeclaredRealQ[ev] := True by pattern matching *)
    CMatrix[
	EigenDef -> cType <> " " <> ev <> ";",
	SetStmt -> "  " <> ev <> " = " <> CExpToCFormString@ToCExp[m, x] <> ";"
    ]
];

SVDToC[m_, f_, IdentityMatrix, IdentityMatrix, scalarType_] := Module[{
	d1, d2, ds, ev
    },
    SetDependenceNode[MassN[f[{__}]], m];
    ds = Min[{d1, d2} = Dimensions[m]];
    ev = ToCMassName[f];
    CMatrix[
	EigenDef -> CEigenArrayType[ds] <> " " <> ev <> ";",
	SetStmt  -> "  " <> ev <> " = " <> ToString@Diagonal[m] <> ";"
    ]
];

SVDToC[m_, f_, u_, u_, scalarType_] := Module[{
	d, ev
    },
    SetDependenceNode[MassN[f[{__}]] | u[_,_], m];
    {d, d} = Dimensions[m];
    ev = ToCMassName[f];
    CMatrix[
	EigenDef ->
	    CEigenArrayType[d] <> " " <> ev <> ";\n" <>
	    CEigenMatrixType["std::complex<double>", d, d] <> " " <>
	    ToValidCSymbolString[u] <> ";",
	SetStmt ->
	    "  {\n" <>
	    CDefTmpMatrix[m, scalarType, "tmpMat"] <> "\n" <>
	    "    fs_diagonalize_symmetric(tmpMat, " <> ev <> ", " <>
		ToValidCSymbolString[u] <> ");\n" <>
	    "  }"
    ]
];

SVDToC[m_, f_, u_, v_, scalarType_] := Module[{
	d1, d2, ds, ev
    },
    SetDependenceNode[MassN[f[{__}]] | (u|v)[_,_], m];
    ds = Min[{d1, d2} = Dimensions[m]];
    ev = ToCMassName[f];
    If[scalarType === "double", DeclaredRealQ[(u|v)[_,_]] := True];
    CMatrix[
	EigenDef ->
	    CEigenArrayType[ds] <> " " <> ev <> ";\n" <>
	    CEigenMatrixType[scalarType, d1, d1] <> " " <>
	    ToValidCSymbolString[u] <> ";\n" <>
	    CEigenMatrixType[scalarType, d2, d2] <> " " <>
	    ToValidCSymbolString[v] <> ";",
	SetStmt ->
	    "  {\n" <>
	    CDefTmpMatrix[m, scalarType, "tmpMat"] <> "\n" <>
	    "    fs_svd(tmpMat, " <> ev <> ", " <>
		ToValidCSymbolString[u] <> ", " <> ToValidCSymbolString[v] <>
		");\n" <>
	    "  }"
    ]
];

HermitianToC[m_, f_, IdentityMatrix, scalarType_] := Module[{
	d, ev
    },
    SetDependenceNode[MassN[f[{__}]], m];
    {d, d} = Dimensions[m];
    ev = ToCMassName[f];
    CMatrix[
	EigenDef -> CEigenArrayType[d] <> " " <> ev <> ";",
	SetStmt  -> "  " <> ev <> " = " <> ToString@Diagonal[m] <> ";"
    ]
];

HermitianToC[m_, f_, z_, scalarType_] := Module[{
	d, ev
    },
    SetDependenceNode[MassN[f[{__}]] | z[_,_], m];
    {d, d} = Dimensions[m];
    ev = ToCMassName[f];
    If[scalarType === "double", DeclaredRealQ[z[_,_]] := True];
    CMatrix[
	EigenDef ->
	    CEigenArrayType[d] <> " " <> ev <> ";\n" <>
	    CEigenMatrixType[scalarType, d, d] <> " " <>
	    ToValidCSymbolString[z] <> ";",
	SetStmt ->
	    "  {\n" <>
	    CDefTmpMatrix[m, scalarType, "tmpMat"] <> "\n" <>
	    "    fs_diagonalize_hermitian(tmpMat, " <> ev <> ", " <>
		ToValidCSymbolString[z] <> ");\n" <>
	    "  }"
    ]
];

CEigenMatrixType[scalarType_String, d1_Integer, d2_Integer] :=
    "Eigen::Matrix<" <> scalarType <> ", " <>
    ToString[d1] <> ", " <> ToString[d2] <> ">";

CEigenArrayType[scalarType_String, len_Integer] :=
    "Eigen::Array<" <> scalarType <> ", " <> ToString[len] <> ", 1>";

CEigenArrayType[len_Integer] := CEigenArrayType["double", len];

ToCMassName[field_Symbol] := ToString @ CForm @ MassN[field];

MassN[field_ /; FieldMassDimension[field] === 3/2] := Lattice`Private`M[field];

MassN[field_ /; FieldMassDimension[field] === 1] := Lattice`Private`M2[field];

Lattice`Private`M [field_ /; FieldMassDimension[field] === 1  ] :=
    Sqrt[Lattice`Private`M2[field]];

Lattice`Private`M2[field_ /; FieldMassDimension[field] === 3/2] :=
    AbsSqr @ Lattice`Private`M[field];

FieldMassDimension[_?IsFermion|_?IsFermion[{__}]] := 3/2;

FieldMassDimension[_] := 1;

CDefTmpMatrix[m_, scalarType_, name_] := Module[{
	d1, d2
    },
    {d1, d2} = Dimensions[m];
    "    " <> CEigenMatrixType[scalarType, d1, d2] <> " " <> name <> ";\n" <>
    CSetMatrix[m, scalarType === "double", name]
];

CSetMatrix[m_, takeReal_, name_] := Module[{
	d1, d2, i1, i2,
	re = If[takeReal, ReCExp, Identity]
    },
    {d1, d2} = Dimensions[m];
    "    " <> name <> " <<\n" <>
    Riffle[Table[
	Riffle[Table[{
	    "    /*", ToString[i1 - 1], ",", ToString[i2 - 1], "*/ ",
	    CExpToCFormString @ re @ ToCExp[m[[i1,i2]], x]}, {i2, d2}], ",\n"],
	{i1, d1}], ",\n"] <> ";"
];


RealQ[z_] := PossibleZeroQ[ConjugateExpand[z - Conjugate[z]]];

RealMatrixQ[m_?MatrixQ] := And @@ (RealQ /@ Flatten[m]);

RealMatrixQ[_] := False;

conjugateExpandDispatch = Dispatch[{
    (* Schwarz reflection principle? *)
    Conjugate[z:(_Plus|_Times|
		 _Power|
		 _Sin|_ArcSin|
		 _Cos|_ArcCos|
		 _Tan|_ArcTan)] :> Conjugate /@ z,
    Conjugate[SARAH`sum[a_, b_, c_, z_]] :> SARAH`sum[a, b, c, Conjugate[z]],
    Conjugate[z:(_SARAH`Delta|_SARAH`ThetaStep|_SARAH`Mass|_SARAH`Mass2)] :> z,
    (* take real part of a loop function *)
    Conjugate[z:(_SARAH`A0|_SARAH`B0|_SARAH`B1|_SARAH`B00|_SARAH`B22|
		 _SARAH`F0|_SARAH`G0|_SARAH`H0)] :> z,
    Conjugate[z_?DeclaredRealQ] :> z
}];

ConjugateExpand[z_] := z //. conjugateExpandDispatch;

WriteMakefile[templateDir_, outputDir_, cppFiles_, templateRules_] :=
    Makefile`ReplaceInMakefiles[{
	{FileNameJoin[{templateDir, "lattice.mk.in"}],
	 FileNameJoin[{outputDir, "lattice.mk"}]}},
	cppFiles, templateRules];

EnumRules[parameters_List] := MapIndexed[
    #1 -> "l" <> ToString[First[#2] - 1] <>
	  ToString[PrivatizeReIm@ToCExp[#1], CForm]&,
    parameters];

EnumParameters[enumRules_List] :=
    StringJoin["enum : size_t { ", {Last[#], ", "}& /@ enumRules,
	       "eftWidth };"];

SetToEnumSymbol[parameter_ -> enum_String] :=
    ToEnumSymbol[PrivatizeReIm@ToCExp[parameter]] = Symbol[enum];

PrivatizeReIm[cexp_] := cexp //. {
    Re[z_] :> Lattice`Private`Re[z],
    Im[z_] :> Lattice`Private`Im[z]
};

PrivatizeParameterReIm[cexp_] := cexp //. {
    Re[z_] /; ValueQ@ToEnumSymbol@Lattice`Private`Re@z :> Lattice`Private`Re@z,
    Im[z_] /; ValueQ@ToEnumSymbol@Lattice`Private`Im@z :> Lattice`Private`Im@z
};

ScaleByA[b:{{(BETA[1, _] -> _)..}, {(BETA[_Integer, _] -> _)...}...},
	 gaugeCouplings_] :=
    Map[ScaleByA[#, gaugeCouplings]&, b, {2}]

ScaleByA[b:BETA[1, p_] -> rhs_, gaugeCouplings_] := (b -> rhs) /;
    MemberQ[gaugeCouplings, p];

ScaleByA[b:BETA[1, _] -> rhs_, _] := b -> a rhs;

ScaleByA[b:BETA[_Integer, _] -> rhs_, _] := b -> a rhs;

RealVariables[parameterRules_] := Module[{
	rhsList = Flatten@parameterRules[[All,2]],
	rvs
    },
    rvs = Variables[rhsList];
    Assert[Length[rvs] === Length@Union[rvs]];
    SortBy[rvs, Position[rhsList, #]&]
];

DifferentiateRules[rules_, parameters_, abbrRules_] :=
    SelfApplyTrivialRules @
    Flatten[DifferentiateRule[#, parameters, abbrRules]& /@ rules]

DifferentiateRule[lhs_ -> rhs_, parameters_, abbrRules_] :=
    (Differentiate[lhs, #, abbrRules] -> Differentiate[rhs, #, abbrRules])& /@
    parameters;

SelfApplyTrivialRules[rules_] := FixedPoint[RewriteTrivialRules, rules];

RewriteTrivialRules[rules_] := Module[{
	trivialRules, nonTrivialRules
    },
    {trivialRules, nonTrivialRules} = SeparateTrivialRules[rules];
    SetRule /@ trivialRules;
    nonTrivialRules
];

SetRule[lhs_ -> rhs_] := lhs = rhs;

SeparateTrivialRules[rules_] := Flatten /@
    Last@Reap[Sow[#, NumericQ@Expand@Last[#]]& /@ rules, {True, False}];

CFxnsToCCode[cfxns_, chunkSize_:Infinity] := Module[{
	decls, defs
    },
    {decls, defs} = Transpose[CFxnToCCode /@ Flatten[cfxns]];
    {StringJoin[decls], StringGroup[StringJoin /@ defs, chunkSize, "\n"]}
];

CFxnToCCode[cfxn_] := Module[{
	returnType, scope, name, args, qualifier, attributes, body,
	argList,
	argsInDecl, argsInDef
    },
    {returnType, scope, name, args, qualifier, attributes, body} =
	{ReturnType, Scope, Name, Args,
	 {" ",Qualifier}, {" ATTR(",Attributes,")"}, Body} /.
	List@@cfxn /. {Scope -> "CLASSNAME::", {___, Qualifier} -> {},
		       {___, Attributes, ___} -> {}};
    SetCFxnScope[name, args, scope];
    argsInDecl = "(" <> Riffle[Riffle[#, " "]& /@ args, ", "] <> ")";
    argsInDef = "(" <> Riffle[
	Riffle[If[Length[#] > 1 &&
	    StringFreeQ[body, RegularExpression["\\b" <> Last[#] <> "\\b"]],
	    Most[#], #], " "]& /@ args,
	", "] <> ")";
    {{returnType, " ", name, argsInDecl, qualifier, attributes, ";\n"},
     {returnType, " ", scope, name, argsInDef, qualifier, "\n", body, "\n"}}
];

StringGroup[strings_List, chunkSize_, separator_:""] := Module[{
	size = 0
    },
    StringJoin@Riffle[#, "\n"]& /@ Last@Reap[
	(Sow[#, Quotient[size, chunkSize]]; size += StringLength[#])& /@
	strings]
];

SetCFxnScope[
    name_String, args_List, scope_String:"CLASSNAME::Interactions::"] :=
Module[{
	symbol = Quiet[Check[Symbol[name], Undefined, {Symbol::symname}],
		       {Symbol::symname}],
	pattern
    },
    If[symbol =!= Undefined,
       pattern = symbol @@ Table[_, {Length[args]}];
       CScope[pattern] = scope]
];

NPointFunctionToC[nPointFunction:_[_, rhs_]] := Module[{
	cType, re,
	name, args
    },
    {cType, re} = If[RealQ[rhs], {"double", ReCExp},
				 {"std::complex<double>", Identity}];
    name = CNPointFunctionName[nPointFunction];
    args = CNPointFunctionArgs[nPointFunction];
    SetDependenceNode[Symbol[name] @@ Table[_, {Length[args]}], rhs];
    CFxn[
	ReturnType -> cType,
	Scope -> "CLASSNAME::Interactions::",
	Name -> name,
	Args -> args,
	Qualifier -> "const",
	Attributes -> "pure",
	Body -> "{\n" <>
	"  return " <> CExpToCFormString[
	    CConversion`oneOver16PiSqr re[ToCExp[rhs, x]]] <> ";\n" <>
	"}\n"
    ]
];

CNPointFunctionName[
    (head:SelfEnergies`FSSelfEnergy|SelfEnergies`FSHeavySelfEnergy|
     SelfEnergies`FSHeavyRotatedSelfEnergy|SelfEnergies`Tadpole)
    [fieldSpec_[lorentzTag:SARAH`PL|SARAH`PR|_Integer], expr_]] :=
    CNPointFunctionName[head[fieldSpec, expr]] <>
    ToValidCSymbolString[lorentzTag];

CNPointFunctionName[
    (head:SelfEnergies`FSSelfEnergy|SelfEnergies`FSHeavySelfEnergy|
     SelfEnergies`FSHeavyRotatedSelfEnergy|SelfEnergies`Tadpole)
    [field_[indices__]|field_Symbol, expr_]] :=
    ToValidCSymbolString[head[field]];

CNPointFunctionArgs[
    (head:SelfEnergies`FSSelfEnergy|SelfEnergies`FSHeavySelfEnergy|
     SelfEnergies`FSHeavyRotatedSelfEnergy|SelfEnergies`Tadpole)
    [fieldSpec_[lorentzTag:SARAH`PL|SARAH`PR|_Integer], expr_]] :=
    CNPointFunctionArgs[head[fieldSpec, expr]];

CNPointFunctionArgs[
    (head:SelfEnergies`FSSelfEnergy|SelfEnergies`FSHeavySelfEnergy|
     SelfEnergies`FSHeavyRotatedSelfEnergy|SelfEnergies`Tadpole)
    [field_Symbol, expr_]] := CNPointFunctionArgs[head[field[], expr]];

CNPointFunctionArgs[
    (head:SelfEnergies`FSSelfEnergy|SelfEnergies`FSHeavySelfEnergy|
     SelfEnergies`FSHeavyRotatedSelfEnergy)
    [field_[indices___], expr_]] :=
    Prepend[{"size_t", ToString[#]}& /@ {indices}, {"double", "p2"}];

CNPointFunctionArgs[
    SelfEnergies`Tadpole
    [field_[indices___], expr_]] :=
    {"size_t", ToString[#]}& /@ {indices};

VertexRuleToC[lhs_ -> rhs_] := Module[{
	cType, re,
	name, args
    },
    name = CVertexFunctionName[lhs];
    args = CVertexFunctionArgs[lhs];
    SetDependenceNode[Symbol[name] @@ Table[_, {Length[args]}], rhs];
    {cType, re} = If[RealQ[rhs],
		     DeclaredRealQ[CVertexFunction[lhs]] := True;
		     {"double", ReCExp},
		     {"std::complex<double>", Identity}];
    CFxn[
	ReturnType -> cType,
	Scope -> "CLASSNAME::Interactions::",
	Name -> name,
	Args -> args,
	Qualifier -> "const",
	Attributes -> "pure",
	Body -> "{\n" <>
	"  return " <> CExpToCFormString @ re @ ToCExp[rhs, x] <> ";\n" <>
	"}\n"
    ]
];

CVertexFunction[cp_] :=
    Symbol[CVertexFunctionName[cp]] @@
    (DecInt /@ Flatten[
	FieldIndexList /@ GetParticleList[cp]]);

CVertexFunctionName[cpPattern_] := Module[{
	fields = GetParticleList[cpPattern]
    },
    ToValidCSymbolString[
	cpPattern /. Thread[fields -> ((# /. h_[_List] :> h)& /@ fields)]]
];

CVertexFunctionArgs[cpPattern_] :=
    Replace[#, {
		index_Pattern :> {"size_t", ToString@First[index]},
		unused_	      :> {"size_t"}
	       }]& /@
    Flatten[FieldIndexList /@ GetParticleList[cpPattern]];

DepNumRuleToC[lhs_ -> rhs_] := (
    SetDependenceNode[lhs[], rhs];
    CFxn[
	ReturnType -> If[RealQ[rhs], DeclaredRealQ[lhs[]] := True; "double",
			 "std::complex<double>"],
	Scope -> "CLASSNAME::Interactions::",
	Name -> CExpToCFormString@ToCExp[lhs],
	Args -> {},
	Qualifier -> "const",
	Attributes -> "pure",
	Body -> "{\n" <>
	"  return " <> CExpToCFormString@ToCExp[rhs, x] <> ";\n" <>
	"}\n"
    ]
);

AbbrRuleToC[lhs_ -> rhs_] :=
CFxn[
    ReturnType -> "double",
    Name -> CExpToCFormString@ToCExp[lhs],
    Args -> {{"const Eigen::VectorXd&", "x"}},
    Qualifier -> "const",
    Attributes -> "pure",
    Body -> "{\n" <>
    "  return " <> CExpToCFormString@ToCExp[rhs, x] <> ";\n" <>
    "}\n"
];

BetaFunctionRulesToC[betanLRules_, enumRules_, abbrRules_] := Flatten[{
Reap[
    CFxn[
	ReturnType -> "void",
	Name -> "dx",
	Args -> {{"double", "a"}, {"const Eigen::VectorXd&", "x"},
		 {"Eigen::VectorXd&", "dx"}, {"size_t", "nloops"}},
	Qualifier -> "const",
	Body -> StringJoin["{\n",
	"  dx.setZero();\n",
	"  dx[l0REt] = 1;\n",
	MapIndexed[{
	"\n  if (nloops < ", ToString@First[#2], ") return;\n",
	Module[{flattened = Flatten[#1], nRules, name},
	nRules = Length[flattened];
	MapIndexed[(
	WriteString["stdout",
		    "[",First[#2],"/",nRules,"] ","BETA"@@First[#1],":"];
	name = "d" <> CExpToCFormString@ToCExp@Last@First[#1] <> "_" <>
	       ToString@First@First[#1] <> "loop";
	Sow[CFxn[
	    ReturnType -> "void",
	    Name -> name,
	    Args -> {{"double", "a"}, {"const Eigen::VectorXd&", "x"},
		     {"Eigen::VectorXd&", "dx"}},
	    Qualifier -> "const",
	    Body -> StringJoin["{\n",
	    "  ", BetaFunctionRuleToCStmt[#1],
	    "}\n"]]
	];
	{"  ", name, "(a, x, dx);\n"})&,
	flattened]]}&, betanLRules],
	"}\n"]
    ]
],
Reap[
    CFxn[
	ReturnType -> "void",
	Name -> "ddx",
	Args -> {{"double", "a"}, {"const Eigen::VectorXd&", "x"},
		 {"Eigen::MatrixXd&", "ddx"}, {"size_t", "nloops"}},
	Qualifier -> "const",
	Body -> StringJoin["{\n",
	"  ddx.setZero();\n",
	MapIndexed[{
	"\n  if (nloops < ", ToString@First[#2], ") return;\n",
	Module[{flattened = Flatten[#1], nRules, name},
	nRules = Length[flattened];
	MapIndexed[(
	name = "dd" <> CExpToCFormString@ToCExp@Last@First[#1] <> "_" <>
	ToString@First@First[#1] <> "loop";
	WriteString["stdout",
		    "[",First[#2],"/",nRules,"] ", "D[BETA"@@First[#1], "]:"];
	Sow[Module[{
		reaped = Reap[
		    BetaFunctionRuleToDerivCStmt[#1, enumRules, abbrRules]]
	    },
	    WriteString[
		"stdout", " diff... ", FormatShortTime[
		Plus@@Cases[Flatten@Last[reaped], tagDiff[t_] :> t]],
		" to CExp... ", FormatShortTime[
		Plus@@Cases[Flatten@Last[reaped], tagCExp[t_] :> t]],
		" to C... ", FormatShortTime[
		Plus@@Cases[Flatten@Last[reaped], tagC   [t_] :> t]], "\n"];
	    CFxn[
		ReturnType -> "void",
		Name -> name,
		Args -> {{"double", "a"}, {"const Eigen::VectorXd&", "x"},
			 {"Eigen::MatrixXd&", "ddx"}},
		Qualifier -> "const",
		Body -> StringJoin["{\n", First[reaped], "}\n"]]]];
	{"  ", name, "(a, x, ddx);\n"})&,
	flattened]]}&, betanLRules],
	"}\n"]
    ]
]}];

BetaFunctionRuleToCStmt[BETA[1, p:(Re|Im)[_]] -> rhs_] :=
    BetaFunctionRuleToAssignment[1, p, rhs, "="];

BetaFunctionRuleToCStmt[BETA[level_Integer, p:(Re|Im)[_]] -> rhs_] :=
    BetaFunctionRuleToAssignment[level, p, rhs, "+="];

BetaFunctionRuleToAssignment[_Integer, _, rhs_, _, _, _] := {} /;
    Expand[rhs] === 0;

BetaFunctionRuleToAssignment[level_Integer, p_, rhs_, op_] :=
    CExpToCFormString[ToCExp[p, dx]] <> " " <> op <> " " <>
    CExpToCFormString[CConversion`oneOver16PiSqr^level
			ToCExp[rhs, x]] <> ";\n";

BetaFunctionRuleToAssignment[level_Integer, p_, rhs_, op_] := Module[{
	rhsCExp = Done[ToCExp[rhs, x], " to CExp... "]
    },
    DoneLn[CExpToCFormString[ToCExp[p, dx]] <> " " <> op <> " " <>
	   CExpToCFormString[CConversion`oneOver16PiSqr^level rhsCExp] <>
	   ";\n",
	   " to C... "]
];

BetaFunctionRuleToDerivCStmt[BETA[1, p:(Re|Im)[_]] -> rhs_,
			     enumRules_, abbrRules_] :=
    BetaFunctionRuleToAssignments[1, p, rhs,
				  enumRules, abbrRules, "="];

BetaFunctionRuleToDerivCStmt[BETA[level_Integer, p:(Re|Im)[_]] -> rhs_,
			     enumRules_, abbrRules_] :=
    BetaFunctionRuleToAssignments[level, p, rhs,
				  enumRules, abbrRules, "+="];

BetaFunctionRuleToAssignments[
    level_Integer, p_, rhs_, enumRules_, abbrRules_, op_] :=
Module[{
	pidx = p /. enumRules
    },
    Module[{
	    q, qidx,
	    deriv, derivCExp, derivC
	},
	{q, qidx} = List @@ #;
	Sow[tagDiff@First@Timing[deriv = Differentiate[rhs, q, abbrRules]]];
	If[Expand[deriv] === 0, {},
	   Sow[tagCExp@First@Timing[derivCExp = ToCExp[deriv, x]]];
	   Sow[tagC@First@Timing[derivC =
	       {"  ddx(", qidx, ",", pidx, ") ", op, " ",
		CExpToCFormString[CConversion`oneOver16PiSqr^level derivCExp],
		";\n"}]];
	   derivC]
    ]& /@ enumRules
];

toCExpDispatch = Dispatch[{
    Re[abbr_Symbol] /; !ValueQ@ToEnumSymbol@Lattice`Private`Re@abbr :>
	Lattice`Private`Re[abbr],
    Im[abbr_Symbol] /; !ValueQ@ToEnumSymbol@Lattice`Private`Im@abbr :>
	Lattice`Private`Im[abbr],
    Conjugate[z_?CRealTypeQ] :> z,
    Conjugate[z_] z_ :> AbsSqr[z],
    SARAH`sum[i_, a_, b_, x_] :> Lattice`Private`SUM[i, a, b, x],
    (* SARAH`Delta[0, _] := 0
       following the principle of greatest astonishment *)
    SARAH`Delta[i__] :> KroneckerDelta[i],
    HoldPattern@Times[a___, SARAH`ThetaStep[i_, j_], b___] :>
	Lattice`Private`ThetaStep[i, j, Times[a, b]],
    HoldPattern[SARAH`Mass [f_]] :> Lattice`Private`M [f],
    HoldPattern[SARAH`Mass2[f_]] :> Lattice`Private`M2[f],
    Complex[p__] :> Lattice`Private`Complex[p]
}];

replaceIndicesDispatch = Dispatch[{
    m:(_KroneckerDelta|_SARAH`ThetaStep|_[__]?HasIndicesQ) :>
	(DecInt /@ m),
    Lattice`Private`SUM[i_, a_, b_, x_] :> Lattice`Private`SUM[
	i, DecInt[a], DecInt[b], x /. replaceIndicesDispatch],
    Lattice`Private`ThetaStep[a_, b_, x_] :> Lattice`Private`ThetaStep[
	DecInt[a], DecInt[b], x /. replaceIndicesDispatch]
}];

DecInt[index_Integer] := index - 1;

DecInt[index_] := index;

ToCExp[parametrization_] := SimplifyParametrization[parametrization] //.
    toCExpDispatch /.
    replaceIndicesDispatch // PrivatizeParameterReIm;

ToCExp[parametrization_, array_Symbol] := ToCExp[parametrization] /.
    d:drv[(Lattice`Private`Re|Lattice`Private`Im)[_],
	  (Lattice`Private`Re|Lattice`Private`Im)[_]] :>
	Symbol[DrvToCFormString[d]][array] /.
    p:(Lattice`Private`Re|Lattice`Private`Im)[_] /; ValueQ@ToEnumSymbol[p] :>
	array@ToEnumSymbol[p] /.
    ap:(Lattice`Private`Re|Lattice`Private`Im)[abbr_Symbol] :> ap[array];

SimplifyParametrization[parametrization_] := Simplify[
    parametrization,
    TimeConstraint -> FlexibleSUSY`FSSimplifyBetaFunctionsTimeConstraint
];

SimplifyParametrization[parametrization_] := Expand[parametrization];

ReCExp[cexp_?CRealTypeQ] := cexp;

ReCExp[cexp_] := (
    WriteString["stdout", "Real expression ", cexp,
		" assumes complex C type,\ntaking its real part\n"];
    Re[cexp]
);

CRealTypeQ[z_Plus|z_Times|
	   z_Power|
	   z_Sin|z_ArcSin|
	   z_Cos|z_ArcCos|
	   z_Tan|z_ArcTan] := And @@ (CRealTypeQ /@ List@@z);

CRealTypeQ[_AbsSqr] := True;

CRealTypeQ[Lattice`Private`SUM[_, _, _, z_] |
	   Lattice`Private`ThetaStep[_, _, z_] |
	   HoldPattern@SARAH`sum[_, _, _, z_]] := CRealTypeQ[z];

CRealTypeQ[_SARAH`Delta|_KroneckerDelta|_SARAH`ThetaStep] := True;

CRealTypeQ[_SARAH`Mass|_SARAH`Mass2|
	   _Lattice`Private`M|_Lattice`Private`M2] := True;

CRealTypeQ[_SARAH`A0|_SARAH`B0|_SARAH`B1|_SARAH`B00|_SARAH`B22|
	   _SARAH`F0|_SARAH`G0|_SARAH`H0] := True;

CRealTypeQ[_Re|_Im|_Lattice`Private`Re|_Lattice`Private`Im] := True;

CRealTypeQ[_Lattice`Private`x] := True;

CRealTypeQ[_?DeclaredRealQ] := True;

CRealTypeQ[z_?NumericQ] := Element[z, Reals];

CRealTypeQ[cexp_] := False;

Differentiate[exp_, x_, abbrRules_] :=
    D[exp, x, NonConstants -> abbrRules[[All,1]]] /.
    HoldPattern@D[f_, y_, ___] -> drv[f, y]

ParametrizeBetaFunctions[betaFunctions_List, sarahAbbrs_, parameterRules_] :=
Module[{
	convertedBetaFunctions = betaFunctions /.
	    sarahOperatorReplacementRules,
	convertedSarahAbbrs = Traces`ConvertSARAHTraces[sarahAbbrs] /.
	    sarahOperatorReplacementRules,
	nestedTraceRules, traceRules,
	symbolizeSarahAbbrs, sarahAbbrRules,
	nBetaFunctions
    },
    symbolizeSarahAbbrRules =
       (# -> ToValidCSymbol[#])& /@ convertedSarahAbbrs[[All,1]];
    nestedTraceRules = AbbreviateTraces[
	{convertedBetaFunctions, convertedSarahAbbrs}];
    traceRules = Flatten[nestedTraceRules];
    convertedBetaFunctions = convertedBetaFunctions /.
       traceRules /. symbolizeSarahAbbrRules;
    sarahAbbrRules = ParametrizeSarahAbbrs[
	convertedSarahAbbrs,
	Join[parameterRules, traceRules, symbolizeSarahAbbrRules]];
    nBetaFunctions = Length[convertedBetaFunctions];
    {MapIndexed[
	(WriteString["stdout",
		     "[", First[#2], "/", nBetaFunctions, "] expanding"];
	 ParametrizeBetaFunction[
	     #1, Join[parameterRules, traceRules, sarahAbbrRules[[All,1]]]])&,
	convertedBetaFunctions],
     Flatten[{ParametrizeTraces[nestedTraceRules, parameterRules],
	      sarahAbbrRules[[All,2]]}]}
];

ParametrizeBetaFunction[
    BetaFunction`BetaFunction[name_, _, betanLs_List], partRules_] :=
Module[{
	result
    },
    result =
       MapIndexed[
	   Flatten@ParametrizeBetanLRules[name, #1, First[#2], partRules]&,
	   betanLs];
    WriteString["stdout", "\n"];
    result
];

ParametrizeBetanLRules[name_, betanL_, n_, partRules_] := Module[{
	equations
    },
Done[
    equations = ParametrizeBetanL[name, betanL, n, partRules];
    EquationsToRules /@ equations,
	" BETA[",n ,", ", name, "]... "
]];

EquationsToRules[equations:HoldPattern@And[(BETA[_Integer, _] == _)..]] :=
    List@@equations /. Equal -> Rule;

EquationsToRules[(lhs:BETA[_Integer, _]) == rhs_] := {lhs -> rhs};

EquationsToRules[True] := {};

EquationsToRules[args__] :=
    (Print["Lattice`EquationsToRules[",args,"] failed."]; Abort[]);

ParametrizeBetanL[name_[i_,j_]?HasIndicesQ, betanL_, n_, partRules_] :=
    ParametrizeMatrixBetanL[name[i,j], betanL, n, partRules];

ParametrizeBetanL[name_[d__]?HasIndicesQ, betanL_, n_, partRules_] :=
    ParametrizeIndexedBetanL[name[d], betanL, n, partRules];

ParametrizeBetanL[name_, betanL_, n_, partRules_] := Module[{
	lBetanL, rBetanL,
	m
    },
    lBetanL = Betaize[n, Hold[name] /. partRules];
    rBetanL = Hold[betanL] /. partRules;
    ReduceBetaEquations[{#}]& /@
    SeparateParts[
	ReleaseHold[rBetanL - lBetanL //. matrixOpRules],
	partRules]
];

ParametrizeMatrixBetanL[name_[i_,j_]?HasIndicesQ, betanL_, n_, partRules_] :=
Module[{
	lBetanL, rBetanL,
	m
    },
    lBetanL = Betaize[n, Hold[name] /. partRules];
    rBetanL = Hold[betanL] /. KroneckerRule[name[i,j]] /. m_[i,j] :> m /.
       partRules;
    ReduceBetaEquations@Flatten[#]& /@
    SeparateParts[
	ReleaseHold[rBetanL - lBetanL //. matrixOpRules],
	partRules]
];

ParametrizeIndexedBetanL[name_[i__]?HasIndicesQ, betanL_, n_, partRules_] :=
Module[{
	dimensions = CouplingDimensions[name],
	lBetanL, rBetanL,
	loopArgs,
	m, l,
	unrollings
    },
    loopArgs = MapIndexed[{l@@#2, #1}&, dimensions];
    lBetanL = Betaize[n, Hold[name[i]] /. m_[i] :> m[[i]] /. partRules];
    rBetanL = Hold[betanL] /. KroneckerRule[name[i]] /. m_[i] :> m[[i]] /.
       partRules;
    (* TODO: unroll sum[]'s if SARAH`NoMatrixMultiplication -> True *)
    unrollings = Flatten[
	Table @@ Prepend[loopArgs, Thread[{i} -> loopArgs[[All, 1]]]],
	Length[dimensions]-1];
    ReduceBetaEquations@Flatten[#]& /@
    SeparateParts[
	ReleaseHold[rBetanL - lBetanL //. matrixOpRules /. #]& /@
	unrollings,
	partRules]
];

KroneckerRule[name_[i__]] := Module[{
	dimensions,
	d
    },
    With[{dimensions = CouplingDimensions[name]},
	Kronecker[d__] :> IdentityMatrix[
	    Flatten[Extract[dimensions, Position[{i}, #]]& /@ {d}]][d]]
];

ReduceBetaEquations[equations_] := Module[{
	variables =
	    Union@Cases[equations, BETA[_Integer, (Re|Im)[__]], {0, Infinity}]
    },
    Reduce[# == 0& /@ equations, variables]
];

SeparateParts[z_, partRules_] := Module[{
	complexes = Union[Variables[partRules[[All, 2]]] /.
			  Re|Im -> Identity]
    },
    ComplexExpand[Through[{Re, Im}[z]], complexes]
];

Betaize[level_, parametrization_] :=
    parametrization /. {
	Re[z_] :> BETA[level, Re[z]],
	Im[z_] :> BETA[level, Im[z]]
    };

ParametrizeSarahAbbrs[sarahAbbrs:{Traces`SARAHTrace[_,_]...}, partRules_] :=
    SarahAbbrToRule /@
    ReleaseHold[Hold[sarahAbbrs] /. partRules //. matrixOpRules]

SarahAbbrToRule[Traces`SARAHTrace[abbr_, rhs_]] :=
    {abbr -> Re[abbr], {Re[abbr] -> rhs}} /; RealQ[rhs];

SarahAbbrToRule[Traces`SARAHTrace[ab_, rhs_]] :=
    {ab -> Re[ab] + I Im[ab], {Re[ab] -> Re[rhs], Im[ab] -> Im[rhs]}};

ParametrizeTraces[nestedTraceRules_List, parameterRules_] :=
    Flatten[ParametrizeTrace[#, parameterRules]& /@ nestedTraceRules]

ParametrizeTrace[{traces_Alternatives -> abbr_, ___}, parameterRules_] :=
    ParametrizeTrace[{First[traces] -> abbr}, parameterRules]

ParametrizeTrace[{trace_SARAH`trace -> Re[ab_] + I Im[ab_], Repeated[_,{0,1}]},
		 parameterRules_] :=
    Thread[{Re[ab], Im[ab]} ->
	   SeparateParts[
	       ReleaseHold[Hold[trace] /. parameterRules //. matrixOpRules],
	       parameterRules]];

ParametrizeTrace[{trace_SARAH`trace -> Re[ab_], Repeated[_,{0,1}]},
		 parameterRules_] :=
Module[{
	re, im,
	z = ReleaseHold[Hold[trace] /. parameterRules //. matrixOpRules]
    },
    {re, im} = SeparateParts[z, parameterRules];
    Assert[RealQ[z]];
    {Re[ab] -> re}
];

ParametrizeTrace[args__] :=
    (Print["Lattice`ParametrizeTrace[",args,"] failed."]; Abort[]);

AbbreviateTraces[exp_] := Module[{
	traces = Union@Flatten@Cases[exp, _SARAH`trace, {0, Infinity}],
	grouped
    },
    grouped =
	Gather[Gather[traces, TraceSameQ],
	       TraceSameQ[Parametrization`cnj@First[#1], First[#2]]&];
    ConjugateTraceAbbrRule /@
    Map[ParametrizeTraceAbbrRule, Map[TraceAbbrRule, grouped, {2}], {2}]
];

ConjugateTraceAbbrRule[rule:{_}] := rule;

ConjugateTraceAbbrRule[{a_ -> b_, c_ -> d_}] :=
    {a -> b, c -> Parametrization`cnj[b]};

ParametrizeTraceAbbrRule[traces_Alternatives -> abbr_] :=
    traces -> Last@ParametrizeTraceAbbrRule[First[traces] -> abbr];

(* Q: can a trace also be purely imaginary? *)

ParametrizeTraceAbbrRule[trace_SARAH`trace?TraceRealQ -> abbr_] :=
    trace -> Re[abbr];

ParametrizeTraceAbbrRule[trace_SARAH`trace -> abbr_] :=
    trace -> Re[abbr] + I Im[abbr];

TraceAbbrRule[{trace_}] := trace -> CConversion`ToValidCSymbol[trace];

TraceAbbrRule[traces:{fst_,__}] :=
    Alternatives@@traces -> CConversion`ToValidCSymbol[fst];

TraceRealQ[trace_SARAH`trace] := TraceSameQ[trace, Parametrization`cnj[trace]];

TraceSameQ[SARAH`trace[a__], SARAH`trace[b__]] :=
    ListSameUptoRotationQ[{a}, {b}] ||
    ListSameUptoRotationQ[{a}, Parametrization`trp /@ Reverse[{b}]];

ListSameUptoRotationQ[a_List, b_List] :=
    Length[a] === Length[b] &&
    Or @@ Table[b === RotateLeft[a, i], {i, 0, Length[a]-1}];

matrixOpRules = Dispatch[{
    SARAH`MatMul :> Dot,
    Parametrization`trp :> Transpose,
    Parametrization`cnj :> Conjugate,
    Parametrization`adj :> ConjugateTranspose,
    SARAH`trace[m__] :> Tr[Dot[m]]
}];

unrollInnermostSumDispatch = Dispatch[{
    Lattice`Private`SUM[i_, a_, b_, x_] /;
	FreeQ[x, _Lattice`Private`SUM | _Lattice`Private`USUM] :>
	Lattice`Private`USUM[i, a, b, x]
}];

cExpToCFormStringDispatch = Dispatch[{
(*
    p:Power[_?NumericQ, _?NumericQ] :> N[p],
*)
    Power[E,z_]			    :> exp[z],
    Power[a_,2]			    :> Global`Sqr[a],
    Power[a_,-2]		    :> 1/Global`Sqr[a],
    Lattice`Private`ThetaStep[a_, b_, x_] :>
	Lattice`Private`LispAnd[HoldForm[a <= b], x]
}];

prefixCScopeDispatch = Dispatch[{
    f_?(ValueQ@CScope[#]&) :> InCScope[CScope[f], f]
}];

CExpToCFormString[expr_] :=
    StringReplace[ToString[expr //. unrollInnermostSumDispatch //.
			   cExpToCFormStringDispatch /.
			   prefixCScopeDispatch, CForm,
			   CharacterEncoding -> "ASCII"],
		  RegularExpression["\\\\\\[(.*?)\\]"] :> "$1"];

SetDependenceNode[form_, expr_, cscope_String:"CLASSNAME::Interactions::"] := (
    DependenceNode[form] = expr /. {
	HoldPattern[SARAH`Mass [f_]] :> Lattice`Private`M [f],
	HoldPattern[SARAH`Mass2[f_]] :> Lattice`Private`M2[f]
    };
    CScope[form] = cscope;
);

FindDependence[expr_] := Module[{
	explicit = Cases[
	    expr,
	    p:Re[_]|Im[_] /; With[{cexp=ToCExp[p]}, ValueQ@ToEnumSymbol@cexp],
	    {0, Infinity}],
	children
    },
    children = Complement[
	Cases[expr, _?(ValueQ@DependenceNode[#]&), {0, Infinity}],
	explicit /. Re|Im -> Identity];
    Union @
    If[children === {},
       explicit,
       Flatten[{explicit, FindDependence /@ Union[DependenceNode/@children]}]]
];

End[] (* `Private` *)

EndPackage[]
