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

BeginPackage["LatticeUtils`", {
    "SARAH`",
    "TreeMasses`"}]

FixDiagonalization::usage;
SingleCase::usage;
Done::usage;
DoneLn::usage;
FormatShortTime::usage;
MajoranaQ::usage;
MajoranaMassMatrixQ::usage;

Begin["`Private`"]

FixDiagonalization[fsMassMatrices_List] := FixDiagonalization/@fsMassMatrices;

(*
   Diagonalization conventions of SARAH:

      SVD: m = u^T diag v,
	where u and v are the 1st and the 2nd mixing matrices from
	DEFINITION[_][MatterSector]

      hermitian: m = z^dagger diag z

   According to the SARAH documentation, the specification of the
   neutraliino mass matrix is indistinguishable from that of a
   hermitian matrix even though it must be diagonalized as

      symmetric: m = u^T diag u

   This leads to the following amendment:
 *)
FixDiagonalization[FSMassMatrix[m_, f_, z_Symbol]?MajoranaMassMatrixQ] :=
    FSMassMatrix[m, f, {z, z}];

FixDiagonalization[
    FSMassMatrix[{m_}?VectorQ, f_?IsFermion /; GetDimension[f] > 1, Null]] :=
    FSMassMatrix[DiagonalMatrix[Table[m, {GetDimension[f]}]], f,
		 {IdentityMatrix, IdentityMatrix}];

FixDiagonalization[
    FSMassMatrix[{m_}?VectorQ, f_ /; GetDimension[f] > 1, Null]] :=
    FSMassMatrix[DiagonalMatrix[Table[m, {GetDimension[f]}]], f,
		 IdentityMatrix];

FixDiagonalization[m_FSMassMatrix] := m;

MajoranaMassMatrixQ[FSMassMatrix[_?MatrixQ, _, _]?MajoranaMassQ] := True;

MajoranaMassMatrixQ[_FSMassMatrix] := False;

MajoranaMassQ[FSMassMatrix[_, _?MajoranaQ, _]] := True;

MajoranaMassQ[_FSMassMatrix] := False;

MajoranaQ[field_] := MemberQ[SARAH`MajoranaPart, field];

SingleCase[args__] := Module[{
	cases = Cases[args]
    },
    Assert[Length[cases] === 1];
    First[cases]
];

SetAttributes[Done, HoldFirst];

Done[exp_, msg__] := Module[{
	result
    },
    WriteString["stdout", msg];
    result = Timing[exp];
    WriteString["stdout", FormatShortTime @ First[result]];
    Last[result]
];

SetAttributes[DoneLn, HoldFirst];

DoneLn[exp_, msg__] := Module[{
	result = Done[exp, msg]
    },
    WriteString["stdout", "\n"];
    result
];

FormatShortTime[seconds_] := Module[{
	time = Round[seconds 1*^3]
    },
    If[time === 0,
       ToString[Round[seconds 1*^6]] <> " us",
       ToString[time] <> " ms"]
];

End[] (* `Private` *)

EndPackage[]
