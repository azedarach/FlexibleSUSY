BeginPackage["Observables`", {"FlexibleSUSY`", "Utils`"}];

(* observables *)
Begin["FlexibleSUSYObservable`"];
FSObservables = { aMuonGM2Calc };
End[];

CalculateObservables::usage="";

Begin["`Private`"];

CalculateObservable[obs_ /; obs === FlexibleSUSYObservable`aMuonGM2Calc, structName_String] :=
    structName <> ".AMUGM2CALC = gm2calc::fs_calculate_amu();\n";

CalculateObservables[something_, structName_String] :=
    Module[{observables},
           observables = Cases[something, a_?(MemberQ[FlexibleSUSYObservable`FSObservables,#]&) :> a, {0, Infinity}];
           Utils`StringJoinWithSeparator[CalculateObservable[#,structName]& /@ observables, "\n"]
          ];

End[];

EndPackage[];