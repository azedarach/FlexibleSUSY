Needs["TestSuite`", "TestSuite.m"];
Needs["Parameters`", "Parameters.m"];

Print["testing IsRealExpression[] ..."];

realPars = {a, b, c};
compPars = {x, y, z};
allPars  = Join[realPars, compPars];

Parameters`SetModelParameters[allPars];
Parameters`AddRealParameter[realPars];
SARAH`RealParameters = {};

TestEquality[Parameters`IsRealParameter[a], True];
TestEquality[Parameters`IsRealParameter[b], True];
TestEquality[Parameters`IsRealParameter[c], True];
TestEquality[Parameters`IsRealParameter[x], False];
TestEquality[Parameters`IsRealParameter[y], False];
TestEquality[Parameters`IsRealParameter[z], False];

TestEquality[Parameters`IsRealExpression[a], True];
TestEquality[Parameters`IsRealExpression[x], False];

TestEquality[Parameters`IsRealExpression[a^2], True];
TestEquality[Parameters`IsRealExpression[a b], True];
TestEquality[Parameters`IsRealExpression[a b c], True];

TestEquality[Parameters`IsRealExpression[Sin[a] b], True];

TestEquality[Parameters`IsRealExpression[a b c x], False];
TestEquality[Parameters`IsRealExpression[x^2], False];
TestEquality[Parameters`IsRealExpression[Conjugate[x] x], True];
TestEquality[Parameters`IsRealExpression[conj[x] x], True];
TestEquality[Parameters`IsRealExpression[Conj[x] x], True];

TestEquality[Parameters`IsRealExpression[trace[a]], True];
TestEquality[Parameters`IsRealExpression[trace[x]], False];

TestEquality[Parameters`IsRealExpression[trace[a b]], True];

TestEquality[Parameters`IsRealExpression[trace[Adj[a],a]], True];
TestEquality[Parameters`IsRealExpression[trace[Adj[x],x]], True];

TestEquality[Parameters`IsRealExpression[trace[Adj[x],x,Adj[y],y]], True];
TestEquality[Parameters`IsRealExpression[trace[Adj[x],x,Adj[y],y,Adj[z],z]], False];

TestEquality[Parameters`IsRealExpression[trace[conj[a],Tp[b]]], True];
TestEquality[Parameters`IsRealExpression[trace[a,Tp[b]]], True];
TestEquality[Parameters`IsRealExpression[trace[a,b]], True];

TestEquality[Parameters`IsRealExpression[trace[conj[x],x]], False];
TestEquality[Parameters`IsRealExpression[trace[conj[x],Tp[x]]], True];
TestEquality[Parameters`IsRealExpression[trace[x,Tp[y]]], False];
TestEquality[Parameters`IsRealExpression[trace[x,y]], False];

TestEquality[Parameters`IsRealExpression[trace[conj[x],Tp[x],Adj[y],y]], True];

(* make a hermitian*)
SARAH`getDimParameters[a] = {1};
TestEquality[Parameters`Private`IsHermitian[a], True];

TestEquality[Parameters`IsRealExpression[trace[Adj[x],x,a]], True];

(* make x hermitian*)
SARAH`getDimParameters[x] = {2,2};
SARAH`ListSoftBreakingScalarMasses = {x};

TestEquality[Parameters`Private`IsHermitian[x], True];
TestEquality[Parameters`IsRealExpression[trace[x]], True];

(* test sum and products *)
TestEquality[Parameters`IsRealExpression[c + b trace[x]], True];
TestEquality[Parameters`IsRealExpression[c - b trace[x] - 2 trace[a]], True];

Print["testing FindAllParameters[] ..."];

modelParameters = { Mu, SARAH`B[Mu], WOp, SARAH`Q[WOp] };
Parameters`SetModelParameters[modelParameters];

expr = 2 * Mu SARAH`B[Mu] + WOp SARAH`Q[WOp] + a;

TestEquality[Sort[Parameters`FindAllParameters[expr]],
             Sort[modelParameters]
            ];

PrintTestSummary[];
