(* ::Package:: *)

BeginPackage["RowCyc`"]
RowCyc::usage="RowCyc[n] gives the Jacobi row-cyclic strategy of order n."
Begin["`Private`"]
RowCyc[n_Integer]:=Flatten[Module[{i,j},Table[{i,j},{i,1,n-1},{j,i+1,n}]],1]
End[]
EndPackage[]
