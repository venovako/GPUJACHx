(* ::Package:: *)

BeginPackage["ListLexOrder`"]
ListLexOrder::usage="ListLexOrder[l1,l2] gives True if List l1 is lexicographically before List l2."
Begin["`Private`"]
ListLexOrder[l1_List,l2_List]:=Module[{i,j1=Length[l1],j2=Length[l2]},j=Min[j1,j2];r=Less[j1,j2];For[i=1,i<=j,++i,If[Equal[l1[[i]],l2[[i]]],Continue[]];r=Less[l1[[i]],l2[[i]]];Break[]];r]
End[]
EndPackage[]
