(* ::Package:: *)

BeginPackage["AllIS`",{"ListLexOrder`"}]
AllIS::usage="AllIS[g,n] returns all independent sets of Graph g of size n in lexicographic order."
Begin["`Private`"]
AllIS[g_Graph,n_Integer]:=Sort[Map[Sort,FindIndependentVertexSet[g,{n},All]],ListLexOrder]
End[]
EndPackage[]
