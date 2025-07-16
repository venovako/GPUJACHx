(* ::Package:: *)

BeginPackage["CollideGraph`",{"PivotsCollideQ`"}]
CollideGraph::usage="CollideGraph[s] creates the collide graph of a Jacobi strategy s."
Begin["`Private`"]
CollideGraph[s_List]:=Module[{i,j,k=Length[s],l={}},For[i=1,i<k,++i,For[j=i+1,j<=k,++j,If[PivotsCollideQ[s[[i]],s[[j]]],AppendTo[l,j\[DirectedEdge]i]]]];Graph[Table[Tooltip[Labeled[i,i,Center],s[[i]]],{i,1,k}],l,GraphLayout->"CircularEmbedding",PerformanceGoal->"Quality",VertexShapeFunction->"Square",VertexSize->Automatic,VertexStyle->White,EdgeStyle->Black]]
End[]
EndPackage[]
