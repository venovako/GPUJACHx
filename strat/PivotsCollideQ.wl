(* ::Package:: *)

BeginPackage["PivotsCollideQ`"]
PivotsCollideQ::usage="PivotsCollideQ[{p,q},{r,s}] returns True iff {p,q}\[Intersection]{r,s}\[NotEqual]\[EmptySet]."
Begin["`Private`"]
PivotsCollideQ[{p_Integer,q_Integer},{r_Integer,s_Integer}]:=(p==r)||(p==s)||(q==r)||(q==s)
End[]
EndPackage[]
