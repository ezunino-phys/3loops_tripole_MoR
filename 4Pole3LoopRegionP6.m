(* ::Package:: *)

AppendTo[$Path,"/home/s2894494"];
$AROptions={TimeConstraint->Infinity,"UseKira"->True,FileBaseName->"4Pole3LoopRegionP6"};
Import["/home/s2894494/ampred/AmpRed/AmpRed.m"];
LaunchKernels[14];
SetOptions[DoKira,{"UserDefinedSystem"->True,Method->2,AlphaBasis->True,"KiraJobFileOptions"->{"run_symmetries: true","run_initiate: true","run_triangular: true","run_back_substitution: true","run_firefly: true"}}];
SetOptions[AlphaDES,{Method->2,AlphaReduce->DoKira}];
SetOptions[QhullNorms, "Executable"->"/home/s2894494/ampred/AmpRed/qhull-2020.2/bin/qhull"];
SP[vI,vJ]=-1/2(aIJ+1/aIJ);
SP[vI,vI]=1;
SP[vJ,vJ]=1;
SP[vK,vK]=1;
SP[vK,vI]=-1/2(aIK+1/aIK)/.aIK->lam/SQRTyIJK;
SP[vK,vJ]=-1/2(aJK+1/aJK)/.aJK->lam SQRTyIJK;
Num1= (-4 SP[q3,vK]^2 SP[vI,vJ]+4 SP[q1,vJ] SP[q3,vK] SP[vI,vK]+SP[q2,vJ] SP[q3,vK] SP[vI,vK]+3 SP[q3,vJ] SP[q3,vK] SP[vI,vK]-SP[q1,vI] SP[q3,vK] SP[vJ,vK]-4 SP[q2,vI] SP[q3,vK] SP[vJ,vK]+3 SP[q3,vI] SP[q3,vK] SP[vJ,vK]+4 SP[q1,q2] SP[vI,vK] SP[vJ,vK]-2 SP[q1,q3] SP[vI,vK] SP[vJ,vK]+2 SP[q2,q3] SP[vI,vK] SP[vJ,vK]-SP[q3,q3] SP[vI,vK] SP[vJ,vK]+SP[q2,vK] (2 SP[q3,vK] SP[vI,vJ]-2 SP[q1,vJ] SP[vI,vK]-SP[q3,vJ] SP[vI,vK]-SP[q1,vI] SP[vJ,vK]+SP[q3,vI] SP[vJ,vK])+SP[q1,vK] (SP[q2,vK] SP[vI,vJ]-2 SP[q3,vK] SP[vI,vJ]-SP[q2,vJ] SP[vI,vK]-SP[q3,vJ] SP[vI,vK]-2 SP[q2,vI] SP[vJ,vK]+SP[q3,vI] SP[vJ,vK])+SP[q1,vI] SP[q2,vJ] SP[vK,vK]-SP[q2,vJ] SP[q3,vI] SP[vK,vK]+SP[q1,vI] SP[q3,vJ] SP[vK,vK]-SP[q3,vI] SP[q3,vJ] SP[vK,vK]);
Num2=(SP[vI, vK]*SP[vJ, vK] - SP[vI, vJ]*SP[vK, vK]);
Dom1=1/(SP[q1,vI]+SP[q3,vI]-1+I\[Epsilon])*1/(-SP[q3,vJ]+SP[q2,vJ]-1+I\[Epsilon])*1/(-SP[q2,vK]+I\[Epsilon])*1/(-SP[q1,vK]+I\[Epsilon])*1/(SP[q1,q1])*1/(SP[q2,q2])*1/(SP[q3,q3])*1/(SP[q1,q1]+2SP[q1,q3]+SP[q3,q3])*1/(SP[q3,q3]-2SP[q3,q2]+SP[q2,q2]);
Dom2=1/(SP[q1,vI]+SP[q3,vI]-1+I\[Epsilon])*1/(-SP[q3,vJ]+SP[q2,vJ]-1+I\[Epsilon])*1/(-SP[q2,vK]+I\[Epsilon])*1/(-SP[q1,vK]+I\[Epsilon])*1/(SP[q1,q1])*1/(SP[q2,q2])*1/(SP[q1,q1]+2SP[q1,q3]+SP[q3,q3])*1/(SP[q3,q3]-2SP[q3,q2]+SP[q2,q2]);
(*int1=Expand[Simplify[ToFeynmanInt[Num1*Dom1,{q1,q2,q3}]]]/.I\[Epsilon]->0;
int2=Expand[Simplify[ToFeynmanInt[Num2*Dom2,{q1,q2,q3}]]]/.I\[Epsilon]->0;
int=int1+int2;*)
int=Expand[Simplify[ToFeynmanInt[Num1*Dom1+Num2*Dom2,{q1,q2,q3}]]]/.I\[Epsilon]->0;
(*intSeries=Series[int,{lam,0,0}];*)
intpp=Sum[Print[i];AlphaParametrize[int[[i]]],{i,1,Length[int]}];
(*CheckDoKira=DoKira[intpp];*)
(*intppSeries=Series[intpp,{lam,0,0}];*)
allgstep=Union[Cases2[intpp,AlphaInt[__]],{}];
repLLstep=Table[Print[i];allgstep[[i]]->AlphaSeries[allgstep[[i]],{lam,0,3}],{i,1,Length[allgstep]}];
intpptemp=Map[Simplify,Collect[intpp/.repLLstep,AlphaInt[__]]];
intppCoLL=Assuming[{lam>0},Expand[PowerExpand[Map[Simplify,Collect[intpptemp,AlphaInt[__]]],{lam}]]];
Regionp6=Sum[Simplify[intppCoLL[[i]]*lam^(-6eps)]/.lam^x_->0,{i,1,Length[intppCoLL]}]
RegionExpp6=Simplify[Normal[Series[Regionp6,{lam,0,0}]]]
Reducedp6=Map[Simplify,Collect[DoKira[RegionExpp6],AlphaInt[__]]]
ReducedEvap6=AlphaIntEvaluate[Reducedp6]
DiffEvap6aIJ=AlphaDES[Cases2[Reducedp6,AlphaInt[__]],aIJ]
DiffEvap6yIJK=AlphaDES[Cases2[Reducedp6,AlphaInt[__]],SQRTyIJK]
ExplicitExpp6=AlphaIntExplicit[Cases2[Reducedp6,AlphaInt[__]]]
SetDirectory["/home/s2894494"];
DumpSave["4Pole3LoopRegionP6.mx",{Regionp6,RegionExpp6,Reducedp6,ReducedEvap6,DiffEvap6aIJ,DiffEvap6yIJK,ExplicitExpp6}];
