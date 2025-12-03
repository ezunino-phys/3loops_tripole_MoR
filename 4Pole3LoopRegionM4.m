(* ::Package:: *)

AppendTo[$Path,"/home/s2894494"];
$AROptions={TimeConstraint->Infinity,"UseKira"->True,FileBaseName->"4Pole3LoopRegionM4"};
Import["/home/s2894494/ampred/AmpRed/AmpRed.m"];
LaunchKernels[32];
SetOptions[DoKira,{"Executable" ->{"/home/s1996839/kira/bin/kira","--parallel=32","--integral_ordering=2"},"UserDefinedSystem"->True,Method->2,AlphaBasis->True,"KiraJobFileOptions"->{"run_symmetries: true","run_initiate: true","run_triangular: true","run_back_substitution: true","run_firefly: true"}}];
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
int=Expand[Simplify[ToFeynmanInt[Num1*Dom1+Num2*Dom2,{q1,q2,q3}]]]/.I\[Epsilon]->0;
intpp=Sum[Print[i];AlphaParametrize[int[[i]]],{i,1,Length[int]}];
allgstep=Union[Cases2[intpp,AlphaInt[__]],{}];
repLLstep=Table[Print[i];allgstep[[i]]->AlphaSeries[allgstep[[i]],{lam,0,3}],{i,1,Length[allgstep]}];
intpptemp=Map[Simplify,Collect[intpp/.repLLstep,AlphaInt[__]]];
intppCoLL=Assuming[{lam>0},Expand[PowerExpand[Map[Simplify,Collect[intpptemp,AlphaInt[__]]],{lam}]]];
Regionm4=Sum[Simplify[intppCoLL[[i]]*lam^(4eps)]/.lam^x_->0,{i,1,Length[intppCoLL]}]
RegionExpm4=Simplify[Normal[Series[Regionm4,{lam,0,0}]]]
Reducedm4=Map[Simplify,Collect[DoKira[RegionExpm4],AlphaInt[__]]]
ReducedEvam4=AlphaIntEvaluate[Reducedm4]
DiffEvam4aIJ=AlphaDES[Cases2[Reducedm4,AlphaInt[__]],aIJ]
DiffEvam4yIJK=AlphaDES[Cases2[Reducedm4,AlphaInt[__]],SQRTyIJK]
ExplicitExpm4=AlphaIntExplicit[Cases2[Reducedm4,AlphaInt[__]]]
SetDirectory["/home/s2894494"];
DumpSave["4Pole3LoopRegionM4.mx",{Regionm4,RegionExpm4,Reducedm4,ReducedEvam4,DiffEvam4aIJ,DiffEvam4yIJK,ExplicitExpm4}];
