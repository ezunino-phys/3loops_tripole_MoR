(* ::Package:: *)

AppendTo[$Path,"/home/s1996839"];
$AROptions={TimeConstraint->Infinity,"UseKira"->True,FileBaseName->"ThreeLoopIRRegionM010"};
Import["/home/s1996839/ampred/AmpRed/AmpRed.m"];
LaunchKernels[32];
SetOptions[DoKira, {"UserDefinedSystem"->True,Method->2,AlphaBasis->True,"KiraJobFileOptions"->{"run_symmetries: true","run_initiate: true","run_triangular: true","run_back_substitution: true","run_firefly: true"}}];
SetOptions[AlphaDES,{Method->2,AlphaReduce->DoKira}];
SetOptions[QhullNorms, "Executable"->"/home/s1996839/ampred/AmpRed/qhull/bin/qhull"];
SP[b1,b2]=2 z12;
SP[b1,b3]=2 z13;
SP[b2,b3]=2 z23;
SP[b1,b4]=(*g14*)-1;
SP[b2,b4]=(*g24*)-1;
SP[b3,b4]=(*g34*)-1;
SP[b1,b1]=2^2lam1 (*g14^2*);
SP[b2,b2]=2^2lam2 (*g24^2*);
SP[b3,b3]=2^2lam3 (*g34^2*);
SP[b4,b4]=1;
term[1]=Simplify[sp[b1,b2]sp[b3,b4]*(2sp[k1,k3]-2sp[k2,k3]+sp[k1,k1]-sp[k2,k2])/.sp->SP];
term[2]=sp[b1,b2](-sp[b3,k3]-2sp[b3,k1]-2sp[b3,k2])*(sp[b4,k1]-sp[b4,k2])/.sp->SP;
term[3]=sp[b1,b2](sp[b4,k1]+sp[b4,k2]-sp[b4,k3])*(sp[b3,k1]-sp[b3,k2])/.sp->SP; 
term[4]=(sp[b1,k1]+2sp[b1,k2])(2sp[b2,k3]+sp[b2,k1]+sp[b2,k2])sp[b3,b4]/.sp->SP;
term[5]=(sp[b1,k1]+2sp[b1,k2])sp[b2,b4](-sp[b3,k3]-2sp[b3,k1]-2sp[b3,k2])/.sp->SP;
term[6]=(sp[b1,k1]+2sp[b1,k2])sp[b2,b3](sp[b4,k1]+sp[b4,k2]-sp[b4,k3])/.sp->SP;
term[7]=(2sp[b1,k3]+sp[b1,k1]+sp[b1,k2])(-sp[b2,k2]-2sp[b2,k1])sp[b3,b4]/.sp->SP;
term[8]=sp[b1,b4](-sp[b2,k2]-2sp[b2,k1])(-sp[b3,k3]-2sp[b3,k1]-2sp[b3,k2])/.sp->SP;
term[9]=sp[b1,b3](-sp[b2,k2]-2sp[b2,k1])(sp[b4,k1]+sp[b4,k2]-sp[b4,k3])/.sp->SP;

NumThreeLoops=Simplify[Sum[term[i],{i,1,9}]];
int3a=Expand[Simplify[ToFeynmanInt[(NumThreeLoops)*1/(SP[k1,k1])*1/(SP[k2,k2])*1/(SP[k3,k3])*1/(SP[k1+k2,k1+k2])*1/(SP[k3+k1+k2,k3+k1+k2])*1/(-SP[b1,k1]+ie)*1/(-SP[k2,b2]+ie)*1/(-SP[k3,b3]+ie)*1/(SP[k3+k1+k2,b4]-1+ie),{k1,k2,k3}]]]/.ie->0;
int3b=Expand[Simplify[ToFeynmanInt[((SP[b1,b3]SP[b2,b4]-SP[b1,b4]SP[b2,b3]))*1/(SP[k1,k1])*1/(SP[k2,k2])*1/(SP[k3,k3])(**1/(SP[k1+k2,k1+k2])*)*1/(SP[k3+k1+k2,k3+k1+k2])*1/(-SP[b1,k1]+ie)*1/(-SP[k2,b2]+ie)*1/(-SP[k3,b3]+ie)*1/(SP[k3+k1+k2,b4]-1+ie),{k1,k2,k3}]]]/.ie->0;
int3=int3a+int3b;
intpp3=Sum[Print[i];AlphaParametrize[int3[[i]]],{i,1,Length[int3]}];
(*Expansion in lam1*)
allgstep1=Union[Cases2[intpp3,AlphaInt[__]],{}];
repLLstep1=Table[Print[i];allgstep1[[i]]->AlphaSeries[allgstep1[[i]],{lam1,0,0}],{i,1,Length[allgstep1]}];
intpp3temp1=Map[Simplify,Collect[intpp3/.repLLstep1,AlphaInt[__]]];
(*Expansion in lam2*)
allgstep2=Union[Cases2[intpp3temp1,AlphaInt[__]],{}];
repLLstep2=Table[Print[i];allgstep2[[i]]->AlphaSeries[allgstep2[[i]],{lam2,0,0}],{i,1,Length[allgstep2]}];
intpp3temp2=Map[Simplify,Collect[intpp3temp1/.repLLstep2,AlphaInt[__]]];
(*Expansion in lam3*)
allgstep3=Union[Cases2[intpp3temp2,AlphaInt[__]],{}];
repLLstep3=Table[Print[i];allgstep3[[i]]->AlphaSeries[allgstep3[[i]],{lam3,0,0}],{i,1,Length[allgstep3]}];
intpp3temp3=Map[Simplify,Collect[intpp3temp2/.repLLstep3,AlphaInt[__]]];
(*Integral after expansion*)
intpp3CoLL=Assuming[{lam1>0,lam2>0,lam3>0},Expand[PowerExpand[Map[Simplify,Collect[intpp3temp3,AlphaInt[__]]],{lam1,lam2,lam3}]]];
(*-2 eps*)
Regionm[0,1,0]=Sum[Coefficient[intpp3CoLL[[i]],lam2^(- eps)]/.lam1^x_->0/.lam2^x_->0/.lam3^x_->0,{i,1,Length[intpp3CoLL]}]
Regionm010Exp=Simplify[Normal[Series[Normal[Series[Normal[Series[Regionm[0,1,0],{lam1,0,0}]],{lam2,0,0}]],{lam3,0,0}]]]
Reducedm010=Map[Simplify,Collect[DoKira[Regionm010Exp],AlphaInt[__]]]
Reducedm010Eva=AlphaIntEvaluate[Reducedm010]
Diffm010Evaz12=AlphaDES[Cases2[Reducedm010,AlphaInt[__]],z12]
Diffm010Evaz13=AlphaDES[Cases2[Reducedm010,AlphaInt[__]],z13]
Diffm010Evaz23=AlphaDES[Cases2[Reducedm010,AlphaInt[__]],z23]
SetDirectory["/home/s1996839"];
DumpSave["ThreeLoopIRm010.mx",{Regionm010Exp,Reducedm010,Reducedm010Eva,Diffm010Evaz12,Diffm010Evaz13,Diffm010Evaz23}]
