f[1]=1-x1-x2
f[2]=x1
f[3]=x2

f[1,1]=D[f[1],{x1}]
f[1,2]=D[f[1],{x2}]
f[2,1]=D[f[2],{x1}]
f[2,2]=D[f[2],{x2}]
f[3,1]=D[f[3],{x1}]
f[3,2]=D[f[3],{x2}]

Print["D=zeros(3,3,2,2);"]
Do[Print["D(",i,",",j,",",a,",",b,")=",Integrate[f[i,a]*f[j,b]*Boole[x1+x2<=1],{x1,0,1},{x2,0,1}],";"],{i,1,3},{j,1,3},{a,1,2},{b,1,2}]

Print[""]

Print["E=zeros(3,3);"]
Do[Print["E(",i,",",j,")=",Integrate[f[i]*f[j]*Boole[x1+x2<=1],{x1,0,1},{x2,0,1}],";"],{i,1,3},{j,1,3}]

Print[""]

Print["C=zeros(3,3,2);"]
Do[Print["C(",i,",",j,",",a,")=",Integrate[f[i]*f[j,a]*Boole[x1+x2<=1],{x1,0,1},{x2,0,1}],";"],{i,1,3},{j,1,3},{a,1,2}]