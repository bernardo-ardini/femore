a=6;
b=6*1.618;
s=1;
h=1/6;

Point(1)={0,0,0,h};
Point(2)={a,0,0,h};
Point(3)={a,b,0,h};
Point(4)={0,b,0,h};

Point(5)={s,s,0,h};
Point(6)={a-s,s,0,h};
Point(7)={a-s,b-s,0,h};
Point(8)={s,b-s,0,h};


For k In {1:3}
   Line(k)={k,k+1};
EndFor
Line(4)={4,1};

For k In {5:7}
   Line(k)={k,k+1};
EndFor
Line(8)={8,5};

Curve Loop(1)={1:4};
Curve Loop(2)={5:8};
Plane Surface(1)={1,2};

Physical Curve(1)={1};
Physical Curve(2)={3};
Physical Surface(800)={1};

Mesh 2;
Save "mesh.m";