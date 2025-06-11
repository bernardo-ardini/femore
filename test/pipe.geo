h=0.1;

L=2;
l=1;
S=2;
s=1;

Point(1)={0,0,0,h};
Point(2)={L,0,0,h};
Point(3)={L,s/2,0,h};
Point(4)={L+l,s/2,0,h};
Point(5)={L+l,0,0,h};
Point(6)={2*L+l,0,0,h};
Point(7)={2*L+l,S,0,h};
Point(8)={L+l,S,0,h};
Point(9)={L+l,S-s/2,0,h};
Point(10)={L,S-s/2,0,h};
Point(11)={L,S,0,h};
Point(12)={0,S,0,h};

For k In {1:11}
   Line(k)={k,k+1};
EndFor

Line(12)={12,1};

Curve Loop(1)={1:12};

Plane Surface(1)={1};

Physical Curve(1)={1:12};
Physical Surface(800)={1};

Mesh 2;
Save "mesh.m";