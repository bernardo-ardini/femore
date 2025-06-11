h=0.04;

Point(1)={0,0,0,h};
Point(2)={1,0,0,h};
Point(3)={1,1,0,h};
Point(4)={0,1,0,h};

For k In {1:3}
   Line(k)={k,k+1};
EndFor

Line(4)={4,1};

Curve Loop(1)={1:4};

Plane Surface(1)={1};

Physical Curve(1)={1:4};
Physical Surface(800)={1};

Mesh 2;
Save "mesh.m";