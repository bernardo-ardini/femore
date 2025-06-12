h=0.1;

Point(1)={0,0,0,h};
Point(2)={3,0,0,h};
Point(3)={3,6,0,h};
Point(4)={2,6,0,h};
Point(5)={2,1,0,h};
Point(6)={1,1,0,h};
Point(7)={1,6,0,h};
Point(8)={0,6,0,h};

For k In {1:7}
   Line(k)={k,k+1};
EndFor

Line(8)={8,1};

Curve Loop(1)={1:8};

Plane Surface(1)={1};

Physical Curve(1)={1:8};
Physical Surface(800)={1};

Mesh 2;
Save "mesh.m";