clear;
close all;

addpath("../scr");

% geometry

%system("gmsh geo.geo");

geo=Geometry();
mesh;
geo.vertices=msh.POS(:,1:2);
geo.triangles=msh.TRIANGLES(:,1:3);
geo.lines=msh.LINES;
geo.initialize();

% function space

V=FunctionSpace(geo,"P1");
V.linesConstraints=1;
V.computeBasisFunctions();


sp=StokesProblem(V);
sp.mu=1;
sp.c=1;

[A,B]=sp.assemble();