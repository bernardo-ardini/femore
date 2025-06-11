clear;
close all;

addpath("../scr");

% geometry

system("gmsh pipe.geo");

geo=Geometry();
mesh;
geo.vertices=msh.POS(:,1:2);
geo.triangles=msh.TRIANGLES(:,1:3);
geo.lines=msh.LINES;
geo.initialize();

% function space

V=FunctionSpace(geo,"P12b");
V.setLinesConstraint(1);

Q=FunctionSpace(geo,"P1");
Q.constrainedVertices=1;

% Stokes problem

sp=StokesProblem(V,Q);
sp.mu=0.1;
sp.c=0;

[A,B]=sp.assemble();

[u,p]=sp.solve(@(x) (2-x(2))*x(2)*((abs(x(1))<=1e-3)+(abs(x(1)-5)<=1e-3))*[1;0]);

% plot

speed=sqrt(u.dof(1:geo.numvertices).^2+u.dof((1:geo.numvertices)+geo.numvertices+geo.numtriangles).^2);

style="None";

subplot(2,2,1);
title("p");
patch("Faces",geo.triangles,"Vertices",geo.vertices,'FaceVertexCData',p.dof,'FaceColor','interp','LineStyle',style);
pbaspect([1,1,1]);
daspect([1,1,1]);
colorbar;

subplot(2,2,2);
quiver(geo.vertices(:,1),geo.vertices(:,2),u.dof(1:geo.numvertices),u.dof((1:geo.numvertices)+geo.numvertices+geo.numtriangles));
title("u");
pbaspect([1,1,1]);
daspect([1,1,1]);

subplot(2,2,3);
title("u_1");
patch("Faces",geo.triangles,"Vertices",geo.vertices,'FaceVertexCData',u.dof(1:geo.numvertices),'FaceColor','interp','LineStyle',style);
pbaspect([1,1,1]);
daspect([1,1,1]);
colorbar;

subplot(2,2,4);
title("u_2");
patch("Faces",geo.triangles,"Vertices",geo.vertices,'FaceVertexCData',u.dof((1:geo.numvertices)+geo.numvertices+geo.numtriangles),'FaceColor','interp','LineStyle',style);
pbaspect([1,1,1]);
daspect([1,1,1]);
colorbar;