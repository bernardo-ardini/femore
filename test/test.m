clear;
close all;

addpath("~/femore/scr/");

% geometry

system("gmsh pipes.geo");

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
sp.delta=0.1;

g=@(x) [0;-1]*((1-x(1,:)).*x(1,:).*(x(2,:)>=6-eps).*(x(1,:)<=1+eps))+[0;1]*((3-x(1,:)).*(x(1,:)-2).*(x(2,:)>=6-eps).*(x(1,:)>=2-eps));
f=@(x) -[0;1]*ones(1,size(x,2));
r=@(x) zeros(1,size(x,2));
[u,p]=sp.solve(g,f,r);

% plot

u1=u.dof(V.vertexComponent2index([1:geo.numvertices;ones(1,geo.numvertices)]'));
u2=u.dof(V.vertexComponent2index([1:geo.numvertices;2*ones(1,geo.numvertices)]'));

style="None";

subplot(2,2,1);
title("p");
patch("Faces",geo.triangles,"Vertices",geo.vertices,'FaceVertexCData',p.dof,'FaceColor','interp','LineStyle',style);
pbaspect([1,1,1]);
daspect([1,1,1]);
colorbar;

subplot(2,2,2);
quiver(geo.vertices(:,1),geo.vertices(:,2),u1,u2);
title("u");
pbaspect([1,1,1]);
daspect([1,1,1]);

subplot(2,2,3);
title("u_1");
patch("Faces",geo.triangles,"Vertices",geo.vertices,'FaceVertexCData',u1,'FaceColor','interp','LineStyle',style);
pbaspect([1,1,1]);
daspect([1,1,1]);
colorbar;

subplot(2,2,4);
title("u_2");
patch("Faces",geo.triangles,"Vertices",geo.vertices,'FaceVertexCData',u2,'FaceColor','interp','LineStyle',style);
pbaspect([1,1,1]);
daspect([1,1,1]);
colorbar;