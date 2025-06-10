clear;
close all;

addpath("../../scr/");

% geometry

system("gmsh geo.geo");

geo=Geometry();
mesh;
geo.vertices=msh.POS(:,1:2);
geo.triangles=msh.TRIANGLES(:,1:3);
geo.lines=msh.LINES;
geo.initialize();

system("rm mesh.m");

% function space

V=FunctionSpace(geo,"P12");
V.addConstrainFromLines(1);
V.computeBasisFunctions();

% reference configuration

u0=Function(V);
u0.computeJacobians();

% koiter analysis

koi=KoiterAnalysis(V);

% constitutive relations

eps=zeros(2);
eps(1,2)=1;
eps(2,1)=-1;
cof=@(F) [F(2,2),-F(2,1);-F(1,2),F(1,1)];

constitutiveModel="KSV";

if constitutiveModel=="NH"
    alpha=10;
    beta=1e4;
    koi.D2W=@(F) alpha*permute(tensorprod(eye(2),eye(2)),[1,3,2,4])+beta*(tensorprod(cof(F),cof(F))+(det(F)-1)*permute(tensorprod(eps,eps),[1,3,2,4]));
    koi.D3W=@(F) beta*(tensorprod(tensorprod(eps,eps),cof(F))+tensorprod(tensorprod(eps,cof(F)),eps)+tensorprod(tensorprod(cof(F),eps),eps));
elseif constitutiveModel=="KSV"
    lmb=120;
    mu=80;
    koi.S=@(F) lmb*(trace(F'*F)/2-1)*F+mu*(F*F'*F-F);
    koi.D2W=@(F) lmb*(tensorprod(F,F)+(trace(F'*F)/2-1)*permute(tensorprod(eye(2),eye(2)),[1,3,2,4]))+mu*(permute(tensorprod(eye(2),F'*F),[1,3,2,4])+permute(tensorprod(F,F),[1,4,3,2])+permute(tensorprod(eye(2),F*F'),[3,1,4,2])-permute(tensorprod(eye(2),eye(2)),[1,3,2,4]));
    koi.D3W=@(F) lmb*(permute(tensorprod(eye(2),tensorprod(eye(2),F)),[1,3,5,6,2,4])+permute(tensorprod(eye(2),tensorprod(eye(2),F)),[1,3,2,4,5,6])+permute(tensorprod(eye(2),tensorprod(eye(2),F)),[5,6,1,3,2,4]))+mu*(permute(tensorprod(eye(2),tensorprod(eye(2),F)),[5,4,1,3,2,6])+permute(tensorprod(eye(2),tensorprod(eye(2),F)),[5,4,1,6,2,3])+permute(tensorprod(eye(2),tensorprod(eye(2),F)),[2,4,1,6,5,3])+permute(tensorprod(eye(2),tensorprod(eye(2),F)),[2,6,1,4,5,3])+permute(tensorprod(eye(2),tensorprod(eye(2),F)),[2,6,5,3,1,4])+permute(tensorprod(eye(2),tensorprod(eye(2),F)),[2,4,5,3,1,6]));
end

% fundamental path

A=koi.assembleFu(u0);
b=koi.assembleFlambda(2,"pressure");

tol=1e-9;
maxit=1e9;
M=ichol(A);
U0hat=pcg(A,-b,tol,maxit,M,M');

u0hat=Function(V);
u0hat.fromFreeDof(U0hat);
u0hat.computeJacobians();

% critical load

B=koi.assembleFuu(u0,u0hat);
[VV,LL]=eigs(A,-B,1,'smallestabs');
Vc=VV(:,end);
lambdac=LL(end,end);

vc=Function(V);
vc.fromFreeDof(Vc);
vc.computeJacobians();

uc=Function(V);
uc.fromFreeDof(lambdac*U0hat);
uc.computeJacobians();

% bifurcated path

C=koi.assembleFuu(uc,vc);
lambda0dot=-0.5*(Vc'*C*Vc)/(Vc'*C*U0hat);

t=40;

lambda=lambdac+lambda0dot*t;
U=lambda*U0hat+Vc*t;

u=Function(V);
u.fromFreeDof(U);
u.computeJacobians();

% plot

figure(1);

func=@(S,F) norm(S,"fro");
style='none';

subplot(1,3,1);
patch("Faces",geo.triangles,"Vertices",geo.vertices,'FaceVertexCData',zeros(geo.numtriangles,1),'FaceColor','flat','LineStyle',style);
pbaspect([1,1,1]);
daspect([1,1,1]);
axis off;

subplot(1,3,2);
P=geo.vertices+reshape(uc.dof,[geo.numvertices,2]);
patch("Faces",geo.triangles,"Vertices",P,'FaceVertexCData',koi.computeTrianglesStress(uc,func),'FaceColor','flat','LineStyle',style);
pbaspect([1,1,1]);
daspect([1,1,1]);
axis off;

subplot(1,3,3);
P=geo.vertices+reshape(u.dof,[geo.numvertices,2]);
patch("Faces",geo.triangles,"Vertices",P,'FaceVertexCData',koi.computeTrianglesStress(u,func),'FaceColor','flat','LineStyle',style);
pbaspect([1,1,1]);
daspect([1,1,1]);
axis off;
