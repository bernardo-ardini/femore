classdef StokesProblem < handle
    properties
        functionSpace
        mu
        c
    end

    methods
        function sp=StokesProblem(V)
            sp.functionSpace=V;
        end

        function [localA,localB]=assembleLocal(sp,e)
            V=sp.functionSpace;
            geo=V.geo;

            % computation of A
            D=[-1,-1;1,0;0,1];
            M=1/24*[2,1,1;1,2,1;1,1,2];
            m=[1/360;1/360;1/360;1/5040];

            localA=sp.mu*geo.areas(e)*D*geo.metricTensor{e}*D'+sp.mu*2*geo.areas(e)*M;
            a=sp.mu*2*geo.areas(e)*m;
            localA=[localA,a(1:3);a(1:3)',a(4)];

            localA=localA(:);

            % computation of B
            L=zeros(3,4,2);
            L(1,1,:)=[-1/6,-1/6];
            L(1,2,:)=[1/6,0];
            L(1,3,:)=[0,1/6];
            L(1,4,:)=[1/120,1/120];
            L(2,1,:)=[-1/6,-1/6];
            L(2,2,:)=[1/6,0];
            L(2,3,:)=[0,1/6];
            L(2,4,:)=[-1/120,0];
            L(3,1,:)=[-1/6,-1/6];
            L(3,2,:)=[1/6,0];
            L(3,3,:)=[0,1/6];
            L(3,4,:)=[0,-1/120];

            localB=-2*geo.areas(e)*tensorprod(L,geo.inverseJacobian{e},3,1);

            localB=localB(:);
        end

        function [A,B]=assemble(sp)
            V=sp.functionSpace;
            geo=V.geo;

            IA=zeros(4*4,geo.numtriangles);
            JA=zeros(4*4,geo.numtriangles);
            valsA=zeros(4*4,geo.numtriangles);

            IB=zeros(4*2*3,geo.numtriangles);
            JB=zeros(4*2*3,geo.numtriangles);
            valsB=zeros(4*2*3,geo.numtriangles);

            for e=1:geo.numtriangles
                [valsA(:,e),valsB(:,e)]=sp.assembleLocal(e);
                IA(:,e)=repmat([geo.triangles(e,1);geo.triangles(e,2);geo.triangles(e,3);geo.numvertices+e],4,1);
                JA(:,e)=[repmat(geo.triangles(e,1),4,1);repmat(geo.triangles(e,2),4,1);repmat(geo.triangles(e,3),4,1);repmat(geo.numvertices+e,4,1)];
                IB(:,e)=repmat([geo.triangles(e,1);geo.triangles(e,2);geo.triangles(e,3)],8,1);
                shift=geo.numvertices+geo.numtriangles;
                JB(:,e)=[repmat(geo.triangles(e,1),3,1);repmat(geo.triangles(e,2),3,1);repmat(geo.triangles(e,3),3,1);repmat(geo.numvertices+e,3,1);...
                    repmat(shift+geo.triangles(e,1),3,1);repmat(shift+geo.triangles(e,2),3,1);repmat(shift+geo.triangles(e,3),3,1);repmat(shift+geo.numvertices+e,3,1)];
            end

            valsA=valsA(:);
            IA=IA(:);
            JA=JA(:);

            valsB=valsB(:);
            IB=IB(:);
            JB=JB(:);
            
            A=sparse(IA,JA,valsA);
            B=sparse(IB,JB,valsB);
        end
    end
end