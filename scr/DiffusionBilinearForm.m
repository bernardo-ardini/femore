classdef DiffusionBilinearForm < handle
    properties
        functionSpace
        diffusionCoefficient        
    end

    methods
        function a=DiffusionBilinearForm(V)
            assert(V.fe=="P1");
            a.functionSpace=V;
        end

        function A=assembleLoc(a,e)
            V=a.functionSpace;
            geo=V.geo;

            X=[ones(3,1),geo.vertices(geo.triangles(e,:),:)];
            M=det(X)*inv(X)';
            Delta=0.5*abs(det(X));
            b=M(:,2);
            c=M(:,3);
            A=1/(4*Delta)*(b*b'+c*c');
            A=A(:);
            A=A*a.diffusionCoefficient;
        end

        function A=assemble(a,varargin)
            V=a.functionSpace;
            geo=V.geo;

            I=zeros(3*3,geo.numtriangles);
            J=zeros(3*3,geo.numtriangles);
            vals=zeros(3*3,geo.numtriangles);

            for e=1:geo.numtriangles
                vals(:,e)=a.assembleLoc(e);
                I(:,e)=repmat([geo.triangles(e,1);geo.triangles(e,2);geo.triangles(e,3)],3,1);
                J(:,e)=[repmat(geo.triangles(e,1),3,1);repmat(geo.triangles(e,2),3,1);repmat(geo.triangles(e,3),3,1)];
            end

            vals=vals(:);
            I=I(:);
            J=J(:);
            
            A=sparse(I,J,vals);

            if ~isempty(varargin) && varargin{1}~="AllDof"
                A=V.toFreeDof(A);
            end
        end
    end
end