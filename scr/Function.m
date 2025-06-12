classdef Function < handle
    properties
        dof
        functionSpace

        jacobian
    end

    methods
        function u=Function(V)
            u.functionSpace=V;
            u.dof=zeros(V.numberDof(),1);
        end

        function fromFunctionHandle(u,f)
            if u.functionSpace.fe=="P1"
                u.dof=f(u.functionSpace.geo.vertices')';
            elseif u.functionSpace.fe=="P12"
                u1=f(u.functionSpace.geo.vertices');
                u2=u1(2,:);
                u1(2,:)=[];
                u.dof=[u1';u2'];
            elseif u.functionSpace.fe=="P12b"
                u1=f(u.functionSpace.geo.vertices');
                u2=u1(2,:);
                u1(2,:)=[];
                ub1=f(u.functionSpace.geo.centroids')-1/3*[sum(u1(u.functionSpace.geo.triangles),2)';sum(u2(u.functionSpace.geo.triangles),2)'];
                ub2=ub1(2,:);
                ub1(2,:)=[];
                u.dof=[u1';ub1';u2';ub2'];
            end
        end

        function w=toFreeDof(u)
            w=u.functionSpace.toFreeDof(u.dof);
        end

        function fromFreeDof(u,w)
            u.dof=u.functionSpace.fromFreeDof(w);
        end

        function fromConstrainedDof(u,w)
            u.dof=u.functionSpace.fromConstrainedDof(w);
        end

        function norm=normL2(u)
            M=u.functionSpace.massMatrix;
            norm=sqrt(u.dof'*M*u.dof);
        end

        function computeJacobians(u)
            V=u.functionSpace;
            geo=V.geo;
            d=geo.d;

            u.jacobian=cell(geo.numtriangles,1);

            if V.fe=="P12"
                for e=1:geo.numtriangles
                    if d==2
                        U=[u.dof(geo.triangles(e,:))';u.dof(geo.numvertices+geo.triangles(e,:))'];
                    elseif d==3
                        U=[u.dof(geo.triangles(e,:))';u.dof(geo.numvertices+geo.triangles(e,:))';u.dof(2*geo.numvertices+geo.triangles(e,:))'];
                    end

                    B=V.B{e};
                    B=B(:,1:d);
                    u.jacobian{e}=U*B;
                end
            end
        end
    end
end