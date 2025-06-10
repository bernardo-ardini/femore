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

        function w=toFreeDof(u)
            w=u.functionSpace.toFreeDof(u.dof);
        end

        function fromFreeDof(u,w)
            V=u.functionSpace;
            geo=V.geo;

            if V.fe=="P12"
                M=eye(2*geo.numvertices);
                M(:,[V.constrainedVertices();geo.numvertices+V.constrainedVertices()])=[];
            elseif V.fe=="P1"
                M=eye(geo.numvertices);
                M(:,V.constrainedVertices())=[];
            end

            u.dof=M*w;
        end

        function fromConstrainedDof(u,w)
            V=u.functionSpace;
            geo=V.geo;

            if V.fe=="P12"
                M=eye(2*geo.numvertices);
                M=M(:,[V.constrainedVertices();geo.numvertices+V.constrainedVertices()]);
            elseif V.fe=="P1"
                M=eye(geo.numvertices);
                M=M(:,V.constrainedVertices());
            end

            u.dof=M*w;
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