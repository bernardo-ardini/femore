classdef FunctionSpace < handle
    properties
        geo
        constrainedVertices
        linesConstraints
        fe
        B
    end

    methods
        function V=FunctionSpace(geo,fe)
            V.geo=geo;
            V.linesConstraints=[];
            V.fe=fe;
        end

        function ndof=numberDof(V)
            if V.fe=="P1"
                ndof=V.geo.numvertices;
            elseif V.fe=="P12"
                ndof=V.geo.d*V.geo.numvertices;
            elseif V.fe=="P12b"
                ndof=V.geo.d*(V.geo.numvertices+V.geo.numtriangles);
            end
        end

        function ndof=numberFreeDof(V)
            if V.fe=="P1"
                ndof=V.geo.numvertices-length(V.constrainedVertices());
            elseif V.fe=="P12"
                ndof=V.geo.d*(V.geo.numvertices-length(V.constrainedVertices()));
            elseif V.fe=="P12b"
                ndof=V.geo.d*(V.geo.numvertices+V.geo.numtriangles-length(V.constrainedVertices()));
            end
        end

        function ndof=numberConstrainedDof(V)
            ndof=V.numberDof()-V.numberFreeDof();
        end

        function setLinesConstraint(V,lines)
            V.linesConstraints=lines;
            V.constrainedVertices=[V.constrainedVertices;V.geo.lines2vertices(lines)];
        end

        function i=vertexComponent2index(V,in)
            if V.fe=="P1"
                i=in;
            elseif V.fe=="P12"
                i=(in(:,2)-1)*V.geo.numvertices+in(:,1);
            elseif V.fe=="P12b"
                i=(in(:,2)-1)*(V.geo.numvertices+V.geo.numtriangles)+in(:,1);
            end
        end

        function i=vertex2index(V,in)
            if V.fe=="P1"
                i=in;
            elseif V.fe=="P12"
                i=[in;V.geo.numvertices+in];
            elseif V.fe=="P12b"
                i=[in;V.geo.numvertices+V.geo.numtriangles+in];
            end
        end

        function out=toFreeDof(V,in)
            if isvector(in)
                out=in;
                out(V.vertex2index(V.constrainedVertices()))=[];
            elseif ismatrix(in)
                out=in;
                out(V.vertex2index(V.constrainedVertices()),:)=[];
                out(:,V.vertex2index(V.constrainedVertices()))=[];
            end
        end

        function out=fromConstrainedDof(V,in)
            M=speye(V.numberDof());
            M=M(:,V.vertex2index(V.constrainedVertices()));

            out=M*in;
        end

        function out=fromFreeDof(V,in)
            M=speye(V.numberDof());
            M(:,V.vertex2index(V.constrainedVertices()))=[];

            out=M*in;
        end

        function computeBasisFunctions(V)
            if ismember(V.fe,["P1","P12"])
                V.B=cell(V.geo.numtriangles,1);
                for e=1:V.geo.numtriangles
                    V.B{e}=inv([V.geo.vertices(V.geo.triangles(e,:),:)';ones(1,V.geo.d+1)]);
                end
            end
        end
    end
end