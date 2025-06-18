classdef StokesProblem < handle
    properties
        V
        Q
        mu
        c
        delta
    end

    methods
        function sp=StokesProblem(V,Q)
            assert(V.fe=="P12b" || V.fe=="P12");
            assert(Q.fe=="P1");
            assert(V.geo==Q.geo);
            sp.V=V;
            sp.Q=Q;
        end

        function [localA,localB,localC]=assembleLocal(sp,e)
            geo=sp.V.geo;
            Finv=geo.inverseJacobian{e};
            area=geo.areas(e);

            if sp.V.fe=="P12b"
                referenceAssemblyP12b;

                localA=2*area*(sp.mu*tensorprod(D,Finv*Finv',[3,4],[1,2])+sp.c*E);
                localB=-2*area*(tensorprod(C,Finv',3,2));

                localA=localA(:);
                localB=localB(:);
                localC=[];
            % elseif sp.V.fe=="P12"
            %     referenceAssemblyP12;
            % 
            %     h=max(geo.lengths(abs(geo.triangles2edges(e,:))));
            % 
            %     localA=2*area*(sp.mu*tensorprod(D,Finv*Finv',[3,4],[1,2])+sp.c*E);
            %     localB=-2*area*(tensorprod(C,Finv',3,2));
            %     localC=-2*area*h^2*sp.delta*tensorprod(D,Finv*Finv',[3,4],[1,2]);
            % 
            %     localA=localA(:);
            %     localB=localB(:);
            %     localC=localC(:);
            % end
            elseif sp.V.fe=="P12"
                referenceAssemblyP12;

                h=max(geo.lengths(abs(geo.triangles2edges(e,:))));

                localA=2*area*(sp.mu*tensorprod(D,Finv*Finv',[3,4],[1,2])+sp.c*E-sp.delta*h^2*sp.c^2*E);
                localB=-2*area*(tensorprod(C,Finv',3,2)+sp.delta*h^2*sp.c*permute(tensorprod(C,Finv',3,2),[2,1,3]));
                localC=-2*area*h^2*sp.delta*tensorprod(D,Finv*Finv',[3,4],[1,2]);

                localA=localA(:);
                localB=localB(:);
                localC=localC(:);
            end
        end

        function [A,B,C]=assemble(sp)
            geo=sp.V.geo;

            if sp.V.fe=="P12b"
                IA=zeros(4*4,geo.numtriangles);
                JA=zeros(4*4,geo.numtriangles);
                valsA=zeros(4*4,geo.numtriangles);
    
                IB=zeros(4*2*3,geo.numtriangles);
                JB=zeros(4*2*3,geo.numtriangles);
                valsB=zeros(4*2*3,geo.numtriangles);
            elseif sp.V.fe=="P12"
                IA=zeros(3*3,geo.numtriangles);
                JA=zeros(3*3,geo.numtriangles);
                valsA=zeros(3*3,geo.numtriangles);
    
                IB=zeros(3*2*3,geo.numtriangles);
                JB=zeros(3*2*3,geo.numtriangles);
                valsB=zeros(3*2*3,geo.numtriangles);

                IC=zeros(3*3,geo.numtriangles);
                JC=zeros(3*3,geo.numtriangles);
                valsC=zeros(3*3,geo.numtriangles);
            end

            for e=1:geo.numtriangles                
                if sp.V.fe=="P12b"
                    [valsA(:,e),valsB(:,e),~]=sp.assembleLocal(e);

                    [X,Y]=meshgrid([geo.triangles(e,1);geo.triangles(e,2);geo.triangles(e,3);geo.numvertices+e],[geo.triangles(e,1);geo.triangles(e,2);geo.triangles(e,3);geo.numvertices+e]);
                    IA(:,e)=Y(:);
                    JA(:,e)=X(:);

                    shift=geo.numvertices+geo.numtriangles;
                    [X,Y]=meshgrid([geo.triangles(e,1);geo.triangles(e,2);geo.triangles(e,3);geo.numvertices+e;shift+geo.triangles(e,1);shift+geo.triangles(e,2);shift+geo.triangles(e,3);shift+geo.numvertices+e],[geo.triangles(e,1);geo.triangles(e,2);geo.triangles(e,3)]);
                    IB(:,e)=Y(:);
                    JB(:,e)=X(:);
                elseif sp.V.fe=="P12"
                    [valsA(:,e),valsB(:,e),valsC(:,e)]=sp.assembleLocal(e);

                    [X,Y]=meshgrid([geo.triangles(e,1);geo.triangles(e,2);geo.triangles(e,3)],[geo.triangles(e,1);geo.triangles(e,2);geo.triangles(e,3)]);
                    IA(:,e)=Y(:);
                    JA(:,e)=X(:);

                    shift=geo.numvertices;
                    [X,Y]=meshgrid([geo.triangles(e,1);geo.triangles(e,2);geo.triangles(e,3);shift+geo.triangles(e,1);shift+geo.triangles(e,2);shift+geo.triangles(e,3)],[geo.triangles(e,1);geo.triangles(e,2);geo.triangles(e,3)]);
                    IB(:,e)=Y(:);
                    JB(:,e)=X(:);

                    [X,Y]=meshgrid([geo.triangles(e,1);geo.triangles(e,2);geo.triangles(e,3)],[geo.triangles(e,1);geo.triangles(e,2);geo.triangles(e,3)]);
                    IC(:,e)=Y(:);
                    JC(:,e)=X(:);
                end
            end

            valsA=valsA(:);
            IA=IA(:);
            JA=JA(:);
            A=sparse(IA,JA,valsA);
            A=blkdiag(A,A);

            valsB=valsB(:);
            IB=IB(:);
            JB=JB(:);
            B=sparse(IB,JB,valsB);
            B(sp.Q.constrainedVertices,:)=[];

            if sp.V.fe=="P12b"
                C=sparse(geo.numvertices-1,geo.numvertices-1);
            elseif sp.V.fe=="P12"
                valsC=valsC(:);
                IC=IC(:);
                JC=JC(:);
                C=sparse(IC,JC,valsC);
                C(sp.Q.constrainedVertices,:)=[];
                C(:,sp.Q.constrainedVertices)=[];
            end
        end

        function [l,m]=assemblerhs(sp,f,r,varargin)
            geo=sp.V.geo;

            p=inputParser;
            addParameter(p,"mode","interpolation");
            parse(p,varargin{:});

            if sp.V.fe=="P1"
                assert(p.Results.mode=="gauss","If GLS stabilization is used the only quadrature mode is gauss");
            end

            if p.Results.mode=="interpolation"
                F=Function(sp.V);
                F.fromFunctionHandle(f);
                l=sp.V.massMatrix*F.dof;

                R=Function(sp.Q);
                R.fromFunctionHandle(r);
                m=-sp.Q.massMatrix*R.dof;
            elseif p.Results.mode=="gauss"
                gaussquad=GaussQuad(5);
    
                if sp.V.fe=="P12b"
                    IV=zeros(4*2,geo.numtriangles);
                    valsV=zeros(4*2,geo.numtriangles);
                elseif sp.V.fe=="P12"
                    IV=zeros(3*2,geo.numtriangles);
                    valsV=zeros(3*2,geo.numtriangles);
                end
    
                IQ=zeros(3,geo.numtriangles);
                valsQ=zeros(3,geo.numtriangles);
    
                for e=1:geo.numtriangles
                    affineTransformation=geo.affineTransformation{e};
                    F=affineTransformation(:,1:2);
                    t=affineTransformation(:,3);
                    h=max(geo.lengths(abs(geo.triangles2edges(e,:))));
                    Finv=geo.inverseJacobian{e};

                    if sp.V.fe=="P12b"
                        int=@(x) [f(F*x+t).*(1-x(1,:)-x(2,:));f(F*x+t).*x(1,:);f(F*x+t).*x(2,:);f(F*x+t).*(27*(1-x(1,:)-x(2,:)).*x(1,:).*x(2,:))];
                        valsV(:,e)=2*geo.areas(e)*gaussquad.integral(int);
                        shift=geo.numvertices+geo.numtriangles;
                        IV(:,e)=[geo.triangles(e,1),geo.triangles(e,1)+shift,geo.triangles(e,2),geo.triangles(e,2)+shift,geo.triangles(e,3),geo.triangles(e,3)+shift,geo.numvertices+e,shift+geo.numvertices+e];

                        int=@(x) [r(F*x+t).*(1-x(1,:)-x(2,:));r(F*x+t).*x(1,:);r(F*x+t).*x(2,:)];
                        valsQ(:,e)=-2*geo.areas(e)*gaussquad.integral(int);
                        IQ(:,e)=[geo.triangles(e,1),geo.triangles(e,2),geo.triangles(e,3)];
                    elseif sp.V.fe=="P12"
                        int=@(x) [f(F*x+t).*(1-x(1,:)-x(2,:));f(F*x+t).*x(1,:);f(F*x+t).*x(2,:)];
                        % valsV(:,e)=2*geo.areas(e)*gaussquad.integral(int);
                        valsV(:,e)=2*geo.areas(e)*(1-sp.delta*h^2*sp.c)*gaussquad.integral(int);
                        shift=geo.numvertices;
                        IV(:,e)=[geo.triangles(e,1),geo.triangles(e,1)+shift,geo.triangles(e,2),geo.triangles(e,2)+shift,geo.triangles(e,3),geo.triangles(e,3)+shift];

                        int=@(x) [r(F*x+t).*(1-x(1,:)-x(2,:))-sp.delta*h^2*diag(f(F*x+t)'*Finv'*([-1;-1]*ones(1,size(x,2))))';r(F*x+t).*x(1,:)-sp.delta*h^2*diag(f(F*x+t)'*Finv'*([1;0]*ones(1,size(x,2))))';r(F*x+t).*x(2,:)-sp.delta*h^2*diag(f(F*x+t)'*Finv'*([0;1]*ones(1,size(x,2))))'];
                        valsQ(:,e)=-2*geo.areas(e)*gaussquad.integral(int);
                        IQ(:,e)=[geo.triangles(e,1),geo.triangles(e,2),geo.triangles(e,3)];
                    end
                end
    
                valsV=valsV(:);
                IV=IV(:);
                l=accumarray(IV,valsV);
    
                valsQ=valsQ(:);
                IQ=IQ(:);
                m=accumarray(IQ,valsQ);
            end
        end

        function [u,p]=solve(sp,g,f,r,varargin)
			parser=inputParser;
			addOptional(parser,'PressurePreconditioner','n',@isstring);
		    parse(parser,varargin{:});

            [A,B,C]=sp.assemble();
     
            G=g(sp.V.geo.vertices(sp.V.constrainedVertices,:)');
            G=[G(1,:)';G(2,:)'];
            G=sp.V.fromConstrainedDof(G);

            [l,m]=sp.assemblerhs(f,r,"mode","gauss");
            l=l-A*G;
            l=sp.V.toFreeDof(l);

            m=sp.Q.toFreeDof(m);
            m=m-B*G;

            A=sp.V.toFreeDof(A);
            B(:,sp.V.vertex2index(sp.V.constrainedVertices))=[];

            opts.type="ict";
            opts.droptol=1e-2;
            L=ichol(A,opts);

            b=B*pcg(A,l,1e-8,500,L,L')-m;

            function v=R(x)
                v=B*pcg(A,B'*x,1e-9,200,L,L')-C*x;
            end

			if parser.Results.PressurePreconditioner=="n"
            	P=pcg(@R,b,1e-7,800);
			elseif parser.Results.PressurePreconditioner=="y"
				P=pcg(@R,b,1e-7,800,B*B');
			end

            bb=l-B'*P;
            U=pcg(A,bb,1e-8,500,L,L');

            p=Function(sp.Q);
            p.fromFreeDof(P);

            u=Function(sp.V);
            u.fromFreeDof(U);
            u.dof=u.dof+G;
        end
    end
end