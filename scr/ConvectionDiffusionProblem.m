classdef ConvectionDiffusionProblem < handle
    properties
        functionSpace
        diffusivity
        velocity
        delta
        xi

        stabilizationParameter
        PechletNumber

        typeStabilization

        theta
        rho
        epsilon

        A
    end

    methods
        function cdp=ConvectionDiffusionProblem(V)
            assert(V.fe=="P1");
            cdp.functionSpace=V;
            cdp.typeStabilization="Assignment";
        end

        function localA=assembleLocal(cdp,e)
            V=cdp.functionSpace;
            geo=V.geo;

            h=max(geo.lengths(abs(geo.triangles2edges(e,:))));

            B=V.B{e};
            B=B(:,1:2);

            alpha=min(eigs(cdp.diffusivity));
            Pe=h*norm(cdp.velocity)/alpha;

            if cdp.typeStabilization=="Assignment"
                gamma=cdp.delta*h/alpha/norm(cdp.velocity);
            elseif cdp.typeStabilization=="Quarteroni"
                gamma=0.5*h/norm(cdp.velocity)*cdp.xi(0.5*Pe);
            elseif cdp.typeStabilization=="My"
                b=norm(cdp.velocity);
                gamma=max([0,h/b*(1-cdp.theta)*(cdp.epsilon+cdp.theta/pi^2/cdp.epsilon-(2-cdp.rho-cdp.theta/cdp.rho)/Pe)]);
            end

            cdp.stabilizationParameter=max([gamma,cdp.stabilizationParameter]);
            cdp.PechletNumber=max([cdp.PechletNumber,Pe]);

            localA=geo.areas(e)*(B*(cdp.diffusivity+gamma*cdp.velocity*cdp.velocity')*B'+1/3*ones(3,1)*cdp.velocity'*B');
            localA=localA(:);
        end

        function assemble(cdp)
            V=cdp.functionSpace;
            geo=V.geo;

            I=zeros(3*3,geo.numtriangles);
            J=zeros(3*3,geo.numtriangles);
            vals=zeros(3*3,geo.numtriangles);

            for e=1:geo.numtriangles
                vals(:,e)=cdp.assembleLocal(e);
                I(:,e)=repmat([geo.triangles(e,1);geo.triangles(e,2);geo.triangles(e,3)],3,1);
                J(:,e)=[repmat(geo.triangles(e,1),3,1);repmat(geo.triangles(e,2),3,1);repmat(geo.triangles(e,3),3,1)];
            end

            vals=vals(:);
            I=I(:);
            J=J(:);
            
            cdp.A=sparse(I,J,vals);
        end

        function u=solve(cdp,dirval)
            % function space
            V=cdp.functionSpace;

            % assemble full stiffness matrix
            cdp.assemble();
            
            % compute all nodal values of lift function
            G=V.fromConstrainedDof(dirval);
            
            % compute rhs
            f=-cdp.A*G;
            f=V.toFreeDof(f);
            
            % compute the stiffness matrix related to free dof only
            AA=V.toFreeDof(cdp.A);
            
            % linear system solution
            restart=10;
            tol=1e-9;
            maxit=20;
            setup.type='nofill';
            setup.milu='off';
            [L,U]=ilu(AA,setup);
            W=gmres(AA,f,restart,tol,maxit,L,U);
            
            % final solution
            u=Function(V);
            u.fromFreeDof(W);
            u.dof=u.dof+G; % add lift function
        end

        function M=massFlux(cdp,u)
            M=sum(-cdp.A*u.dof);

            geo=cdp.functionSpace.geo;

            geo.boundaryEdges();

            for k=geo.boundaryEdges()'
                e=geo.edges2triangles(k,2);
                kk=1:3;
                kk=kk(abs(geo.triangles2edges(e,:))==k);
                nu=geo.normals(k,:)'*sign(geo.triangles2edges(e,kk));

                M=M+cdp.velocity'*nu*geo.lengths(k)/2*sum(u.dof(geo.edges2verticesList(k)));
            end

            M=M+cdp.velocity'*nu*geo.lengths(k)/2*sum(u.dof(geo.edges2verticesList(k)));
        end
    end
end