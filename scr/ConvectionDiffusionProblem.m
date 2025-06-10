classdef ConvectionDiffusionProblem < handle
    properties
        functionSpace
        diffusivity
        velocity
        tau
    end

    methods
        function a=ConvectionDiffusionProblem(V)
            assert(V.fe=="P1");
            a.functionSpace=V;
        end

        function A=assembleLocalStiffness(prob,e)
            V=prob.functionSpace;
            geo=V.geo;

            x=geo.vertices(geo.triangles(e,:),:);
            F=[x(2,:)-x(1,:);x(3,:)-x(1,:)]';
            c=x(1,:)';

            centroid=geo.centroids(e,:)';

            B=V.B{e};
            b=B(:,end);
            B=B(:,1:2);
            
            diff=prob.diffusivity(centroid);
            vel=prob.velocity(centroid);

            h=max(geo.lengths(abs(geo.triangles2edges(e,:))));

            A=geo.areas(e)*(B*diff*B'+(B*(1/3*F*[1;1]+c)+b)*vel'*B');

            if norm(vel)~=0
                A=A+geo.areas(e)*prob.tau*h/(norm(vel)*norm(diff(:),Inf))*B*(vel*vel')*B';
            end

            A=A(:);
        end

        function A=assembleStiffness(a)
            V=a.functionSpace;
            geo=V.geo;

            I=zeros(3*3,geo.numtriangles);
            J=zeros(3*3,geo.numtriangles);
            vals=zeros(3*3,geo.numtriangles);

            for e=1:geo.numtriangles
                vals(:,e)=a.assembleLocalStiffness(e);
                I(:,e)=repmat([geo.triangles(e,1);geo.triangles(e,2);geo.triangles(e,3)],3,1);
                J(:,e)=[repmat(geo.triangles(e,1),3,1);repmat(geo.triangles(e,2),3,1);repmat(geo.triangles(e,3),3,1)];
            end

            vals=vals(:);
            I=I(:);
            J=J(:);
            
            A=sparse(I,J,vals);
        end
    end
end