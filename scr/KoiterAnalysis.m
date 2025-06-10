classdef KoiterAnalysis < handle
    properties
        functionSpace % fornire un function space P12 con le condizioni di Dirichlet

        W
        S
        D2W
        D3W
    end

    methods
        function koi=KoiterAnalysis(V)
            assert(V.fe=="P12");
            koi.functionSpace=V;
        end

        function stress=computeTrianglesStress(koi,u,func)
            V=koi.functionSpace;
            geo=V.geo;

            stress=zeros(geo.numtriangles,1);
            
            for e=1:geo.numtriangles
                Du=u.jacobian{e};
                stress(e)=stress(e)+func(koi.S(eye(2)+Du),eye(2)+Du);
            end
        end

        function stress=computeVerticesStress(koi,u,func)
            V=koi.functionSpace;
            geo=V.geo;

            normalization=zeros(geo.numvertices,1);
            stress=zeros(geo.numvertices,1);
            
            for e=1:geo.numtriangles
                Du=u.jacobian{e};
                
                for a=1:3
                    aa=geo.triangles(e,a);

                    normalization(aa)=normalization(aa)+geo.areas(e);
                    stress(aa)=stress(aa)+geo.areas(e)*func(koi.S(eye(2)+Du),eye(2)+Du);
                end
            end

            stress=stress./normalization;
        end

        function Flambda=assembleFlambda(koi,indices,load)
            V=koi.functionSpace;
            geo=V.geo;
            d=geo.d;

            edg=geo.lines2edges(indices);

            if size(edg,1)==0
                Flambda=zeros(V.numberFreeDof(),1);
                return
            end

            I=zeros(d*(d+1),size(edg,1));
            vals=zeros(d*(d+1),size(edg,1));

            cc=1;

            for b=edg'
                e=geo.edges2triangles(b,2);
                mask=abs(geo.triangles2edges(e,:))==b;
                bb=1:3;
                bb=bb(mask);
                
                nu=sign(geo.triangles2edges(e,bb))*geo.normals(b,:)';

                l=geo.lengths(b);
                v=(geo.vertices(geo.edges(b,1),:)+geo.vertices(geo.edges(b,2),:))/2;
                B=V.B{e};

                if isstring(load) && load=="pressure"
                    val=l*nu*[v,1]*B';
                else
                    val=-l*load*[v,1]*B';
                end
                
                val=val';
                vals(:,cc)=val(:);

                c=1;
                for i=1:d
                    for a=1:d+1
                        I(c,cc)=(i-1)*geo.numvertices+geo.triangles(e,a);
                        c=c+1;
                    end
                end

                cc=cc+1;
            end

            vals=vals(:);
            I=I(:);

            vals=[vals;zeros(2*geo.numvertices,1)];
            I=[I;(1:2*geo.numvertices)'];

            Flambda=accumarray(I,vals);
            Flambda=V.toFreeDof(Flambda);
        end

        function Fu=assembleFu(koi,u)
            % matrice associata a v|->Fu(u,l\lambda)v

            V=koi.functionSpace;
            geo=V.geo;
            d=geo.d;

            I=zeros(d^2*(d+1)^2,geo.numtriangles);
            J=zeros(d^2*(d+1)^2,geo.numtriangles);
            vals=zeros(d^2*(d+1)^2,geo.numtriangles);

            for e=1:geo.numtriangles
                Du=u.jacobian{e};

                B=koi.functionSpace.B{e};
                B=B(:,1:d);

                val=geo.areas(e)*permute(tensorprod(tensorprod(koi.D2W(eye(d)+Du),B,2,2),B,3,2),[3,1,4,2]);
                vals(:,e)=val(:);

                c=1;
                for j=1:d
                    for b=1:d+1
                        for i=1:d
                            for a=1:d+1
                                I(c,e)=(i-1)*geo.numvertices+geo.triangles(e,a);
                                J(c,e)=(j-1)*geo.numvertices+geo.triangles(e,b);
                                c=c+1;
                            end
                        end
                    end
                end
            end

            vals=vals(:);
            I=I(:);
            J=J(:);

            Fu=sparse(I,J,vals);
            Fu=V.toFreeDof(Fu);
        end

        function Fuu=assembleFuu(koi,u,w)
            % matrice associata a v|->Fuu(u,l\lambda)wv

            V=koi.functionSpace;
            geo=V.geo;
            d=geo.d;

            I=zeros(d^2*(d+1)^2,geo.numtriangles);
            J=zeros(d^2*(d+1)^2,geo.numtriangles);
            vals=zeros(d^2*(d+1)^2,geo.numtriangles);

            for e=1:geo.numtriangles
                Du=u.jacobian{e};
                Dw=w.jacobian{e};

                B=koi.functionSpace.B{e};
                B=B(:,1:d);

                val=geo.areas(e)*permute(tensorprod(tensorprod(tensorprod(koi.D3W(eye(d)+Du),Dw,[1,2],[1,2]),B,2,2),B,3,2),[3,1,4,2]);
                vals(:,e)=val(:);

                c=1;
                for j=1:d
                    for b=1:d+1
                        for i=1:d
                            for a=1:d+1
                                I(c,e)=(i-1)*geo.numvertices+geo.triangles(e,a);
                                J(c,e)=(j-1)*geo.numvertices+geo.triangles(e,b);
                                c=c+1;
                            end
                        end
                    end
                end
            end

            vals=vals(:);
            I=I(:);
            J=J(:);

            Fuu=sparse(I,J,vals);
            Fuu=V.toFreeDof(Fuu);
        end
    end
end