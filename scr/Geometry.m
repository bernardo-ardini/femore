classdef Geometry < handle
    properties
        % descrizione della triangolazione
        vertices % coordinate dei vertici
        edges % vertici dei latiindices
        triangles % vertici dei triangoli
        triangles2edges % lati dei triangoli (con segno)
        edges2triangles % triangoli che condividono un vertice (ordinati, 0 è l'esterno)
        normals % normali ai lati (sono uscenti da un triangolo se il segno corrispondente in trangles2edges è positivo)
        vertices2edges % lati corrispondenti ad una coppia di vertici (con segno, 0 se i vertici non sono collegati)
        numvertices % numero vertici
        numedges % numero lati
        numtriangles % numero triangoli

        areas % area degli elementi
        lengths % lunghezze dei lati
        centroids % baricentri degli elementi

        lines % matrice contenente in ogni riga gli indici di due vertici e un indice che dice a quale
        % delle linee quella coppia di vertici appartiene

        d % dimensione

        affineTransformation
        inverseJacobian
    end

    methods
        function geo=Geometry()
            geo.d=2;
        end

        function edg=boundaryEdges(geo)
            mask=geo.edges2triangles(:,1)==0;
            edg=1:geo.numedges;
            edg=edg(mask)';
        end

        function edg=lines2edges(geo,indices)
            pcv=sort(geo.lines(ismember(geo.lines(:,3),indices),1:2),2);
            edg=full(diag(geo.vertices2edges(pcv(:,1),pcv(:,2))));
        end

        function vert=lines2vertices(geo,indices)
            % restituisce il vettore con i vertici che appartengono alle
            % linee con indici indices
            vert=geo.lines(ismember(geo.lines(:,3),indices),1:2);
            vert=vert(:);
            vert=unique(vert);
        end

        function initialize(geo)
            % calcola tutto quello che c'è da calcolare a partire da solo
            % le coordinate dei vertici e i triangoli

            % salva le dimensioni
            geo.numvertices=size(geo.vertices,1);
            geo.numtriangles=size(geo.triangles,1);
            geo.numedges=geo.numvertices+geo.numtriangles-1;

            % prealloca
            geo.areas=zeros(geo.numtriangles,1);
            geo.vertices2edges=spalloc(geo.numvertices,geo.numvertices,2*geo.numedges);
            geo.edges=zeros(geo.numedges,2);
            geo.triangles2edges=zeros(geo.numtriangles,3);
            geo.edges2triangles=zeros(geo.numedges,2);
            geo.normals=zeros(geo.numedges,2);
            geo.lengths=zeros(geo.numedges,1);
            geo.centroids=zeros(geo.numtriangles,2);
            geo.affineTransformation=cell(geo.numtriangles,1);
            geo.inverseJacobian=cell(geo.numtriangles,1);
            
            % conteggio dei lati trovati
            index=0;

            for e=1:geo.numtriangles % faccio il giro dei triangoli
                % calcola la trasformazione affine di ciascun triangolo rispetto al triangolo di riferimento
                x=geo.vertices(geo.triangles(e,:),:);
                F=[x(2,:)-x(1,:);x(3,:)-x(1,:)]';
                geo.affineTransformation{e}=[F,x(1,:)'];
                geo.inverseJacobian{e}=inv(F);

                % calcolo l'area dei triangoli e i baricentri
                geo.areas(e)=0.5*abs(det(F));
                geo.centroids(e,:)=sum(x)/3;

                for alpha=1:3 % faccio il giro dei vertici dentro il triangolo
                    % trovo i due vertici a1 e a2 di un lato del triangolo
                    % e anche il terzo vertice a3
                    a1=geo.triangles(e,alpha);
                    a2=[2,3,1];
                    a2=a2(alpha);
                    a2=geo.triangles(e,a2);
                    a3=[3,1,2];
                    a3=a3(alpha);
                    a3=geo.triangles(e,a3);

                    % ordino i due vertici
                    aa1=min([a1,a2]);
                    aa2=max([a1,a2]);

                    if geo.vertices2edges(aa1,aa2)==0 % ho trovato un nuovo lato
                        index=index+1; % incremento il conteggio

                        geo.edges(index,:)=[aa1,aa2]; % aggiungo alla lista dei lati

                        % salvo nella matrice per trovare i lati dai
                        % vertici
                        geo.vertices2edges(aa1,aa2)=index;
                        geo.vertices2edges(aa2,aa1)=-index;

                        % calcolo i versori normali
                        orth=@(v) [v(2),-v(1)]; % ortogonale di un vettore
                        geo.normals(index,:)=orth(geo.vertices(aa2,:)-geo.vertices(aa1,:));
                        geo.lengths(index)=norm(geo.normals(index,:));
                        geo.normals(index,:)=geo.normals(index,:)/geo.lengths(index); % normalizzo

                        geo.edges2triangles(index,2)=e; % salvo nella matrice per trovare i triangoli dai lati
                    else % lato che avevo già trovato
                        % salvo nella matrice per trovare i triangoli dai
                        % lati ma li inverto per avere prima il triangolo
                        % minore
                        geo.edges2triangles(geo.vertices2edges(aa1,aa2),1)=geo.edges2triangles(geo.vertices2edges(aa1,aa2),2);
                        geo.edges2triangles(geo.vertices2edges(aa1,aa2),2)=e;
                    end

                    sigma=sign(geo.normals(abs(geo.vertices2edges(a1,a2)),:)*(0.5*(geo.vertices(a1,:)+geo.vertices(a2,:))-geo.vertices(a3,:))'); % calcolo se il versore è entrante o uscente
                    geo.triangles2edges(e,alpha)=sigma*abs(geo.vertices2edges(a1,a2)); % salvo nella matrice per trovare i lati di un elemento
                end
            end

            % ricalcolo il numero di edges se la formula di Eulero non
            % aveva dato il corretto numero perchè il dominio non era
            % semplicemente connesso
            geo.numedges=fix(nnz(geo.vertices2edges)/2);
        end
    end
end