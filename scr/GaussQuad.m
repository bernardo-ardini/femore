classdef GaussQuad < handle
    properties
        points
        weights
    end

    methods
        function gaussquad=GaussQuad(order)
            gaussquad.points=readintegrationdata(order);
            gaussquad.weights=gaussquad.points(:,3);
            gaussquad.points=gaussquad.points(:,1:2);
        end

        function I=integral(gaussquad,f,varargin)
            p=inputParser;
            addOptional(p,"typeIntegrand","handle");
            parse(p,varargin{:});
        
            if p.Results.typeIntegrand=="constant"
                I=0.5*f;
                return;
            end

            values=f(gaussquad.points');
            I=1/4*values*gaussquad.weights;
        end
    end
end