classdef BasisNodes1d
    
    properties
        N;
        V;
        invV;
        
        r;
        
        basisEvaluator@BasisEvaluator1d;
        
        iVertices;
        iEdges;
        
        numNodes;
    end
    
    methods
        function obj = BasisNodes1d(N)
            obj.N = N;
            obj.basisEvaluator = BasisEvaluator1d(N);
            
            obj.r = support.nodes1d(N);
            %obj.V = support.vandermonde(N, obj.r);
            obj.V = obj.basisEvaluator.vandermonde(obj.r);
            obj.invV = inv(obj.V);
            
            obj.iVertices = [1, N];
            obj.iEdges = 2:(N-1);
            
            obj.numNodes = N;
        end
        
        
        % ---- DIFFERENTIATION
        
        function Dr = gradientMatrix(obj, rr)
            
            %dVdr = support.gradVandermonde(obj.N, rr);
            dVdr = obj.basisEvaluator.vandermondeDerivative(rr);
            
            Dr = dVdr * obj.invV;
        end
        
        % ---- INTERPOLATION
        
        function M = interpolationMatrix(obj, rr)
            %V2 = support.vandermonde(obj.N, rr);
            V2 = obj.basisEvaluator.vandermonde(rr);
            M = V2 * obj.invV;
        end
        
        % ---- NODE COORDINATES
        
        function rs = getNodes(obj, varargin)
            % Return row vector of nodes on the edge
            
            if nargin < 2
                orientation = 1;
            else
                orientation = varargin{1};
            end
            
            if orientation > 0
                rs = obj.r;
            else
                rs = obj.r(end:-1:1);
            end
        end
        
        function rs = getInteriorNodes(obj, varargin)
            % Return row vector of interior nodes on the edge, or empty
            % array for N <= 2.
            
            if nargin < 2
                orientation = 1;
            else
                orientation = varargin{1};
            end
            
            if obj.N > 2
                if orientation > 0
                    rs = obj.r(2:(obj.N-1));
                else
                    rs = obj.r((obj.N-1):-1:2);
                end
            else
                rs = [];
            end
        end
    end
    
    
end