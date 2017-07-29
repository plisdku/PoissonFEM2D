classdef BasisNodes1d
    
    properties
        N;
        V;
        invV;
        
        r;
        
        iVertices;
        iEdges;
        
        numNodes;
    end
    
    methods
        function obj = BasisNodes1d(N)
            obj.N = N;
            
            obj.V = support.vandermonde(N);
            obj.invV = inv(obj.V);
            
            obj.r = support.nodes1d(N);
            
            obj.iVertices = [1, N];
            obj.iEdges = 2:(N-1);
            
            obj.numNodes = N;
        end
        
        function rs = getNodes(obj)
            % Return row vector of nodes on the edge
            rs = obj.r;
        end
        
        function rs = getInteriorNodes(obj)
            % Return row vector of interior nodes on the edge, or empty
            % array for N <= 2.
            if obj.N > 2
                rs = obj.r(2:(obj.N-1));
            else
                rs = [];
            end
        end
    end
    
    
end