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
            
            obj.r = support.nodes1d(N);
            obj.V = support.vandermonde(N, obj.r);
            obj.invV = inv(obj.V);
            
            
            obj.iVertices = [1, N];
            obj.iEdges = 2:(N-1);
            
            obj.numNodes = N;
        end
        
        
        % ---- INTERPOLATION
        
        function M = interpolationMatrix_r(obj, rr)
            V2 = support.vandermonde(obj.N, rr);
            M = V2 * obj.invV;
        end
        
        %function outVals = interpolate(obj, vals, rr)
        %%%%function outVals = interpolate(obj, vals, edgeVerts, xx)
        function outVals = interpolate(obj, vals, varargin)
            
            if nargin == 3
                rr = varargin{1};
                M = obj.interpolationMatrix_r(rr);
                outVals = M*vals;
            elseif nargin == 4
                error('unimplemented');
%                 triVerts = varargin{1};
%                 xx = varargin{2};
%                 
%                 M = obj.interpolationMatrix_xy(triVerts, xx, yy);
%                 outVals = M*vals;
            else
                error('shit');
            end
                
        end
        
        % ---- NODE COORDINATES
        
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