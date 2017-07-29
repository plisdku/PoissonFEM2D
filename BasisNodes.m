classdef BasisNodes
    
    properties
        N;   % number of nodes on one edge
        V;   % vandermonde matrix
        invV; % inverse of vandermonde matrix
        r;  % r-coordinates
        s;  % s-coordinates
        
        iVertices;   % node numbers at the 3 corners
        iEdges;      % node numbers along the 3 edges, oriented.  Size [3, m]
        iCenter;     % node numbers in the center
        
        numNodes;
    end
    
    
    methods
        function obj = BasisNodes(N)
            obj.N = N;
            
            obj.V = support2d.vandermonde(N);
            obj.invV = inv(obj.V);
            
            [obj.r, obj.s] = support2d.nodes2d(N);
            
            [obj.iVertices, iEdges, obj.iCenter] = support2d.classifyNodes(N);
            obj.iEdges = [iEdges{1}; iEdges{2}; iEdges{3}];
            
            obj.numNodes = N*(N+1)/2;
        end
        
        % ---- INTERPOLATION
        
        function M = interpolationMatrix(obj, rr, ss)
            V2 = support2d.vandermonde(obj.N, rr, ss);
            M = V2 * obj.invV;
        end
        
        %function outVals = interpolate(obj, vals, rr, ss)
        %function outVals = interpolate(obj, vals, triVerts, xx, yy)
        function outVals = interpolate(obj, vals, varargin)
            
            if nargin == 4
                rr = varargin{1};
                ss = varargin{2};
                M = obj.interpolationMatrix(rr,ss);
                outVals = M*vals;
            elseif nargin == 5
                triVerts = varargin{1};
                
                xx = varargin{2};
                yy = varargin{3};
                
                xy = [reshape(xx,[],1), reshape(yy,[],1)];
                rs = support2d.xy2rs(triVerts, xy);
                
                M = obj.interpolationMatrix(rs(:,1), rs(:,2));
                outVals = M*vals;
            else
                error('shit');
            end
                
        end
        
        function outVals = scatteredInterp(obj, vals, rr, ss)
            
            interp = scatteredInterpolant(obj.r, obj.s, vals);
            outVals = interp(rr, ss);
        end
        
        % ---- NODE COORDINATES
        
        function rs = getNodes(obj)
            rs = [obj.r, obj.s];
        end
        
        function rs = getInteriorNodes(obj)
            rs = [obj.r(obj.iCenter), obj.s(obj.iCenter)];
        end
        
        function rs = getEdgeInteriorNodes(obj, iEdge)
            rs = [obj.r(obj.iEdges(iEdge,:)), obj.s(obj.iEdges(iEdge,:))];
        end
        
    end
    
    
end