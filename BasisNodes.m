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
        
        function M = interpolationMatrix_rs(obj, rr, ss)
            V2 = support2d.vandermonde(obj.N, rr, ss);
            M = V2 * obj.invV;
        end
        
        function M = interpolationMatrix_xy(obj, triVerts, xx, yy)
            xy = [reshape(xx,1,[]); reshape(yy,1,[])];
            rs = support2d.xy2rs(triVerts, xy);
            V2 = support2d.vandermonde(obj.N, rs(1,:), rs(2,:));
            M = V2 * obj.invV;
        end
        
        function DM = interpolationMatrixSensivitity_nodal_xy(obj, triVerts, DtriVerts, xx, yy)
            % Works for interpolating nodal fields (values stuck to nodes)
            % but not for interpolating fields in xy space.  For that,
            % need the full derivative D(V2*inv(V)) instead of this which
            % is just D(V2)*inv(V).
            
            xy = [reshape(xx,1,[]); reshape(yy,1,[])];
            rs_query = support2d.xy2rs(triVerts, xy);
            Drs = support2d.xy2rsSensitivity(triVerts, DtriVerts, xy);
            [~, dVdr, dVds] = support2d.vandermonde(obj.N, rs_query(1,:), rs_query(2,:));
            DV = bsxfun(@times, dVdr, Drs(1,:)') + bsxfun(@times, dVds, Drs(2,:)');
            DM = DV*obj.invV;
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
                
                M = obj.interpolationMatrix_xy(triVerts, xx, yy);
                outVals = M*vals;
                
                %xy = [reshape(xx,1,[]); reshape(yy,1,[])];
                %rs = support2d.xy2rs(triVerts, xy);
                %M = obj.interpolationMatrix(rs(1,:), rs(2,:));
                %outVals = M*vals;
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