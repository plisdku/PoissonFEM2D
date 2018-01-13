classdef ParameterizedGeometry2D < handle
    
    properties
        contours;
        %geometry;
        %dvx_dp;
        %dvy_dp;
    end
    
    
    methods
        
        function obj = ParameterizedGeometry2D()
            obj.contours = [];
        end
        
        function addContour(obj, x, y, meshSize, dirichletOrNeumann, boundaryFunction)
            
            assert(isa(x, 'function_handle'), 'x must be a function handle returning vertex x coordinates');
            assert(isa(y, 'function_handle'), 'y must be a function handle returning vertex y coordinates');
            assert(isa(meshSize, 'numeric'), 'meshSize must be a scalar or an array of one mesh size per vertex');
            assert(isa(boundaryFunction, 'function_handle'), 'boundaryFunction must be a function handle returning a scalar for each vertex');
            
            c = struct('xFunc', x, 'yFunc', y, 'meshSize', meshSize, ...
                'type', dirichletOrNeumann, 'boundaryFunc', boundaryFunction);
            
            if isempty(obj.contours)
                obj.contours = c;
            else
                obj.contours(end+1) = c;
            end
        end
        
        function instance = instantiate(obj, p)
            geometry = obj.evaluateGeometry(p);
            [dvx_dp, dvy_dp] = obj.evaluateGeometryJacobian(p);
            
            instance = InstantiatedGeometry2D(geometry, dvx_dp, dvy_dp);
        end
        
        function g = evaluateGeometry(obj, p)
            
            col = @(A) reshape(A, [], 1);
            
            outContourVertexIndices = {};
            outMeshSizes = {};
            
            outGeometryVertices = [];
            outGeometryLines = [];
            idxGeomVert = 1;
            
            for cc = 1:length(obj.contours)
                x = obj.contours(cc).xFunc(p); %local
                y = obj.contours(cc).yFunc(p); %local
                contourLength = length(x); %local
                meshSize = obj.contours(cc).meshSize; %local
                
                if length(meshSize) == 1
                    meshSize = repmat(meshSize, length(x), 1);
                end
                
                %outContourVertices{cc} = [col(x), col(y)];
                outMeshSizes{cc} = meshSize;
                
                outGeometryVertices = [outGeometryVertices; [col(x), col(y)]];
                
                arange = (0:(contourLength-1))'; %local
                
                outContourVertexIndices{cc} = idxGeomVert + arange;
                
                outGeometryLines = [outGeometryLines; idxGeomVert + [arange, mod(arange+1, contourLength)]];
                idxGeomVert = idxGeomVert + contourLength; % local
            end
            
            g = struct('contourVertexIndices', {outContourVertexIndices}, ...
                'contourMeshSizes', {outMeshSizes}, ...
                'vertices', outGeometryVertices, ...
                'lines', outGeometryLines);
        end
        
        
        function [dvxdp, dvydp] = evaluateGeometryJacobian(obj, p)
            
            numParams = length(p);
            
            rows = [];
            cols = [];
            vx_vals = [];
            vy_vals = [];
            
            delta = 1e-8;
            for nn = 1:numParams
                pLow = p;
                pLow(nn) = pLow(nn) - delta;
                
                pHigh = p;
                pHigh(nn) = pHigh(nn) + delta;
                
                vLow = obj.evaluateGeometry(pLow).vertices;
                vHigh = obj.evaluateGeometry(pHigh).vertices;
                
                numVertices = size(vLow,1);
                rows = [rows; (1:numVertices)'];
                cols = [cols; repmat(nn, numVertices, 1)];
                vx_vals = [vx_vals; (vHigh(:,1)-vLow(:,1))/(2*delta)];
                vy_vals = [vy_vals; (vHigh(:,2)-vLow(:,2))/(2*delta)];
            end
            
            dvxdp = sparse(rows, cols, vx_vals, numVertices, numParams);
            dvydp = sparse(rows, cols, vy_vals, numVertices, numParams);
        end
        
    end
    
    
    
    
end