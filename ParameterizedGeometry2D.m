classdef ParameterizedGeometry2D < handle
    
    properties
        contours;
        %geometry;
        %dvx_dp;
        %dvy_dp;
    end
    
    
    methods (Static)
        function bc = dirichlet(value, lineIds)
            bc = ParameterizedGeometry2D.boundaryCondition('dirichlet', value, lineIds);
        end
        
        function bc = neumann(value, lineIds)
            bc = ParameterizedGeometry2D.boundaryCondition('neumann', value, lineIds);
        end
        
        function bc = boundaryCondition(bcType, value, lineIds)
            assert(isa(value, 'function_handle') || isscalar(value), ...
                'Boundary value must be a scalar or a function handle of (p,x,y)');
            if isscalar(value)
                value = @(p,x,y) value;
            end
            bc = struct('type', bcType, 'func', value, 'relativeLineIds', lineIds);
        end
    end
    
    methods
        
        function obj = ParameterizedGeometry2D()
            obj.contours = [];
        end
        
        function addContour(obj, x, y, meshSize, varargin)
            
            assert(isa(x, 'function_handle'), 'x must be a function handle returning vertex x coordinates');
            assert(isa(y, 'function_handle'), 'y must be a function handle returning vertex y coordinates');
            assert(isa(meshSize, 'numeric'), 'meshSize must be a scalar or an array of one mesh size per vertex');
            
            lineLabels = [];
            for nn = 1:2:length(varargin)
                labelId = varargin{nn};
                lineIdsInContour = varargin{nn+1};
                maxLineId = max(lineIdsInContour);
                
                if length(lineLabels) < maxLineId
                    tmp = zeros(1, maxLineId);
                    tmp(1:length(lineLabels)) = lineLabels;
                    lineLabels = tmp;
                end
                
                lineLabels(lineIdsInContour) = labelId;
            end
            
            % Duh, I should put the labels HERE
            contourStruct = struct('xFunc', x, 'yFunc', y, 'meshSize', meshSize, 'lineLabels', lineLabels);
            
            if isempty(obj.contours)
                obj.contours = contourStruct;
            else
                obj.contours(end+1) = contourStruct;
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
            outGeometryLineLabels = [];
            idxGeomVert = 1;
            
            contourLengths = zeros(size(obj.contours));
            
            for cc = 1:length(obj.contours)
                x = obj.contours(cc).xFunc(p); %local
                y = obj.contours(cc).yFunc(p); %local
                contourLength = length(x); %local
                numLineLabels = length(obj.contours(cc).lineLabels);
                
                if numLineLabels ~= contourLength
                    error(sprintf('Contour %i has %i line labels for %i lines\n', cc, numLineLabels, contourLength));
                end
                
                contourLengths(cc) = contourLength;
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
                outGeometryLineLabels = [outGeometryLineLabels; reshape(obj.contours(cc).lineLabels, [], 1)];
                idxGeomVert = idxGeomVert + contourLength; % local
            end
            
            g = struct('contourVertexIndices', {outContourVertexIndices}, ...
                'contourMeshSizes', {outMeshSizes}, ...
                'vertices', outGeometryVertices, ...
                'lines', outGeometryLines, ...
                'lineLabels', outGeometryLineLabels);
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