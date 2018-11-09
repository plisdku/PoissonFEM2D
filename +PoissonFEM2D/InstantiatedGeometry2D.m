classdef InstantiatedGeometry2D < handle
    
    properties
        parameterizedGeometry;
        
        geometry;
        
        dvx_dp;
        dvy_dp;
        
        meshStruct;
        geom2mesh; % a matrix
        mesh2mesh; % a matrix
        
        tnMesh%@TriNodalMesh;
        geomNodeJacobian;
        
        N_geom;
        N_field;
        N_quad;
        
        isAxisymmetric;
    end
    
    methods
        
        function obj = InstantiatedGeometry2D(inParameterizedGeometry, inN_field, inN_geom, inN_quad, inIsAxisymmetric)
            obj.parameterizedGeometry = inParameterizedGeometry;
            obj.geometry = [];
            
            obj.meshStruct = [];
            obj.geom2mesh = []; % it's a Jacobian actually
            obj.mesh2mesh = [];
            
            obj.geomNodeJacobian = [];
            
            obj.N_field = inN_field;
            obj.N_geom = inN_geom;
            obj.N_quad = inN_quad;
            
            if nargin < 5
                obj.isAxisymmetric = 0;
            else
                obj.isAxisymmetric = inIsAxisymmetric;
            end
        end
        
        function instantiateMesh(obj, paramVector)
            
            obj.geometry = obj.parameterizedGeometry.evaluateGeometry(paramVector);
            %obj.plotGeometry();
            if ~isempty(paramVector)
                [obj.dvx_dp, obj.dvy_dp] = obj.parameterizedGeometry.evaluateGeometryJacobian(paramVector);
            else
                obj.dvx_dp = sparse(length(obj.geometry.vertices), 0);
                obj.dvy_dp = sparse(length(obj.geometry.vertices), 0);
            end
            
            obj.createMesh();
            obj.geom2mesh = obj.getGeomToMeshMatrix();
            obj.mesh2mesh = obj.getMeshInteriorMatrix();
        end
        
        function adjustMesh(obj, paramVector)
            % Tweak meshStruct assuming that the topology is the same now.
            
            oldVertices = obj.geometry.vertices;
            obj.geometry = obj.parameterizedGeometry.evaluateGeometry(paramVector);
            [obj.dvx_dp, obj.dvy_dp] = obj.parameterizedGeometry.evaluateGeometryJacobian(paramVector);
            
            if ~isempty(obj.meshStruct)
                idxBdyVerts = obj.meshStruct.boundaryEdges(:,1);
                idxInteriorVerts = setdiff(1:size(obj.meshStruct.vertices,1), idxBdyVerts);
                
                obj.meshStruct.vertices(idxBdyVerts,:) = obj.geom2mesh * obj.geometry.vertices; % sole difference from instantiateMesh
                
                obj.tnMesh.xyNodes = obj.tnMesh.xyNodes + obj.geomNodeJacobian * (obj.geometry.vertices - oldVertices);
                
                doInterior = 0;
                if doInterior
                    % Interior-interior matrix
                    A_ii = speye(length(idxInteriorVerts)) - obj.mesh2mesh(idxInteriorVerts, idxInteriorVerts);

                    % Interior-boundary matrix
                    A_ib = obj.mesh2mesh(idxInteriorVerts, idxBdyVerts);

                    obj.meshStruct.vertices(idxInteriorVerts,:) = A_ii \ (A_ib * obj.meshStruct.vertices(idxBdyVerts,:));                
                end
            end
        end
        
        function jac = getGeomToMeshMatrix(obj)
            % Decompose each vertex in obj.meshStruct as a linear
            % combination of geometry vertices.
            %
            % that is, boundaryMeshVertices = A * geomVertices
            
            numBoundaryEdges = size(obj.meshStruct.boundaryEdges, 1);
            jac = sparse(numBoundaryEdges, size(obj.geometry.vertices,1));
            
            for ee = 1:numBoundaryEdges
                idxMeshVert0 = obj.meshStruct.boundaryEdges(ee,1); % to avoid redundancy, each edge contributes only one vertex
                idxEdgeLine = obj.meshStruct.edgeGeometryLines(ee);
                idxGeomVerts = obj.geometry.lines(idxEdgeLine,:);
                
                xy0 = obj.geometry.vertices(idxGeomVerts(1),:);
                xy1 = obj.geometry.vertices(idxGeomVerts(2),:);
                alpha = [xy1' - xy0'] \ (obj.meshStruct.vertices(idxMeshVert0,:)' - xy0');
                
                jac(ee, idxGeomVerts(1)) = 1 - alpha;
                jac(ee, idxGeomVerts(2)) = alpha;
            end
        end
        
        function outMatrix = getMeshInteriorMatrix(obj)
            % meshVertices = A * meshVertices.  uhhhh yeah.  Weird huh.
            
            numVertices = size(obj.meshStruct.vertices, 1);
            outMatrix = sparse(numVertices, numVertices);
            
            vv = VVMesh.fv2vv(obj.meshStruct.faces);
            vertAdjMatrix = VVMesh.vv2adjacency(vv);
            
            for idxVert = 1:numVertices
                
                idxNeighbors = find(vertAdjMatrix(idxVert,:));
                
                thisVertex = obj.meshStruct.vertices(idxVert,:);
                neighborVertices = obj.meshStruct.vertices(idxNeighbors,:);
                
                neighborDisplacements = bsxfun(@minus, neighborVertices, thisVertex);
                neighborDistances = sqrt(sum(neighborDisplacements.^2, 2));
                neighborWeights = neighborDistances / sum(neighborDistances);
                
                outMatrix(idxVert, idxNeighbors) = neighborWeights;
            end
            
        end
        
        function outFieldNodeLabels = getFieldNodeLabels(obj)
            % Last-added labels win when a node is on two lines!!!
            
            %fieldNodeLines = obj.getFieldNodeLines();
            %[row,col,val] = find(fieldNodeLines);
            
            labelVals = obj.geometry.lineLabels(val);
            
            numFieldNodes = obj.tnMesh.hFieldNodes.getNumNodes();
            outFieldNodeLine = sparse(numFieldNodes, 1);
            numLines = numel(unique(obj.meshStruct.edgeGeometryLines));
            for ll = 1:numLines
                iLineVertices = obj.meshStruct.boundaryEdges(obj.meshStruct.edgeGeometryLines == ll,:);
                iLineEdges = obj.tnMesh.hMesh.getVertexEdgesExclusive(unique(iLineVertices(:)));
                iLineNodes = obj.tnMesh.hFieldNodes.getEdgeNodes(iLineEdges);
                outFieldNodeLine(iLineNodes) = ll;
            end
            
            outFieldNodeLabels = sparse(row, col, labelVals);
        end
        
        function outFieldNodeContour = getFieldNodeContours(obj)
            % Get the contour index for each field node.  Only boundary
            % nodes will have a contour.
            
            numFieldNodes = obj.tnMesh.hFieldNodes.getNumNodes(); % local
            outFieldNodeContour = sparse(numFieldNodes, 1); % Needed for boundary conditions!!
            for cc = 1:length(obj.geometry.contourVertexIndices)
                iContourVertices = obj.meshStruct.boundaryEdges(obj.meshStruct.edgeGeometryContours == cc,:); % local
                iContourEdges = obj.tnMesh.hMesh.getVertexEdgesExclusive(unique(iContourVertices(:))); % local
                iContourNodes = obj.tnMesh.hFieldNodes.getEdgeNodes(iContourEdges); % local
                outFieldNodeContour(iContourNodes) = cc;
            end
        end
        
        function iNodes = getLabeledFieldNodes(obj, label)
            iLines = find(obj.geometry.lineLabels == label);
            iNodes = obj.getLineFieldNodes(iLines);
        end
        
        function iNodes = getLineFieldNodes(obj, iLines)
            iLineVertices = obj.meshStruct.boundaryEdges(ismember(obj.meshStruct.edgeGeometryLines, iLines),:);
            iLineEdges = obj.tnMesh.hMesh.getVertexEdgesExclusive(unique(iLineVertices(:)));
            iNodes = obj.tnMesh.hFieldNodes.getEdgeNodes(iLineEdges);
        end
        
        function outFieldNodeLine = getFieldNodeLines(obj)
            % Get the line index for each field node.  Only boundary nodes
            % will have a line.
            
            numFieldNodes = obj.tnMesh.hFieldNodes.getNumNodes();
            outFieldNodeLine = sparse(numFieldNodes, 1);
            numLines = numel(unique(obj.meshStruct.edgeGeometryLines));
            for ll = 1:numLines
                iLineVertices = obj.meshStruct.boundaryEdges(obj.meshStruct.edgeGeometryLines == ll,:);
                iLineEdges = obj.tnMesh.hMesh.getVertexEdgesExclusive(unique(iLineVertices(:)));
                iLineNodes = obj.tnMesh.hFieldNodes.getEdgeNodes(iLineEdges);
                outFieldNodeLine(iLineNodes) = ll;
            end
        end
        
        function outGeomNodeLine = getGeomNodeLines(obj)
            % Get the line index for each geometry node.  Only boundary
            % nodes will have a line.
            
            numGeomNodes = obj.tnMesh.hGeomNodes.getNumNodes();
            outGeomNodeLine = sparse(numGeomNodes, 1);
            numLines = numel(unique(obj.meshStruct.edgeGeometryLines));
            for ll = 1:numLines
                iLineVertices = obj.meshStruct.boundaryEdges(obj.meshStruct.edgeGeometryLines == ll,:);
                iLineEdges = obj.tnMesh.hMesh.getVertexEdgesExclusive(unique(iLineVertices(:)));
                iLineNodes = obj.tnMesh.hGeomNodes.getEdgeNodes(iLineEdges);
                outGeomNodeLine(iLineNodes) = ll;
            end
            
        end
        
        function jac = getMeshGeometryJacobians(obj)
            
            jac = sparse(obj.tnMesh.hGeomNodes.getNumNodes(), size(obj.geometry.vertices, 1));
            
            geomNodeLine = obj.getGeomNodeLines();
            
            for ll = 1:size(obj.geometry.lines,1)
                
                iGeomNode = find(geomNodeLine == ll);  % indices of geometry nodes on this line
                %plot(xyGeomNodes(ii,1), xyGeomNodes(ii,2), 'o');
                
                iLineVert0 = obj.geometry.lines(ll,1);
                iLineVert1 = obj.geometry.lines(ll,2);
                xy0 = obj.geometry.vertices(iLineVert0,:);
                xy1 = obj.geometry.vertices(iLineVert1,:);
                
                % Node positions = xy0 + alpha*xy1
                % Let's do this a lazy way
                alpha = [xy1' - xy0'] \ [ obj.tnMesh.xyNodes(iGeomNode,1)' - xy0(1); obj.tnMesh.xyNodes(iGeomNode,2)' - xy0(2) ];
                
                jac(iGeomNode, iLineVert0) = 1 - alpha;
                jac(iGeomNode, iLineVert1) = alpha;
            end
            
        end
        
        
        function createMesh(obj)
            
            contourVertices = {};
            for cc = 1:length(obj.geometry.contourVertexIndices)
                contourVertices{cc} = obj.geometry.vertices(obj.geometry.contourVertexIndices{cc},:);
            end
            i_error = 0;
             while(1)
%                  try
                    PoissonFEM2D.writeGEO('fromMatlab.geo', contourVertices, obj.geometry.contourMeshSizes);
                    gmshPath = PoissonFEM2D.getGmshPath();
                    cmd = sprintf('%s -2 fromMatlab.geo > gmshOut.txt', gmshPath);
%                     path = pwd;
%                     path = strrep(path,'\','/');
%                     path = strrep(path,' ','\ ');
%                     file = sprintf('%s/fromMatlab.geo',path);
%                     system(sprintf('C:\\cygwin64\\bin\\bash --login -c " /cygdrive/c/Users/larstn/Desktop/gmsh-4.0.2-Windows64/gmsh-4.0.2-Windows64/gmsh -2 -format msh2 %s > gmshOut.txt"', file));
%                     %[status, result] = unix('/usr/local/bin/gmsh -2 fromMatlab.geo > gmshOut.txt');
                    [status, result] = unix(cmd);
                    if status
%                        contourMeshSizes = obj.geometry.contourMeshSizes;
%                        save('GmshWritingErrorData','contourVertices', 'contourMeshSizes')
%                         
%                        warning('An error occured while WRITING the gmsh file')
%                        warning('please refer to file GmshwritingErrorData.m for the data causing this error')
                    %while status  
                        warning(result);
                        PoissonFEM2D.writeGEO('fromMatlab.geo', contourVertices, obj.geometry.contourMeshSizes);
                        [status, result] = unix(cmd);                        
                    else
                        
                    %end
                        break
%                 catch e
%                     pause(0.5)
%                     fprintf('%s \n', e.identifier)
%                     warning('Matlab fopen error occured, cleaning up...')
%                     if exist('fromMatlab.geo','file') == 2
%                         delete fromMatlab.geo
%                     end
%                     pause(0.2)
%                     i_error = i_error + 1;
%                     if i_error > 20
%                         keyboard
%                     end
                    end
                
            end
            
            [mshFaceVertices, mshEdgeVertices, mshVerts, mshEdgeContour, mshEdgeLine] = PoissonFEM2D.readMSH('fromMatlab.msh');
            
            obj.meshStruct = struct('faces', mshFaceVertices,...
                'boundaryEdges', mshEdgeVertices, ...
                'vertices', mshVerts, ...
                'edgeGeometryContours', mshEdgeContour, ...
                'edgeGeometryLines', mshEdgeLine);
            
            lng = PoissonFEM2D.LinearNodalGeometry(obj.meshStruct.faces, obj.meshStruct.vertices, obj.N_geom);
            obj.tnMesh = PoissonFEM2D.TriNodalMesh(obj.meshStruct.faces, lng.getNodeCoordinates(), obj.N_field, obj.N_geom, obj.N_quad, obj.isAxisymmetric);
            
            obj.geomNodeJacobian = obj.getMeshGeometryJacobians();
        end
        
        
        function plotMesh(obj, varargin)
            if isempty(obj.meshStruct)
                error('Mesh has not been created');
            end
            
            patch('Faces', obj.meshStruct.faces, 'Vertices', obj.meshStruct.vertices, varargin{:});
        end
        
        function plotGeometry(obj)
            geometry = obj.geometry();
            plot([geometry.vertices(geometry.lines(:,1),1), geometry.vertices(geometry.lines(:,2),1)]', ...
                [geometry.vertices(geometry.lines(:,1),2), geometry.vertices(geometry.lines(:,2),2)]', 'color', [0.8 0.8 0.8], 'linewidth', 2)
        end
        
    end
    
end


