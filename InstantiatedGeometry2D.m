classdef InstantiatedGeometry2D < handle
    
    properties
        parameterizedGeometry;
        
        contourVertexIndices;
        contourMeshSizes;
        vertices;
        lines;
        
        dvx_dp;
        dvy_dp;
        
        meshStruct;
        geom2mesh; % a matrix
        mesh2mesh; % a matrix
        
        tnMesh@TriNodalMesh;
        geomNodeJacobian;
        
        N_geom;
        N_field;
        N_quad;
    end
    
    methods
        
        function obj = InstantiatedGeometry2D(inParameterizedGeometry, inN_field, inN_geom, inN_quad)
            obj.parameterizedGeometry = inParameterizedGeometry;
            
            obj.meshStruct = [];
            obj.geom2mesh = []; % it's a Jacobian actually
            obj.mesh2mesh = [];
            
            obj.geomNodeJacobian = [];
            
            obj.N_field = inN_field;
            obj.N_geom = inN_geom;
            obj.N_quad = inN_quad;
        end
        
        function instantiateMesh(obj, paramVector)
            
            geometry = obj.parameterizedGeometry.evaluateGeometry(paramVector);
            [obj.dvx_dp, obj.dvy_dp] = obj.parameterizedGeometry.evaluateGeometryJacobian(paramVector);
            
            obj.contourVertexIndices = geometry.contourVertexIndices;
            obj.contourMeshSizes = geometry.contourMeshSizes;
            obj.vertices = geometry.vertices;
            obj.lines = geometry.lines;
            
            obj.createMesh();
            obj.geom2mesh = obj.getGeomToMeshMatrix();
            obj.mesh2mesh = obj.getMeshInteriorMatrix();
        end
        
        function adjustMesh(obj, paramVector)
            % Tweak meshStruct assuming that the topology is the same now.
            
            geometry = obj.parameterizedGeometry.evaluateGeometry(paramVector);
            [obj.dvx_dp, obj.dvy_dp] = obj.parameterizedGeometry.evaluateGeometryJacobian(paramVector);
            
            obj.contourVertexIndices = geometry.contourVertexIndices;
            obj.contourMeshSizes = geometry.contourMeshSizes;
            oldVertices = obj.vertices;
            obj.vertices = geometry.vertices;
            obj.lines = geometry.lines;
            
            if ~isempty(obj.meshStruct)
                idxBdyVerts = obj.meshStruct.boundaryEdges(:,1);
                idxInteriorVerts = setdiff(1:size(obj.meshStruct.vertices,1), idxBdyVerts);
                
                obj.meshStruct.vertices(idxBdyVerts,:) = obj.geom2mesh * obj.vertices; % sole difference from instantiateMesh
                
                obj.tnMesh.xyNodes = obj.tnMesh.xyNodes + obj.geomNodeJacobian * (obj.vertices - oldVertices);
                
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
            jac = sparse(numBoundaryEdges, size(obj.vertices,1));
            
            for ee = 1:numBoundaryEdges
                idxMeshVert0 = obj.meshStruct.boundaryEdges(ee,1); % to avoid redundancy, each edge contributes only one vertex
                idxEdgeLine = obj.meshStruct.edgeGeometryLines(ee);
                idxGeomVerts = obj.lines(idxEdgeLine,:);
                
                xy0 = obj.vertices(idxGeomVerts(1),:);
                xy1 = obj.vertices(idxGeomVerts(2),:);
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
        
        
        function outFieldNodeContour = getFieldNodeContours(obj)
            % Get the contour index for each field node.  Only boundary
            % nodes will have a contour.
            
            numFieldNodes = obj.tnMesh.hFieldNodes.getNumNodes(); % local
            outFieldNodeContour = sparse(numFieldNodes, 1); % Needed for boundary conditions!!
            for cc = 1:length(obj.contourVertexIndices)
                iContourVertices = obj.meshStruct.boundaryEdges(obj.meshStruct.edgeGeometryContours == cc,:); % local
                iContourEdges = obj.tnMesh.hMesh.getVertexEdgesExclusive(unique(iContourVertices(:))); % local
                iContourNodes = obj.tnMesh.hFieldNodes.getEdgeNodes(iContourEdges); % local
                outFieldNodeContour(iContourNodes) = cc;
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
                iContourNodes = obj.tnMesh.hGeomNodes.getEdgeNodes(iLineEdges);
                outGeomNodeLine(iContourNodes) = ll;
            end
            
        end
        
        function jac = getMeshGeometryJacobians(obj)
            
            jac = sparse(obj.tnMesh.hGeomNodes.getNumNodes(), size(obj.vertices, 1));
            
            geomNodeLine = obj.getGeomNodeLines();
            
            for ll = 1:size(obj.lines,1)
                
                iGeomNode = find(geomNodeLine == ll);  % indices of geometry nodes on this line
                %plot(xyGeomNodes(ii,1), xyGeomNodes(ii,2), 'o');
                
                iLineVert0 = obj.lines(ll,1);
                iLineVert1 = obj.lines(ll,2);
                xy0 = obj.vertices(iLineVert0,:);
                xy1 = obj.vertices(iLineVert1,:);
                
                % Node positions = xy0 + alpha*xy1
                % Let's do this a lazy way
                alpha = [xy1' - xy0'] \ [ obj.tnMesh.xyNodes(iGeomNode,1)' - xy0(1); obj.tnMesh.xyNodes(iGeomNode,2)' - xy0(2) ];
                
                jac(iGeomNode, iLineVert0) = 1 - alpha;
                jac(iGeomNode, iLineVert1) = alpha;
            end
            
        end
        
        
        function createMesh(obj)
            
            contourVertices = {};
            for cc = 1:length(obj.contourVertexIndices)
                contourVertices{cc} = obj.vertices(obj.contourVertexIndices{cc},:);
            end
            
            writeGEO('fromMatlab.geo', contourVertices, obj.contourMeshSizes);
            [status, result] = unix('/usr/local/bin/gmsh -2 fromMatlab.geo > gmshOut.txt');
            [mshFaceVertices, mshEdgeVertices, mshVerts, mshEdgeContour, mshEdgeLine] = readMSH('fromMatlab.msh');
            
            obj.meshStruct = struct('faces', mshFaceVertices,...
                'boundaryEdges', mshEdgeVertices, ...
                'vertices', mshVerts, ...
                'edgeGeometryContours', mshEdgeContour, ...
                'edgeGeometryLines', mshEdgeLine);
            
            lng = LinearNodalGeometry(obj.meshStruct.faces, obj.meshStruct.vertices, obj.N_geom);
            obj.tnMesh = TriNodalMesh(obj.meshStruct.faces, lng.getNodeCoordinates(), obj.N_field, obj.N_geom, obj.N_quad);
            
            obj.geomNodeJacobian = obj.getMeshGeometryJacobians();
        end
        
        
        function plotMesh(obj, varargin)
            if isempty(obj.meshStruct)
                error('Mesh has not been created');
            end
            
            patch('Faces', obj.meshStruct.faces, 'Vertices', obj.meshStruct.vertices, varargin{:});
        end
            
            
        
    end
    
end


