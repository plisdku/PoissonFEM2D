classdef MeshTopology < handle
% MeshTopology Topology of faces, edges, vertices and FEM nodes in 2D
% triangulated mesh.
%
% The mesh geometry is not handled here.
    
    properties (Access = private)
        N;  % nodes per edge of element
        
        basis;  % BasisNodes object
        
        interiorNodesPerEdge;  % number of interior nodes per edge (vertices excluded)
        interiorNodesPerFace;  % number of interior nodes per face (edges and corners excluded)
        
        % Mesh topology and geometry.
        faceVertices;
        faceEdges;
        faceEdgeOrientations; % +1 or -1
        edgeVertices;  % not oriented.
        
        numVertices; % inferred from the mesh
        
        % Nodes are numbered consecutively in three groups: vertex nodes,
        % edge nodes and face nodes.
        %   iVertexNode0 is the index beginning the vertex nodes
        %   iEdgeNode0 is the index beginning the edge nodes
        %   iFaceNode0 is the index beginning the face nodes
        iVertexNode0;
        iEdgeNode0;
        iFaceNode0;
    end

    methods
        
        % ---- CONSTRUCTOR
        
        function obj = MeshTopology(faceVertices, N)
            obj.N = N;
            
            obj.basis = BasisNodes(N);
            
            obj.interiorNodesPerEdge = max(0, N-2);
            obj.interiorNodesPerFace = max(0, (N-3)*(N-2)/2);
            
            
            obj.faceVertices = faceVertices;
            [obj.faceEdges, obj.faceEdgeOrientations, obj.edgeVertices] = VVMesh.getEdges(obj.faceVertices);
            obj.numVertices = numel(unique(faceVertices(:)));
            
            
            obj.iVertexNode0 = 0;
            obj.iEdgeNode0 = obj.getNumVertices();
            obj.iFaceNode0 = obj.iEdgeNode0 + obj.interiorNodesPerEdge*obj.getNumEdges();
        end
        
        % ---- NUMBERS OF THINGS
        
        function n = getNodalOrder(obj)
            n = obj.N;
        end
        
        function n = getNumVertices(obj)
            n = obj.numVertices;
        end
        
        function n = getNumEdges(obj)
            n = size(obj.edgeVertices, 1);
        end

        function n = getNumFaces(obj)
            n = size(obj.faceVertices, 1);
        end

        function n = getNumNodes(obj)
            n = obj.getNumVertices() + ...
                obj.interiorNodesPerEdge*obj.getNumEdges() + ...
                obj.interiorNodesPerFace*obj.getNumFaces();
        end
        
        % ---- BASIC TOPOLOGY
        
        function [iEdges, orientations] = getFaceEdges(obj, iFaces)
            % Return ordered edge indices and orientations for each face
            %
            % getFaceEdges(iFaces)
            
            if nargin == 1
                iEdges = obj.faceEdges;
                orientations = obj.faceEdgeOrientations;
                return
            end
            
            iEdges = obj.faceEdges(iFaces,:);
            orientations = obj.faceEdgeOrientations(iFaces,:);
        end
        
        function iVertices = getFaceVertices(obj, iFaces)
            % Return ordered vertex indices for each face
            %
            % getFaceVertices(iFaces)
            
            if nargin == 1
                iVertices = obj.faceVertices;
                return
            end
            
            iVertices = obj.faceVertices(iFaces,:);
        end
        
        function iVertices = getEdgeVertices(obj, iEdges, varargin)
            % Return start and end vertex indices for each edge
            %
            % getEdgeVertices(iEdges)
            % getEdgeVertices(iEdges, orientations)
            
            if nargin == 1
                iVertices = obj.edgeVertices;
                return
            end
            
            iVertices = obj.edgeVertices(iEdges,:);
            
            if length(varargin) >= 1
                iFlippedEdges = varargin{1} < 0;
                iVertices(iFlippedEdges,:) = iVertices(iFlippedEdges, [2, 1]);
            end
        end
        
        % ---- ADJACENCY MATRICES
        
        function A = getFaceEdgeAdjacency(obj)
            % getFaceEdgeAdjacency()     indexed (face, edge)
            
            nFaces = obj.getNumFaces();
            nEdges = obj.getNumEdges();
            
            iFaces = repmat( (1:nFaces)', 1, 3);
            
            A = sparse(iFaces(:), obj.faceEdges(:), ones(3*nFaces,1), nFaces, nEdges);
        end
        
        function O = getFaceEdgeOrientation(obj)
            % getFaceEdgeOrientation()    indexed (face, edge)
            
            nFaces = obj.getNumFaces();
            nEdges = obj.getNumEdges();
            
            iFaces = repmat( (1:nFaces)', 1, 3);
            
            O = sparse(iFaces(:), obj.faceEdges(:), obj.faceEdgeOrientations(:), nFaces, nEdges);
        end
        
        function A = getFaceVertexAdjacency(obj)
            
            nFaces = obj.getNumFaces();
            nVertices = obj.getNumVertices();
            
            iFaces = repmat( (1:nFaces)', 1, 3);
            
            A = sparse(iFaces(:), obj.faceVertices(:), ones(3*nFaces,1), nFaces, nVertices);
        end
        
        function A = getEdgeVertexAdjacency(obj)
            
            nEdges = obj.getNumEdges();
            nVertices = obj.getNumVertices();
            
            iEdges = repmat( (1:nEdges)', 1, 2);
            
            A = sparse(iEdges(:), obj.edgeVertices(:), ones(2*nEdges,1), nEdges, nVertices);
        end
        
        % ---- NODE ACCESSORS
        
        
        function iNodes = getFaceNodes(obj, iFace)
            % getFaceNodes(iFace)   Get node indices for all nodes in a face.
            
            iVertexNodes = obj.getVertexNodes(obj.faceVertices(iFace,:));
            iEdgeNodes1 = obj.getFaceEdgeInteriorNodes(iFace, 1);
            iEdgeNodes2 = obj.getFaceEdgeInteriorNodes(iFace, 2);
            iEdgeNodes3 = obj.getFaceEdgeInteriorNodes(iFace, 3);
            iFaceNodes = obj.getFaceInteriorNodes(iFace);
            
            % Now we somehow need to slam them together in the right order.
            
            iNodes = zeros(obj.basis.numNodes, 1);
            iNodes(obj.basis.iVertices) = iVertexNodes;
            iNodes(obj.basis.iEdges(1,:)) = iEdgeNodes1;
            iNodes(obj.basis.iEdges(2,:)) = iEdgeNodes2;
            iNodes(obj.basis.iEdges(3,:)) = iEdgeNodes3;
            iNodes(obj.basis.iCenter) = iFaceNodes;
        end
        
        function iNodes = getFaceInteriorNodes(obj, iFace)
            assert(numel(iFace) == 1, 'Only one face at a time');
            iNodes = obj.iFaceNode0 + (iFace-1)*obj.interiorNodesPerFace + ...
                (1:obj.interiorNodesPerFace);
        end
        
        function iNodes = getFaceEdgeNodes(obj, iFace, iEdgeLocal)
            iEdgeGlobal = obj.faceEdges(iFace, iEdgeLocal);
            orientation = obj.faceEdgeOrientations(iFace, iEdgeLocal);
            iInteriorNodes = obj.getEdgeInteriorNodes(iEdgeGlobal, orientation);
            
            edgeVertsLocal = [iEdgeLocal, 1 + mod(iEdgeLocal,3)];
            iVertexNodes = obj.getVertexNodes(obj.faces(iFace, edgeVertsLocal));
            
            iNodes = [iVertexNodes(1), iInteriorNodes, iVertexNodes(2)];
        end
        
        function iNodes = getFaceEdgeInteriorNodes(obj, iFace, iEdgeLocal)
            iEdgeGlobal = obj.faceEdges(iFace, iEdgeLocal);
            orientation = obj.faceEdgeOrientations(iFace, iEdgeLocal);
            iNodes = obj.getEdgeInteriorNodes(iEdgeGlobal, orientation);
        end
        
        function iNodes = getFaceVertexNodes(obj, iFace, iVerticesLocal)
            iVertices = obj.faces(iFace, iVerticesLocal);
            iNodes = obj.getVertexNodes(iVertices);
        end
        
        function iNodes = getEdgeNodes(obj, iEdge, varargin)
            % getEdgeNodes(iEdge)
            % getEdgeNodes(iEdge, orientation)
            assert(numel(iEdge) == 1, 'Only one edge at a time');
            
            iInteriorNodes = obj.iEdgeNode0 + (iEdge-1)*obj.interiorNodesPerEdge + ...
                (1:obj.interiorNodesPerEdge);
            iVertexNodes = obj.getVertexNodes(obj.edgeVertices(iEdge,:));
            
            iNodes = [iVertexNodes(1), iInteriorNodes, iVertexNodes(2)];
            
            if ~isempty(varargin)
                if varargin{1} < 0
                    iNodes = iNodes(end:-1:1);
                end
            end
        end
        
        function iNodes = getEdgeInteriorNodes(obj, iEdge, varargin)
            % getEdgeInteriorNodes(iEdge)
            % getEdgeInteriorNodes(iEdge, orientation)
            assert(numel(iEdge) == 1, 'Only one edge at a time');
            iNodes = obj.iEdgeNode0 + (iEdge-1)*obj.interiorNodesPerEdge + ...
                (1:obj.interiorNodesPerEdge);
            
            if ~isempty(varargin)
               if varargin{1} < 0
                   iNodes = iNodes(end:-1:1);
               end
            end
        end
        
        function iNodes = getVertexNodes(obj, iVertices)
            iNodes = obj.iVertexNode0 + iVertices;
        end
        
        
        % ---- MESH TOPOLOGY
        
        function [edgeIndices, orientations] = getBoundaryEdges(obj)
            edgeUses = obj.faceEdges(:);
            edgeUseCounts = hist(edgeUses, 1:obj.getNumEdges());
            
            edgeIndices = find(edgeUseCounts == 1);
            
            % Get the faces for these edges, and the orientations they
            % use the edges with.
            
            O = obj.getFaceEdgeOrientation();
            [~, ii, ori] = find(O(:,edgeIndices));
            
            % Select the orientations corresponding to the boundary edges.
            orientations = 0*edgeIndices;
            orientations(ii) = ori;
        end
        
        % Get an arbitrarily-ordered list of the nodes on the outer boundary
        function iNodes = getBoundaryNodes(obj)
            
            iEdgeIndices = obj.getBoundaryEdges();
            iEdgeVertexIndices = obj.edgeVertices(iEdgeIndices,:);
            
            iNodes = zeros(numel(iEdgeIndices), obj.interiorNodesPerEdge);
            for ii = 1:length(iEdgeIndices)
                iNodes(ii,:) = obj.getEdgeInteriorNodes(iEdgeIndices(ii));
            end
            iNodes = unique(iNodes(:));

            % throw in the edge vertices
            iEdgeCornerNodes = unique(obj.getVertexNodes(iEdgeVertexIndices(:)));

            iNodes = [iNodes; iEdgeCornerNodes];
            
        end
        
        function iCenterNodes = getInteriorNodes(obj)
            
            iCenterNodes = setdiff(1:obj.getNumNodes(), obj.getBoundaryNodes());
            
        end
        
        
    end % methods


end