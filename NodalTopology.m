classdef NodalTopology < handle
    
    properties
        N; % num nodes per edge; edge polynomial order is N-1.
        
        basis@BasisNodes;
        basis1d@BasisNodes1d;
        
        numInteriorNodesPerEdge;
        numInteriorNodesPerFace;
        
        % Nodes are numbered consecutively in three groups: vertex nodes,
        % edge nodes and face nodes.
        %   iVertexNode0 is the index beginning the vertex nodes
        %   iEdgeNode0 is the index beginning the edge nodes
        %   iFaceNode0 is the index beginning the face nodes
        iVertexNode0;
        iEdgeNode0;
        iFaceNode0;
        
        hMesh@MeshTopology;  % handle to a MeshTopology which can be shared
    end
    
    
    methods
        % ---- CONSTRUCTOR
        
        function obj = NodalTopology(hMeshTopology, N)
            obj.hMesh = hMeshTopology;
            
            obj.N = N;
            obj.basis = BasisNodes(N);
            obj.basis1d = BasisNodes1d(N);
            
            obj.numInteriorNodesPerEdge = max(0, N-2);
            obj.numInteriorNodesPerFace = max(0, (N-3)*(N-2)/2);
            
            obj.iVertexNode0 = 0;
            obj.iEdgeNode0 = obj.hMesh.getNumVertices();
            obj.iFaceNode0 = obj.iEdgeNode0 + ...
                obj.numInteriorNodesPerEdge * obj.hMesh.getNumEdges();
        end
        
        % ---- NUMBERS OF THINGS
        
        function n = getNumNodes(obj)
            n = obj.hMesh.getNumVertices() + ...
                obj.numInteriorNodesPerEdge*obj.hMesh.getNumEdges() + ...
                obj.numInteriorNodesPerFace*obj.hMesh.getNumFaces();
        end
        
        
        % ---- NODE ACCESSORS
        
        function iNodes = getFaceNodes(obj, iFace)
            % getFaceNodes(iFace)   Get node indices for all nodes in a face.       
            
            iNodes = zeros(obj.basis.numNodes, 1);
            iNodes_idx = [obj.basis.iVertices, ...
                obj.basis.iEdges(1,:), ...
                obj.basis.iEdges(2,:), ...
                obj.basis.iEdges(3,:), ...
                obj.basis.iCenter];
            iNodes_val = ...
                [obj.getVertexNodes(obj.hMesh.faceVertices(iFace,:)),...
                obj.getFaceEdgeInteriorNodes(iFace, 1),...
                obj.getFaceEdgeInteriorNodes(iFace, 2),...
                obj.getFaceEdgeInteriorNodes(iFace, 3),...
                obj.getFaceInteriorNodes(iFace)];
            iNodes(iNodes_idx) = iNodes_val;
           
            
        end
        
        function iNodes = getFaceInteriorNodes(obj, iFace)
            assert(numel(iFace) == 1, 'Only one face at a time');
            iNodes = obj.iFaceNode0 + (iFace-1)*obj.numInteriorNodesPerFace + ...
                (1:obj.numInteriorNodesPerFace);
        end
        
        function iNodes = getFaceEdgeNodes(obj, iFace, iEdgeLocal)
            iEdgeGlobal = obj.hMesh.faceEdges(iFace, iEdgeLocal);
            orientation = obj.hMesh.faceEdgeOrientations(iFace, iEdgeLocal);
            iInteriorNodes = obj.getEdgeInteriorNodes(iEdgeGlobal, orientation);
            
            edgeVertsLocal = [iEdgeLocal, 1 + mod(iEdgeLocal,3)];
            iVertexNodes = obj.getVertexNodes(obj.hMesh.faceVertices(iFace, edgeVertsLocal));
            
            iNodes = [iVertexNodes(1), iInteriorNodes, iVertexNodes(2)];
        end
        
        function iNodes = getFaceEdgeInteriorNodes(obj, iFace, iEdgeLocal)
            iEdgeGlobal = obj.hMesh.faceEdges(iFace, iEdgeLocal);
            orientation = obj.hMesh.faceEdgeOrientations(iFace, iEdgeLocal);
            iNodes = obj.getEdgeInteriorNodes(iEdgeGlobal, orientation);
        end
        
        function iNodes = getFaceVertexNodes(obj, iFace, iVerticesLocal)
            if nargin < 3
                iVerticesLocal = [1;2;3];
            end
            iVertices = obj.hMesh.faceVertices(iFace, iVerticesLocal);
            iNodes = obj.getVertexNodes(iVertices);
        end
        
        function iNodes = getEdgeNodes(obj, iEdge, varargin)
            % getEdgeNodes(iEdge)
            % getEdgeNodes(iEdge, orientation)
            
            
            % MULTIPLE EDGE CASE
            if numel(iEdge) > 1
                numEdges = length(iEdge);
                
                allNodes = zeros(1, numEdges*obj.N); % overestimate
                iNextNode = 1;
                
                if ~isempty(varargin)
                    orientations = varargin{1};
                    if numel(orientations) ~= numel(iEdge)
                        orientations = repmat(orientations(1), numel(iEdge), 1);
                    end
                    %assert(numel(varargin{1}) == numel(iEdge), 'Must supply one orientation per edge');
                else
                    orientations = ones(numEdges,1);
                end
                
                for ii = 1:numEdges
                    allNodes(iNextNode:iNextNode+obj.N-1) = obj.getEdgeNodes(iEdge(ii),orientations(ii));
                    iNextNode = iNextNode + obj.N;
                end
                
                iNodes = unique(allNodes);
                
                return
            end
            
            
            % SINGLE EDGE CASE
            %assert(numel(iEdge) == 1, 'Only one edge at a time');
            
            iInteriorNodes = obj.iEdgeNode0 + (iEdge-1)*obj.numInteriorNodesPerEdge + ...
                (1:obj.numInteriorNodesPerEdge);
            iVertexNodes = obj.getVertexNodes(obj.hMesh.edgeVertices(iEdge,:));
            
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
            iNodes = obj.iEdgeNode0 + (iEdge-1)*obj.numInteriorNodesPerEdge + ...
                (1:obj.numInteriorNodesPerEdge);
            
            if ~isempty(varargin)
               if varargin{1} < 0
                   iNodes = iNodes(end:-1:1);
               end
            end
        end
        
        function iNodes = getVertexNodes(obj, iVertices)
            iNodes = obj.iVertexNode0 + iVertices;
        end
        
        % ---- BOUNDARY AND INTERIOR NODES
        
        % Get an arbitrarily-ordered list of the nodes on the outer boundary
        function iNodes = getBoundaryNodes(obj)
            
            iEdgeIndices = obj.hMesh.getBoundaryEdges();
            iEdgeVertexIndices = obj.hMesh.edgeVertices(iEdgeIndices,:);
            
            iNodes = zeros(numel(iEdgeIndices), obj.numInteriorNodesPerEdge);
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
        
        % ---- NODAL TRIANGULATED MESH
        % Make a denser mesh by triangulating nodes.
        
        function globalFaces = getNodalMesh(obj)
            
            numFaces = obj.hMesh.getNumFaces();
            
            % Create a subtriangulation of a single triangle.
            % We'll reuse this over the whole mesh.
            
            numTrisPerFace = (obj.N-1)^2;
            localFaces = zeros(numTrisPerFace,3);
            
            % Make a helper grid of node indices.
            % We'll just use a portion of it.
            nodeGrid = bsxfun(@plus, 1:obj.N, cumsum([0, obj.N:-1:2])');
            
            ff = 1;
            for jj = 1:(obj.N-1)
                iMax = obj.N-jj;
                for ii = 1:iMax
                    idx1 = nodeGrid(jj,ii);
                    idx2 = nodeGrid(jj,ii+1);
                    idx3 = nodeGrid(jj+1,ii);
                    idx4 = nodeGrid(jj+1,ii+1);
                    
                    localFaces(ff,:) = [idx1, idx2, idx3];
                    ff = ff+1;
                    
                    if ii < iMax
                        localFaces(ff,:) = [idx3, idx2, idx4];
                        ff = ff+1;
                    end
                end
            end
            
            globalFaces = zeros(numTrisPerFace*numFaces, 3);
            gg = 1;
            
            for ff = 1:numFaces
                
                iFaceNodes = obj.getFaceNodes(ff);
                
                subtriangleNodes = iFaceNodes(localFaces);
                globalFaces(gg:gg+numTrisPerFace-1,:) = subtriangleNodes;
                gg = gg + numTrisPerFace;
                
            end
            
        end
        
    end
    
    
end