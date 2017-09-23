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
            
            iVertexNodes = obj.getVertexNodes(obj.hMesh.faceVertices(iFace,:));
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
            iNodes = obj.iFaceNode0 + (iFace-1)*obj.numInteriorNodesPerFace + ...
                (1:obj.numInteriorNodesPerFace);
        end
        
        function iNodes = getFaceEdgeNodes(obj, iFace, iEdgeLocal)
            iEdgeGlobal = obj.hMesh.faceEdges(iFace, iEdgeLocal);
            orientation = obj.hMesh.faceEdgeOrientations(iFace, iEdgeLocal);
            iInteriorNodes = obj.getEdgeInteriorNodes(iEdgeGlobal, orientation);
            
            edgeVertsLocal = [iEdgeLocal, 1 + mod(iEdgeLocal,3)];
            iVertexNodes = obj.getVertexNodes(obj.hMesh.faces(iFace, edgeVertsLocal));
            
            iNodes = [iVertexNodes(1), iInteriorNodes, iVertexNodes(2)];
        end
        
        function iNodes = getFaceEdgeInteriorNodes(obj, iFace, iEdgeLocal)
            iEdgeGlobal = obj.hMesh.faceEdges(iFace, iEdgeLocal);
            orientation = obj.hMesh.faceEdgeOrientations(iFace, iEdgeLocal);
            iNodes = obj.getEdgeInteriorNodes(iEdgeGlobal, orientation);
        end
        
        function iNodes = getFaceVertexNodes(obj, iFace, iVerticesLocal)
            iVertices = obj.hMesh.faces(iFace, iVerticesLocal);
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
                    assert(numel(varargin{1}) == numel(iEdge), 'Must supply one orientation per edge');
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
        
        
    end
    
    
end