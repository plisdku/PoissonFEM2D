classdef MeshTopology < handle
% MeshTopology Topology of faces, edges, vertices and FEM nodes in 2D
% triangulated mesh.
%
% The mesh geometry is not handled here.
    
    properties % (Access = private)
        % Mesh topology.
        faceVertices;
        faceEdges;
        faceEdgeOrientations; % +1 or -1
        edgeVertices;  % not oriented.
        
        numVertices; % inferred from the mesh
    end

    methods
        
        % ---- CONSTRUCTOR
        
        function obj = MeshTopology(faceVertices)
            obj.faceVertices = faceVertices;
            [obj.faceEdges, obj.faceEdgeOrientations, obj.edgeVertices] = VVMesh.getEdges(obj.faceVertices);
            obj.numVertices = numel(unique(faceVertices(:)));
        end
        
        % ---- NUMBERS OF THINGS
        
        function n = getNumVertices(obj)
            n = obj.numVertices;
        end
        
        function n = getNumEdges(obj)
            n = size(obj.edgeVertices, 1);
        end

        function n = getNumFaces(obj)
            n = size(obj.faceVertices, 1);
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
        
        function iVertices = getFaceEdgeVertices(obj, iFace)
            % Return ordered vertex indices for an edge of a face
            %
            % getFaceEdgeVertices(iFace, iEdge)
            %
            % Convenience function.
            
            [iEdges, orientations] = obj.getFaceEdges(iFace);
            iVertices = obj.getEdgeVertices(iEdges, orientations);
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
        
    end % methods


end