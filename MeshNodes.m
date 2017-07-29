classdef MeshNodes < handle

    properties
        N;  % nodes per edge of element
        
        % Mesh topology and geometry.
        faces;
        vertices;
        faceEdges;
        faceEdgeOrientations;
        edgeVertices;
        
        % Nodes are numbered consecutively in three groups: vertex nodes,
        % edge nodes and face nodes.
        %   iVertexNode0 is the index beginning the vertex nodes
        %   iEdgeNode0 is the index beginning the edge nodes
        %   iFaceNode0 is the index beginning the face nodes
        iVertexNode0;
        iEdgeNode0;
        iFaceNode0;
        
        nodesPerEdge;  % number of interior nodes per edge (vertices excluded)
        nodesPerFace;  % number of interior nodes per face (edges and corners excluded)
        
        nodesLocal_rs; % cached coordinates of basis nodes
        nodesInteriorLocal_rs; % cached coordinates of interior basis nodes (no edges or vertices)
        
        nodesLocalEdge; % cached positions of nodes on basis interval
        nodesLocalEdgeInterior; % interior nodes on basis interval
    end

    methods
        
        % ---- CONSTRUCTOR
        
        function obj = MeshNodes(faces, vertices, N)
            obj.N = N;
            obj.faces = faces;
            obj.vertices = vertices;
            [obj.faceEdges, obj.faceEdgeOrientations, obj.edgeVertices] = VVMesh.getEdges(obj.faces);
            
            obj.nodesPerEdge = max(0, N-2);
            obj.nodesPerFace = max(0, (N-3)*(N-2)/2);
            
            obj.iVertexNode0 = 0;
            obj.iEdgeNode0 = obj.getNumVertices();
            obj.iFaceNode0 = obj.iEdgeNode0 + obj.nodesPerEdge*obj.getNumEdges();
            
            obj.nodesLocalEdge = support.gaussLobatto(N);
            obj.nodesLocalEdgeInterior = obj.nodesLocalEdge(2:end-1);
            
            [rr, ss] = support2d.nodes2d(N);
            obj.nodesLocal_rs = [rr,ss];
            [iVertices,iEdges,iCenters] = support2d.classifyNodes(N);
            obj.nodesInteriorLocal_rs = obj.nodesLocal_rs(iCenters,:);
            
            
        end
        
        % ---- NUMBERS OF THINGS

        function n = getNumVertices(obj)
            n = size(obj.vertices, 1);
        end
        
        function n = getNumEdges(obj)
            n = size(obj.edgeVertices, 1);
        end

        function n = getNumFaces(obj)
            n = size(obj.faces, 1);
        end

        function n = getNumNodes(obj)
            n = obj.getNumVertices() + obj.nodesPerEdge*obj.getNumEdges() + obj.nodesPerFace*obj.getNumFaces();
        end
        
        % ---- TOPOLOGY
        
        function [edgeIndices] = getBoundaryEdges(obj)
            edgeUses = obj.faceEdges(:);
            edgeUseCounts = hist(edgeUses, unique(edgeUses));
            
            edgeIndices = find(edgeUseCounts == 1);
        end
        
        % Get an arbitrarily-ordered list of the nodes on the outer boundary
        function iNodes = getBoundaryNodes(obj)
            
            iEdgeIndices = obj.getBoundaryEdges();
            iEdgeVertexIndices = obj.edgeVertices(iEdgeIndices,:);

            iNodes = zeros(numel(iEdgeIndices), obj.nodesPerEdge);
            for ii = 1:length(iEdgeIndices)
                iNodes(ii,:) = obj.getEdgeNodeIndices(iEdgeIndices(ii));
            end
            iNodes = unique(iNodes(:));

            % throw in the edge vertices
            iEdgeCornerNodes = unique(obj.getCornerNodeIndices(iEdgeVertexIndices(:)));

            iNodes = [iNodes; iEdgeCornerNodes];
            
        end
        
        function iCenterNodes = getInteriorNodes(obj)
            
            iCenterNodes = setdiff(1:obj.getNumNodes(), obj.getBoundaryNodes());
            
        end
        
        
        % ---- LOCAL TO GLOBAL MAP
        
        function iNodeGlobal = local2global_vertex(obj, iFaceGlobal, iCornerLocal)
            % Return the global index of a node on a vertex of a face.
            iNodeGlobal = obj.getCornerNodeIndices(obj.faces(iFaceGlobal, iCornerLocal));
        end

        function iNodesGlobal = local2global_edgeInterior(obj, iFaceGlobal, iEdgeLocal)
            % Return global indices of nodes interior to one edge of a face, in
            % counterclockwise order with respect to the face.
            iEdgeGlobal = obj.faceEdges(iFaceGlobal, iEdgeLocal);
            iNodesGlobal = obj.getEdgeNodeIndices(iEdgeGlobal);
            
            if obj.faceEdgeOrientations(iFaceGlobal, iEdgeLocal) < 0
                iNodesGlobal = iNodesGlobal(end:-1:1);
            end
        end
        
        function iNodesGlobal = local2global_edge(obj, iFaceGlobal, iEdgeLocal)
            % Return global indices of all edges on one edge of a face,
            % including endpoints (vertices), in counterclockwise order
            % with respect to the face.
            
            iEdgeGlobal = obj.faceEdges(iFaceGlobal, iEdgeLocal);
            iEdgeNodesGlobal = obj.getEdgeNodeIndices(iEdgeGlobal);
            
            iVertexNodesGlobal = obj.getCornerNodeIndices(obj.edgeVertices(obj.faceEdges(iFaceGlobal, iEdgeLocal),:));
            
            if obj.faceEdgeOrientations(iFaceGlobal, iEdgeLocal) < 0
                iNodesGlobal = [iVertexNodesGlobal(2), iEdgeNodesGlobal(end:-1:1), iVertexNodesGlobal(1)];
            else
                iNodesGlobal = [iVertexNodesGlobal(1), iEdgeNodesGlobal, iVertexNodesGlobal(2)];
            end
        end

        function iNodesGlobal = local2global_faceInterior(obj, iFaceGlobal)
            % Return global indices of interior nodes of a face.
            iNodesGlobal = obj.getFaceNodeIndices(iFaceGlobal);
        end

        function iGlobal = local2global_face(obj, iFace)
            % Return global indices of all nodes associated with a face,
            % oriented as on the basis element.
            [iCorners, iEdges, iCenter] = support2d.classifyNodes(obj.N);
            
            iGlobal = zeros(obj.N*(obj.N+1)/2, 1);
            
            iGlobal(iCorners) = obj.local2global_vertex(iFace, 1:3);
            iGlobal(iEdges{1}) = obj.local2global_edgeInterior(iFace, 1);
            iGlobal(iEdges{2}) = obj.local2global_edgeInterior(iFace, 2);
            iGlobal(iEdges{3}) = obj.local2global_edgeInterior(iFace, 3);
            iGlobal(iCenter) = obj.local2global_faceInterior(iFace);
        end
        
        
        % ---- NODE INDEX CALCULATORS
        % getCornerNodeIndices, getEdgeNodeIndices and getFaceNodeIndices return
        % disjoint node lists... these are for storage of nodes and matrix
        % indexing.
        
        % Get indices of nodes located at vertices of the mesh
        function iNodes = getCornerNodeIndices(obj, iVertices)
            iNodes = obj.iVertexNode0 + iVertices;
        end

        % Get indices of nodes in interior of edge
        function iNodes = getEdgeNodeIndices(obj, iEdge)
            iNodes = obj.iEdgeNode0 + (iEdge-1)*obj.nodesPerEdge + ...
                (1:obj.nodesPerEdge);
        end
        
        % Get indices of nodes in interior of faces
        function iNodes = getFaceNodeIndices(obj, iFace)
            iNodes = obj.iFaceNode0 + (iFace-1)*obj.nodesPerFace + ...
                (1:obj.nodesPerFace);
        end
        
        
        % ---- NODE COORDINATE CALCULATORS

        function xyz = getVertexNodeCoordinates(obj, iVertex)
            xyz = obj.vertices(iVertex,:);
        end

        function xyz = getEdgeNodeCoordinates(obj, iEdge)
            if obj.N < 3
                xyz = zeros(0,3);
                return
            end
            
            v1 = obj.vertices(obj.edgeVertices(iEdge,1),:);
            v2 = obj.vertices(obj.edgeVertices(iEdge,2),:);
            
            %d = linspace(0, 1, obj.N)';
            %d = d(2:end-1);
            
            d = 0.5 + 0.5*obj.nodesLocalEdgeInterior';
            
            x = v1(1) + (v2(1)-v1(1))*d;
            y = v1(2) + (v2(2)-v1(2))*d;
            xyz = [x, y];
        end
        
        function xyz = getFaceNodeCoordinates(obj, iFace)
            if obj.N < 4
                xyz = zeros(0,3);
                return
            end
            threeVertices = obj.vertices(obj.faces(iFace,:),:);
            xyz = support2d.rs2xy(threeVertices', obj.nodesInteriorLocal_rs')';
        end
        
        function threeVertices = getFaceVertices(obj, iFace)
            threeVertices = obj.vertices(obj.faces(iFace,:),:);
        end
        
        
        function xyz = getNodeCoordinates(obj)
            numNodes = obj.getNumNodes();
            xyz = zeros(numNodes,2);

            % Nodes, section 1/3: Vertices
            xyz(1:obj.getNumVertices(),:) = obj.vertices;
            
            % Nodes, section 2/3: Edge-centers
            for iEdge = 1:obj.getNumEdges()
                xyz(obj.getEdgeNodeIndices(iEdge),:) = obj.getEdgeNodeCoordinates(iEdge);
            end

            % Nodes, section 3/3: Face-centers
            for iFace = 1:obj.getNumFaces()
                xyz(obj.getFaceNodeIndices(iFace),:) = obj.getFaceNodeCoordinates(iFace);
            end
        end
            
            
        
        % ---- JACOBIAN CALCULATOR
        
        
        function jac = getLinearJacobian(obj, iFace)
            % This is the Jacobian of the mapping from (r,s) to (x,y).
            % 
            % [T, v0] = support2d.rs2xy_affineParameters(xyTri);
            % xy = bsxfun(@plus, v0, T*rs);
            % so the Jacobian is just T.
            
            threeVertices = obj.vertices(obj.faces(iFace,:),:);
            jac = support2d.rs2xy_affineParameters(threeVertices');
            
        end
        
        
        function jac1d = getLinearJacobian1d(obj, iFace, iEdge)
            % This is the Jacobian of the mapping from r to (x,y).
            %
            % xy = (v1+v2)/2 + (v2-v1)/2 * r;
            %
            % so d(xy)/dr = (v2-v1)/2.
            
            i0 = iEdge;
            i1 = 1 + mod(iEdge,3);
            twoVertices = obj.vertices(obj.faces(iFace,[i0,i1]), :);
            
            jac1d = 0.5*(twoVertices(2,:) - twoVertices(1,:))';
            
            
        end
        
        % ---- SENSITIVITY ANALYSIS
        
        
        function dJdv = getLinearJacobianSensitivity(obj, iFace)
            % dJdv = getLinearJacobianSensitivity(obj, iFace)
            %
            % dJdv is a 4D array indexed by (row, col, vertex, xy).
            
            dJdv = zeros(2,2,3,2); % row, col, vertex, vertex xy
            
            threeVertices = obj.vertices(obj.faces(iFace,:),:);
            
            for iVert = 1:3
                for iXY = 1:2
                    dThreeV = 0*threeVertices;
                    dThreeV(iVert,iXY) = 1.0;
                    
                    dJdv(:,:,iVert,iXY) = support2d.rs2xy_affineParameterSensitivities(threeVertices', dThreeV');
                end
            end
        end
        
        
        % ---- OPERATORS
        
        function [outDx, outDy, count] = getGradientOperators(obj)
            
            numNodes = obj.getNumNodes();
            outDx = sparse(numNodes, numNodes);
            outDy = sparse(numNodes, numNodes);
            
            count = zeros(numNodes, 1);
            
            numFaces = obj.getNumFaces();
            
            [Dr, Ds] = support2d.gradients(obj.N, obj.nodesLocal_rs(:,1), obj.nodesLocal_rs(:,2));
            
            for ff = 1:numFaces
            %for ff = fff
                jac = obj.getLinearJacobian(ff);
                invJac = inv(jac);
                
                Dx = Dr*invJac(1,1) + Ds*invJac(2,1);
                Dy = Dr*invJac(1,2) + Ds*invJac(2,2);
                
                
                % first silly approach: on edges between faces we will
                % repeatedly overwrite the matrix elements.  Bogus!!
                iGlobal = obj.local2global(ff);
                outDx(iGlobal, iGlobal) = outDx(iGlobal, iGlobal) + Dx;
                outDy(iGlobal, iGlobal) = outDy(iGlobal, iGlobal) + Dy;
                
                count(iGlobal) = count(iGlobal) + 1;
            end
            
            % this handles averaging on boundaries
            normalizer = spdiags(1./count, 0, numNodes, numNodes);
            outDx = normalizer * outDx;
            outDy = normalizer * outDy;
            
            
        end
        
        function outI = getInterpolationOperator(obj, xs, ys)
            
            numPts = length(xs);
            numNodes = obj.getNumNodes();
            outI = sparse(numPts, numNodes);
            count = zeros(numPts, 1);
            
            nodalBasis = NodalBasis(obj.N);
            
            tr = triangulation(obj.faces, obj.vertices);
            iFaces = tr.pointLocation(xs(:), ys(:));
            numFaces = obj.getNumFaces();
            
            for ff = 1:numFaces
                ii = find(iFaces == ff);
                xy = [xs(ii)'; ys(ii)'];
                
                xyTri = obj.getFaceVertices(ff)';
                rs = support2d.xy2rs(xyTri, xy);
                
                M = nodalBasis.interpolationMatrix(rs(1,:), rs(2,:));
                
                iGlobal = obj.local2global(ff);
                outI(ii, iGlobal) = outI(ii, iGlobal) + M;
                count(ii) = count(ii)+1;
                
            end

            % this handles averaging on boundaries
            normalizer = spdiags(1./count, 0, numPts, numPts);
            outI = normalizer * outI;
            
        end
        
        
        % ---- VISUALIZATION
        %
        
        function plotMatrix(obj, A, varargin)
            % Draw every FEM node affected by this matrix.
            
            assert(size(A,1) == obj.getNumNodes());
            
            [ii,~,~] = find(A);
            
            xy = obj.getNodeCoordinates();
            
            iii = unique(ii);
            plot(xy(iii,1), xy(iii,2), varargin{:});
        end
        
        function plotNodeIndices(obj)
            
            xy = obj.getNodeCoordinates();
            
            for ii = 1:size(xy,1)
                text(xy(ii,1), xy(ii,2), num2str(ii));
            end
            
        end
        
        
    end % methods


end