classdef MeshNodes < handle

    properties
        N;  % nodes per edge of element
        faces;
        vertices;
        faceEdges;
        faceEdgeOrientations;
        edgeVertices;
        
        Dvertices; % Jacobian of vertices
        
        iVertexNode0;
        iEdgeNode0;
        iFaceNode0;
        
        nodesPerEdge;
        nodesPerFace;
        
        nodesLocal_rs;
        nodesInteriorLocal_rs;
        
        nodesLocalEdge;
        nodesLocalEdgeInterior;
    end

    methods
        
        % ---- CONSTRUCTOR
        
        function obj = MeshNodes(faces, vertices, N)
            obj.N = N;
            obj.faces = faces;
            obj.vertices = vertices;
            obj.Dvertices = 0*vertices;
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
            [~,~,iCenters] = support2d.classifyNodes(N);
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
        
        function edgeIndices = getBoundaryEdges(obj)
            edgeUses = obj.faceEdges(:);
            edgeUseCounts = hist(edgeUses, unique(edgeUses));
            
            edgeIndices = find(edgeUseCounts == 1);
        end
        
        function iNodes = getBoundaryNodes(obj)
            
            iEdgeIndices = obj.getBoundaryEdges();
            iEdgeVertexIndices = obj.edgeVertices(iEdgeIndices,:);

            iNodes = zeros(numel(iEdgeIndices), obj.nodesPerEdge);
            for ii = 1:length(iEdgeIndices)
                iNodes(ii,:) = obj.getEdgeNodes(iEdgeIndices(ii));
            end
            iNodes = unique(iNodes(:));

            % throw in the edge vertices
            iEdgeCornerNodes = unique(obj.getCornerNodes(iEdgeVertexIndices(:)));

            iNodes = [iNodes; iEdgeCornerNodes];
            
        end
        
        function iCenterNodes = getInteriorNodes(obj)
            
            iCenterNodes = setdiff(1:obj.getNumNodes(), obj.getBoundaryNodes());
            
        end
        
        
        % ---- LOCAL TO GLOBAL MAP
        
        function iNodeGlobal = local2global_corner(obj, iFaceGlobal, iCornerLocal)
            iNodeGlobal = obj.getCornerNodes(obj.faces(iFaceGlobal, iCornerLocal));
        end

        function iNodesGlobal = local2global_edge(obj, iFaceGlobal, iEdgeLocal)
            iEdgeGlobal = obj.faceEdges(iFaceGlobal, iEdgeLocal);
            iNodesGlobal = obj.getEdgeNodes(iEdgeGlobal);
            
            if obj.faceEdgeOrientations(iFaceGlobal, iEdgeLocal) < 0
                iNodesGlobal = iNodesGlobal(end:-1:1);
            end
        end

        function iNodesGlobal = local2global_face(obj, iFaceGlobal)
            iNodesGlobal = obj.getFaceNodes(iFaceGlobal);
        end

        function iGlobal = local2global(obj, iFace)
            
            [iCorners, iEdges, iCenter] = support2d.classifyNodes(obj.N);
            
            iGlobal = zeros(obj.N*(obj.N+1)/2, 1);
            
            iGlobal(iCorners) = obj.local2global_corner(iFace, 1:3);
            iGlobal(iEdges{1}) = obj.local2global_edge(iFace, 1);
            iGlobal(iEdges{2}) = obj.local2global_edge(iFace, 2);
            iGlobal(iEdges{3}) = obj.local2global_edge(iFace, 3);
            iGlobal(iCenter) = obj.local2global_face(iFace);
        end
        
        
        % ---- NODE INDEX CALCULATORS
        
        function iNodes = getCornerNodes(obj, iVertices)
            iNodes = obj.iVertexNode0 + iVertices;
        end

        function iNodes = getEdgeNodes(obj, iEdge)
            iNodes = obj.iEdgeNode0 + (iEdge-1)*obj.nodesPerEdge + ...
                (1:obj.nodesPerEdge);
        end

        function iNodes = getFaceNodes(obj, iFace)
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
                xyz(obj.getEdgeNodes(iEdge),:) = obj.getEdgeNodeCoordinates(iEdge);
            end

            % Nodes, section 3/3: Face-centers
            for iFace = 1:obj.getNumFaces()
                xyz(obj.getFaceNodes(iFace),:) = obj.getFaceNodeCoordinates(iFace);
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
        
        % ---- SENSITIVITY ANALYSIS
        
        function setTestPerturbation(obj, iVert, iXY)
            % Set the vertex Jacobian to represent a unit perturbation
            % of one vertex in one direction.
            assert(iVert >= 1 && iVert <= obj.getNumVertices());
            assert(iXY >= 1 && iXY <= 2, 'iXY must be 1 or 2');
            
            obj.Dvertices = 0*obj.vertices;
            obj.Dvertices(iVert, iXY) = 1.0;
        end
        
%         function Djac = getLinearJacobianSensitivity(obj, iFace)
%             % This is the sensitivity of the Jacobian of the mapping from
%             % (r,s) to (x,y).
%             
%             threeVertices = obj.vertices(obj.faces(iFace,:),:);
%             dThreeV = obj.Dvertices(obj.faces(iFace,:),:);
%             Djac = support2d.rs2xy_affineParameterSensitivities(threeVertices', dThreeV');
%         end
        
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
        
        function dJdv = getAllLinearJacobianSensitivities(obj)
            
            numFaces = obj.getNumFaces();
            
            % Initialize the Jacobian Jacobian.  :-D
            % It should be sparse.
            dJdv = cell(2,2);
            for nn = 1:4
                dJdv{nn} = sparse(size(obj.vertices,1), size(obj.vertices,2));
            end
            
            for ff = 1:numFaces
                
                dJ = obj.getLinearJacobianSensitivity(ff);
                
                for ii = 1:2
                    for jj = 1:2
                        dJ_sparse = sparse(size(obj.vertices,1), size(obj.vertices,2));
                        dJ_sparse(obj.faces(ff,:), :) = dJ(ii,jj,:,:);
                
                        disp('----')
                        disp(obj.faces(ff,:))
                        disp(dJ_sparse)
                        dJdv{ii,jj} = dJdv{ii,jj} + dJ_sparse;
                    end
                end
                
            end
            
        end % getAllLinearJacobianSensitivities
        
%         function Djac = getLinearJacobianSensitivity(obj, iFace, iVertInFace, iDirection)
%             % This is the sensitivity of the Jacobian of the mapping from
%             % (r,s) to (x,y).
%             
%             threeVertices = obj.vertices(obj.faces(iFace,:),:);
%             Dvertices = 0*threeVertices;
%             Dvertices(iVertInFace, iDirection) = 1;
%             Djac = support2d.rs2xy_affineParameterSensitivities(threeVertices', Dvertices');
%         end
        
        
        
    end


end