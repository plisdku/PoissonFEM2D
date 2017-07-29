classdef TriNodalMesh < handle
% TriNodalMesh Geometry and topology of a triangulated FEM mesh with nodes
    
    properties
        N;
        topology;   % MeshTopology
        basis;      % BasisNodes
        basis1d;
        
        vertices;
    end
    
    
    
    methods
        
        % ---- CONSTRUCTOR
        
        
        function obj = TriNodalMesh(N, faces, vertices)
            
            obj.N = N;
            obj.topology = MeshTopology(faces, N);
            obj.basis = BasisNodes(N);
            obj.basis1d = BasisNodes1d(N);
            
            obj.vertices = vertices;
            
        end
        
        % ---- JACOBIANS
        
        function jac = getLinearJacobian(obj, iFace)
            % Calculate the Jacobian of the mapping from (r,s) to (x,y).
            % 
            % [T, v0] = support2d.rs2xy_affineParameters(xyTri);
            % xy = bsxfun(@plus, v0, T*rs);
            % so the Jacobian is just T.
            
            threeVertices = obj.vertices(obj.topology.getFaceVertices(iFace),:);
            jac = support2d.rs2xy_affineParameters(threeVertices');
        end
        
        function jac1d = getLinearJacobian1d(obj, iEdge, orientation)
            % This is the Jacobian of the mapping from r to (x,y).
            %
            % xy = (v1+v2)/2 + (v2-v1)/2 * r;
            %
            % so d(xy)/dr = (v2-v1)/2.
            
            twoVertices = obj.vertices(obj.topology.getEdgeVertices(iEdge, orientation), :);
            jac1d = 0.5*(twoVertices(2,:) - twoVertices(1,:))';
        end
%         
%         function jac1d = getLinearJacobian1d(obj, iFace, iLocalEdge)
%             % This is the Jacobian of the mapping from r to (x,y).
%             %
%             % xy = (v1+v2)/2 + (v2-v1)/2 * r;
%             %
%             % so d(xy)/dr = (v2-v1)/2.
%             
%             twoVertices = obj.vertices(obj.topology.getFaceEdgeVertices(iFace, iLocalEdge), :);
%             jac1d = 0.5*(twoVertices(2,:) - twoVertices(1,:))';
%         end
        
        
        function dJdv = getLinearJacobianSensitivity(obj, iFace)
            % dJdv = getLinearJacobianSensitivity(iFace)
            %
            % dJdv is a 4D array indexed by (row, col, vertex, xy).
            %
            % Affine transformations are simple and actually the face
            % argument is unnecessary.
            
            dJdv = zeros(2,2,3,2); % row, col, vertex, vertex xy
            
            threeVertices = obj.vertices(obj.topology.getFaceVertices(iFace),:);
            
            % I could build this without the affineParameter... abstraction
            % if I did it like getLinearJacobianSensitivity1d.  Think about
            % it.  (Long term, what about non-affine transformations?)
            for iVert = 1:3
                for iXY = 1:2
                    dThreeV = 0*threeVertices;
                    dThreeV(iVert,iXY) = 1.0;
                    
                    dJdv(:,:,iVert,iXY) = support2d.rs2xy_affineParameterSensitivities(threeVertices', dThreeV');
                end
            end
        end
        
        
        function dJdv1d = getLinearJacobianSensitivity1d(obj, iFace, iLocalEdge)
            % dJdv1d = getLinearJacobianSensitivity1d(iFace, iLocalEdge)
            %
            % dJdv is a 3D array indexed by (row, vertex, xy).
            %
            % Affine transformations are simple and actually the face and
            % edge arguments are unnecessary.
            
            % Transformation maps (r) -> (x, y)
            % Jacobian is thus [dx/dr; dy/dr]
            % Jacobian dimensions are [2,1].  We'll represent it as a
            % column vector.
            % 
            dJdv1d = zeros(2,2,2); % row, vertex, vertex xy
            
            %twoVertices = obj.vertices(obj.getFaceEdgeVertices(iFace, iLocalEdge),:);
            
            %jac1d = 0.5*(v2 - v1)    (a column vector)
            %
            % d/dv1x = 0.5*( [-1; 0] )
            % d/dv1y = 0.5*( [0; -1] )
            % d/dv2x = 0.5*( [1; 0] )
            % d/dv2y = 0.5*( [0; 1] )
            
            dJdv1d(:,1,1) = 0.5*[-1; 0];
            dJdv1d(:,1,2) = 0.5*[0; -1];
            dJdv1d(:,2,1) = 0.5*[1; 0];
            dJdv1d(:,2,2) = 0.5*[0; 1];
        end
        
        
        
        % ---- OPERATORS
        
        function [outDx, outDy, count] = getGradientOperators(obj)
            
            numNodes = obj.topology.getNumNodes();
            outDx = sparse(numNodes, numNodes);
            outDy = sparse(numNodes, numNodes);
            
            count = zeros(numNodes, 1);
            
            numFaces = obj.topology.getNumFaces();
            
            rs = obj.basis.getNodes();
            [Dr, Ds] = support2d.gradients(obj.N, rs(:,1), rs(:,2));
            
            for ff = 1:numFaces
            %for ff = fff
                jac = obj.getLinearJacobian(ff);
                invJac = inv(jac);
                
                Dx = Dr*invJac(1,1) + Ds*invJac(2,1);
                Dy = Dr*invJac(1,2) + Ds*invJac(2,2);
                
                
                % first silly approach: on edges between faces we will
                % repeatedly overwrite the matrix elements.  Bogus!!
                iGlobal = obj.topology.getFaceNodes(ff);
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
            numNodes = obj.topology.getNumNodes();
            outI = sparse(numPts, numNodes);
            count = zeros(numPts, 1);
            
            if numPts == 0
                return
            end
            
            tr = triangulation(obj.topology.getFaceVertices(), obj.vertices);
            iFaces = tr.pointLocation(xs(:), ys(:));
            numFaces = obj.topology.getNumFaces();
            
            for ff = 1:numFaces
                ii = find(iFaces == ff);
                xy = [xs(ii)'; ys(ii)'];
                
                xyTri = obj.vertices(obj.topology.getFaceVertices(ff), :)';
                rs = support2d.xy2rs(xyTri, xy);
                
                M = obj.basis.interpolationMatrix(rs(1,:), rs(2,:));
                
                iGlobal = obj.topology.getFaceNodes(ff);
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
            
            assert(size(A,1) == obj.topology.getNumNodes());
            
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
        
        % ---- NODE COORDINATES
        
        function xyz = getVertexNodeCoordinates(obj, iVertex)
            % Get [x,y] coordinates of the node on a given vertex
            %
            % getVertexNodeCoordinates(iVertex)
            
            xyz = obj.vertices(iVertex,:);
        end

        function xyz = getEdgeInteriorNodeCoordinates(obj, iEdge, varargin)
            % Get ordered [x,y] coordinates of interior nodes on a given edge
            %
            % getEdgeInteriorNodeCoordinates(iEdge)
            % getEdgeInteriorNodeCoordinates(iEdge, orientation)
            
            if obj.N < 3
                xyz = zeros(0,2);
                return
            end
            
            verts = obj.vertices(obj.topology.getEdgeVertices(iEdge),:);
            %v2 = obj.vertices(obj.topology.getEdgeVertices(iEdge),:);
            
            %d = linspace(0, 1, obj.N)';
            %d = d(2:end-1);
            
            d = 0.5 + 0.5*transpose(obj.basis1d.getInteriorNodes());
            
            x = verts(1,1) + (verts(2,1)-verts(1,1))*d;
            y = verts(1,2) + (verts(2,2)-verts(1,2))*d;
            xyz = [x, y];
            
            if ~isempty(varargin)
                orientation = varargin{1};
                if orientation < 0
                    xyz = xyz(end:-1:1,:);
                end
            end
        end
        
        function xy = getEdgeNodeCoordinates(obj, iEdge, varargin)
            % Get ordered [x,y] coordinates of nodes on a given edge
            %
            % getEdgeNodeCoordinates(iEdge)
            % getEdgeNodeCoordinates(iEdge, orientation)
            
            verts = obj.vertices(obj.topology.getEdgeVertices(iEdge),:);
            
            d = 0.5 + 0.5*transpose(obj.basis1d.getNodes());
            
            x = verts(1,1) + (verts(2,1)-verts(1,1))*d;
            y = verts(1,2) + (verts(2,2)-verts(1,2))*d;
            
            xy = [x, y];
            
            if ~isempty(varargin)
                orientation = varargin{1};
                if orientation < 0
                    xy = xy(end:-1:1,:);
                end
            end
        end
        
        
        function xy = getFaceInteriorNodeCoordinates(obj, iFace)
            % Get ordered [x,y] coordinates of interior nodes in a given face
            %
            % getFaceInteriorNodeCoordinates(iFace)
            
            if obj.N < 4
                xy = zeros(0,2);
                return
            end
            threeVertices = obj.vertices(obj.topology.getFaceVertices(iFace),:);
            xy = support2d.rs2xy(threeVertices', obj.basis.getInteriorNodes()')';
        end
        
        
        function xy = getFaceNodeCoordinates(obj, iFace)
            % Get ordered [x,y] coordinates of all nodes in a given face
            %
            % getFaceNodeCoordinates(iFace)
            
            if obj.N < 4
                xy = zeros(0,2);
                return
            end
            threeVertices = obj.vertices(obj.topology.getFaceVertices(iFace),:);
            xy = support2d.rs2xy(threeVertices', obj.basis.getNodes()')';
        end
        
        
        function xyz = getNodeCoordinates(obj)
            % Get ordered [x,y] coordinates of all nodes in the mesh
            %
            % getNodeCoordinates()
            %
            % The global node ordering is defined in MeshTopology:
            % nodes = [vertex nodes;
            %          edge interior nodes;
            %          face interior nodes]
            % with vertex nodes ordered by vertex number, edge nodes
            % ordered by edge number and face nodes ordered by face number.
            
            numNodes = obj.topology.getNumNodes();
            xyz = zeros(numNodes,2);

            % Nodes, section 1/3: Vertices
            xyz(1:obj.topology.getNumVertices(),:) = obj.vertices;
            
            % Nodes, section 2/3: Edge-centers
            for iEdge = 1:obj.topology.getNumEdges()
                xyz(obj.topology.getEdgeInteriorNodes(iEdge),:) = obj.getEdgeInteriorNodeCoordinates(iEdge);
            end

            % Nodes, section 3/3: Face-centers
            for iFace = 1:obj.topology.getNumFaces()
                xyz(obj.topology.getFaceNodes(iFace),:) = obj.getFaceNodeCoordinates(iFace);
            end
        end
        
        function xy = getBoundaryNodeCoordinates(obj)
            xy = obj.getNodeCoordinates();
            xy = xy(obj.topology.getBoundaryNodes(),:);
        end
        
        function xy = getInteriorNodeCoordinates(obj)
            xy = obj.getNodeCoordinates();
            xy = xy(obj.topology.getInteriorNodes(),:);
        end
        
    end % methods
    
    
end