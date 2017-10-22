classdef TriNodalMesh < handle
% TriNodalMesh Geometry and topology of a triangulated FEM mesh with nodes
    
    properties
        xyNodes; % the vertex positions come first followed by other nodes
        
        hMesh@MeshTopology;
        hFieldNodes@NodalTopology;
        hGeomNodes@NodalTopology;
        hQuadNodes@NodalTopology;
        
    end
    
    
    
    methods
        
        % ---- CONSTRUCTOR
        
        
        function obj = TriNodalMesh(faces, xyNodes, N_field, N_geom, N_quad)
            obj.hMesh = MeshTopology(faces);
            obj.hFieldNodes = NodalTopology(obj.hMesh, N_field);
            obj.hGeomNodes = NodalTopology(obj.hMesh, N_geom);
            obj.hQuadNodes = NodalTopology(obj.hMesh, N_quad);
            
            assert(size(xyNodes,2) == 2, 'Vertices must be Nx2'); % test because 3d verts are common
            obj.xyNodes = xyNodes;
        end
        
        % ---- NODAL COORDINATES
        % a.k.a. coordinate transformation
        
        function xyz = getVertexNodeCoordinates(obj, iVertex)
            % Get [x,y] coordinates of the node on a given vertex
            %
            % getVertexNodeCoordinates(iVertex)
            
            xyz = obj.xyNodes(iVertex,:);
        end

        function xy = getEdgeCoordinates(obj, iEdge, rs, varargin)
            % Get ordered [x,y] coordinates at positions on a given edge.
            % getEdgeNodeCoordinates(obj, iEdge, rs)
            % getEdgeNodeCoordinates(obj, iEdge, rs, orientation)
            
            if nargin < 4
                orientation = 1;
            else
                orientation = varargin{1};
            end
            
            M = obj.hGeomNodes.basis1d.interpolationMatrix_r(rs);
            xy = M*obj.xyNodes(obj.hGeomNodes.getEdgeNodes(iEdge, orientation),:);
            
        end
        
        function xy = getFaceCoordinates(obj, iFace, rr, ss)
            % Get ordered [x,y] coordinates at positions on a given face.
            % getFaceCoordinates(obj, iFace, rr, ss)
            
            M = obj.hGeomNodes.basis.interpolationMatrix_rs(rr, ss);
            xy = M*obj.xyNodes(obj.hGeomNodes.getFaceNodes(iFace),:);
            
        end
            
        function xyz = getEdgeInteriorNodeCoordinates(obj, iEdge, varargin)
            % Get ordered [x,y] coordinates of interior nodes on a given edge
            %
            % getEdgeInteriorNodeCoordinates(iEdge)
            % getEdgeInteriorNodeCoordinates(iEdge, orientation)
            
            if obj.hFieldNodes.N < 3
                xyz = zeros(0,2);
                return
            end
            
            rField = obj.hFieldNodes.basis1d.getInteriorNodes();
            xyz = obj.getEdgeCoordinates(iEdge, rField, varargin{:});
        end
        
        function xy = getEdgeNodeCoordinates(obj, iEdge, varargin)
            % Get ordered [x,y] coordinates of nodes on a given edge
            %
            % getEdgeNodeCoordinates(iEdge)
            % getEdgeNodeCoordinates(iEdge, orientation)
            
            rField = obj.hFieldNodes.basis1d.getNodes();
            xy = obj.getEdgeCoordinates(iEdge, rField, varargin{:});
        end
        
        
        
        function xy = getFaceInteriorNodeCoordinates(obj, iFace)
            % Get ordered [x,y] coordinates of interior nodes in a given face
            %
            % getFaceInteriorNodeCoordinates(iFace)
            
            if obj.hFieldNodes.N < 4
                xy = zeros(0,2);
                return
            end
            
            rs = obj.hFieldNodes.basis.getInteriorNodes();
            xy = obj.getFaceCoordinates(iFace, rs(:,1), rs(:,2));
            
        end
        
        
        function xy = getFaceNodeCoordinates(obj, iFace)
            % Get ordered [x,y] coordinates of all nodes in a given face
            %
            % getFaceNodeCoordinates(iFace)
            
            rs = obj.hFieldNodes.basis.getNodes();
            xy = obj.getFaceCoordinates(iFace, rs(:,1), rs(:,2));
        end
        
        
        function xy = getNodeCoordinates(obj)
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
            
            numNodes = obj.hFieldNodes.getNumNodes();
            xy = zeros(numNodes,2);

            % Nodes, section 1/3: Vertices
            xy(1:obj.hMesh.getNumVertices(),:) = obj.xyNodes(1:obj.hMesh.getNumVertices(),:);
            
            % Nodes, section 2/3: Edge-centers
            for iEdge = 1:obj.hMesh.getNumEdges()
                xy(obj.hFieldNodes.getEdgeInteriorNodes(iEdge),:) = obj.getEdgeInteriorNodeCoordinates(iEdge);
            end

            % Nodes, section 3/3: Face-centers
            for iFace = 1:obj.hMesh.getNumFaces()
                xy(obj.hFieldNodes.getFaceNodes(iFace),:) = obj.getFaceNodeCoordinates(iFace);
            end
        end
        
        function xy = getBoundaryNodeCoordinates(obj)
            xy = obj.getNodeCoordinates();
            xy = xy(obj.hFieldNodes.getBoundaryNodes(),:);
        end
        
        function xy = getInteriorNodeCoordinates(obj)
            xy = obj.getNodeCoordinates();
            xy = xy(obj.hFieldNodes.getInteriorNodes(),:);
        end
        
        % ---- JACOBIANS
        
        function [dxy_dr, dxy_ds] = getJacobian(obj, iFace, rr, ss)
            % Calculate the Jacobian of the mapping from (r,s) to (x,y).
            %
            % For a single point (r,s), the Jacobian is
            %   [Dr*x, Ds*x; Dr*y, Ds*y].
            % where (x,y) are the geometry node coordinates.
            
            % Multiply geometry nodal (x,y).
            xy = obj.xyNodes(obj.hGeomNodes.getFaceNodes(iFace),:);
            
            % Get gradient matrices for geom nodes
            [Dr, Ds] = obj.hGeomNodes.basis.gradientMatrix_rs(rr,ss);
            
            dxy_dr = Dr*xy;
            dxy_ds = Ds*xy;
        end
        
        function dxy_dr = getEdgeJacobian(obj, iEdge, rr, varargin)
            % Calculate the Jacobian of the mapping from r to (x,y).
            %
            % dxy_dr = getEdgeJacobian(obj, iEdge, rr)
            % dxy_dr = getEdgeJacobian(obj, iEdge, rr, orientation)
            
            if nargin < 4
                orientation = 1;
            else
                orientation = varargin{1};
            end
            
            xy = obj.xyNodes(obj.hGeomNodes.getEdgeNodes(iEdge, orientation), :);
            
            % Get gradient matrix for geom nodes
            Dr = obj.hGeomNodes.basis1d.gradientMatrix_rs(rr);
            
            dxy_dr = Dr*xy;
        end
        
        function [dxy_dr, dxy_ds] = getFieldJacobian(obj, iFace)
            
            % Evaluate AT field nodal (r,s).
            rr = obj.hFieldNodes.basis.r;
            ss = obj.hFieldNodes.basis.s;
            
            [dxy_dr, dxy_ds] = obj.getJacobian(iFace, rr, ss);
        end
        
        function [dxy_dr, dxy_ds] = getQuadratureJacobian(obj, iFace)
            
            % Evaluate AT field nodal (r,s).
            rr = obj.hQuadNodes.basis.r;
            ss = obj.hQuadNodes.basis.s;
            
            [dxy_dr, dxy_ds] = obj.getJacobian(iFace, rr, ss);
        end
        
        
        
        % ---- OPERATORS
        
        
        
        function [outDx, outDy, count] = getGradientOperators(obj)
            
            numNodes = obj.hNodes.getNumNodes();
            outDx = sparse(numNodes, numNodes);
            outDy = sparse(numNodes, numNodes);
            
            count = zeros(numNodes, 1);
            
            numFaces = obj.hMesh.getNumFaces();
            
            rs = obj.hNodes.basis.getNodes();
            [Dr, Ds] = support2d.gradients(obj.hNodes.N, rs(:,1), rs(:,2));
            
            for ff = 1:numFaces
                jac = obj.getLinearJacobian(ff);
                
                invJac = inv(jac);
                
                Dx = Dr*invJac(1,1) + Ds*invJac(2,1);
                Dy = Dr*invJac(1,2) + Ds*invJac(2,2);
                
                
                % first silly approach: on edges between faces we will
                % repeatedly overwrite the matrix elements.  Bogus!!
                iGlobal = obj.hNodes.getFaceNodes(ff);
                outDx(iGlobal, iGlobal) = outDx(iGlobal, iGlobal) + Dx;
                outDy(iGlobal, iGlobal) = outDy(iGlobal, iGlobal) + Dy;
                
                count(iGlobal) = count(iGlobal) + 1;
            end
            
            % this handles averaging on boundaries
            normalizer = spdiags(1./count, 0, numNodes, numNodes);
            outDx = normalizer * outDx;
            outDy = normalizer * outDy;
        end
        
        function [outI] = getInterpolationOperator(obj, xs, ys)
            
            numPts = length(xs);
            numNodes = obj.hNodes.getNumNodes();
            outI = sparse(numPts, numNodes);
            count = zeros(numPts, 1);
            
            if numPts == 0
                return
            end
            
            tr = triangulation(obj.hMesh.getFaceVertices(), obj.vertices);
            iEnclosingFaces = tr.pointLocation(xs(:), ys(:));
            numFaces = obj.hMesh.getNumFaces();
            
            for ff = 1:numFaces
                iPoint = find(iEnclosingFaces == ff);
                
                if isempty(iPoint)
                    continue
                end
                
                xy = [xs(iPoint)'; ys(iPoint)'];
                
                xyTri = obj.vertices(obj.hMesh.getFaceVertices(ff), :)';
                %rs = support2d.xy2rs(xyTri, xy);
                M = obj.hNodes.basis.interpolationMatrix_xy(xyTri, xy(1,:), xy(2,:));
                %M = obj.basis.interpolationMatrix_rs(rs(1,:), rs(2,:));
                
                iGlobal = obj.hNodes.getFaceNodes(ff);
                
                outI(iPoint, iGlobal) = M;
                % Original idea: sum and then keep a count to handle the
                % averaging on boundaries.  I want to be simpler.
                %outI(iPoint, iGlobal) = outI(iPoint, iGlobal) + M;
                %count(iPoint) = count(iPoint)+1;
                
            end

            % This handles averaging on boundaries.
            % Usually each output point will inside only one face.
            % However in case a point is on a boundary between faces
            % ...
            % Instead of summing up above let's just overwrite.
            %normalizer = spdiags(1./count, 0, numPts, numPts);
            %outI = normalizer * outI;
            
        end
        
        
        function dIdv = getInterpolationOperatorSensitivity(obj, xs, ys)
            % dIdv = getInterpolationOperatorSensitivity(xs, ys)
            %
            
            numPts = length(xs);
            numNodes = obj.hNodes.getNumNodes();
            numVertices = obj.hMesh.getNumVertices();
            numFaces = obj.hMesh.getNumFaces();
            
            if numPts == 0
                return
            end
            
            dIdv = cell(numVertices, 2);
            for nn = 1:numel(dIdv)
                dIdv{nn} = sparse(numPts, numNodes);
            end
            
            %count = cell(numVertices,2);
            %for nn = 1:numel(count)
            %    count{nn} = zeros(numPts, 1);
            %end
            
            tr = triangulation(obj.hMesh.getFaceVertices(), obj.vertices);
            iEnclosingFaces = tr.pointLocation(xs(:), ys(:));
            
            for ff = 1:numFaces
                iPoint = find(iEnclosingFaces == ff);
                if isempty(iPoint)
                    continue
                end
                iGlobal = obj.hNodes.getFaceNodes(ff);
                iGlobalVertices = obj.hMesh.getFaceVertices(ff);
                
                xy = [xs(iPoint)'; ys(iPoint)'];
                xyTri = obj.vertices(obj.hMesh.getFaceVertices(ff),:)';
                
                %M = obj.basis.interpolationMatrix_xy(xyTri, xy(1,:), xy(2,:));
                
                for iVert = 1:3
                    iGlobalVert = iGlobalVertices(iVert);
                    for iXY = 1:2
                        DxyTri = zeros(size(xyTri));
                        DxyTri(iXY, iVert) = 1;
                        
                        DM = obj.hNodes.basis.interpolationMatrixSensitivity_nodal_xy(xyTri, DxyTri, xy(1,:), xy(2,:));
                        dIdv{iGlobalVert,iXY}(iPoint, iGlobal) = DM;
                        %dIdv{iVert,iXY}(iPoint, iGlobal) = dIdv{iVert,iXY}(iPoint, iGlobal) + DM;
                        %count{iVert,iXY}(iPoint) = count{iVert,iXY}(iPoint) + 1;
                    end
                end
                
            end
            
            %for nn = 1:numel(dIdv)
            %    % this handles averaging on boundaries
            %    normalizer = spdiags(1./count{nn}, 0, numPts, numPts);
            %    dIdv{nn} = normalizer * dIdv{nn};
            %end
            
        end
        
        
        % ---- VISUALIZATION
        %
        
        function plotMesh(obj, varargin)
            
            numEdges = obj.hMesh.getNumEdges();
            rr = linspace(-1, 1, 15);
            
            for ee = 1:numEdges
                xy = obj.getEdgeCoordinates(ee, rr);
                
                line(xy(:,1), xy(:,2), varargin{:});
            end
            
        end
        
        
        function plotMatrix(obj, A, varargin)
            % Draw every FEM node affected by this matrix.
            
            assert(size(A,1) == obj.hNodes.getNumNodes());
            
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
        
        function plotEdgeIndices(obj)
            
            edges = obj.hMesh.getEdgeVertices();
            
            for ii = 1:size(edges,1)
                
                v0 = obj.vertices(edges(ii,1),:);
                v1 = obj.vertices(edges(ii,2),:);
                
                vCenter = 0.5*(v0+v1);
                
                text(vCenter(1), vCenter(2), num2str(ii));
                
            end
        end
        
        function plotFaceIndices(obj)
            
            for ff = 1:obj.hMesh.getNumFaces()
                vCenter = mean(obj.vertices(obj.hMesh.getFaceVertices(ff),:), 1);
                text(vCenter(1), vCenter(2), num2str(ff));
            end
        end
        
        
    end % methods
    
    
end