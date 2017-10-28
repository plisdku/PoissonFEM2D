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
        
        function invJacs = getInverseJacobian(obj, iFace, rr, ss)
            % Get the inverse Jacobian in a face at each desired point.
            % invJacs = obj.getInverseJacobian(obj, iFace, rr, ss)
            %
            % The Jacobians are indexed invJacs(i,j,iNode).
            
            [dxy_dr, dxy_ds] = obj.getJacobian(iFace, rr, ss);
            
            numPositions = length(rr);
            invJacs = zeros(2, 2, numPositions);
            
            for nn = 1:numPositions
                invJacs(:,:,nn) = inv([dxy_dr(nn,:)', dxy_ds(nn,:)']);
            end
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
        
        % ---- INVERSE COORDINATE TRANSFORM
        
        function rs = inverseCoordinateTransform(obj, iFace, xx, yy)
            
            % Run Newton's method for a while.
            numIters = 10;
            
            xyGoal = [xx, yy]';
            
            rs = 0*xyGoal - 0.5;
            
            %figure(1); clf
            %obj.plotMesh();
            %hold on
            %plot(xyGoal(1,:), xyGoal(2,:), 'rx');
            
            for ii = 1:numIters
                
                xy = obj.getFaceCoordinates(iFace, rs(1,:), rs(2,:))';
                
                %plot(xy(1,:), xy(2,:), 'go');
                
                xyResidual = xy - xyGoal;
                
                invJac = obj.getInverseJacobian(iFace, rs(1,:), rs(2,:));
                
                rsOld = rs;
                for dd = 1:size(invJac,3)
                    rs(:,dd) = rs(:,dd) - invJac(:,:,dd)*xyResidual(:,dd);
                    
                    while rs(1,dd) < -1.0
                        rs(:,dd) = 0.5*( rs(:,dd) + rsOld(:,dd) );
                    end
                    
                    while rs(2,dd) < -1.0
                        rs(:,dd) = 0.5*( rs(:,dd) + rsOld(:,dd) );
                    end
                    
                    while rs(2,dd) > -rs(1,dd)
                        rs(:,dd) = 0.5*( rs(:,dd) + rsOld(:,dd) );
                    end
                end
                
                
                %pause
            end
            %xyResidual = xyGoal - xy;
            
            %disp(norm(xyResidual));
            
        end
        
        % ---- QUADRATURE
        
        function Q = getQuadratureMatrix(obj, iFace)
            % I need Ifq, Vq, and det(J)q.
            % For the Jacobian and interpolation matrix I need the nodes in (r,s).
            
            rsQuad = obj.hQuadNodes.basis.getNodes();
            
            % Interpolation matrix from field nodes to quadrature nodes
            Ifq = obj.hFieldNodes.basis.interpolationMatrix_rs(rsQuad(:,1), rsQuad(:,2));
            
            % Integration kernel on quadrature nodes
            invVq = obj.hQuadNodes.basis.invV;
            Qq = invVq' * invVq;
            
            % Jacobian on quadrature nodes
            [dxy_dr, dxy_ds] = obj.getJacobian(iFace, rsQuad(:,1), rsQuad(:,2));
            detJq = dxy_dr(:,1).*dxy_ds(:,2) - dxy_dr(:,2).*dxy_ds(:,1);
            assert(all(detJq > 0));
            
            Q = Ifq' * Qq * diag(detJq) * Ifq;
        end
        
        
        function Q = getQuadratureMatrix1d(obj, iEdge, varargin)
            if nargin < 3
                orientation = 1;
            else
                orientation = varargin{1};
            end
            
            rQuad = obj.hQuadNodes.basis1d.getNodes();
            if orientation < 0
                rQuad = rQuad(end:-1:1);
            end
            
            Ifq = obj.hFieldNodes.basis1d.interpolationMatrix_r(rQuad);
            
            invVq = obj.hQuadNodes.basis1d.invV;
            Qq = invVq' * invVq;
            
            dxy_dr = obj.getEdgeJacobian(iEdge, rQuad, orientation);
            
            % The determinant should be sqrt(dxy_dr' * dxy_dr).  Right?
            % But I have to sort of do it by hand here.
            detJq = sqrt(dxy_dr(:,1).^2 + dxy_dr(:,2).^2);
            
            Q = Ifq' * Qq * diag(detJq) * Ifq;
            
        end
        
        
        % ---- FULL-MESH OPERATORS
        
        function outQ = getQuadratureOperator(obj)
            
            numFieldNodes = obj.hFieldNodes.getNumNodes();
            outQ = sparse(numFieldNodes, numFieldNodes);
            numFaces = obj.hMesh.getNumFaces();
            
            for ff = 1:numFaces
                %jac = obj.getLinearJacobian(ff);
                
                Q = obj.getQuadratureMatrix(ff);
                
                % first silly approach: on edges between faces we will
                % repeatedly overwrite the matrix elements.  Bogus!!
                iGlobal = obj.hFieldNodes.getFaceNodes(ff);
                outQ(iGlobal, iGlobal) = outQ(iGlobal, iGlobal) + Q;
            end
        end
        
        function [outDx, outDy, count] = getGradientOperators(obj)
            
            numFieldNodes = obj.hFieldNodes.getNumNodes();
            outDx = sparse(numFieldNodes, numFieldNodes);
            outDy = sparse(numFieldNodes, numFieldNodes);
            
            count = zeros(numFieldNodes, 1); % used for averaging on edges.
            
            numFaces = obj.hMesh.getNumFaces();
            
            rs = obj.hFieldNodes.basis.getNodes();
            [Dr, Ds] = obj.hFieldNodes.basis.gradientMatrix_rs(rs(:,1), rs(:,2));
            %[Dr, Ds] = support2d.gradients(obj.hNodes.N, rs(:,1), rs(:,2));
            
            for ff = 1:numFaces
                %jac = obj.getLinearJacobian(ff);
                
                %[dxy_dr, dxy_ds] = obj.getJacobian(ff, rr, ss);
                invJacs = obj.getInverseJacobian(ff, rs(:,1), rs(:,2));
                
                Dx = diag(squish(invJacs(1,1,:)))*Dr + diag(squish(invJacs(2,1,:)))*Ds;
                Dy = diag(squish(invJacs(1,2,:)))*Dr + diag(squish(invJacs(2,2,:)))*Ds;
                
                % first silly approach: on edges between faces we will
                % repeatedly overwrite the matrix elements.  Bogus!!
                iGlobal = obj.hFieldNodes.getFaceNodes(ff);
                outDx(iGlobal, iGlobal) = outDx(iGlobal, iGlobal) + Dx;
                outDy(iGlobal, iGlobal) = outDy(iGlobal, iGlobal) + Dy;
                
                count(iGlobal) = count(iGlobal) + 1;
            end
            
            % this handles averaging on boundaries
            normalizer = spdiags(1./count, 0, numFieldNodes, numFieldNodes);
            outDx = normalizer * outDx;
            outDy = normalizer * outDy;
        end
        
        
        % To implement this function I will need to invert the coordinate
        % transformation.
        function [outI] = getInterpolationOperator(obj, xs, ys)
            
            numPts = length(xs);
            numNodes = obj.hFieldNodes.getNumNodes();
            outI = sparse(numPts, numNodes);
            count = zeros(numPts, 1);
            
            if numPts == 0
                return
            end
            
            tr = triangulation(obj.hGeomNodes.getNodalMesh(), obj.xyNodes);
            numSubTris = (obj.hGeomNodes.N-1)^2;
            iEnclosingFaces = ceil(tr.pointLocation(xs(:), ys(:)) / numSubTris);
            
            %tr = triangulation(obj.hMesh.getFaceVertices(), obj.xyNodes);
            %iEnclosingFaces = tr.pointLocation(xs(:), ys(:));
            numFaces = obj.hMesh.getNumFaces();
            
            for ff = 1:numFaces
                iPoint = find(iEnclosingFaces == ff);
                
                if isempty(iPoint)
                    continue
                end
                
                xy = [xs(iPoint)'; ys(iPoint)'];
                
                rs = obj.inverseCoordinateTransform(ff, xs(iPoint), ys(iPoint));
                M = obj.hFieldNodes.basis.interpolationMatrix_rs(rs(1,:), rs(2,:));
                
                %xyTri = obj.vertices(obj.hMesh.getFaceVertices(ff), :)';
                %M = obj.hNodes.basis.interpolationMatrix_xy(xyTri, xy(1,:), xy(2,:));
                
                iGlobal = obj.hFieldNodes.getFaceNodes(ff);
                
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
            
            if obj.hGeomNodes.N == 2
                numPointsPerEdge = 2;
            else
                numPointsPerEdge = 4*obj.hGeomNodes.N;
            end
            
            rr = linspace(-1, 1, numPointsPerEdge);
            
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