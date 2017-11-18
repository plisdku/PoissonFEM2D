classdef TriNodalMesh < handle
% TriNodalMesh Geometry and topology of a triangulated FEM mesh with nodes
    
    properties
        xyNodes; % the vertex positions come first followed by other nodes
        
        hMesh@MeshTopology;
        hFieldNodes@NodalTopology;
        hGeomNodes@NodalTopology;
        hQuadNodes@NodalTopology;
        hMeshNodes@NodalTopology; % this is the N=2 version
        
    end
    
    
    
    methods
        
        % ---- CONSTRUCTOR
        
        
        function obj = TriNodalMesh(faces, xyNodes, N_field, N_geom, N_quad)
            
            assert(N_quad >= N_field);
            
            obj.hMesh = MeshTopology(faces);
            obj.hFieldNodes = NodalTopology(obj.hMesh, N_field);
            obj.hGeomNodes = NodalTopology(obj.hMesh, N_geom);
            obj.hQuadNodes = NodalTopology(obj.hMesh, N_quad);
            obj.hMeshNodes = NodalTopology(obj.hMesh, 2);
            
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
        
        function xy = getLinearFaceCoordinates(obj, iFace, rr, ss)
            % Get ordered [x,y] coordinates at positions on a given face.
            % Treats all elements as triangles.
            % getLinearFaceCoordinates(obj, iFace, rr, ss)
            
            M = obj.hMeshNodes.basis.interpolationMatrix_rs(rr, ss);
            xy = M*obj.xyNodes(obj.hMeshNodes.getFaceNodes(iFace),:);
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
        
        function [dxy_dr, dxy_ds] = getLinearJacobian(obj, iFace, rr, ss)
            
            % Multiply geometry nodal (x,y).
            xy = obj.xyNodes(obj.hMeshNodes.getFaceNodes(iFace),:);
            
            % Get gradient matrices for geom nodes
            [Dr, Ds] = obj.hMeshNodes.basis.gradientMatrix_rs(rr,ss);
            
            dxy_dr = Dr*xy;
            dxy_ds = Ds*xy;
        end
        
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
        
        
        function invJacs = getInverseLinearJacobian(obj, iFace, rr, ss)
            % Get the inverse Jacobian in a face at each desired point.
            % invJacs = obj.getInverseJacobian(obj, iFace, rr, ss)
            %
            % The Jacobians are indexed invJacs(i,j,iNode).
            
            [dxy_dr, dxy_ds] = obj.getLinearJacobian(iFace, rr, ss);
            
            numPositions = length(rr);
            invJacs = zeros(2, 2, numPositions);
            
            for nn = 1:numPositions
                invJacs(:,:,nn) = inv([dxy_dr(nn,:)', dxy_ds(nn,:)']);
            end
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
        
        function rs = linearInverseCoordinateTransform(obj, iFace, xx, yy)
            % rs = linearInverseCoordinateTransform(obj, iFace, xx, yy)
            % Compute inverse coordinate transform as if the element were a
            % straight-sided triangle.  Initial guess for the real
            % coordinate transform.
            
            xyVerts = obj.getVertexNodeCoordinates(obj.hMesh.getFaceVertices(iFace));
            xyOrigin = 0.5*(xyVerts(2,:) + xyVerts(3,:));
            xyEdge1 = 0.5*(xyVerts(2,:) - xyVerts(1,:));
            xyEdge2 = 0.5*(xyVerts(3,:) - xyVerts(1,:));
            
            xyGoal = [reshape(xx, 1, []); reshape(yy, 1, [])];
            
            M = [xyEdge1', xyEdge2'];
            b = bsxfun(@plus, xyGoal, -xyOrigin');

            rsGuess = M \ b;
            rs = rsGuess;
        end
        
        function [rs, bad, outOfBounds, bigSteps] = inverseCoordinateTransform(obj, iFace, xx, yy)
            % 
            xyGoal = [reshape(xx,1,[]); reshape(yy,1,[])];
            
            
            doPlots = 0;
            
            % A note about step sizes.
            % When the algorithm is nearing convergence, I see no reason
            % to not take a full Newton step.  When it's still far from
            % convergence, a smaller step would be more prudent.  I think 
            % Levenberg-Marquardt is trying to do this balancing act on
            % the fly.  Since I know roughly how far I am from the solution
            % I can probably just guess the right max step size and fix
            % it.  So I should set maxStep to something good but let
            % stepScalar = 1.0 all the time.
            
            doLinearize = 0;
            if doLinearize
                rs = obj.linearInverseCoordinateTransform(iFace, xx, yy);
                %rs = 0.5*(rs + 0.5) - 0.5;
                numIters = 100;
                linearitySchedule = linspace(1.0, 0.0, numIters).^2;
                stepScalar = 1.0;
                maxStep = 0.1;
            else
                rs = obj.linearInverseCoordinateTransform(iFace, xx, yy);
                rs = 0.5*(rs + 0.5) - 0.5;
                numIters = 20;
                linearitySchedule = zeros(1, numIters);
                stepScalar = 1.0;
                maxStep = 0.08;
            end
            
            
            if doPlots
                figure(3); clf
                subplot(1,2,1);
                plot([-1,-1], [1,-1], 'b-');
                hold on
                plot([1,-1], [-1,1], 'b-');
                plot([-1,1], [-1,-1], 'b-');
                xlim([-1.5, 1.5]);
                ylim([-1.5, 1.5]);
                
                subplot(1,2,2);
                obj.plotMesh('linewidth', 3);
                hold on
                xlim([min(xx(:)), max(xx(:))]);
                ylim([min(yy(:)), max(yy(:))]);
                %plot(xyGoal(1,:), xyGoal(2,:), 'k.');
            end
            
            % Newton's method.
            % We will try to successively deform the triangle more and more
            % towards its final shape and pull the coordinate solutions
            % with it.
            
            for ii = 1:numIters
                ll = linearitySchedule(ii);
                
                xyNonlinear = obj.getFaceCoordinates(iFace, rs(1,:), rs(2,:))';
                xyLinear = obj.getLinearFaceCoordinates(iFace, rs(1,:), rs(2,:))';
                
                xy = ll*xyLinear + (1-ll)*xyNonlinear;

                if doPlots
                    pause(0.0001);
                    if exist('p1', 'var')
                        delete(p1);
                        delete(p2);
                        delete(p3);
                        %delete(p4);
                    end
                end
                
                if doPlots
                    
                    [dxy_dr, dxy_ds] = obj.getJacobian(iFace, rs(1,:), rs(2,:));
                    detJ = dxy_dr(:,1).*dxy_ds(:,2) - dxy_dr(:,2).*dxy_ds(:,1);
                    
                    iInverted = detJ <= 0;
                    iNotInverted = ~iInverted;
                    
                    iInsideOk = iNotInverted' & rs(1,:) > -1 & rs(2,:) > -1 & rs(1,:) + rs(2,:) < 0;
                    
                    figure(3);
                    subplot(1,2,1)
                    p1 = plot(rs(1,iNotInverted), rs(2,iNotInverted), 'b.', 'markersize', 7);
                    p2 = plot(rs(1,iInverted), rs(2,iInverted), 'r.', 'markersize', 7);
                    
                    subplot(1,2,2);
                    %p1 = plot(xyLinear(1,:), xyLinear(2,:), 'b.');
                    %p2 = plot(xyNonlinear(1,:), xyNonlinear(2,:), 'r.');
                    
                    p3 = plot(xyNonlinear(1,iInsideOk), xyNonlinear(2,iInsideOk), 'b.', 'markersize', 7);
                    %p3 = plot(xy(1,iNotInverted), xy(2,iNotInverted), 'b.', 'markersize', 7);
                    %p4 = plot(xy(1,iInverted), xy(2,iInverted), 'ro');
                end
                
                xyResidual = xy - xyGoal;
                
                [dxy_dr_nl, dxy_ds_nl] = obj.getJacobian(iFace, rs(1,:), rs(2,:));
                [dxy_dr_l, dxy_ds_l] = obj.getLinearJacobian(iFace, rs(1,:), rs(2,:));
                
                dxy_dr = ll*dxy_dr_l + (1-ll)*dxy_dr_nl;
                dxy_ds = ll*dxy_ds_l + (1-ll)*dxy_ds_nl;
                
                rsOld = rs;
                for dd = 1:numel(xx)
                    
                    jac = [dxy_dr(dd,:)', dxy_ds(dd,:)'];
                
                    rsStep = stepScalar*max(min(-jac \ xyResidual(:,dd), maxStep), -maxStep);
                    rs(:,dd) = rs(:,dd) + rsStep;
                end
            end
            
            xy = obj.getFaceCoordinates(iFace, rs(1,:), rs(2,:))';
            %plot(xy(1,:), xy(2,:), 'go');
            
            % Now test which ones are in/out and which ones converged.
            delta = 1e-4;
            
            lastStep = rs - rsOld;
            bigSteps = abs(lastStep(1,:)) > delta | abs(lastStep(2,:)) > delta;
            outOfBounds = rs(1,:) < -1.0-delta | rs(2,:) < -1.0-delta | rs(2,:)+rs(1,:) > delta;
            bad = outOfBounds | bigSteps;
        end
        
        
        % ---- DIFFERENTIATION
        
        function [outDx, outDy] = getFaceGradientMatrices(obj, iFace)
            
            rBasis = obj.hFieldNodes.basis.r;
            sBasis = obj.hFieldNodes.basis.s;
            
            [Dr, Ds] = obj.hFieldNodes.basis.gradientMatrix_rs(rBasis, sBasis);
            
            invJacs = obj.getInverseJacobian(iFace, rBasis, sBasis);
            
            outDx = diag(squish(invJacs(1,1,:)))*Dr + diag(squish(invJacs(2,1,:)))*Ds;
            outDy = diag(squish(invJacs(1,2,:)))*Dr + diag(squish(invJacs(2,2,:)))*Ds;
            
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
            
            for ff = 1:numFaces
                [Dx,Dy] = obj.getFaceGradientMatrices(ff);
                
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
        
        function yesNo = pointInElement(obj, iFace, xs, ys)
            [rs, outOfBounds] = obj.inverseCoordinateTransform(iFace, xs, ys);
            
            yesNo = ~outOfBounds;
        end
        
        function iFace = getEnclosingFaces(obj, xs, ys)
            
            % We'll take a coarse first pass and then correct it.
            
            tr = triangulation(obj.hGeomNodes.getNodalMesh(), obj.xyNodes);
            numSubTris = (obj.hGeomNodes.N-1)^2;
            iFace0 = ceil(tr.pointLocation(xs(:), ys(:)) / numSubTris); % first guess.
            
            if obj.hGeomNodes.N == 2
                % No mistakes to be made if the sides are straight.
                iFace = iFace0;
                return;
            else
                iFace = nan(size(iFace0));
            end
            
            % We have curved sides.  Need to double-check all points.
            
            ffa = obj.hMesh.getFaceFaceAdjacency();
            
            numFaces = obj.hMesh.getNumFaces();
            
            for ff = 1:numFaces
                
                iNeighborFaces = find(ffa(ff,:));
                iNeighborFaces = [];
                
                iPoint = find(iFace0 == ff);
                
                if isempty(iPoint)
                    continue
                end
                %fprintf('Face %i\n', ff);
                xy = [xs(iPoint)'; ys(iPoint)'];
                [rs, bad] = obj.inverseCoordinateTransform(ff, xs(iPoint), ys(iPoint));
                iFace(iPoint(~bad)) = ff;
                
                iOutOfBounds = find(bad);
                
                for nf = iNeighborFaces
                    [rs2, newBad] = obj.inverseCoordinateTransform(nf, xs(iPoint(isBad)), ys(iPoint(isBad)));
                    
                    bad(~newBad) = 0;
                    %fprintf('Fixed %i good from %i bad\n', nnz(isGood), length(isBad));
                    
                    iFace(iPoint(iOutOfBounds(~newBad))) = nf;
                end
                
                % rs(x or y, point index) should be inside the unit
                % triangle.  For all the bad points, we will try to
                % locate them in other triangles.
                %
                % A test of the residual should suffice though.
                
            end
            
            %iFace = iFace0;
            
        end % getEnclosingFaces
        
        % To implement this function I will need to invert the coordinate
        % transformation.
        function [outI] = getInterpolationOperator(obj, xs, ys)
            
            assert(isvector(xs));
            assert(isvector(ys));
            
            numPts = length(xs);
            numNodes = obj.hFieldNodes.getNumNodes();
            outI = sparse(numPts, numNodes);
            count = zeros(numPts, 1);
            
            if numPts == 0
                return
            end
            
            iEnclosingFaces = obj.getEnclosingFaces(xs,ys);
            
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
        
        % ---- POINT-IN-ELEMENT
        
        
        function xy = getFaceBoundary(obj, iFace)
            
            if obj.hGeomNodes.N == 2
                numPointsPerEdge = 2;
            else
                numPointsPerEdge = 3*obj.hGeomNodes.N;
            end
            
            [iFaceEdges, ori] = obj.hMesh.getFaceEdges(iFace);
            rr = linspace(-1, 1, numPointsPerEdge);
            rr = rr(1:end-1);
            
            xy = [];
            
            for ee = 1:3
                
                xyEdge = obj.getEdgeCoordinates(iFaceEdges(ee), rr, ori(ee));
                xy = [xy; xyEdge];
                
            end
            
        end
        
        
        
        % ---- VISUALIZATION
        %
        
        function plotMesh(obj, varargin)
            
            numEdges = obj.hMesh.getNumEdges();
            
            if obj.hGeomNodes.N == 2
                numPointsPerEdge = 2;
            else
                numPointsPerEdge = 10*obj.hGeomNodes.N;
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