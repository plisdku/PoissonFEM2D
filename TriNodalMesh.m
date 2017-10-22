classdef TriNodalMesh < handle
% TriNodalMesh Geometry and topology of a triangulated FEM mesh with nodes
    
    properties
        vertices;
        
        hMesh@MeshTopology;
        hFieldNodes@NodalTopology;
        hGeomNodes@NodalTopology;
        hQuadNodes@NodalTopology;
        
    end
    
    
    
    methods
        
        % ---- CONSTRUCTOR
        
        
        function obj = TriNodalMesh(N_field, N_geom, N_quad, faces, vertices)
            obj.hMesh = MeshTopology(faces);
            obj.hFieldNodes = NodalTopology(obj.hMesh, N_field);
            obj.hGeomNodes = NodalTopology(obj.hMesh, N_geom);
            obj.hQuadNodes = NodalTopology(obj.hMesh, N_quad);
            
            assert(size(vertices,2) == 2, 'Vertices must be Nx2'); % test because 3d verts are common
            obj.vertices = vertices;
        end
        
        % ---- JACOBIANS
        
        function [dxy_dr, dxy_ds] = getJacobian(obj, iFace, rr, ss)
            % Calculate the Jacobian of the mapping from (r,s) to (x,y).
            %
            % For a single point (r,s), the Jacobian is
            %   [Dr*x, Ds*x; Dr*y, Ds*y].
            % where (x,y) are the geometry node coordinates.
            
            % Multiply geometry nodal (x,y).
            xy = obj.vertices(obj.hGeomNodes.getFaceVertexNodes(iFace),:);
            
            % Get gradient matrices for geom nodes
            [Dr, Ds] = obj.hGeomNodes.basis.gradientMatrix_rs(rr,ss);
            
            dxy_dr = Dr*xy;
            dxy_ds = Ds*xy;
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
        
        function jac = getLinearJacobian(obj, iFace)
            % Calculate the Jacobian of the mapping from (r,s) to (x,y).
            % 
            % [T, v0] = support2d.rs2xy_affineParameters(xyTri);
            % xy = bsxfun(@plus, v0, T*rs);
            % so the Jacobian is just T.
            
            threeVertices = obj.vertices(obj.hMesh.getFaceVertices(iFace),:);
            jac = support2d.rs2xy_affineParameters(threeVertices');
        end
        
        function jac1d = getLinearJacobian1d(obj, iEdge, orientation)
            % This is the Jacobian of the mapping from r to (x,y).
            %
            % xy = (v1+v2)/2 + (v2-v1)/2 * r;
            %
            % so d(xy)/dr = (v2-v1)/2.
            
            twoVertices = obj.vertices(obj.hMesh.getEdgeVertices(iEdge, orientation), :);
            jac1d = 0.5*(twoVertices(2,:) - twoVertices(1,:))';
        end
        
        
        function dJdv = getLinearJacobianSensitivity(obj, iFace)
            % dJdv = getLinearJacobianSensitivity(iFace)
            %
            % dJdv is a 4D array indexed by (row, col, vertex, xy).
            %
            % Affine transformations are simple and actually the face
            % argument is unnecessary.
            
            dJdv = zeros(2,2,3,2); % row, col, vertex, vertex xy
            
            threeVertices = obj.vertices(obj.hMesh.getFaceVertices(iFace),:);
            
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
        
        
        function dJdv1d = getLinearJacobianSensitivity1d(obj, iEdge, orientation)
            % dJdv1d = getLinearJacobianSensitivity1d(iEdge, orientation)
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
            
            % uh... i don't think we need to flip the sign because of
            % the orientation...
            %if orientation < 0
            %    dJdv1d = -dJdv1d;
            %end
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
            
            if obj.hNodes.N < 3
                xyz = zeros(0,2);
                return
            end
            
            verts = obj.vertices(obj.hMesh.getEdgeVertices(iEdge),:);
            %v2 = obj.vertices(obj.getEdgeVertices(iEdge),:);
            
            %d = linspace(0, 1, obj.N)';
            %d = d(2:end-1);
            
            d = 0.5 + 0.5*transpose(obj.hNodes.basis1d.getInteriorNodes());
            
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
            
            verts = obj.vertices(obj.hMesh.getEdgeVertices(iEdge),:);
            
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
            threeVertices = obj.vertices(obj.hMesh.getFaceVertices(iFace),:);
            xy = support2d.rs2xy(threeVertices', obj.hNodes.basis.getInteriorNodes()')';
        end
        
        
        function xy = getFaceNodeCoordinates(obj, iFace)
            % Get ordered [x,y] coordinates of all nodes in a given face
            %
            % getFaceNodeCoordinates(iFace)
            
            threeVertices = obj.vertices(obj.hMesh.getFaceVertices(iFace),:);
            xy = support2d.rs2xy(threeVertices', obj.hNodes.basis.getNodes()')';
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
            
            numNodes = obj.hNodes.getNumNodes();
            xy = zeros(numNodes,2);

            % Nodes, section 1/3: Vertices
            xy(1:obj.hMesh.getNumVertices(),:) = obj.vertices;
            
            % Nodes, section 2/3: Edge-centers
            for iEdge = 1:obj.hMesh.getNumEdges()
                xy(obj.hNodes.getEdgeInteriorNodes(iEdge),:) = obj.getEdgeInteriorNodeCoordinates(iEdge);
            end

            % Nodes, section 3/3: Face-centers
            for iFace = 1:obj.hMesh.getNumFaces()
                xy(obj.hNodes.getFaceNodes(iFace),:) = obj.getFaceNodeCoordinates(iFace);
            end
        end
        
        function xy = getBoundaryNodeCoordinates(obj)
            xy = obj.getNodeCoordinates();
            xy = xy(obj.getBoundaryNodes(),:);
        end
        
        function xy = getInteriorNodeCoordinates(obj)
            xy = obj.getNodeCoordinates();
            xy = xy(obj.hNodes.getInteriorNodes(),:);
        end
        
        % ---- Node coordinate sensitivities
        
        function dxy_dv = getNodeCoordinateSensitivities(obj)
            
            numNodes = obj.hNodes.getNumNodes();
            numVertices = obj.hMesh.getNumVertices();
            numEdges = obj.hMesh.getNumEdges();
            numFaces = obj.hMesh.getNumFaces();
            
            dxy_dv = cell(numVertices,2);
            for nn = 1:numel(dxy_dv)
                dxy_dv{nn} = sparse(numNodes, 2);
            end
            
            % Nodes, section 1/3: Vertices
            for iVert = 1:numVertices
                for iXY = 1:2
                    dxy_dv{iVert,iXY}(iVert, iXY) = 1;
                end
            end
            
            % Nodes, section 2/3: Edge-centers
            
            if obj.hNodes.N > 2
                r_1d = obj.hNodes.basis1d.getInteriorNodes();
                for iEdge = 1:numEdges
                    iEdgeVertices = obj.hMesh.getEdgeVertices(iEdge);
                    iv0 = iEdgeVertices(1);
                    iv1 = iEdgeVertices(2);

                    iNodes = obj.hNodes.getEdgeInteriorNodes(iEdge);

                    % Get fractional distance from one end to the other
                    d = 0.5 + 0.5*r_1d;

                    % ... so node = v0(1-d) + v1(d)
                    dxy_dv0 = 1-d;
                    dxy_dv1 = d;

                    dxy_dv{iv0,1}(iNodes,1) = dxy_dv0;
                    dxy_dv{iv0,2}(iNodes,2) = dxy_dv0;
                    dxy_dv{iv1,1}(iNodes,1) = dxy_dv1;
                    dxy_dv{iv1,2}(iNodes,2) = dxy_dv1;
                end
            end % obj.N > 2
            
            % Nodes, section 3/3: Face-centers
            
            if obj.hNodes.N > 3
                rs = obj.hNodes.basis.getInteriorNodes()';
                for iFace = 1:numFaces
                    iFaceVertices = obj.hMesh.getFaceVertices(iFace);

                    xyTri = obj.vertices(iFaceVertices, :)';

                    iNodes = obj.hNodes.getFaceInteriorNodes(iFace);

                    for iLocalVert = 1:3
                        iGlobalVert = iFaceVertices(iLocalVert);
                        for iXY = 1:2
                            DxyTri = zeros(2,3);
                            DxyTri(iXY, iLocalVert) = 1.0;

                            [DT, Dx0] = support2d.rs2xy_affineParameterSensitivities(xyTri, DxyTri);

                            Dxy = bsxfun(@plus, DT*rs, Dx0);
                            dxy_dv{iGlobalVert,iXY}(iNodes,:) = Dxy';
                        end
                    end
                end
            end % obj.N > 3
            
        end
        
    end % methods
    
    
end