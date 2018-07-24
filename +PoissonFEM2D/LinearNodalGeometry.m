classdef LinearNodalGeometry < handle
    
    properties
        vertices;
        
        hMesh%@MeshTopology;
        hNodes%@NodalTopology;
    end
    
    
    methods
        
        % ---- CONSTRUCTOR

        function obj = LinearNodalGeometry(faceVertices, vertices, N)
            obj.hMesh = PoissonFEM2D.MeshTopology(faceVertices);
            obj.hNodes = PoissonFEM2D.NodalTopology(obj.hMesh, N);
            obj.vertices = vertices;
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
            xy = PoissonFEM2D.support2d.rs2xy(threeVertices', obj.hNodes.basis.getInteriorNodes()')';
        end
        
        
        function xy = getFaceNodeCoordinates(obj, iFace)
            % Get ordered [x,y] coordinates of all nodes in a given face
            %
            % getFaceNodeCoordinates(iFace)
            
            threeVertices = obj.vertices(obj.hMesh.getFaceVertices(iFace),:);
            xy = PoissonFEM2D.support2d.rs2xy(threeVertices', obj.hNodes.basis.getNodes()')';
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
            xy = xy(obj.hNodes.getBoundaryNodes(),:);
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

                            [DT, Dx0] = PoissonFEM2D.support2d.rs2xy_affineParameterSensitivities(xyTri, DxyTri);

                            Dxy = bsxfun(@plus, DT*rs, Dx0);
                            dxy_dv{iGlobalVert,iXY}(iNodes,:) = Dxy';
                        end
                    end
                end
            end % obj.N > 3
            
        end
        
    end % methods
    
end

