classdef FEMInterface < handle
    
    properties
        contours;
        freeChargeFunction;
        
        N_field;
        N_geom;
        N_quad;
    end
    
    methods
        
        function obj = FEMInterface(N_field, N_geom, N_quad)
            obj.N_field = N_field;
            obj.N_geom = N_geom;
            obj.N_quad = N_quad;
            obj.contours = [];
            obj.freeChargeFunction = [];
        end
        
        function addContour(obj, x, y, meshSize, dirichletOrNeumann, boundaryFunction)
            
            assert(isa(x, 'function_handle'));
            assert(isa(y, 'function_handle'));
            assert(isa(meshSize, 'function_handle'));
            assert(isa(boundaryFunction, 'function_handle'));
            
            c = struct('xFunc', x, 'yFunc', y, 'meshSizeFunc', meshSize, ...
                'type', dirichletOrNeumann, 'boundaryFunc', boundaryFunction);
            
            if isempty(obj.contours)
                obj.contours = c;
            else
                obj.contours(end+1) = c;
            end
        end
        
        function setFreeCharge(obj, freeChargeFunc)
            assert(isa(freeChargeFunc, 'function_handle'));
            
            obj.freeChargeFunction = freeChargeFunc;
        end
        
        function g = evaluateGeometry(obj, p)
            
            col = @(A) reshape(A, [], 1);
            
            outContourVertices = {};
            outMeshSizes = {};
            
            outGeometryVertices = [];
            outGeometryLines = [];
            idxGeomVert = 1;
            
            for cc = 1:length(obj.contours)
                x = obj.contours(cc).xFunc(p); %local
                y = obj.contours(cc).yFunc(p); %local
                contourLength = length(x); %local
                meshSize = obj.contours(cc).meshSizeFunc(p); %local
                
                if length(meshSize) == 1
                    meshSize = repmat(meshSize, length(x), 1);
                end
                
                outContourVertices{cc} = [col(x), col(y)];
                outMeshSizes{cc} = meshSize;
                
                outGeometryVertices = [outGeometryVertices; outContourVertices{cc}];
                
                arange = (0:(contourLength-1))'; %local
                
                outGeometryLines = [outGeometryLines; idxGeomVert + [arange, mod(arange+1, contourLength)]];
                idxGeomVert = idxGeomVert + contourLength; % local
            end
            
            g = struct('contourVertices', {outContourVertices}, ...
                'contourMeshSizes', {outMeshSizes}, ...
                'vertices', outGeometryVertices, ...
                'lines', outGeometryLines);
        end
        
        function [dvx_dp, dvy_dp] = evaluateGeometryJacobian(obj, p)
            
            numParams = length(p);
            
            rows = [];
            cols = [];
            vx_vals = [];
            vy_vals = [];
            
            delta = 1e-8;
            for nn = 1:numParams
                pLow = p;
                pLow(nn) = pLow(nn) - delta;
                
                pHigh = p;
                pHigh(nn) = pHigh(nn) + delta;
                
                vLow = obj.evaluateGeometry(pLow).vertices;
                vHigh = obj.evaluateGeometry(pHigh).vertices;
                
                numVertices = size(vLow,1);
                rows = [rows; (1:numVertices)'];
                cols = [cols; repmat(nn, numVertices, 1)];
                vx_vals = [vx_vals; (vHigh(:,1)-vLow(:,1))/(2*delta)];
                vy_vals = [vy_vals; (vHigh(:,2)-vLow(:,2))/(2*delta)];
            end
            
            dvx_dp = sparse(rows, cols, vx_vals, numVertices, numParams);
            dvy_dp = sparse(rows, cols, vy_vals, numVertices, numParams);
        end
        
        function gmshObj = meshGeometry(obj, contourVertices, contourMeshSizes)
            
            writeGEO('fromMatlab.geo', contourVertices, contourMeshSizes);
            !/usr/local/bin/gmsh -2 fromMatlab.geo > gmshOut.txt
            [mshFaceVertices, mshEdgeVertices, mshVerts, mshEdgeContour, mshEdgeLine] = readMSH('fromMatlab.msh');
            
            gmshObj = struct('faces', mshFaceVertices,...
                'boundaryEdges', mshEdgeVertices, ...
                'vertices', mshVerts, ...
                'edgeGeometryContours', mshEdgeContour, ...
                'edgeGeometryLines', mshEdgeLine);
            
        end
        
        function outFieldNodeContour = getFieldNodeContours(obj, tnMesh, boundaryEdges, edgeGeometryContours)
            % Get the contour index for each field node.  Only boundary
            % nodes will have a contour.
            
            numFieldNodes = tnMesh.hFieldNodes.getNumNodes(); % local
            outFieldNodeContour = sparse(numFieldNodes, 1); % Needed for boundary conditions!!
            for cc = 1:length(obj.contours)
                iContourVertices = boundaryEdges(edgeGeometryContours == cc,:); % local
                iContourEdges = tnMesh.hMesh.getVertexEdgesExclusive(unique(iContourVertices(:))); % local
                iContourNodes = tnMesh.hFieldNodes.getEdgeNodes(iContourEdges); % local
                outFieldNodeContour(iContourNodes) = cc;
            end
        end
        
        function outGeomNodeLine = getGeomNodeLines(obj, tnMesh, boundaryEdges, edgeGeometryLines)
            
            % GEOMETRY NODE LINE
            xyGeomNodes = tnMesh.xyNodes;
            numGeomNodes = tnMesh.hGeomNodes.getNumNodes();
            outGeomNodeLine = sparse(numGeomNodes, 1);
            numLines = numel(unique(edgeGeometryLines));
            for ll = 1:numLines
                iLineVertices = boundaryEdges(edgeGeometryLines == ll,:);
                iLineEdges = tnMesh.hMesh.getVertexEdgesExclusive(unique(iLineVertices(:)));
                iContourNodes = tnMesh.hGeomNodes.getEdgeNodes(iLineEdges);
                outGeomNodeLine(iContourNodes) = ll;
            end
            
        end
        
        function jac = getGeometryNodeJacobians(obj, geometry, tnMesh, geomNodeLine)
            
            jac = sparse(tnMesh.hGeomNodes.getNumNodes(), size(geometry.vertices, 1));
            
            for ll = 1:size(geometry.lines,1)
                
                ii = find(geomNodeLine == ll);
                %plot(xyGeomNodes(ii,1), xyGeomNodes(ii,2), 'o');
                
                iLineVert0 = geometry.lines(ll,1);
                iLineVert1 = geometry.lines(ll,2);
                xy0 = geometry.vertices(iLineVert0,:);
                xy1 = geometry.vertices(iLineVert1,:);
                
                % Node positions = xy0 + alpha*xy1
                % Let's do this a lazy way
                alpha = [xy1' - xy0'] \ [ tnMesh.xyNodes(ii,1)' - xy0(1); tnMesh.xyNodes(ii,2)' - xy0(2) ];
                
                jac(ii, iLineVert0) = 1 - alpha;
                jac(ii, iLineVert1) = alpha;
            end
            
        end
        
        function [funcVals, jac] = getFunctionJacobian(obj, p, xx, yy, func)
            % Evaluate Jacobian of parameterized source function e.g.
            % dirichlet(p,x,y), neumann(p,x,y), freeCharge(p,x,y).
            
            jac = sparse(length(xx), length(p));
            delta = 1e-8;
            for nn = 1:length(xx)
                
                funcVals(nn) = func(p,xx,yy);
                
                for mm = 1:length(p)
                    pLow = p;
                    pLow(mm) = pLow(mm) - delta;
                    pHigh = p;
                    pHigh(mm) = pHigh(mm) + delta;
                    
                    jac(nn,mm) = (func(pHigh,xx(nn),yy(nn)) - func(pLow,xx(nn),yy(nn))) / (2*delta);
                end
            end
        end
        
        function [iDirichlet, iNeumann, dirichletVals, neumannVals, dirichletJacobian, neumannJacobian] = getBoundaryConditions(obj, p, tnMesh, fieldNodeContour)
            
            % Set Dirichlet and Neumann boundary conditions
            iDirichlet = [];
            iNeumann = [];
            
            contourValues = {};
            contourJacobians = {};
            contourNodes = {};
            
            xyFieldNodes = tnMesh.getNodeCoordinates();
            for cc = 1:length(obj.contours)
                iContourNodes = find(fieldNodeContour == cc);
                
                [bdyVals, bdyJac] = obj.getFunctionJacobian(p, xyFieldNodes(iContourNodes,1), xyFieldNodes(iContourNodes,2), obj.contours(cc).boundaryFunc);
                contourValues{cc} = bdyVals;
                contourJacobians{cc} = bdyJac;
                contourNodes{cc} = iContourNodes;
                
                if strcmpi(obj.contours(cc).type, 'dirichlet')
                    iDirichlet = [iDirichlet; iContourNodes];
                elseif strcmpi(obj.contours(cc).type, 'neumann')
                    iNeumann = [iNeumann; iContourNodes];
                else
                    error('Invalid boundary type %s', obj.contours(cc).type);
                end
            end
            
            dirichletVals = zeros(length(iDirichlet),1);
            neumannVals = zeros(length(iNeumann),1);
            
            dirichletJacobian = sparse(length(iDirichlet), length(p));
            neumannJacobian = sparse(length(iNeumann), length(p));
            
            dIdx = 1;
            nIdx = 1;
            for cc = 1:length(obj.contours)
                contourLength = length(contourNodes{cc});
                
                if strcmpi(obj.contours(cc).type, 'dirichlet')
                    indices = dIdx:(dIdx+contourLength-1);
                    dirichletJacobian(indices,:) = contourJacobians{cc};
                    dirichletVals(indices) = contourValues{cc};
                    dIdx = dIdx + contourLength;
                elseif strcmpi(obj.contours(cc).type, 'neumann')
                    indices = dIdx:(dIdx+contourLength-1);
                    neumannJacobian(indices,:) = contourJacobians{cc};
                    neumannVals(indices) = contourValues{cc};
                    nIdx = nIdx + contourLength;
                else
                    error('Invalid boundary type %s', obj.contours(cc).type);
                end
            end
        end
        
        
        function [femp, geometry, dDirichlet_dp, dnx_dp, dny_dp] = instantiateProblem(obj, p)
            
            geometry = obj.evaluateGeometry(p);
            
            gmsh = obj.meshGeometry(geometry.contourVertices, geometry.contourMeshSizes);
            
            lng = LinearNodalGeometry(gmsh.faces, gmsh.vertices, obj.N_geom);
            xyGeomNodes = lng.getNodeCoordinates();
            tnMesh = TriNodalMesh(gmsh.faces, xyGeomNodes, obj.N_field, obj.N_geom, obj.N_quad);
            
            % Assign contour indices to all field nodes on boundaries.
            % We'll use this to impose Dirichlet and Neumann conditions.
            fieldNodeContour = obj.getFieldNodeContours(tnMesh, gmsh.boundaryEdges, gmsh.edgeGeometryContours);
            
            % Assign geometry line indices to all geometry nodes on
            % boundaries.  We'll use this for sensitivity calculations.
            geomNodeLine = obj.getGeomNodeLines(tnMesh, gmsh.boundaryEdges, gmsh.edgeGeometryLines);
            
            % Derivative of geometry nodes with respect to geometry
            % vertices (user vertices)
            [dvx_dp, dvy_dp] = obj.evaluateGeometryJacobian(p);
            dNode_dv = obj.getGeometryNodeJacobians(geometry, tnMesh, geomNodeLine);
            
            dnx_dp = dNode_dv * dvx_dp;
            dny_dp = dNode_dv * dvy_dp;
            
            [iDirichlet, iNeumann, dirichletVals, neumannVals, dDirichlet_dp, dNeumann_dp] = obj.getBoundaryConditions(p, tnMesh, fieldNodeContour);
            
            % Assemble the FEM problem
            poi = PoissonFEM2D(tnMesh);
            femp = FEMProblem(poi);
            
            femp.setDirichlet(iDirichlet, dirichletVals);
            femp.setNeumann(iNeumann, neumannVals);
            femp.setFreeCharge(@(x,y) obj.freeChargeFunction(p, x, y));
        end
        
    end % methods
    
end