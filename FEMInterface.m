classdef FEMInterface < handle
    
    properties
        contours;
        freeChargeFunction;
        
        femProblem;
    end
    
    methods
        
        function obj = FEMInterface()
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
        
        function [mshFaceVertices, mshVerts] = solve(obj, p, objectiveFunc)
            
            col = @(A) reshape(A, [], 1);
            
            contourVertices = {};
            meshSizes = {};
            
            for cc = 1:length(obj.contours)
                x = obj.contours(cc).xFunc(p);
                y = obj.contours(cc).yFunc(p);
                meshSize = obj.contours(cc).meshSizeFunc(p);
                
                if length(meshSize) == 1
                    meshSize = repmat(meshSize, length(x), 1);
                end
                
                contourVertices{cc} = [col(x), col(y)];
                meshSizes{cc} = meshSize;
            end
            
            writeGEO('fromMatlab.geo', contourVertices, meshSizes);
            !/usr/local/bin/gmsh -2 fromMatlab.geo > gmshOut.txt
            [mshFaceVertices, mshEdgeVertices, mshVerts, mshEdgeContour, mshEdgeLine] = readMSH('fromMatlab.msh');
            
            %figure(101); clf
            %patch('Faces', mshFaceVertices, 'Vertices', mshVerts, 'FaceColor', 'r');
            %axis xy image
            %hold on
            
            N_field = 4;
            N_geom = 2;
            N_quad = N_field;
            
            lng = LinearNodalGeometry(mshFaceVertices, mshVerts, N_geom);
            xyGeomNodes = lng.getNodeCoordinates();
            tnMesh = TriNodalMesh(mshFaceVertices, xyGeomNodes, N_field, N_geom, N_quad);
            poi = PoissonFEM2D(tnMesh);
            femp = FEMProblem(poi);
            
            % Now I need to make my nodal mesh and determine which nodes
            % are Dirichlet and Neumann and apply the functions to them.
            
            % I need maps:
            % nodeContour(ii) = contour index of field node ii
            %
            % lineVertices(ii) = geometry vertices of line ii (lines must
            % be inferred from the contours as I go)
            %
            % nodeLine(ii) = line index of geometry node ii.  I can
            % overwrite, it will be ok.
            
            % FIELD NODE CONTOUR
            numFieldNodes = tnMesh.hFieldNodes.getNumNodes();
            fieldNodeContour = zeros(numFieldNodes, 1);
            for cc = 1:length(obj.contours)
                iContourVertices = mshEdgeVertices(mshEdgeContour == cc,:);
                iContourEdges = tnMesh.hMesh.getVertexEdgesExclusive(unique(iContourVertices(:)));
                iContourNodes = tnMesh.hFieldNodes.getEdgeNodes(iContourEdges);
                fieldNodeContour(iContourNodes) = cc;
            end
            
            xyFieldNodes = tnMesh.getNodeCoordinates();
            
            % Plot field nodes.
            %for cc = 1:length(obj.contours)
            %    ii = find(fieldNodeContour == cc);
            %    
            %    plot(xyFieldNodes(ii,1), xyFieldNodes(ii,2), 'o');
            %end
            
            % Set Dirichlet and Neumann boundary conditions
            iDirichlet = [];
            iNeumann = [];
            
            dirichletVals = [];
            neumannVals = [];
            
            for cc = 1:length(obj.contours)
                iContourNodes = find(fieldNodeContour == cc);
                xx = xyFieldNodes(iContourNodes,1);
                yy = xyFieldNodes(iContourNodes,2);
                bdyVals = zeros(length(iContourNodes),1);
                
                for nn = 1:length(bdyVals)
                    bdyVals(nn) = obj.contours(cc).boundaryFunc(p, xx(nn), yy(nn));
                end
                
                if strcmpi(obj.contours(cc).type, 'dirichlet')
                    iDirichlet = [iDirichlet; iContourNodes];
                    dirichletVals = [dirichletVals; bdyVals];
                elseif strcmpi(obj.contours(cc).type, 'neumann')
                    iNeumann = [iNeumann; iContourNodes];
                    neumannVals = [neumannVals; bdyVals];
                else
                    error('Invalid boundary type %s', obj.contours(cc).type);
                end
                
            end
            
            femp.setDirichlet(iDirichlet, dirichletVals);
            femp.setNeumann(iNeumann, neumannVals);
            femp.setFreeCharge(@(x,y) obj.freeChargeFunction(p, x, y));
            
            femp.solve(objectiveFunc);
            
            obj.femProblem = femp;
        end
        
        function solveAdjoint(obj)
            
        end
        
    end % methods
    
end