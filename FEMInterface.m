classdef FEMInterface < handle
    
    properties
        %contours;
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
%            obj.contours = [];
            obj.freeChargeFunction = [];
        end
        
        function setFreeCharge(obj, freeChargeFunc)
            assert(isa(freeChargeFunc, 'function_handle'));
            
            obj.freeChargeFunction = freeChargeFunc;
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
        
        function [iDirichlet, iNeumann, dirichletVals, neumannVals, dirichletJacobian, neumannJacobian] = getBoundaryConditions(obj, p, tnMesh, paramGeom, fieldNodeContour)
            
            % Set Dirichlet and Neumann boundary conditions
            iDirichlet = [];
            iNeumann = [];
            
            contourValues = {};
            contourJacobians = {};
            contourNodes = {};
            
            xyFieldNodes = tnMesh.getNodeCoordinates();
            for cc = 1:length(paramGeom.contours)
                iContourNodes = find(fieldNodeContour == cc);
                
                [bdyVals, bdyJac] = obj.getFunctionJacobian(p, xyFieldNodes(iContourNodes,1), xyFieldNodes(iContourNodes,2), paramGeom.contours(cc).boundaryFunc);
                contourValues{cc} = bdyVals;
                contourJacobians{cc} = bdyJac;
                contourNodes{cc} = iContourNodes;
                
                if strcmpi(paramGeom.contours(cc).type, 'dirichlet')
                    iDirichlet = [iDirichlet; iContourNodes];
                elseif strcmpi(paramGeom.contours(cc).type, 'neumann')
                    iNeumann = [iNeumann; iContourNodes];
                else
                    error('Invalid boundary type %s', paramGeom.contours(cc).type);
                end
            end
            
            dirichletVals = zeros(length(iDirichlet),1);
            neumannVals = zeros(length(iNeumann),1);
            
            dirichletJacobian = sparse(length(iDirichlet), length(p));
            neumannJacobian = sparse(length(iNeumann), length(p));
            
            dIdx = 1;
            nIdx = 1;
            for cc = 1:length(paramGeom.contours)
                contourLength = length(contourNodes{cc});
                
                if strcmpi(paramGeom.contours(cc).type, 'dirichlet')
                    indices = dIdx:(dIdx+contourLength-1);
                    dirichletJacobian(indices,:) = contourJacobians{cc};
                    dirichletVals(indices) = contourValues{cc};
                    dIdx = dIdx + contourLength;
                elseif strcmpi(paramGeom.contours(cc).type, 'neumann')
                    indices = dIdx:(dIdx+contourLength-1);
                    neumannJacobian(indices,:) = contourJacobians{cc};
                    neumannVals(indices) = contourValues{cc};
                    nIdx = nIdx + contourLength;
                else
                    error('Invalid boundary type %s', paramGeom.contours(cc).type);
                end
            end
        end
        
        
        function [femp, geometry, dDirichlet_dp, dnx_dp, dny_dp] = instantiateProblem(obj, p)
            
            geometry = obj.evaluateGeometry(p);
            
            gmsh = obj.meshGeometry(geometry.contourVertices, geometry.contourMeshSizes);
            
            [femp, geometry, dDirichlet_dp, dnx_dp, dny_dp] = obj.instantiateProblemWithMesh(p, geometry, gmsh);
        
        end
        
        
        
        function [femp, dDirichlet_dp, dnx_dp, dny_dp] = instantiateProblemNew(obj, p, ig, paramGeom)  % gmsh is needed for topology and labels.
            
            % Geometry node parameter sensitivity
            dnx_dp = ig.geomNodeJacobian * ig.dvx_dp;
            dny_dp = ig.geomNodeJacobian * ig.dvy_dp;
            
            % Assign contour indices to all field nodes on boundaries.
            % We'll use this to impose Dirichlet and Neumann conditions.
            fieldNodeContour = ig.getFieldNodeContours();
            
            [iDirichlet, iNeumann, dirichletVals, neumannVals, dDirichlet_dp, dNeumann_dp] = obj.getBoundaryConditions(p, ig.tnMesh, paramGeom, fieldNodeContour);
            
            
            % Assemble the FEM problem
            poi = PoissonFEM2D(ig.tnMesh);
            femp = FEMProblem(poi);
            
            femp.setDirichlet(iDirichlet, dirichletVals);
            femp.setNeumann(iNeumann, neumannVals);
            femp.setFreeCharge(@(x,y) obj.freeChargeFunction(p, x, y));
        end
        
        %function [femp, dDirichlet_dp, dnx_dp, dny_dp] = reinstantiateProblem(obj, p, ig, paramGeom)
        %    
        %end
        
        
        function [femp, geometry, dDirichlet_dp, dnx_dp, dny_dp] = instantiateProblemWithMesh(obj, p, ig)  % gmsh is needed for topology and labels.
            lng = LinearNodalGeometry(gmsh.faces, gmsh.vertices, obj.N_geom);
            xyGeomNodes = lng.getNodeCoordinates();
            tnMesh = TriNodalMesh(gmsh.faces, xyGeomNodes, obj.N_field, obj.N_geom, obj.N_quad);
            
            
            % Assign geometry line indices to all geometry nodes on
            % boundaries.  We'll use this for sensitivity calculations.
            geomNodeLine = obj.getGeomNodeLines(tnMesh, gmsh.boundaryEdges, gmsh.edgeGeometryLines);
            
            % Derivative of geometry nodes with respect to geometry
            % vertices (user vertices).
            [dvx_dp, dvy_dp] = obj.evaluateGeometryJacobian(p); % ig.dvx_dp, ig.dvy_dp
            dNode_dv = obj.getGeometryNodeJacobians(geometry, tnMesh, geomNodeLine); % ig.geomNodeJacobian
            
            dnx_dp = dNode_dv * dvx_dp;
            dny_dp = dNode_dv * dvy_dp;
            
            
            % Assign contour indices to all field nodes on boundaries.
            % We'll use this to impose Dirichlet and Neumann conditions.
            fieldNodeContour = obj.getFieldNodeContours(tnMesh, gmsh.boundaryEdges, gmsh.edgeGeometryContours);
            
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