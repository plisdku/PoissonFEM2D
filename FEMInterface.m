classdef FEMInterface < handle
    
    properties
        freeChargeFunction;
        instantiatedGeom@InstantiatedGeometry2D;
        boundaryConditions;
    end
    
    methods
        
        function obj = FEMInterface(parameterizedGeom, N_field, N_geom, N_quad, isAxisymmetric)
            obj.instantiatedGeom = InstantiatedGeometry2D(parameterizedGeom, N_field, N_geom, N_quad, isAxisymmetric);
            obj.freeChargeFunction = [];
            obj.boundaryConditions = {};
        end
        
        function setFreeCharge(obj, freeChargeFunc)
            assert(isa(freeChargeFunc, 'function_handle'));
            obj.freeChargeFunction = freeChargeFunc;
        end
        
        function setDirichlet(obj, inLabel, inFunc)
            obj.boundaryConditions{end+1} = struct('label', inLabel, 'func', inFunc, 'type', 'dirichlet');
        end
        
        function setNeumann(obj, inLabel, inFunc)
            obj.boundaryConditions{end+1} = struct('label', inLabel, 'func', inFunc, 'type', 'neumann');
        end
        
        function [funcVals, jac] = getFunctionJacobian(obj, p, xx, yy, func)
            % Evaluate Jacobian of parameterized source function e.g.
            % dirichlet(p,x,y), neumann(p,x,y), freeCharge(p,x,y).
            
            jac = sparse(length(xx), length(p));
            delta = 1e-8;
            for nn = 1:length(xx)
                
                funcVals(nn) = func(p,xx(nn),yy(nn));
                
                for mm = 1:length(p)
                    pLow = p;
                    pLow(mm) = pLow(mm) - delta;
                    pHigh = p;
                    pHigh(mm) = pHigh(mm) + delta;
                    
                    jac(nn,mm) = (func(pHigh,xx(nn),yy(nn)) - func(pLow,xx(nn),yy(nn))) / (2*delta);
                end
            end
        end
        
        function [iDirichlet, iNeumann, iBoundaryNodes, boundaryNodeLabels] = getBoundaryFieldNodeIndices(obj)
            
            %[iBoundaryNodes,~,nodeLabels] = find(obj.instantiatedGeom.getFieldNodeLabels());
            %boundaryDorN = zeros(size(nodeLabels));
            
            numFieldNodes = obj.instantiatedGeom.tnMesh.hFieldNodes.getNumNodes();
            boundaryDorN = sparse(numFieldNodes, 1);
            nodeLabels = sparse(numFieldNodes, 1);
            
            for ii = 1:length(obj.boundaryConditions)
                bcType = obj.boundaryConditions{ii}.type;
                bcLabel = obj.boundaryConditions{ii}.label;
                
                iNodes = obj.instantiatedGeom.getLabeledFieldNodes(bcLabel);
                nodeLabels(iNodes) = bcLabel;
                
                if strcmpi(bcType, 'dirichlet')
                    boundaryDorN(iNodes) = 'd';
                    %boundaryDorN(nodeLabels == bcLabel) = 'd';
                else
                    boundaryDorN(iNodes) = 'n';
                    %boundaryDorN(nodeLabels == bcLabel) = 'n';
                end
            end
            [iBoundaryNodes,~,bcType] = find(boundaryDorN);
            iDirichlet = iBoundaryNodes(bcType == 'd');
            iNeumann = iBoundaryNodes(bcType == 'n');
            [iBoundaryNodes,~,boundaryNodeLabels] = find(nodeLabels);
        end
        
        function [iDirichlet, iNeumann, dirichletVals, neumannVals, dirichletJacobian, neumannJacobian] = ...
                getBoundaryConditions(obj, p)
            
            % Copy the handles out here for convenience.
            paramGeom = obj.instantiatedGeom.parameterizedGeometry;
            tnMesh = obj.instantiatedGeom.tnMesh;
            xyFieldNodes = tnMesh.getNodeCoordinates();
            numFieldNodes = size(xyFieldNodes, 1);
            numParameters = length(p);
            
            % Indices of field nodes that are Dirichlet or Neumann.
            % iBoundary is the union of iDirichlet and iNeumann.
            [iDirichlet, iNeumann, iBoundary, nodeLabels] = obj.getBoundaryFieldNodeIndices();
            
            boundaryVals = sparse(numFieldNodes, 1);
            boundaryJacobian = sparse(numFieldNodes, numParameters);
            
            for ii = 1:length(obj.boundaryConditions)
                bc = obj.boundaryConditions{ii};
                iNodes = iBoundary(nodeLabels == bc.label);
                
                [bdyVals, bdyJac] = obj.getFunctionJacobian(p, xyFieldNodes(iNodes,1), xyFieldNodes(iNodes,2), bc.func);
                
                boundaryVals(iNodes) = bdyVals;
                boundaryJacobian(iNodes,:) = bdyJac;
            end
            
            dirichletVals = boundaryVals(iDirichlet,:);
            dirichletJacobian = boundaryJacobian(iDirichlet,:);
            
            neumannVals = boundaryVals(iNeumann,:);
            neumannJacobian = boundaryJacobian(iNeumann,:);
        end
        
        
        function [femp, dDirichlet_dp, dnx_dp, dny_dp] = instantiateProblem(obj, p)
            obj.instantiatedGeom.instantiateMesh(p);
            
            [femp, dDirichlet_dp, dnx_dp, dny_dp] = obj.assembleProblemWithMesh(p);
     
        end
        
        function [femp, dDirichlet_dp, dnx_dp, dny_dp] = adjustProblem(obj, p)
            obj.instantiatedGeom.adjustMesh(p);
            [femp, dDirichlet_dp, dnx_dp, dny_dp] = obj.assembleProblemWithMesh(p);
        end
        
        function [femp, dDirichlet_dp, dnx_dp, dny_dp] = assembleProblemWithMesh(obj, p)
            
            % Geometry node parameter sensitivity
            dnx_dp = obj.instantiatedGeom.geomNodeJacobian * obj.instantiatedGeom.dvx_dp;
            dny_dp = obj.instantiatedGeom.geomNodeJacobian * obj.instantiatedGeom.dvy_dp;
            
            [iDirichlet, iNeumann, dirichletVals, neumannVals, dDirichlet_dp, dNeumann_dp] = obj.getBoundaryConditions(p);
            
            % Assemble the FEM problem
            poi = PoissonFEM2D(obj.instantiatedGeom.tnMesh);
            
            femp = FEMProblem(poi);
           
            
            femp.setDirichlet(iDirichlet, dirichletVals);
            femp.setNeumann(iNeumann, neumannVals);
            femp.setFreeCharge(@(x,y) obj.freeChargeFunction(p, x, y));
        end
        
        
    end % methods
    
end