classdef FEMProblem < handle
    
    properties
        poi@PoissonFEM2D;
        
        iDirichlet
        iNeumann
        iCenter
        
        dirichletFunc % unused
        neumannFunc % unused
        chargeFunc % USED (damn)
        
        u;
        v;
        
        F;
        
        u0_dirichlet;
        en_neumann;
        
        freeCharge;
        dFreeCharge_dx;
        dFreeCharge_dy;
        
        A, dA;
        B, dB;
        NM, dNM;
        
        % For Cartesian mode ("only mode")
        uCartesian;
        interpolationOperator;
        dInterpolationOperator; % for input/output
        
        dF_dCharge;
        dF_dDirichlet;
        dF_dNeumann;
        dF_dxy;
        
    end
    
    
    methods
        
        function obj = FEMProblem(poissonFEM)
            obj.poi = poissonFEM;
            
            numNodes = obj.poi.tnMesh.hFieldNodes.getNumNodes();
            
            obj.iDirichlet = [];
            obj.iNeumann = [];
            obj.iCenter = 1:numNodes;
        end
        
        function other = copyModel(obj)
            other = FEMProblem(obj.poi.copy());
            other.iDirichlet = obj.iDirichlet;
            other.iNeumann = obj.iNeumann;
            other.iCenter = obj.iCenter;
            
            other.u0_dirichlet = obj.u0_dirichlet;
            other.en_neumann = obj.en_neumann;
            other.setFreeCharge(obj.chargeFunc);
            
            %other.setSources(obj.chargeFunc, obj.dirichletFunc, obj.neumannFunc);
        end
        
        function other = perturbedDirichlet(obj, dirichletIdx, delta)
            other = obj.copyModel();
            other.u0_dirichlet(dirichletIdx) = ...
                other.u0_dirichlet(dirichletIdx) + delta;
        end
        
        function other = perturbedNeumann(obj, neumannIdx, delta)
            other = obj.copyModel();
            other.en_neumann(neumannIdx) = ...
                other.en_neumann(neumannIdx) + delta;
        end
        
        function other = perturbedFreeCharge(obj, chargeIdx, delta)
            other = obj.copyModel();
            other.freeCharge(chargeIdx) = ...
                other.freeCharge(chargeIdx) + delta;
        end
        
        function other = perturbedMesh(obj, nodeIdx, dirIdx, delta)
            other = obj.copyModel();
            other.poi.tnMesh = obj.poi.tnMesh.perturbed(nodeIdx, dirIdx, delta);
            
            other.setFreeCharge(obj.chargeFunc); % need to rerun
        end
        
        function setDirichlet(obj, iNodes, nodeVals)
            obj.iDirichlet = iNodes;
            
            if isa(nodeVals, 'function_handle')
                xyFieldNodes = obj.poi.tnMesh.getNodeCoordinates();
                obj.u0_dirichlet = arrayfun(nodeVals, xyFieldNodes(iNodes,1), xyFieldNodes(iNodes, 2));
            else
                obj.u0_dirichlet = nodeVals;
            end
            
            numNodes = obj.poi.tnMesh.hFieldNodes.getNumNodes();
            obj.iCenter = setdiff(1:numNodes, obj.iDirichlet);
        end
        
        function setNeumann(obj, iNodes, nodeVals)
            obj.iNeumann = iNodes;
            
            if isa(nodeVals, 'function_handle')
                xyFieldNodes = obj.poi.tnMesh.getNodeCoordinates();
                obj.en_neumann = arrayfun(nodeVals, xyFieldNodes(iNodes,1), xyFieldNodes(iNodes, 2));
            else
                obj.en_neumann = nodeVals;
            end
            
        end
        
        function setFreeCharge(obj, freeChargeFunction)
            obj.chargeFunc = freeChargeFunction;
            [obj.freeCharge, obj.dFreeCharge_dx, obj.dFreeCharge_dy] = obj.poi.evaluateOnNodes(obj.chargeFunc);
        end
        
        function solveForwardCartesian(obj, xy0, xy1, Nxy)
            
            obj.solve();
            
            obj.interpolationOperator = obj.poi.tnMesh.getRasterInterpolationOperator(xy0, xy1, Nxy);
            
            obj.uCartesian = reshape(obj.interpolationOperator*obj.u, Nxy);
            
        end
        
        function solveCartesian(obj, xy0, xy1, Nxy)
            obj.solve();
            
            [obj.dInterpolationOperator, obj.interpolationOperator] = ...
                obj.poi.tnMesh.getRasterInterpolationOperatorSensitivity(xy0, xy1, Nxy);
            
            obj.uCartesian = reshape(obj.interpolationOperator*obj.u, Nxy);
        end
        
        function solveAdjointCartesian(obj, dF_du_cartesian)
            
            dF_du_rowVector = reshape(dF_du_cartesian, 1, []);
            
            Df_val = dF_du_rowVector * obj.interpolationOperator;
            %Df_val = objFunDerivative(obj.u);
            assert(isrow(Df_val));
            
            Df_center = Df_val(obj.iCenter);
            
            A_center = obj.A(obj.iCenter, obj.iCenter);
            B_center = obj.B(obj.iCenter, :);
            NM_neumann = obj.NM(obj.iCenter, obj.iNeumann);
            A_dirichlet = obj.A(obj.iCenter, obj.iDirichlet);
            
            % Solve for the adjoint variable
            
            v_center = A_center' \ Df_center';
            obj.v = 0*obj.u;
            obj.v(obj.iCenter) = v_center;
            
            % Matrix sensitivities
            
            
            % Get all faces that have a boundary vertex.
            % This is more faces than just those that have a boundary edge.
            [boundaryEdges, ~] = obj.poi.tnMesh.hMesh.getBoundaryEdges();
            evMat = obj.poi.tnMesh.hMesh.getEdgeVertexAdjacency();
            fvMat = obj.poi.tnMesh.hMesh.getFaceVertexAdjacency();
            feMat = fvMat * evMat';
            [iBoundaryFaces,~,~] = find(feMat(:,boundaryEdges));
            iBoundaryFaces = unique(iBoundaryFaces);
            %[boundaryFaces, = find(adjMat(:,boundaryEdges));
            
            obj.dA = obj.poi.getSystemMatrixSensitivity();
            obj.dB = obj.poi.getRhsMatrixSensitivity();
            obj.dNM = obj.poi.getNeumannMatrixSensitivity();
            
            % Sensitivity to free charge (1 x N)
            
            obj.dF_dCharge = v_center' * B_center;
            
            % Sensitivity to Dirichlet boundary value (1 x N_dirichlet)
            
            obj.dF_dDirichlet = -v_center' * A_dirichlet + Df_val(obj.iDirichlet);
            
            % Sensitivity to Neumann boundary value (1 x N_neumann)
            
            obj.dF_dNeumann = v_center' * NM_neumann;
            
            % Sensitivity to perturbing geometry nodes
            
            numGeomNodes = obj.poi.tnMesh.hGeomNodes.getNumNodes();
            
            obj.dF_dxy = zeros(numGeomNodes,2);
            
            for mm = 1:numGeomNodes
                
                wx = -obj.dA{1,mm}(obj.iCenter, obj.iCenter)*obj.u(obj.iCenter)...
                    - obj.dA{1,mm}(obj.iCenter, obj.iDirichlet)*obj.u0_dirichlet ...
                    + obj.dNM{1,mm}(obj.iCenter, obj.iNeumann)*obj.en_neumann ...
                    + obj.dB{1,mm}(obj.iCenter, :)*obj.freeCharge ...
                    + obj.B(obj.iCenter,:)*obj.dFreeCharge_dx(:,mm);
                
                wy = -obj.dA{2,mm}(obj.iCenter, obj.iCenter)*obj.u(obj.iCenter)...
                    - obj.dA{2,mm}(obj.iCenter, obj.iDirichlet)*obj.u0_dirichlet ...
                    + obj.dNM{2,mm}(obj.iCenter, obj.iNeumann)*obj.en_neumann ...
                    + obj.dB{2,mm}(obj.iCenter, :)*obj.freeCharge ...
                    + obj.B(obj.iCenter,:)*obj.dFreeCharge_dy(:,mm);
                
                dFdx = v_center'*wx + dF_du_rowVector * obj.dInterpolationOperator{1,mm} * obj.u;
                dFdy = v_center'*wy + dF_du_rowVector * obj.dInterpolationOperator{2,mm} * obj.u;
                
                obj.dF_dxy(mm,1) = dFdx;
                obj.dF_dxy(mm,2) = dFdy;
            end
            
        end
        
        
        function solve(obj)
            
            obj.A = obj.poi.getSystemMatrix();
            obj.B = obj.poi.getRhsMatrix();
            obj.NM = obj.poi.getNeumannMatrix();
            
            A_center = obj.A(obj.iCenter, obj.iCenter);
            B_center = obj.B(obj.iCenter, :);
            NM_neumann = obj.NM(obj.iCenter, obj.iNeumann);
            A_dirichlet = obj.A(obj.iCenter, obj.iDirichlet);
            
            if numel(obj.iNeumann) > 0
                u_center = A_center \ (B_center*obj.freeCharge - A_dirichlet*obj.u0_dirichlet - NM_neumann*obj.en_neumann);
            else
                u_center = A_center \ (B_center*obj.freeCharge - A_dirichlet*obj.u0_dirichlet);
            end
            
            % Make a bigger version of u including boundaries
            numNodes = obj.poi.tnMesh.hFieldNodes.getNumNodes();
            
            obj.u = zeros(numNodes,1);
            obj.u(obj.iCenter) = u_center;
            obj.u(obj.iDirichlet) = obj.u0_dirichlet;
            
            % Evaluate the functional
            %obj.F = objFun(obj.u);
        end
        
        
        
        function solveAdjoint(obj, objFunDerivative)
            
            Df_val = objFunDerivative(obj.u);
            assert(isrow(Df_val));
            
            Df_center = Df_val(obj.iCenter);
            
            A_center = obj.A(obj.iCenter, obj.iCenter);
            B_center = obj.B(obj.iCenter, :);
            NM_neumann = obj.NM(obj.iCenter, obj.iNeumann);
            A_dirichlet = obj.A(obj.iCenter, obj.iDirichlet);
            
            % Solve for the adjoint variable
            
            v_center = A_center' \ Df_center';
            obj.v = 0*obj.u;
            obj.v(obj.iCenter) = v_center;
            
            % Matrix sensitivities
            
            obj.dA = obj.poi.getSystemMatrixSensitivity();
            obj.dB = obj.poi.getRhsMatrixSensitivity();
            obj.dNM = obj.poi.getNeumannMatrixSensitivity();
            
            % Sensitivity to free charge (1 x N)
            
            obj.dF_dCharge = v_center' * B_center;
            
            % Sensitivity to Dirichlet boundary value (1 x N_dirichlet)
            
            obj.dF_dDirichlet = -v_center' * A_dirichlet + Df_val(obj.iDirichlet);
            
            % Sensitivity to Neumann boundary value (1 x N_neumann)
            
            obj.dF_dNeumann = v_center' * NM_neumann;
            
            % Sensitivity to perturbing geometry nodes
            
            numGeomNodes = obj.poi.tnMesh.hGeomNodes.getNumNodes();
            
            obj.dF_dxy = zeros(numGeomNodes,2);
            
            for mm = 1:numGeomNodes
                
                wx = -obj.dA{1,mm}(obj.iCenter, obj.iCenter)*obj.u(obj.iCenter)...
                    - obj.dA{1,mm}(obj.iCenter, obj.iDirichlet)*obj.u0_dirichlet ...
                    + obj.dNM{1,mm}(obj.iCenter, obj.iNeumann)*obj.en_neumann ...
                    + obj.dB{1,mm}(obj.iCenter, :)*obj.freeCharge ...
                    + obj.B(obj.iCenter,:)*obj.dFreeCharge_dx(:,mm);
                
                wy = -obj.dA{2,mm}(obj.iCenter, obj.iCenter)*obj.u(obj.iCenter)...
                    - obj.dA{2,mm}(obj.iCenter, obj.iDirichlet)*obj.u0_dirichlet ...
                    + obj.dNM{2,mm}(obj.iCenter, obj.iNeumann)*obj.en_neumann ...
                    + obj.dB{2,mm}(obj.iCenter, :)*obj.freeCharge ...
                    + obj.B(obj.iCenter,:)*obj.dFreeCharge_dy(:,mm);
                
                dFdx = v_center'*wx; % + Df_val * DoutI{1,mm} * obj.u;
                dFdy = v_center'*wy; % + Df_val * DoutI{2,mm} * obj.u;
                
                obj.dF_dxy(mm,1) = dFdx;
                obj.dF_dxy(mm,2) = dFdy;
            end
        end
        
        % dirichletPredicate([v1x,v1y], [v2x,v2y])
        function [iDirichlet, iNeumann, iBoundary] = classifyBoundary(obj, dirichletPredicate)
            
            iBoundary = obj.poi.tnMesh.hFieldNodes.getBoundaryNodes();
            xyBoundary = obj.poi.tnMesh.getBoundaryNodeCoordinates();
            
            flagDirichlet = 0*iBoundary;
            for ii = 1:length(iBoundary)
                xy = xyBoundary(ii,:);
                flagDirichlet(ii) = dirichletPredicate(xy(1), xy(2));
            end
            
            iDirichlet = iBoundary(flagDirichlet ~= 0);
            iNeumann = setdiff(iBoundary, iDirichlet);
        end
        
    end
    
end