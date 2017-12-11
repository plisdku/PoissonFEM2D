classdef FEMProblem < handle
    
    properties
        poi@PoissonFEM2D;
        
        iDirichlet
        iNeumann
        iCenter
        
        dirichletFunc
        neumannFunc
        chargeFunc
        
        u;
        v;
        
        F, u0_dirichlet, dF_dud, freeCharge, dFreeCharge_dx, dFreeCharge_dy, dF_df, en_neumann, dF_den, dFdv_total;
        
        dF_dCharge;
        dF_dDirichlet;
        dF_dNeumann;
        dF_dxy;
        
        A, dA;
        B, dB;
        NM, dNM;
        
    end
    
    
    methods
        
        function obj = FEMProblem(poissonFEM, dirichletPredicate)
            obj.poi = poissonFEM;
            
            numNodes = obj.poi.tnMesh.hFieldNodes.getNumNodes();
            
            [obj.iDirichlet, obj.iNeumann] = obj.classifyBoundary(dirichletPredicate);
            obj.iCenter = setdiff(1:numNodes, obj.iDirichlet);
        end
        
        function other = copyModel(obj)
            other = FEMProblem(obj.poi.copy(), @(x,y) 1);
            other.iDirichlet = obj.iDirichlet;
            other.iNeumann = obj.iNeumann;
            other.iCenter = obj.iCenter;
            other.setSources(obj.chargeFunc, obj.dirichletFunc, obj.neumannFunc);
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
            other.setSources(obj.chargeFunc, obj.dirichletFunc, obj.neumannFunc); % need to rerun
        end
        
        function setSources(obj, freeChargeFunction, dirichletFunction, neumannFunction)
            
            obj.chargeFunc = freeChargeFunction;
            obj.dirichletFunc = dirichletFunction;
            obj.neumannFunc = neumannFunction;
            
            xy = obj.poi.tnMesh.getNodeCoordinates();
            
            numNodes = size(xy, 1);
            
            obj.freeCharge = zeros(numNodes,1);
            obj.u0_dirichlet = zeros(length(obj.iDirichlet),1);
            obj.en_neumann = zeros(length(obj.iNeumann),1);
            
            [obj.freeCharge, obj.dFreeCharge_dx, obj.dFreeCharge_dy] = obj.poi.evaluateOnNodes(obj.chargeFunc);
            
            for nn = 1:length(obj.iDirichlet)
                jj = obj.iDirichlet(nn);
                obj.u0_dirichlet(nn) = obj.dirichletFunc(xy(jj,1), xy(jj,2));
            end
            
            for nn = 1:length(obj.iNeumann)
                jj = obj.iNeumann(nn);
                obj.en_neumann(nn) = obj.neumannFunc(xy(jj,1), xy(jj,2));
            end
        end
        
        function solve(obj, objFun)
            
            obj.A = obj.poi.getSystemMatrix();
            obj.B = obj.poi.getRhsMatrix();
            obj.NM = obj.poi.getNeumannMatrix();
            
            A_center = obj.A(obj.iCenter, obj.iCenter);
            B_center = obj.B(obj.iCenter, :);
            NM_center = obj.NM(obj.iCenter, obj.iCenter);
            NM_neumann = obj.NM(obj.iCenter, obj.iNeumann);
            A_dirichlet = obj.A(obj.iCenter, obj.iDirichlet);
            
            freeCharge_center = obj.freeCharge(obj.iCenter);
            
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
            obj.F = objFun(obj.u);
        end
        
        function solveAdjoint(obj, objFunDerivative)
            
            Df_val = objFunDerivative(obj.u);
            Df_center = Df_val(obj.iCenter);
            
            A_center = obj.A(obj.iCenter, obj.iCenter);
            B_center = obj.B(obj.iCenter, :);
            NM_center = obj.NM(obj.iCenter, obj.iCenter);
            NM_neumann = obj.NM(obj.iCenter, obj.iNeumann);
            A_dirichlet = obj.A(obj.iCenter, obj.iDirichlet);
            
            %freeCharge_center = obj.freeCharge(obj.iCenter);
            
            % Solve for the adjoint variable
            
            v_center = A_center' \ Df_center';
            obj.v = 0*obj.u;
            obj.v(obj.iCenter) = v_center;
            
            % Matrix sensitivities
            
            obj.dA = obj.poi.getSystemMatrixSensitivity();
            obj.dB = obj.poi.getRhsMatrixSensitivity();
            obj.dNM = obj.poi.getNeumannMatrixSensitivity();
            
            % Sensitivity to free charge
            
            obj.dF_dCharge = v_center' * B_center;
            
            % Sensitivity to Dirichlet boundary value
            
            obj.dF_dDirichlet = -v_center' * A_dirichlet;
            
            % Sensitivity to Neumann boundary value
            
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
                
                dFdx = v_center'*wx;
                dFdy = v_center'*wy;
                
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