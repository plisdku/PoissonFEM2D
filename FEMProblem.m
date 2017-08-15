classdef FEMProblem < handle
    
    properties
        fem@PoissonFEM2D;
        
        iNodesDirichlet
        iNodesNeumann
        iNodesCenter
        
        u;
        v;
        
        F, u0_dirichlet, dF_dud, freeCharge, dFreeCharge_dv, dF_df, en_neumann, dF_den, dFdv_total;
        
        A, dA_dv;
        B, dB_dv;
        NM, dNM_dv;
        
    end
    
    
    methods
        
        function obj = FEMProblem(N, faces, vertices, dirichletPredicate)
            meshNodes = TriNodalMesh(N, faces, vertices);
            obj.fem = PoissonFEM2D(meshNodes);
            
            numNodes = obj.fem.meshNodes.getNumNodes();
            [obj.iNodesDirichlet, obj.iNodesNeumann] = obj.classifyEdges(dirichletPredicate);
            obj.iNodesCenter = setdiff(1:numNodes, obj.iNodesDirichlet);
        end
        
        function setSources(obj, freeChargeFunction, dirichletFunction, neumannFunction)
            
            xyNodes = obj.fem.meshNodes.getNodeCoordinates();
            
            % Free charge
            [obj.freeCharge, obj.dFreeCharge_dv] = obj.fem.evaluateOnNodes(freeChargeFunction);
            
            % Dirichlet boundary condition
            obj.u0_dirichlet = dirichletFunction(xyNodes(obj.iNodesDirichlet,:));
            
            % Neumann boundary condition
            obj.en_neumann = neumannFunction(xyNodes(obj.iNodesNeumann,:));
            
            
        end
        
        function solve(obj, evalPt)

            numNodes = obj.fem.meshNodes.getNumNodes();
            
            [A, dA_dv] = obj.fem.systemMatrix();
            [B, dB_dv] = obj.fem.rhsMatrix();
            [NM, dNM_dv] = obj.fem.neumannMatrix();
            
            obj.A = A;
            obj.dA_dv = dA_dv;
            obj.B = B;
            obj.dB_dv = dB_dv;
            obj.NM = NM;
            obj.dNM_dv = dNM_dv;
            
            
            
            % Partition the system matrices to separate out Dirichlet nodes
            iDirichlet = obj.iNodesDirichlet;
            iNeumann = obj.iNodesNeumann;
            iCenter = obj.iNodesCenter;
            
            A_center = A(iCenter, iCenter);
            B_center = B(iCenter, iCenter);
            NM_center = NM(iCenter, iCenter);
            NM_neumann = NM(iCenter, iNeumann);
            A_dirichlet = A(iCenter, iDirichlet);
            
            freeCharge_center = obj.freeCharge(iCenter);
            
            if numel(iNeumann) > 0
                u_center = A_center \ (B_center*freeCharge_center - A_dirichlet*obj.u0_dirichlet - NM_neumann*obj.en_neumann);
            else
                u_center = A_center \ (B_center*freeCharge_center - A_dirichlet*obj.u0_dirichlet);
            end
            
            % Make a bigger version of u including boundaries
            u = zeros(numNodes,1);
            u(iCenter) = u_center;
            u(iDirichlet) = obj.u0_dirichlet;
            obj.u = u;
            
            % Evaluate the functional
            [F, dFdv, dFdu] = obj.fem.pointEvaluationFunctional(@multiplyByOne, evalPt, u);
            obj.F = F;
            
            % Solve the adjoint system
            
            dFdu_center = dFdu(iCenter);
            
            v_center = A_center' \ dFdu_center';
            v = 0*u;
            v(iCenter) = v_center;
            obj.v = v;
            
            % Get some sensitivities
            
            % Sensitivity with respect to Dirichlet boundary values
            dF_dud = -v_center'*A_dirichlet;
            dF_dud_big = 0*u;
            dF_dud_big(iDirichlet) = dF_dud;
            
            %obj.u0_dirichlet = u0_dirichlet;
            obj.dF_dud = dF_dud;
            
            % Sensitivity with respect to free charge
            dF_df = v_center'*B_center;
            dF_df_big = 0*u;
            dF_df_big(iCenter) = dF_df;
            
            %obj.freeCharge = freeCharge;
            obj.dF_df = dF_df_big;
            
            % Sensitivity with respect to normal E on Neumann boundary nodes
            dF_den = -v_center'*NM_neumann;
            dF_den_big = 0*u;
            dF_den_big(iNeumann) = dF_den;
            
            %obj.en_neumann = en_neumann;
            obj.dF_den = dF_den;
            
            % Perturb the mesh
            numVerts = obj.fem.meshNodes.getNumVertices();

            dFdv_total = zeros(numVerts,2);

            % u_center = A_center \ (B_center*freeCharge_center - A_dirichlet*obj.u0_dirichlet - NM_neumann*obj.en_neumann);
            % A*u = B*f - Ad*ud - NM*en
            % 
            % A*du = dB*f + B*df - dAd*ud - Ad*dud - dNM*en - NM*den - dA*u
            % 
            for vv = 1:numVerts

                wx = -dA_dv{vv,1}(iCenter, iCenter)*u_center ...
                    - dA_dv{vv,1}(iCenter, iDirichlet)*obj.u0_dirichlet ...
                    + dNM_dv{vv,1}(iCenter, iNeumann)*obj.en_neumann ...
                    + dB_dv{vv,1}(iCenter, iCenter)*freeCharge_center ...
                    + B_center*obj.dFreeCharge_dv{vv,1}(iCenter);
                
                wy = -dA_dv{vv,2}(iCenter, iCenter)*u_center ...
                    - dA_dv{vv,2}(iCenter, iDirichlet)*obj.u0_dirichlet ...
                    + dNM_dv{vv,2}(iCenter, iNeumann)*obj.en_neumann ...
                    + dB_dv{vv,2}(iCenter, iCenter)*freeCharge_center...
                    + B_center*obj.dFreeCharge_dv{vv,2}(iCenter);
                
                dFdvx = v_center'*wx;
                dFdvy = v_center'*wy;

                dFdv_total(vv,1) = dFdvx;
                dFdv_total(vv,2) = dFdvy;
            end

            dFdv_total = dFdv_total + dFdv;
            obj.dFdv_total = dFdv_total;
            
        end
        
        
        % dirichletPredicate([v1x,v1y], [v2x,v2y])
        function [iNodesDirichlet, iNodesNeumann, iNodesBoundary] = classifyEdges(obj, dirichletPredicate)
            
            iBoundaryEdges = obj.fem.meshNodes.getBoundaryEdges();
            
            flagDirichlet = false(length(iBoundaryEdges),1);
            for ii = 1:length(iBoundaryEdges)
                iEdge = iBoundaryEdges(ii);
                verts = obj.fem.meshNodes.getVertexNodeCoordinates(...
                    obj.fem.meshNodes.getEdgeVertices(iEdge));
                
                flagDirichlet(ii) = dirichletPredicate(verts(1,:), verts(2,:));
            end
            
            iDirichletEdges = iBoundaryEdges(flagDirichlet);
            
            iNodesBoundary = obj.fem.meshNodes.getBoundaryNodes();
            iNodesDirichlet = obj.fem.meshNodes.getEdgeNodes(iDirichletEdges);
            iNodesNeumann = setdiff(iNodesBoundary, iNodesDirichlet);
        end
        
    end
    
end