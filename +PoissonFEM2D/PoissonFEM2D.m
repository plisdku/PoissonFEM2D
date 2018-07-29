classdef PoissonFEM2D < handle
    
    properties
        tnMesh%@TriNodalMesh
        Dr;  % differentiation matrix on basis element
        Ds;  % differentiation matrix on basis element
        Q;   % quadrature matrix on basis element
        Q1d; % quadrature matrix on basis edge
        %L;   % cell array of three edge selection matrices (why)
    end
    
    
    methods
        
        % ---- CONSTRUCTOR
        function obj = PoissonFEM2D(tnMesh)
            obj.tnMesh = tnMesh;
        end
        
        function other = copy(obj)
            other = PoissonFEM2D(obj.tnMesh.copy());
        end
        
        % ---- Elemental matrices
        
        function A = getElementPotentialMatrix(obj, iFace)
            
            [Dx,Dy] = obj.tnMesh.getFaceGradientMatrices(iFace);
            QQ = obj.tnMesh.getQuadratureMatrix(iFace);
            
            % Matrix that multiplies the potential
            A = -(Dx'*QQ*Dx + Dy'*QQ*Dy);
            
        end % getElementPotentialMatrix()
        
        function DA = getElementPotentialMatrixSensitivity(obj, iFace)
            
            [Dx,Dy] = obj.tnMesh.getFaceGradientMatrices(iFace);
            [DDx, DDy] = obj.tnMesh.getFaceGradientMatrixSensitivities(iFace);
            QQ = obj.tnMesh.getQuadratureMatrix(iFace);
            DQQ = obj.tnMesh.getQuadratureMatrixSensitivity(iFace);
            
            numGeomNodes = obj.tnMesh.hGeomNodes.basis.numNodes;
            
            DA = zeros([size(Dx), 2, numGeomNodes]);
            
            for m = 1:numGeomNodes
                for k = 1:2
                    DA(:,:,k,m) = -( ...
                        DDx(:,:,k,m)'*QQ*Dx + Dx'*DQQ(:,:,k,m)*Dx + Dx'*QQ*DDx(:,:,k,m) + ...
                        DDy(:,:,k,m)'*QQ*Dy + Dy'*DQQ(:,:,k,m)*Dy + Dy'*QQ*DDy(:,:,k,m) ...
                        );
                end
            end
            
        end % getElementPotentialMatrixSensitivity()
        
        
        function B = getElementChargeMatrix(obj, iFace)
            B = obj.tnMesh.getQuadratureMatrix(iFace);
        end % getElementChargeMatrix()
        
        function DB = getElementChargeMatrixSensitivity(obj, iFace)
            DB = obj.tnMesh.getQuadratureMatrixSensitivity(iFace);
        end % getElementChargeMatrixSensitivity()
        
        
        function C = getElementNeumannMatrix(obj, iEdge, orientation)
            C = obj.tnMesh.getQuadratureMatrix1d(iEdge, orientation);
        end % getElementNeumannMatrix()
        
        function DC = getElementNeumannMatrixSensitivity(obj, iEdge, orientation)
            DC = obj.tnMesh.getQuadratureMatrixSensitivity1d(iEdge, orientation);
        end % getElementNeumannMatrixSensitivity()
        
        
        % ---- System Matrices
        
        function NM = getNeumannMatrix(obj)
            numNodes = obj.tnMesh.hFieldNodes.getNumNodes();
            NM = sparse(numNodes, numNodes);
            
            [boundaryEdges, orientations] = obj.tnMesh.hMesh.getBoundaryEdges();
            
            numEdges = length(boundaryEdges);
            for ii = 1:numEdges
                ee = boundaryEdges(ii);
                oo = orientations(ii);
                
                iGlobal = obj.tnMesh.hFieldNodes.getEdgeNodes(ee, oo);
                
                edgeNeumannMatrix = obj.getElementNeumannMatrix(ee,oo);
                
                NM(iGlobal,iGlobal) = NM(iGlobal,iGlobal) + edgeNeumannMatrix;
                
            end
            
        end % getNeumannMatrix
        
        function DNM = getNeumannMatrixSensitivity(obj, iEdges)
            numFieldNodes = obj.tnMesh.hFieldNodes.getNumNodes();
            numGeomNodes = obj.tnMesh.hGeomNodes.getNumNodes();
            
            DNM = cell(2,numGeomNodes);
            for nn = 1:numel(DNM)
                DNM{nn} = sparse(numFieldNodes,numFieldNodes);
            end
            
            [boundaryEdges, orientations] = obj.tnMesh.hMesh.getBoundaryEdges();
            
            if nargin < 2
                numEdges = length(boundaryEdges);
                iEdges = 1:numEdges;
            end
                
            
            for ii = 1:reshape(iEdges, 1, [])
                ee = boundaryEdges(ii);
                oo = orientations(ii);
                
                iGeomGlobal = obj.tnMesh.hGeomNodes.getEdgeNodes(ee,oo);
                iFieldGlobal = obj.tnMesh.hFieldNodes.getEdgeNodes(ee,oo);
                
                DNM_edge = obj.getElementNeumannMatrixSensitivity(ee,oo);
                
                for mm = 1:length(iGeomGlobal)
                    for kk = 1:2
                        DNM{kk,iGeomGlobal(mm)}(iFieldGlobal,iFieldGlobal) = ...
                            DNM{kk,iGeomGlobal(mm)}(iFieldGlobal,iFieldGlobal) + DNM_edge(:,:,kk,mm);
                    end
                end
            end
            
        end % getNeumannMatrixSensitivity
        
        
        
        function M = getRhsMatrix(obj)
            numNodes = obj.tnMesh.hFieldNodes.getNumNodes();
           % M = sparse(numNodes, numNodes);
            
            numFaces = obj.tnMesh.hMesh.getNumFaces();
            
            iGlobals1 = {};%[];
            iGlobals2 = {};%[];
            faceMatrices = {};%[];
            
            for ff = 1:numFaces
                faceMatrix = obj.getElementChargeMatrix(ff);
                
                iGlobal = obj.tnMesh.hFieldNodes.getFaceNodes(ff);
                
                [repGlobals1, repGlobals2] = ndgrid(iGlobal, iGlobal); 
                iGlobals1{ff} = repGlobals1(:);%[iGlobals1; repGlobals1(:)];
                iGlobals2{ff} = repGlobals2(:);%[iGlobals2; repGlobals2(:)];
                faceMatrices{ff} = faceMatrix(:);%[faceMatrices; faceMatrix(:)];
                %iGlobals2 = 
               % M(iGlobal, iGlobal) = M(iGlobal,iGlobal) + faceMatrix;
            end
            
            
            M = sparse(cell2mat(iGlobals1)',cell2mat(iGlobals2)',...
                cell2mat(faceMatrices)',numNodes,numNodes);%sparse(iGlobals1, iGlobals2, faceMatrices, numNodes, numNodes);
            
            
        end % getRhsMatrix
        
        function DM = getRhsMatrixSensitivity(obj, iFaces)
            numFieldNodes = obj.tnMesh.hFieldNodes.getNumNodes();
            numGeomNodes = obj.tnMesh.hFieldNodes.getNumNodes();
            DM = cell(2, numGeomNodes);
            for nn = 1:numel(DM)
                DM{nn} = sparse(numFieldNodes,numFieldNodes);
            end
            
            
            if nargin < 2
                numFaces = obj.tnMesh.hMesh.getNumFaces();
                iFaces = 1:numFaces;
            end
            
            for ff = reshape(iFaces, 1, [])
                DM_face = obj.getElementChargeMatrixSensitivity(ff);
                
                iFieldGlobal = obj.tnMesh.hFieldNodes.getFaceNodes(ff);
                iGeomGlobal = obj.tnMesh.hGeomNodes.getFaceNodes(ff);
                
                for mm = 1:length(iGeomGlobal)
                    for kk = 1:2
                        DM{kk,iGeomGlobal(mm)}(iFieldGlobal,iFieldGlobal) = ...
                            DM{kk,iGeomGlobal(mm)}(iFieldGlobal,iFieldGlobal) + DM_face(:,:,kk,mm);
                    end
                end
            end
        end % getRhsMatrixSensitivity
        
        
        
        function M = getSystemMatrix(obj)
            
            numNodes = obj.tnMesh.hFieldNodes.getNumNodes();
            
           % M = sparse(numNodes, numNodes);
           
            iGlobals1 = {};%[];
            iGlobals2 = {};%[];
            facePotentials = {};%[];
            
            numFaces = obj.tnMesh.hMesh.getNumFaces();
            
            for ff = 1:numFaces
                
                facePotentialM = obj.getElementPotentialMatrix(ff);
                
                iGlobal = obj.tnMesh.hFieldNodes.getFaceNodes(ff);
                [repGlobals1, repGlobals2] = ndgrid(iGlobal, iGlobal); 
                iGlobals1{ff} = repGlobals1(:);%[iGlobals1; repGlobals1(:)];
                iGlobals2{ff} = repGlobals2(:);%[iGlobals2; repGlobals2(:)];
                facePotentials{ff} = facePotentialM(:);%[facePotentials; facePotentialM(:)];
                %M(iGlobal,iGlobal) = M(iGlobal, iGlobal) + facePotentialM;
                
            end
            
            M = sparse(cell2mat(iGlobals1)',cell2mat(iGlobals2)',...
                cell2mat(facePotentials)',numNodes,numNodes);%sparse(iGlobals1, iGlobals2, facePotentials, numNodes, numNodes);
            
        end % getSystemMatrix
        
        
        function DM = getSystemMatrixSensitivity(obj, iFaces)
            numFieldNodes = obj.tnMesh.hFieldNodes.getNumNodes();
            numGeomNodes = obj.tnMesh.hGeomNodes.getNumNodes();
             
           
%             DM = cell(2, numGeomNodes);
%             for nn = 1:numel(DM)
%                 DM{nn} = sparse(numFieldNodes, numFieldNodes);
%             end

            DM{2,numGeomNodes} = sparse(numFieldNodes, numFieldNodes);
            DM(:) = {sparse(numFieldNodes, numFieldNodes)};
            
            numFaces = obj.tnMesh.hMesh.getNumFaces();
            if nargin < 2
                iFaces = 1:numFaces;
            end
            
            for ff = reshape(iFaces, 1, [])
                DM_face = obj.getElementPotentialMatrixSensitivity(ff);
                iFieldGlobal = obj.tnMesh.hFieldNodes.getFaceNodes(ff);
                iGeomGlobal = obj.tnMesh.hGeomNodes.getFaceNodes(ff);
                
                
                for mm = 1:length(iGeomGlobal)
                    for kk = 1:2
                        DM{kk,iGeomGlobal(mm)}(iFieldGlobal,iFieldGlobal) = ...
                            DM{kk,iGeomGlobal(mm)}(iFieldGlobal,iFieldGlobal) + DM_face(:,:,kk,mm);
                    end
                end
            end
          
            
%             tic
%             numFaces = obj.tnMesh.hMesh.getNumFaces();
%             if nargin < 2
%                 iFaces = 1:numFaces;
%             end
%             iFieldGlobals = [];
%             iGeomGlobals = [];
%             DM_faces{2, length(iGeomGlobal)} = [];
%             disp('done initialisation')
%             tic
%             toc
%             for ff = reshape(iFaces, 1, [])
%                 
%                 [iFeldGlobalrep, iGeomGlobalrep] = ndgrid(obj.tnMesh.hFieldNodes.getFaceNodes(ff), obj.tnMesh.hGeomNodes.getFaceNodes(ff));
% 
%                 %iFieldGlobals = [iFieldGlobals; iFeldGlobalrep(:)];
%                 %iGeomGlobals = [iGeomGlobals; iGeomGlobalrep(:)];
%                 keyboard
%                 DM_face = obj.getElementPotentialMatrixSensitivity(ff);
%                 %toc
%                 %disp('done getting indexes (ff)')
%                 %disp(ff)
%                 for mm = 1:length(iGeomGlobal)
%                     %DM_faces{1,mm} = [DM_faces{1,mm}; reshape(DM_face(:,:,1,mm),[],1)];
%                     %DM_faces{2,mm} = [DM_faces{2,mm}; reshape(DM_face(:,:,2,mm),[],1)];
%                     %iFieldGlobals = [iFieldGlobals; iFeldGlobalrep(:)];
%                     %iGeomGlobals = [iGeomGlobals; iGeomGlobalrep(:)];
%                     
%                     %toc
%                     %disp('done creating one round of values')
%                     %disp(mm)
%                 end 
% 
%             end
%             toc
%             DM{2,numGeomNodes} = sparse(numFieldNodes, numFieldNodes);
%             toc
%             disp('done init DM')
%             for mm = 1:length(iGeomGlobal)
%                 DM{1,mm} = sparse(iFieldGlobals, iGeomGlobals, DM_faces{1,mm}, numFieldNodes, numFieldNodes);
%                 DM{2,mm} = sparse(iFieldGlobals, iGeomGlobals, DM_faces{2,mm}, numFieldNodes, numFieldNodes);
%                 toc
%                 disp('done replacing one matrix')
%                 disp(mm)
%             end 
%             
%             toc
%             
            
            
            
            
        end % getSystemMatrixSensitivity
        
        
        
        
        % ---- Function evaluation
        
        function [f, df_dxg, df_dyg] = evaluateOnNodes(obj,func)
            % [f, df_dxg, df_dyg] = evaluateOnNodes(func)
            %
            % Evaluate the given function on the field nodes.
            % Calculate the sensitivity of the evaluated values to changing
            % the geometry nodes.
            
            xy = obj.tnMesh.getNodeCoordinates();
            
            numNodes = size(xy,1);
            f = zeros(numNodes,1);
            
            for nn = 1:numNodes
                f(nn) = func(xy(nn,1), xy(nn,2));
            end
            
            if nargout == 1
                return
            end
            
            dfdxy = zeros(numNodes,2);
            delta = 1e-8;
            for nn = 1:numNodes
                dfdxy(nn,1) = (func(xy(nn,1)+delta, xy(nn,2)) - func(xy(nn,1)-delta, xy(nn,2)))/(2*delta);
                dfdxy(nn,2) = (func(xy(nn,1), xy(nn,2)+delta) - func(xy(nn,1), xy(nn,2)-delta))/(2*delta);
            end
            
            % Sensitivity of the function on field nodes to changing
            % the geometry nodes
            
            %numGeomNodes = obj.tnMesh.hGeomNodes.getNumNodes();
            %df_dxg = sparse(numNodes, numGeomNodes);
            %df_dyg = sparse(numNodes, numGeomNodes);
            
            % Some nodes belong to one element, some to two, some to many.
            % This makes it a little harder for me to get the field node
            % position sensitivities out.
            
            dxdx = obj.tnMesh.getNodeCoordinateSensitivities();
            
            df_dxg = spdiags(dfdxy(:,1),0,size(dfdxy(:,1),1),size(dfdxy(:,1),1)) * dxdx;

            df_dyg = spdiags(dfdxy(:,2),0,size(dfdxy(:,2),1),size(dfdxy(:,2),1)) * dxdx;
            
            
            
        end % evaluateOnNodess
        
        function [f, dfdv] = oldEvaluateOnNodes(obj, func)
            % [f, dfdv] = evaluateOnNodes(func)
            %
            % Evaluate the given function on the field nodes.
            % Calculate the sensitivity of the evaluated values to changing
            % the mesh vertices.
            
            % I guess it will work like... uhhhhhh...
            
            xy = obj.tnMesh.getNodeCoordinates();
            dxy_dv = obj.tnMesh.getNodeCoordinateSensitivities();
            
            % Calculate f and its gradient components
            f = func(xy(:,1), xy(:,2));
            
            delta = 1e-6;
            dfdx = (func(xy(:,1)+delta, xy(:,2)) - func(xy(:,1)-delta, xy(:,2)))/2/delta;
            dfdy = (func(xy(:,1), xy(:,2)+delta) - func(xy(:,1), xy(:,2)-delta))/2/delta;
            
            % Somehow we need to get dfdv now...
            
            dfdv = cell(obj.tnMesh.hMesh.getNumVertices(), 2);
            for nn = 1:numel(dfdv)
                dfdv{nn} = sparse(obj.tnMesh.hNodes.getNumNodes(), 1);
            end
            
            % Chain rule to get dfdv
            
            for iVert = 1:size(dxy_dv,1)
                for iXY = 1:2
                    dfdv{iVert,iXY} = dfdx .* dxy_dv{iVert,iXY}(:,1) + ...
                        dfdy .* dxy_dv{iVert,iXY}(:,2);
                end
            end
            
        end
        
        % ---- Functionals
        
        function [I, dIdJ, dIdu] = elementIntegralFunctional(obj, f, dfdu, jacobian)
            % [I, dIdJ, dIdu] = elementIntegralFunctional(f, dfdu, jacobian)
            %
            % Evaluate integral of f(u) in a single element.
            
            assert(iscolumn(f), 'f must be a column vector');
            assert(iscolumn(dfdu), 'dfdu must be a column vector, pointwise sensitivities df/du');
            
            [Qv, dQvdJ] = obj.elementQuadratureVector(jacobian);
            % dQvdJ indices are (q1, jacobian1, jacobian2)
            
            I = Qv'*f;
            dIdJ = multiplyTensors.txt(dQvdJ, 3, f, 1, 1, 1);
            
            dIdu = Qv'.*dfdu';
        end
        
        function [F, dFdv, dFdu] = surfaceIntegralFunctional(obj, integrandFunction, elementIndices, u)
            % [F, dFdp, dFdu] = surfaceIntegralFunctional(integrandFunction, elementIndices, u)
            %
            % Evaluate the integral of a function of u over selected
            % elements.
            %
            % The sensitivity dFdp is not calculated yet because I don't
            % know if I want to jump to d/dvertex or not.
            
            
            % TODO: implement dFdv
                
            numNodes = obj.tnMesh.hNodes.getNumNodes();
            dFdu = sparse(1, numNodes);
            dFdv = [];
            
            numIntegrandFaces = length(elementIndices);
            
            F = 0;
            for ii = 1:numIntegrandFaces
                ff = elementIndices(ii);
                
                iGlobal = obj.tnMesh.hNodes.getFaceNodes(ff);
                [func, dfunc_du] = integrandFunction(u(iGlobal));
                jacobian = obj.tnMesh.getLinearJacobian(ff);
                dJdv = obj.tnMesh.getLinearJacobianSensitivity(ff);
                
                [I, dIdJ, dIdu] = obj.elementIntegralFunctional(func, dfunc_du, jacobian);
                dIdv = multiplyTensors.txt(dIdJ, 2, dJdv, 4, 1:2, 1:2);
                
                F = F + I;
                dFdu(iGlobal) = dFdu(iGlobal) + dIdu;
                %dFdv(iGlobal) = dFdv(iGlobal) + dIdv; % this seems to not
                %be implemented yet
                
            end
            
        end
        
        function [F, dFdv, dFdu] = pointEvaluationFunctional(obj, pointFunction, xy, u)
            % [F, dFdp, dFdu] = pointEvaluationFunctional(pointFunction, xy, u)
            %
            % Calculate f(u) and df/du at a single point.
            
            assert(length(xy) == 2, 'Functional must be evaluated at a single point');
            
            interpMatrix = obj.tnMesh.getInterpolationOperator(xy(1), xy(2));
            
            uEval = interpMatrix*u;
            
            [F, dfunc_du] = pointFunction(uEval);
            dFdu = interpMatrix * dfunc_du;
            
            assert(isrow(dFdu));
            
            % now the sensitivity.
            
            dIdv = obj.tnMesh.getInterpolationOperatorSensitivity(xy(1), xy(2));
            dFdv = zeros(size(dIdv));
            
            for nn = 1:numel(dFdv)
                dFdv(nn) = dIdv{nn}*u;
            end
            
        end
        
        
    end % methods
    
    
    
end