classdef PoissonFEM2D < handle
    
    properties
        meshNodes@TriNodalMesh
        Dr;  % differentiation matrix on basis element
        Ds;  % differentiation matrix on basis element
        Q;   % quadrature matrix on basis element
        Q1d; % quadrature matrix on basis edge
        L;   % cell array of three lift matrices
    end
    
    
    methods
        
        % ---- CONSTRUCTOR
        function obj = PoissonFEM2D(meshNodes)
            
            obj.meshNodes = meshNodes;
            
            rs = obj.meshNodes.basis.r;
            ss = obj.meshNodes.basis.s;
            
            [obj.Dr, obj.Ds] = support2d.gradients(obj.meshNodes.N, rs, ss);
            obj.Q = support2d.quadratureKernel(obj.meshNodes.N, rs, ss);
            
            obj.Q1d = support.quadratureKernel(obj.meshNodes.N, obj.meshNodes.basis1d.r);

            % make the lift matrices
            [iCorners, iEdgeCenters, ~] = support2d.classifyNodes(obj.meshNodes.N);
            
            iEdge1 = [iCorners(1), iEdgeCenters{1}, iCorners(2)];
            iEdge2 = [iCorners(2), iEdgeCenters{2}, iCorners(3)];
            iEdge3 = [iCorners(3), iEdgeCenters{3}, iCorners(1)];
            nEdge = length(iEdge1);
            nNodes = size(obj.Q,1);
            unos = ones(size(iEdge1));
            obj.L{1} = sparse(1:nEdge, iEdge1, unos, nEdge, nNodes);
            obj.L{2} = sparse(1:nEdge, iEdge2, unos, nEdge, nNodes);
            obj.L{3} = sparse(1:nEdge, iEdge3, unos, nEdge, nNodes);
        end
        
        % ---- Helper matrices
        
        function [Dx, Dy, dDxdJ, dDydJ] = elementGradientMatrix(obj, jacobian)
            invJac = inv(jacobian);
            
            Dx = obj.Dr*invJac(1,1) + obj.Ds*invJac(2,1);
            Dy = obj.Dr*invJac(1,2) + obj.Ds*invJac(2,2);
            
            dDxdJ = zeros([size(Dx), 2, 2]);
            dDydJ = zeros([size(Dy), 2, 2]);
            
            for ii = 1:2
                for jj = 1:2
                    dDxdJ(:,:,ii,jj) = -invJac(1,ii)*invJac(jj,1)*obj.Dr - invJac(2,ii)*invJac(jj,1)*obj.Ds;
                    dDydJ(:,:,ii,jj) = -invJac(1,ii)*invJac(jj,2)*obj.Dr - invJac(2,ii)*invJac(jj,2)*obj.Ds;
                end
            end
            
%             dDxdJ = cell(2,2);
%             dDydJ = cell(2,2);
%             
%             for ii = 1:2
%                 for jj = 1:2
%                     dDxdJ{ii,jj} = -invJac(1,ii)*invJac(jj,1)*obj.Dr - invJac(2,ii)*invJac(jj,1)*obj.Ds;
%                     dDydJ{ii,jj} = -invJac(1,ii)*invJac(jj,2)*obj.Dr - invJac(2,ii)*invJac(jj,2)*obj.Ds;
%                 end
%             end
        end % elementGradientMatrix
        
        
        function [Q, dQdJ] = elementQuadratureMatrix(obj, jacobian)
            % [Q, dQdJ] = elementQuadratureMatrix(obj, jacobian)
            %
            % Sensitivities to jacobian(ii,jj) are indexed dQdJ(iNode1, iNode2, ii, jj).
            %
            
            invJac = inv(jacobian);
            detJac = det(jacobian);
            Q = obj.Q*detJac;
            
            dDetJ_dJ = detJac*transpose(invJac); % this is all 4 sensitivities in a matrix
            
            dQdJ = zeros([size(Q), 2, 2]);
            for ii = 1:2
                for jj = 1:2
                    dQdJ(:,:,ii,jj) = obj.Q*dDetJ_dJ(ii,jj);
                end
            end
        end % elementQuadratureMatrix
        
        
        function [Qv, dQvdJ] = elementQuadratureVector(obj, jacobian)
            % [Qv, dQvdJ] = elementQuadratureVector(obj, jacobian)
            %
            % Qv is a column vector indexed Qv(node)
            % dQvdJ is indexed dQvdJ(node, jacobian1, jacobian2)
            
            [Q, dQdJ] = obj.elementQuadratureMatrix(jacobian);
            
            sz_dQdJ = size(dQdJ);
            
            Qv = reshape(sum(Q, 1), [], 1);
            dQvdJ = reshape(sum(dQdJ, 1), sz_dQdJ(2:4));
            
        end
            
        
        % ---- Elemental matrices
        function [A, dAdJ] = elementPotentialMatrix(obj, jacobian)
            % [A, dA_dJij] = elementPotentialMatrix(obj, jacobian)
            %
            % Sensitivities are indexed dA_dJij{ii}{jj}.
            %
            
            [Dx, Dy, dDxdJ, dDydJ] = obj.elementGradientMatrix(jacobian);
            [QQ, dQdJ] = obj.elementQuadratureMatrix(jacobian);
            
            % Matrix that multiplies the potential
            A = - (Dx'*QQ*Dx + Dy'*QQ*Dy);
            
            if nargout == 1
                return
            end
            
            % Its sensitivities with respect to Jacobian elements
           
            
            dAdJ = zeros([size(A), 2,2]);
            for ii = 1:2
                for jj = 1:2
                    
                    dAdJ(:,:,ii,jj) = -(dDxdJ(:,:,ii,jj)'*QQ*Dx + Dx'*dQdJ(:,:,ii,jj)*Dx + Dx'*QQ*dDxdJ(:,:,ii,jj)) ...
                        - (dDydJ(:,:,ii,jj)'*QQ*Dy + Dy'*dQdJ(:,:,ii,jj)*Dy + Dy'*QQ*dDydJ(:,:,ii,jj));
                end
            end
            
%             dAdJ = cell(2,2);
%             for ii = 1:2
%                 for jj = 1:2
%                     
%                     dAdJ{ii,jj} = -(dDxdJ{ii,jj}'*QQ*Dx + Dx'*dQdJ{ii,jj}*Dx + Dx'*QQ*dDxdJ{ii,jj}) ...
%                         - (dDydJ{ii,jj}'*QQ*Dy + Dy'*dQdJ{ii,jj}*Dy + Dy'*QQ*dDydJ{ii,jj});
%                 end
%             end
        end % elementPotentialMatrix()
        
        function [B, dBdJ] = elementChargeMatrix(obj, jacobian)
            % [B, dBdJ_ij] = elementChargeMatrix(obj, jacobian)
            %
            % Sensitivities are indexed dB_dJ_ij{ii}{jj}.
            %
            
            [B, dBdJ] = obj.elementQuadratureMatrix(jacobian);
            
        end % elementChargeMatrix()
        
        
        function [C, dCdJ] = elementNeumannMatrix(obj, jacobian1d, nodeFactors)
            % [C, dCdJ] = elementNeumannMatrix(jacobian1d)
            % [C, dCdJ] = elementNeumannMatrix(jacobian1d, nodeFactors)
            %
            % Calculate Neumann matrix for a single edge, with or without
            % factors.
            
            detJ = sqrt(det(jacobian1d'*jacobian1d));
            
            if nargin < 3
                C = obj.Q1d * detJ;
            else
                C = obj.Q1d * detJ * diag(nodeFactors);
            end
            
            %dDetJ_dJ = detJ*transpose(invJac); % this is all 4 sensitivities in a matrix
            dDetJ_dJ = detJ*jacobian1d*inv(jacobian1d'*jacobian1d);
            
            dCdJ = zeros([size(C), 2]);
            
            for ii = 1:2
                if nargin < 3
                    dCdJ(:,:,ii) = obj.Q1d*dDetJ_dJ(ii);
                else
                    dCdJ(:,:,ii) = obj.Q1d*dDetJ_dJ(ii) * diag(nodeFactors);
                end
            end
            
        end % elementNeumannMatrix
        
        % ---- System Matrices
        
        function [NM, DNM_dv] = neumannMatrix(obj, nodeFactors)
            
            % The Neumann matrix is sum( [get edge nodes] * Q * detJ ).
            % Its raw dimensions are nNodes x nNodes and on the outside
            % we can bake it down to the size of only the Neumann boundary
            % nodes.
            
            numNodes = obj.meshNodes.getNumNodes();
            numVertices = obj.meshNodes.getNumVertices();
            
            NM = sparse(numNodes, numNodes);
            DNM_dv = cell(numVertices, 2);
            for nn = 1:numel(DNM_dv)
                DNM_dv{nn} = sparse(numNodes, numNodes);
            end
            
            [boundaryEdges, orientations] = obj.meshNodes.getBoundaryEdges();
            
            numEdges = length(boundaryEdges);
            
            for ii = 1:numEdges
                ee = boundaryEdges(ii);
                oo = orientations(ii);
                
                dJdv = obj.meshNodes.getLinearJacobianSensitivity1d(ee, oo);
                
                iVertices = obj.meshNodes.getEdgeVertices(ee,oo);
                iGlobal = obj.meshNodes.getEdgeNodes(ee, oo);
                
                jac1d = obj.meshNodes.getLinearJacobian1d(ee, oo);
                
                if nargin < 2
                    [neuMat, dNdJ] = obj.elementNeumannMatrix(jac1d);
                else
                    [neuMat, dNdJ] = obj.elementNeumannMatrix(jac1d, nodeFactors(iGlobal));
                end
                dNdv_local = multiplyTensors.txt(dNdJ, 3, dJdv, 3, 3, 1);
                
                NM(iGlobal, iGlobal) = NM(iGlobal, iGlobal) + neuMat;
                
                for iVert = 1:2
                    iVertGlobal = iVertices(iVert);
                    for iXY = 1:2
                        DNM_dv{iVertGlobal, iXY}(iGlobal, iGlobal) = ...
                            DNM_dv{iVertGlobal, iXY}(iGlobal, iGlobal) + ...
                            dNdv_local(:,:,iVert,iXY);
                    end
                end
                
            end
            
        end
        
        function [rhsMatrix, DrhsMatrix_dv] = rhsMatrix(obj)
            
            numNodes = obj.meshNodes.getNumNodes();
            numVertices = obj.meshNodes.getNumVertices();
            rhsMatrix = sparse(numNodes, numNodes);
            
            DrhsMatrix_dv = cell(numVertices,2);
            
            for nn = 1:numel(DrhsMatrix_dv)
                DrhsMatrix_dv{nn} = sparse(numNodes, numNodes);
            end
            
            
            % Fill in the matrix!
            numFaces = obj.meshNodes.getNumFaces();

            for ff = 1:numFaces
                jacobian = obj.meshNodes.getLinearJacobian(ff);
                dJdv = obj.meshNodes.getLinearJacobianSensitivity(ff);
                
                [M2, dM2dJ] = obj.elementChargeMatrix(jacobian);
                dM2dv = multiplyTensors.txt(dM2dJ, 4, dJdv, 4, 3:4, 1:2);
                
                iVertices = obj.meshNodes.getFaceVertices(ff);
                iGlobal = obj.meshNodes.getFaceNodes(ff);

                rhsMatrix(iGlobal, iGlobal) = rhsMatrix(iGlobal, iGlobal) + M2;
                
                for iVert = 1:3
                    iVertGlobal = iVertices(iVert);
                    for iXY = 1:2
                        DrhsMatrix_dv{iVertGlobal, iXY}(iGlobal, iGlobal) = ...
                            DrhsMatrix_dv{iVertGlobal, iXY}(iGlobal, iGlobal) + ...
                            dM2dv(:,:,iVert,iXY);
                    end
                end
            end
        end
        
        function [systemMatrix, DsystemMatrix_dv] = systemMatrix(obj)
            % [A, B, dAdJ, dBdJ] = systemMatrix(obj)
            %
            % The jacobian override is for sensitivity testing.
            
            numNodes = obj.meshNodes.getNumNodes();
            numVertices = obj.meshNodes.getNumVertices();
            
            systemMatrix = sparse(numNodes, numNodes);

            % We will accumulate the system matrix sensitivities with respect to
            % elemental Jacobians, first.
            %
            % DsystemMatrix_dJ{ii,jj} is the matrix of sensitivity to J(ii,jj).

            DsystemMatrix_dv = cell(numVertices,2);
            for nn = 1:numel(DsystemMatrix_dv)
                DsystemMatrix_dv{nn} = sparse(numNodes, numNodes);
            end
            
            % Fill in the matrices!
            numFaces = obj.meshNodes.getNumFaces();

            for ff = 1:numFaces
                
                jacobian = obj.meshNodes.getLinearJacobian(ff);
                dJdv = obj.meshNodes.getLinearJacobianSensitivity(ff);
                
                [M1, dM1dJ] = obj.elementPotentialMatrix(jacobian);
                dM1dv = multiplyTensors.txt(dM1dJ, 4, dJdv, 4, 3:4, 1:2);
                
                iVertices = obj.meshNodes.getFaceVertices(ff);
                iGlobal = obj.meshNodes.getFaceNodes(ff);
                
                systemMatrix(iGlobal,iGlobal) = systemMatrix(iGlobal, iGlobal) + M1;
                
                for iVert = 1:3
                    for iXY = 1:2
                        iVertGlobal = iVertices(iVert);
                        DsystemMatrix_dv{iVertGlobal, iXY}(iGlobal, iGlobal) = ...
                            DsystemMatrix_dv{iVertGlobal, iXY}(iGlobal, iGlobal) + ...
                            dM1dv(:,:,iVert,iXY);
                    end
                end

            end

            
        end % systemMatrix()
        
        
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
            
            numNodes = obj.meshNodes.getNumNodes();
            dFdu = sparse(1, numNodes);
            dFdv = [];
            
            numIntegrandFaces = length(elementIndices);
            
            F = 0;
            for ii = 1:numIntegrandFaces
                ff = elementIndices(ii);
                
                iGlobal = obj.meshNodes.getFaceNodes(ff);
                [func, dfunc_du] = integrandFunction(u(iGlobal));
                jacobian = obj.meshNodes.getLinearJacobian(ff);
                dJdv = obj.meshNodes.getLinearJacobianSensitivity(ff);
                
                [I, dIdJ, dIdu] = obj.elementIntegralFunctional(func, dfunc_du, jacobian);
                dIdv = multiplyTensors.txt(dIdJ, 2, dJdv, 4, 1:2, 1:2);
                
                F = F + I;
                dFdu(iGlobal) = dFdu(iGlobal) + dIdu;
                dFdv(iGlobal) = dFdv(iGlobal) + dIdv;
            end
            
        end
        
        function [F, dFdv, dFdu] = pointEvaluationFunctional(obj, pointFunction, xy, u)
            % [F, dFdp, dFdu] = pointEvaluationFunctional(pointFunction, xy, u)
            %
            % Calculate f(u) and df/du at a single point.
            
            assert(length(xy) == 2, 'Functional must be evaluated at a single point');
            
            interpMatrix = obj.meshNodes.getInterpolationOperator(xy(1), xy(2));
            
            uEval = interpMatrix*u;
            
            [F, dfunc_du] = pointFunction(uEval);
            dFdu = interpMatrix * dfunc_du;
            
            assert(isrow(dFdu));
            
            dFdv = [];
            
        end
        
        
    end % methods
    
    
    
end