classdef PoissonFEM2D < handle
    
    properties
        meshNodes;
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
            
            %dQdJ = cell(2,2);
%             for ii = 1:2
%                 for jj = 1:2
%                     dQdJ{ii,jj} = obj.Q*dDetJ_dJ(ii,jj);
%                 end
%             end
        end % elementQuadratureMatrix
            
            
        
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
        
        
        function [C, dCdJ] = elementNeumannMatrix(obj, jacobian1d)
            
            detJ = sqrt(det(jacobian1d'*jacobian1d));
            C = obj.Q1d * detJ;
            
            %dDetJ_dJ = detJ*transpose(invJac); % this is all 4 sensitivities in a matrix
            dDetJ_dJ = detJ*jacobian1d*inv(jacobian1d'*jacobian1d);
            
            dCdJ = zeros([size(C), 2]);
            
            for ii = 1:2
                dCdJ(:,:,ii) = obj.Q1d*dDetJ_dJ(ii);
            end
            
        end % elementNeumannMatrix
        
        % ---- System Matrices
        
        function [NM, DNM_dv] = neumannMatrix(obj)
            
            % The Neumann matrix is sum( [get edge nodes] * Q * detJ ).
            % Its raw dimensions are nNodes x nNodes and on the outside
            % we can bake it down to the size of only the Neumann boundary
            % nodes.
            
            numNodes = obj.meshNodes.topology.getNumNodes();
            numVertices = obj.meshNodes.topology.getNumVertices();
            
            NM = sparse(numNodes, numNodes);
            DNM_dv = cell(numVertices, 2);
            for nn = 1:numel(DNM_dv)
                DNM_dv{nn} = sparse(numNodes, numNodes);
            end
            
            [boundaryEdges, orientations] = obj.meshNodes.topology.getBoundaryEdges();
            
            numEdges = length(boundaryEdges);
            
            for ii = 1:numEdges
                ee = boundaryEdges(ii);
                oo = orientations(ii);
                
                dJdv = obj.meshNodes.getLinearJacobianSensitivity1d(ee, oo);
                
                iVertices = obj.meshNodes.topology.getEdgeVertices(ee,oo);
                iGlobal = obj.meshNodes.topology.getEdgeNodes(ee, oo);
                
                jac1d = obj.meshNodes.getLinearJacobian1d(ee, oo);
                [neuMat, dNdJ] = obj.elementNeumannMatrix(jac1d);
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
            
            numNodes = obj.meshNodes.topology.getNumNodes();
            numVertices = obj.meshNodes.topology.getNumVertices();
            rhsMatrix = sparse(numNodes, numNodes);
            
            DrhsMatrix_dv = cell(numVertices,2);
            
            for nn = 1:numel(DrhsMatrix_dv)
                DrhsMatrix_dv{nn} = sparse(numNodes, numNodes);
            end
            
            
            % Fill in the matrix!
            numFaces = obj.meshNodes.topology.getNumFaces();

            for ff = 1:numFaces
                jacobian = obj.meshNodes.getLinearJacobian(ff);
                dJdv = obj.meshNodes.getLinearJacobianSensitivity(ff);
                
                [M2, dM2dJ] = obj.elementChargeMatrix(jacobian);
                dM2dv = multiplyTensors.txt(dM2dJ, 4, dJdv, 4, 3:4, 1:2);
                
                iVertices = obj.meshNodes.topology.getFaceVertices(ff);
                iGlobal = obj.meshNodes.topology.getFaceNodes(ff);

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
            
            numNodes = obj.meshNodes.topology.getNumNodes();
            numVertices = obj.meshNodes.topology.getNumVertices();
            
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
            numFaces = obj.meshNodes.topology.getNumFaces();

            for ff = 1:numFaces
                
                jacobian = obj.meshNodes.getLinearJacobian(ff);
                dJdv = obj.meshNodes.getLinearJacobianSensitivity(ff);
                
                [M1, dM1dJ] = obj.elementPotentialMatrix(jacobian);
                dM1dv = multiplyTensors.txt(dM1dJ, 4, dJdv, 4, 3:4, 1:2);
                
                iVertices = obj.meshNodes.topology.getFaceVertices(ff);
                iGlobal = obj.meshNodes.topology.getFaceNodes(ff);
                
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
        
        
        
    end % methods
    
    
    
end