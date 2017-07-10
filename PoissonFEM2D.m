classdef PoissonFEM2D < handle
    
    properties
        meshNodes;
        Dr;  % differentiation matrix on basis element
        Ds;  % differentiation matrix on basis element
        Q;   % quadrature matrix on basis element
    end
    
    
    methods
        
        % ---- CONSTRUCTOR
        function obj = PoissonFEM2D(meshNodes)
            
            obj.meshNodes = meshNodes;
            
            rs = obj.meshNodes.nodesLocal_rs(:,1);
            ss = obj.meshNodes.nodesLocal_rs(:,2);
            
            [obj.Dr, obj.Ds] = support2d.gradients(obj.meshNodes.N, rs, ss);
            obj.Q = support2d.quadratureKernel(obj.meshNodes.N, rs, ss);
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
        
        
        function [C, dC_dJ_ij] = elementNeumannMatrix(obj, jacobian1d)
            
            error('Unimplemented');
            
            % C = L'*Q1d*det(J1d).  Must be done per edge.
            
        end % elementNeumannMatrix
        
        
        % ---- System Matrices
        
        function [systemMatrix, rhsMatrix, DsystemMatrix_dv, DrhsMatrix_dv] = systemMatrix(obj, varargin)
            % [A, B, dAdJ, dBdJ] = systemMatrix(obj)
            % [A, B, dAdJ, dBdJ] = systemMatrix(obj, jacobianOverride)
            %
            % The jacobian override is for sensitivity testing.
            
            if nargin == 2
                jacobianOverride = varargin{1};
            else
                jacobianOverride = [];
            end
            
            numNodes = obj.meshNodes.getNumNodes();
            systemMatrix = sparse(numNodes, numNodes);
            rhsMatrix = sparse(numNodes, numNodes);

            % We will accumulate the system matrix sensitivities with respect to
            % elemental Jacobians, first.
            %
            % DsystemMatrix_dJ{ii,jj} is the matrix of sensitivity to J(ii,jj).

            DsystemMatrix_dJ = cell(2,2);
            DrhsMatrix_dJ = cell(2,2);
            for nn = 1:4
                DsystemMatrix_dJ{nn} = sparse(numNodes, numNodes);
                DrhsMatrix_dJ{nn} = sparse(numNodes, numNodes);
            end
            
            % Fill in the matrices!
            numFaces = size(obj.meshNodes.faces,1);

            for ff = 1:numFaces
                % Pieces of the system matrix
                if isempty(jacobianOverride)
                    jac = obj.meshNodes.getLinearJacobian(ff);
                else
                    jac = jacobianOverride;
                    %disp('JACKED')
                end

                [M1, dM1dJ] = obj.elementPotentialMatrix(jac);
                [M2, dM2dJ] = obj.elementChargeMatrix(jac);

                iGlobal = obj.meshNodes.local2global(ff);

                systemMatrix(iGlobal,iGlobal) = systemMatrix(iGlobal, iGlobal) + M1;
                rhsMatrix(iGlobal, iGlobal) = rhsMatrix(iGlobal, iGlobal) + M2;

                for ii = 1:2
                    for jj = 1:2
                        DsystemMatrix_dJ{ii,jj}(iGlobal,iGlobal) = DsystemMatrix_dJ{ii,jj}(iGlobal,iGlobal) + dM1dJ{ii,jj};
                        DrhsMatrix_dJ{ii,jj}(iGlobal,iGlobal) = DrhsMatrix_dJ{ii,jj}(iGlobal, iGlobal) + dM2dJ{ii,jj};
                    end
                end

            end

            
        end % systemMatrix()
        
        
    end % methods
    
    
    
end