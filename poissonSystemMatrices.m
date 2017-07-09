function [systemMatrix, rhsMatrix, DsystemMatrix, DrhsMatrix] = poissonSystemMatrices(meshNodes)
% [M1, M2] = poissonSystemMatrices(meshNodes)
% [M1, M2, DM1, DM2] = poissonSystemMatrices(meshNodes)
%
% Calculate the LHS and RHS matrix for the Poisson equation.  Optionally
% also obtain the sensitivity of the matrices.

%%
N = meshNodes.N;

[rr,ss] = support2d.nodes2d(N);
V = support2d.vandermonde(N, rr, ss);
Q = support2d.quadratureKernel(N, rr, ss);

%%
% How can I get the matrices I need?
% I need a jacobian for each triangle.  Ask MeshNodes for that.
% 

numNodes = meshNodes.getNumNodes();
systemMatrix = sparse(numNodes, numNodes);
rhsMatrix = sparse(numNodes, numNodes);

DsystemMatrix = sparse(numNodes, numNodes);
DrhsMatrix = sparse(numNodes, numNodes);

%% Sum up the elemental matrices: system matrix w/o boundary conditions

numFaces = size(meshNodes.faces,1);

for ff = 1:numFaces
    % Pieces of the system matrix
    jac = meshNodes.getLinearJacobian(ff);
    invJac = inv(jac);
    [Dx,Dy] = support2d.gradients_xy(N, rr, ss, invJac);
    
    Q_xy = Q*det(jac);
    
    M1 = (Dx'*Q_xy*Dx + Dy'*Q_xy*Dy);
    M2 = Q_xy;
    
    % Pieces of the primal matrix
    Djac = meshNodes.getLinearJacobianSensitivity(ff);
    DinvJac = -invJac * Djac * invJac;
    [dDx, dDy] = support2d.gradientSensitivities_xy(N, rr, ss, DinvJac);
    
    DdetJ = det(jac)*trace(invJac * Djac);
    
    DQ_xy = Q*DdetJ;
    
    DM1 = dDx'*Q_xy*Dx + Dx'*DQ_xy*Dx + Dx'*Q_xy*dDx + ...
        dDy'*Q_xy*Dy + Dy'*DQ_xy*Dy + Dy'*Q_xy*dDy;
    DM2 = DQ_xy;
    
    % Accumulate the global matrices
    
    iGlobal = meshNodes.local2global(ff);
    
    systemMatrix(iGlobal,iGlobal) = systemMatrix(iGlobal, iGlobal) + M1;
    rhsMatrix(iGlobal, iGlobal) = rhsMatrix(iGlobal, iGlobal) + M2;
    
    DsystemMatrix(iGlobal,iGlobal) = DsystemMatrix(iGlobal,iGlobal) + DM1;
    DrhsMatrix(iGlobal,iGlobal) = DrhsMatrix(iGlobal, iGlobal) + DM2;
end
