%% Unit tests for MeshNodes!!

checkClose = @(a,b) assert( norm(a-b) < 1e-6 ); %, sprintf('(%g) is not close to (%g)', a, b));

%%

triVerts = [0,0; 2,0; 0,1];
triFaces = [1, 2, 3];
N = 3;

%% constructor

triMeshNodes = MeshNodes(triFaces, triVerts, N);

% face edges


% getBoundaryEdges
% getBoundaryNodes
% getInteriorNodes
% local2global_corner, _edge, _face
% local2global
% getCornerNodes
% getEdgeNodes
% getFaceNodes
% getVertexNodeCoordinates
% getEdgeNodeCoordinates
% getFaceNodeCoordinates
% getFaceVertices
% getNodeCoordinates
% getLinearJacobian
%% getLinearJacobian1d

% In general for an edge along direction d, the Jacobian is d/2.
% first line direction is d = [2;0]: Jacobian is [1;0].

jac = triMeshNodes.getLinearJacobian1d(1, 1);
checkClose(jac, [1;0]);

jac = triMeshNodes.getLinearJacobian1d(1, 2);
checkClose(jac, [-1;0.5]);

jac = triMeshNodes.getLinearJacobian1d(1, 3);
checkClose(jac, [0; -0.5]);

% getLinearJacobianSensitivity
% getGradientOperators
% getInterpolationOperator


