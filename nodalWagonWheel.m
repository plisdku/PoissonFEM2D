function [xyNodes, faces] = nodalWagonWheel(nOuter, N)
% [xyNodes, faces] = nodalWagonWheel(nOuter, N)

assert(N >= 2, 'Nodes per edge must be at least 2');

[vertices, faces] = VVMesh.wagonWheel(nOuter);

lng = LinearNodalGeometry(faces, vertices(:,1:2), N);
xyNodes = lng.getNodeCoordinates();