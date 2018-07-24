function [xyNodes, faces] = nodalMesh(faces, vertices, N)
% [xyNodes, faces] = nodalMesh(faces, vertices, N)

lng = PoissonFEM2D.LinearNodalGeometry(faces, vertices, N);
xyNodes = lng.getNodeCoordinates();