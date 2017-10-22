function [xyNodes, faces] = nodalMesh(faces, vertices, N)

lng = LinearNodalGeometry(faces, vertices, N);
xyNodes = lng.getNodeCoordinates();