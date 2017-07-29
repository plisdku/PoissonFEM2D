function [iCorners, iEdgeCenters, iCenters] = classifyNodes(N)
% [iCorners, iEdgeCenters, iCenters] = classifyNodes(N)
%
% Classify the nodes of an element of order N.  Return indices of
% - three corners, [v0 v1 v2]
% - edges excluding corners, { [v00 v01 v02...], [v10 v11 v12...], ... }
% - remaining nodes in the middle, [v0 v1 v2 ...]
%
% The edge nodes are returned in counter-clockwise order (positive
% orientation).

allNodeIndices = 1:( N*(N+1)/2 );

[r,s] = support2d.nodes2d(N);

FUDGE = 1e-6;

isEdge1 = s < -1+FUDGE;
isEdge2 = s+r > -FUDGE;
isEdge3 = r < -1+FUDGE;
isEdge = isEdge1 | isEdge2 | isEdge3;

isCorner1 = isEdge1 & isEdge3;
isCorner2 = isEdge1 & isEdge2;
isCorner3 = isEdge2 & isEdge3;
isCorner = isCorner1 | isCorner2 | isCorner3;

isCenter = ~isEdge;

isEdgeCenter = ~(isCenter | isCorner);

%%

iCorners = allNodeIndices(isCorner);
%iEdgeCenters = allNodeIndices(isEdgeCenter);
iCenters = allNodeIndices(isCenter);

iEdgeCenters = { allNodeIndices(isEdgeCenter & isEdge1),
    allNodeIndices(isEdgeCenter & isEdge2),
    fliplr(allNodeIndices(isEdgeCenter & isEdge3)) };  % Note orientation!!

