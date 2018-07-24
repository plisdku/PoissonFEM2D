function [T, x0] = rs2xy_affineParameters(xyTri)
% [T, x0] = rs2xy_affineParameters(xyTri)

% Triangle vertices are three column vectors
assert(size(xyTri,1) == 2);
assert(size(xyTri,2) == 3);

x0 = 0.5*( xyTri(:,2) + xyTri(:,3) );
T = 0.5 * [ xyTri(:,2)-xyTri(:,1), xyTri(:,3)-xyTri(:,1) ];
