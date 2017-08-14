function [DT, Dx0] = rs2xy_affineParameterSensitivities(xyTri, DxyTri)
% [DT, Dx0] = rs2xy_affineParameterSensitivities(xyTri, DxyTri)

% Triangle vertices are three column vectors
assert(size(xyTri,1) == 2, 'xyTri must be 2x3 (xy x numVertices)');
assert(size(xyTri,2) == 3, 'xyTri must be 2x3 (xy x numVertices)');
assert(size(DxyTri,1) == 2, 'DxyTri must be 2x3 (xy x numVertices)');
assert(size(DxyTri,2) == 3, 'DxyTri must be 2x3 (xy x numVertices)');

Dx0 = 0.5*( DxyTri(:,2) + DxyTri(:,3) );
DT = 0.5 * [ DxyTri(:,2)-DxyTri(:,1), DxyTri(:,3)-DxyTri(:,1) ];
