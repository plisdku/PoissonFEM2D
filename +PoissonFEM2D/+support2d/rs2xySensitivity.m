function Dxy = rs2xySensitivity(xyTri, DxyTri, rs)
% Dxy = rs2xySensitivity(xyTri, DxyTri, rs)
import PoissonFEM2D.*
% Triangle vertices are three column vectors
assert(size(xyTri,1) == 2);
assert(size(xyTri,2) == 3);

% Input points are column vectors
assert(size(rs,1) == 2);

[T, x0] = support2d.rs2xy_affineParameters(xyTri);
[DT, Dx0] = support2d.rs2xy_affineParameterSensitivities(xyTri, DxyTri);

Dxy = bsxfun(@plus, Dx0, DT*rs);