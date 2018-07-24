function Drs = xy2rsSensitivity(xyTri, DxyTri, xy)
% Drs = xy2rs_sensitivity(xyTri, DxyTri, xy)
import PoissonFEM2D.*
% Triangle vertices are three column vectors
assert(size(xyTri,1) == 2, 'Triangle vertices must be three column vectors');
assert(size(xyTri,2) == 3, 'Triangle vertices must be three column vectors');

assert(isequal(size(DxyTri), size(xyTri)), 'DxyTri must be the same size as xyTri');

% Input points are column vectors
assert(size(xy,1) == 2, 'Query points must be column vectors');

[T,x0] = support2d.rs2xy_affineParameters(xyTri);
[DT, Dx0] = support2d.rs2xy_affineParameterSensitivities(xyTri, DxyTri);

invT = inv(T);
DinvT = -invT*DT*invT;

Drs = bsxfun(@plus, -invT*Dx0, DinvT*bsxfun(@plus, xy, -x0));
