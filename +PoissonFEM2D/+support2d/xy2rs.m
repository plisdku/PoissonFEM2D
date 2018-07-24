function rs = xy2rs(xyTri, xy)
% xy = xy2rs(rsTri, rs)
import PoissonFEM2D.*
% Triangle vertices are three column vectors
assert(size(xyTri,1) == 2, 'Triangle vertices must be three column vectors');
assert(size(xyTri,2) == 3, 'Triangle vertices must be three column vectors');

% Input points are column vectors
assert(size(xy,1) == 2, 'Query points must be column vectors');

[T, x0] = support2d.rs2xy_affineParameters(xyTri);

%xy = bsxfun(@plus, v0, T*rs);
rs = T \ bsxfun(@minus, xy, x0);

