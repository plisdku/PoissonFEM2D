function rs = xy2rs(xyTri, xy)
% xy = xy2rs(rsTri, rs)

% Triangle vertices are three column vectors
assert(size(xyTri,1) == 2);
assert(size(xyTri,2) == 3);

% Input points are column vectors
assert(size(xy,1) == 2);

[T, v0] = support2d.rs2xy_affineParameters(xyTri);

%xy = bsxfun(@plus, v0, T*rs);
rs = T \ bsxfun(@minus, xy, v0);

