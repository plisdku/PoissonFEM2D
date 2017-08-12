function xy = rs2xy(xyTri, rs)
% xy = rs2xy(rsTri, rs)

% Triangle vertices are three column vectors
assert(size(xyTri,1) == 2);
assert(size(xyTri,2) == 3);

% Input points are column vectors
assert(size(rs,1) == 2);

[T, x0] = support2d.rs2xy_affineParameters(xyTri);

xy = bsxfun(@plus, x0, T*rs);

