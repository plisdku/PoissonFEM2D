function dr = warp(gaussLobattoPoints, r)
% dr = warp(gaussLobattoPoints, r)
% Return the displacement of points on [-1,1] that yields Gauss-Lobatto
% points.
% 

N = numel(gaussLobattoPoints);

rr = linspace(-1, 1, N);
dr = spline(rr, gaussLobattoPoints, r) - r;

if any( r + dr > 1.0 | r + dr < -1.0)
    warning('wha')
end