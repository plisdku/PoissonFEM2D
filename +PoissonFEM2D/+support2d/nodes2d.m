function [r,s] = nodes2d(N)
% [r,s] = nodes2d(N)
% Get warped nodes-on-triangle grid
import PoissonFEM2D.*

[is,js,ks] = support2d.orders2d(N);

r = 2*(is/(N-1)) - 1;
s = 2*(js/(N-1)) - 1;

% Triangle corners:
triVerts = [-1, 1, -1; -1, -1, 1];
baryBasis = [triVerts(:,1)-triVerts(:,3), triVerts(:,2)-triVerts(:,3)];

% Barycentric coords
uv = baryBasis \ bsxfun(@plus, [r'; s'], -triVerts(:,3));

uu = uv(1,:)';
vv = uv(2,:)';
ww = 1 - uu - vv;

%uu = is/(N-1);
%vv = js/(N-1);
%ww = ks/(N-1);

% Blend and warp
xGL = support.gaussLobatto(N);
b1 = support2d.blend(uu,vv,ww).*support.warp(xGL, ww-vv);
b2 = support2d.blend(vv,ww,uu).*support.warp(xGL, uu-ww);
b3 = support2d.blend(ww,uu,vv).*support.warp(xGL, vv-uu);

% Warp direction
edge1 = triVerts(:,3)'-triVerts(:,2)';
edge2 = triVerts(:,1)'-triVerts(:,3)';
edge3 = triVerts(:,2)'-triVerts(:,1)';

warpVec = b1*edge1 + b2*edge2 + b3*edge3;

%figure(10); clf
%quiver(r, s, warpVec(:,1), warpVec(:,2))
%hold on
%plot(r, s, 'o')

% Final positions
r = r + 0.5*warpVec(:,1);
s = s + 0.5*warpVec(:,2);

