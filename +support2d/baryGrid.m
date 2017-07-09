function b = baryGrid(N)
% [b1 b2 b3] = baryGrid(N)
% 
% Return an evenly-spaced array of barycentric coords spanning a triangle

[is,js,ks] = support2d.orders2d(N);

b1 = is/(N-1);
b2 = js/(N-1);
b3 = ks/(N-1);

b = [b1 b2 b3];