function f_eval = interpolate(f, xy_eval, xy_tri, V)
% f_eval = interpolate(f, xy_eval, xy_tri, V)

assert(size(xy_tri,1) == 2, 'xy_tri should be three column vectors of length 2');
assert(size(xy_tri,2) == 3, 'xy_tri should be three column vectors of length 2');
assert(size(xy_eval,1) == 2, 'xy_eval should be 2xN');

numNodes = size(V,1);
N = 0.5*(-1 + sqrt(1 + 8*numNodes));
%fprintf('N = %d\n', N);

rs_eval = support2d.xy2rs(xy_tri, xy_eval);

V_eval = support2d.vandermonde(N, rs_eval(1,:), rs_eval(2,:));

f_eval = V_eval * (V \ f);