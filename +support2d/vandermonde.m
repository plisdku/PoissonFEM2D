function [V, dVdr, dVds] = vandermonde(N, rr, ss)
% Create vandermonde matrix at given points

[is,js] = ndgrid(0:N-1,0:N-1);
pickThese = (is+js < N);
is = is(pickThese);
js = js(pickThese);

if nargin < 3
    [rr,ss] = support2d.nodes2d(N);
end

for nn = 1:length(is)
    [z, dzdr, dzds] = support2d.jacobiNormalized2D(N, is(nn), js(nn), rr, ss);
    V(:,nn) = z;
    dVdr(:,nn) = dzdr;
    dVds(:,nn) = dzds;
end

