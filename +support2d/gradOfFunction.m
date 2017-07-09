function [gx, gy] = gradOfFunction(g, x0, y0, varargin)
%function [gx, gy] = gradHack(g, x0, y0, delta)

if nargin < 4
    delta = 1e-3;
else
    delta = varargin{1};
end

gx = (g(x0+delta,y0) - g(x0-delta,y0))/delta/2;
gy = (g(x0, y0+delta) - g(x0, y0-delta))/delta/2;