function [f, dfdu] = multiplyByOne(u)
% [f, dfdu] = multiplyByOne(u)
%
% Returns [u, ones(size(u))].

f = u;
dfdu = ones(size(u));