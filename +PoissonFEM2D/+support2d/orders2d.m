function [is,js,ks] = orders2d(N)

[is,js] = ndgrid(0:N-1, 0:N-1);
is = is(:); js = js(:);
pickThese = (is+js < N);
is = is(pickThese);
js = js(pickThese);
ks = (N-1) - is - js;