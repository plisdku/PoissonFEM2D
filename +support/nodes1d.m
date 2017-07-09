function rr = nodes1d(N)
% [r,s] = nodes2d(N)
% Get warped nodes-on-triangle grid

rr = support.gaussLobatto(N);
