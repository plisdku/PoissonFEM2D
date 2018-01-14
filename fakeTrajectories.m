function [xOut,yOut] = fakeTrajectories(x0, y0, x1, y1, numT, numParticles)

% numParticles = 10;
% numT = 100;
% x0 = 0;
% y0 = 0;
% x1 = 5;
% y1 = 2;

xOut = zeros(numT, numParticles);
yOut = zeros(numT, numParticles);

ts = linspace(0, 1, numT);

for nn = 1:numParticles
    p = [randn([1,4]), 0] ./ [6,3,2,1,1];
    xs = polyval(p, ts) + linspace(x0,x1,numT);
    p = [randn([1,4]), 0] ./ [6,3,2,1,1];
    ys = polyval(p, ts) + linspace(y0,y1,numT);
    
    xOut(:,nn) = xs;
    yOut(:,nn) = ys;
end

%figure(10); clf
%plot(xOut, yOut);
