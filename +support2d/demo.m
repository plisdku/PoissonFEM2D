%% Demo of support2d tools

import support2d.*

%% Basis nodes

N = 10;
[r, s] = nodes2d(N);

figure(1); clf
plot(r, s, 'o')
title(sprintf('Nodes for N=%i', N));
xlabel('r')
ylabel('s')
grid on

%% Basis modes

V = vandermonde(N, r, s);

for nn = 1:size(V,2)
    figure(2); clf
    plot3(r, s, V(:,nn), 'o');
    xlabel('r')
    ylabel('s')
    title(sprintf('Mode %i for N=%i', nn, N));
    pause(0.1);
end

%% Vandermonde derivatives

[V, dVdr, dVds] = vandermonde(N, r, s);

modeNum = 3;
[ii,jj] = orders2d(N);


for nn = 1:length(r)
    r0 = r(nn);
    s0 = s(nn);
    
    g = @(x,y) jacobiNormalized2D(N, ii(modeNum), jj(modeNum), x, y);
    
    [dPdr, dPds] = gradOfFunction(g, r0, s0);
    
    fprintf('\nPosition %i (%f, %f):\n', nn, r0, s0)
    fprintf('dPdr: Vandermonde %f vs finite diff %f\n', dPdr, dVdr(nn,modeNum));
    fprintf('dPds: Vandermonde %f vs finite diff %f\n', dPds, dVds(nn,modeNum));
    
end

%% Derivative operators applied to test function

testFunc = @(x,y) sin(pi*x).*y + y.^3/2;
dtestFunc_dr = @(x,y) cos(pi*x).*y*pi;
dtestFunc_ds = @(x,y) sin(pi*x) + 3*y.^2/2;

%testFunc = @(x,y) x.*y;
%dtestFunc_dr = @(x,y) y;
%dtestFunc_ds = @(x,y) x;

%testFunc = @(x,y) sin(pi*x);
%dtestFunc_dr = @(x,y) cos(pi*x)*pi;
%dtestFunc_ds = @(x,y) 0*x;

figure(1); clf
plot3(r, s, testFunc(r,s), 'o')
grid on
xlabel('r'); ylabel('s');
title('Test function')

[Dr, Ds] = gradients(N, r, s);

dfdr = Dr*testFunc(r,s);
dfds = Ds*testFunc(r,s);

dfdr_exact = dtestFunc_dr(r,s);
dfds_exact = dtestFunc_ds(r,s);

relErr_r = abs( (dfdr - dfdr_exact)./dfdr_exact );
relErr_s = abs( (dfds - dfds_exact)./dfds_exact );

relErrNorm_r = norm(dfdr-dfdr_exact)/norm(dfdr_exact);
relErrNorm_s = norm(dfds-dfds_exact)/norm(dfds_exact);

fprintf('Relative error of df/dr and df/ds (%%):\n');
fprintf('  %2.2f%%   %2.2f%%\n', ([relErr_r relErr_s]*100))
fprintf('\nOverall relative error in the norm:\n');
fprintf('  %2.2f%%   %2.2f%%\n', 100*[relErrNorm_r relErrNorm_s])

%% Jacobian operators

x0 = [0, 0]';
x1 = [3.1, 1.1]';
x2 = [0.5, 2]';

triVerts = [x0 x1 x2];
jacobian = support2d.rs2xy_jacobian(triVerts);

r0 = -0.1;
s0 = -0.2;

xyFunc = @(r,s) support2d.rs2xy(triVerts, [r; s]);

[dxy_dr, dxy_ds] = support2d.gradOfFunction(xyFunc, r0, s0);

fprintf('Expected Jacobian:\n');
disp(jacobian)
fprintf('Measured Jacobian:\n');
disp([dxy_dr, dxy_ds])

%% Quadrature in (r,s) space: normalization of the basis elements

Q = support2d.quadratureKernel(N, r, s);

for modeNum = 1:4
    basisFunc = V(:,modeNum);
    integralOfSquare = basisFunc' * Q * basisFunc;
    fprintf('Integral of square of basis function %i = %f\n', ...
        modeNum, integralOfSquare);
end

%% Quadrature in (r,s) space: area of the basis triangle

uno = ones(size(r));

triArea = uno' * Q * uno;

fprintf('Area of triangle should be 2, and quadrature says %f\n', ...
    triArea);

%% Quadrature in (r,s) space: integral of some functions

int_of_r = -2/3;  % by hand
val = ones(size(r))'*Q*(r);

fprintf('Integral of r: expected %f got %f\n', int_of_r, val);

int_of_s = -2/3; % by hand
val = ones(size(r))'*Q*(s);

fprintf('Integral of s: expected %f got %f\n', int_of_s, val);

int_of_rs = 0;
val = ones(size(r))'*Q*(r.*s);

fprintf('Integral of r*s: expected %f got %f\n', int_of_rs, val);

%% Quadrature in (x,y) space: integral of some functions

% Unit simplex.  Why not, it's easy!!!
x0 = [0, 0]';
x1 = [1, 0]';
x2 = [0, 1]';

triVerts = [x0 x1 x2];
jacobian = support2d.rs2xy_jacobian(triVerts);

xy = support2d.rs2xy(triVerts, [r,s]');
x = xy(1,:)';
y = xy(2,:)';

Q = support2d.quadratureKernel(N, r, s, jacobian);

% Area of unit simplex should be 1/2

shouldBeHalf = ones(size(r))'*Q*ones(size(r));

fprintf('Area of unit simplex: expected %f got %f\n', 0.5, shouldBeHalf);

% Integral of some functions

int_of_x2 = 1/12;
val = ones(size(r))'*Q*x.^2;

fprintf('Integral of x on unit simplex: expected %f got %f\n', ...
    int_of_x2, val);

%% Partial derivatives in (x,y) space

% Pick any old thing!
x0 = [0, 0]';
x1 = [1, 0.5]';
x2 = [0, 1]';

triVerts = [x0 x1 x2];
jacobian = support2d.rs2xy_jacobian(triVerts);

xy = support2d.rs2xy(triVerts, [r,s]');
x = xy(1,:)';
y = xy(2,:)';

Q = support2d.quadratureKernel(N, r, s, jacobian);

invJac = inv(jacobian);

%Dx = Dr*invJac(1,1) + Ds*invJac(2,1);
%Dy = Dr*invJac(1,2) + Ds*invJac(2,2);

[Dx, Dy] = support2d.gradients_xy(N, r, s, invJac);

fprintf('Dx*x (1) and Dx*y (0):\n');
disp([Dx*x, Dx*y])

fprintf('Dy*x (0) and Dy*y (1):\n');
disp([Dy*x, Dy*y])

%% Sensitivity of inverse of Jacobian
% to elements of the Jacobian.

jacobian = [1, 2; 3, 4];
invJac = inv(jacobian);

delta = 1e-6;
for ii = 1:2
    for jj = 1:2
        
        dinvJ_dJij = -invJac(:,ii)*invJac(jj,:);
        
        % Test them...
        jac2 = jacobian;
        jac2(ii,jj) = jac2(ii,jj) + delta;
        invJac2 = inv(jac2);

        dInvJac_meas = (invJac2-invJac)/delta;

        relErr = norm(dInvJac_meas - dinvJ_dJij)/norm(dInvJac_meas);
        fprintf('Relative error(%g,%g) = %g\n', ii, jj, relErr);
    end
end


%% Sensitivity of determinant of Jacobian
% Test with a single perturbed Jacobians.

detJ = det(jacobian);
detJ2 = det(jacobian2);

ddetJ = detJ * trace(invJac * dJacobian);
ddetJ_meas = (detJ2-detJ)/delta;

fprintf('d|J| relative error in the norm:\n');
disp(abs(ddetJ_meas-ddetJ)/abs(ddetJ_meas));

%% Sensitivity of the determinant of the Jacobian
% to elements of the Jacobian

jacobian = [1, 2; 3, 4];
detJac = det(jacobian);
invJac = inv(jacobian);

% Put the 4 gradient elements into a matrix.
dinvJac_dJ = det(jacobian)*transpose(inv(jacobian));

delta = 1e-6;
for ii = 1:2
    for jj = 1:2
        
        % Test them...
        jac2 = jacobian;
        jac2(ii,jj) = jac2(ii,jj) + delta;
        
        ddetJ_meas = (det(jac2)-det(jacobian))/delta;
        
        relErr = abs(ddetJ_meas - dinvJac_dJ(ii,jj))/abs(ddetJ_meas);
        fprintf('Relative error(%g,%g) = %g\n', ii, jj, relErr);
    end
end


%% Sensitivity of partial derivatives to change of Jacobian
% Test with a single perturbed Jacobian.

dJacobian = [1.1, 2.2; -3.3, -4.4];

dInvJacobian = -invJac*dJacobian*invJac;

% Calculate new Dx and Dy matrices for a perturbed jacobian
delta = 1e-6;

jacobian2 = jacobian + delta * dJacobian;

[Dx2, Dy2] = support2d.gradients_xy(N, r, s, inv(jacobian2));

[dDx, dDy] = support2d.gradientSensitivities_xy(N, r, s, dInvJacobian);

dDx_meas = (Dx2-Dx)/delta;
dDy_meas = (Dy2-Dy)/delta;

fprintf('dDx relative error in the norm:\n');
disp(norm(dDx_meas-dDx)/norm(dDx))
fprintf('dDy relative error in the norm:\n');
disp(norm(dDy_meas-dDy)/norm(dDy))

%% Sensitivity of partial derivatives to change of each Jacobian element

jacobian = [1, 2; 3, 4];
invJac = inv(jacobian);
[Dx, Dy] = support2d.gradients_xy(N, r, s, inv(jacobian));

delta = 1e-6;
for ii = 1:2
    for jj = 1:2
        
        jac2 = jacobian;
        jac2(ii,jj) = jac2(ii,jj) + delta;
        
        [Dx2, Dy2] = support2d.gradients_xy(N, r, s, inv(jac2));
        
        dDx_meas = (Dx2-Dx)/delta;
        dDy_meas = (Dy2-Dy)/delta;
        
        % Use crazy formulae
        dDx = -invJac(1,ii)*invJac(jj,1)*Dr - invJac(2,ii)*invJac(jj,1)*Ds;
        dDy = -invJac(1,ii)*invJac(jj,2)*Dr - invJac(2,ii)*invJac(jj,2)*Ds;
        
        relErr_Dx = norm(dDx - dDx_meas)/norm(dDx_meas);
        relErr_Dy = norm(dDy - dDy_meas)/norm(dDy_meas);
        
        fprintf('Relative error in dDx/dJ(%i,%i) = %g\n', ii, jj, relErr_Dx);
        fprintf('Relative error in dDy/dJ(%i,%i) = %g\n', ii, jj, relErr_Dy);
    end
end

%% Sensitivity of affine parameters in one triangle

% Unit simplex, in [iVert xy] orientation like I use in data structures
xy = [0, 0; 1, 0; 0, 1];
Dxy = [0, 1; 0, 1; 1, 0];  % just... something.  try several.
delta = 1e-6;
xy2 = xy + delta*Dxy;

[T, v0] = support2d.rs2xy_affineParameters(xy');
[T2,v02] = support2d.rs2xy_affineParameters(xy2');




