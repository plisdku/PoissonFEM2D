%% Test a whole adjoint thing

lx = [0, 1, 1, 0];
ly = [0, 0, 1, 1];

in_lx = 0.5 + 0.05*[-1, -1, 1, 1];
in_ly = 0.5 + 0.05*[-1, 1, 1, -1];

density = 4;
[domainV,domainF] = meshPolygon(lx, ly, density, in_lx, in_ly);

figure(1); clf
VVMesh.plotFV(domainF, domainV, 'k-');
patch('Faces', domainF, 'Vertices', domainV, 'FaceColor', 'r', 'FaceAlpha', 0.1, 'EdgeAlpha', 0);

%% Make an FEM object

dirichletPredicate = @(v1,v2) norm(v1) < 0.25 || norm(v2) < 0.25;
femp = FEMProblem(N, domainF, domainV, dirichletPredicate);

%% Define some shit

x0 = 0.25;
y0 = 0.5;
sigma = 0.05;
freeCharge = @(x,y) exp( (-(x-x0).^2 - (y-y0).^2)/(2*sigma^2));

dirichletVal = @(xy) zeros(size(xy,1),1);

neumannVal = @(xy) zeros(size(xy,1),1);

measPt = [0.2, 0.22];

femp.setSources(freeCharge, dirichletVal, neumannVal);
femp.solve(measPt);

%% Perturb (Dirichlet)

delta = 1e-6;
femp2 = FEMProblem(N, domainF, domainV, dirichletPredicate);
femp2.setSources(freeCharge, dirichletVal, neumannVal);
femp2.u0_dirichlet(1) = femp2.u0_dirichlet(1) + delta;
femp2.solve(measPt);

dF_meas = (femp2.F - femp.F)/delta;
dF_calc = femp.dF_dud(1);

%% Perturb (Neumann)

femp2.setSources(freeCharge, dirichletVal, neumannVal);
femp2.en_neumann(1) = femp2.en_neumann(1) + delta;
femp2.solve(measPt);

dF_meas = (femp2.F - femp.F)/delta;
dF_calc = femp.dF_den(1);

%% Perturb (free charge)

femp2.setSources(freeCharge, dirichletVal, neumannVal);
femp2.freeCharge(10) = femp2.freeCharge(1) + delta;
femp2.solve(measPt);

dF_meas = (femp2.F - femp.F)/delta;
dF_calc = femp.dF_df(10);

%% Perturb (vertices)

iXY = 1;
for iVert = 1

    femp2.fem.meshNodes.vertices(iVert,iXY) = femp2.fem.meshNodes.vertices(iVert,iXY) + delta;
    femp2.setSources(freeCharge, dirichletVal, neumannVal);
    femp2.solve(measPt);

    dF_meas = (femp2.F - femp.F)/delta;
    dF_calc = femp.dFdv_total(iVert,iXY);
    
    
end


% TODO: check all the pieces


