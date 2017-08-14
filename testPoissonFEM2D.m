%% Interpolation operator sensitivity

[vertices,faces] = VVMesh.wagonWheel(5);
vertices = vertices(:,1:2);

delta = 1e-6;
iVert = 2;
iXY = 2;
Dvertices = 0*vertices;
Dvertices(iVert, iXY) = 1.0;
vertices2 = vertices + delta*Dvertices;

N = 4;
tnMesh = TriNodalMesh(N, faces, vertices);
tnMesh2 = TriNodalMesh(N, faces, vertices2);

%figure(1); clf
%VVMesh.plotFV(faces, vertices, 'k');

fem = PoissonFEM2D(tnMesh);
fem2 = PoissonFEM2D(tnMesh2);

%% 

myFunc = @(x,y) sin(x) + cos(y);
[f, dfdv] = fem.evaluateOnNodes(myFunc);
f2 = fem2.evaluateOnNodes(myFunc);

Df_meas = (f2-f)/delta;
Df_calc = full(dfdv{iVert,iXY});

fprintf('Norm of Df_meas = %f\n', norm(Df_meas));
fprintf('Norm of error = %g\n', norm(Df_meas - Df_calc));
