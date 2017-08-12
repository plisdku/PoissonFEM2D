%% Test element quadrature tiringly

hack_verts = [0,0; 1,0; 0,1];
hack_faces = [1,2,3];
hack_N = 8;

hack_mesh = TriNodalMesh(hack_N, hack_faces, hack_verts);
hack_fem = PoissonFEM2D(hack_mesh);
jacobian = hack_fem.meshNodes.getLinearJacobian(1);
qv = hack_fem.elementQuadratureVector(jacobian);

%% Now let's do some integrals!
% 

xy = hack_mesh.getFaceNodeCoordinates(1);
integral_oracle = @(nn,mm) factorial(mm+1)*factorial(nn)/factorial(mm+nn+2)/(mm+1);

%%

ns = 0:6;
ms = 0:6;

rel_errs = [];


for n = ns
    for m = ms

        fprintf('For n = %d, m = %d:\n', n, m);

        integrand = xy(:,1).^n .* xy(:,2).^m;
        quadrature_result = qv'*integrand;
        exact_result = integral_oracle(n, m);
        err = quadrature_result - exact_result;
        rel_err = abs(err / exact_result);

        fprintf('Expected %e\n', integral_oracle(n,m));
        fprintf('Received %e\n\n', quadrature_result);
        fprintf('REL ERR  %e\n\n', rel_err);
        
        rel_errs(n+1,m+1) = rel_err;
    end
end

figure(3); clf
imagesc(ns, ms, log10(rel_errs)');
xlabel('n');
ylabel('m');
colorbar
title(sprintf('log10 relative error, order=%d', hack_N));