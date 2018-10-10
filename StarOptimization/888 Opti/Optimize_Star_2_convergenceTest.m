import PoissonFEM2D.*
%%
mesh_factors = [50 20 10 5 1 0.9 0.75 0.5];
box_factors  = [1 1.5 2 2.5 3 4];
measNxy = [380, 18500];
cartesianField = {};
gradients_xy = {};
dnxs = {};
dnys = {};
gradients_p = {};

for ii = 1:length(mesh_factors)
    
    for jj = 1:length(box_factors)
    
    
        Lx_outer = 60e-3*box_factors(jj);
        isAxisymmetric = 1;

        N_field = 5;
        N_geom = 2;
        N_quad = N_field + isAxisymmetric;

        s = mesh_factors(ii)*1.28e-3;%4e-3; %1.3e-3; % mesh scale
        ratio = 0.3;
        geom2d = PoissonFEM2D.ParameterizedGeometry2D();


        Vb = 29e3;
        D1 = 3e-3;

        Ly = 45e-3*box_factors(jj);

        D12 = 5e-3;
        D23 = D12;
        L1 = 2e-3;
        L2 = 5e-3;
        L3 = L1;

        r_i = 3e-3;
        r_o = 11e-3;

        Ncontour_vec = [3 3 3 3];%[100 100 100 100];%[30 30 30 30];

        E2_V1 = [-L2/2, r_i];
        E2_V2 = [L2/2, r_i];
        E2_V3 = [L2/2, r_o];
        E2_V4 = [-L2/2, r_o];
        E2 = [E2_V1; E2_V2; E2_V3; E2_V4];
        [E2_x, E2_y, s2_vec] = BoxVector(E2, Ncontour_vec, @(N) s_vec_lin(N,1,2));
        l2 = length(E2_x);
        l22 = 2*l2;

        E1_V1 = [-L2/2, r_i] - [D12+L1 0];
        E1_V2 = [-L2/2, r_i] - [D12 0];
        E1_V3 = [-L2/2, r_o] - [D12 0];
        E1_V4 = [-L2/2, r_o] - [D12+L1 0];
        E1 = [E1_V1; E1_V2; E1_V3; E1_V4];
        [E1_x, E1_y, s1_vec] = BoxVector(E1, Ncontour_vec, @(N) s_vec_lin(N,1,2));
        l1 = length(E1_x);
        l12 = 2*l1;
        %
        E3_V1 = [L2/2, r_i] + [D23 0];
        E3_V2 = [L2/2, r_i] + [D23+L3 0];
        E3_V3 = [L2/2, r_o] + [D23+L3 0];
        E3_V4 = [L2/2, r_o] + [D23 0];
        E3 = [E3_V1; E3_V2; E3_V3; E3_V4];
        [E3_x, E3_y, s3_vec] = BoxVector(E3, Ncontour_vec, @(N) s_vec_lin(N,1,2));
        l3 = length(E3_x);
        l32 = 2*l3;
        %
        plotElectrodePoints(E1_x, E1_y, E2_x, E2_y, E3_x, E3_y)

        [E1x_func, E1y_func, end_p11] = createGeomFunction(E1_x, E1_y, 0);
        [E2x_func, E2y_func, end_p12] = createGeomFunction(E2_x, E2_y, end_p11);
        [E3x_func, E3y_func, end_p13] = createGeomFunction(E3_x, E3_y, end_p12);

        r_alpha1 = sqrt( (0.5*L1)^2 + ((r_o-r_i)*0.5)^2);
        r_alpha2 = sqrt( (0.5*L2)^2 + ((r_o-r_i)*0.5)^2);

        r1_start = [r_alpha1 ((r_o-r_i)*0.5) r_alpha1 L1/2 r_alpha1 ((r_o-r_i)*0.5) r_alpha1 L1/2];
        r2_start = [r_alpha2 ((r_o-r_i)*0.5) r_alpha2 L2/2 r_alpha2 ((r_o-r_i)*0.5) r_alpha2 L2/2];
        r3_start = r1_start;

        center1_start = [-L2/2-D12-L1/2 ((r_o-r_i)*0.5)+r_i];
        center2_start = [0 ((r_o-r_i)*0.5)+r_i];
        center3_start = [L2/2+D23+L3/2  ((r_o-r_i)*0.5)+r_i];

        alpha = atan(((r_o-r_i)*0.5)/(0.5*L1));
        alpha2 = atan(((r_o-r_i)*0.5)/(0.5*L2));
        angles1_start = [pi+alpha 3/2*pi -alpha 0 alpha 1/2*pi pi-alpha pi];
        angles2_start = [pi+alpha2 3/2*pi -alpha2 0 alpha2 1/2*pi pi-alpha2 pi];
        angles3_start = angles1_start;

        x_vec = r1_start .* cos(angles1_start) + center1_start(1);
        y_vec = r1_start .* sin(angles1_start) + center1_start(2);


        %%

        [x1_star, y1_star, l1, end_p1] = xy_star(length(angles1_start), 0, r1_start,...
            center1_start, angles1_start);
        [x2_star, y2_star, l2, end_p2] = xy_star(length(angles1_start), end_p1, r2_start,...
            center2_start', angles2_start);
        [x3_star, y3_star, l3, end_p3] = xy_star(length(angles1_start), end_p2, r3_start,...
            center3_start, angles3_start);


        geom2d.addContour(@(p) [-Lx_outer -L2/2-D12-L1 0 L2/2+D23+L3 Lx_outer Lx_outer 0 -Lx_outer], @(p) [0, 0, 0, 0, 0, Ly, Ly, Ly], s*[15, ratio, ratio, ratio, 15, 25, 25, 25], 1, 1:8);


        geom2d.addContour(x1_star, y1_star, s*ratio*[1 1 1 1.5 2 2 2 1.5], 2, 1:l1);

        geom2d.addContour(x2_star, y2_star, s*ratio*[1 1 1 1.5 2 2 2 1.5], 3, 1:l2);

        geom2d.addContour(x3_star, y3_star, s*ratio*[1 1 1 1.5 2 2 2 1.5], 4, 1:l3);

        safetyplotGeom(geom2d, zeros(1,end_p3))

        fem = PoissonFEM2D.FEMInterface(geom2d, N_field, N_geom, N_quad, isAxisymmetric);
        fem.setNeumann(1, @(p,x,y) 0.0);

        fem.setDirichlet(2, @(p,x,y) 0 );
        fem.setDirichlet(3, @(p,x,y) Vb );
        fem.setDirichlet(4, @(p,x,y) 0 );

        fem.setFreeCharge(@(p,x,y) 0.0);


        %%


        p = zeros(1,end_p3)';

        [~] = fem.instantiateProblem(p);

        [femProblem, dDirichlet_dp, dnx_dp, dny_dp] = fem.adjustProblem(p);

        [idxN, ~,~ ] = find(abs(dnx_dp) + abs(dny_dp));
        idxN = unique(idxN(:));
        xyGeomNodes = femProblem.poi.tnMesh.xyNodes;
        geometry = fem.instantiatedGeom.geometry;

        PoissonFEM2D.plotGeom(107, geometry)
        PoissonFEM2D.plotGeomNodes(325, xyGeomNodes, idxN)

        measBox = [-55e-3, 0, 55e-3, 2.5e-4];
        measNxy = [380, 18500];

        fprintf('Forward solution... ');

        tic
        femProblem.solveCartesian(measBox(1:2), measBox(3:4), measNxy);
        toc
        fprintf('Solved Forward ')

        [particles, hit_objective] = SetupParticles_star();
        [VV] = ElectronSetup_star(femProblem.uCartesian, measBox, measNxy, particles, hit_objective); 

        fprintf('Adjoint solution... ');
        tic
        femProblem.solveAdjointCartesian(VV.dFdV, idxN);
        toc
        fprintf('complete.\n');
        dFdp = femProblem.dF_dxy(:,1)' * dnx_dp + femProblem.dF_dxy(:,2)' * dny_dp ...
            + femProblem.dF_dDirichlet * dDirichlet_dp;

        F = VV.Fval;


        cartesianField{ii, jj} = femProblem.uCartesian;
        gradients_xy{ii,jj} = femProblem.dF_dxy;

        dnxs{ii,jj} = dnx_dp;
        dnys{ii,jj} = dny_dp;
        
        gradients_p{ii,jj} = dFdp;

        disp('iteration done')
        save('intermediateResult', 'cartesianField', 'gradients_xy', 'dnys', 'dnxs', 'gradients_p')
    
    end
    
    
end
    

%%



finField = cartesianField{end,end};

for kk = 1:length(mesh_factors)
    for ll = 1:length(box_factors)
        caField = cartesianField{kk,ll};
        diff_cartesianField(kk,ll) = norm(caField(:) - finField(:),2) ./ norm(finField(:),2);
        
    end 
end

figure(641)
yvalues = [50 20 10 5 1 0.9 0.75 0.5];
xvalues = [1 1.5 2 2.5 3 4];
h = heatmap(xvalues, yvalues, log10(diff_cartesianField));

h.XLabel = 'BoxExpansionFactor';
h.YLabel = 'MeshExpansionFactor';


finGrad = gradients_p{end,end};


for kk = 1:length(mesh_factors)
    for ll = 1:length(box_factors)
        gradi = gradients_p{kk,ll};
        diff_gradient(kk,ll) = norm(gradi(:) - finGrad(:),2) ./norm(finGrad(:),2); 
        
    end
end

figure(642)
yvalues = [50 20 10 5 1 0.9 0.75 0.5];
xvalues = [1 1.5 2 2.5 3 4];
h = heatmap(xvalues, yvalues, log10(diff_gradient));

h.XLabel = 'BoxExpansionFactor';
h.YLabel = 'MeshExpansionFactor';





%%

figure(652)
imagesc(log10(xvalues), log10(yvalues), log10(diff_cartesianField))


figure(6662)
clf
hold on
for iii = 1:length(box_factors)
    loglog(yvalues, diff_cartesianField(:,iii))
end
figure(6663)
clf
hold on
for jjj = 1:length(mesh_factors)
    loglog(xvalues, diff_cartesianField(jjj,:))
end

figure(662)
loglog(xvalues, diff_cartesianField')
xlabel('Box factor')
title('cartesian Field convergance')
figure(663)
loglog(xvalues, diff_gradient')
xlabel('Box factor')
title('gradient convergance')


figure(672)
loglog(yvalues, diff_cartesianField')
xlabel('Meshfactor')
title('gradient convergance')

figure(673)
loglog(yvalues, diff_gradient')
xlabel('MeshFactor')
title('cartesian Field convergance')







































































    
    