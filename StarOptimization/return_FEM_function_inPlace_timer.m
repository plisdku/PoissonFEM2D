function [F, dFdp] = return_FEM_function_inPlace_timer(fem, p, Lx, Ly)
    tic
    fprintf('Forward solution... ');
    
    [~] = fem.instantiateProblem(p);
   
    [femProblem, dDirichlet_dp, dnx_dp, dny_dp] = fem.adjustProblem(p);
    
    [idxN, ~,~ ] = find(abs(dnx_dp) + abs(dny_dp));
    idxN = unique(idxN(:));
    xyGeomNodes = femProblem.poi.tnMesh.xyNodes;
    geometry = fem.instantiatedGeom.geometry;
   
    PoissonFEM2D.plotGeom(107, geometry)
    PoissonFEM2D.plotGeomNodes(325, xyGeomNodes, idxN)
    
    %measBox = [-55e-3, 0, 55e-3, 1e-3];
    measBox = [-55e-3, 0, 55e-3, 2.5e-4];
    measNxy = [380, 18500];%[100 350];%[350, 18000];
    
    
    
    femProblem.solveForwardCartesian(measBox(1:2), measBox(3:4), measNxy);
    %disp('just forward')
    toc
    %tic
    %femProblem.solveCartesian(measBox(1:2), measBox(3:4), measNxy);
    %toc
    fprintf('Solved Forward \n')
    
    fprintf('Start VV \n')
    tic
    [particles, hit_objective] = SetupParticles_Array();
    [VV] = ElectronSetup_Array(femProblem.uCartesian, measBox, measNxy, particles, hit_objective); 
    fprintf('Total VV (forward/dual)')
    toc
    
    fprintf('Adjoint solution... ');
    tic
    femProblem.solveAdjointCartesian_new(VV.dFdV, measBox(1:2), measBox(3:4), measNxy,idxN)
    %femProblem.solveAdjointCartesian(VV.dFdV, idxN);
    
    fprintf('complete.\n');
    dFdp = femProblem.dF_dxy(:,1)' * dnx_dp + femProblem.dF_dxy(:,2)' * dny_dp ...
        + femProblem.dF_dDirichlet * dDirichlet_dp;

    F = VV.Fval;
    toc
    PoissonFEM2D.FullResultsPlot(789, femProblem, VV, geometry, idxN,...
    measBox, xyGeomNodes, Lx, Ly)

%     figure(); clf
%     hold on

    
    disp('time not for gradient')
    disp(femProblem.timer_variable)
%     
%     for i = 1:length(VV.ParticleArray)
%         plot(VV.ParticleArray(i).xx,VV.ParticleArray(i).yy,'k', 'LineWidth', 1)
%     %plot(VV.ParticleArray(i).xv(ix_x(1)),VV.ParticleArray(i).xv(ix_y(1)),'rx')
%     end
% 
%     figure(324)
%     for i = 1:length(VV.ParticleArray)
%         plot(VV.ParticleArray(i).xx,VV.ParticleArray(i).yy,'k', 'LineWidth', 1)
%     %plot(VV.ParticleArray(i).xv(ix_x(1)),VV.ParticleArray(i).xv(ix_y(1)),'rx')
%     end
%     
end 