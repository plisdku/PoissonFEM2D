%% Input parameters for this test case
%clc 
%%
% Setting natural constants
% in case that they are needed
function [VV] = ElectronSetup_star(u, measBox, measNxy, particles, hit_objective)



    xyz = [measBox(1), measBox(3); measBox(2),measBox(4); 0,2];
    Nxy = [measNxy(1), measNxy(2), 2];
    VV = cylindricalVelocity_Verletdebug(xyz, Nxy, particles, hit_objective, u);
    VV = VV.calculateF_dF;

    fprintf('Trajectory and Dual Trajectory Calculations done')



end