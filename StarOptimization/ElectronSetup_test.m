%% Input parameters for this test case
%clc 
%%
% Setting natural constants
% in case that they are needed
function [VV] = ElectronSetup_test(u, measBox, measNxy, particles, particles2, hit_objective, hit_objective2)



    xyz = [measBox(1), measBox(3); measBox(2),measBox(4); 0,2];
    Nxy = [measNxy(1), measNxy(2), 2];
    VV = cylindricalVelocity_Verletdebug2(xyz, Nxy, particles, hit_objective, u);
    VV = VV.calculateF_dF;
    VV2 = cylindricalVelocity_VerletArray(xyz, Nxy, particles2, hit_objective2, u);
%     for ii = 1:15
%         
%         xv3(:,ii) = VV.ParticleArray(ii).xv;
%         %accel3(:,:,ii) = VV.ParticleArray
%         
%     end 
%     VV2.ParticleArray.xv = xv3;
    %VV2 = VV2.calculateF;
    VV2 = VV2.calculateF_dF;
%     nP = 15;
%     Traj1 = VV2.ParticleArray.xv;
%     Traj2 = zeros(9000,nP);
%     %figure()
%     %hold on
%     for i = 1:nP
%         figure()
%         hold on
%         Traj2(:,i) = VV.ParticleArray(i).xv;
%         plot(VV.ParticleArray(i).xx)
%         plot(VV2.ParticleArray.xx(:,i))
%     end
%     
%     Acc2 = VV2.ParticleArray.accelerations;
%     
%     Acc1 = zeros(1500,3,nP);
%     for i = 1:nP
%         Acc1(:,:,i) = VV.ParticleArray(i).accelerations;
%     end

%     for ii = 1:15
%         
%         xv3(:,ii) = VV.ParticleArray.xv;
%         
%     end 
    [G3, DG3] = hit_objective2(xv3);
    for ii = 1:15
        
        Diff(:,ii) =  abs((xv3(:,ii) - VV2.ParticleArray.xv(:,ii)) ./ xv3(:,ii));
        Diff(isnan(Diff)) = 0;
        figure()
        plot(Diff(:,ii))
        k = sum(Diff(:,ii) ~= 0); 
        title(sprintf('%d , %d', ii, k))
        max(abs(Diff(:,ii)))
        
    end 
    for ii = 1:15
        disp(ii)
        isequal(VV.D_I_1{ii},VV2.D_I_1{ii})
        isequal(VV.D_I_2{ii},VV2.D_I_2{ii})
       isequal(VV.D_I_3{ii},VV2.D_I_3{ii})
       isequal(VV.Ix{ii}, VV2.Ix{ii})
       isequal(VV.Iz{ii}, VV2.Iz{ii})
        %[G4 DG4] = hit_objective(VV2.ParticleArray.xv(:,ii));
        
    end 
    
    dGda_1 = VV.dGda(:);
    
    
    C = (VV.DG - VV2.DG) ./ VV.DG;
    figure()
    %imagesc(C)
    C( isnan(C) ) = 0;
	imagesc(C)
    
    C2 = (DG3 - VV.DG) ./ VV.DG;
    C2( isnan(C2) ) = 0;
    figure()
    imagesc(C2)
    
    D = (VV.dFdEx - VV2.dFdEx) ./ VV.dFdEx;
    D (isnan(D)) = 0;
    figure()
    imagesc(D)
    max(max(abs(D)))
    
    D = (VV.dFdEr - VV2.dFdEr) ./ VV.dFdEr;
    D (isnan(D)) = 0;
    figure()
    imagesc(D)
    max(max(abs(D)))
    
    D = (VV.dFdEz - VV2.dFdEz) ./ VV.dFdEz;
    D (isnan(D)) = 0;
    figure()
    imagesc(D)
    max(max(abs(D)))
    
    D = (VV.dFdV - VV2.dFdV) ./ VV.dFdV;
    D (isnan(D)) = 0;
    figure()
    imagesc(D)
    max(max(abs(D)))
    
    fprintf('Trajectory and Dual Trajectory Calculations done')



end