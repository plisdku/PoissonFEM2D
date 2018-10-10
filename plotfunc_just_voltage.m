function plotfunc_just_voltage(xHistory, fHistory, DfHistory, max_x,alpha)
    iter = length(fHistory);
    figure(10100)
    subplot(3,1,1)
    semilogy(1:iter, fHistory)
    xlabel('iteration')
    ylabel('F')
%     subplot(3,1,2)
%     plot(1:iter, xHistory(1:16,1:iter))
%     xlabel('iteration')
%     ylabel('change in geom. parameters')
%     N_idx = N_vec*4+4;
%     idx_r = 1:1:N_idx(1);
%     idx_c = N_idx(1)+1:1:(N_idx(1)+2);
%     idx_a = (N_idx(1)+3):1:(2*N_idx(1)+2);
%     
%     idx_radii = idx_r;
%     idx_center = idx_c; 
%     idx_angles = idx_a; 
%     
%     for i=1:(length(N_idx)-1)
%         
%         idx_r = idx_r + N_idx(i)*2+2;
%         idx_c = idx_c + N_idx(i)*2+2;
%         idx_a = idx_a + N_idx(i)*2+2;
%         
%         idx_radii = [idx_radii, idx_r];
%         idx_center = [idx_center, idx_c];
%         idx_angles = [idx_angles, idx_a];
%         
%         
%     end
% 
%     idx_voltage = (max_p+1):1:(max_p+N_voltage);
%     
%     
%     subplot(4,2,3)
%     plot(1:iter, xHistory(idx_radii,1:iter))
%     xlabel('iteration')
%     ylabel('change in radii')
%     
%     subplot(4,2,4)
%     plot(1:iter, xHistory(idx_angles,1:iter))
%     xlabel('iteration')
%     ylabel('change in angles')
%     
%     subplot(4,2,5)
%     plot(1:iter, xHistory(idx_center,1:iter))
%     xlabel('iteration')
%     ylabel('change in center point')
%     
%     subplot(4,2,6)
%     plot(1:iter, xHistory(idx_voltage,1:iter))
%     xlabel('iteration')
%     ylabel('change in voltage')
%     
    subplot(3,1,2)
    plot(1:iter, xHistory(:,1:iter))

    subplot(3,1,3)
    plot(1:iter, alpha)
    xlabel('iteration')
    ylabel('Stepsize')
end