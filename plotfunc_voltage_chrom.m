function plotfunc_voltage_chrom(xHistory, fHistory, DfHistory)
    iter = length(fHistory);
    figure(20170)
    subplot(2,1,1)
    semilogy(1:iter, fHistory)
    xlabel('iteration')
    ylabel('log F')
%     subplot(3,1,2)
%     plot(1:iter, xHistory(1:16,1:iter))
%     xlabel('iteration')
%     ylabel('change in geom. parameters')
    subplot(2,1,2)
    plot(1:iter, xHistory(1:3,1:iter))
    xlabel('iteration')
    ylabel('change in voltage parameters')
end