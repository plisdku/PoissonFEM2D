%% Plot function for Geometries

function plotElectrodePoints(varargin)

    N_elec = nargin;
    
    figure(914162120)
    clf
    hold on
    for i=1:2:N_elec
        Ex = varargin{i};
        Ey = varargin{i+1};
        Ex_ex = [Ex Ex(1,:)];
        Ey_ex = [Ey Ey(1,:)];
        plot(Ex_ex, Ey_ex,'-o')
    end
    xlim(xlim+[-10e-3 10e-3]);
    ylim(ylim+[-10e-3 10e-3]);
    
end