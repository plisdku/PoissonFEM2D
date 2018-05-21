figure(); hold on
       
[ix_x, ix_y, ix_z, ix_vx, ix_vy, ix_vz] =...
            get_Index3D(Nt);
        for ii = 1:size(xv,2)
            subplot(2,1,1)
            plot(ts,xv(ix_vx,ii),'b', 'LineWidth', 1)
            title('v_x')
            subplot(2,1,2)
            plot(ts,xv(ix_vy,ii),'b', 'LineWidth', 1)
            title('v_y')

            %plot(xv(ix_vx(1),ii),xv(ix_y(1),ii),'rx')
        end
        %plot(x_p,y_p,'bx')
        
        