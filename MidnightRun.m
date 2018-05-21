   a = zeros(2,size(E_r,1));
for ii = 1:size(E_r,1)
        
        figure(71)
        plot(y_grid(1:111),E_r(ii,length(y_grid)+1:(length(y_grid) + 111)))
        title('E_y field')
        xlabel('y')
        ylabel('E_y')
        title(sprintf('Slice of E_y at x = %0.4f', x_grid(ii)))
        disp(ii)
        disp('done')
        %pause
        [a(:,ii)]=polyfit(y_grid(1:111),E_r(ii,length(y_grid)+1:(length(y_grid)+ 111)),1);

        
end
    
    

%%
x_intercepts = -a(2,:)./a(1,:);
figure(77)
plot(x_grid,x_intercepts,'LineWidth',2)
set(gca,'FontSize', 18)
ylabel('y intercept')
xlabel('x - pos')

mean_intercept = mean(x_intercepts);
std_intercept = sqrt(var(x_intercepts));
title(sprintf('x-Intercepts (zero crossings) at different pos., mean = %0.2e, std. dev. = %0.2e',mean_intercept, std_intercept))



%%


        uu(:,:,1) = femProblem.uCartesian;  
        uu(:,:,2) = femProblem.uCartesian;  
        U = [uu(:,end:-1:2,:) uu(:,:,:)];

        r_grid = [-y_grid(end:-1:2) y_grid];
        x_grid = linspace(measBox(1),measBox(3),measNxy(1));
        y_grid = linspace(measBox(2),measBox(4),measNxy(2));
        z_grid = [-2, 0];
        d_x = x_grid(2) - x_grid(1);
        d_r = y_grid(2) - y_grid(1);
        d_z = z_grid(2) - z_grid(1);


        E_x         = -centeredDiff(U, 1) / d_x;
        E_r4        = -centeredDiff(U, 2) / d_r;
        E_z         = -centeredDiff(U, 3) / d_z;
        
%%
   a2 = zeros(2,size(E_r4,1));
for ii = 1:size(E_r4,1)
        
        figure(72)
        plot(y_grid(2:round(0.2*length(y_grid))), E_r4(ii,(length(y_grid)+1):( length(y_grid)  -1 + round(0.2*length(y_grid)) )))
        title('E_y field')
        xlabel('y')
        ylabel('E_y')
        title(sprintf('Slice of E_y at x = %0.4f', x_grid(ii)))
        disp(ii)
        disp('done')
        %pause
        [a2(:,ii)]=polyfit(y_grid(2:end),E_r4(ii,length(y_grid)+1:size(E_r4,2)),1);

        
end
    

x_intercepts = -a2(2,:)./a2(1,:);
figure(77)
plot(x_grid,x_intercepts,'LineWidth',2)
set(gca,'FontSize', 18)
ylabel('y intercept')
xlabel('x - pos')

mean_intercept = mean(x_intercepts);
std_intercept = sqrt(var(x_intercepts));
title(sprintf('x-Intercepts (zero crossings) at different pos., mean = %0.2e, std. dev. = %0.2e, N=4 low res',mean_intercept, std_intercept))
%%
a3 = zeros(2,Nt);
for ii = 1:Nt
        
        figure(72)
        plot(1:(size(xv_matrix,2)-1),xv_matrix(ix_y(ii),2:end))
        %title('E_y field')
        xlabel('n_Particle')
        ylabel('y_pos')
        title(sprintf('y_pos of an particle at Nt = %i', ii))
        disp(ii)
        disp('done')
        %pause
        [a3(:,ii)]=polyfit(1:(size(xv_matrix,2)-1),xv_matrix(ix_y(ii),2:end),1);

        
end 
    
%%
x_intercepts = -a3(2,:)./a3(1,:);
figure(77)
plot(1:Nt,x_intercepts,'LineWidth',2)
set(gca,'FontSize', 18)
ylabel('y intercept')
xlabel('x - pos')

mean_intercept = mean(x_intercepts);
std_intercept = sqrt(var(x_intercepts));
title(sprintf('r-Intercepts (zero crossings, interpolated) at different Time Steps, mean = %0.2e, std. dev. = %0.2e, N=4 low res',mean_intercept, std_intercept))


%%

 [ix_x, ix_y, ix_z, iv_x, iv_y, iv_z] =...
            get_Index3D(Nt);

a4 = zeros(2,Nt);
for ii = 1:Nt
        
        figure(72)
        plot(xv_matrix(ix_y(ii),1:end),xv_matrix(iv_y(ii),1:end))
        %title('E_y field')
        xlabel('n_Particle')
        ylabel('y_pos')
        title(sprintf('y_pos of an particle at Nt = %i', ii))
        disp(ii)
        disp('done')
        %pause
        [a4(:,ii)]=polyfit(xv_matrix(ix_y(ii),1:end),xv_matrix(iv_y(ii),1:end),1);

        
end 
    
%%
x_intercepts = -a4(2,:)./a4(1,:);
figure(77)
clf
hold on
subplot(2,1,1)
plot(xv_matrix(ix_x,1),x_intercepts,'LineWidth',2)
grid on 
%set(gca,'FontSize', 18)
ylabel('r intercept')
xlabel('x pos.')

%mean_intercept = mean((x_intercepts));
%std_intercept = sqrt(var(isx_intercepts));
title(sprintf('r-Intercepts (zero crossings, interpolated) of v_r vs. r at different Time Steps N=5 low res'))
ylim([-0.5e-7 1.5e-7])
%xlim([xv_matrix(ix_x(1),1) xv_matrix(ix_x(end),1)])
xlim([-0.06 0.06])
 subplot(2,1,2)
 plot(xv_matrix(ix_x,1),a4(1,:), 'LineWidth', 2)
 xlabel('x pos.')
 ylabel('v_r/r slopes')
title('Slopes')
ylim([-0.3e8 0.3e8])
xlim([-0.06 0.06])
xlim([xv_matrix(ix_x(1),1) xv_matrix(ix_x(end),1)])

