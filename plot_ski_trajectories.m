 figure(1510+i)
    clf
    %imagesc(x_grid,r_grid,U(:,:,2)')
    axis xy
    axis square
    hold on 
    femProblem.poi.tnMesh.plotMesh();
    
    [xxv, yyv] = ndgrid(x_grid, r_grid);
    %U_matrix = U(:,:,2)';

    for ii = 1:size(xv_matrix,2)
        plot(xv_matrix(ix_x,ii),xv_matrix(ix_y,ii),'Color','w', 'LineWidth', 1)
        plot(xv_matrix(ix_x(1),ii),xv_matrix(ix_y(1),ii),'rx')
    end
    colormap orangecrush
    
    plot(xxv(:),yyv(:),'wo'),
%%

        
figure(15100+i)
        clf
  imagesc(x_grid,r_grid,U(:,:,2)')
  %imagesc(U(:,:,2)')
  axis xy
  axis square
        hold on 
        %femProblem.poi.tnMesh.plotMesh();

for ii = 1:size(xv_matrix,2)
            plot(xv_matrix(ix_x,ii),xv_matrix(ix_y,ii),'Color','w', 'LineWidth', 1)
            plot(xv_matrix(ix_x(1),ii),xv_matrix(ix_y(1),ii),'rx')
            disp(ii)
end
        colormap orangecrush
        
        
        %plot(xxv(:),yyv(:),'wo'),


    [xxv, yyv] = ndgrid(x_grid, r_grid);

    %plot(xxv(:), yyv(:), 'g.')
    
    
    
%%
    
    
    for ii = 1:size(E_r,1)
        
        figure(71)
        plot(y_grid(1:111),E_r(ii,length(y_grid)+1:(length(y_grid) + 111)))
        title('E_y field')
        xlabel('y')
        ylabel('E_y')
        title(sprintf('Slice of E_y at x = %0.4f', x_grid(ii)))
        disp(ii)
        disp('done')
        pause
        
    end