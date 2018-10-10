function [Ex, Ey, ExSpline, EySpline]= createBoxSpline(Ncontrol, Ex_func, Ey_func, p, Nspine)

% Function to create a Spline representation for a 4-edge geometry element.
% Input are the number of control points per side (non-exclusive, verticies
% are counted twice, a function to describe the x-coordinates and the
% y-coordinates, both assume a counter-clockwise direction. The function
% gets evalued with the vector p. The number of returned points on the
% spine is Nspine. 
% The function returns an x-Vector and any-Vector with all the control 
% points as well as each with all the interpolated spine points. 


    Ex = Ex_func(p);
    Ey = Ey_func(p);
    
    

    xx1 = linspace(Ex(1), Ex(Ncontrol(1)), Nspine);
    yy1 = spline(Ex(1:(Ncontrol(1))), Ey(1:(Ncontrol(1))), xx1);
    
    yy2 = linspace(Ey(Ncontrol(1)), Ey(Ncontrol(1)+Ncontrol(2)-1), ...
        Nspine);
    xx2 = spline(Ey(Ncontrol(1):(Ncontrol(1)+Ncontrol(2)-1)), ...
        Ex(Ncontrol(1):(Ncontrol(1)+Ncontrol(2)-1)), yy2);
    
    xx3 = linspace(Ex((Ncontrol(1)+Ncontrol(2)-1)),...
        Ex((Ncontrol(1)+Ncontrol(2)+Ncontrol(3)-2)), Nspine);
    yy3 = spline(Ex((Ncontrol(1)+Ncontrol(2)-1):...
        (Ncontrol(1)+Ncontrol(2)+Ncontrol(3)-2)),...
        Ey((Ncontrol(1)+Ncontrol(2)-1):...
        (Ncontrol(1)+Ncontrol(2)+Ncontrol(3)-2)), xx3);
    
    yy4 = linspace(Ey((Ncontrol(1)+Ncontrol(2)+Ncontrol(3)-2)), ...
        Ey(1), Nspine);
    xx4 = spline([Ey((Ncontrol(1)+Ncontrol(2)+Ncontrol(3)-2):...
        (Ncontrol(1)+Ncontrol(2)+Ncontrol(3)+Ncontrol(4)-4))...
        Ey(1)], [Ex((Ncontrol(1)+Ncontrol(2)+Ncontrol(3)-2):...
        (Ncontrol(1)+Ncontrol(2)+Ncontrol(3)+Ncontrol(4)-4)) ...
        Ex(1)], yy4);
    
    
   
% 
%     figure(566); clf
%     hold on
% 
%     plot(Ex, Ey, 'ro')
%     plot(xx1,yy1,'r-')
%     plot(xx2,yy2,'r-')
%     plot(xx3,yy3,'r-')
%     plot(xx4,yy4,'r-')
    
    
    ExSpline = [xx1(1:(end-1)) xx2(1:(end-1)) xx3(1:(end-1)) ...
        xx4(1:(end-1))];
    EySpline = [yy1(1:(end-1)) yy2(1:(end-1)) yy3(1:(end-1)) ...
        yy4(1:(end-1))];
    
%     figure(567)
%     plot(ExSpline,EySpline,'rx')
    
    
    
end 