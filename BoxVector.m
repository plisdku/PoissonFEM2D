%% Function for arbitrary number of points per line for a box
% Input: 4 Vertex Points, vector with number of points for ech of the four
% lines
% Output: Vector with points in the style required by inhouse 2D FEM Code
% (starting in the lower left corner, counter-clockwise until lower left
% corner is reached again) 

function varargout = BoxVector(Vertices, NlinePoints, varargin) 

    assert(size(Vertices,1) == size(NlinePoints,2), 'Number of Line Length elements must be the same as the number of vertices')

    Box_x = linspace(Vertices(1,1),Vertices(2,1),NlinePoints(1));
    Box_y = linspace(Vertices(1,2),Vertices(2,2),NlinePoints(1));

        for ii = 2:size(NlinePoints,2)-1
            x = linspace(Vertices(ii,1),Vertices(ii+1,1),...
                NlinePoints(ii));
            Box_x = ...
                [Box_x x(2:end)];
            y = linspace(Vertices(ii,2),Vertices(ii+1,2),...
                NlinePoints(ii));
            Box_y = ...
                [Box_y y(2:end)];

        end 
    x = linspace(Vertices(ii+1,1),Vertices(1,1),NlinePoints(ii+1));
    y = linspace(Vertices(ii+1,2),Vertices(1,2),NlinePoints(ii+1));
    Box_x = [Box_x x(2:(end-1))];
    Box_y = [Box_y y(2:(end-1))];
    
    varargout{1} = Box_x;
    varargout{2} = Box_y;
    
    if nargin > 2
        
        s_func = varargin{1};
        assert(isa(s_func,'function_handle'), 'Must use a function handle to define a function for the s vector')
       
        s_vec = s_func(NlinePoints);
        
        varargout{3} = s_vec;
    end

end