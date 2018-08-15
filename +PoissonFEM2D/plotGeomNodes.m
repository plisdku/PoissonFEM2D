function plotGeomNodes(nFigure, xyGeomNodes, varargin)

% nFigure, geometry, varargin of the form femProblem and a list of indices
    
    figure(nFigure); clf
    hold on 
    plot(xyGeomNodes(:,1),xyGeomNodes(:,2),'o')
    if nargin > 2
       n_idx = nargin - 2;
       
       for i = 1:n_idx
           
           idxN = varargin{i};
           plot(xyGeomNodes(idxN,1),xyGeomNodes(idxN,2),'x', 'Color', 'r')
            
       end
    

end