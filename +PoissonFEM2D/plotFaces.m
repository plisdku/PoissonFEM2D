function plotFaces(nFigure, xyGeomNodes, hMesh, varargin)

        faceVerts = hMesh.faceVertices;

        if nargin == 3
            face_X = [xyGeomNodes(faceVerts(:,1),1)';xyGeomNodes(faceVerts(:,2),1)'; xyGeomNodes(faceVerts(:,3),1)'];
            face_Y = [xyGeomNodes(faceVerts(:,1),2)';xyGeomNodes(faceVerts(:,2),2)'; xyGeomNodes(faceVerts(:,3),2)'];


            figure(nFigure); clf    
            hold on
            patch(face_X, face_Y,'r')
        else
            iF = varargin{1};
            
            face_X_idxN = [xyGeomNodes(faceVerts(iF,1),1)';xyGeomNodes(faceVerts(iF,2),1)'; xyGeomNodes(faceVerts(iF,3),1)'];
            face_Y_idxN = [xyGeomNodes(faceVerts(iF,1),2)';xyGeomNodes(faceVerts(iF,2),2)'; xyGeomNodes(faceVerts(iF,3),2)'];
            figure(nFigure); clf
            
            patch(face_X_idxN, face_Y_idxN, 'r')
            
        end 
            
                
                
end