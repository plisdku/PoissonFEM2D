function [faceNodes, edgeNodes, xy, edgePhysicalEntity, edgeElementaryEntity] = readMSH(fname)

fh = fopen(fname);

done = false;

while false == done
    
    line = fgetl(fh);
    token = strtok(line);
    
    switch token
        case -1
            done = true;
        case '$MeshFormat'
            blockIsDone = false;
            while false == blockIsDone
                line = fgetl(fh);
                
                if strcmpi(line, '$EndMeshFormat')
                    blockIsDone = true;
                else
                    %fprintf('Version string: %s\n', line);
                end
            end
            
        case '$PhysicalNames'
            line = fgetl(fh);
            numNames = str2num(line);
            %fprintf('There are %i names\n', numNames);
            
            for nn = 1:numNames
                line = fgetl(fh);
                [dat, count] = sscanf(line, '%f %f %s');
                dimension = dat(1); 
                physicalNumber = dat(2);
                physicalName = char(dat(4:end-1)');
                %fprintf('#%i %id (%s)\n', physicalNumber, dimension, physicalName);
            end
            fgetl(fh); % read the $EndPhysicalNames tag.  don't have to.
            
        case '$Nodes'
            line = fgetl(fh);
            numNodes = str2num(line);
            
            xy = zeros(numNodes, 2);
            
            for nn = 1:numNodes
                line = fgetl(fh);
                [dat, count] = sscanf(line, '%f %f %f %f');
                nodeId = dat(1);
                xy(nn,1) = dat(2);
                xy(nn,2) = dat(3);
            end
            
            %fprintf('Read %i nodes\n', numNodes);
            
        case '$Elements'
            line = fgetl(fh);
            numElements = str2num(line);
            
            % Don't know in advance how many edges and how many faces
            % so let's over-allocate and trim the arrays later.
            edgeNodes = zeros(numElements, 2);
            faceNodes = zeros(numElements, 3);
            
            edgeIds = zeros(numElements, 1);
            faceIds = zeros(numElements, 1);
            
            edgePhysicalEntity = zeros(numElements, 1);
            facePhysicalEntity = zeros(numElements, 1);
            edgeElementaryEntity = zeros(numElements, 1);
            faceElementaryEntity = zeros(numElements, 1);
            
            edgeIdx = 1;
            faceIdx = 1;
            
            for ee = 1:numElements
                line = fgetl(fh);
                
                [dat, count] = sscanf(line, '%f %f %f %f %f %f %f %f');
                if count == 7  % edge
                    edgeIds(edgeIdx) = dat(1);
                    edgePhysicalEntity(edgeIdx) = dat(4);
                    edgeElementaryEntity(edgeIdx) = dat(5);
                    edgeNodes(edgeIdx,:) = dat(6:7);
                    edgeIdx = edgeIdx + 1;
                elseif count == 8 % face
                    faceIds(faceIdx) = dat(1);
                    facePhysicalEntity(faceIdx) = dat(4);
                    faceElementaryEntity(faceIdx) = dat(5);
                    faceNodes(faceIdx,:) = dat(6:8);
                    faceIdx = faceIdx + 1;
                else
                    error('this did not happen');
                end
            end
            
            edgeNodes = edgeNodes(1:edgeIdx-1,:);
            edgeIds = edgeIds(1:edgeIdx-1,:);
            edgePhysicalEntity = edgePhysicalEntity(1:edgeIdx-1);
            edgeElementaryEntity = edgeElementaryEntity(1:edgeIdx-1);
            
            faceNodes = faceNodes(1:faceIdx-1,:);
            faceIds = faceIds(1:faceIdx-1,:);
            facePhysicalEntity = facePhysicalEntity(1:faceIdx-1);
            faceElementaryEntity = faceElementaryEntity(1:faceIdx-1);

        otherwise
    end
end