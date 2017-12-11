%% Make a test mesh.

vertices = [0,0; 1,0; 0,1];
faces = [1,2,3];

% Node orders
N_field = 4;
N_geom = 3;
N_quad = N_field;

%N_field=2;
%N_geom=2;
%N_quad=2;

lng = LinearNodalGeometry(faces, vertices, N_geom);
xyNodes = lng.getNodeCoordinates();
tnMesh = TriNodalMesh(faces, xyNodes, N_field, N_geom, N_quad);

poi = PoissonFEM2D(tnMesh);

delta = 1e-8;


%% Element potential matrix

A = poi.getElementPotentialMatrix(1);
DA = poi.getElementPotentialMatrixSensitivity(1);

iGeomGlobal = tnMesh.hGeomNodes.getFaceNodes(1);

% Iterate over geometry nodes, perturb, test
for mm = 1:tnMesh.hGeomNodes.getNumNodes()
    for dirIdx = 1:2
        A2 = PoissonFEM2D(tnMesh.perturbed(iGeomGlobal(mm), dirIdx, delta)).getElementPotentialMatrix(1);
        DA_meas = (A2-A)/delta;
        DA_calc = DA(:,:,dirIdx,mm);
        
        [same, relErr, normDiff, normExact] = compareNorms(DA_calc, DA_meas);
        if ~same
            fprintf('DJ error %0.4e out of %0.4e ***\n', normDiff, normExact)
        end
    end
end

%% Element charge matrix

A = poi.getElementChargeMatrix(1);
DA = poi.getElementChargeMatrixSensitivity(1);

iGeomGlobal = tnMesh.hGeomNodes.getFaceNodes(1);

% Iterate over geometry nodes, perturb, test
for mm = 1:tnMesh.hGeomNodes.getNumNodes()
    for dirIdx = 1:2
        A2 = PoissonFEM2D(tnMesh.perturbed(iGeomGlobal(mm), dirIdx, delta)).getElementChargeMatrix(1);
        DA_meas = (A2-A)/delta;
        DA_calc = DA(:,:,dirIdx,mm);
        
        [same, relErr, normDiff, normExact] = compareNorms(DA_calc, DA_meas);
        if ~same
            fprintf('DJ error %0.4e out of %0.4e ***\n', normDiff, normExact)
        end
    end
end

%% Element Neumann matrix

iEdge = 1;
orientation = 1;
A = poi.getElementNeumannMatrix(iEdge, orientation);
DA = poi.getElementNeumannMatrixSensitivity(iEdge, orientation);

iGeomGlobal = tnMesh.hGeomNodes.getEdgeNodes(iEdge, orientation);

% Iterate over geometry nodes, perturb, test
for mm = 1:length(iGeomGlobal)
    for dirIdx = 1:2
        A2 = PoissonFEM2D(tnMesh.perturbed(iGeomGlobal(mm), dirIdx, delta)).getElementNeumannMatrix(iEdge, orientation);
        DA_meas = (A2-A)/delta;
        DA_calc = DA(:,:,dirIdx,mm);
        
        [same, relErr, normDiff, normExact] = compareNorms(DA_calc, DA_meas, 1e-6, 1e-6); % tolerance fail lol
        if ~same
            fprintf('DJ error %0.4e out of %0.4e ***\n', normDiff, normExact)
        end
    end
end



%% Neumann matrix

A = poi.getNeumannMatrix();
DA = poi.getNeumannMatrixSensitivity();

for mm = 1:tnMesh.hGeomNodes.getNumNodes()
    for dirIdx = 1:2
        A2 = PoissonFEM2D(tnMesh.perturbed(mm, dirIdx, delta)).getNeumannMatrix();
        DA_meas = (A2-A)/delta;
        DA_calc = DA{dirIdx,mm};
        
        [same, relErr, normDiff, normExact] = compareNorms(DA_calc, DA_meas, 1e-6, 1e-6); % tolerance fail lol
        if ~same
            fprintf('DJ error %0.4e out of %0.4e ***\n', normDiff, normExact);
        end
    end
end

%% Rhs Matrix

A = poi.getRhsMatrix();
DA = poi.getRhsMatrixSensitivity();

for mm = 1:tnMesh.hGeomNodes.getNumNodes()
    for dirIdx = 1:2
        A2 = PoissonFEM2D(tnMesh.perturbed(mm, dirIdx, delta)).getRhsMatrix();
        DA_meas = (A2-A)/delta;
        DA_calc = DA{dirIdx,mm};
        
        [same, relErr, normDiff, normExact] = compareNorms(DA_calc, DA_meas);
        if ~same
            fprintf('DJ error %0.4e out of %0.4e ***\n', normDiff, normExact);
        end
    end
end

%% System matrix

A = poi.getSystemMatrix();
DA = poi.getSystemMatrixSensitivity();

for mm = 1:tnMesh.hGeomNodes.getNumNodes()
    for dirIdx = 1:2
        A2 = PoissonFEM2D(tnMesh.perturbed(mm, dirIdx, delta)).getSystemMatrix();
        DA_meas = (A2-A)/delta;
        DA_calc = DA{dirIdx,mm};
        
        [same, relErr, normDiff, normExact] = compareNorms(DA_calc, DA_meas, 1e-6, 1e-6); % tolerance fail lol
        if ~same
            fprintf('DJ error %0.4e out of %0.4e ***\n', normDiff, normExact);
        end
    end
end

%% Function evaluation on field nodes

func_x = @(x,y) x; % Try lots of things
[f, dfdxg, dfdyg] = poi.evaluateOnNodes(func_x);

for mm = 1:tnMesh.hGeomNodes.getNumNodes()
    for dirIdx = 1:2
        f2 = PoissonFEM2D(tnMesh.perturbed(mm, dirIdx, delta)).evaluateOnNodes(func_x);
        df_meas = (f2-f)/delta;
        
        if dirIdx == 1
            df_calc = dfdxg(:,mm);
        else
            df_calc = dfdyg(:,mm);
        end
        
        [same, relErr, normDiff, normExact] = compareNorms(df_calc, df_meas);
        if ~same
            fprintf('df error %0.4e out of %0.4e ***\n', normDiff, normExact);
        end
    end
end






