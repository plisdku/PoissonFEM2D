
checkClose = @(a,b) assert( norm(a-b) < 1e-6 ); %, sprintf('(%g) is not close to (%g)', a, b));

vertices = [0, 0; 2, 0; 0, 1];
faces = [1, 2, 3];

m = TriNodalMesh(4, faces, vertices);