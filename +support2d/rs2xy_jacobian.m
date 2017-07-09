function J = rs2xy_jacobian(rsTri)

[T, ~] = support2d.rs2xy_affineParameters(rsTri);

J = T;