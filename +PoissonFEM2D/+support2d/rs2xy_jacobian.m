function J = rs2xy_jacobian(rsTri)
import PoissonFEM2D.*
[T, ~] = support2d.rs2xy_affineParameters(rsTri);

J = T;