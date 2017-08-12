classdef FEMProblem < handle
    
    properties
        poisson@PoissonFEM2D;
    end
    
    
    methods
        
        function obj = FEMProblem(N, faces, vertices )
            meshNodes = TriNodalMesh(N, faces, vertices);
            obj.poisson = PoissonFEM2D(meshNodes);
        end
        
        function solve(obj)
            
            
        end
        
    end
    
end