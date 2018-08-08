%% Collective Tests for PoissonFEM2D Tests

import PoissonFEM2D.*

%% I
try 
    
    PoissonFEM2D.testFEMFunctionals
    
    fprintf('testFEMFunctionals has run \n')
    
catch 
    
    fprintf('testFEMFunctionals failed \n')
    
    % testFEMFunctionals uses a reference to the property hNodes (line 15)
    % this is not a property of TriNodalMesh (anymore) and might have 
    % moved into a different object. 
    % Tests could not run through, there might be more errors. 
    
end 

%% II

try 
    
    PoissonFEM2D.testFEMSensitivity 
    
    fprintf('testFEMSensitivity has run \n')
    
catch 
    
    fprintf('testFEMSensitivity failed \n')
    
    % testFEMFunctionals uses an objective function to call the solve
    % method (line 60)
    % the solve method does not take any input parameters (anymore) This
    % feature has probably moved into a different function (see comments 
    % in the class) 
    % Tests could not run through, there might be more errors. 
    
end 

%% III

try 
    
    PoissonFEM2D.testJacobiEvaluator 
    
    fprintf('testJacobiEvaluator  has run \n')
    
catch 
    
    fprintf('testJacobiEvaluator failed \n')
    
    % Is expected to work 
    
end 

%% IV 

try 
    
    PoissonFEM2D.testPoissonFEM2DSensitivities
    
    fprintf('testPoissonFEM2DSensitivities has run \n')
    
catch 
    
    fprintf('testPoissonFEM2DSensitivities  failed \n')
    
    % Is expected to work 
    
end 

%% V 

try 
    
    PoissonFEM2D.testTriNodalMesh
    
    fprintf('testTriNodalMesh has run \n')
    % While this test runs through, there are errors in the test starting
    % in line 239. Moreover, there is a note that several other tests need
    % to be rewritten. 
    
catch 
    
    fprintf('testTriNodalMesh  failed \n')
    
    % See above 
    
end 

%% VI 

try 
    
    PoissonFEM2D.testTriNodalMeshSensitivities
    
    fprintf('testTriNodalMeshSensitivites has run \n')
    % While this test runs through, this test creates an error in the 
    % 'Face gradient matrix sensitivity' test on all machines except 
    % Paul's laptop 
    % *** DDy error 1.2515e-07 out of 1.2515e-07 ***
    
catch 
    
    fprintf('testTriNodalMeshSensitivites  failed \n')
    
    % See above 
    
end 
