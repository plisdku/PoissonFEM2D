[particles3, hit_objective3] = SetupParticles_Array();
[particles, hit_objective, particles2, hit_objective2] = SetupParticles_test();
[VV] = ElectronSetup_test(femProblem.uCartesian, measBox, measNxy, particles, particles3, hit_objective, hit_objective3); 
    