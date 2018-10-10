function [particles, hit_objective] = SetupParticles_chrom()

    elementary_charge   = 1.60217662e-19;
    ion_mass = 5.1477e-26;
    electron_mass       = 1.6605e-27;

    E_center = 30*1e3*elementary_charge;
    v_center = sqrt(2*E_center / ion_mass);
    
    Ne = 10;
    v_center_vec = linspace(v_center*0.98,v_center*1.02,Ne);

    Nr = 10;

    y_vec = linspace(-0.2e-3,0.2e-3,Nr);

    x_pos = -50e-3;
    
    [Y_pos, V_center] = ndgrid(y_vec,v_center_vec);

    x_pos = -50e-3;
    
    xv0 = [x_pos*ones(1,Nr*Ne);Y_pos(:)';zeros(1,Nr*Ne); V_center(:)'; zeros(1,Nr*Ne); zeros(1,Nr*Ne)];
    
    Nt = 1500;
    end_time = (120e-3 ./ (v_center))*2.5;
    end_time = end_time - 0.6*end_time;

    n_charges = -1;
    n_masses = ion_mass/electron_mass;


    particles(size(xv0,2)) = Particle(xv0(4:6,size(xv0,2)),'Velocity',n_masses,n_charges,xv0(1:3,size(xv0,2)), 0, end_time, Nt);
    for zz = 1:size(xv0,2)
         particles(zz) = Particle(xv0(4:6,zz),'Velocity',n_masses,n_charges,xv0(1:3,zz), 0, end_time, Nt);
    end
    

    x_p = 50e-3;%39.5e-3;
    y_p  = 0;
    z_p  = 0;
    vx_p = 1e-3;
    vy_p = 0;
    vz_p = 0;

    obj_weights = [1, 1, 0, 0, 0, 0];
    hit_objective = @(x_v) hitObjective3D_wrap(...
        x_v, x_p, y_p, z_p, vx_p, vy_p, vz_p, obj_weights);
end