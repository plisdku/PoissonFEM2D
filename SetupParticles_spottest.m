function [particles, hit_objective] = SetupParticles_spottest()

    elementary_charge   = 1.60217662e-19;
    ion_mass = 5.1477e-26;
    electron_mass       = 1.6605e-27;

    E_center = 30*1e3*elementary_charge;
    v_center = sqrt(2*E_center / ion_mass);

    Nr = 100;
    %angle = linspace(-2.5e-3,2.5e-3,Nr);

    y_vec = linspace(0e-6,5e-6,Nr);
    x_pos = -173.5e-3;
    
    velocities = v_center;
    
    %[xx, yy, angles, vv] = ndgrid(x_pos, y_vec, angle, velocities);

    %angles = angles(:);
    %vv = vv(:);
    
    xv0 = [x_pos*ones(1,Nr);y_vec ;zeros(1,Nr); v_center*ones(1,Nr); zeros(1,Nr); zeros(1,Nr)];
    
    Nt = 1500;
    end_time = (120e-3 ./ (v_center))*2.5;
    end_time = end_time - 0.2*end_time;

    n_charges = -1;
    n_masses = ion_mass/electron_mass;



    particles(size(xv0,2)) = Particle(xv0(4:6,size(xv0,2)),'Velocity',n_masses,n_charges,xv0(1:3,size(xv0,2)), 0, end_time, Nt);
    for zz = 1:size(xv0,2)
         particles(zz) = Particle(xv0(4:6,zz),'Velocity',n_masses,n_charges,xv0(1:3,zz), 0, end_time, Nt);
    end
    
    %%
     x_p = 52.5e-3;
    y_p  = 0;
    z_p  = 0;
    vx_p = 1e-3;
    vy_p = 0;
    vz_p = 0;

    obj_weights = [1, 1, 0, 0, 0, 0];
    hit_objective = @(x_v) hitObjective3D_wrap(...
        x_v, x_p, y_p, z_p, vx_p, vy_p, vz_p, obj_weights);
end