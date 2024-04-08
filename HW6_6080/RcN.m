function [RcN] = RcN(t)
    
    rmag_gmo = 20424.2; % [km] radius of gmo orbit
    r_mars = 3396.19; %[km] radius of mars
    h_lmo = 400; %[km] altitude of lmo above mars' surface
    rmag_lmo = r_mars + h_lmo; % [km] radius of lmo orbit in hill fram
    mu = 42828.3; % [km^3/s^2] mars gravitational co
    n1 = [1;0;0];
    n2 = [0;1;0];
    n3 = [0;0;1];
    % Initial orbit fram orientation angles
    OMEGA_lmo = 20*pi/180; % [rad]
    i_lmo = 30*pi/180; % [rad]
    theta_lmo = 60*pi/180; % [rad]
%     thetadot_lmo = 0.000884797; % [rad/s]
    thetadot_lmo = sqrt(mu/rmag_lmo^3); % [rad/s]

    OMEGA_gmo = 0; % [deg]
    i_gmo = 0; % [deg]
    theta_gmo = 250*pi/180; % [rad]
%     thetadot_gmo = 0.0000709003; % [rad/s]
    thetadot_gmo = sqrt(mu/rmag_gmo^3); % [rad/s]

    Xlmo = [rmag_lmo; OMEGA_lmo; i_lmo; theta_lmo+thetadot_lmo*t];
    r_lmo = rrdot(Xlmo);
    Xgmo = [rmag_gmo; OMEGA_gmo; i_gmo; theta_gmo+thetadot_gmo*t];
    r_gmo = rrdot(Xgmo);
    
    deltar = r_gmo-r_lmo;
    
    r1 = -deltar/norm(deltar);
    r2 = cross(deltar, n3)/norm(cross(deltar,n3));
    r3 = cross(r1,r2);
    RcN = [r1'; r2'; r3'];
end

