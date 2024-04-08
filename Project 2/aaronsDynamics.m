global m2km sec2hr sec2day R_E mu_E mu_S A_m AU w_dot P_phi epoch_0
    %conversions
m2km = (1/1000);
sec2hr = (1/3600);
sec2day = sec2hr/24;
    %other
R_E = 6378.1363; % km
R_S = 696e3; % km
mu_E = 3.986004415E5; % km^3/s^2
mu_S = 132712440017.987; % km^3/s^2
AU = 149597870.7; % km
A_m = m2km*m2km*0.01; % km^2/kg
w_dot = 7.29211585275553e-5; % rad/sec
phi_flux = 1357; % W/m^2 at 1 AU
c_sol = 2.99792458e8; % m/s %%%%%%% maybe change to just 3e8???
P_phi = AU^2*(phi_flux/c_sol)/m2km; %kg/(km*s^2)
epoch_0 = 2456296.25; % JD (Julian Day)

% function x_dot = ode_3bp(t,x)

% global mu_E mu_S P_phi A_m epoch_0 sec2day

%find sun position at current time (ECI frame)
[E_r, ~, ~] = Ephem(epoch_0+(t*sec2day),3,'EME2000');
r_S_E = -E_r';
r_S_E_mag = norm(r_S_E);

%r_dot = v
x_dot(1) = x(4);
x_dot(2) = x(5);
x_dot(3) = x(6);

%calcs for acceleration equation
r_sc_E = x(1:3);
r_sc_E_mag = norm(r_sc_E);
r_S_sc = r_S_E - r_sc_E';
r_S_sc_mag = norm(r_S_sc);
r_S_sc_unit = r_S_sc*(1/r_S_sc_mag);
a_Earth = -mu_E/r_sc_E_mag^3 * [x(1),x(2),x(3)];
a_SRP = -x(7)*(P_phi/r_S_sc_mag^2)*A_m*r_S_sc_unit; 
a_Sun = mu_S*( (r_S_sc/r_S_sc_mag^3) - (r_S_E/r_S_E_mag^3));

%acceration equations
x_dot(4) = a_Earth(1) + a_SRP(1) + a_Sun(1);
x_dot(5) = a_Earth(2) + a_SRP(2) + a_Sun(2);
x_dot(6) = a_Earth(3) + a_SRP(3) + a_Sun(3);
x_dot(7) = 0;

% r_Sun2Earth = -r_S_E
% r_SC2Sun = r_S_sc
% a_Earth = a_Earth
% a_Sun = a_Sun
% a_SRP = a_SRP

x_dot = x_dot';