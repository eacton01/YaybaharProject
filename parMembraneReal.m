function[par] = parMembraneReal
%parameters for membrane

%% PHYSICAL PARAMETERS %%%%%%%%%%
% membrane parameters
par.Rad = 0.25;     % Radius
par.h = 0.2e-3;     % Thickness
par.sig = 0.27;     % surface density
par.E = 3.5e9;      % Young's modulus
par.v = 0.2;        % Poisson ratio
par.dp1 = 1.25;     % Freq ind damping
par.dp3 = 5e-4;     % Freq dep damping
par.Tau = 780;      % Tension
par.c = sqrt(par.Tau/par.sig);  % wave speed
par.D = (par.E*par.h^3) / (12*(1-par.v^2)); % Bending stiffness
par.area = pi*par.Rad^2;    % Area
par.mass = par.area*par.sig;    % Mass


