function[par] = parSlinkyReal
% parameters for slinky

%% PHYSICAL PARAMETERS %%%%%%%%%%
par.L = 2.1;                % length (m)
par.varrho = 0.03;          % mass per unit length
par.Tau = 0;                % tension   
par.c = sqrt(par.Tau/par.varrho);   % wave velocity (m/s)    
par.kappa = 11.5e-3;        % stiffness
par.b1 = 5e-1;              % freq independant damping coef
par.b2 = 12e-7;             % freq dependant damping coef
par.GL = 0;                 % boundary admittance on left side of string


