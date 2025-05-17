%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Yaybahar with Sound Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;

% NUMERICAL TIME CONSTANTS AND PARAMETER LOADING %%%%%%%%%%%%%%%%%%%%%%%%
parString = parStringReal;            % Load string constants
parSlinky = parSlinkyReal;            % Load slinky constants
parMembrane = parMembraneReal;        % Load membrane constants

M_mem = 55;             % Number of m-modes for membrane (number of n-modes = M+1, total modes = N*(M+1))
Fs = 44100;             % sampling frequency
T = 1/Fs;               % temporal resolution
dur = 6;                % duration
Ns = floor(dur*Fs);     % number of samples
t = (0:Ns-1)*T;         % times vector

% EXCITATION POSITIONS AND SIGNALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% relative input position along slinky 
reEx_slinky = 0.23;

% relative coupling, excitation & pickup positions on string
reCoup_str = 0.015;
reEx_str = 0.4;
rePick1_str = 0.5;
rePick2_str = 0.6;

% relative coupling, excitation, and pickup positions on membrane radius
reCoup_mem = 0.3;
reEx_mem = reCoup_mem;
rePick1_mem = reCoup_mem + 0.03;
phiPick1 = -pi/30;      % angle of pickup position 1
phiPick2 = pi/20;       % angle of pickup position 2
rePick2_mem = reCoup_mem - 0.03;
if rePick2_mem < 0.03
    rePick2_mem = 0.03;     % prevent pickup point from being too near to centre
end

% input signals
inpSlinky = zeros(1,Ns);
hannLenSlinky = 7;
inpSlinky(1:hannLenSlinky) = 1e7*hanning(hannLenSlinky);

inpString = zeros(1,Ns);
hannLenString = 10;
%inpString(1:hannLenString) = 1e3*hanning(hannLenString);

inpMem = zeros(1,Ns);
hannLenMem = 40;
%inpMem(1:hannLenMem) = 1e1*hanning(hannLenMem);

% output signals
outp1 = zeros(1,Ns);
outp2 = zeros(1,Ns);
outp3 = zeros(1,Ns);
outp4 = zeros(1,Ns);

% PHYSICAL PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% slinky (u) parameters 
L_u = parSlinky.L;              % length
varrho_u = parSlinky.varrho;    % mass per unit length
Tau_u = parSlinky.Tau;          % tension (0)  
c_u = parSlinky.c;              % wave velocity (0)   
kappa_u = parSlinky.kappa;      % stiffness
b1_u = parSlinky.b1;            % damping coef
b2_u = parSlinky.b2;            % damping coef
GL_u = parSlinky.GL;            % N/A

% string parameters
L_str = parString.L;                % string length
rhoA_str = parString.rhoA;          % string mass per unit length
Tau_str = parString.T;              % string tension
EI_str = parString.EI;              % string stiffness
sig0_str = parString.eta0;          % string damping constant
sig2_str = parString.eta2;          % string damping constant
m_str = 0.5*rhoA_str*L_str;         % modal mass, kg
c_str = sqrt(Tau_str/(rhoA_str));   % Wave velocity, m/s

% membrane parameters
Rad_mem = parMembrane.Rad;      % Radius
h_mem = parMembrane.h;          % Thickness
sig_mem = parMembrane.sig;      % surface density
E_mem = parMembrane.E;          % Young's modulus
v_mem = parMembrane.v;          % Poisson ratio
dp1_mem = parMembrane.dp1;      % Freq ind damping
dp3_mem = parMembrane.dp3;      % Freq dep damping
Tau_mem = parMembrane.Tau;      % Tension
c_mem = parMembrane.c;          % wave speed
D_mem = parMembrane.D;          % Bending stiffness
area_mem = parMembrane.area;    % Area
mass_mem = parMembrane.mass;    % Mass

% SLINKY FD PRE-CALCUALTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[A_u,BL_u,BR_u,lam_u,mu_u,nu1_u,nu2_u,M_u,X_u] = string2FD(L_u,varrho_u,Tau_u,kappa_u,b1_u,b2_u,T,GL_u,0);
Np_u = M_u + 1;                     % number of nodes along slinky length
ke = round(reEx_slinky*M_u) + 1;    % excitation node
% initiation of variables & constants 
U0 = zeros(Np_u,1);                 % slinky displacement at time (n+1)
U1 = zeros(Np_u,1);                 % slinky displacement at time (n)    
U2 = zeros(Np_u,1);                 % slinky displacement at time (n-1)
F1 = zeros(Np_u,1);                 % force per unit slinky length at time (n)


% STRING MODAL PRE-CALCUALTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FINDING NUMBER OF MODES
Binh_str = (EI_str*pi^2)/(Tau_str*L_str^2);     % Inharmonicity factor
LAM_str = 0.9;                                  % Lambda
om1_str = (pi/L_str)*sqrt(Tau_str/(rhoA_str));  % Fundamental freq
N_str = 2000;                                   % Number of points on string

% CALCULATE NUMBER OF STRING MODES
eqC1 = Binh_str*om1_str^2 - (1/16)*pi^4*sig2_str^2;
eqC2 = om1_str^2 - 2*(pi/L_str)^2*sig0_str*sig2_str;
eqC3 = sig0_str^2 + (pi/T)^2;
M_str = floor(LAM_str*sqrt((sqrt(eqC2^2+4*eqC1*eqC3) - eqC2)/(2*eqC1)));

% MODE SHAPES MATRIX
x_str = linspace(0,L_str,N_str);
iM_str = 1:M_str;
W_str = sin(iM_str'*x_str*pi/L_str);

% modal oscillator coefficients
bet_str = iM_str*pi/L_str;
k = 0.5*L_str*(EI_str*bet_str.^4 + Tau_str*bet_str.^2);
alp_str = sig0_str + sig2_str*bet_str.^2;   % Modal decay rates
w_str = sqrt(k/m_str - alp_str.^2);         % Modal frequencies

R_str = exp(-alp_str*T);
Om_str = cos(w_str*T);

a_hat_str = (1 - 2*R_str.*Om_str + R_str.^2)./(1 +2*R_str.*Om_str + R_str.^2);
b_hat_str = 2*(1 - R_str.^2)./(1 + 2*R_str.*Om_str + R_str.^2);
xi_str = (T^2)/m_str;

A_str = 2*(1-a_hat_str) ./ (1+a_hat_str+b_hat_str);
B_str = -(1+a_hat_str-b_hat_str) ./ (1+a_hat_str+b_hat_str);
C_str = xi_str ./ (1+a_hat_str+b_hat_str);

coupPos_str = floor(1 + reCoup_str*(N_str-1));   % Coupling position along string
hitPos_str = floor(1 + reEx_str*(N_str-1));      % Excitation position along string
p1Pos_str = floor(1 + rePick1_str*(N_str-1));    % Pickup position 1
p2Pos_str = floor(1 + rePick2_str*(N_str-1));    % Pickup position 2

% initiation of variables & constants %
Wc_str = W_str(:,coupPos_str);                  % modal coupling weight vector (column)
Wex_str = W_str(:,hitPos_str);                  % modal coupling weight vector (column)
Wp1_str = W_str(:,p1Pos_str);                   % modal pickup weight vector (column) 
Wp2_str = W_str(:,p2Pos_str);                   % modal pickup weight vector (column)

Nmta_str = length(Wc_str);      % total number of modes;

A1_str = -A_str';               % modal oscillator coefficients
A2_str = -B_str';               % modal oscillator coefficients
B1_str = (C_str'.*Wc_str);      % modal oscillator coefficients
B2_str = (C_str'.*Wex_str);     % modal oscillator coefficients

% variables initialisation
Y0_str = zeros(Nmta_str,1);         % modal state (n+1)
Y1_str = zeros(Nmta_str,1);         % modal state (n)
Y2_str = zeros(Nmta_str,1);         % modal state (n-1)
V1_str = 0;                         % velocity at output 1
V2_str = 0;                         % velocity at output 2
Fc_str = 0;                         % coupling force
VC0_str = 0;                        % Displacement (n+1)
VC1_str = 0;                        % Displacement (n)
VC2_str = 0;                        % Displacement (n-1)

% MEMBRANE MODAL PRE-CALCUALTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mth zero of nth order bessel function
N_mem = M_mem;
M_mem = N_mem;
n_mem = (0:N_mem)';
m_mem = 1:M_mem;
bessZeros = besselzero(n_mem,M_mem);

% Loss factor alpha, center freq omega
alp_2_mem = (1/sig_mem)*(dp1_mem+dp3_mem*(bessZeros/Rad_mem).^2);
alp_mem = alp_2_mem / 2;        % Modal decay rates
om_squar_mem = ((bessZeros/Rad_mem).^2) .* ((D_mem/sig_mem) * ((bessZeros/Rad_mem).^2) + c_mem^2);
om_mem = sqrt(om_squar_mem);    % Modal frequencies

% modal oscillator coefficients
Om_mem = cos(om_mem*T);
R_mem = exp(-alp_mem*T);

a_hat_mem = (1-2*R_mem.*Om_mem+R_mem.^2)./(1+2*R_mem.*Om_mem+R_mem.^2);
b_hat_mem = 2*(1-R_mem.^2)./(1+2*R_mem.*Om_mem+R_mem.^2);
xi_mem = (T^2)/sig_mem;

A_mem = 2*(1-a_hat_mem) ./ (1+a_hat_mem+b_hat_mem);
B_mem = (b_hat_mem-a_hat_mem-1) ./ (1+a_hat_mem+b_hat_mem);
C_mem = xi_mem ./ (1+a_hat_mem+b_hat_mem);

% Plotting mode shapes
phi_mem = linspace(0,2*pi,100)';
phi_len_mem = length(phi_mem);
r_mem = linspace(0,Rad_mem,120);
r_len_mem = length(r_mem);
[X_mem,Y_mem] = pol2cart(phi_mem,r_mem);

% Create Matrix of mode shapes
K_mem = zeros(N_mem+1,M_mem,phi_len_mem,r_len_mem);
for i=0:N_mem
    for j=1:M_mem
        K_mem(i+1,j,:,:) = cos(i*phi_mem) * besselj(i,(r_mem/Rad_mem)*bessZeros(i+1,j));
    end
end

% Set exitation and pickup positions
r_coup_mem = reCoup_mem*Rad_mem;   % Coupling pos
phi_coup_mem = 0;
r_hit_mem = reEx_mem*Rad_mem;      % Excitation pos
phi_hit_mem = 0;
r_p1_mem = rePick1_mem*Rad_mem;    % Pickup pos 1
phi_p1_mem = phiPick1;
r_p2_mem = rePick2_mem*Rad_mem;    % Pickup pos 2
phi_p2_mem = phiPick2;

% Modal shapes/weight at excitation and pickup positions
K_coup_mem = zeros(N_mem+1,M_mem);
K_ex_mem = zeros(N_mem+1,M_mem);
K_p1_mem = zeros(N_mem+1,M_mem);
K_p2_mem = zeros(N_mem+1,M_mem);
for i=0:N_mem
    K_coup_mem(i+1,:) = cos(i*phi_coup_mem) * besselj(i,(r_coup_mem/Rad_mem)*bessZeros(i+1,:));
    K_ex_mem(i+1,:) = cos(i*phi_hit_mem) * besselj(i,(r_hit_mem/Rad_mem)*bessZeros(i+1,:));
    K_p1_mem(i+1,:) = cos(i*phi_p1_mem) * besselj(i,(r_p1_mem/Rad_mem)*bessZeros(i+1,:));
    K_p2_mem(i+1,:) = cos(i*phi_p2_mem) * besselj(i,(r_p2_mem/Rad_mem)*bessZeros(i+1,:));
end

% Reshape matrices for loop
A_mem = reshape(A_mem,1,[]);
B_mem = reshape(B_mem,1,[]);
C_mem = reshape(C_mem,1,[]);
K_coup_mem = (reshape(K_coup_mem,1,[]));
K_ex_mem = reshape(K_ex_mem,1,[]);
K_p1_mem = reshape(K_p1_mem,1,[]);
K_p2_mem = reshape(K_p2_mem,1,[]);

% initiation of variables & constants
Wc_mem = K_coup_mem';           % modal coupling weights
W_mem = K_ex_mem';              % modal excitation weights
Wp1_mem = K_p1_mem';            % modal pickup weights
Wp2_mem = K_p2_mem';            % modal pickup weights

Nmta_mem = length(K_coup_mem);  % total number of modes;

A1_mem = -A_mem';               % modal oscillator coefficients
A2_mem = -B_mem';               % modal oscillator coefficients
B1_mem = (C_mem.*K_coup_mem)';  % modal oscillator coefficients
B2_mem = (C_mem.*K_ex_mem)';    % modal oscillator coefficients

% variables initialisation
Y0_mem = zeros(Nmta_mem,1);     % modal state (n+1)
Y1_mem = zeros(Nmta_mem,1);     % modal state (n)
Y2_mem = zeros(Nmta_mem,1);     % modal state (n-1)
V1_mem = 0;                     % velocity at output 1
V2_mem = 0;                     % velocity at output 2
Fc_mem = 0;                     % coupling force
VC0_mem = 0;                    % Displacement (n+1)
VC1_mem = 0;                    % Displacement (n)
VC2_mem = 0;                    % Displacement (n-1)

% COUPLING COEFFICIENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta = (varrho_u*X_u)/(2*T*T);

beta_mem = sum((Wc_mem).*B1_mem);
d1 = -2*mu_u*mu_u - 2*lam_u*lam_u + 2;
d2 = 2*(lam_u*lam_u + 2*mu_u*mu_u);
d3 = -2*mu_u*mu_u;
d4 = -1 - nu1_u;
d5 = -1 + nu1_u;
d6 = (T*T)/varrho_u;

% membrane
dg_mem =  1 - beta_mem*theta*d4;
g1_mem = (beta_mem*theta*d1)/dg_mem;
g2_mem = (beta_mem*theta*d2)/dg_mem;
g3_mem = (beta_mem*theta*d3)/dg_mem;
g4_mem = (beta_mem*theta*d5)/dg_mem;
g5_mem = (beta_mem*theta*d6)/dg_mem;
g6_mem = -1/dg_mem;

% string
beta_str = sum((Wc_str).*B1_str);
dg_str =  1 - beta_str*theta*d4;
g1_str = (beta_str*theta*d1)/dg_str;
g2_str = (beta_str*theta*d2)/dg_str;
g3_str = (beta_str*theta*d3)/dg_str;
g4_str = (beta_str*theta*d5)/dg_str;
g5_str = (beta_str*theta*d6)/dg_str;
g6_str = -1/dg_str;


% MAIN LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n=1:Ns
    % FORCE INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    F1(ke) = inpSlinky(n);
    
    % slinky %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    k = 2;        % near left boundary
    U0(k) = A_u(1)*(U1(k+2) + U1(k-1) + U1(k-1) - U1(k)) + A_u(2)*(U1(k+1) + U1(k-1)) + A_u(3)*U1(k) + A_u(4)*U2(k) + A_u(5)*(U2(k+1) + U2(k-1)) + A_u(6)*F1(k);    
    k = 3:Np_u-2;   % interior     
    U0(k) = A_u(1)*(U1(k+2) + U1(k-2)) + A_u(2)*(U1(k + 1) + U1(k-1)) + A_u(3)*U1(k)  +  A_u(4)*U2(k) + A_u(5)*(U2(k+1) + U2(k-1)) + A_u(6)*F1(k);    
    k = Np_u-1;     % near right boundary
    U0(k) = A_u(1)*(U1(k+1) + U1(k+1) - U1(k) + U1(k-2)) + A_u(2)*(U1(k+1) + U1(k-1)) + A_u(3)*U1(k) + A_u(4)*U2(k) + A_u(5)*(U2(k+1) + U2(k-1)) + A_u(6)*F1(k);    
     
    % STRING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % connection point
    k = 1;
    H1 = sum( Wc_str.*(A1_str.*Y1_str + A2_str.*Y2_str - B2_str*inpString(n)));
    VC0_str = g1_str*VC1_str + g2_str*U1(k+1) + g3_str*U1(k+2) + g4_str*VC2_str  + g5_str*F1(k) + g6_str*H1;
    U0(k) = VC0_str;
    % modal oscillators
    Fc_str = theta*(d1*VC1_str + d2*U1(k+1) + d3*U1(k+2) + d4*VC0_str + d5*VC2_str + d6*F1(k));
    Y0_str = B1_str*Fc_str - A1_str.*Y1_str - A2_str.*Y2_str;
    % pick-up points
    V1_str = (Wp1_str')*(Y0_str - Y2_str)/(2*T);                
    V2_str = (Wp2_str')*(Y0_str - Y2_str)/(2*T);   

    % MEMBRANE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % connection point
    k = Np_u;
    H1 = sum( Wc_mem.*(A1_mem.*Y1_mem + A2_mem.*Y2_mem - B2_mem*inpMem(n)));
    VC0_mem = g1_mem*VC1_mem + g2_mem*U1(k-1) + g3_mem*U1(k-2) + g4_mem*VC2_mem  + g5_mem*F1(k) + g6_mem*H1;
    U0(k) = VC0_mem;
    % modal oscillators
    Fc_mem = theta*(d1*VC1_mem + d2*U1(k-1) + d3*U1(k-2) + d4*VC0_mem + d5*VC2_mem + d6*F1(k));
    Y0_mem = B1_mem*Fc_mem - A1_mem.*Y1_mem - A2_mem.*Y2_mem;
    % pick-up points
    V1_mem = (Wp1_mem')*(Y0_mem - Y2_mem)/(2*T);                
    V2_mem = (Wp2_mem')*(Y0_mem - Y2_mem)/(2*T);                
        
    %%% record output signals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    outp1(n) = V1_mem;
    outp2(n) = V2_mem;
    outp3(n) = Fc_mem;
    outp4(n) = VC0_mem;
           
    %%% memorise variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    U2 = U1;
    U1 = U0;
    Y2_str = Y1_str;
    Y1_str = Y0_str;
    VC2_str = VC1_str;
    VC1_str = VC0_str;
    Y2_mem = Y1_mem;
    Y1_mem = Y0_mem;
    VC2_mem = VC1_mem;
    VC1_mem = VC0_mem;

end

% SOUND
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Release Ramp
ramp = ones(1,Ns);  % Initialise ramp with ones
onsetDur = 0.005;   % Assign ramp onset length
NsOnset = floor(onsetDur*Fs);   % Length of onset/offset in samples
ramp((end-NsOnset+1):end) = (NsOnset-1:-1:0) / (NsOnset-1); % Release of ramp

% Scale output signal
m1 = max(abs(outp1));
m2 = max(abs(outp2));
maxAmp = max([m1 m2]);
outp1 = (0.99/maxAmp)*outp1.*ramp;
outp2 = (0.99/maxAmp)*outp2.*ramp;
outp = [outp1' outp2'];

% Output sound
soundsc(outp,Fs);
audiowrite('realisticYaybahar.wav',outp,Fs);
                                                                                                                                                                                                                                                                                                                                                                                  