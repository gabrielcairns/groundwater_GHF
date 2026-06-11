% Function to generate default parameter values, scalings, and bed
% geometry
function [scales, params, bedf, bedxf] = default_parameter_values()

% Porosity
phi = 0.2;

% Water latent heat, gravity, water viscosity
L = 3.3e5;
g = 9.81;
eta = 1e-3;

% Densities
rho_w = 1000;
rho_r = 2000;
rho_i = 900;

% Specific heat capacities
C_w = 4200;
C_r = 700;
rhoC_sb = rho_r*C_r*(1-phi) + rho_w*C_w*phi;

% Thermal conductivity
k_w = 0.56;
k_r = 2.8;
k_sb = k_r*(1-phi) + k_w*phi;

K_sb = 1e-13*rho_w*g/eta;

% Ice sheet model
betaW = 0.1*7.624e6;
A = 4.227e-25;
n = 3;
params.n = n;
mW = 1/3;
params.mW = mW;

% Specific storage, Skempton coeff.
Ss = 1e-5;
scales.xi = 0.2;

% Seconds per year, useful for conversions
spery = 31557600;

% Horizontal lengthscale, accumulation
xd = 5e5;   
scales.xdkm = xd/1e3;
ad = 0.3/spery; 

% Vertical lengthscale, timescale
uid = (rho_i*g*ad^2*xd/betaW)^(1/(mW+2));    
scales.uidmyr = uid*spery;
Hd = (ad*xd)/uid;
scales.Hd = Hd;
td = xd/uid;
scales.tdyr = td/spery;
scales.tkyr = scales.tdyr/1e3;

% Used in the ice sheet model
params.eps = (uid/(xd*A))^(1/n) / (2*rho_i*g*Hd);

% Density difference between ice and water
params.delta = 1 - rho_i/rho_w ;

% GHF and temperature scale
scales.Gd = 0.06;
DeltaT = scales.Gd*Hd/k_sb;
scales.DeltaT = DeltaT;

% Melt rate
md = k_sb*DeltaT/(rho_w*L*Hd);
scales.md = md;
scales.mdmmyr = md*1000*spery; 

% 'Stefan number': = Pe_sb / kappa 
params.St = C_w*DeltaT/L;

% Kappa and Sigma, and some conversion factors
% We have an additional zeta := kappa / sigma
usbd = K_sb*Hd/xd;
params.kappa = (Hd/xd)*usbd/md;
scales.Ksb_to_kappa = (Hd/xd).^2/md;
scales.Ss_to_zeta = (1-scales.xi)*xd^2/(K_sb*td);
scales.Ss_to_Sigma = (1-scales.xi)*Hd^2/(md*td);
params.Sigma = Ss*scales.Ss_to_Sigma;
params.zeta = Ss*scales.Ss_to_zeta;

% Other parameters
params.lambda_sb = rhoC_sb*uid/(rho_w*C_w*usbd);
params.Pe_sb = Hd^2*usbd*rho_w*C_w/(xd*k_sb);
params.Xi_sb = g*Hd/(C_w*DeltaT);

% Bed geometry
bedf = @(x) -100/Hd -0.2*x;
bedxf = @(x) -0.2*ones(size(x));


end