% Calculation of modified temperature and heat flux for Figures 3, 4, 7, 9
% and 10

function [Tsb,F_sb, converged] = test_case(K_sb,Ss,case_to_test)
options = optimoptions('fsolve','Display','off');
% Loading default parameter values
[scales, params, bedf, bedxf] = default_parameter_values();
Hd = scales.Hd;
delta = params.delta; 

% Updating the parameter values
kappa = scales.Ksb_to_kappa*K_sb;
Sigma = scales.Ss_to_Sigma*Ss;
Pe_sb = kappa*params.St;
zeta = Sigma/kappa;
params.Pe_sb = Pe_sb;
params.kappa = kappa;
params.zeta= zeta;
params.Sigma = Sigma;

% Loading variable grid and ice sheet model
load("variable_grid_300.mat","XH");
if isfile('default_weertman_ice.mat')
    load('default_weertman_ice.mat','Hi','xg');
else
    % If ice thickness data isn't there, generate it
    [Hi, xg] = make_default_weertman_ice();
    save('default_weertman_ice.mat','Hi','xg');
end

% Case-dependent: sedimentary basin thickness and GHF
if isequal(case_to_test,'default')
    Hsbf= @(x) (2000/Hd)*exp( - ((x-.7)/.3).^2) ;
    Gf = @(x) (0.08/scales.Gd)*ones(size(x));
elseif isequal(case_to_test,'alt')
    Hsbf= @(x) (2000/Hd)*(.5 - .5*tanh( (x -1.1 )/.1) -.8*exp(-((x-.4)/.1).^2 ) );
    Gf = @(x) (0.06/scales.Gd)*(1 + .5*exp(- ((x-.9)/.1).^2 ));
elseif isequal(case_to_test,'CS')
    Hsbf= @(x) (2000/Hd)*(.5 + .5*tanh( (x-.05)/.05).*tanh( (.45-x)/.05)) ;
    Gf = @(x) (0.08/scales.Gd)*ones(size(x));
else
    error("Test case must be one of the following: 'default' 'alt 'CS'");
end

% Cross-section case (Figure 7) is a little different to the others
if isequal(case_to_test,'CS')
    params.Sigma = 0;
    bedf = @(x) -(1000/Hd)*sin(pi*x/.5);
    bedxf = @(x) -(pi/.4)*(1000/Hd)*cos(pi*x/.5);
    % Uniform grid
    XH = linspace(0,1,301)';
    xg = 0.5;
    Hi = (1000/Hd)-bedf(xg*XH);
end

x = xg*XH;

psit = -(1/0.3)*ones(size(x)); % Just uniform, 1 m/yr

b = bedf(x); bx = bedxf(x);
Hsb = Hsbf(x); 
G = Gf(x);

% Vertical meshgrid
Lz = 20;
Z = linspace(0,1,Lz);
[ZZ, xx] = meshgrid(Z,x);

% Hydraulic potential
psi = (0.9*(1-delta)*Hi+b);

% Inital guess is unmodified case
Tsbi = -G.*(ZZ-b);

% Nonlinear solve for Tsb
[Tsb, ~, ex]  = fsolve(@(t) OBJ_FUN_heat_sb_SS(t,xx,ZZ,psi,zeta.*psit,G, bx, gradient(Hsb,x), Hsb,Pe_sb,0),Tsbi,options);

% Calculate heat flux at top
F_sb = -(Tsb(:,Lz)-Tsb(:,Lz-1))./(Hsb.*(ZZ(:,Lz)-ZZ(:,Lz-1)));

converged = (ex>0);

end