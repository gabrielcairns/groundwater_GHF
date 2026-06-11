% Makes the default weertman ice sheet
function [Hi, xg] = make_default_weertman_ice()
options = optimoptions('fsolve','Display','off');

% Loading parameters
[~, params, bedf, bedxf] = default_parameter_values();
delta = params.delta;   n = params.n;    mW = params.mW;    eps = params.eps;  
params.Sigma = 0;

% Variable meshgrid
load("variable_grid_300.mat","XA","XH","XU");
J = 300;

% Sliding law
fricf = @(u) abs(u + 1e-6).^(mW-1).*u;

% Accumulation
accf = @(x) ones(size(x));
intaccf = @(x) x;

% Initial rough guess for xg (see Schoof, 2007)
fxg = @(x) ( (delta/(8*eps))^n / (1-delta)^(mW+n+3) )^(1/(mW+1)) * (-bedf(x)).^((n+mW+3)/(mW+1)) - intaccf(x);
xgR = fsolve(@(x) fxg(x), 2,options); 

% Initial guess for H and U
HR = ones(size(XH));
HR(end) = -bedf(xgR)/(1-delta);

for k = length(XH):-1:2
    HxT = -(intaccf(xgR*XH(k))).^mW./HR(k).^(mW+1) - bedxf(xgR*XH(k)) ;
    HR(k-1) = HR(k)- xgR*(XH(k)-XH(k-1))*HxT;
end
UR = [0; intaccf(xgR*XU(2:end))./HR];
HUx = [HR; UR; xgR];

% Nonlinear solve for H, U and xg
obj = @(hux) OBJ_FUNC_HUx(hux, hux([1:J+1 2*J+4]), XA, bedf, accf, fricf, 100, params);
HUx = fsolve(obj, HUx);

xg = HUx(2*J+4);
Hi = HUx(1:J+1);

end