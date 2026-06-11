% Nonlinear solver for ice sheet model with a general sliding law in terms
% of u
function fun = OBJ_FUNC_HUx(HUx, HxI, XA, bedf, accf, fricfu, dt, params)

delta = params.delta;   n = params.n;    
eps = params.eps; 

J = (length(XA)-3)/2;

% Mesh grids
XH = XA(2:2:2*J+2);
XU = XA(1:2:2*J+3);

% Current values to solve for
H = HUx(1:J+1);
U = HUx(J+2:2*J+3);
xg =  HUx(2*J+4);

% Previous value
HI = HxI(1:J+1);
xgI = HxI(J+2);

% Time derivatives
xgt = (xg - xgI)/dt;

a = accf(xg*XH);
bed = bedf(xg*XH);

% MOMENTUM EQUATION:
Ux = (xg^(-1))*(U(2:J+2)-U(1:J+1))./(XU(2:J+2)-XU(1:J+1));

exstr = H.* (Ux.^2 + 1e-4).^((1-n)/(2*n)).*Ux ;

% Extensional stress
tauex = 4*eps* (xg^(-1))* (exstr(2:J+1)-exstr(1:J))./(XH(2:J+1)-XH(1:J)) ;
% Driving stress
taudr = .5*(H(2:J+1)+H(1:J) ) .*  (xg^(-1)).* (H(2:J+1)-H(1:J) + bed(2:J+1)-bed(1:J))./(XH(2:J+1)-XH(1:J));

% Basal stress (general sliding law)
tausl = fricfu(U(2:J+1));

% Flux condition
fun(1) = U(1)+U(2); 

% Momentum balance
fun(2:J+1) = tauex - tausl - taudr;
% Extensional stress condition
fun(J+2) = U(J+2)-U(J+1) - xg*(XU(J+2)-XU(J+1))*(delta*H(J+1)/(8*eps)).^n;

% MASS EQUATION:
Hp = [H(1)+bed(2)-bed(1); H]; % minus-oneth grid cell
XHp = [-XH(2); XH];
% Advection (downwind method works best, by experience)
HU = Hp.*U(1:J+2); 
divflux = (xg^(-1))*(HU(2:J+2)-HU(1:J+1))./(XU(2:J+2)-XU(1:J+1)); 
% Accounting for moving grid
GLdrift = (xgt/xg).*XHp(1:J+1).*( (Hp(2:J+2)-Hp(1:J+1))./(XHp(2:J+2)-XHp(1:J+1)));

% Mass balance
dHdt = GLdrift + a - divflux;

% Mass equation
fun(J+3:2*J+3) = H-HI-dt*dHdt;
% Flotation condition
fun(2*J+4) = H(J+1)+bed(J+1)/(1-delta);
end