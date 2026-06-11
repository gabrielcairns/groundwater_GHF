% Nonlinear solver of the heat equation in the sedimentary basin

function fun = OBJ_FUN_heat_sb_SS(T,xx,ZZ,psi,zeta_psit,G,bx,Hsbx,Hsb,Pe,Xi)
% ZZ = (z-(b-Hsb))/Hsb is a transformed coordinate to make the grid
% rectangular

% Meshgrids have x increasing down the columns and z increasing along rows
Lx = size(xx,1);    Lz = size(xx,2);
% psi known at whole x midpoints

% Finding u_sb, ensuring smoothness into final cells
u1 = -gradient(psi,xx(:,Lz));
u = u1;
u(1) =  3*u1(2)-3*u1(3)+u1(4);
u(Lx) =  3*u1(Lx-1)-3*u1(Lx-2)+u1(Lx-3);

% Half-points in x
xh = [xx(1,Lz)-0.5*(xx(2,Lz)-xx(1,Lz)); 0.5*(xx(2:Lx,Lz)+xx(1:Lx-1,Lz)); xx(Lx,Lz)+0.5*(xx(Lx,Lz)-xx(Lx-1,Lz))];

% w at whole x and half Z
psixx = -(u(2:Lx)-u(1:Lx-1))./(xh(2:Lx)-xh(1:Lx-1));
psixx = [2*psixx(1)-psixx(2); psixx];
w = (psixx - zeta_psit).*Hsb.*(ZZ(1,1:Lz-1)-ZZ(1,1)) + u(1:Lx).*(bx-Hsbx);

% Half-points in Z
Zh = [ZZ(:,1)-0.5*(ZZ(:,2)-ZZ(:,1)) 0.5*(ZZ(:,2:Lz)+ZZ(:,1:Lz-1)) ZZ(:,Lz)+0.5*(ZZ(:,Lz)-ZZ(:,Lz-1))];

% Horizontal advection: upwind method
advx = [zeros(1,Lz); (T(1:Lx-1,:).*(u(1:Lx-1,:)>=0) + T(2:Lx,:).*(u(1:Lx-1,:)<0)).*u(1:Lx-1,:);  T(Lx,:).*u(Lx,:)];
advx(1,:) = 3*advx(2,:)-3*advx(3,:)+advx(4,:);

% Coordinate change variable
dZdx = -(bx+Hsbx.*(ZZ-1))./Hsb;

% Extra advection in Z direction due to coordinate change
advxZ = 0.5*(advx(1:Lx,:)+advx(2:Lx+1,:));

% Vertical advection: upwind method
advz = [(bx-Hsbx).*advxZ (T(:,1:Lz-1).*(w>=0) + T(:,2:Lz).*(w<0)).*w];

% Diffusive heat flux (including prescribed GHF)
diffz = [-G (T(:,2:Lz)-T(:,1:Lz-1))./( Hsb.*(ZZ(:,2:Lz)-ZZ(:,1:Lz-1)) )];

% Dissipation
dissip = Xi*u(1:Lx).^2;

% Final objective function
fun(1:Lx,1:Lz-1) = Pe*(advx(2:Lx+1,1:Lz-1)-advx(1:Lx,1:Lz-1))./(xh(2:Lx+1)-xh(1:Lx)) + ...
     Pe*dZdx(:,1:Lz-1).*(advxZ(:,2:Lz) -advxZ(:,1:Lz-1))./(Zh(:,2:Lz)-Zh(:,1:Lz-1)) + ...
     (Pe*(advz(:,2:Lz) -advz(:,1:Lz-1)) -(diffz(:,2:Lz) -diffz(:,1:Lz-1) )  )./( Hsb.*(Zh(:,2:Lz)-Zh(:,1:Lz-1))) ...
     - Pe*dissip + Pe*zeta_psit.*T(:,1:Lz-1);

% Dirichlet boundary condition at top
fun(1:Lx,Lz) = T(:,Lz);
end