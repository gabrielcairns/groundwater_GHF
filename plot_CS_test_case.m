% Produces Figure 7 of the manuscript "Modification of Antarctic geothermal heat flux by groundwater flow"
% Requires cmocean package (https://mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps)

% Horizontal meshgrid
xg = 0.5;
XH = linspace(0,1,301)';
x = xg*XH;
Lx = 301;

% Vertical meshgrid
Lz = 20;
Z = linspace(0,1,Lz);
[ZZ, xx] = meshgrid(Z,x);

% Loading parameter values
[scales, params, ~, ~] = default_parameter_values();
Hd = scales.Hd; xdkm = scales.xdkm; Gd = scales.Gd;
delta = params.delta;  
mdmmyr = scales.mdmmyr;

% Bed, sedimentary basin, heat flux, thinning rate
bedf = @(x) -(1000/Hd)*sin(pi*x/.5);
bedxf = @(x) -(pi/.4)*(1000/Hd)*cos(pi*x/.5);
Hsbf= @(x) (2000/Hd)*(.5 + .5*tanh( (x-.05)/.05).*tanh( (.45-x)/.05)) ;
Gf = @(x) (0.08/scales.Gd)*ones(size(x));
psit = zeros(size(x));

G = Gf(x);
Hsb = Hsbf(x);
b = bedf(x);
bx = bedxf(x);

% Ice thickness
Hi = (1000/Hd)-b;
% Shreve potential
psi_Shr = (1-delta)*Hi + b;

% Creation of unmodified case
if ~isfile('test_case/test_CS_Ksb_0.mat')
    F_sb = G;
    Tsb = -G.*Hsb.*(ZZ-1);
    save('test_case/test_CS_Ksb_0.mat',"F_sb","Tsb");
end

% Values of Ksb to calculate
Ksbs = [0 1e-6 3e-6];
Ksbp = 1e-6; % Particular value to plot
% Text for use in legend
Ksb_text = legend_text_K_Pe(Ksbp, params.St*scales.Ksb_to_kappa*Ksbp);

Pes = params.St*scales.Ksb_to_kappa*Ksbs;

for K_sb = Ksbs(2:end)

    % If the solution isn't already saved, we generate and save it
    filename = sprintf('test_case/test_CS_Ksb_%g.mat',K_sb);
    if ~isfile(filename)
        [Tsb,F_sb,converged] = test_case(K_sb,0,'CS');
        save(filename,"Tsb","F_sb");
    end

end

% Initialisation of figure
figure(1); clf; set(gcf(),'Position',[100  100  1075  775])
tiledlayout(3,1)
nexttile([2 1]);
aT = make_nice_axes(gca(), '', '$z$ / m');  nexttile;
aF = make_nice_axes(gca(), '$y$ / km', '$-k_{sb}\partial T/ \partial z$ / mW m$^{-2}$');  %nexttile;

xlim(aT,[0 xdkm*xg]);   colormap(aT,cmocean('thermal'));    aT.XTickLabels = {};
xlim(aF,[0 xdkm*xg]); 

L = length(Ksbs);
clrs = cmocean('haline',L+2);
clrs = clrs([1 3:L+1],:);

% Plot of bed
plot(aT,xdkm*xg*XH, Hd* b, 'k', 'LineWidth',2);

% Plot of heat flux in each case
for l = 1:L
    load(sprintf('test_case/test_CS_Ksb_%g.mat',Ksbs(l)),"F_sb");    
    plot(aF,xdkm*xg*XH, 1000*Gd*F_sb, 'LineWidth',2, 'Color',clrs(l,:));
end

% Hydraulic potential
psi = psi_Shr- .1*(1-delta)*Hi;
plot(aT,xdkm*xg*XH, Hd*psi, 'k-' , 'LineWidth',2);

% Ice surface
plot(aT,xdkm*xg*XH, Hd*(Hi+ b), 'k--', 'LineWidth',2);
ylim(aF,[30 130]);

% Particular solution for temperature 
load(sprintf('test_case/test_CS_Ksb_%g.mat',3e-6),"Tsb");
contourf(aT,xdkm*xx,Hd*(Hsb.*(ZZ-1)+b),scales.DeltaT*Tsb,0:5:100,'EdgeColor','none'); 

ylim(aT, [-3100 1100])

% Colorbar
cb = make_nice_colorbar(aT,'east','$T$ / $^o$C');
cb.Position =[0.8800 0.4200 0.0200 0.2000];

% Velocity streamlines
u1 = -gradient(psi,xx(:,Lz));
u = u1;
u(1) =  3*u1(2)-3*u1(3)+u1(4);
u(Lx) =  3*u1(Lx-1)-3*u1(Lx-2)+u1(Lx-3);
contour(aT, xdkm*xx,Hd*(Hsb.*(ZZ-1)+b),u.*Hsb.*ZZ,20,'w')

% Adding legends
legtext = legend_text_K_Pe(Ksbs,Pes);
legend(aT, aT.Children([3 4]), {'$H_i+b$','$\psi$'}, 'FontSize',14, 'Interpreter','latex', 'NumColumns',1 ,'Position',[0.13 0.8 0.1 0.07]);
legend(aF,legtext, 'FontSize',14, 'Interpreter','latex','Location', 'NorthOutside','NumColumns',3);

text(aT,5,-2900,legend_text_K_Pe(3e-6, params.St*scales.Ksb_to_kappa*3e-6),'FontSize',14, 'Interpreter','latex');
% Labelling subplots
text(aT,240,500,'(a)','BackgroundColor','w','FontSize',14,'Interpreter','latex');
text(aF,240,120,'(b)','BackgroundColor','w','FontSize',14,'Interpreter','latex');