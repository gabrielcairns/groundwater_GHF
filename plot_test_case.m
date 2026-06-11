% Produces Figure 3, 4, 9 and 10 of the manuscript "Modification of Antarctic geothermal heat flux by groundwater flow"
% Requires cmocean package (https://mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps)

% Prompt user to choose between 'default' and 'alt' via a dropdown dialog
% 'default' produces Figures 3 and 4, 'alt' produces 9 and 10
[idx,tf] = listdlg('PromptString','Select test case to plot:',...
                   'SelectionMode','single',...
                   'ListString',{'default','alt'},...
                   'ListSize',[160 60]);
if ~tf || idx==1
    case_to_test = 'default';
else
    case_to_test = 'alt';
end

if isequal(case_to_test,'default')
    name_part = 'test_def';
elseif isequal(case_to_test,'alt')
    name_part = 'test_alt';
end

% Loading grid and ice model (creating this if it doesn't exist already)
load("variable_grid_300.mat","XH");
if isfile('default_weertman_ice.mat')
    load('default_weertman_ice.mat','Hi','xg');
else
    [Hi, xg] = make_default_weertman_ice();
    save('default_weertman_ice.mat','Hi','xg');
end

% Meshgrid
x = xg*XH;
Lx = 301;

% Vertical meshgrid
Lz = 20;
Z = linspace(0,1,Lz);
[ZZ, xx] = meshgrid(Z,x);

% Loading parameter values and bed geometry
[scales, params, bedf, bedxf] = default_parameter_values();
Hd = scales.Hd; xdkm = scales.xdkm; Gd = scales.Gd;
delta = params.delta;
mdmmyr = scales.mdmmyr;

% Sedimentary basin thickness and GHF
if isequal(case_to_test,'default')
    Hsbf= @(x) (2000/Hd)*exp( - ((x-.7)/.3).^2) ;
    Gf = @(x) (0.08/scales.Gd)*ones(size(x));
elseif isequal(case_to_test,'alt')
    Hsbf= @(x) (2000/Hd)*(.5 - .5*tanh( (x -1.1 )/.1) -.8*exp(-((x-.4)/.1).^2 ) );
    Gf = @(x) (0.06/scales.Gd)*(1 + .5*exp(- ((x-.9)/.1).^2 ));
end
% Thinning rate
psit = -(1/0.3)*ones(size(x));

G = Gf(x);    Hsb = Hsbf(x);
b = bedf(x);
bx = bedxf(x);

psi_Shr = (1-delta)*Hi + b;

% Generation of unmodified case
if ~isfile(strcat('test_case/',name_part,'_Ksb_0_Ss_0.mat'))
    F_sb = G;
    Tsb = -G.*Hsb.*(ZZ-1);
    qE = zeros(size(x));
    save(strcat('test_case/',name_part,'_Ksb_0_Ss_0.mat'),"F_sb","Tsb");
end

% List of Ksb to find
Ksbs = [0 1e-6 3e-6 5e-6];
Ksbp = 1e-6; % Particular one to vary Ss for
% Text for legend
Ksb_text = legend_text_K_Pe(Ksbp, params.St*scales.Ksb_to_kappa*Ksbp);

Pes = params.St*scales.Ksb_to_kappa*Ksbs;
% For each Ksb, generate the solution if it doesn't exist already
for K_sb = Ksbs(2:end)

    filename = sprintf('test_case/%s_Ksb_%g_Ss_0.mat',name_part,K_sb);
    if ~isfile(filename)
        [Tsb,F_sb,converged] = test_case(K_sb,0,case_to_test);
        save(filename,"Tsb","F_sb");
    end
end

% Likewise for different values of Ss
Sss = [0 3e-7 1e-6 3e-6];
for Ss = Sss
    filename = sprintf('test_case/%s_Ksb_%g_Ss_%g.mat',name_part,Ksbp,Ss);
    if ~isfile(filename)
        [Tsb,F_sb,converged] = test_case(Ksbp,Ss,case_to_test);
        save(filename,"Tsb","F_sb");
    end 
end

% Figure 3 / 9 showing setup, temperature, and perturbed heat flux
% Axes initialisation
figure(1); clf; set(gcf(),'Position',[370 130 1130 1000])
tiledlayout(4,1,'TileSpacing','Compact','Padding','compact');
nexttile(1,[2 1]);
aT = make_nice_axes(gca(), '', '$z$ / m');  nexttile;
aF1 = make_nice_axes(gca(), '', '$-k_{sb}\partial T/ \partial z$ / mW m$^{-2}$');  nexttile;
aF2 = make_nice_axes(gca(), '$x$ / km', '$-k_{sb}\partial T/ \partial z$ / mW m$^{-2}$'); 
xlim(aT,[0 xdkm*xg]);   colormap(aT,cmocean('thermal'));    aT.XTickLabels = {};
xlim(aF1,[0 xdkm*xg]);   aF1.XTickLabels = {};
xlim(aF2,[0 xdkm*xg]);

L = length(Ksbs);

clrs = cmocean('haline',L+2);
clrs = clrs([1 3:L+1],:);

clrs2 = cmocean('ice',length(Sss)+2);

% Plotting various heat fluxes
for l = 1:L
    load(sprintf('test_case/%s_Ksb_%g_Ss_0.mat',name_part,Ksbs(l)),"F_sb");    
    plot(aF1,xdkm*xg*XH, 1000*Gd*F_sb, 'LineWidth',2, 'Color',clrs(l,:));
end

for l = 1:length(Sss)
    load(sprintf('test_case/%s_Ksb_1e-06_Ss_%g.mat',name_part,Sss(l)),"F_sb");
    plot(aF2,xdkm*xg*XH, 1000*Gd*F_sb, 'LineWidth',2, 'Color',clrs2(l,:));
end

% Hydraulic potential
psi = psi_Shr-0.1*(1-delta)*Hi;
plot(aT,xdkm*xg*XH, Hd*psi, 'k-' , 'LineWidth',2);

% Ice surface
plot(aT,xdkm*xg*XH, Hd*(Hi+ b), 'k--', 'LineWidth',2);

if isequal(case_to_test,'default')
    ylim(aF1,[60 120]);
    ylim(aF2,[60 120]);
elseif isequal(case_to_test,'alt')
    ylim(aF1,[30 130]);
    ylim(aF2,[30 130]);
end

% Temperature contour plot
load(sprintf('test_case/%s_Ksb_%g_Ss_0.mat',name_part,5e-6),"Tsb");
contourf(aT,xdkm*xx,Hd*(Hsb.*(ZZ-1)+b),scales.DeltaT*Tsb,0:5:100,'EdgeColor','none'); 

% Colorbar
cb = make_nice_colorbar(aT,'east','$T$ / $^o$C');
cb.Position = [.91 .55 .02 .19];
clim(aT,[0 min(scales.DeltaT*max(Tsb,[],'all'),100)])

ylim(aT, [-2300 1200])

% Streamlines of groundwater
u1 = -gradient(psi,xx(:,Lz));
u = u1;
u(1) =  3*u1(2)-3*u1(3)+u1(4);
u(Lx) =  3*u1(Lx-1)-3*u1(Lx-2)+u1(Lx-3);
contour(aT, xdkm*xx,Hd*(Hsb.*(ZZ-1)+b),u.*Hsb.*ZZ,20,'w')

% Legend
legtext1 = legend_text_K_Pe(Ksbs,Pes);
legtext2 = legend_text_Ss_Sigma(Sss,Sss*scales.Ss_to_Sigma);

legend(aT, aT.Children([3 4]), {'$H_i+b$','$\psi$'}, 'FontSize',14, 'Interpreter','latex','Location','northeast', 'NumColumns',2)
legend(aF1,legtext1, 'FontSize',14, 'Interpreter','latex','Location', 'NorthWest','NumColumns',2);
legend(aF2,legtext2, 'FontSize',14, 'Interpreter','latex','Location', 'NorthWest','NumColumns',1);
  
% Text describing particular case
text(aF2,510,100,Ksb_text{1},'FontSize',14, 'Interpreter','latex');
text(aT,10,300,legend_text_K_Pe(5e-6, params.St*scales.Ksb_to_kappa*5e-6),'FontSize',14, 'Interpreter','latex');

% Subplot labels
text(aT,680,600,'(a)','BackgroundColor','w','FontSize',14,'Interpreter','latex');
text(aF1,680,110,'(b)','BackgroundColor','w','FontSize',14,'Interpreter','latex');
text(aF2,680,110,'(c)','BackgroundColor','w','FontSize',14,'Interpreter','latex');

% Figure 4 / 10 showing numerical solution compared to asymptotic
% Axes initalisation
figure(2); clf;  set(gcf(),'Position', [100 100 1020 700])
tiledlayout(2,1,'TileSpacing','Compact','Padding','compact'); nexttile;
aG1 = make_nice_axes(gca(),'$x$ / km', '$-k_{sb}\partial T/ \partial z$ / mW m$^{-2}$');
xlim(aG1,xdkm*[0 xg]);
nexttile;
aG2 = make_nice_axes(gca(),'$x$ / km', '$-k_{sb}\partial T/ \partial z$ / mW m$^{-2}$');
xlim(aG2,xdkm*[0 xg]);

clrs = cmocean('haline',L+2);
clrs = clrs([1 3:L+1],:);

% For each Ksb, calculate asymptotic solution and plot the two together
psix = gradient(psi,x);
for l = 1:length(Ksbs)
    Pe = params.St*scales.Ksb_to_kappa*Ksbs(l);

    load(sprintf('test_case/%s_Ksb_%g_Ss_0.mat',name_part,Ksbs(l)),"F_sb");

    corrflux =  0.5*Pe*G.*Hsb.^2.*psix;

    F_sb_asymp = G+gradient(corrflux,x);

    plot(aG1,xdkm*xg*XH,1000*Gd*F_sb,'-','Color',clrs(l,:),'LineWidth',2);
    plot(aG1,xdkm*xg*XH,1000*Gd*F_sb_asymp,'--','Color',clrs(l,:),'LineWidth',2);

end

Pe = params.St*scales.Ksb_to_kappa*Ksbp;

% Likewise for Ss
for l = 1:length(Sss)
    zeta = scales.Ss_to_Sigma*Sss(l)./(scales.Ksb_to_kappa*Ksbp);

    load(sprintf('test_case/%s_Ksb_%g_Ss_%g.mat',name_part,Ksbp,Sss(l)),"F_sb");


    corrflux =  0.5*Pe*G.*Hsb.^2.*psix;

    F_sb_asymp = G+gradient(corrflux,x) - 0.5*Pe*zeta*G.*Hsb.^2.*psit;

    plot(aG2,xdkm*xg*XH,1000*Gd*F_sb,'-','Color',clrs2(l,:),'LineWidth',2);
    plot(aG2,xdkm*xg*XH,1000*Gd*F_sb_asymp,'--','Color',clrs2(l,:),'LineWidth',2);
end
% Legends
legtext1 = legend_text_K_Pe(Ksbs,Pes);
legend(aG1,aG1.Children(end:-2:2),legtext1,'fontsize',14,'interpreter','latex','Location','northwest');
legtext2 = legend_text_Ss_Sigma(Sss,scales.Ss_to_Sigma*Sss);
legend(aG2,aG2.Children(end:-2:2),legtext2,'fontsize',14,'interpreter','latex','Location','northwest');
% Label particular case
text(aG2,510,100,Ksb_text{1},'FontSize',14, 'Interpreter','latex');
% Subplot labels
text(aG1,680,92,'(a)','BackgroundColor','w','FontSize',14,'Interpreter','latex');
text(aG2,680,105,'(b)','BackgroundColor','w','FontSize',14,'Interpreter','latex');