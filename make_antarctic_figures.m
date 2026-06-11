% Generates Figures 1, 5, 6, 8, 11 and 12 of the manuscript "Modification
% of Antarctic geothermal heat flux by groundwater flow" 

% NOTE: This code requires several packages and datasets in order to run,
% which are listed as follows:

% Packages:
% Antarctic Mapping Tools (https://mathworks.com/matlabcentral/fileexchange/47638-antarctic-mapping-tools)
% MEaSUREs (https://mathworks.com/matlabcentral/fileexchange/47329-measures)
% Antarctic boundaries, grounding line and masks from InSAR (https://mathworks.com/matlabcentral/fileexchange/60246-antarctic-boundaries-grounding-line-and-masks-from-insar)
% Hatchfill2 (https://mathworks.com/matlabcentral/fileexchange/53593-hatchfill2)
% cmocean perceptually-uniform colormaps (https://mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps)

% Data:
% "/Bedmap3/" 
% Link: https://data.bas.ac.uk/full-record.php?id=GB/NERC/BAS/PDC/01615
% Pritchard et al. (2025), Bedmap3 updated ice bed, surface and thickness gridded datasets for Antarctica. Scientific data, 12(1), 414.

% "AIS_BaTh_v1.nc"
% Link: https://zenodo.org/records/15556691
% Seiner, O. et al. (2025). A synthesis of the basal thermal state of the Antarctic ice sheet. Journal of Glaciology, 1-30.

% "Ant_Crust.mat"
% Link: https://zenodo.org/records/10242299
% Li, L., & Aitken, A. R. A. (2024). Crustal heterogeneity of Antarctica signals spatially variable radiogenic heat production. Geophysical Research Letters, 51(2), e2023GL106201.

% "aq1_01_20.nc" 
% Link: https://doi.pangaea.de/10.1594/PANGAEA.924857
% Stål, T., et al. (2021). Antarctic geothermal heat flow model: Aq1. Geochemistry, Geophysics, Geosystems, 22(2), e2020GC009428.

% "/ICESat1_ICESat2_mass_change_updated_2_2021/" 
% Link: http://hdl.handle.net/1773/45388
% Smith, B., et al. (2020). Pervasive ice sheet mass loss reflects competing ocean and atmosphere processes. Science, 368(6496), 1239-1242.

function make_antarctic_figures()

% Load parameter values (in order to non-dimensionalise)
[scales, params] = default_parameter_values();
xdkm = scales.xdkm;

[x_g, y_g] = antbounds_data('gl','xy');
x_g = x_g/1000;
y_g = y_g/1000;

% Values of storage / permeability
Ss_L = 1e-5;
K_sb_L = 1e-6;
kappa = scales.Ksb_to_kappa*K_sb_L;
Sigma = scales.Ss_to_Sigma*Ss_L;

St = params.St;

% 20km grid (for sedimentary basin thickness data from Li et al.)
XL = linspace(-3.33e6,3.33e6,334);
YL = linspace(-3.33e6,3.33e6,334);
[XXL, YYL] = meshgrid(XL,YL);

% Grid for GHF model (from Stål et al.) 
XG = linspace(-2.79e6,2.79e6,280);
YG = linspace(-2.79e6,2.79e6,280);
[XXG, YYG] = meshgrid(XG,YG);

% Grid for frozen bed mask (from Seiner et al.)
XS = linspace(-3.04e6,3.04e6,761);
YS = linspace(-3.04e6,3.04e6,761);
[XXS, YYS] = meshgrid(XS,YS);
frozen_mask = ncread("AIS_BaTh_v1.nc",'binaryIsmip');

XX = XXG;
YY = YYG;


% Reading, smoothing and regridding basin thickness 
load("Ant_Crust.mat","MeanSSB");
Hsb = smoothdata2(MeanSSB,"gaussian",5);
Hsb = interp2(XXL,YYL,Hsb,XX,YY);

% Reading frozen bed mask
frozen_mask = 1-round(interp2(XXS,YYS,frozen_mask,XX,YY))';

% Reading GHF model
G = ncread('aq1_01_20.nc','Q')';

% Mask to focus on grounded ice
grounded_mask = isgrounded(XX,YY);

% Reading Bedmap3 ice bed and thickness data (from Pritchard et al., 2025)
[b,~] = readgeoraster('Bedmap3/bm3_bed.tif');
[Hi,~] = readgeoraster('Bedmap3/bm3_thickness.tif');

% Regridding and smoothing data
Hi = flip(interp2(-3333250:500:3333250,-3333250:500:3333250,double(Hi),XX,YY),1);
b = flip(interp2(-3333250:500:3333250,-3333250:500:3333250,double(b),XX,YY),1);
Hi(Hi<0)= NaN;
b = smoothdata2(b,"gaussian",5);
Hi = smoothdata2(Hi,"gaussian",5);

% Reading and regridding mass change data
BB = geotiffinfo('ais_dmdt_grounded_filt.tif').BoundingBox;
dHidt = flip(readgeoraster('ais_dmdt_grounded_filt.tif'));
[meshX, meshY] = meshgrid(linspace(BB(1,1), BB(2,1), size(dHidt,2)),linspace(BB(1,2), BB(2,2), size(dHidt,1)));
dHidt = interp2(meshX,meshY,dHidt,XX,YY);


% Non-dimensionalising variables
G_ = G/scales.Gd; 
XX_ = XX/(1000*scales.xdkm);
YY_ = YY/(1000*scales.xdkm); 
Hsb_ = Hsb./scales.Hd;
b_ = b./scales.Hd; 
dHidt_ = dHidt/(scales.Hd/scales.tdyr); 
dHidt_(isnan(dHidt_))=0;
dHidt_(~grounded_mask) = nan;

% Boundary limits
xlims = xdkm*[-5.5 5.5];
ylims = xdkm*[-4.5 4.5];

% Figure 1: (a) GHF model and (b) sedimentary basin thickness data
figure(1); clf; set(gcf(),'color','w','Position',[370 620 1400 500]); clf; 
tiledlayout(1,2,'TileSpacing','Compact','Padding','compact');

% 1(a)
nexttile;
aG = set_up_axes(gca(),'thermal','(a)');
title(aG,'GHF $G$','FontSize',14,'Interpreter','latex')
make_nice_colorbar(aG,'EastOutside','$G$ / mW m$^{-2}$');
contourf(aG,xdkm*XX_,xdkm*YY_,1000*scales.Gd*G_,20:5:180,'EdgeColor','none');
clim(aG,[0 180]);
[~,fr] = contourf(aG,xdkm*XX_,xdkm*YY_,frozen_mask,[1 1],'FaceColor','none','EdgeColor','w','LineWidth',0.5);
hatchfill2(fr,'HatchDensity',120);

% 1(b)
nexttile; gca();
asb = set_up_axes(gca(),'tempo','(b)');
title(asb,'Sedimentary basin thickness $H_{sb}$','FontSize',14,'Interpreter','latex')
make_nice_colorbar(asb,'EastOutside','$H_{sb}$ / m');
contourf(asb,xdkm*XX_,xdkm*YY_,scales.Hd*Hsb_,0:250:6000,'EdgeColor','none');
text(asb,2400,2000,'(b)','BackgroundColor','w','FontSize',14,'Interpreter','latex')
[psix_, psiy_] = gradient(((1-params.delta)*Hi*0.9+b)/scales.Hd,20/scales.xdkm);
plot(asb,x_g, y_g, 'k', 'linewidth',1);
hatch_frozen(asb);

% Calculation of heat flux modification and qE
Qx_ = 0.5*G_.*Hsb_.^2.*psix_;
Qy_ = 0.5*G_.*Hsb_.^2.*psiy_;

Fx_ = Hsb_.*psix_;
Fy_ = Hsb_.*psiy_;

[Qxx_, ~] = gradient(Qx_,20/scales.xdkm);
[~, Qyy_] = gradient(Qy_,20/scales.xdkm);

[Fxx_, ~] = gradient(Fx_,20/scales.xdkm);
[~, Fyy_] = gradient(Fy_,20/scales.xdkm);

% Topographic
mod_T_ = (Qxx_ + Qyy_);
mod_T_(~grounded_mask)=nan;

qET_ = Fxx_ + Fyy_;
qET_(~grounded_mask)=nan;

% Compaction
qEC_ = -Hsb_.*dHidt_ ;
mod_C_ = -0.5*Hsb_.^2.*G_.*dHidt_ ;

% Figure 5: (a) Total heat flux modification compared with (b) G (c)
% topographic term (d) compaction term
figure(5); clf; set(gcf(),'color','w','Position',[120 160 1200 1070]);
tiledlayout(2,2,'TileSpacing','Compact','Padding','compact');
nexttile; ab1 = set_up_axes(gca(),'balance','(a)');
nexttile; ab4 = set_up_axes(gca(),'thermal','(b)');
nexttile; ab2 = set_up_axes(gca(),'balance','(c)');
nexttile; ab3 = set_up_axes(gca(),'balance','(d)');

pcolor(ab1,xdkm*XX_,xdkm*YY_,1000*scales.Gd*St*(kappa*mod_T_+Sigma*mod_C_),'EdgeColor','none','FaceColor','interp');
pcolor(ab2,xdkm*XX_,xdkm*YY_,1000*scales.Gd*(St*kappa*mod_T_),'EdgeColor','none','FaceColor','interp');
pcolor(ab3,xdkm*XX_,xdkm*YY_,1000*scales.Gd*(St*Sigma*mod_C_),'EdgeColor','none','FaceColor','interp');

max_G = 60;
clim(ab1,[-max_G max_G]);    hatch_frozen(ab1);
clim(ab2,[-max_G max_G]);    hatch_frozen(ab2);
clim(ab3,[-max_G max_G]);    hatch_frozen(ab3);
make_nice_colorbar(ab1,'NorthOutside','mW m$^{-2}$');

add_location_labels(ab1,xlims,ylims);
add_location_labels(ab2,xlims,ylims);
add_location_labels(ab3,xlims,ylims);

title(ab1,'Total','FontSize',16,'Interpreter','latex','FontWeight','bold')
title(ab2,'Topographic','FontSize',16,'Interpreter','latex','FontWeight','bold')
title(ab3,'Compaction','FontSize',16,'Interpreter','latex','FontWeight','bold')

% Drawing circle around ICESAT data gap
rectangle(ab3,'Position',xdkm*[-.86 -.86 1.76 1.76],'Curvature',1,'LineStyle','--','linewidth',1.5)
rectangle(ab1,'Position',xdkm*[-.86 -.86 1.76 1.76],'Curvature',1,'LineStyle','--','linewidth',1.5)

% 5(d)
title(ab4,'Unmodified GHF $G$','FontSize',14,'Interpreter','latex')
make_nice_colorbar(ab4,'NorthOutside','$G$ / mW m$^{-2}$');
contourf(ab4,xdkm*XX_,xdkm*YY_,1000*scales.Gd*G_,20:5:180,'EdgeColor','none');
clim(ab4,[0 180]);
[~,fr] = contourf(ab4,xdkm*XX_,xdkm*YY_,frozen_mask,[1 1],'FaceColor','none','EdgeColor','w','LineWidth',0.5);
hatchfill2(fr,'HatchDensity',120);
add_location_labels(ab4,xlims,ylims);

% Figure 8: (a) Melt modification compared to (b) qE
figure(8); clf; set(gcf(),'color','w','Position',[70 590 1200 590]);
tiledlayout(1,2,'TileSpacing','Compact','Padding','compact');
nexttile; bm = set_up_axes(gca(),'balance','(a)');
nexttile; bqe = set_up_axes(gca(),'balance','(b)');

pcolor(bm,xdkm*XX_,xdkm*YY_,scales.mdmmyr*(St*kappa*mod_T_ + St*Sigma*mod_C_),'EdgeColor','none','FaceColor','interp');
pcolor(bqe,xdkm*XX_,xdkm*YY_,scales.mdmmyr*(kappa*qET_ + Sigma*qEC_),'EdgeColor','none','FaceColor','interp');

clim(bm,[-5 5]);    hatch_frozen(bm);
clim(bqe,[-15 15]);    hatch_frozen(bqe);
make_nice_colorbar(bm,'SouthOutside','mm yr$^{-1}$');
make_nice_colorbar(bqe,'SouthOutside','mm yr$^{-1}$');

title(bm,'Modification to melt rate $\Delta m$','FontSize',16,'Interpreter','latex','FontWeight','bold');
title(bqe,'Exfiltration $q_E$','FontSize',16,'Interpreter','latex','FontWeight','bold');

% Drawing circle around ICESAT data gap
rectangle(bm,'Position',xdkm*[-.86 -.86 1.76 1.76],'Curvature',1,'LineStyle','--','linewidth',1.5)
rectangle(bqe,'Position',xdkm*[-.86 -.86 1.76 1.76],'Curvature',1,'LineStyle','--','linewidth',1.5)

add_location_labels(bm,xlims,ylims);
add_location_labels(bqe,xlims,ylims);

% Figure 12: Lambda = c_p * G * Hsb / k_sb * L
figure(12); clf; set(gcf(),'color','w','Position',[400 600 700 500]); clf;
aLam = set_up_axes(gca(),'tempo','');

pcolor(aLam,xdkm*XX_,xdkm*YY_,params.St*G_.*Hsb_,'EdgeColor','none','FaceColor','interp');
clim(aLam,[0 1.5]);    
make_nice_colorbar(aLam,'EastOutside','');
hatch_frozen(aLam);
title(aLam,'$\Lambda = c_p G H_{sb} / k_{sb} L$','FontSize',16,'Interpreter','latex','FontWeight','bold');

% Reading MEaSUREs velocity data
skip = 5;
xv = -2800000:450*skip:2800000; 
yv = (-2800000:450*skip:2800000)'; 

[xxv, yyv] = meshgrid(xv/1000,yv/1000);

vx = ncread('antarctica_ice_velocity_450m_v2.nc','VX');
vx = rot90(vx(1:skip:end,1:skip:end)); 
vy = ncread('antarctica_ice_velocity_450m_v2.nc','VY');
vy = rot90(vy(1:skip:end,1:skip:end)); 

% Figure 6: regional plots
% 6(c) WAIS
figure(2); clf; set(gcf(),'color','w','Position',[1 1 1600 600]);     
regional_plots([-1800 -200],[-1300 300],'(c)');

% 6(b) Coats land
figure(3); clf; set(gcf(),'color','w','Position',[1 1 1600 600]);     
regional_plots([-700 500],[400 1600],'(b)');

% 6(a) Wilkes / Aurora 
figure(4); clf; set(gcf(),'color','w','Position',[1 1 1600 600]);     
regional_plots([800 2400],[-2000 -400],'(a)');

% Figure 11: Ice thinning rate data
figure(11); clf; set(gcf(),'color','w','Position',[400 600 700 500]); clf;
aHt = set_up_axes(gca(),'-balance','');
title(aHt,'Ice thickness rate of change $\partial H_i / \partial t$','FontSize',14,'Interpreter','latex')
make_nice_colorbar(aHt,'EastOutside','$\partial H_i / \partial t$ / m yr$^{-1}$');
pcolor(aHt,xdkm*XX_,xdkm*YY_,dHidt,'EdgeColor','none');
clim(aHt,[-5 5]);

[~,fr] = contourf(aHt,xdkm*XX_,xdkm*YY_,frozen_mask,[1 1],'FaceColor','none','EdgeColor','k','LineWidth',0.5);
hatchfill2(fr,'HatchDensity',120);

rectangle(aHt,'Position',xdkm*[-.86 -.86 1.76 1.76],'Curvature',1,'LineStyle','--','linewidth',1.5)

% Function to initialise axes for map plots
    function ax = set_up_axes(ax,clrmap,subplot_label)
        
        % Custom function for consistent style
        make_nice_axes(ax,'','');
        set(ax,'FontSize',12);
        
        % Make geographically accurate
        axis equal 

        % Axis limits
        xlim(ax,xlims);
        ylim(ax,ylims);

        % Axis label
        xlabel(ax,'$x$ / km','FontSize',12,'Interpreter','latex');
        ylabel(ax,'$y$ / km','FontSize',12,'Interpreter','latex');

        % Specified colormap
        colormap(ax,cmocean(clrmap))

        % Add subplot label
        text(ax,2400,2000,subplot_label,'BackgroundColor','w','FontSize',14,'Interpreter','latex')
    end

% Function to hatch regions of frozen bed in black
    function hatch_frozen(ax)
        [~,fr] = contourf(ax,xdkm*XX_,xdkm*YY_,frozen_mask,[1 1],'FaceColor','none','EdgeColor','k','LineWidth',0.5);
        hatchfill2(fr,'HatchDensity',120);
    end

% Function to set up regional plots (Figure 6)
    function regional_plots(xlimsR,ylimsR,letter)

        % Creation and initalisation of axes
        tiledlayout(2,4,'TileSpacing','Compact','Padding','compact');
        nexttile(1,[2 2]); amr = make_nice_axes(gca(),'$x$ / km', '$y$ / km', 'GHF modification'); axis equal;
        nexttile; aGr = make_nice_axes(gca(),'','', '$G$'); axis equal;
        nexttile; aHsbr = make_nice_axes(gca(),'','', '$H_{sb}$'); axis equal;
        nexttile; abr = make_nice_axes(gca(),'','',  '$b$'); axis equal;
        nexttile; avr = make_nice_axes(gca(),'','', '$|u_i|$'); axis equal;

        set(amr,'FontSize',14);
        set(aGr,'FontSize',14);
        set(aHsbr,'FontSize',14);
        set(abr,'FontSize',14);
        set(avr,'FontSize',14);
    
        xlim(amr,xlimsR);    ylim(amr,ylimsR);
        xlim(aGr,xlimsR);    ylim(aGr,ylimsR);
        xlim(aHsbr,xlimsR);    ylim(aHsbr,ylimsR);
        xlim(abr,xlimsR);    ylim(abr,ylimsR);
        xlim(avr,xlimsR);    ylim(avr,ylimsR);

        colormap(amr,cmocean('balance'));
        colormap(aGr,cmocean('thermal'));
        colormap(aHsbr,cmocean('tempo'));
        colormap(abr,cmocean('-tempo'));
        colormap(avr,cmocean('thermal'));

        % Plotting various data
        pcolor(amr,xdkm*XX_,xdkm*YY_,1000*scales.Gd*(St*kappa*mod_T_),'EdgeColor','none','FaceColor','interp');
        pcolor(aGr,xdkm*XX_,xdkm*YY_,1000*scales.Gd*G_,'EdgeColor','none','FaceColor','interp');
        pcolor(aHsbr,xdkm*XX_,xdkm*YY_,scales.Hd*Hsb_,'EdgeColor','none','FaceColor','interp');
        pcolor(abr,xdkm*XX_,xdkm*YY_,scales.Hd*b_,'EdgeColor','none','FaceColor','interp');
        pcolor(avr,xxv,yyv,sqrt(vx.^2 + vy.^2),'EdgeColor','none','FaceColor','interp');
        set(avr,'Colorscale','log');
 
        % Adding colorbars
        make_nice_colorbar(amr,'EastOutside','mW m$^{-2}$');
        make_nice_colorbar(aGr,'EastOutside','mW m$^{-2}$');
        make_nice_colorbar(aHsbr,'EastOutside','m');
        make_nice_colorbar(abr,'EastOutside','m');
        make_nice_colorbar(avr,'EastOutside','m yr$^{-1}$');

        % Adding grounding line
        plot(aHsbr,x_g, y_g, 'k', 'linewidth',1);
        plot(avr,x_g, y_g, 'k', 'linewidth',1);
        plot(amr,x_g, y_g, 'k', 'linewidth',1);
        plot(aGr,x_g, y_g, 'k', 'linewidth',1);
        plot(abr,x_g, y_g, 'k', 'linewidth',1);

        % Colorbar saturation and hatching frozen region
        clim(amr,[-max_G max_G]);    hatch_frozen(amr);
        clim(aHsbr,[0 5000]);       hatch_frozen(aHsbr);

        clim(abr,[-1500 1000]);    
        hatch_frozen(aGr);
        hatch_frozen(abr);
        hatch_frozen(avr);

        % Adding subplot labels
        lettx = xlimsR(1) + 0.06*(xlimsR(2)-xlimsR(1));
        letty = ylimsR(1) + 0.9*(ylimsR(2)-ylimsR(1));

        text(amr,(xlimsR(1) + 0.03*(xlimsR(2)-xlimsR(1))),(ylimsR(1) + 0.95*(ylimsR(2)-ylimsR(1))),[letter '(i)'],'BackgroundColor','w','FontSize',14,'Interpreter','latex');
        text(aGr,lettx,letty,[letter '(ii)'],'BackgroundColor','w','FontSize',14,'Interpreter','latex');
        text(aHsbr,lettx,letty,[letter '(iii)'],'BackgroundColor','w','FontSize',14,'Interpreter','latex');
        text(abr,lettx,letty,[letter '(iv)'],'BackgroundColor','w','FontSize',14,'Interpreter','latex');
        text(avr,lettx,letty,[letter '(v)'],'BackgroundColor','w','FontSize',14,'Interpreter','latex');
    end

% Function to label selected locations on a map
    function add_location_labels(ax,xlims,ylims)
        % Names to lookup locations with 'scarloc'
        locnames = {'recovery glacier','wilkes subglacial basin',...
            'aurora subglacial basin','siple coast','pine island glacier',...
            'thwaites glacier','berkner island','lake vostok'};
        % Abbreviated names to plot
        locabbrv = {'RG','WSB','ASB','SC','PIG','TG','BI','LV'};
        
        % Locations for text so they don't clutter map 
        textloc = [-1200 1600; 2000 -2000; 2200 -1300; -400 -1000; -2200 -1100 ; -1600 -1500;-1300 1200; 2200 600];

        % Add text label and draw a line from label to location
        for l = 1:length(locnames)
            loc = scarloc(locnames{l},'xy','km');
            if isbetween(loc(1),xlims(1),xlims(2)) && isbetween(loc(2),ylims(1),ylims(2))
                plot(ax, [loc(1) textloc(l,1)], [loc(2) textloc(l,2)] ,'k','linewidth',1);
                text(ax,textloc(l,1),textloc(l,2),locabbrv{l},'FontSize',16,'Interpreter','latex'); 
            end
        end
    end

end