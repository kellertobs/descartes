
% prepare for plotting
TX = {'Interpreter','Latex'}; FS = {'FontSize',12};
TL = {'TickLabelInterpreter','Latex'}; TS = {'FontSize',10};
UN = {'Units','Centimeters'};
CL = {'Color',[0.0 0.0 0.0],[0.80 0.15 0.10],[0.10 0.15 0.65],[0.45 0.60 0.95]};
LW = {'LineWidth',2};
if plot_op
    VIS = {'Visible','on'};
else
    VIS = {'Visible','off'};
end

% adjust scales and units for intuitive visualisation
if time < hr
    TimeScale = 1;
    TimeUnits = 'sec';
elseif time >= hr && time < 1e3*hr
    TimeScale = hr;
    TimeUnits = 'hr';
elseif time >= 1e3*hr && time < 1e2*yr
    TimeScale = yr;
    TimeUnits = 'yr';
elseif time >= 1e2*yr
    TimeScale = 1e3*yr;
    TimeUnits = 'kyr';
end
if D < 1e3
    SpaceScale = 1;
    SpaceUnits = 'm';
elseif D >= 1e3
    SpaceScale = 1e3;
    SpaceUnits = 'km';
end
if max(Vel(:)) < 1000/yr
    SpeedScale = 1/yr;
    SpeedUnits = 'm/yr';
elseif max(Vel(:)) >= 1000/yr && max(Vel(:)) < 1000/hr
    SpeedScale = 1/hr;
    SpeedUnits = 'm/hr';
elseif max(Vel(:)) >= 1000/hr
    SpeedScale = 1;
    SpeedUnits = 'm/s';
end
Xsc  = Xc./SpaceScale;
Zsc  = Zc./SpaceScale;

cmap = ocean([20,195,160,60],:).*1.1;
cmap = [cmap;[0.95 0.90 0.85]];

% set axis and border dimensions
axh = 8.00*sqrt(D/L); axw = 8.00*sqrt(L/D)+1.50;
ahs = 0.80; avs = 0.80;
axb = 1.25; axt = 1.50;
axl = 1.75; axr = 0.50;

% initialize figures and axes
if ~exist('fh1','var'); fh1 = figure(VIS{:});
else; set(0, 'CurrentFigure', fh1); clf;
end
colormap(ocean);
fh = axb + 2*axh + 1*avs + axt;
fw = axl + 2*axw + 1*ahs + axr;
set(fh1,UN{:},'Position',[5 5 fw fh]);
set(fh1,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(fh1,'Color','w','InvertHardcopy','off','Resize','off');
ax(11) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
ax(12) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
ax(13) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
ax(14) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);

% plot velocity-pressure solution in Fig. 1
set(0,'CurrentFigure',fh1)
set(fh1,'CurrentAxes',ax(11));
imagesc(Xsc,Zsc,-W(:      ,2:end-1)./SpeedScale); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$W$ [',SpeedUnits,']'],TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel(['Depth [',SpaceUnits,']'],TX{:},FS{:});
set(fh1,'CurrentAxes',ax(12));
imagesc(Xsc,Zsc, U(2:end-1,:      )./SpeedScale); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$U$ [',SpeedUnits,']'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
text(-0.1,1.1,['time = ',num2str(time/TimeScale,3),' [',TimeUnits,']'],TX{:},FS{:},'Color','k','HorizontalAlignment','center','Units','normalized');
set(fh1,'CurrentAxes',ax(13));
imagesc(Xsc,Zsc, P(2:end-1,2:end-1)); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$P$ [Pa]'],TX{:},FS{:}); ylabel(['Depth [',SpaceUnits,']'],TX{:},FS{:}); xlabel(['Width [',SpaceUnits,']'],TX{:},FS{:});
set(fh1,'CurrentAxes',ax(14));
imagesc(Xsc,Zsc,rho(2:end-1,2:end-1)); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\rho$ [kg/m$^3$]'],TX{:},FS{:}); xlabel(['Width [',SpaceUnits,']'],TX{:},FS{:}); set(gca,'YTickLabel',[]);

if ~exist('fh2','var'); fh2 = figure(VIS{:});
else; set(0, 'CurrentFigure', fh2); clf;
end
colormap(ocean);
fh = axb + 2.0*axh + 0*avs + axt/2;
fw = axl + 2.0*axw + 0*ahs + axr;
set(fh2,UN{:},'Position',[10 10 fw fh]);
set(fh2,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(fh2,'Color','w','InvertHardcopy','off','Resize','off');
ax(21) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs 2*axw 2*axh]);

set(fh2,'CurrentAxes',ax(21));
CRGB = zeros(Nz,Nx,3);
for it=1:Nt
    CRGB = CRGB +    C(:,:,it) .*permute(repmat(cmap(it,:).',1,Nz,Nx),[2 3 1]);
    colmap(it,:) = cmap(it,:);
end
CRGB = CRGB + (1-sum(C,3)).*permute(repmat(cmap(end,:).',1,Nz,Nx),[2 3 1]);
colmap(it+1,:) = cmap(end,:);
imagesc(Xsc,Zsc,CRGB); axis ij equal tight; box on; cb = colorbar; colormap(ax(21),colmap);
cols = ['k','r','b','g'];
hold on;
rv = D./sqrt(Np)/2;
for it = 1:Nt
    idx = find(tp == it);
    viscircles([xp(idx)   zp(idx)], rp(it),'Color',cmap(it,:),'LineWidth',1.0);
    viscircles([xp(idx)   zp(idx)], rv(it),'Color',[0 0 0],'LineStyle',':','LineWidth',0.5,'EnhanceVisibility',0);
    viscircles([xp(idx)-L zp(idx)], rp(it),'Color',cmap(it,:),'LineWidth',1.0);
    viscircles([xp(idx)-L zp(idx)], rv(it),'Color',[0 0 0],'LineStyle',':','LineWidth',0.5,'EnhanceVisibility',0);
    viscircles([xp(idx)+L zp(idx)], rp(it),'Color',cmap(it,:),'LineWidth',1.0);
    viscircles([xp(idx)+L zp(idx)], rv(it),'Color',[0 0 0],'LineStyle',':','LineWidth',0.5,'EnhanceVisibility',0);
    viscircles([xp(idx) zp(idx)-D], rp(it),'Color',cmap(it,:),'LineWidth',1.0);
    viscircles([xp(idx) zp(idx)-D], rv(it),'Color',[0 0 0],'LineStyle',':','LineWidth',0.5,'EnhanceVisibility',0);
    viscircles([xp(idx) zp(idx)+D], rp(it),'Color',cmap(it,:),'LineWidth',1.0);
    viscircles([xp(idx) zp(idx)+D], rv(it),'Color',[0 0 0],'LineStyle',':','LineWidth',0.5,'EnhanceVisibility',0);
end
xlim([0 L]);
ylim([0 D]);
set(cb,TL{:},TS{:},'Ticks',linspace(1/2/(Nt+1),1-1/2/(Nt+1),Nt+1),'TickLabels',strp); set(gca,TL{:},TS{:}); 
title(['Phases'],TX{:},FS{:}); ylabel(['Depth [',SpaceUnits,']'],TX{:},FS{:}); xlabel(['Width [',SpaceUnits,']'],TX{:},FS{:});
text(0.8,1.02,['time = ',num2str(time/TimeScale,3),' [',TimeUnits,']'],TX{:},FS{:},'Color','k','HorizontalAlignment','center','Units','normalized');

if step>0

% initialize figures and axes
if ~exist('fh3','var'); fh3 = figure(VIS{:});
else; set(0, 'CurrentFigure', fh3); clf;
end
colormap(ocean);
fh = axb + 1.5*axh + 0*avs + axt/2;
fw = axl + 2.0*axw + 0*ahs + axr;
set(fh3,UN{:},'Position',[10 10 fw fh]);
set(fh3,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(fh3,'Color','w','InvertHardcopy','off','Resize','off');
ax(31) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs 2*axw 1.5*axh]);

% plot segregation speed history in Fig. 3
set(0,'CurrentFigure',fh3)
set(fh3,'CurrentAxes',ax(31));
ph(1) = plot(HST.time./TimeScale, HST.Wc_mean./SpeedScale,'-' ,'Color','k','LineWidth',2); axis tight; box on; hold on
ph(2) = plot(HST.time./TimeScale,(HST.Wc_mean+HST.Wc_std)./SpeedScale,':' ,'Color','k','LineWidth',1.5);
ph(2) = plot(HST.time./TimeScale,(HST.Wc_mean-HST.Wc_std)./SpeedScale,':' ,'Color','k','LineWidth',1.5);
ph(3) = plot(HST.time./TimeScale, HST.Wc_tavg./SpeedScale,'-.' ,'Color','k','LineWidth',1.5);
for it = 1:Nt
pph(it) = plot(HST.time./TimeScale,- HST.DWp_mean(:,it)./SpeedScale,'-' ,'Color',cmap(it,:),'LineWidth',2); axis tight; box on; hold on
plot(HST.time./TimeScale,-(HST.DWp_mean(:,it)+HST.DWp_std(:,it))./SpeedScale,':' ,'Color',cmap(it,:),'LineWidth',1.5);
plot(HST.time./TimeScale,-(HST.DWp_mean(:,it)-HST.DWp_std(:,it))./SpeedScale,':' ,'Color',cmap(it,:),'LineWidth',1.5);
plot(HST.time./TimeScale,- HST.DWp_tavg(:,it)./SpeedScale,'-.' ,'Color',cmap(it,:),'LineWidth',1.5);
end
maxspd = max([-(HST.DWp_mean(:)+HST.DWp_std(:));-(HST.DWp_mean(:)-HST.DWp_std(:));])./SpeedScale;
minspd = min([-(HST.DWp_mean(:)+HST.DWp_std(:));-(HST.DWp_mean(:)-HST.DWp_std(:));])./SpeedScale;
line([time/2/TimeScale,time/2/TimeScale],[minspd,maxspd],'Color','k','LineStyle','--','LineWidth',1)
set(gca,TL{:},TS{:}); 
legend([ph(1:3),pph(1:Nt)],[{'mean'},{'std'},{'time avg.'},strp(1:Nt)],TX{:},FS{:},'Location','southwest');
title(['Convection \& Segregation Speeds'],TX{:},FS{:}); 
xlabel(['Time [',TimeUnits,']'],TX{:},FS{:});
ylabel(['Speed [',SpeedUnits,']'],TX{:},FS{:});

% initialize figures and axes
if ~exist('fh4','var'); fh4 = figure(VIS{:});
else; set(0, 'CurrentFigure', fh4); clf;
end
colormap(ocean);
fh = axb + Nt*axh + (Nt-1)*avs + axt/2;
fw = axl + 1.5*axw + 0*ahs + axr;
set(fh4,UN{:},'Position',[15 15 fw fh]);
set(fh4,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(fh4,'Color','w','InvertHardcopy','off','Resize','off');
for it=1:Nt
ax(40+it) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+(it-1)*axh+(it-1)*avs 1.5*axw axh]);
end

% plot phase speed histograms in Fig. 4
set(0,'CurrentFigure',fh4)
wbin = (max([Wp(:);Wm(:)])-min([Wp(:);Wm(:)]))/SpeedScale/25;
limx = [min(-[Wp(:);Wm(:)]),max(-[Wp(:);Wm(:)])]/SpeedScale;
for it=1:Nt
    set(fh4,'CurrentAxes',ax(40+it));
    histogram(-Wp(tp==it)./SpeedScale,'BinWidth',wbin,'FaceColor',cmap(it,:),'FaceAlpha',0.7); axis tight; box on; hold on
    histogram(-Wm(tp==it)./SpeedScale,'BinWidth',wbin,'FaceColor',cmap(it,:),'FaceAlpha',0.2); axis tight; box on; hold on
    set(gca,TL{:},TS{:},'xlim',limx);
    if it==1; xlabel(['Speed [',SpeedUnits,']'],TX{:},FS{:}); end
    ylabel('Bin Count [1]',TX{:},FS{:});
end
sgtitle(['Phase Speed Distributions'],TX{:},FS{:}); 
text(0.8,1.03,['time = ',num2str(time/TimeScale,3),' [',TimeUnits,']'],TX{:},FS{:},'Color','k','HorizontalAlignment','center','Units','normalized');

end

drawnow

% save output to file
if save_op && ~restart
    name = [outdir,'/',runID,'/',runID,'_sol_',num2str(floor(step/nop))];
    print(fh1,name,'-dpng','-r300','-image');
    name = [outdir,'/',runID,'/',runID,'_phs_',num2str(floor(step/nop))];
    print(fh2,name,'-dpng','-r300','-image');
    if step>0
    name = [outdir,'/',runID,'/',runID,'_sgr'];
    print(fh3,name,'-dpng','-r300','-image');
    name = [outdir,'/',runID,'/',runID,'_hst_',num2str(floor(step/nop))];
    print(fh4,name,'-dpng','-r300','-image');
    end

    name = [outdir,'/',runID,'/',runID,'_',num2str(floor(step/nop))];
    save(name,'U','W','P','C','rho','eta','eII','tII','dt','time','step');
    name = [outdir,'/',runID,'/',runID,'_cont'];
    save(name,'U','W','P','C','rho','eta','eII','tII','dt','time','step');
    name = [outdir,'/',runID,'/',runID,'_HST'];
    save(name,'HST');

end

if save_op && (step==0 || restart)
    logfile = [outdir,'/',runID,'/',runID,'.log'];
    if exist(logfile,'file') && step==0; delete(logfile); end
    diary(logfile)
end
    