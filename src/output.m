
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

cmap = [colororder;[0 0 0]];

% set axis and border dimensions
axh = 8.00*sqrt(D/L); axw = 8.00*sqrt(L/D)+1.50;
ahs = 0.80; avs = 0.80;
axb = 1.25; axt = 1.50;
axl = 1.50; axr = 0.50;

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
imagesc(Xsc,Zsc,rho); axis ij equal tight; box on; cb = colorbar;
cols = ['k','r','b','g'];
hold on;
for it = 1:Nt
    idx = find(tp == it);
    viscircles([xp(idx)   zp(idx)], rp(it)  ,'Color',cmap(it,:),'LineWidth',1.0);
    viscircles([xp(idx)   zp(idx)], rp(it)*4,'Color',[0 0 0],'LineStyle',':','LineWidth',0.5,'EnhanceVisibility',0);
    viscircles([xp(idx)-L zp(idx)], rp(it)  ,'Color',cmap(it,:),'LineWidth',1.0);
    viscircles([xp(idx)-L zp(idx)], rp(it)*4,'Color',[0 0 0],'LineStyle',':','LineWidth',0.5,'EnhanceVisibility',0);
    viscircles([xp(idx)+L zp(idx)], rp(it)  ,'Color',cmap(it,:),'LineWidth',1.0);
    viscircles([xp(idx)+L zp(idx)], rp(it)*4,'Color',[0 0 0],'LineStyle',':','LineWidth',0.5,'EnhanceVisibility',0);
    viscircles([xp(idx) zp(idx)-D], rp(it)  ,'Color',cmap(it,:),'LineWidth',1.0);
    viscircles([xp(idx) zp(idx)-D], rp(it)*4,'Color',[0 0 0],'LineStyle',':','LineWidth',0.5,'EnhanceVisibility',0);
    viscircles([xp(idx) zp(idx)+D], rp(it)  ,'Color',cmap(it,:),'LineWidth',1.0);
    viscircles([xp(idx) zp(idx)+D], rp(it)*4,'Color',[0 0 0],'LineStyle',':','LineWidth',0.5,'EnhanceVisibility',0);
end
xlim([0 L]);
ylim([0 D]);
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\rho$ [kg/m$^3$]'],TX{:},FS{:}); xlabel(['Width [',SpaceUnits,']'],TX{:},FS{:}); set(gca,'YTickLabel',[]);

if step>0

% initialize figures and axes
if ~exist('fh2','var'); fh2 = figure(VIS{:});
else; set(0, 'CurrentFigure', fh2); clf;
end
colormap(ocean);
fh = axb + 1.5*axh + 0*avs + axt;
fw = axl + 2.0*axw + 0*ahs + axr;
set(fh2,UN{:},'Position',[10 10 fw fh]);
set(fh2,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(fh2,'Color','w','InvertHardcopy','off','Resize','off');
ax(21) = axes(UN{:},'position',[axl+0*2*axw+0*ahs axb+0*1.5*axh+0*avs 2*axw 1.5*axh]);

% plot velocity-pressure solution in Fig. 1
set(0,'CurrentFigure',fh2)
set(fh2,'CurrentAxes',ax(21));
for it = 1:Nt
plot(HST.time./TimeScale,-HST.DWp_NM(:,it)./SpeedScale,'-' ,'Color',cmap(it,:),'LineWidth',1.5); axis tight; box on; hold on
plot(HST.time./TimeScale,-(HST.DWp_NM(:,it)+[2,-2].*HST.DWp_std(:,it))./SpeedScale,':' ,'Color',cmap(it,:),'LineWidth',1.5);
plot(HST.time./TimeScale,-HST.DWp_EM(:,it)./SpeedScale,'--','Color',cmap(it,:),'LineWidth',1.5);
end
set(gca,TL{:},TS{:}); 
xlabel(['Time [',TimeUnits,']'],TX{:},FS{:});
ylabel(['Speed [',SpeedUnits,']'],TX{:},FS{:});

% initialize figures and axes
if ~exist('fh3','var'); fh3 = figure(VIS{:});
else; set(0, 'CurrentFigure', fh3); clf;
end
colormap(ocean);
fh = axb + 1.5*axh + 0*avs + axt;
fw = axl + 2.0*axw + 0*ahs + axr;
set(fh3,UN{:},'Position',[15 15 fw fh]);
set(fh3,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(fh3,'Color','w','InvertHardcopy','off','Resize','off');
ax(31) = axes(UN{:},'position',[axl+0*2*axw+0*ahs axb+0*1.5*axh+0*avs 2*axw 1.5*axh]);

% plot velocity-pressure solution in Fig. 1
set(0,'CurrentFigure',fh3)
set(fh3,'CurrentAxes',ax(31));
wbin = (max(Wp(:))-min(Wp(:)))/SpeedScale/30;
for it=1:Nt
    histogram(-Wp(tp==it)./SpeedScale,'BinWidth',wbin,'FaceColor',cmap(it,:),'FaceAlpha',0.75); axis tight; box on; hold on
    histogram(-Wm(tp==it)./SpeedScale,'BinWidth',wbin,'FaceColor',cmap(it,:),'FaceAlpha',0.25); axis tight; box on; hold on
end
set(gca,TL{:},TS{:}); 
xlabel(['Speed [',SpeedUnits,']'],TX{:},FS{:});
ylabel('Bin Count [1]',TX{:},FS{:});

end

drawnow

% save output to file
if save_op && ~restart
    name = [outdir,'/',runID,'/',runID,'_sol_',num2str(floor(step/nop))];
    print(fh1,name,'-dpng','-r300','-image');
    if step>0
    name = [outdir,'/',runID,'/',runID,'_sgr_',num2str(floor(step/nop))];
    print(fh2,name,'-dpng','-r300','-image');
    name = [outdir,'/',runID,'/',runID,'_hst_',num2str(floor(step/nop))];
    print(fh3,name,'-dpng','-r300','-image');
    end

    name = [outdir,'/',runID,'/',runID,'_',num2str(floor(step/nop))];
    save(name,'U','W','P','C','rho','eta','dCdt','eII','tII','dt','time','step');
    name = [outdir,'/',runID,'/',runID,'_cont'];
    save(name,'U','W','P','C','rho','eta','dCdt','eII','tII','dt','time','step');
    name = [outdir,'/',runID,'/',runID,'_HST'];
    save(name,'HST');

end

if save_op && (step==0 || restart)
    logfile = [outdir,'/',runID,'/',runID,'.log'];
    if exist(logfile,'file') && step==0; delete(logfile); end
    diary(logfile)
end
    