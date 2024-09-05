
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
if time < 1e3*hr
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
Xfsc = Xcf./SpaceScale;
Zfsc = Zcf./SpaceScale;

% set axis and border dimensions
axh = 6.00*sqrt(D/L); axw = 6.00*sqrt(L/D)+1.50;
ahs = 0.60; avs = 0.80;
axb = 1.20; axt = 1.50;
axl = 1.50; axr = 0.50;

% initialize figures and axes
if ~exist('fh1','var'); fh1 = figure(VIS{:});
else; set(0, 'CurrentFigure', fh1); clf;
end
colormap(ocean);
fh = axb + 2*axh + 1*avs + axt;
fw = axl + 2*axw + 1*ahs + axr;
set(fh1,UN{:},'Position',[1 1 fw fh]);
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
imagesc(Xfsc,Zfsc,rho); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\nabla \cdot \mathbf{v}$ [1/',TimeUnits,']'],TX{:},FS{:}); xlabel(['Width [',SpaceUnits,']'],TX{:},FS{:}); set(gca,'YTickLabel',[]);

drawnow

% % save output to file
% if save_op && ~restart
%     if Nx <= 10 && Nz <= 10  % print 0D plots
%         name = [outdir,'/',runID,'/',runID,'_tch_',num2str(floor(step/nop))];
%         print(fh1,name,'-dpng','-r300','-image');
%         name = [outdir,'/',runID,'/',runID,'_aux_',num2str(floor(step/nop))];
%         print(fh2,name,'-dpng','-r300','-image');
%         name = [outdir,'/',runID,'/',runID,'_cmp',num2str(floor(step/nop))];
%         print(fh3,name,'-dpng','-r300','-image');
%         name = [outdir,'/',runID,'/',runID,'_TAS',num2str(floor(step/nop))];
%         print(fh11,name,'-dpng','-r300','-image');
%         name = [outdir,'/',runID,'/',runID,'_AFM',num2str(floor(step/nop))];
%         print(fh12,name,'-dpng','-r300','-image');
%     elseif Nx <= 10  % create 1D plots
%         name = [outdir,'/',runID,'/',runID,'_sol_',num2str(floor(step/nop))];
%         print(fh1,name,'-dpng','-r300','-image');
%         name = [outdir,'/',runID,'/',runID,'_aux_',num2str(floor(step/nop))];
%         print(fh2,name,'-dpng','-r300','-image');
%         name = [outdir,'/',runID,'/',runID,'_cmp_',num2str(floor(step/nop))];
%         print(fh3,name,'-dpng','-r300','-image');
%         name = [outdir,'/',runID,'/',runID,'_TAS',num2str(floor(step/nop))];
%         print(fh11,name,'-dpng','-r300','-image');
%         name = [outdir,'/',runID,'/',runID,'_AFM',num2str(floor(step/nop))];
%         print(fh12,name,'-dpng','-r300','-image');
%     else
%         name = [outdir,'/',runID,'/',runID,'_vep_',num2str(floor(step/nop))];
%         print(fh1,name,'-dpng','-r300','-image');
%         name = [outdir,'/',runID,'/',runID,'_tch_',num2str(floor(step/nop))];
%         print(fh2,name,'-dpng','-r300','-image');
%         name = [outdir,'/',runID,'/',runID,'_phs_',num2str(floor(step/nop))];
%         print(fh3,name,'-dpng','-r300','-image');
%         name = [outdir,'/',runID,'/',runID,'_sgr_',num2str(floor(step/nop))];
%         print(fh4,name,'-dpng','-r300','-image');
%         name = [outdir,'/',runID,'/',runID,'_cmp',num2str(floor(step/nop))];
%         print(fh5,name,'-dpng','-r300','-image');
%         name = [outdir,'/',runID,'/',runID,'_oxd',num2str(floor(step/nop))];
%         print(fh6,name,'-dpng','-r300','-image');
%         name = [outdir,'/',runID,'/',runID,'_mnr',num2str(floor(step/nop))];
%         print(fh7,name,'-dpng','-r300','-image');
%         name = [outdir,'/',runID,'/',runID,'_gch',num2str(floor(step/nop))];
%         print(fh8,name,'-dpng','-r300','-image');
%         name = [outdir,'/',runID,'/',runID,'_TAS',num2str(floor(step/nop))];
%         print(fh11,name,'-dpng','-r300','-image');
%         name = [outdir,'/',runID,'/',runID,'_AFM',num2str(floor(step/nop))];
%         print(fh12,name,'-dpng','-r300','-image');
%     end
% 
%     name = [outdir,'/',runID,'/',runID,'_',num2str(floor(step/nop))];
%     save(name,'U','W','P','Pt','Pchmb','f','x','m','fq','xq','mq','phi','chi','mu','X','F','M','S','C','T','Tp','c','cm','cx','cf','TRC','trc','dSdt','dCdt','dFdt','dXdt','dMdt','drhodt','dTRCdt','Gf','Gx','Gm','rho','eta','eII','tII','dt','time','step','VolSrc','wf','wx','wm','cal');
%     name = [outdir,'/',runID,'/',runID,'_cont'];
%     save(name,'U','W','P','Pt','Pchmb','f','x','m','fq','xq','mq','phi','chi','mu','X','F','M','S','C','T','Tp','c','cm','cx','cf','TRC','trc','dSdt','dCdt','dFdt','dXdt','dMdt','drhodt','dTRCdt','Gf','Gx','Gm','rho','eta','eII','tII','dt','time','step','VolSrc','wf','wx','wm','cal');
%     name = [outdir,'/',runID,'/',runID,'_hist'];
%     save(name,'hist');
% 
% end
% 
% if save_op && (step==0 || restart)
%     logfile = [outdir,'/',runID,'/',runID,'.log'];
%     if exist(logfile,'file') && step==0; delete(logfile); end
%     diary(logfile)
% end
    