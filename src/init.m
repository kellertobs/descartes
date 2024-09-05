% create output directory
if ~isfolder([outdir,'/',runID])
    mkdir([outdir,'/',runID]);
end

% save input parameters and runtime options (unless restarting)
if restart == 0 && save_op == 1
    parfile = [outdir,'/',runID,'/',runID,'_par'];
    save(parfile);
end

fprintf('\n\n')
fprintf('*************************************************************\n');
fprintf('*****  RUN DESCARTES MODEL | %s  *************\n',datetime('now'));
fprintf('*************************************************************\n');
fprintf('\n   run ID: %s \n\n',runID);

load ocean;                  % load custom colormap
if periodic % periodic sides
    BCA     =  {'periodic','periodic'};  % boundary condition on advection (top/bot, sides)
    BCD     =  {'periodic','periodic'};  % boundary condition on advection (top/bot, sides)
else % closed sides
    BCA     =  {'closed','closed'};  % boundary condition on advection (top/bot, sides)
    BCD     =  {'closed','closed'};  % boundary condition on advection (top/bot, sides) 
end

% get coordinate arrays
Xc        = -h/2:h:L+h/2;
Zc        = -h/2:h:D+h/2;
Xf        = (Xc(1:end-1)+Xc(2:end))./2;
Zf        = (Zc(1:end-1)+Zc(2:end))./2;
[XXu,ZZu] = meshgrid(Xf,Zc);
[XXw,ZZw] = meshgrid(Xc,Zf);
Xc        = Xc(2:end-1);
Zc        = Zc(2:end-1);
[XX,ZZ]   = meshgrid(Xc,Zc);

Nx = length(Xc);
Nz = length(Zc);

% get smoothed initialisation field
rng(seed);
% smth = smth*Nx*Nz*1e-4;
% rp   = randn(Nz,Nx);
% for i = 1:round(smth)
%     rp = rp + diffus(rp,1/8*ones(size(rp)),1,[1,2],BCD);
%     rp = rp - mean(mean(rp));
% end
% rp = rp./max(abs(rp(:)));

% get mapping arrays
NP = (Nz+2) * (Nx+2);
NW = (Nz+1) * (Nx+2);
NU = (Nz+2) * (Nx+1);
MapP = reshape(1:NP,Nz+2,Nx+2);
MapW = reshape(1:NW,Nz+1,Nx+2);
MapU = reshape(1:NU,Nz+2,Nx+1) + NW;

% set ghosted index arrays
if periodic  % periodic side boundaries
    icx = [Nx,1:Nx,1];
    icz = [1,1:Nz,Nz];
    ifx = [Nx,1:Nx+1,2];
    ifz = [2,1:Nz+1,Nz];
else         % closed side boundaries
    icx = [1,1:Nx,Nx];
    icz = [1,1:Nz,Nz];
    ifx = [2,1:Nx+1,Nx];
    ifz = [2,1:Nz+1,Nz];
end

% initialise fluid mechanics solution fields
U   =  zeros(Nz+2,Nx+1);  UBG = U; Ui = U; upd_U = 0*U;
W   =  zeros(Nz+1,Nx+2);  WBG = W; Wi = W; wf = 0.*W; wx = 0.*W; wm = 0.*W; upd_W = 0*W;
P   =  zeros(Nz+2,Nx+2);  Vel = 0.*P; upd_P = 0*P;
SOL = [W(:);U(:);P(:)];

% initialise particle fields
np = round(D*L*fp./(pi.*rp.^2));
nr = 1; ir = 1; jr = 1; d = 0;
for ip = 1:Np
    xp = rand(np(ip),1).*L;
    zp = rand(np(ip),1).*D;
    xp(xp>L) = xp(xp>L)-L;  xp(xp<0) = xp(xp<0)+L;
    zp(zp>D) = zp(zp>D)-D;  zp(zp<0) = zp(zp<0)+D;
    figure(1); clf; plot(xp,zp,'ko');
    while nr>0
        d = (sqrt((xp-xp').^2 + (zp-zp').^2));
        Fxp = mean((xp-xp')./max(1e-3,(xp-xp').^2)).'/1000; Fxp = Fxp-mean(Fxp(:));
        Fzp = mean((zp-zp')./max(1e-3,(zp-zp').^2)).'/1000; Fzp = Fzp-mean(Fzp(:));
        xp = xp - Fxp;
        zp = zp - Fzp;
        xp(xp>L) = xp(xp>L)-L;  xp(xp<0) = xp(xp<0)+L;
        zp(zp>D) = zp(zp>D)-D;  zp(zp<0) = zp(zp<0)+D;
        % d  = triu(sqrt((xp-xp').^2 + (zp-zp').^2));
        % [ir,jr] = find(d>0 & (d<=mean(d(:))/5 | d>=mean(d(:))*5));
        % nr = length(ir);
        plot(xp,zp,'ko');
    end

end

% initialise auxiliary variables 
dCdt   = 0.*C;  dCdto  = dCdt;
upd_C  = 0.*C;

% initialise timing and iterative parameters
frst    = 1;
step    = 0;
time    = 0;
iter    = 0;
hist    = [];
dsumMdt = 0; dsumMdto = 0;
dsumCdt = 0; dsumCdto = 0;

% overwrite fields from file if restarting run
if restart
    if     restart < 0  % restart from last continuation frame
        name = [outdir,'/',runID,'/',runID,'_cont.mat'];
    elseif restart > 0  % restart from specified continuation frame
        name = [outdir,'/',runID,'/',runID,'_',num2str(restart),'.mat'];
    end
    if exist(name,'file')
        fprintf('\n   restart from %s \n\n',name);
        load(name,'U','W','P','Pt','C','dCdt','rho','eta','eII','tII','dt','time','step');
        name = [outdir,'/',runID,'/',runID,'_hist'];
        load(name,'hist');

        SOL = [W(:);U(:);P(:)];

        update; 
        
        Co    = C;
        dCdto = dCdt;
        dto   = dt;

        fluidmech;
        update;
        output;

        time    = time+dt;
        step    = step+1;

    else % continuation file does not exist, start from scratch
        fprintf('\n   !!! restart file does not exist !!! \n   => starting run from scratch %s \n\n',runID);
        fluidmech;
        update;
        history;
        output;
    end
else
    % complete, plot, and save initial condition
    fluidmech;
    update;
    history;
    output;
    step = step+1;
end

restart = 0;
