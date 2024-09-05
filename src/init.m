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
[XXco,ZZco] = meshgrid(Xf,Zf);
Xc        = Xc(2:end-1);
Zc        = Zc(2:end-1);
[XX,ZZ]   = meshgrid(Xc,Zc);

hf        = h/4;
Xcf       = hf/2:hf:L-hf/2;
Zcf       = hf/2:hf:D-hf/2;
[XXf,ZZf] = meshgrid(Xcf,Zcf);
Xcfg      = -hf/2:hf:L+hf/2;
Zcfg      = -hf/2:hf:D+hf/2;
[XXfg,ZZfg] = meshgrid(Xcfg,Zcfg);
Xff        = (Xcfg(1:end-1)+Xcfg(2:end))./2;
Zff        = (Zcfg(1:end-1)+Zcfg(2:end))./2;
[XXuf,ZZuf] = meshgrid(Xff,Zcfg);
[XXwf,ZZwf] = meshgrid(Xcfg,Zff);

Nx = length(Xc);  Nxf = length(Xcf);
Nz = length(Zc);  Nzf = length(Zcf);

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
NC =  Nzf   *  Nxf  ;
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
[xp, zp, tp] = generate_particles(Np, rp, fp, D, L, ptol);

C = zeros(Nzf,Nxf,Np);
for ip = 1:length(xp)
    C(:,:,tp(ip)) = min(1,C(:,:,tp(ip)) + double(sqrt((XXf-xp(ip)).^2 + (ZZf-zp(ip)).^2) < rp(tp(ip))));
end

% initialise auxiliary parameters 
Co = C;  Coo = Co;
dCdt   = 0.*C;  dCdto  = dCdt;  dCdtoo = dCdto;
upd_C  = 0.*C;
dto = dt;
a1 = 1; a2 = 1; a3 = 0;
b1 = 1; b2 = 0; b3 = 0;

% initialise coefficient fields
update;

rhoW = rhofz.*W(2:end-1,2:end-1); rhoWo = rhoW; rhoWoo = rhoWo;
rhoU = rhofx.*U(2:end-1,2:end-1); rhoUo = rhoU; rhoUoo = rhoUo;

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
    % history;
    output;
    step = step+1;
end

restart = 0;

function [px, pz, pt] = generate_particles(Np, rp, fp, D, L, tol)
    % Inputs:
    % Np  - number of particle types (e.g. 2 for two particle types)
    % rp  - vector of particle radii [rp1, rp2, ...]
    % fp  - vector of area fractions [fp1, fp2, ...] for each particle type
    % D   - width of the rectangular domain
    % L   - length of the rectangular domain
    % tol - tolerance for the minimum allowed spacing between particles (percentage of particle radius)
    
    % Outputs:
    % x_coords, y_coords - arrays of particle coordinates
    % particle_type - array denoting the type of each particle
    
    % Calculate total domain area
    domain_area = D * L;
    
    % Initialize variables for particle storage
    px = [];
    pz = [];
    pt = [];
    
    for t = 1:Np
        % Area to be occupied by particle type t
        target_area = fp(t) * domain_area;
        % Calculate number of particles of type t based on area fraction
        num_particles = round(target_area / (pi * rp(t)^2));
        
        for i = 1:num_particles
            placed = false;
            while ~placed
                % Randomly place particle inside the domain
                x = rand() * L;
                y = rand() * D;
                
                % Check if the new particle overlaps with any existing particles
                min_dist = tol * (rp(t)+rp(pt));  % Minimum distance between particles (scaled by tolerance)
                overlaps = false;
                
                for j = 1:length(px)
                    dist = sqrt((x - px(j))^2 + (y - pz(j))^2);
                    if dist < min_dist
                        overlaps = true;
                        break;
                    end
                end
                
                if ~overlaps
                    % If no overlap, store the particle
                    px = [px; x];
                    pz = [pz; y];
                    pt = [pt; t];
                    placed = true;
                end
            end
        end
    end
    
    px = [px;px(px<0+max(rp))+L]; pz = [pz;pz(px<0+max(rp))]; pt = [pt;pt(px<0+max(rp))];
    px = [px;px(px>L-max(rp))-L]; pz = [pz;pz(px>L-max(rp))]; pt = [pt;pt(px>L-max(rp))];

    pz = [pz;pz(pz<0+max(rp))+D]; px = [px;px(pz<0+max(rp))]; pt = [pt;pt(pz<0+max(rp))];
    pz = [pz;pz(pz>D-max(rp))-D]; px = [px;px(pz>D-max(rp))]; pt = [pt;pt(pz>L-max(rp))];

    % Plot particles for visualization
    figure;
    cols = ['k','r','b','g'];
    hold on;
    for t = 1:Np
        idx = find(pt == t);
        viscircles([px(idx) pz(idx)], rp(t) * ones(length(idx), 1),'Color',cols(t));
    end
    axis ij equal tight; box on;
    xlim([0 L]);
    ylim([0 D]);
    title('Randomized Particle Distribution');
    hold off;

end
