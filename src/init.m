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
BCA     =  {'periodic','periodic'};  % boundary condition on advection (top/bot, sides)
BCD     =  {'periodic','periodic'};  % boundary condition on advection (top/bot, sides)

% get coordinate arrays
Xc        = -h/2:h:L+h/2;
Zc        = -h/2:h:D+h/2;
Xf        = (Xc(1:end-1)+Xc(2:end))./2;
Zf        = (Zc(1:end-1)+Zc(2:end))./2;
[XXu,ZZu]   = meshgrid(Xf,Zc);
[XXw,ZZw]   = meshgrid(Xc,Zf);
[XXco,ZZco] = meshgrid(Xf,Zf);
Xc        = Xc(2:end-1);
Zc        = Zc(2:end-1);
[XX,ZZ]   = meshgrid(Xc,Zc);

Nx = length(Xc);
Nz = length(Zc);

rng(seed);

% get mapping arrays
NP = (Nz+2) * (Nx+2);
NW = (Nz+1) * (Nx+2);
NU = (Nz+2) * (Nx+1);
NC =  Nz    *  Nx   ;
MapP = reshape(1:NP,Nz+2,Nx+2);
MapW = reshape(1:NW,Nz+1,Nx+2);
MapU = reshape(1:NU,Nz+2,Nx+1) + NW;

% set ghosted index arrays
icx = [Nx,1:Nx,1];
icz = [Nz,1:Nz,1];
ifx = [Nx,1:Nx+1,2];
ifz = [Nz,1:Nz+1,2];

% initialise fluid mechanics solution fields
U   =  zeros(Nz+2,Nx+1);  UBG = U; Ui = U; upd_U = 0*U;
W   =  zeros(Nz+1,Nx+2);  WBG = W; Wi = W; wf = 0.*W; wx = 0.*W; wm = 0.*W; upd_W = 0*W;
P   =  zeros(Nz+2,Nx+2);  Vel = 0.*P(2:end-1,2:end-1);  upd_P = 0*P;
SOL = [W(:);U(:);P(:)];

% initialise particle fields
[xp, zp, tp, Np] = generate_particles(Nt, rp, fp, D, L, ptol);

xpo = xp;
zpo = zp;
Wp = zeros(Np,1);  Wpo = Wp;  Wm = Wp;
Up = zeros(Np,1);  Upo = Up;  Um = Up;

C = zeros(Nz,Nx,Nt);
for ip = 1:Np
    C(:,:,tp(ip)) = min(1,C(:,:,tp(ip)) + double(sqrt((XX-xp(ip)).^2 + (ZZ-zp(ip)).^2) < rp(tp(ip))));
    C(:,:,tp(ip)) = min(1,C(:,:,tp(ip)) + double(sqrt((XX-xp(ip)-L).^2 + (ZZ-zp(ip)+0).^2) < rp(tp(ip))));
    C(:,:,tp(ip)) = min(1,C(:,:,tp(ip)) + double(sqrt((XX-xp(ip)+L).^2 + (ZZ-zp(ip)+0).^2) < rp(tp(ip))));
    C(:,:,tp(ip)) = min(1,C(:,:,tp(ip)) + double(sqrt((XX-xp(ip)+0).^2 + (ZZ-zp(ip)-D).^2) < rp(tp(ip))));
    C(:,:,tp(ip)) = min(1,C(:,:,tp(ip)) + double(sqrt((XX-xp(ip)+0).^2 + (ZZ-zp(ip)+D).^2) < rp(tp(ip))));
end

for i=1:2*min(rp/h)
    C = C + (diff(C(icz,:,:),2,1) + diff(C(:,icx,:),2,2))/8;
    C = C./max(max(max(C,[],1)),max(max(C,[],2)));
end
C = 2*C-1;

C0 = C;

% initialise auxiliary parameters 
Co = C;  Coo = Co;
dCdt  = 0.*C;  dCdto = dCdt;  dCdtoo = dCdto;
upd_C = 0.*C;
upd_xp = xp;
upd_zp = zp;
dto = dt;
a1 = 1; a2 = 1; a3 = 0;
b1 = 1; b2 = 0; b3 = 0;
frst    = 1;
step    = 0;
time    = 0;
iter    = 0;
dsumMdt = 0; dsumMdto = 0;
dsumCdt = 0; dsumCdto = 0;
HST     = [];

% initialise coefficient fields
update;

rhoWo = rhoW;
rhoUo = rhoU;

store;

% overwrite fields from file if restarting run
if restart
    if     restart < 0  % restart from last continuation frame
        name = [outdir,'/',runID,'/',runID,'_cont.mat'];
    elseif restart > 0  % restart from specified continuation frame
        name = [outdir,'/',runID,'/',runID,'_',num2str(restart),'.mat'];
    end
    if exist(name,'file')
        fprintf('\n   restart from %s \n\n',name);
        load(name,'U','W','P','C','rho','eta','HST','dCdt','eII','tII','dt','time','step');
        name = [outdir,'/',runID,'/',runID,'_HST'];
        load(name,'HST');

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
end

restart = 0;

function [px, pz, pt, Np] = generate_particles(Nt, rp, fp, D, L, tol)
    % Inputs:
    % Nt  - number of particle types (e.g. 2 for two particle types)
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

    for t = 1:Nt
        % Area to be occupied by particle type t
        target_area = fp(t) * domain_area;
        % Calculate number of particles of type t based on area fraction
        np = round(target_area / (pi * rp(t)^2));
        
        for i = 1:np
            placed = false;
            while ~placed
                % Randomly place particle inside the domain
                if np==1
                    x = L/2;
                    z = D/2;
                else
                    x = rand() * L;
                    z = rand() * D;
                end
                
                % Check if the new particle overlaps with any existing particles
                overlaps = false;
                
                for j = 1:length(px)
                    min_dist = tol * (rp(t) + rp(pt(j)));  % Minimum distance between particles (scaled by tolerance)
                    
                    % Calculate distance with periodic boundaries in both x and z directions
                    dx = abs(x - px(j));
                    dz = abs(z - pz(j));
                    
                    % Apply periodic boundary conditions
                    dx = min(dx, L - dx);
                    dz = min(dz, D - dz);
                    
                    % Calculate the Euclidean distance considering periodicity
                    dist = sqrt(dx^2 + dz^2);

                    if any(dist(:) < min_dist)
                        overlaps = true;
                        break;
                    end
                end
                
                if ~overlaps
                    % If no overlap, store the particle
                    px = [px; x(1)];
                    pz = [pz; z(1)];
                    pt = [pt; t   ];
                    placed = true;
                end
            end
        end
    end

    % px = [px;px(px<0+max(rp))+L]; pz = [pz;pz(px<0+max(rp))]; pt = [pt;pt(px<0+max(rp))];
    % px = [px;px(px>L-max(rp))-L]; pz = [pz;pz(px>L-max(rp))]; pt = [pt;pt(px>L-max(rp))];
    % 
    % pz = [pz;pz(pz<0+max(rp))+D]; px = [px;px(pz<0+max(rp))]; pt = [pt;pt(pz<0+max(rp))];
    % pz = [pz;pz(pz>D-max(rp))-D]; px = [px;px(pz>D-max(rp))]; pt = [pt;pt(pz>L-max(rp))];

    Np = length(px);

    % % Plot particles for visualization
    % figure;
    % cols = ['k','r','b','g'];
    % hold on;
    % for t = 1:Nt
    %     idx = find(pt == t);
    %     viscircles([px(idx)   pz(idx)], rp(t) ,'Color',cols(t));
    %     viscircles([px(idx)-L pz(idx)], rp(t) ,'Color',cols(t));
    %     viscircles([px(idx)+L pz(idx)], rp(t) ,'Color',cols(t));
    %     viscircles([px(idx) pz(idx)-D], rp(t) ,'Color',cols(t));
    %     viscircles([px(idx) pz(idx)+D], rp(t) ,'Color',cols(t));
    % end
    % axis ij equal tight; box on;
    % xlim([0 L]);
    % ylim([0 D]);
    % title('Randomized Particle Distribution');
    % hold off;

end
