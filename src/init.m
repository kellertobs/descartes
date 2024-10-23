% Model initialization script for Descartes numerical modelling code
% This script sets up the computational grid, initializes variables,
% and handles file IO for saving and restarting.

% Create output directory if it doesn't exist
output_dir = fullfile(outdir, runID);
if ~isfolder(output_dir)
    mkdir(output_dir);
end

% Save input parameters and runtime options (unless restarting)
if restart == 0 && save_op == 1
    parfile = fullfile(output_dir, [runID '_par']);
    save(parfile);
end

% Print initialization message with current timestamp and runID
fprintf('\n\n');
fprintf('*************************************************************\n');
fprintf('*****  RUN DESCARTES MODEL | %s  **********\n', datetime('now'));
fprintf('*************************************************************\n');
fprintf('\n   run ID: %s \n\n', runID);

% Load custom colormap for visualization
load ocean;

% Set boundary conditions for advection (periodic by default)
BCA = {'periodic', 'periodic'};  % Boundary condition on advection (top/bot, sides)
BCD = {'periodic', 'periodic'};  % Boundary condition for diffusion (top/bot, sides)

% Generate coordinate arrays for both cell-centered and face-centered grids
h  = D/N;        % grid spacing (equal in both dimensions) [m]
Xc = -h/2 : h : L + h/2;  % Cell-centered in X
Zc = -h/2 : h : D + h/2;  % Cell-centered in Z
Xf = (Xc(1:end-1) + Xc(2:end)) / 2;  % Face-centered in X
Zf = (Zc(1:end-1) + Zc(2:end)) / 2;  % Face-centered in Z

% Remove ghost points from grid for internal calculations
Xc = Xc(2:end-1);
Zc = Zc(2:end-1);
[XX, ZZ] = meshgrid(Xc, Zc);

% Grid sizes
Nx = length(Xc); % grid size in x-direction
Nz = length(Zc); % grid size in z-direction

% Set time unit conversion factors
hr       =  3600;                % seconds per hour
yr       =  hr*24*365.25;        % seconds per year

% Set random seed for reproducibility
rng(seed);

% Generate mapping arrays for various fields
NP = (Nz + 0) * (Nx + 0);  % Pressure field size (with ghost cells)
NW = (Nz + 1) * (Nx + 0);  % W-velocity field size
NU = (Nz + 0) * (Nx + 1);  % U-velocity field size
NC = Nz * Nx;              % Cell-centered field size
MapP = reshape(1:NP, Nz, Nx);  % Mapping for Pressure
MapW = reshape(1:NW, Nz+1, Nx);  % Mapping for W-velocity
MapU = reshape(1:NU, Nz, Nx+1) + NW;  % Mapping for U-velocity

% Set ghosted index arrays for periodic boundary conditions
icx = [Nx, 1:Nx, 1];  % Index for periodic X boundary
icz = [Nz, 1:Nz, 1];  % Index for periodic Z boundary
ifx = [Nx, 1:Nx+1, 2];  % Index for face-centered X boundary
ifz = [Nz, 1:Nz+1, 2];  % Index for face-centered Z boundary

MapSten = zeros(Nz*Nx,9);
k = 0;
for i=2:Nz+1
    for j=2:Nx+1
        k=k+1;
        MapSten(k,:) = [MapP(icz(i  ),icx(j  )), ...
                        MapP(icz(i-1),icx(j  )), MapP(icz(i  ),icx(j-1)), ...
                        MapP(icz(i+1),icx(j  )), MapP(icz(i  ),icx(j+1)), ...
                        MapP(icz(i-1),icx(j-1)), MapP(icz(i-1),icx(j+1)), ...
                        MapP(icz(i+1),icx(j-1)), MapP(icz(i+1),icx(j+1))];
    end
end

% Initialize fluid mechanics solution fields (U, W, P)
U = zeros(Nz, Nx+1); upd_U = zeros(size(U));
W = zeros(Nz+1, Nx); upd_W = zeros(size(W));
P = zeros(Nz, Nx); Vel = zeros(Nz, Nx); upd_P = zeros(size(P));  res_P = 0*P;
SOL = [W(:); U(:); P(:)];

% Initialize particle fields
[xp, zp, tp, Np] = generate_particles(Nt, rp, fp, D, L, pord);

% Initialize particle velocity arrays
Wp = zeros(Np, 1); Wpo = Wp; Wm = Wp;
Up = zeros(Np, 1); Upo = Up; Um = Up;

% Initialize auxiliary parameters
xpo = xp; zpo = zp;
upd_xp = xp; upd_zp = zp;
dto = dt;
a1 = 1; a2 = 1; a3 = 0;
b1 = 1; b2 = 0; b3 = 0;
frst = 1;
step = 0;
time = 0;
iter = 0;
dsumMdt = 0; dsumMdto = 0;
dsumCdt = 0; dsumCdto = 0;
HST = [];

% Initialize coefficient fields
eta = 0*P + etam;
update;
C     = Cq;
dCdt  = zeros(size(C)); dCdto = dCdt; dCdtoo = dCdto;
upd_C = zeros(size(C));

% Store the initial fields
Co = C; Coo = Co;
rhoWo = rhoW;
rhoUo = rhoU;
store;

% Store initial concentration and particle fields
C0 = C;
xp0 = xp;
zp0 = zp;
indp0 = indp;

% Overwrite fields from file if restarting the run
if restart
    handle_restart();
else
    % Complete, plot, and save the initial condition
    fluidmech;
    update;
    output;
end

restart = 0;


%% Helper Functions

function handle_restart()
    % Handle restart logic based on the restart flag
    if restart < 0  % Restart from the last continuation frame
        name = fullfile(outdir, runID, [runID '_cont.mat']);
    elseif restart > 0  % Restart from a specified continuation frame
        name = fullfile(outdir, runID, [runID '_' num2str(restart) '.mat']);
    end
    
    if exist(name, 'file')
        fprintf('\n   Restarting from %s \n\n', name);
        load(name, 'U', 'W', 'P', 'C', 'rho', 'eta', 'HST', 'dCdt', 'eII', 'tII', 'dt', 'time', 'step');
        name = fullfile(outdir, runID, [runID '_HST']);
        load(name, 'HST');

        SOL = [W(:); U(:); P(:)];

        update;
        Co = C;
        dCdto = dCdt;
        dto = dt;

        fluidmech;
        update;
        output;

        time = time + dt;
        step = step + 1;
    else
        fprintf('\n   !!! Restart file does not exist !!!\n   Starting run from scratch: %s \n\n', runID);
        fluidmech;
        update;
        history;
        output;
    end
end

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
                    min_dist = tol * (rp(t) + rp(pt(j))/2);  % Minimum distance between particles (scaled by tolerance)

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
                        tol      = (1-1e-6)*tol;
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

    Np = length(px);

end