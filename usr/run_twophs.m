% Prepare workspace
clear; close all;

% Set Run control and I/O parameters
runID    =  'twophs';            % run identifier
outdir   = '../out/';            % output directory
srcdir   = '../src/';            % source directory
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      =  10;                  % output frame plotted/saved every 'nop' time steps
plot_op  =  1;                   % switch on to live plot results
plot_crc =  0;                   % switch on particle and melt halo circles on plots
save_op  =  1;                   % switch on to save output to file
colourmap = 'lapaz';             % choose colourmap

load ../src/colmap/lapaz.mat;    % load colormap
typeclrs = lapaz([24,120,240],:);% select phase colours

% Set model domain parameters
D        =  1;                   % chamber depth [m]
L        =  D;                   % chamber width [m]
N        =  200;                 % number of grid points in z-direction

% Set model timing parameters
Nstep    =  1e5;                 % number of time steps to take
tend     =  1e3;                 % end time for simulation [s]
dt       =  1e-3;                % initial time step [s]
CFL      =  1.00;                % (physical) time stepping Courant-Friedrich-Levy number (multiplies stable step) [0,1]

% Set physical model parameters
grav     =  10;                  % gravity [m/s2]
ptol     =  0.9;                 % particle ordering tolerance (0-1, larger for more regular distribution)
seed     =  15;                  % random number generator seed for reproducibility
Nt       =  2;                   % number of particle types
rp       =  [0.01,0.01];         % particle radius [m]
fp       =  [0.01,0.01];         % particle fraction [vol]
rhop     =  [3200,2600];         % particle density [kg/m3]
strp     =  {'olv','plg','mlt'}; % name strings for particle and fluid phases
rhom     =  2900;                % melt density [kg/m3]
etap     =  1e6;                 % particle-melt viscosity contrast [Pas]
etam     =  1e1;                 % melt viscosity [Pas]

% Set numerical model parameters
TINT     =  'bd2im';             % time integration scheme ('be1im','bd2im','cn2si','bd2si')
ADVN     =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
rtol     =  1e-3;                % relative nonlinear residual tolerance
atol     =  1e-7;                % absolute nonlinear residual tolerance
maxit    =  6;                   % maximum outer its
alpha    =  0.95;                % iterative update step size parameter

%*****  RUN DESCARTES MODEL  **********************************************
run([srcdir,'main'])
%**************************************************************************

