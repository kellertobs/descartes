% Prepare workspace
clear; close all;

% Set Run control and I/O parameters
runID    =  'demo';              % run identifier
outdir   = '../out/';            % output directory
srcdir   = '../src/';            % source directory
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      =  1;                   % output frame plotted/saved every 'nop' time steps
plot_op  =  1;                   % switch on to live plot results
save_op  =  0;                   % switch on to save output to file

% Set model domain parameters
D        =  1;                   % chamber depth [m]
L        =  D;                   % chamber width [m]
N        =  300;                 % number of grid points in z-direction

% Set model timing parameters
Nstep    =  1e5;                 % number of time steps to take
tend     =  1e3;                 % end time for simulation [s]
dt       =  1e-3;                % initial time step [s]
CFL      =  1.00;                % (physical) time stepping Courant-Friedrich-Levy number (multiplies stable step) [0,1]

% Set physical model parameters
grav     =  10;                              % gravity [m/s2]
pord     =  10;                              % particle ordering parameter (larger for more regular distribution)
seed     =  15;                              % random number generator seed for reproducibility
Nt       =  3;                               % number of particle types
rp       =  [0.03,0.02,0.01];                % particle radius [m]
fp       =  [0.02,0.03,0.01];                % particle fraction [vol]
rhop     =  [3200,2600,4200];                % particle density [kg/m3]
strp     =  {'pxn','plg','spn','mlt'};       % name strings for particle and fluid phases
rhom     =  2800;                            % melt density [kg/m3]
etap     =  1e6;                             % particle-melt viscosity contrast [Pas]
etam     =  1e1;                             % melt viscosity [Pas]
DW0      =  (rhop-rhom).*grav.*rp.^2./etam;  % Stokes particle settling speeds
Rep      =  DW0.*(rhop-rhom).*rp./etam;      % particle Reynolds numbers

% Set numerical model parameters
TINT     =  'bd2im';             % time integration scheme ('be1im','bd2im','cn2si','bd2si')
ADVN     =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
rtol     =  1e-4;                % relative nonlinear residual tolerance
atol     =  1e-7;                % absolute nonlinear residual tolerance
maxit    =  6;                   % maximum outer its
alpha    =  0.95;                % iterative update step size parameter
beta     =  0.00;                % iterative update damping parameter

%*****  RUN DESCARTES MODEL  **********************************************
run([srcdir,'main'])
%**************************************************************************

