% prepare workspace
clear; close all;

% load default parameters
% run('./par_default')

% set run parameters
runID    =  'demo';              % run identifier
outdir   = '../out/';
srcdir   = '../src/';
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      =  10;                 % output frame plotted/saved every 'nop' time steps
plot_op  =  1;                   % switch on to live plot results
save_op  =  1;                   % switch on to save output to file

% set model domain parameters
D        =  1;                   % chamber depth [m]
N        =  200;                 % number of grid points in z-direction
h        =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]
L        =  D;                   % chamber width (equal to h for 1-D mode) [m]

% set model timing parameters
Nt       =  1e5;                 % number of time steps to take
tend     =  1e6;                   % end time for simulation [s]
dt       =  0.1;                 % initial time step [s]
hr       =  3600;
yr       =  hr*24*365.25;

% set physical model parameters
seed     =  24;
ptol     =  1.6;
Np       =  2;                   % number of particle types
rp       =  [5*h,5*h];         % particle radius [m]
fp       =  [0.1,0.1];           % target particle fraction [vol]
rhop     =  [3300,2500];         % particle density [kg/m3]
rhom     =  2800;                % melt density [kg/m3]
etap     =  1e6;                 % particle viscosity contrast rel. to melt [Pas]
etam     =  1e0;                 % melt viscosity
grav     =  10;                  % gravity [m/s2]

% set numerical model parameters
TINT     =  'bd2si';             % time integration scheme ('be1im','bd2im','cn2si','bd2si')
ADVN     =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
CFL      =  0.50;                % (physical) time stepping courant number (multiplies stable step) [0,1]
rtol     =  1e-3;                % relative nonlinear residual tolerance
atol     =  1e-6;                % absolute nonlinear residual tolerance
maxit    =  6;                  % maximum outer its
alpha    =  0.95;
beta     =  0.0;

%*****  RUN NAKHLA MODEL  *************************************************
run([srcdir,'main'])
%**************************************************************************

