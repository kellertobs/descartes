% prepare workspace
clear; close all;

% load default parameters
% run('./par_default')

% set run parameters
runID    =  'demo';              % run identifier
outdir   = '../out/';
srcdir   = '../src/';
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      =  100;                 % output frame plotted/saved every 'nop' time steps
plot_op  =  1;                   % switch on to live plot results
save_op  =  1;                   % switch on to save output to file
periodic =  1;

% set model domain parameters
D        =  1;                   % chamber depth [m]
N        =  300;                 % number of grid points in z-direction
h        =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]
L        =  D;                   % chamber width (equal to h for 1-D mode) [m]

% set model timing parameters
Nt       =  5e5;                 % number of time steps to take
tend     =  1;                   % end time for simulation [s]
dt       =  36;                  % initial time step [s]

% set physical model parameters
seed     =  15;
Np       =  2;                   % number of particle types
rp       =  [1e-2,1e-2];         % particle radius [m]
fp       =  [0.1,0.1];           % target particle fraction [vol]
rhop     =  [2600,3300];         % particle density [kg/m3]
rhom     =  2700;                % melt density [kg/m3]
etap     =  1e6;                 % particle viscosity contrast rel. to melt [Pas]
etam     =  1e2;                 % melt viscosity
grav     =  10;                  % gravity [m/s2]

% set numerical model parameters
TINT     =  'bd2im';             % time integration scheme ('be1im','bd2im','cn2si','bd2si')
ADVN     =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
CFL      =  0.50;                % (physical) time stepping courant number (multiplies stable step) [0,1]
rtol     =  1e-4;                % outer its relative tolerance
atol     =  1e-8;                % outer its absolute tolerance
maxit    =  50;                  % maximum outer its

%*****  RUN NAKHLA MODEL  *************************************************
run([srcdir,'main'])
%**************************************************************************

