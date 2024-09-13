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
N        =  400;                 % number of grid points in z-direction
h        =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]
L        =  D;                   % chamber width (equal to h for 1-D mode) [m]

% set model timing parameters
Nstep    =  1e5;                 % number of time steps to take
tend     =  1e6;                 % end time for simulation [s]
dt       =  0.01;                % initial time step [s]
hr       =  3600;
yr       =  hr*24*365.25;

% set physical model parameters
seed     =  24;
ptol     =  2.1;
Ns       =  20;                  % size of sampling volume for settling speed evaluation
Nt       =  2;                   % number of particle types
rp       =  [0.02,0.02];           % particle radius [m]
fp       =  [0.05,0.05];         % target particle fraction [vol]
rhop     =  [2500,3300];         % particle density [kg/m3]
rhom     =  2800;                % melt density [kg/m3]
etap     =  1e6;                 % particle viscosity contrast rel. to melt [Pas]
etam     =  1e2;                 % melt viscosity
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

