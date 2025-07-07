% Prepare workspace
clear; close all;

% Set Run control and I/O parameters
runID    =  'jupiter_trace';     % run identifier
outdir   = '../out/';            % output directory
srcdir   = '../src/';            % source directory
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      =  20;                  % output frame plotted/saved every 'nop' time steps
plot_op  =  1;                   % switch on to live plot results
plot_crc =  0;
save_op  =  1;                   % switch on to save output to file
colourmap = 'lapaz';
% typeclrs = [ 32  26 116; ...
%              50  96 154; ...
%             114 152 160;...
%             251 232 225]./255;
typeclrs = [ 20  30  85; ...
             40  80 120; ...
            100 140 160; ...
            175 150 140; ...
            170 120  95; ...
            155  75  60; ...
            215 210 200]./255;
% load ../src/colmap/ocean.mat;
% typeclrs = ocean([20,195,160,60],:).*1.1;
% typeclrs = [typeclrs;[0.95 0.90 0.85]];


% Set model domain parameters
D        =  1;                   % chamber depth [m]
L        =  D;                   % chamber width [m]
N        =  600;                 % number of grid points in z-direction

% Set model timing parameters
Nstep    =  1e5;                 % number of time steps to take
tend     =  1e3;                 % end time for simulation [s]
dt       =  1e-3;                % initial time step [s]
CFL      =  1.00;                % (physical) time stepping Courant-Friedrich-Levy number (multiplies stable step) [0,1]

% Set physical model parameters
grav     =  10;                              % gravity [m/s2]
pord     =  24;                              % particle ordering parameter (larger for more regular distribution)
seed     =  15;                               % random number generator seed for reproducibility
Nt       =  6;                               % number of particle types
rp       =  [0.015,0.02,0.01,0.01,0.02,0.015]; % particle radius [m]
fp       =  [0.02,0.03,0.01,0.01,0.03,0.02]; % particle fraction [vol]
rhop     =  [3400,3200,3000,2600,2400,2200]; % particle density [kg/m3]
strp     =  {'dbl','blu','trq','bei','amb','rst','gry'};       % name strings for particle and fluid phases
rhom     =  2800;                            % melt density [kg/m3]
etap     =  1e6;                             % particle-melt viscosity contrast [Pas]
etam     =  1e-2;                            % melt viscosity [Pas]
DW0      =  (rhop-rhom).*grav.*rp.^2./etam;  % Stokes particle settling speeds
Rep      =  DW0.*(rhop-rhom).*rp./etam;      % particle Reynolds numbers

% Set numerical model parameters
TINT     =  'bd2im';             % time integration scheme ('be1im','bd2im','cn2si','bd2si')
ADVN     =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
rtol     =  1e-4;                % relative nonlinear residual tolerance
atol     =  1e-7;                % absolute nonlinear residual tolerance
maxit    =  3;                   % maximum outer its
alpha    =  0.95;                % iterative update step size parameter
beta     =  0.00;                % iterative update damping parameter

%*****  RUN DESCARTES MODEL  **********************************************
run([srcdir,'main'])
%**************************************************************************

