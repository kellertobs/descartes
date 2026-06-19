% *****  UPDATE PARAMETERS & AUXILIARY FIELDS  ****************************

if step > 0

    % ---------------------------------------------------------------------
    % ADVECT CONCENTRATION FIELD
    % ---------------------------------------------------------------------

    % Advect particle fields
    advn = -advect(C, U(:, :), W(:, :), h, {ADVN, ''}, [1, 2], BCA);
   
    src  = max(0,Cq-C)/2/dt + min(0,Cq-C)/1000/dt;

    advn = advn + src;

    % Residual of particle field update
    res_C = (a1 * C - a2 * Co - a3 * Coo) / dt - (b1 * advn + b2 * dCdto + b3 * dCdtoo);

    % Semi-implicit update of major component density
    upd_C = - alpha * res_C * dt / a1;

    C = max(0,min(1, C + upd_C ));

    % ---------------------------------------------------------------------
    % UPDATE PARTICLE POSITIONS
    % ---------------------------------------------------------------------
    % Residuals for particle positions
    res_zp = (a1*zp - a2*zpo - a3*zpoo)/dt - (b1*Wp + b2*Wpo + b3*Wpoo);
    res_xp = (a1*xp - a2*xpo - a3*xpoo)/dt - (b1*Up + b2*Upo + b3*Upoo);

    % Semi-implicit update of particle positions
    upd_xp = - alpha*res_xp*dt/a1;
    xp     = xp + upd_xp;

    upd_zp = - alpha*res_zp*dt/a1;
    zp     = zp + upd_zp;

    % Apply periodic boundary conditions
    [xpoo(xp>L), xpo(xp>L), xp(xp>L)] = deal(xpoo(xp>L)-L, xpo(xp>L)-L, xp(xp>L)-L);
    [xpoo(xp<0), xpo(xp<0), xp(xp<0)] = deal(xpoo(xp<0)+L, xpo(xp<0)+L, xp(xp<0)+L);
    [zpoo(zp>D), zpo(zp>D), zp(zp>D)] = deal(zpoo(zp>D)-D, zpo(zp>D)-D, zp(zp>D)-D);
    [zpoo(zp<0), zpo(zp<0), zp(zp<0)] = deal(zpoo(zp<0)+D, zpo(zp<0)+D, zp(zp<0)+D);

end

% ---------------------------------------------------------------------
% PROCESS PARTICLES AT NEW POSITIONS
% ---------------------------------------------------------------------

[fp,Wp,Up,Vp,Pp,Wm,Um,Vm,Pm,DWp,DUp,DVp,DPp,indp,indm,Cq,chi] = process_part(Wc,Uc,P,xp,zp,rp,rm,tp,XX,ZZ,L,D,kp2);


% -------------------------------------------------------------------------
% UPDATE COEFFICIENT FIELDS (ETA, RHO)
% -------------------------------------------------------------------------

% Update eta (viscosity field)
eta = etam .* ones(Nz, Nx);
for it = 1:Nt
    eta = eta .* etap.^indp(:,:,it);
end

% Update rho (density field)
rho = rhom .* ones(Nz, Nx);
for it = 1:Nt
    rho = rho + (rhop(it) - rhom) .* indp(:,:,it);
end

% Compute averaged rho values at cell faces
rhofz = (rho(icz(1:end-1), :) + rho(icz(2:end), :)) / 2;
rhofx = (rho(:, icx(1:end-1)) + rho(:, icx(2:end))) / 2;

% Compute rho*W and rho*U for inertial terms
rhoW = rhofz .* W;
rhoU = rhofx .* U;

% Update eta at cell centers and edges
etacc = eta;
etaco = (eta(icz(1:end-1), icx(1:end-1)) .* eta(icz(1:end-1), icx(2:end)) ...
      .* eta(icz(2:end  ), icx(1:end-1)) .* eta(icz(2:end  ), icx(2:end))).^0.25;


% -------------------------------------------------------------------------
% UPDATE TIME STEP
% -------------------------------------------------------------------------
% Compute time step size based on CFL condition
dta = h / (2 * max(abs([U(:); W(:)])));
dt  = min([1.1 * dto, CFL * dta]);


% -------------------------------------------------------------------------
% UPDATE STRAIN RATES AND STRESSES
% -------------------------------------------------------------------------

% Strain rates
ups = ddz(W(:, :), h) + ddx(U(:, :), h);  % Velocity divergence

exx = ddx(U(:, :),h) - ups / 3;  % x-normal strain rate
ezz = ddz(W(:, :),h) - ups / 3;  % z-normal strain rate
exz = 0.5 * (ddz(U(icz,:),h) + ddx(W(:,icx),h));  % Shear strain rate

% Effective strain rate
eII = sqrt(0.5 * (exx.^2 + ezz.^2 + 2 * mean([exz(1:end-1, 1:end-1).^2, ...
                                              exz(2:end  , 1:end-1).^2, ...
                                              exz(1:end-1, 2:end  ).^2, ...
                                              exz(2:end  , 2:end  ).^2], 'all'))) + eps;

txx = etacc .* exx;  % x-normal stress
tzz = etacc .* ezz;  % z-normal stress
txz = etaco .* exz;  % xz-shear stress

% Effective stress
tII = sqrt(0.5 * (txx.^2 + tzz.^2 + 2 * mean([txz(1:end-1, 1:end-1).^2, ...
                                              txz(2:end  , 1:end-1).^2, ...
                                              txz(1:end-1, 2:end  ).^2, ...
                                              txz(2:end  , 2:end  ).^2], 'all'))) + eps;

% -------------------------------------------------------------------------
% UPDATE REYNOLDS NUMBER
% -------------------------------------------------------------------------

Re = Vel .* rho .* D ./ eta;  % Reynolds number calculation


% -------------------------------------------------------------------------
% HELPER FUNCTIONS
% -------------------------------------------------------------------------

function [fp,Wp,Up,Vp,Pp,Wm,Um,Vm,Pm,DWp,DUp,DVp,DPp,indp,indm,Cq,chi] = process_part(Wc,Uc,P,xp,zp,rp,rm,tp,XX,ZZ,L,D,kp2)

% % update particle indicator fields
% [Nz,Nx] = size(Wc);
% h       = D/Nz;
% Nt      = length(Np);
% indp    = zeros(Nz, Nx, Np);
% indm    = zeros(Nz, Nx, Np);
% for ip = 1:Np
%     dx = mod(XX - xp(ip) + L/2, L) - L/2;
%     dz = mod(ZZ - zp(ip) + D/2, D) - D/2;
% 
%     d = sqrt(dx.^2 + dz.^2);
% 
%     phi = rp(tp(ip)) - d;    % positive inside particle
%     indp(:,:,ip) = 0.5*(1 + erf(phi/(sqrt(2)*h/2)));
%     phi = rv(tp(ip)) - d;    % positive inside particle
%     indm(:,:,ip) = 0.5*(1 + erf(phi/(sqrt(2)*h/2)));
% end
% indm = max(0,indm - indp);

% update particle indicator fields
Nt      = length(rp);
Np      = length(xp);
[Nz,Nx] = size(Wc);
h       = D/Nz;
sgm     = h/4;

% precompute particle-specific radii
rpart = rp(tp);
rhalo = rm(tp);

indp = zeros(Nz,Nx,Nt);
indm = zeros(Nz,Nx,Nt);
Cq   = zeros(Nz,Nx,Nt);

Wp   = zeros(Np, 1);  Up = zeros(Np, 1);  Vp = zeros(Np, 1);  Pp = zeros(Np, 1);
Wm   = zeros(Np, 1);  Um = zeros(Np, 1);  Vm = zeros(Np, 1);  Pm = zeros(Np, 1);
fp   = zeros(Np, 1);

for ip = 1:Np

    % truncate stencil where erf contribution is negligible
    Rsupp = rhalo(ip) + 5*sgm;

    n = ceil(Rsupp/h);

    % nearest grid indices to particle centre
    ix0 = floor(xp(ip)/h) + 1;
    iz0 = floor(zp(ip)/h) + 1;

    % periodic stencil indices
    ix = mod((ix0-n):(ix0+n)-1,Nx) + 1;
    iz = mod((iz0-n):(iz0+n)-1,Nz) + 1;

    % local coordinates
    Xloc = XX(iz,ix);
    Zloc = ZZ(iz,ix);

    % exact periodic distances
    dx = mod(Xloc - xp(ip) + L/2, L) - L/2;
    dz = mod(Zloc - zp(ip) + D/2, D) - D/2;

    d = hypot(dx,dz);

    % smooth particle indicator
    wp = 0.5*(1 + erf((rpart(ip) - d)/(sqrt(2)*sgm)));

    % smooth halo indicator
    wm = 0.5*(1 + erf((rhalo(ip) - d)/(sqrt(2)*sgm)));
    wm = max(0, wm - wp);  % remove particle from melt halo

    % smooth concentration field
    wc = (1 + erf(min(0,rpart(ip) - d)/(sqrt(2)*sgm*4)));

    % get particle and melt halo metrics
    Wp(ip) = sum(Wc(iz,ix).*wp,'all') ./ sum(wp,'all');
    Up(ip) = sum(Uc(iz,ix).*wp,'all') ./ sum(wp,'all');
    Vp(ip) = sqrt(Wp(ip).^2 + Up(ip).^2);
    Pp(ip) = sum(P (iz,ix).*wp,'all') ./ sum(wp,'all');
    fp(ip) = sum(wp,'all')./sum(ones(size(Xloc)),'all');
    Wm(ip) = sum(Wc(iz,ix).*wm,'all') ./ sum(wm,'all');
    Um(ip) = sum(Uc(iz,ix).*wm,'all') ./ sum(wm,'all');
    Vm(ip) = sqrt(Wm(ip).^2 + Um(ip).^2);
    Pm(ip) = sum(P (iz,ix).*wm,'all') ./ sum(wm,'all');

    % accumulate global union fields
    indp(iz,ix,tp(ip)) = max(indp(iz,ix,tp(ip)), wp);
    indm(iz,ix,tp(ip)) = max(indm(iz,ix,tp(ip)), wm);
    Cq  (iz,ix,tp(ip)) = max(Cq  (iz,ix,tp(ip)), wc);

end

indm = max(0, indm - sum(indp,3));
Cq   = Cq.*(1-sum(indp,3)+indp);

% Particle segregation speeds
DWp = Wp - Wm;
DUp = Up - Um;
DVp = sqrt(DWp.^2 + DUp.^2);
DPp = Pp - Pm;

% Smoothed particle fraction
chi = zeros(Nz,Nx,Nt);

for it=1:Nt
    rr = 10*rp(it);
    Gk = exp(-rr^2/2 * kp2);
    chi(:,:,it) = real(ifft2(Gk .* fft2(indp(:,:,it))));
end

end
