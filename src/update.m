% *****  UPDATE PARAMETERS & AUXILIARY FIELDS  ****************************

if step > 0

    % ---------------------------------------------------------------------
    % ADVECT CONCENTRATION FIELD
    % ---------------------------------------------------------------------

    % Advect particle fields
    advn = -advect(C, U(:, :), W(:, :), h, {ADVN, ''}, [1, 2], BCA);
   
    src  = max(0,Cq-C)/2/dt + min(0,Cq-C)/2000/dt;

    advn = advn + src;

    % Residual of particle field update
    res_C = (a1 * C - a2 * Co - a3 * Coo) / dt - (b1 * advn + b2 * dCdto + b3 * dCdtoo);

    % Semi-implicit update of major component density
    upd_C = - alpha * res_C * dt / a1 + beta * upd_C;

    C = max(0,min(1, C + upd_C ));

    % ---------------------------------------------------------------------
    % UPDATE PARTICLE POSITIONS
    % ---------------------------------------------------------------------
    % Residuals for particle positions
    res_zp = (a1*zp - a2*zpo - a3*zpoo)/dt - (b1*Wp + b2*Wpo + b3*Wpoo);
    res_xp = (a1*xp - a2*xpo - a3*xpoo)/dt - (b1*Up + b2*Upo + b3*Upoo);

    % Semi-implicit update of particle positions
    upd_xp = - alpha*res_xp*dt/a1 + beta*upd_xp;
    xp     = xp + upd_xp;

    upd_zp = - alpha*res_zp*dt/a1 + beta*upd_zp;
    zp     = zp + upd_zp;

    % Apply periodic boundary conditions
    [xpoo(xp>L), xpo(xp>L), xp(xp>L)] = deal(xpoo(xp>L)-L, xpo(xp>L)-L, xp(xp>L)-L);
    [xpoo(xp<0), xpo(xp<0), xp(xp<0)] = deal(xpoo(xp<0)+L, xpo(xp<0)+L, xp(xp<0)+L);
    [zpoo(zp>D), zpo(zp>D), zp(zp>D)] = deal(zpoo(zp>D)-D, zpo(zp>D)-D, zp(zp>D)-D);
    [zpoo(zp<0), zpo(zp<0), zp(zp<0)] = deal(zpoo(zp<0)+D, zpo(zp<0)+D, zp(zp<0)+D);

end

% ---------------------------------------------------------------------
% UPDATE PARTICLE VELOCITIES (Wp, Up) AND MELT VELOCITIES (Wm, Um)
% ---------------------------------------------------------------------

[fp,Wp,Up,Vp,Wm,Um,Vm,DWp,DUp,DVp,indp,indm] = phase_vel(Wc,Uc,xp,zp,rp,tp,Np,XX,ZZ,L,D);


% ---------------------------------------------------------------------
% UPDATE EQUILIBRIUM CONCENTRATIONS (Cq)
% ---------------------------------------------------------------------

Cq = get_conc(XX,ZZ,xp,zp,tp,rp,indp,L,D,Nz,Nx,Nt,Np);


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
dt = min([1.1 * dto, CFL * dta]);
% 
% dtauW = (h/2)^2 ./ max( eta(icz(1:end-1),:), eta(icz(2:end),:));
% dtauU = (h/2)^2 ./ max( eta(:,icx(1:end-1)), eta(:,icx(2:end)));
% dtauP = eta;

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

function [fp,Wp,Up,Vp,Wm,Um,Vm,DWp,DUp,DVp,indpt,indmt] = phase_vel(Wc,Uc,xp,zp,rp,tp,Np,XX,ZZ,L,D)

Wp = zeros(sum(Np), 1);  Up = zeros(sum(Np), 1);
Wm = zeros(sum(Np), 1);  Um = zeros(sum(Np), 1);
fp = zeros(sum(Np), 1);
rv = rp + D./sqrt(sum(Np))/2;

% update particle indicator fields
[Nz,Nx] = size(Wc);
Nt      = length(Np);
indp    = zeros(Nz, Nx, sum(Np));
indm    = zeros(Nz, Nx, sum(Np));
for ip = 1:sum(Np)
    % Update the particle indicator field for a given particle
    indp(:,:,ip) = min(1, indp(:,:,ip) + double(sqrt((XX - xp(ip)    ).^2 + (ZZ - zp(ip)    ).^2) < rp(tp(ip))));
    % Apply periodic boundary conditions in X and Z
    indp(:,:,ip) = min(1, indp(:,:,ip) + double(sqrt((XX - xp(ip) - L).^2 + (ZZ - zp(ip)    ).^2) < rp(tp(ip))));
    indp(:,:,ip) = min(1, indp(:,:,ip) + double(sqrt((XX - xp(ip) + L).^2 + (ZZ - zp(ip)    ).^2) < rp(tp(ip))));
    indp(:,:,ip) = min(1, indp(:,:,ip) + double(sqrt((XX - xp(ip)    ).^2 + (ZZ - zp(ip) - D).^2) < rp(tp(ip))));
    indp(:,:,ip) = min(1, indp(:,:,ip) + double(sqrt((XX - xp(ip)    ).^2 + (ZZ - zp(ip) + D).^2) < rp(tp(ip))));
    indp(:,:,ip) = min(1, indp(:,:,ip) + double(sqrt((XX - xp(ip) - L).^2 + (ZZ - zp(ip) - D).^2) < rp(tp(ip))));
    indp(:,:,ip) = min(1, indp(:,:,ip) + double(sqrt((XX - xp(ip) + L).^2 + (ZZ - zp(ip) - D).^2) < rp(tp(ip))));
    indp(:,:,ip) = min(1, indp(:,:,ip) + double(sqrt((XX - xp(ip) - L).^2 + (ZZ - zp(ip) + D).^2) < rp(tp(ip))));
    indp(:,:,ip) = min(1, indp(:,:,ip) + double(sqrt((XX - xp(ip) + L).^2 + (ZZ - zp(ip) + D).^2) < rp(tp(ip))));

    % Update the particle indicator field for a given particle
    indm(:,:,ip) = min(1, indm(:,:,ip) + double(sqrt((XX - xp(ip)    ).^2 + (ZZ - zp(ip)    ).^2) < rv(tp(ip))));
    % Apply periodic boundary conditions in X and Z
    indm(:,:,ip) = min(1, indm(:,:,ip) + double(sqrt((XX - xp(ip) - L).^2 + (ZZ - zp(ip)    ).^2) < rv(tp(ip))));
    indm(:,:,ip) = min(1, indm(:,:,ip) + double(sqrt((XX - xp(ip) + L).^2 + (ZZ - zp(ip)    ).^2) < rv(tp(ip))));
    indm(:,:,ip) = min(1, indm(:,:,ip) + double(sqrt((XX - xp(ip)    ).^2 + (ZZ - zp(ip) - D).^2) < rv(tp(ip))));
    indm(:,:,ip) = min(1, indm(:,:,ip) + double(sqrt((XX - xp(ip)    ).^2 + (ZZ - zp(ip) + D).^2) < rv(tp(ip))));
    indm(:,:,ip) = min(1, indm(:,:,ip) + double(sqrt((XX - xp(ip) - L).^2 + (ZZ - zp(ip) - D).^2) < rv(tp(ip))));
    indm(:,:,ip) = min(1, indm(:,:,ip) + double(sqrt((XX - xp(ip) + L).^2 + (ZZ - zp(ip) - D).^2) < rv(tp(ip))));
    indm(:,:,ip) = min(1, indm(:,:,ip) + double(sqrt((XX - xp(ip) - L).^2 + (ZZ - zp(ip) + D).^2) < rv(tp(ip))));
    indm(:,:,ip) = min(1, indm(:,:,ip) + double(sqrt((XX - xp(ip) + L).^2 + (ZZ - zp(ip) + D).^2) < rv(tp(ip))));
end
indm = max(0,indm - indp - sum(indp,3).*(1-indp));
indp = indp + (diff(indp(:,[end-1,1:end,2],:),2,2) + diff(indp([end-1,1:end,2],:,:),2,1))/8;
indm = indm + (diff(indm(:,[end-1,1:end,2],:),2,2) + diff(indm([end-1,1:end,2],:,:),2,1))/8;

for ip = 1:sum(Np)
    Wp(ip) = sum(Wc.*indp(:,:,ip),'all') ./ sum(indp(:,:,ip),'all');
    Up(ip) = sum(Uc.*indp(:,:,ip),'all') ./ sum(indp(:,:,ip),'all');
    Vp(ip) = sqrt(Wp(ip).^2 + Up(ip).^2);
    fp(ip) = sum(indp(:,:,ip),'all')./sum(ones(size(XX)),'all');
    Wm(ip) = sum(Wc.*indm(:,:,ip),'all') ./ sum(indm(:,:,ip),'all');
    Um(ip) = sum(Uc.*indm(:,:,ip),'all') ./ sum(indm(:,:,ip),'all');
    Vm(ip) = sqrt(Wm(ip).^2 + Um(ip).^2);
end

% Particle segregation speeds
DWp = Wp - Wm;
DUp = Up - Um;
DVp = sqrt(DWp.^2 + DUp.^2);

% Sum up particle indicator function by type
indpt = zeros(Nz,Nx,Nt);
indmt = zeros(Nz,Nx,Nt);
for it=1:Nt
    indpt(:,:,it) = min(1,sum(indp(:,:,tp==it),3));
    indmt(:,:,it) = min(1,sum(indm(:,:,tp==it),3));
end

end

function Cq = get_conc(XX,ZZ,xp,zp,tp,rp,indp,L,D,Nz,Nx,Nt,Np)

Cq = zeros(Nz,Nx,Nt);
for ip=1:sum(Np)
    ff = exp(1/6);
    rr = (1.25*rp(tp(ip))).^2;
    Cq(:,:,tp(ip)) = min(1,Cq(:,:,tp(ip)) + ff.*exp(-(XX-xp(ip)  ).^2./rr) .* exp(-(ZZ-zp(ip)  ).^2./rr));
    Cq(:,:,tp(ip)) = min(1,Cq(:,:,tp(ip)) + ff.*exp(-(XX-xp(ip)-L).^2./rr) .* exp(-(ZZ-zp(ip)  ).^2./rr));
    Cq(:,:,tp(ip)) = min(1,Cq(:,:,tp(ip)) + ff.*exp(-(XX-xp(ip)+L).^2./rr) .* exp(-(ZZ-zp(ip)  ).^2./rr));
    Cq(:,:,tp(ip)) = min(1,Cq(:,:,tp(ip)) + ff.*exp(-(XX-xp(ip)  ).^2./rr) .* exp(-(ZZ-zp(ip)-D).^2./rr));
    Cq(:,:,tp(ip)) = min(1,Cq(:,:,tp(ip)) + ff.*exp(-(XX-xp(ip)  ).^2./rr) .* exp(-(ZZ-zp(ip)+D).^2./rr));
    Cq(:,:,tp(ip)) = min(1,Cq(:,:,tp(ip)) + ff.*exp(-(XX-xp(ip)-L).^2./rr) .* exp(-(ZZ-zp(ip)-D).^2./rr));
    Cq(:,:,tp(ip)) = min(1,Cq(:,:,tp(ip)) + ff.*exp(-(XX-xp(ip)+L).^2./rr) .* exp(-(ZZ-zp(ip)-D).^2./rr));
    Cq(:,:,tp(ip)) = min(1,Cq(:,:,tp(ip)) + ff.*exp(-(XX-xp(ip)-L).^2./rr) .* exp(-(ZZ-zp(ip)+D).^2./rr));
    Cq(:,:,tp(ip)) = min(1,Cq(:,:,tp(ip)) + ff.*exp(-(XX-xp(ip)+L).^2./rr) .* exp(-(ZZ-zp(ip)+D).^2./rr));
end
for it=1:length(Np)
    Cq(:,:,it) = Cq(:,:,it).*(1-sum(indp,3)+indp(:,:,it));
end
Cq = Cq./max(1,sum(Cq,3));
end