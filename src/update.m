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
    % UPDATE PARTICLE VELOCITIES (Wp, Up) AND MELT VELOCITIES (Wm, Um)
    % ---------------------------------------------------------------------

    [Wp,Up,Wm,Um,DWp,DUp] = phase_vel(Wc,Uc,xp,zp,rp,tp,Np,XX,ZZ,L,D);


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


% update particle indicator fields
indp = zeros(Nz, Nx, Nt);
for ip = 1:Np
    % Update the particle indicator field for a given particle
    indp(:,:,tp(ip)) = min(1, indp(:,:,tp(ip)) + double(sqrt((XX - xp(ip)).^2 + (ZZ - zp(ip)).^2) < rp(tp(ip))));
    % Apply periodic boundary conditions in X and Z
    indp(:,:,tp(ip)) = min(1, indp(:,:,tp(ip)) + double(sqrt((XX - xp(ip) - L).^2 + (ZZ - zp(ip)).^2) < rp(tp(ip))));
    indp(:,:,tp(ip)) = min(1, indp(:,:,tp(ip)) + double(sqrt((XX - xp(ip) + L).^2 + (ZZ - zp(ip)).^2) < rp(tp(ip))));
    indp(:,:,tp(ip)) = min(1, indp(:,:,tp(ip)) + double(sqrt((XX - xp(ip)).^2 + (ZZ - zp(ip) - D).^2) < rp(tp(ip))));
    indp(:,:,tp(ip)) = min(1, indp(:,:,tp(ip)) + double(sqrt((XX - xp(ip)).^2 + (ZZ - zp(ip) + D).^2) < rp(tp(ip))));
end

Cq = get_conc(XX,ZZ,xp,zp,tp,rp,indp,L,D,Nz,Nx,Nt,Np);


% -------------------------------------------------------------------------
% UPDATE COEFFICIENT FIELDS (ETA, RHO)
% -------------------------------------------------------------------------

% Update eta (viscosity field)
eta = etam .* ones(Nz, Nx);
for it = 1:Nt
    eta = eta .* etap.^(indp(:,:,it) > 0);
end

% Update rho (density field)
rho = rhom .* ones(Nz, Nx);
for it = 1:Nt
    rho = rho + (rhop(it) - rhom) .* (indp(:,:,it) > 0);
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

function [Wp,Up,Wm,Um,DWp,DUp] = phase_vel(Wc,Uc,xp,zp,rp,tp,Np,XX,ZZ,L,D)

Wp = zeros(Np, 1);  Up = zeros(Np, 1);
Wm = zeros(Np, 1);  Um = zeros(Np, 1);

for ip = 1:Np
    indp   =  sqrt((XX-xp(ip)  ).^2 + (ZZ-zp(ip)  ).^2) < rp(tp(ip)) ...
           |  sqrt((XX-xp(ip)-L).^2 + (ZZ-zp(ip)  ).^2) < rp(tp(ip)) ...
           |  sqrt((XX-xp(ip)+L).^2 + (ZZ-zp(ip)  ).^2) < rp(tp(ip)) ...
           |  sqrt((XX-xp(ip)  ).^2 + (ZZ-zp(ip)-D).^2) < rp(tp(ip)) ...
           |  sqrt((XX-xp(ip)  ).^2 + (ZZ-zp(ip)+D).^2) < rp(tp(ip)) ...
           |  sqrt((XX-xp(ip)-L).^2 + (ZZ-zp(ip)-D).^2) < rp(tp(ip)) ...
           |  sqrt((XX-xp(ip)+L).^2 + (ZZ-zp(ip)-D).^2) < rp(tp(ip)) ...
           |  sqrt((XX-xp(ip)-L).^2 + (ZZ-zp(ip)+D).^2) < rp(tp(ip)) ...
           |  sqrt((XX-xp(ip)+L).^2 + (ZZ-zp(ip)+D).^2) < rp(tp(ip));
    Wp(ip) = sum(Wc.*indp,'all') ./ sum(indp,'all');
    Up(ip) = sum(Uc.*indp,'all') ./ sum(indp,'all');

    indm   = (sqrt((XX-xp(ip)  ).^2 + (ZZ-zp(ip)).^2  ) >= rp(tp(ip))    ...
           |  sqrt((XX-xp(ip)-L).^2 + (ZZ-zp(ip)).^2  ) >= rp(tp(ip))    ...
           |  sqrt((XX-xp(ip)+L).^2 + (ZZ-zp(ip)).^2  ) >= rp(tp(ip))    ...
           |  sqrt((XX-xp(ip)  ).^2 + (ZZ-zp(ip)-D).^2) >= rp(tp(ip))    ...
           |  sqrt((XX-xp(ip)  ).^2 + (ZZ-zp(ip)+D).^2) >= rp(tp(ip))    ...
           |  sqrt((XX-xp(ip)-L).^2 + (ZZ-zp(ip)-D).^2) >= rp(tp(ip))    ...
           |  sqrt((XX-xp(ip)+L).^2 + (ZZ-zp(ip)-D).^2) >= rp(tp(ip))    ...
           |  sqrt((XX-xp(ip)-L).^2 + (ZZ-zp(ip)+D).^2) >= rp(tp(ip))    ...
           |  sqrt((XX-xp(ip)+L).^2 + (ZZ-zp(ip)+D).^2) >= rp(tp(ip)))   ...
           & (sqrt((XX-xp(ip)  ).^2 + (ZZ-zp(ip)).^2  ) <  rp(tp(ip))*4  ...
           |  sqrt((XX-xp(ip)-L).^2 + (ZZ-zp(ip)).^2  ) <  rp(tp(ip))*4  ...
           |  sqrt((XX-xp(ip)+L).^2 + (ZZ-zp(ip)).^2  ) <  rp(tp(ip))*4  ...
           |  sqrt((XX-xp(ip)  ).^2 + (ZZ-zp(ip)-D).^2) <  rp(tp(ip))*4  ...
           |  sqrt((XX-xp(ip)  ).^2 + (ZZ-zp(ip)+D).^2) <  rp(tp(ip))*4  ...
           |  sqrt((XX-xp(ip)-L).^2 + (ZZ-zp(ip)-D).^2) <  rp(tp(ip))*4  ...
           |  sqrt((XX-xp(ip)+L).^2 + (ZZ-zp(ip)-D).^2) <  rp(tp(ip))*4  ...
           |  sqrt((XX-xp(ip)-L).^2 + (ZZ-zp(ip)+D).^2) <  rp(tp(ip))*4  ...
           |  sqrt((XX-xp(ip)+L).^2 + (ZZ-zp(ip)+D).^2) <  rp(tp(ip))*4) ...
           &  all(indp==0,3);
    Wm(ip) = sum(Wc.*indm,'all') ./ sum(indm,'all');
    Um(ip) = sum(Uc.*indm,'all') ./ sum(indm,'all');
end

% Particle segregation speeds
DWp = Wp - Wm;
DUp = Up - Um;

end

function Cq = get_conc(XX,ZZ,xp,zp,tp,rp,indp,L,D,Nz,Nx,Nt,Np)

Cq = zeros(Nz,Nx,Nt);
for ip=1:Np
    ff = exp(1/6);
    rr = (1.5*rp(tp(ip))).^2;
    Cq(:,:,tp(ip)) = min(1,Cq(:,:,tp(ip)) + ff.*exp(-(XX-xp(ip)  ).^2./rr) .* exp(-(ZZ-zp(ip)  ).^2./rr));
    Cq(:,:,tp(ip)) = min(1,Cq(:,:,tp(ip)) + ff.*exp(-(XX-xp(ip)-L).^2./rr) .* exp(-(ZZ-zp(ip)  ).^2./rr));
    Cq(:,:,tp(ip)) = min(1,Cq(:,:,tp(ip)) + ff.*exp(-(XX-xp(ip)+L).^2./rr) .* exp(-(ZZ-zp(ip)  ).^2./rr));
    Cq(:,:,tp(ip)) = min(1,Cq(:,:,tp(ip)) + ff.*exp(-(XX-xp(ip)  ).^2./rr) .* exp(-(ZZ-zp(ip)-D).^2./rr));
    Cq(:,:,tp(ip)) = min(1,Cq(:,:,tp(ip)) + ff.*exp(-(XX-xp(ip)  ).^2./rr) .* exp(-(ZZ-zp(ip)+D).^2./rr));
    Cq(:,:,tp(ip)) = min(1,Cq(:,:,tp(ip)) + ff.*exp(-(XX-xp(ip)-L).^2./rr) .* exp(-(ZZ-zp(ip)-D).^2./rr));
    Cq(:,:,tp(ip)) = min(1,Cq(:,:,tp(ip)) + ff.*exp(-(XX-xp(ip)+L).^2./rr) .* exp(-(ZZ-zp(ip)-D).^2./rr));
    Cq(:,:,tp(ip)) = min(1,Cq(:,:,tp(ip)) + ff.*exp(-(XX-xp(ip)-L).^2./rr) .* exp(-(ZZ-zp(ip)+D).^2./rr));
    Cq(:,:,tp(ip)) = min(1,Cq(:,:,tp(ip)) + ff.*exp(-(XX-xp(ip)+L).^2./rr) .* exp(-(ZZ-zp(ip)+D).^2./rr));
    Cq(:,:,tp(ip)) = Cq(:,:,tp(ip)).*(1-sum(indp,3)+indp(:,:,tp(ip)));
end

end