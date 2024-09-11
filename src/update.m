%%*****  UPDATE PARAMETERS & AUXILIARY FIELDS  ****************************

if step>0

% advect particle fields
dCdt = - advect(C,U(2:end-1,:),W(:,2:end-1),h,{ADVN,''},[1,2],BCA);

% residual of particle field update
res_C = (a1*C-a2*Co-a3*Coo)/dt - (b1*dCdt + b2*dCdto + b3*dCdtoo);

% % update particle fields
% C = (a2*Co + a3*Coo + (b1*dCdt + b2*dCdto + b3*dCdtoo)*dt)/a1;

% semi-implicit update of major component density
upd_C = - alpha*res_C*dt/a1 + beta*upd_C;
C     = C + upd_C;

end

% update coefficient fields
eta = etam.*ones(Nz,Nx);
for it=1:Np
    eta = eta .* etap.^(C(:,:,it)>=0);
end

rho = rhom.*ones(Nz,Nx);
for it=1:Np
    rho = rho + (rhop(it)-rhom).*(C(:,:,it)>=0);
end

rhofz = (rho(icz(1:end-1),:) + rho(icz(2:end),:))/2;
rhofx = (rho(:,icx(1:end-1)) + rho(:,icx(2:end)))/2;
etacc = eta;
etaco = (eta(icz(1:end-1),icx(1:end-1)) .* eta(icz(1:end-1),icx(2:end-0)) ...
      .* eta(icz(2:end-0),icx(1:end-1)) .* eta(icz(2:end-0),icx(2:end-0))).^0.25;

% update time step
dta = h/2/max(abs([U(:);W(:)]));
dt  = min([1.1*dto,CFL*dta]);                                       % time step size

% update velocity divergence
Div_V = ddz(W(:,2:end-1),h) + ddx(U(2:end-1,:),h);                         % get velocity divergence

% update strain rates
exx = diff(U(2:end-1,:),1,2)./h - Div_V./3;                                % x-normal strain rate
ezz = diff(W(:,2:end-1),1,1)./h - Div_V./3;                                % z-normal strain rate
exz = 1/2.*(diff(U,1,1)./h+diff(W,1,2)./h);                                % shear strain rate

eII = (0.5.*(exx.^2 + ezz.^2 ...
       + 2.*(exz(1:end-1,1:end-1).^2+exz(2:end,1:end-1).^2 ...
       +     exz(1:end-1,2:end  ).^2+exz(2:end,2:end  ).^2)/4)).^0.5 + eps;

% update stresses
txx = etacc .* exx;                                                        % x-normal stress
tzz = etacc .* ezz;                                                        % z-normal stress
txz = etaco .* exz;                                                        % xz-shear stress

tII = (0.5.*(txx.^2 + tzz.^2 ...
       + 2.*(txz(1:end-1,1:end-1).^2+txz(2:end,1:end-1).^2 ...
       +     txz(1:end-1,2:end  ).^2+txz(2:end,2:end  ).^2)/4)).^0.5 + eps;