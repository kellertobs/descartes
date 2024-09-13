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

for ip = 1:length(xp)
    indp   = sqrt((XX-xp(ip)).^2 + (ZZ-zp(ip)).^2) < rp(tp(ip));
    Wp(ip) = sum(Wc.*indp,'all') ./ sum(indp,'all');
    Up(ip) = sum(Uc.*indp,'all') ./ sum(indp,'all');
    indm   = sqrt((XX-xp(ip)).^2 + (ZZ-zp(ip)).^2) >= rp(tp(ip)) ...
           & sqrt((XX-xp(ip)).^2 + (ZZ-zp(ip)).^2) <  rp(tp(ip))*4 ...
           & all(C<0,3);
    Wm(ip) = sum(Wc.*indm,'all') ./ sum(indm,'all');
    Um(ip) = sum(Uc.*indm,'all') ./ sum(indm,'all');
end
DWp = Wp - Wm;

res_zp = (a1*zp-a2*zpo-a3*zpoo)/dt - (b1*Wp + b2*Wpo + b3*Wpoo);
res_xp = (a1*xp-a2*xpo-a3*xpoo)/dt - (b1*Up + b2*Upo + b3*Upoo);

upd_xp = - alpha*res_xp*dt/a1 + beta*upd_xp;
xp     = xp + upd_xp;

upd_zp = - alpha*res_zp*dt/a1 + beta*upd_zp;
zp     = zp + upd_zp;

xpo(xp>L) = xpo(xp>L)-L; xpoo(xp>L) = xpoo(xp>L)-L; xp(xp>L) = xp(xp>L)-L;
xpo(xp<0) = xpo(xp<0)+L; xpoo(xp<0) = xpoo(xp<0)+L; xp(xp<0) = xp(xp<0)+L;
zpo(zp>D) = zpo(zp>D)-D; zpoo(zp>D) = zpoo(zp>D)-D; zp(zp>D) = zp(zp>D)-D;
zpo(zp<0) = zpo(zp<0)+D; zpoo(zp<0) = zpoo(zp<0)+D; zp(zp<0) = zp(zp<0)+D;

end

% update coefficient fields
eta = etam.*ones(Nz,Nx);
for it=1:Nt
    eta = eta .* etap.^(C(:,:,it)>=0);
end

rho = rhom.*ones(Nz,Nx);
for it=1:Nt
    rho = rho + (rhop(it)-rhom).*(C(:,:,it)>=0);
end

rhofz = (rho(icz(1:end-1),:) + rho(icz(2:end),:))/2;
rhofx = (rho(:,icx(1:end-1)) + rho(:,icx(2:end)))/2;

rhoW = rhofz.*W(:,2:end-1); 
rhoU = rhofx.*U(2:end-1,:); 

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

% update Reynolds number
Re = Vel.*rho.*D./eta;
