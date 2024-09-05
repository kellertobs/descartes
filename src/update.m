%%*****  UPDATE PARAMETERS & AUXILIARY FIELDS  ****************************

Uf = interp2(XXu,ZZu,U,XXuf,ZZuf,'linear',0);
Wf = interp2(XXw,ZZw,W,XXwf,ZZwf,'linear',0);

% advance particle fields
dCdt = - advect(C,Uf(2:end-1,:),Wf(:,2:end-1),h,{ADVN,''},[1,2],BCA);

% residual of major component evolution
res_C = (a1*C-a2*Co-a3*Coo)/dt - (b1*dCdt + b2*dCdto + b3*dCdtoo);

% semi-implicit update of major component density
upd_C = max(-C, - alpha*res_C*dt/a1 + beta*upd_C );
C     = C + upd_C;

% update coefficient fields
eta = etam.*ones(Nzf,Nxf);
for it=1:Np
    eta = eta .* etap.^C(:,:,it);
end

rho = rhom.*ones(Nzf,Nxf);
for it=1:Np
    rho = rho + (rhop(it)-rhom).*C(:,:,it);
end

rhofz = interp2(XXf,ZZf,rho,XXw(2:end-1,2:end-1),ZZw(2:end-1,2:end-1));
rhofx = interp2(XXf,ZZf,rho,XXu(2:end-1,2:end-1),ZZu(2:end-1,2:end-1));
rhocc =     interp2(XXf,ZZf,      rho ,XX  ,ZZ  );
etacc = 10.^interp2(XXf,ZZf,log10(eta),XX  ,ZZ  );
etaco = 10.^interp2(XXfg,ZZfg,log10(eta([end,1:end,1],[end,1:end,1])),XXco,ZZco,'linear');

% update time step
dta = h/2   /max(abs([U(:);W(:)]));
dt  = min([1.01*dto,CFL*dta]);                                       % time step size

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