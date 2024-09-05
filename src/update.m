%%*****  UPDATE PARAMETERS & AUXILIARY FIELDS  ****************************
tic;

% update time step
dtk = (h/2)^2/max([kc(:);kwm(:);kwx(:);kwf(:);(kT(:)+ks(:).*T(:))./rho(:)./cP(:)])*0.9; % diffusive time step size
dta =  h/2   /max(abs([Um(:).* mux(:);Wm(:).* muz(:); ...  % advective time step size
    Ux(:).*chix(:);Wx(:).*chiz(:); ...
    Uf(:).*phix(:);Wf(:).*phiz(:)]+eps));
dtc = 0.01./max(abs([advn_X(:)./rho(:);advn_M(:)./rho(:);advn_F(:)./rho(:)]));
dt  = min([1.01*dto,min([dtk,CFL*dta,dtc]),dtmax,tau_T/100]);                         % time step size

% update density and viscosity fields
rho  = rhoc.*c;
eta  = etac.^c;

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
txx = eta   .* exx;                                                        % x-normal stress
tzz = eta   .* ezz;                                                        % z-normal stress
txz = etaco .* exz;                                                        % xz-shear stress

tII = (0.5.*(txx.^2 + tzz.^2 ...
       + 2.*(txz(1:end-1,1:end-1).^2+txz(2:end,1:end-1).^2 ...
       +     txz(1:end-1,2:end  ).^2+txz(2:end,2:end  ).^2)/4)).^0.5 + eps;

UDtime = UDtime + toc;
