stp = max(1,step);

% record time
HST.time(stp) = time;

% record mixture and melt velocities
HST.Wc_rms(stp) = rms(W  ,  'all');     % convection z-speed mean
HST.Wc_std(stp) = std(W  ,1,'all');     % convection z-speed std

HST.Uc_rms(stp) = rms(U  ,  'all');     % convection z-speed mean
HST.Uc_std(stp) = std(U  ,1,'all');     % convection z-speed std

HST.Vc_rms(stp) = rms(Vel,  'all');     % convection speed mean
HST.Vc_std(stp) = std(Vel,1,'all');     % convection speed std

HST.Wm_rms(stp) = rms(Wm ,  'all');     % melt flow z-speed mean
HST.Wm_std(stp) = std(Wm ,1,'all');     % melt flow z-speed std

HST.Um_rms(stp) = rms(Um ,  'all');     % melt flow x-speed mean
HST.Um_std(stp) = std(Um ,1,'all');     % melt flow x-speed std

HST.Vm_rms(stp) = rms(Vm ,  'all');     % melt flow speed mean
HST.Vm_std(stp) = std(Vm ,1,'all');     % melt flow speed std

HST.Vm_tavg(stp) = mean(HST.Vm_rms(end-floor(stp/2):end));
HST.Vc_tavg(stp) = mean(HST.Vc_rms(end-floor(stp/2):end));

% record particle velocities
k = 0;
for it = 1:Nt
    HST.fp_sum(stp,it) = sum(fp(tp==it),  'all');     % phase fractions of sum of particles
    HST.fp_std(stp,it) = std(fp(tp==it),1,'all');     % phase fractions standard deviations (std)

    HST.Wp_rms(stp,it) = rms(Wp(tp==it),  'all');     % particle settling z-speed mean
    HST.Wp_std(stp,it) = std(Wp(tp==it),1,'all');     % particle settling z-speed std

    HST.Up_rms(stp,it) = rms(Up(tp==it),  'all');     % particle settling x-speed mean
    HST.Up_std(stp,it) = std(Up(tp==it),1,'all');     % particle settling x-speed std

    HST.Vp_rms(stp,it) = rms(Vp(tp==it),  'all');     % particle settling speed mean
    HST.Vp_std(stp,it) = std(Vp(tp==it),1,'all');     % particle settling speed std

    dWp = Wp - HST.Wp_rms(stp,it);
    dUp = Up - 0;
    dVp = sqrt(dWp.^2 + dUp.^2);

    HST.dWp_rms(stp,it) = rms(abs(dWp(tp==it)),'all'); % mean particle z-speed fluctuation
    HST.dUp_rms(stp,it) = rms(abs(dUp(tp==it)),'all'); % mean particle x-speed fluctuation
    HST.dVp_rms(stp,it) = rms(abs(dVp(tp==it)),'all'); % mean particle   speed fluctuation

    HST.DWp_rms(stp,it) = rms(DWp(tp==it),  'all');   % particle-melt z-speed difference mean
    HST.DWp_std(stp,it) = std(DWp(tp==it),1,'all');   % particle-melt z-speed difference std 

    HST.DVp_rms(stp,it) = rms(DVp(tp==it),  'all');   % particle-melt speed difference mean
    HST.DVp_std(stp,it) = std(DVp(tp==it),1,'all');   % particle-melt speed difference mean

    % hindered settling diagnostics averaged over latter half of simulation
    HST.fp_tavg (stp,it) = mean(HST.fp_sum (end-floor(stp/2):end,it));
    HST.Vp_tavg (stp,it) = mean(HST.Vp_rms (end-floor(stp/2):end,it));
    HST.DWp_tavg(stp,it) = mean(HST.DWp_rms(end-floor(stp/2):end,it));
    HST.DVp_tavg(stp,it) = mean(HST.DVp_rms(end-floor(stp/2):end,it));

    % particle positions
    for ip = 1:Np(it)
        k = k+1;
        HST.pos(stp,k,1) = zp(k);
        HST.pos(stp,k,2) = xp(k);
    end
end

