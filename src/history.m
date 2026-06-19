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

HST.P_mean(stp) = mean(P ,  'all');     % dynamic pressure mean
HST.P_std (stp) = std (P ,1,'all');     % dynamic pressure std

HST.Pm_mean(stp) = mean(Vm ,  'all');     % melt flow speed mean
HST.Pm_std (stp) = std (Vm ,1,'all');     % melt flow speed std

% get time average quantities over latter half of simulation
HST.Vm_tavg(stp) = mean(HST.Vm_rms(end-floor(stp/2):end));
HST.Vc_tavg(stp) = mean(HST.Vc_rms(end-floor(stp/2):end));
HST.P_tavg (stp) = mean(HST.P_mean(end-floor(stp/2):end));
HST.Pm_tavg(stp) = mean(HST.Pm_mean(end-floor(stp/2):end));

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

    HST.Pp_mean(stp,it) = mean(Pp(tp==it),  'all');   % particle pressure diff mean
    HST.Pp_std (stp,it) = std (Pp(tp==it),1,'all');   % particle pressure diff std

    dWp = Wp - HST.Wp_rms(stp,it);                    % particle deviation from rms z-speed
    dUp = Up - HST.Up_rms(stp,it);                    % particle deviation from rms x-speed
    dVp = sqrt(dWp.^2 + dUp.^2);                      % particle deviation from rms   speed

    HST.dWp_rms(stp,it) = rms(abs(dWp(tp==it)),'all'); % mean particle z-speed fluctuation
    HST.dUp_rms(stp,it) = rms(abs(dUp(tp==it)),'all'); % mean particle x-speed fluctuation
    HST.dVp_rms(stp,it) = rms(abs(dVp(tp==it)),'all'); % mean particle   speed fluctuation

    HST.DWp_rms(stp,it) = rms(DWp(tp==it),  'all');   % particle-melt z-speed difference mean
    HST.DWp_std(stp,it) = std(DWp(tp==it),1,'all');   % particle-melt z-speed difference std 

    HST.DVp_rms(stp,it) = rms(DVp(tp==it),  'all');   % particle-melt speed difference mean
    HST.DVp_std(stp,it) = std(DVp(tp==it),1,'all');   % particle-melt speed difference mean

    HST.DPp_mean(stp,it) = rms(DPp(tp==it),  'all');  % particle-melt pressure difference mean
    HST.DPp_std (stp,it) = std(DPp(tp==it),1,'all');  % particle-melt pressure difference mean

    % hindered settling diagnostics averaged over latter half of simulation
    HST.fp_tavg (stp,it) = sum(HST.fp_sum (end-floor(stp/2):end,it).*HST.time(end-floor(stp/2):end).') ./ sum(HST.time(end-floor(stp/2):end));
    HST.Wp_tavg (stp,it) = sum(HST.Wp_rms (end-floor(stp/2):end,it).*HST.time(end-floor(stp/2):end).') ./ sum(HST.time(end-floor(stp/2):end));
    HST.Vp_tavg (stp,it) = sum(HST.Vp_rms (end-floor(stp/2):end,it).*HST.time(end-floor(stp/2):end).') ./ sum(HST.time(end-floor(stp/2):end));
    HST.DWp_tavg(stp,it) = sum(HST.DWp_rms(end-floor(stp/2):end,it).*HST.time(end-floor(stp/2):end).') ./ sum(HST.time(end-floor(stp/2):end));
    HST.DVp_tavg(stp,it) = sum(HST.DVp_rms(end-floor(stp/2):end,it).*HST.time(end-floor(stp/2):end).') ./ sum(HST.time(end-floor(stp/2):end));
    HST.DPp_tavg(stp,it) = sum(HST.DPp_mean(end-floor(stp/2):end,it).*HST.time(end-floor(stp/2):end).') ./ sum(HST.time(end-floor(stp/2):end));

    % particle positions
    for ip = 1:Np(it)
        k = k+1;
        HST.pos(stp,k,1) = zp(k);
        HST.pos(stp,k,2) = xp(k);
    end
end

% record smoothed particle fractions
HST.chi_mean(stp,:) = mean(reshape(chi,Nz*Nx,Nt));
HST.chi_std (stp,:) = std (reshape(chi,Nz*Nx,Nt));
HST.chi_tavg(stp  ) = mean(HST.chi_mean(end-floor(stp/2):end));
